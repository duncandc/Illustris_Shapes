"""
script to measure galaxy angular momentum and circularity
"""

from __future__ import print_function, division
import numpy as np
import sys

# I/O
from astropy.table import Table, Column
from astropy.io import ascii

# illustris python functions
from illustris_python.snapshot import loadHalo, snapPath, loadSubhalo
from illustris_python.groupcat import gcPath, loadHalos, loadSubhalos, loadHeader

# utilities
from halotools.utils import normalized_vectors

# Illustris simulation properties
from simulation_props import sim_prop_dict

# gravitational constant
from astropy.constants import G
G = G.to('Mpc km^2/(Msun s^2)').value

# progress bar
from tqdm import tqdm


def specific_angular_momentum(x, v, m):
    """
    specific angular momentum of a group of particles

    Parameters
    ----------
    x : array_like
        array particle positions of shape (Nptcl, ndim)

    v : array_like
        array of particle velcities wth shape (Nptcl, ndim)

    m : array_like
        array of particle masses of shape (Nptcl,)

    Returns
    -------
    L : nump.array
        specific angular momentum vector
    """
    return (m[:,np.newaxis]*np.cross(x,v)).sum(axis=0)/m.sum()


def disk_fraction(coords, vels, masses, L, e_thresh=0.7):
    """
    fraction of stars oin disk
    """

    # circulatiry parameter
    epsilon = circularity(coords, vels, masses, L)

    mask = (epsilon>=e_thresh)

    return np.sum(masses[mask])/np.sum(masses)


def circularity(coords, vels, masses, L):
    """
    circularity parameter
    """

    r = np.sqrt(np.sum(coords*coords, axis=-1))
    r_sort_inds = np.argsort(r)

    r = r[r_sort_inds]
    coords = coords[r_sort_inds]
    vels = vels[r_sort_inds]
    masses = masses[r_sort_inds]

    # mass internal to r
    m_within_r = np.cumsum(masses)

    # specific angular momentum of a circular orbit at r
    j_circ = r*np.sqrt(G*m_within_r/r)

    # specific angular momomentum of particles
    j = np.cross(coords,vels)
    # component along system angular momentum vector
    j_z = np.dot(j,L/np.sqrt(np.sum(L*L,axis=-1)))

    # circulatiry parameter
    epsilon = j_z/j_circ

    return epsilon


def format_particles(center, coords, Lbox):
    """
    center the partivcle coordinates on (0,0,0) and account for PBCs

    Parameters
    ----------
    center : array_like
        array of shape (3,) storing the coordinates of the
        center of the particle distribution

    coords :  array_like
        array of shape (nptcls, 3) storing the coordinates of the
        particles

    Lbox : array_like
        length 3-array giving the box size in each dimension

    Returns
    -------
    coords : numpy.array
        array of shape (nptcls, 3) storing the centered particle coordinates
    """

    dx = coords[:,0] - center[0]
    dy = coords[:,1] - center[1]
    dz = coords[:,2] - center[2]

    # x-coordinate
    mask = (dx > Lbox[0]/2.0)
    dx[mask] = dx[mask] - Lbox[0]
    mask = (dx < -Lbox[0]/2.0)
    dx[mask] = dx[mask] + Lbox[0]

    # y-coordinate
    mask = (dy > Lbox[1]/2.0)
    dy[mask] = dy[mask] - Lbox[1]
    mask = (dy < -Lbox[1]/2.0)
    dy[mask] = dy[mask] + Lbox[1]

    # z-coordinate
    mask = (dz > Lbox[2]/2.0)
    dz[mask] = dz[mask] - Lbox[2]
    mask = (dz < -Lbox[2]/2.0)
    dz[mask] = dz[mask] + Lbox[2]

    # format coordinates
    coords = np.vstack((dx,dy,dz)).T

    return coords


def particle_selection(gal_id, ptcl_coords, basePath, snapNum, radial_mask=True, num_r_half=10):
    """
    Apply the selection criteria to galaxy particles

    Returns
    -------
    mask : numpy.array
       boolean array of shape (nptcls,)
    """

    # don't use wind particles
    sf_time = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['GFM_StellarFormationTime'])
    star_mask = (sf_time>=0.0) # don't use wind particles

    # get the half mass radius
    gal_rhalfs = loadSubhalos(basePath, snapNum, fields=['SubhaloHalfmassRadType'])[:,4]/1000.0
    gal_rhalf = gal_rhalfs[gal_id]

    # use only particles within 2 * R_half
    r = np.sqrt(np.sum(ptcl_coords**2, axis=1))/gal_rhalf
    if radial_mask:
        radial_mask = (r<=num_r_half)
    else:
        radial_mask = np.array([True]*len(star_mask))

    return (radial_mask) & (star_mask)


def galaxy_center(gal_id, basePath, snapNum):
    """
    Return the coordinates of the center of the galaxy
    """

    # load position of the most bound particle (of any type)
    coords = loadSubhalos(basePath, snapNum, fields=['SubhaloPos'])/1000.0
    coord = coords[gal_id]

    return coord


def galaxy_selection(min_mstar, basePath, snapNum):
    """
    make a cut on galaxy properties
    """

    # make selection
    galaxy_table = loadSubhalos(basePath, snapNum, fields=['SubhaloGrNr', 'SubhaloMassInRadType'])

    gal_ids = np.arange(0,len(galaxy_table['SubhaloGrNr']))

    # mass of stellar particles within 2*R_half
    mstar = galaxy_table['SubhaloMassInRadType'][:,4]
    mstar = mstar*10**10

    mask = (mstar >= min_mstar)

    return mask, gal_ids[mask]


def format_velocities(ptcl_vels, ptcl_masses, basePath, snapNum):
    """
    """

    # scale to get peculiar velocity
    header = loadHeader(basePath, snapNum)
    scale_factor = header['Time']
    ptcl_vels = ptcl_vels*np.sqrt(scale_factor)

    # remove center of mass velocity
    vel0 = np.sum(ptcl_masses[:,np.newaxis]*ptcl_vels, axis=0)/np.sum(ptcl_masses)
    ptcl_vels = ptcl_vels - vel0

    return ptcl_vels


def galaxy_circularity(gal_id, basePath, snapNum, Lbox, num_r_half):
    """
    Parameters
    ----------
    gal_id : int

    basepath : string

    snapNum : int

    Lbox : array_like

    reduced : bool

    Returns
    -------
    eig_vals, eig_vecs
    """

    # choose a 'center' for each galaxy
    gal_position = galaxy_center(gal_id, basePath, snapNum)

    # load stellar particle positions and masses
    ptcl_coords = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['Coordinates'])/1000.0
    ptcl_vels = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['Velocities'])
    ptcl_masses = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['Masses'])*10.0**10

    # center and account for PBCs
    ptcl_coords = format_particles(gal_position, ptcl_coords, Lbox)

    # use center of mass velocity to subtract bulk velocity
    ptcl_vels = format_velocities(ptcl_vels, ptcl_masses, basePath, snapNum)

    ptcl_mask = particle_selection(gal_id, ptcl_coords, basePath, snapNum, radial_mask=True, num_r_half=num_r_half)

    # specific angular momentum of the system
    L = specific_angular_momentum(ptcl_coords[ptcl_mask], ptcl_vels[ptcl_mask], ptcl_masses[ptcl_mask])

    # fraction of stellar mass with circularity > a threshold
    f_disk = disk_fraction(ptcl_coords[ptcl_mask], ptcl_vels[ptcl_mask], ptcl_masses[ptcl_mask], L, e_thresh=0.7)

    return f_disk, L


def main():

    if len(sys.argv)>1:
        sim_name = sys.argv[1]
        snapNum = int(sys.argv[2])
        num_r_half = float(sys.argv[3])
    else:
        sim_name = 'Illustris-1' # full physics high-res run
        snapNum = 135  # z=0
        num_r_half = 10.0

    # get simulation properties
    d = sim_prop_dict[sim_name]
    basePath = d['basePath']
    m_dm = d['m_dm']
    litte_h = d['litte_h']
    Lbox = d['Lbox']

    # make galaxy selection
    min_mstar = litte_h*10.0**9.0
    mask, gal_ids = galaxy_selection(min_mstar, basePath, snapNum)

    # number of galaxies in selection
    Ngals = len(gal_ids)
    print("number of galaxies in selection: {0}".format(Ngals))

    Lx = np.zeros(Ngals)
    Ly = np.zeros(Ngals)
    Lz = np.zeros(Ngals)
    f_disk = np.zeros(Ngals)

    for i in tqdm(range(Ngals)):
        gal_id = gal_ids[i]
        f, vec = galaxy_circularity(gal_id, basePath, snapNum, Lbox, num_r_half)
        f_disk[i] = f
        Lx[i] = vec[0]
        Ly[i] = vec[1]
        Lz[i] = vec[2]

    # save measurements
    fpath = './data/shape_catalogs/'
    fname = sim_name + '_' + str(snapNum) +'_galaxy_circularities_' + "{:.1f}".format(num_r_half) + '.dat'
    ascii.write([gal_ids, Lx, Ly, Lz, f_disk],
                 fpath+fname,
                 names=['gal_id', 'Lx', 'Ly', 'Lz', 'f_disk'],
                 overwrite=True)






if __name__ == '__main__':
    main()
