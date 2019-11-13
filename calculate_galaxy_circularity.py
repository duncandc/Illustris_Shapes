"""
script to measure galaxy angular momentum and circularity
"""

from __future__ import print_function, division
import numpy as np
import sys

# I/O
from astropy.io import ascii

# illustris python functions
from illustris_python.snapshot import loadSubhalo
from illustris_python.groupcat import loadSubhalos, loadHeader

# utilities
from scipy import interpolate
from rotations.vector_utilities import normalized_vectors

# Illustris simulation properties
from simulation_props import sim_prop_dict

# gravitational constant
from astropy.constants import G
G = G.to('Mpc km^2/(Msun s^2)').value

# progress bar
from tqdm import tqdm

def particle_cos_alpha(coords, vels, L):
    """
    cosine of the angle between the orbital angular momentum
    of the particle and that of the galaxy system

    Parameters
    ----------
    coords :
        normalized particle positions

    vels : array_like
        normalized particle velocities

    L : array_like
        array of 3D angular momentum vector

    Returns
    -------
    cos_alpha : array_like
       array of cosine(alpha) values
    """

    # specific angular momomentum of particles
    j = normalized_vectors(np.cross(coords, vels))
    l = normalized_vectors(L)[0]

    return np.dot(j, l)


def particle_cos_beta(coords, L):
    """
    cosine of the angle between the radial vector
    and the spin axis of the galaxy

    Parameters
    ----------
    coords :
        normalized particle positions

    L : array_like
        array of 3D angular momentum vector

    Returns
    -------
    cos_alpha : array_like
       array of cosine(alpha) values
    """

    # specific angular momomentum of particles
    r = normalized_vectors(coords)
    l = normalized_vectors(L)[0]

    return np.dot(r, l)


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


def disk_fraction(masses, epsilons, disk_threshold):
    """
    mass fraction of particles in disk

    Parameters
    ----------
    masses : array_like
        array of particle masses

    epsilons : array_like
        array of particle circularities

    disk_threshold : float
        circulalrity threshold for disk

    Returns
    -------
    f_disk : float
        fraction of mass in disk
    """

    # define disk particles as those where epsilon > e_thresh
    mask = (epsilons >= disk_threshold)

    return np.sum(masses[mask])/np.sum(masses)


def circular_velocity(r, coords, masses):
    """
    circular velocity at radial position r

    Parmaters
    ---------
    r : array_like
        radial positions

    coords : array_like
        array of massive particle coordinates (centered on the galaxy)

    masses : array_like
        array of massive particle masses

    Returns
    -------
    v_circ : array_like
        array of circular velocities
    """

    # radial position of massive particles
    ptcl_r = np.sqrt(np.sum(coords*coords, axis=-1))

    # sort by radial position
    r_sort_inds = np.argsort(ptcl_r)
    ptcl_r = ptcl_r[r_sort_inds]
    masses = masses[r_sort_inds]

    # mass internal to r
    m_within_r = np.cumsum(masses)

    # circular velocity
    vc = np.sqrt(G*m_within_r/ptcl_r)

    # interpolate circular velocity curve
    f = interpolate.interp1d(ptcl_r, vc, kind='linear', fill_value='extrapolate')

    return f(r)


def loadSubhaloAll(basePath, snapNum, gal_id, field, m_dm=None):
    """
    Load particles/cells (of all types) for a single field.

    Parameters
    ----------
    basePath :

    snapNum :

    gal_id : int

    field : string

    m_dm : float
        dark matter particle mass

    Returns
    -------
    result : numpy.array

    Notes
    -----
    This is a simple wrapper around the illustris_python.loadSubhalo() function.
    """

    ptypes = [0,1,4,5]

    # must hard code in various behaviour for different fields
    fields = ['Coordinates', 'Masses']
    if field not in fields:
        msg = ('field must be one of:', fields)
        raise ValueError(msg)

    # since all dm particle masses are the same,
    # this quantity isn't stored in the table.
    if (field=='Masses') & (m_dm is None):
        msg = ('dm particle mass must be passed if field==Masses')
        raise ValueError(msg)

    data = []
    for ptype in ptypes:

        # create dm mass array
        if (ptype==1) & (field=='Masses'):
            x = loadSubhalo(basePath, snapNum, gal_id, ptype, fields=['ParticleIDs'])
            # if no particles are found, loadSubhalo returns {'count': 0}
            if type(x) is dict:
                x = []
            x = np.array([m_dm]*len(x), dtype='float32')/10.0**10
        else:
            x = loadSubhalo(basePath, snapNum, gal_id, ptype, fields=[field])

        # if no particles are found, loadSubhalo returns {'count': 0}
        if type(x) is dict:
            pass
        else:
            data.append(x)

    if field in [fields[1]]:
        return np.hstack(tuple(data))
    else:
        return np.vstack(tuple(data))


def circularity(coords, vels, masses, v_circs, L):
    """
    circularity parameter
    """

    r = np.sqrt(np.sum(coords*coords, axis=-1))

    # specific angular momentum of a circular orbit at r
    j_circ = r*v_circs

    # specific angular momomentum of particles
    j = np.cross(coords, vels)

    # component along system angular momentum vector
    j_z = np.dot(j, L/np.sqrt(np.sum(L*L, axis=-1)))

    # circulatiry parameter
    epsilon = j_z/j_circ

    return epsilon


def format_particles(center, coords, Lbox):
    """
    center the particle coordinates on (0,0,0) and account for PBCs

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


def particle_selection(gal_id, ptcl_coords, galaxy_table, basePath, snapNum, radial_mask=True, num_r_half=10):
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
    gal_rhalfs = galaxy_table['SubhaloHalfmassRadType'][:,4]/1000.0
    gal_rhalf = gal_rhalfs[gal_id]

    # use only particles within 2 * R_half
    r = np.sqrt(np.sum(ptcl_coords**2, axis=1))/gal_rhalf
    if radial_mask:
        radial_mask = (r<=num_r_half)
    else:
        radial_mask = np.array([True]*len(star_mask))

    return (radial_mask) & (star_mask)


def galaxy_center(gal_id, galaxy_table):
    """
    Return the coordinates of the center of the galaxy
    """

    # load position of the most bound particle (of any type)
    coords = galaxy_table['SubhaloPos']/1000.0
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


def galaxy_circularity(gal_id, galaxy_table, basePath, snapNum, Lbox, num_r_half, m_dm):
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
    gal_position = galaxy_center(gal_id, galaxy_table)

    # load stellar particle positions and masses
    ptcl_coords = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['Coordinates'])/1000.0
    ptcl_vels = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['Velocities'])
    ptcl_masses = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['Masses'])*10.0**10

    # load all particles positions and masses
    all_ptcl_coords = loadSubhaloAll(basePath, snapNum, gal_id, field='Coordinates')/1000.0
    all_ptcl_masses = loadSubhaloAll(basePath, snapNum, gal_id, field='Masses', m_dm=m_dm)*10.0**10

    # center and account for PBCs
    ptcl_coords = format_particles(gal_position, ptcl_coords, Lbox)
    all_ptcl_coords = format_particles(gal_position, all_ptcl_coords, Lbox)

    # use center of mass velocity to subtract bulk velocity
    ptcl_vels = format_velocities(ptcl_vels, ptcl_masses, basePath, snapNum)

    # use a subset of particles for stellar properties
    ptcl_mask = particle_selection(gal_id, ptcl_coords, galaxy_table, basePath, snapNum, radial_mask=True, num_r_half=num_r_half)

    # specific angular momentum of the stellar system (use selected disk stars)
    L = specific_angular_momentum(ptcl_coords[ptcl_mask], ptcl_vels[ptcl_mask], ptcl_masses[ptcl_mask])

    # radial position of stellar particles
    r = np.sqrt(np.sum(ptcl_coords*ptcl_coords, axis=-1))

    # circular velocity at the radial position of stellar particles
    v_circs = circular_velocity(r, all_ptcl_coords, all_ptcl_masses)

    # circularity parameter for stellar particles
    epsilons = circularity(ptcl_coords, ptcl_vels, ptcl_masses, v_circs, L)

    # fraction of stellar mass with circularity > a threshold (use selected disk stars)
    f_disk = disk_fraction(ptcl_masses[ptcl_mask], epsilons[ptcl_mask], disk_threshold=0.7)

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

    # load galaxy table
    fields = ['SubhaloGrNr', 'SubhaloMassInRadType', 'SubhaloPos', 'SubhaloHalfmassRadType']
    galaxy_table = loadSubhalos(basePath, snapNum, fields=fields)

    Lx = np.zeros(Ngals)
    Ly = np.zeros(Ngals)
    Lz = np.zeros(Ngals)
    f_disk = np.zeros(Ngals)

    for i in tqdm(range(Ngals)):
        gal_id = gal_ids[i]
        f, vec = galaxy_circularity(gal_id, galaxy_table, basePath, snapNum, Lbox, num_r_half, m_dm)
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
