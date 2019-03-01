"""
script to measure galaxy shapes in Illustris
"""

from __future__ import print_function, division
import numpy as np
from astropy.table import Table, Column
import time
from astropy.io import ascii
import sys

from illustris_python.snapshot import loadHalo, snapPath, loadSubhalo
from illustris_python.groupcat import gcPath, loadHalos, loadSubhalos
from inertia_tensors import inertia_tensors, reduced_inertia_tensors

from halotools.utils import normalized_vectors


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


def galaxy_ang_mom(gal_id, basePath, snapNum, Lbox, reduced=True):
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

    # load galaxy position (most bound particle)
    gal_positions = loadSubhalos(basePath, snapNum, fields=['SubhaloPos'])/1000.0
    gal_position = gal_positions[gal_id]

    # half mass radius
    gal_rhalfs = loadSubhalos(basePath, snapNum, fields=['SubhaloHalfmassRadType'])[:,4]/1000.0
    gal_rhalf = gal_rhalfs[gal_id]

    # load stellar particles
    ptcl_coords = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['Coordinates'])/1000.0
    ptcl_masses = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['Masses'])*10.0**10
    ptcl_vels = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['Velocities'])
    sf_time = loadSubhalo(basePath, snapNum, gal_id, 4, fields=['GFM_StellarFormationTime'])
    is_a_star = (sf_time>=0.0) # don't use wind particles

    # account for PBCs
    dx = ptcl_coords[:,0] - gal_position[0]
    dy = ptcl_coords[:,1] - gal_position[1]
    dz = ptcl_coords[:,2] - gal_position[2]

    mask = (dx > Lbox[0]/2.0)
    dx[mask] = dx[mask] - Lbox[0]
    mask = (dx < -Lbox[0]/2.0)
    dx[mask] = dx[mask] + Lbox[0]

    mask = (dy > Lbox[1]/2.0)
    dy[mask] = dy[mask] - Lbox[1]
    mask = (dy < -Lbox[1]/2.0)
    dy[mask] = dy[mask] + Lbox[1]

    mask = (dz > Lbox[2]/2.0)
    dz[mask] = dz[mask] - Lbox[2]
    mask = (dz < -Lbox[2]/2.0)
    dz[mask] = dz[mask] + Lbox[2]

    ptcl_coords = np.vstack((dx,dy,dz)).T

    r = np.sqrt(np.sum(ptcl_coords**2, axis=1))/gal_rhalf
    mask = (r<=10.0) & (is_a_star)

    L = specific_angular_momentum(ptcl_coords[mask], ptcl_vels[mask], ptcl_masses[mask])

    mag_L = np.sqrt(np.sum(L**2,axis=-1))

    return mag_L, L/mag_L


def main():

    basePath = '/Volumes/G-RAID/simulations/unprocessed/Illustris/Illustris-1'
    snapNum = 135
    m_dm = 6.3*10**6.0 #dark matter particle mass
    litte_h = 0.704
    Lbox = np.array([75.0]*3)

    # make selection
    x = loadSubhalos(basePath, snapNum)

    gal_ids = np.arange(0,len(x['SubhaloGrNr']))
    gal_stellar_mass = x['SubhaloMassInRadType'][:,4] # mass within 2*R_half

    mask = np.log10(gal_stellar_mass*10**(10)/litte_h)>=9.0

    gal_ids = gal_ids[mask]
    Ngals = len(gal_ids)
    print("number of galaxies in selection: {0}".format(Ngals))

    L = np.zeros(Ngals)
    Lx = np.zeros(Ngals)
    Ly = np.zeros(Ngals)
    Lz = np.zeros(Ngals)
    for i, gal_id in enumerate(gal_ids):
        print(1.0*i/Ngals)
        m, vec = galaxy_ang_mom(gal_id, basePath, snapNum, Lbox, reduced=True)
        L[i] = m
        Lx[i] = vec[0]
        Ly[i] = vec[1]
        Lz[i] = vec[2]

    # save measurements
    fpath = './data/'
    fname = 'galaxy_angular_momentum.dat'
    ascii.write([gal_ids, L, Lx, Ly, Lz],
                fpath+fname,
                names=['gal_id', 'L', 'Lx', 'Ly','Lz'],
                overwrite=True)






if __name__ == '__main__':
    main()
