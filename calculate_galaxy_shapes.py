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


def galaxy_shape(gal_id, basePath, snapNum, Lbox, reduced=True):
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
    mask = (r<=2.0) & (is_a_star)

    if reduced:
        I = reduced_inertia_tensors(ptcl_coords[mask], ptcl_masses[mask])
    else:
        I = inertia_tensors(ptcl_coords[mask], ptcl_masses[mask])

    evals, evecs = np.linalg.eigh(I)
    evals = np.sqrt(evals)

    return evals[0], evecs[0]


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

    a = np.zeros(Ngals)
    b = np.zeros(Ngals)
    c = np.zeros(Ngals)
    av = np.zeros((Ngals,3))
    bv = np.zeros((Ngals,3))
    cv = np.zeros((Ngals,3))
    for i, gal_id in enumerate(gal_ids):
        print(1.0*i/Ngals)
        evals, evecs = galaxy_shape(gal_id, basePath, snapNum, Lbox, reduced=True)
        a[i] = evals[2]
        b[i] = evals[1]
        c[i] = evals[0]
        av[i,:] = evecs[2]
        bv[i,:] = evecs[1]
        cv[i,:] = evecs[0]

    # save measurements
    fpath = './data/'
    fname = 'galaxy_shapes_1.dat'
    ascii.write([gal_ids, a, b, c,
                 av[:,0], av[:,1], av[:,2],
                 bv[:,0], bv[:,1], bv[:,2],
                 cv[:,0], cv[:,1], cv[:,2]],
                fpath+fname,
                names=['gal_id', 'a', 'b', 'c',
                       'av_x', 'av_y','av_z',
                       'bv_x', 'bv_y','bv_z',
                       'cv_x', 'cv_y','cv_z'],
                overwrite=True)






if __name__ == '__main__':
    main()
