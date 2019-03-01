"""
script to make avalue added galaxy catalog for Illustris 
"""

from __future__ import print_function, division
import numpy as np
from astropy.table import Table, Column
import time
from astropy.io import ascii
import sys

from illustris_python.snapshot import loadHalo, snapPath, loadSubhalo
from illustris_python.groupcat import gcPath, loadHalos, loadSubhalos

def main():
    
    basePath = '/Volumes/G-RAID/simulations/unprocessed/Illustris/Illustris-1'
    snapNum = 135
    m_dm = 6.3*10**6.0 #dark matter particle mass
    little_h = 0.704
    Lbox = np.array([75.0]*3)

    ptcl_type = {'gas':0,
                 'dark-matter':1,
                 'tracers':3,
                 'stars/wind':4,
                 'black holes':5}

    from astropy.table import Table
    t_1 = Table.read('./data/galaxy_shapes_1.dat', format='ascii')
    t_2 = Table.read('./data/halo_shapes_1.dat', format='ascii')

    gal_ids = t_1['gal_id']

    central_ids = loadHalos(basePath, snapNum, fields=['GroupFirstSub'])  
    centrals_mask = np.in1d(gal_ids, central_ids)

    host_ids = loadSubhalos(basePath, snapNum, fields=['SubhaloGrNr'])
    host_ids = host_ids[gal_ids]

    galaxy_positions = loadSubhalos(basePath, snapNum, fields=['SubhaloPos'])/1000.0
    galaxy_positions = galaxy_positions[gal_ids]

    t = Table([gal_ids],names=['gal_ids'])
    t['central'] = 0
    t['central'][centrals_mask] = 1

    print(t)
    


if __name__ == '__main__':
    main()