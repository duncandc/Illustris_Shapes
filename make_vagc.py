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

    import h5py
    f = h5py.File(basePath + '/snapdir_135/stellar_circs.hdf5', 'r')
    circ_data = f.get('Snapshot_135')
    f_disk = np.arrayt(circ_data['CircAbove07Frac'].value)
    f_bulge = 1.0-f_disk

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

    t['host_halo_id'] = host_ids

    t['x'] = galaxy_positions[:,0]
    t['y'] = galaxy_positions[:,1]
    t['z'] = galaxy_positions[:,2]

    t['galaxy_axisA_x'] = t_1['av_x']
    t['galaxy_axisA_y'] = t_1['av_y']
    t['galaxy_axisA_z'] = t_1['av_z']

    t['galaxy_axisB_x'] = t_1['bv_x']
    t['galaxy_axisB_y'] = t_1['bv_y']
    t['galaxy_axisB_z'] = t_1['bv_z']

    t['galaxy_axisC_x'] = t_1['cv_x']
    t['galaxy_axisC_y'] = t_1['cv_y']
    t['galaxy_axisC_z'] = t_1['cv_z']

    t['galaxy_b_to_a'] = t_1['b']/t_1['a']
    t['galaxy_c_to_a'] = t_1['c']/t_1['a']

    from halotools.utils import crossmatch
    idx, idy = crossmatch(t['gal_ids'], circ_data['SubfindID'].value)
    t['galaxy_f_bulge'] = -99.0
    t['galaxy_f_bulge'][idx] = f_bulge[idy]

    # some halo properties
    host_halo_sizes = loadHalos(basePath, snapNum, fields=['Group_R_Mean200'])/1000.0
    host_halo_sizes = host_halo_sizes[host_ids]

    host_halo_m200 = loadHalos(basePath, snapNum, fields=['Group_M_Mean200'])*10**10
    host_halo_m200 = host_halo_m200[host_ids]

    t['host_halo_m200b'] = host_halo_sizes
    t['host_halo_r200b'] = host_halo_m200

    # halo shapes
    halo_mask = np.in1d(t_2['halo_id'], host_ids[centrals_mask])
    
    t['halo_axisA_x'] = -99.0
    t['halo_axisA_x'][centrals_mask] = t_2['av_x'][halo_mask]
    t['halo_axisA_y'] = -99.0
    t['halo_axisA_y'][centrals_mask] = t_2['av_y'][halo_mask]
    t['halo_axisA_z'] = -99.0
    t['halo_axisA_z'][centrals_mask] = t_2['av_z'][halo_mask]

    t['halo_axisB_x'] = -99.0
    t['halo_axisB_x'][centrals_mask] = t_2['bv_x'][halo_mask]
    t['halo_axisB_y'] = -99.0
    t['halo_axisB_y'][centrals_mask] = t_2['bv_y'][halo_mask]
    t['halo_axisB_z'] = -99.0
    t['halo_axisB_z'][centrals_mask] = t_2['bv_z'][halo_mask]

    t['halo_axisC_x'] = -99.0
    t['halo_axisC_x'][centrals_mask] = t_2['cv_x'][halo_mask]
    t['halo_axisC_y'] = -99.0
    t['halo_axisC_y'][centrals_mask] = t_2['cv_y'][halo_mask]
    t['halo_axisC_z'] = -99.0
    t['halo_axisC_z'][centrals_mask] = t_2['cv_z'][halo_mask]


    t['halo_b_to_a'] = -99.0
    t['halo_b_to_a'][centrals_mask] = t_2['b'][halo_mask]/t_2['a'][halo_mask]
    
    t['halo_c_to_a'] = -99.0
    t['halo_c_to_a'][centrals_mask] = t_2['c'][halo_mask]/t_2['a'][halo_mask]

    print(t)

    # save catalog
    fpath = './data/'
    fname = 'illustris_shapes_vagc_1.dat'
    t.write(fpath + fname, format='ascii', overwrite=True)


if __name__ == '__main__':
    main()