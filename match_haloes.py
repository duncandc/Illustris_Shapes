"""
For Illuistris simulations, match DMO host haloes to full physics host haloes.
"""


from __future__ import print_function, division
import numpy as np
from astropy.table import Table, Column
from astropy.io import ascii
import sys

from simulation_props import sim_prop_dict

from illustris_python.groupcat import gcPath, loadHalos, loadSubhalos

from scipy.spatial import cKDTree as KDTree


def main():

    if len(sys.argv)>1:
        sim_name = sys.argv[1]
        snapNum = int(sys.argv[2])
    else:
        sim_name = 'Illustris-1' # full physics high-res run
        snapNum = 135  # z=0

    # load simulation properties
    d = sim_prop_dict[sim_name]
    basePath_1 = d['basePath']
    m_dm_1 = d['m_dm']

    d = sim_prop_dict[sim_name + '-Dark']
    basePath_2 = d['basePath']
    m_dm_2 = d['m_dm']

    litte_h = d['litte_h']
    Lbox = d['Lbox']

    # load full physics catalog
    fields = ['GroupPos', 'GroupMass', 'Group_R_Mean200']
    halo_table_1 = loadHalos(basePath_1, snapNum, fields=fields)
    host_ids_1 = np.arange(0,len(halo_table_1['GroupMass'])).astype('int')

    coords_1 = halo_table_1['GroupPos']/1000.0

    # load DMO catalog
    fields = ['GroupPos', 'GroupMass']
    halo_table_2 = loadHalos(basePath_2, snapNum, fields=fields)
    host_ids_2 = np.arange(0,len(halo_table_2['GroupMass'])).astype('int')

    # build KD tree
    coords_2 = halo_table_2['GroupPos']/1000.0
    tree = KDTree(coords_2, boxsize=Lbox)

    # query tree for nearest neighbors
    r_max = 1.0
    result = tree.query(coords_1, k=1, distance_upper_bound=r_max)

    idx = result[1]
    r = result[0]

    no_match = (idx==len(halo_table_2['GroupMass']))
    idx[no_match] = -1

    # save table
    fpath = './data/value_added_catalogs/'
    fname = sim_name + '_' + str(snapNum) + '_host_halo_matches.dat'
    ascii.write([host_ids_1, idx, r],
                fpath+fname,
                names=['host_halo_id', 'dmo_host_halo_id', 'r'],
                overwrite=True)






if __name__ == '__main__':
    main()
