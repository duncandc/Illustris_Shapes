"""
create DMO only halo catalogs for Illustris Simulations
"""

from __future__ import print_function, division
import numpy as np
from astropy.table import Table, Column
from astropy.io import ascii
import sys

from simulation_props import sim_prop_dict

from illustris_python.groupcat import gcPath, loadHalos, loadSubhalos, loadHeader


def main():

    if len(sys.argv)>1:
        sim_name = sys.argv[1]
        snapNum = int(sys.argv[2])
    else:
        sim_name = 'Illustris-1' # full physics high-res run
        snapNum = 135  # z=0

    # load simulation properties
    d = sim_prop_dict[sim_name]
    basePath = d['basePath']
    m_dm = d['m_dm']
    litte_h = d['litte_h']
    Lbox = d['Lbox']

    # scale to get peculiar velocity
    header = loadHeader(basePath, snapNum)
    scale_factor = header['Time']

    # load host halo catalog
    fields = ['GroupPos', 'GroupVel', 'GroupMass', 'Group_M_Mean200', 'Group_R_Mean200', 'GroupFirstSub']
    halo_table = loadHalos(basePath, snapNum, fields=fields)
    host_ids = np.arange(0,len(halo_table['GroupMass'])).astype('int')

    coords = halo_table['GroupPos']/1000.0
    x = coords[:,0]
    y = coords[:,1]
    z = coords[:,2]

    vels = halo_table['GroupVel']*np.sqrt(scale_factor)
    vx = vels[:,0]
    vy = vels[:,1]
    vz = vels[:,2]

    halo_mass_200m = halo_table['Group_M_Mean200']*10**10
    host_halo_radius_200m = halo_table['Group_R_Mean200']/1000.0
    host_halo_fof_mass = halo_table['GroupMass']*10**10

    central_id = halo_table['GroupFirstSub']

    # save table
    fpath = './data/value_added_catalogs/'
    fname = sim_name + '_' + str(snapNum) + '_host_halo_catalog.dat'
    ascii.write([host_ids, central_id,
    	         x, y, z,
    	         vx, vy, vz,
    	         halo_mass_200m, host_halo_radius_200m,
    	         host_halo_fof_mass],
                fpath+fname,
                names=['host_halo_id', 'central_id',
                       'x', 'y', 'z',
                       'vx', 'vy', 'vz',
                       'host_halo_mass_200m', 'host_halo_radius_200m',
                       'host_halo_fof_mass'],
                overwrite=True)


if __name__ == '__main__':
    main()
