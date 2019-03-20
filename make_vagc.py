"""
script to measure galaxy shapes in the Illustris simmulations
"""

from __future__ import print_function, division
import numpy as np
from astropy.table import Table, Column
from astropy.io import ascii
import sys

from simulation_props import sim_prop_dict

from illustris_python.groupcat import gcPath, loadHalos, loadSubhalos


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


def main():

    if len(sys.argv)>1:
        sim_name = sys.argv[1]
        snapNum = int(sys.argv[2])
    else:
        sim_name = 'Illustris-1' # full physics high-res run
        snapNum = 135  # z=0

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

    fields = ['SubhaloGrNr', 'SubhaloMassInRadType', 'SubhaloMassType', 'SubhaloPos', 'SubhaloVel']

    galaxy_table = loadSubhalos(basePath, snapNum, fields=fields)
    
    # FoF host halo ID
    host_halo_ids = galaxy_table['SubhaloGrNr'][gal_ids]
    
    x = galaxy_table['SubhaloPos'][:,0][gal_ids]
    y = galaxy_table['SubhaloPos'][:,1][gal_ids]
    z = galaxy_table['SubhaloPos'][:,2][gal_ids]

    vx = galaxy_table['SubhaloVel'][:,0][gal_ids]
    vy = galaxy_table['SubhaloVel'][:,1][gal_ids]
    vz = galaxy_table['SubhaloVel'][:,2][gal_ids]

    mstar_in_twice_halfrad = galaxy_table['SubhaloMassInRadType'][:,4][gal_ids]
    mstar_in_twice_halfrad = mstar_in_twice_halfrad*10**10

    mstar_all = galaxy_table['SubhaloMassType'][:,4][gal_ids]
    mstar_all = mstar_all*10**10


    # central galaxy ID of group
    fields = ['GroupFirstSub', 'Group_M_Mean200', 'Group_R_Mean200', 'GroupMass']
    host_halo_table = loadHalos(basePath, snapNum, fields=fields)
    central_id = host_halo_table['GroupFirstSub'][host_halo_ids]

    # host halo properties
    host_halo_mass_200m = host_halo_table['Group_M_Mean200'][host_halo_ids]
    host_halo_radius_200m = host_halo_table['Group_R_Mean200'][host_halo_ids]
    host_halo_fof_mass = host_halo_table['GroupMass'][host_halo_ids]

    # save table
    fpath = './data/value_added_catalogs/'
    fname = sim_name + '_' + str(snapNum) + '_vagc.dat'
    ascii.write([gal_ids, host_halo_ids, central_id,
    	         x, y, z,
    	         vx, vy, vz,
    	         mstar_in_twice_halfrad, mstar_all,
    	         host_halo_mass_200m, host_halo_radius_200m, host_halo_fof_mass],
                fpath+fname,
                names=['gal_id', 'host_halo_id', 'central_id',
                 'x', 'y', 'z',
                 'vx', 'vy', 'vz',
                 'stellar_mass_in_twice_halfrad', 'stellar_mass_all',
                 'host_halo_mass_200m', 'host_halo_radius_200m', 'host_halo_fof_mass'],
                overwrite=True)



if __name__ == '__main__':
    main()