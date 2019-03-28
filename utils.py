"""
utility functions
"""

from astropy.table import Table
from project_settings import PROJECT_DIRECTORY

import numpy as np
from halotools.utils import crossmatch


def match_centrals(simname='Illustris-1', snapnum=135, data_path=PROJECT_DIRECTORY+'data/value_added_catalogs/'):
    """
    match central galaxies and primary subhaloes between full physics and dmo simulations.

    Returns
    -------
    fp_central_id, dmo_central_id : np.arrays
        arrays of matching subfind halo IDs

    Notes
    -----
    This function requires that host halo matches and halo catalogs have been precomputed.
    The location of the precomputed files is sepcified by ``data_path``.
    """

    # open host halo matching file
    fname = simname + '_' + str(snapnum) + '_host_halo_matches.dat'
    host_halo_match_table = Table.read(data_path + fname, format='ascii')

    # open full physics host halo catalog
    fname = simname + '_' + str(snapnum) + '_host_halo_catalog.dat'
    fp_halo_table = Table.read(data_path + fname, format='ascii')

    # open dmo host halo catalog
    fname = simname + '-Dark' + '_' + str(snapnum) + '_host_halo_catalog.dat'
    dmo_halo_table = Table.read(data_path + fname, format='ascii')

    # match halo catalogs to matches
    idx_1, idy_1 = crossmatch(fp_halo_table['host_halo_id'], host_halo_match_table['host_halo_id'])

    fp_central_id = np.array(fp_halo_table['central_id'][idx_1])
    dmo_host_id = host_halo_match_table['dmo_host_halo_id'][idy_1]


    idx_2, idy_2 = crossmatch(dmo_host_id, dmo_halo_table['host_halo_id'])

    dmo_central_id = np.zeros(len(dmo_host_id)).astype('int')-1

    dmo_central_id[idx_2] = dmo_halo_table['central_id'][idy_2]


    return fp_central_id, dmo_central_id


