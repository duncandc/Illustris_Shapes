"""
"""

from __future__ import print_function
import sys
from utils import get
import time

base_savepath = "/Volumes/G-RAID/simulations/unprocessed/Illustris/"

def main():

    if len(sys.argv)>1:
        sim = sys.argv[1]
        snapnum = str(sys.argv[2])
    else:
        sim = 'Illustris-1'
        snapnum = str(135)

    # restart download from file number xxx
    if len(sys.argv)>3:
        start_file_num = int(sys.argv[3])
    else:
        start_file_num = 0

    # check to see if it is a dm only simulation
    if sim[-4:]=='Dark':
        dm_only=True
    else:
        dm_only=False

    if sim[:3] == 'TNG':
        base_url = "http://www.tng-project.org/api/" + sim + "/ 
    else:
        base_url = "http://www.illustris-project.org/api/" + sim + "/"
    snapnum = snapnum.zfill(3)

    # location to save data
    if sim[:3] == 'TNG':
        savepath = base_savepath + sim + "/output/" + "snapdir_" + snapnum + "/"
    else:
        savepath = base_savepath + sim + "/" + "snapdir_" + snapnum + "/"
    print("saving files to: "+savepath)
    print("     ")

    sim_metadata = get(base_url)

    # select which columns to download
    if dm_only:
        params = {'dm':'Coordinates,Velocities'}
    else:
        params = {'stars':'Coordinates,Velocities,Masses,GFM_StellarFormationTime',
                  'dm':'Coordinates,Velocities'}

    outer_start = time.time()
    n_snaps = sim_metadata['num_files_snapshot']
    for i in range(start_file_num, n_snaps):
        inner_start = time.time()

        file_url = base_url + "files/snapshot-"+str(snapnum).zfill(3)+"." + str(i) + ".hdf5"
        saved_filename = get(file_url, params, savepath)

        dt = time.time()-inner_start
        n_remaining = n_snaps-i-1  # number of files left to download
        print(i, "estimated time remaining: {0} hours".format(n_remaining*dt/60.0/60.0))

    print("     ")
    print("total time ellapsed: {0} hours".format(time.time()/60.0/60.0-outer_start/60.0/60.0))


if __name__ == '__main__':
    main()
