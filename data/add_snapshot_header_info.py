"""
adds information to header of snapshot files
"""

import h5py
import sys
import os
import numpy as np

base_savepath = "/Volumes/G-RAID/simulations/unprocessed/Illustris/"

def main():

    if len(sys.argv)>1:
        sim_name = sys.argv[1]
        snapnum = sys.argv[2]
    else:
        sim_name = 'TNG100-1'
        snapnum = str(99)

    snapnum = snapnum.zfill(3)  

    # location of data
    if sim_name[:3] == 'TNG':
        savepath = base_savepath + sim_name + "/output/" + "snapdir_" + snapnum + "/"
    else:
        savepath = base_savepath + sim_name + "/" + "snapdir_" + snapnum + "/"
    
    # get list of snapshot files
    filenames = os.listdir(savepath)

    part_types = [0,1,2,3,4,5]

    for filename in filenames:
    	if filename[-4:]=='hdf5':
    		f = h5py.File(savepath + filename, 'r')

    		num_ptcls_type = np.zeros(5, dtype='int32')
    		for part_type in part_types:
    			try:
    			    x = f['PartType'+str(part_type)]
    			    keys = x.keys()
    			except KeyError:
    				keys = []
    			
    			if len(keys)>0:
    				num_ptcls_type[part_type] = len(x[keys[0]])

    		print(filename, num_ptcls_type)

    		# add/replce header entry for 'NumPart_ThisFile' 
    		#f['Header']['NumPart_ThisFile'] = num_ptcls_type






if __name__ == '__main__':
    main()