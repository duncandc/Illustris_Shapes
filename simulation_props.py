"""
dictionaries of Illustris simulation properties
"""

import numpy as np

d_1 = {'basePath':'/Volumes/G-RAID/simulations/unprocessed/Illustris/Illustris-1',
       'm_dm': 6.3*10**6.0,
       'litte_h': 0.704,
       'Lbox': np.array([75.0]*3)
       }

d_1_dmo = {'basePath': '/Volumes/G-RAID/simulations/unprocessed/Illustris/Illustris-1-Dark',
           'm_dm': 6.3*10**6.0,
           'litte_h': 0.704,
           'Lbox': np.array([75.0]*3)
           }

sim_prop_dict = {'Illustris-1': d_1,
                 'Illustris-1-Dark': d_1_dmo}
