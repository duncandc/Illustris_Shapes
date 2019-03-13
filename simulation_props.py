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

d_2 = {'basePath': '/Volumes/G-RAID/simulations/unprocessed/Illustris/TNG300-1',
       'm_dm': 0.00398342749867548*10**6.0,
       'litte_h': 0.6774,
       'Lbox': np.array([205.00]*3)
       }

sim_prop_dict = {'Illustris-1': d_1,
                 'Illustris-1-Dark': d_1_dmo,
                 'TNG300-1': d_2}
