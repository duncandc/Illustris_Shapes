"""
dictionaries of Illustris simulation properties
"""

from project_settings import base_savepath

import numpy as np

d_1 = {'basePath': base_savepath + 'Illustris-1',
       'm_dm': 6.3*10**6.0,
       'litte_h': 0.704,
       'Lbox': np.array([75.0]*3)
       }

d_1_dmo = {'basePath': base_savepath + 'Illustris-1-Dark',
           'm_dm': 6.3*10**6.0,
           'litte_h': 0.704,
           'Lbox': np.array([75.0]*3)
           }

d_2 = {'basePath': base_savepath + 'TNG300-1/output',
       'm_dm': 0.00398342749867548*10**10,
       'litte_h': 0.6774,
       'Lbox': np.array([205.00]*3)
       }

d_2_dmo = {'basePath': base_savepath + 'TNG300-1-Dark/output',
       'm_dm': 0.0047271638660809*10**10,
       'litte_h': 0.6774,
       'Lbox': np.array([205.00]*3)
       }

d_3 = {'basePath': base_savepath + 'TNG100-1/output',
       'm_dm': 0.000505574296436975*10**10,
       'litte_h': 0.6774,
       'Lbox': np.array([75.0]*3)
       }

d_3_dmo = {'basePath': base_savepath + 'TNG100-1-Dark/output',
       'm_dm': 0.000599968882709879*10**10,
       'litte_h': 0.6774,
       'Lbox': np.array([75.0]*3)
       }

sim_prop_dict = {'Illustris-1': d_1,
                 'Illustris-1-Dark': d_1_dmo,
                 'TNG300-1': d_2,
                 'TNG300-1-Dark': d_2_dmo,
                 'TNG100-1': d_3,
                 'TNG100-1-Dark': d_3_dmo}
