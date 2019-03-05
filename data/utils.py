"""
helper function modified from http://www.illustris-project.org/data/docs/api/
"""

import requests
from credentials import headers

def get(path, params=None, savepath='./'):
    """
    path : string

    params : dict

    savepath : string
    """

    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        filename = savepath+filename
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r
