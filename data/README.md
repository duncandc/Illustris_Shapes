# Data

This directory stores data associated with this project and scripts to retreive external data files.


## Downloading Data

In order to download particle data, you must have an account on the Illustris project [website](http://www.illustris-project.org/data/).  You will need to modify the API key in `credentials.py` before running any of the download scripts. 


### Particle Data

Warning, the particle data from Illustris-1 for this project is ~200 Gb per snapshot.

The location to save particle data is set in `download_ptcl_data.py` by modifying the `base_savepath` variable.

Individual data sets can be downloaded.  For example, to download Illustris-1 z=0 (snapshot 135) execute the following command: 

```
$user: python download_ptcl_data.py Illustris-1 135
```

All particle data sets used in this project can be downloaded by running the `download_all_ptcl_data.sh` script. 

The Illustris TNG snapshots downloaded using the `download_ptcl_data.py` script are missing some header information.  In this case, you need to run the `add_snapshot_header_info.py` script for each snapshot.  For example, this script can be run as:

```
$user: python add_snapshot_header_info.py TNG100-1 99
``` 


### Halo Catalogs

Halo catalogs are downloaded in their entirety.  


### Genel et al. (2015) Shape Catalog

There is a publicaly available shape catalog for Illustris-1 from [Genel et al. (2015)](http://adsabs.harvard.edu/abs/2015arXiv150301117G) described [here](http://www.illustris-project.org/data/docs/specifications/#sec4c).  This only includes the eigenvalues and not the eigenvectors, i.e. no orientations.  I have also confirmed that both the `MassTensorEigenVals` and `ReducedMassTensorEigenVals` use the same eigenvectors in the calculation, so there is a slight descrepancy between my calculation of the reduced mass tensor eigenvalues and those in the Genel catalog (private communication).

## Directory Structure

The Illustris data should be stored within a specific directory structure in order for the convenience functions in [illustris_python](https://bitbucket.org/illustris/illustris_python) to function.  The structures are slightly different for OG and TNG Illustris.

For TNG simulations:

- [base_savepath]
	- [simname]
		- output
			- snapdir_[snapnum]
			- groups_[snapnum]
		- postprocessing
			- offsets

For OG simulations:

- [base_savepath]
	- [simname]
		- snapdir_[snapnum]
		- groups_[snapnum]

Note that the offsets are downloaded as a seperate file for Illustris TNG.  

## Notes

checksums can be run on the data on Macs using the default shasum program by executing the following command in the directory storing the files.  

```
shasum -a 256 -c checksums.txt
```


Git will ignore changes to `credentials.py`.  This was accomplished with the following command:

```
git update-index --assume-unchanged credentials.py
```