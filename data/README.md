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


### Halo Catalogs

Halo catalogs are downloaded in their entirety.  


## Notes

checksums can be run on the data on Macs using the default shasum program by executing the following command in the directoruy storing the files.  

```
shasum -a 256 -c checksums.txt
```


Git will ignore changes to `credentials.py`.  This was accomplished with the following command:

```
git update-index --assume-unchanged path/to/file.txt
```