# Value Added Catalogs

This directory value added galaxy catalogs (VAGC) for various Illusrtris Simulations.  See the `README.md` in the main project directory for a full description of the galaxies that are in each catalog and the method used to calculate derived quantities quantities.

Many quantities are taken directly from the FoF halo catalogs and Subfind group cagtalogs provided by Illustris, described [here](http://www.illustris-project.org/data/docs/specifications/#sec2a). 


## Description of Data

### VAGC

The file names are of the form: [simulation name]\_[snapshot number]\_vagc.dat

Note that all values in the catalog are scaled to h=1.
  
Each catalog is saved in ascii format as space seperated values.  The first row contains the column names:

* `gal_id` : subfind galaxy/subhalo ID 
* `x`,`y`,`z` : position [h^-1 Mpc]
* `vx`, `vy`, `vz` : peculiar velocity [km/s]
* `stellar_mass_in_twice_halfrad` : stellar mass within twice the half-mass radius (exclusive of substructure) [h^-1 Msol]
* `stellar_mass_all` : stellar mass associated with subhalo (exclusive of substructure) [h^-1 Msol]
* `host_halo_mass_200m` : 200m halo mass [h^-1 Msol]
* `host_halo_radius_200m` : 200m halo radius [h^-1 Mpc] 
* `host_halo_fof_mass` : total fof halo mass [h^-1 Msol]


### Halo Matching

crossmatched catalogs of FoF host haloes between full physics and DMO runs

The file names are of the form: [simulation name]\_[snapshot number]\_halo_matches.dat

Each catalog is saved in ascii format as space seperated values.  The first row contains the column names:
 
* `host_halo_id` : full physics FoF halo ID
* `dmo_host_halo_id` : DMO FoF halo ID 
* `r` : distance between halo positions [h^-1 Mpc]


### Halo Catalogs

FoF host halo catalogs.

The file names are of the form: [simulation name]\_[snapshot number]\_halo_catalog.dat

Each catalog is saved in ascii format as space seperated values.  The first row contains the column names:
 
* `host_halo_id` : fof halo ID
* `central_id` : subfind galaxy/halo ID
* `x`, `y`, `z` : position [h^-1 Mpc]
* `vx`, `vy`, `vz` : peculiar velocity [km/s]
* `host_halo_mass_200m` : 200m halo mass [h^-1 Msol]
* `host_halo_radius_200m` : 200m halo radius [h^-1 Mpc] 
* `host_halo_fof_mass` : total fof halo mass [h^-1 Msol]


## Loading Illustris VAGC

Catalogs can be loaded using the astropy.ascii module.
For example, to open the z=0 Illustris-1 vagc catalog:

```
from astropy.table import Table
t = Table.read('Illustris-1_135_vagc.dat', format='ascii')
```
