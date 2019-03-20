# Value Added Catalogs

This directory value added galaxy catalogs (VAGC) for various Illusrtris Simulations.  See the `README.md` in the main project directory for a full description of the galaxies that are in each catalog and the method used to calculate derived quantities quantities.

Many quantities are taken directly from the FoF group catalogs and Subfind group cagtalogs provided by Illustris, described [here](http://www.illustris-project.org/data/docs/specifications/#sec2a). 


## Description of Data

The file names are of the form: [simulation name]\_[snapshot number]\_vagc.dat

Note that all values in the catalog are scaled to h=1.
  
Each catalog is saved in ascii format as space seperated values.  The first row contains the column names:

* gal_id : subfind galaxy/subhalo ID 


## Loading Illustris VAGC

Catalogs can be loaded using the astropy.ascii module.
For example, to open the z=0 Illustris-1 galaxy shapes catalog:

```
from astropy.table import Table
t = Table.read('Illustris-1_135_vagc.dat', format='ascii')
```
