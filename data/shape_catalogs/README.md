# Shape Catalogs

This directory stores shape catalogs for various Illusrtris Simulations.

## Description of Data
  
Each shape catalog is saved in ascii as space seperated values.  The first row contains the column names:

* gal_id : galaxy ID 
* a : primary eigenvalue
* b : intermediate eigenvalue
* c : minor eigenvalue
* av_x : x-component of the primary eigenvector
* av_y : y-component of the primary eigenvector
* av_z : z-component of the primary eigenvector
* bv_x : x-component of the intermediate eigenvector
* bv_y : y-component of the intermediate eigenvector
* bv_z : z-component of the intermediate eigenvector
* cv_x : x-component of the minor eigenvector
* cv_y : y-component of the minor eigenvector
* cv_z : z-component of the minor eigenvector

## Loading Shape Catalog Tables

Catalogs can be loaded using the astropy.ascii module.
For example, to open the z=0 Illustris-1 galaxy shapes catalog:

```
from astropy.table import Table
t = Table.read('Illustris-1_135_reduced_galaxy_shapes.dat', format='ascii')
```


