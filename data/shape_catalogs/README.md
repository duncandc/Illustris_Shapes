# Shape Catalogs

This directory stores shape catalogs for various Illusrtris Simulations.  See the `README.md` in the main project directory for a full description of the galaxies/haloes that are in each catalog and the method used to calculate these quantities. 

## Description of Data

The file names are of the form:[simulation name]\_[snapshot number]\_[inertia tensor type]\_[galaxy/halo]_shapes.dat
  
Each shape catalog is saved in ascii as space seperated values.  The first row contains the column names:

* gal_id : galaxy ID 
* a : primary half-axis length
* b : intermediate half-axis length
* c : minor half-axis length
* av_x : x-component of the primary eigenvector
* av_y : y-component of the primary eigenvector
* av_z : z-component of the primary eigenvector
* bv_x : x-component of the intermediate eigenvector
* bv_y : y-component of the intermediate eigenvector
* bv_z : z-component of the intermediate eigenvector
* cv_x : x-component of the minor eigenvector
* cv_y : y-component of the minor eigenvector
* cv_z : z-component of the minor eigenvector

Note that the half-axis lengths are the square root of the eigenvalues.

## Loading Shape Catalog Tables

Catalogs can be loaded using the astropy.ascii module.
For example, to open the z=0 Illustris-1 galaxy shapes catalog:

```
from astropy.table import Table
t = Table.read('Illustris-1_135_reduced_galaxy_shapes.dat', format='ascii')
```


