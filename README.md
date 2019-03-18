# Illustris Galaxy Shapes

This is a project to measure galaxy and halo shapes/orientations in the Illustris Simulations.

![](./notebooks/figures/demo_shapes.png)

## Requirements

This project requires the following Python packages installed:

* [numpy](http://www.numpy.org)
* [astropy](http://www.astropy.org)
* [inertia_tensors](https://github.com/duncandc/inertia_tensors/edit/master/README.md)
* [illustris_python](https://bitbucket.org/illustris/illustris_python)


## Galaxy & Halo Shapes

Currently, I calculate galaxy and halo shapes for the following simulations:

* Illustris-1
	* z=0.0 (snapshot 135)
	* z=0.6 (snapshot 099)
	* z=1.0 (snapshot 085)
* Illustris-1-Dark
	* z=0.0 (snapshot 135)
	* z=0.6 (snapshot 099)
	* z=1.0 (snapshot 085)

Galaxy shapes are calculated for galaxies with a stellar mass within two times the stellar half-mass radius of at least log(Mstar)>= 9 + log(h), using all stellar particles (excluding wind particles) that are within two times the stellar half-mass radius.  The center of a galaxy is taken to be the most bound particle within the subfind subhalo, regardless of particle type.

(Sub-)halo shapes are calculated using all dark matter particles that belong to a subfind 'subhalo' with at least 1000 particles.  The center of (sub-)haloes is taken to be the most bound particle in the subfind subhalo, regardless of particle type.  Note that subfind sub-haloes exclude particles that belong to substructures, i.e. central subhaloes do not include particles from satellite subhaloes.  

Galaxy and halo shapes/orientations are calculated using the following scripts:

* `calculate_galaxy_shapes.py`
* `calculate_halo_shapes.py`

Each of these scripts can be run to create a shape catalog for a given simulation, snapshot, and shape calculation method.  For example, to calculate galaxy shapes for the Illustris-1 simulation at z=0 (snapshot 135) using the reduced inertia tensor method, you would execute the following command:

```
$user: python calculate_galaxy_shapes.py Illustris-1 135 reduced
```

The resulting shape catalogs are saved in the `./data/shape_catalogs/` directory.

All galaxy and halo shapes are determined by calculating an inertia tensor for a particle distribution.  The code to calculate inertia tensors is part of the [inertia_tensors](https://github.com/duncandc/inertia_tensors/edit/master/README.md) package.  There are three methods implemented: non-reduced, reduced, and iterative.  See the docs in the [inertia_tensors](https://github.com/duncandc/inertia_tensors/edit/master/README.md) package for details on how these quantities are calculated.

Note that in order to run these scripts, you must have the required particle data downloaded on to your disk.  The location of this data for each simulation is set in the `simulation_props.py` file as the `basePath` key in each dictionary.  Scripts to download particle data are stored in the `./data/` directory.


## Galaxy Circularity

For each sample for which galaxy shapes are calculated, I also calculate the specific angular momentum and circularity for galaxies.  This is calcuated using the `calculate_galaxy_circularity.py` script using the same positional arguments as the shape scripts described above.  

The specific angular momentum for galaxies is calculated using all star particles (excluding wind particles) within two times the stellar half-mass radius.  Again, the center of a galaxy is taken to be the most bound particle within the subfind subhalo, regardless of particle type.  I record the magnitude and direction of each galaxies specific angular momentum.  

In addition, the circularity is calculated for each stellar particle.  See [Scannapieco et al. (2009)](https://arxiv.org/abs/0812.0976) for details of the calculation.  Briefly, I calculate the component of a galaxy's a stellar particle's angular momentum aligned with the global angular momentum axis.  This is compared to that of a particle on a circular orbit.  Particles with high circularity (e>0.7) are considered disk stars.  I record the fraction of stellar mass in a disk component.  

The resulting circularity catalogs are saved in the `./data/shape_catalogs/` directory.  


## Data

This project uses a large amount of particle data available on the Illustris data access webpage.  Scripts to efficiently download the required data are available in the `./data` directory.

The data products created by this project are stored in the `./data` directory. 

contact:
duncanc@andrew.cmu.edu
