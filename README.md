# SPECFEM2D

SPECFEM2D allows users to perform 2D and 2.5D (i.e., axisymmetric) simulations
of acoustic, elastic, viscoelastic, and poroelastic seismic wave propagation.
The package can also be used for full waveform imaging (FWI) or adjoint tomography.


Main "historical" developers: Dimitri Komatitsch and Jeroen Tromp
  (there are currently many more!)


## Installation

Instructions on how to install and use `SPECFEM2D` are
available in the

- PDF manual located in directory: [doc/USER_MANUAL](doc/USER_MANUAL).

- HTML manual (latest version): [specfem2d.readthedocs.io](http://specfem2d.readthedocs.io/)


For a quick test, run the default example with these commands:

  ./configure FC=gfortran
  make all
  ./bin/xmeshfem2D
  ./bin/xspecfem2D

and check the output files in ./OUTPUT_FILES/


## Development

[![Build Status](https://travis-ci.org/geodynamics/specfem2d.svg?branch=devel)](https://travis-ci.org/geodynamics/specfem2d) [![codecov](https://codecov.io/gh/geodynamics/specfem2d/branch/devel/graph/badge.svg)](https://codecov.io/gh/geodynamics/specfem2d)

Development is hosted on GitHub in the
[geodynamics/specfem2d repository](https://github.com/geodynamics/specfem2d).


To contribute, please follow the guidelines in the SPECFEM3D github wiki:
[specfem3d wiki](https://github.com/geodynamics/specfem3d/wiki)


## Computational Infrastructure for Geodynamics (CIG)

Seismology software repository: [SPECFEM2D](https://geodynamics.org/cig/software/specfem2d/)

