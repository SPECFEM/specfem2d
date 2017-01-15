#!/bin/bash

rm -f *.o *.mod *genmod* xadd_CPML_layers_to_an_existing_mesh xconvert_external_layers_of_a_given_mesh_to_CPML_layers

gfortran -O3 -o xadd_CPML_layers_to_an_existing_mesh add_CPML_layers_to_an_existing_mesh.f90

gfortran -O3 -o xconvert_external_layers_of_a_given_mesh_to_CPML_layers convert_external_layers_of_a_given_mesh_to_CPML_layers.f90

rm -f *.o *.mod *genmod*

