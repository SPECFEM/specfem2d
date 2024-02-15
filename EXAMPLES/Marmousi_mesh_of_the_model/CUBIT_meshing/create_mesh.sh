#!/bin/bash

mkdir -p MESH

# create all surfaces from horizons
./1_create_surfaces.py
if [[ $? -ne 0 ]]; then exit 1; fi

# mesh surfaces with quads
./2_mesh_marmousi.py
if [[ $? -ne 0 ]]; then exit 1; fi

# improve mesh quality
./3_smooth_mesh.py
if [[ $? -ne 0 ]]; then exit 1; fi

# export to SPECFEM2D external mesh format
./4_export_mesh.py
if [[ $? -ne 0 ]]; then exit 1; fi

# create material properties file nummaterial_velocity_file_marmousi2
./5_convert_surface_rock_to_velocities.py
#or
#./5_convert_surface_rock_to_velocities.py --without-water
if [[ $? -ne 0 ]]; then exit 1; fi


echo
echo "all done"
echo

