#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take a few minutes)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p OUTPUT_FILES
mkdir -p DATA

# sets up local DATA/ directory
cd DATA/
cp ../Par_file Par_file
cp ../SOURCE SOURCE
cp ../STATIONS STATIONS
cp ../modelY1_absorbing_surface_file modelY1_absorbing_surface_file
cp ../modelY1_elements_CPML_list modelY1_elements_CPML_list
cp ../modelY1_free_surface_file modelY1_free_surface_file
cp ../modelY1_materials_file modelY1_materials_file
cp ../modelY1_mesh_file modelY1_mesh_file
cp ../modelY1_nodes_coords_file modelY1_nodes_coords_file
#cp ../modelY1_geometry_ACIS_ascii_format.sat modelY1_geometry_ACIS_ascii_format.sat
#cp ../modelY1_geometry_IGES_ascii_format.iges modelY1_geometry_IGES_ascii_format.iges
#cp ../modelY1_geometry_STL_ascii_format.txt modelY1_geometry_STL_ascii_format.txt
#cp ../modelY1_PML_inner_edge_for_plane_wave_file modelY1_PML_inner_edge_for_plane_wave_file
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

cd $currentdir

# links executables
rm -f xmeshfem2D xspecfem2D
ln -s ../../../bin/xmeshfem2D
ln -s ../../../bin/xspecfem2D

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# runs database generation
echo
echo "  running mesher..."
echo
./xmeshfem2D

# runs simulation
echo
echo "  running solver..."
echo
./xspecfem2D

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
