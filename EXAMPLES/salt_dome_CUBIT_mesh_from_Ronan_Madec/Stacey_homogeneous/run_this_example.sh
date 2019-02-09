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
cp ../modelY1_absorbing_surface_file modelY1_absorbing_surface_file
cp ../modelY1_free_surface_file modelY1_free_surface_file
cp ../modelY1_materials_file modelY1_materials_file
cp ../modelY1_mesh_file modelY1_mesh_file
cp ../modelY1_nodes_coords_file modelY1_nodes_coords_file
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

cd $currentdir

# links executables
mkdir -p bin
cd bin/
rm -f xmeshfem2D xspecfem2D
ln -s ../../../bin/xmeshfem2D
ln -s ../../../bin/xspecfem2D
cd ../

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# runs database generation
echo
echo "  running mesher..."
echo
./bin/xmeshfem2D

# runs simulation
echo
echo "  running solver..."
echo
./bin/xspecfem2D

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
