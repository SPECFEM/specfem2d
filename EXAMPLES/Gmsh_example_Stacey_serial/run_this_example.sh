#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take just a few seconds)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p OUTPUT_FILES
mkdir -p DATA

# sets up local DATA/ directory
cd DATA/
cp ../Par_file_Gmsh_SqrCircles.in Par_file
cp ../SOURCE.SqrCirc SOURCE
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

cd $currentdir

# links executables
rm -f xmeshfem2D xspecfem2D
ln -s ../../bin/xmeshfem2D
ln -s ../../bin/xspecfem2D

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
cp DATA/SOURCE_xz.dat OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/STATIONS_target OUTPUT_FILES/

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
