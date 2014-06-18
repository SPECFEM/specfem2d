#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 2 minutes)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p OUTPUT_FILES
mkdir -p DATA

# sets up local DATA/ directory
cd DATA/
ln -s ../Par_file_Abel_Balanche_bathy_source_solid Par_file
ln -s ../SOURCE_Abel_Balanche_bathy_source_solid SOURCE
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

# compiles executables in root directory
cd ../../
make
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
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`
