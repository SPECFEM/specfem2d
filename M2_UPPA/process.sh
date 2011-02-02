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
ln -s ../Par_file_M2_UPPA Par_file
ln -s ../SOURCE_M2_UPPA SOURCE
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

# compiles executables in root directory
cd ../../
make > tmp.log
cd $currentdir

# links executables
ln -s ../../xmeshfem2D
ln -s ../../xspecfem2D

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# runs database generation
echo
echo "  running mesher..."
echo
./xmeshfem2D > OUTPUT_FILES/output_mesher.txt

# runs simulation
echo
echo "  running solver..."
echo
./xspecfem2D > OUTPUT_FILES/output_solver.txt

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`
