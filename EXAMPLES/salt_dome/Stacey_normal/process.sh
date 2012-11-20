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
#cp ../interfaces_elastic_analytic.dat .
cp ../SOURCE SOURCE
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

# compiles executables in root directory
cd ../../../
make > tmp.log
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
./xmeshfem2D > OUTPUT_FILES/output_mesher.txt

# runs simulation
echo
echo "  running solver..."
echo
./xspecfem2D > OUTPUT_FILES/output_solver.txt

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`
