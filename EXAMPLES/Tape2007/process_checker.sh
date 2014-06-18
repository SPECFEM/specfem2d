#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 1 minute)"
echo

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

\rm -rf OUTPUT_FILES DATA
mkdir OUTPUT_FILES
mkdir DATA

# sets up local DATA/ directory
cd DATA/
ln -s ../Par_file_Tape2007_132rec_checker Par_file
ln -s ../SOURCE_005 SOURCE
ln -s ../model_velocity.dat_checker model_velocity.dat_input 
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

# compiles executables in root directory
cd ../../
make > tmp.log
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
./xmeshfem2D > OUTPUT_FILES/output_mesher.txt

# runs simulation
echo
echo "  running solver..."
echo
./xspecfem2D > OUTPUT_FILES/output_solver.txt

# stores output
cp DATA/SOURCE_xz.dat OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/STATIONS_target OUTPUT_FILES/

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`
