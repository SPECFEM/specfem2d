#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`

# compiles executables in root directory
cd ../../
make > tmp.log
cd $currentdir

# links executables
rm -f xmeshfem2D xspecfem2D
ln -s ../../bin/xmeshfem2D
ln -s ../../bin/xspecfem2D

# stores setup
cp DATA/Par_file OUTPUT_FILES/Par_file_kernel

# runs database generation
echo
echo "  running mesher..."
echo
./xmeshfem2D > OUTPUT_FILES/output_mesher.txt

# runs simulation
echo
echo "  running solver..."
echo
./xspecfem2D > OUTPUT_FILES/output_kernel.txt

# stores output
cp DATA/SOURCE_xz.dat OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
cp DATA/STATIONS_target OUTPUT_FILES/

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
