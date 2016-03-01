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

cd $currentdir

# runs database generation
# In principle we do not need rerun xmeshfem2D in the adjoint run.
# However, in the current structure of the code, the xspecfem2D program can not detect the
# the parameter change in Par_file. Since in the adjoint run we change the SIMULATION_TYPE and the save_forward
echo
echo "  running mesher..."
echo
./xmeshfem2D

# runs simulation
echo
echo "  running adjoint and backward solver..."
echo
./xspecfem2D

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

# plot the traveltime kernel
cp ../../utils/visualization/plot_wavefield.pl ./
./plot_wavefield.pl 400/3000/400 400/2800/800 0/480/0/480 120/20/120/20 -8/-8/-8 1/1/1 -48.0/0.06 0/1/0/0/0 2.0/1/0 1/0/1/120 Tape2007_kernel onerec_homo


echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
