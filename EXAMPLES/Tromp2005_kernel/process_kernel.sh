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

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
