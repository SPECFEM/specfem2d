#!/bin/bash
#
# script runs solver in serial
#

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 1 minute)"
echo

# runs simulation
echo
echo "  running mesher and solver..."
echo
./xmeshfem2D > OUTPUT_FILES/output_mesher_kernel.txt
./xspecfem2D > OUTPUT_FILES/output_solver_kernel.txt

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`
