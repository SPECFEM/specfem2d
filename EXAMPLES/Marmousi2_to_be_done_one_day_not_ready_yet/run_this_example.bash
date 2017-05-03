#!/bin/bash -eu
#################################
# USER PARAMETER

NPROC=20

#################################

echo "setting up example: `pwd`"
echo

DIR_SEM="../.."

cp -v $DIR_SEM/bin/xspecfem2D $DIR_SEM/bin/xmeshfem2D ./

mkdir -p OUTPUT_FILES
# cleanup
rm -rf OUTPUT_FILES/*

echo
echo "running mesher..."
echo
./bin/xmeshfem2D

echo
echo "running solver..."
echo
mpirun -np $NPROC ./bin/xspecfem2D

echo
echo "done"
echo

