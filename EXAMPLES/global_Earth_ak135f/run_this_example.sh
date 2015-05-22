#!/bin/bash

./make_specific_mesher_for_2D_Earth.csh # Build the mesh

rm -rf OUTPUT_FILES
mkdir -p OUTPUT_FILES # Will contains the results of the simulation

../../bin/xmeshfem2D # Runs the mesher (to create the database)

# Get the number of processors from Par_file
NPROC=`grep nproc DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

if [ "$NPROC" -eq 1 ]; then # If nproc == 1 this is a serial simulation
  echo
  echo " Running solver..."
  echo
  ../../bin/xspecfem2D
else # Else this is a MPI simulation
  echo
  echo " Running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ../../bin/xspecfem2D
fi
