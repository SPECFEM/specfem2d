#!/bin/bash
#
# script runs mesher and solver
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directoy
echo
echo "setting up example..."
echo

mkdir -p OUTPUT_FILES
mkdir -p DATA

# Sets up local DATA/ directory
cd DATA/
cp ../Par_file Par_file
cp ../interfaces.dat interfaces.dat
cp ../SOURCE SOURCE
cd ../

# Clean output files
rm -rf OUTPUT_FILES/*

# Store setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running mesher..."
  echo
  xmeshfem2D > OUTPUT_FILES/output_mesher.txt
else
  # This is a MPI simulation
  echo
  echo "running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC xmeshfem2D > OUTPUT_FILES/output_mesher.txt
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver..."
  echo
  xspecfem2D > OUTPUT_FILES/output_solver.txt
else
  # This is a MPI simulation
  echo
  echo "running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC xspecfem2D > OUTPUT_FILES/output_solver.txt
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`

