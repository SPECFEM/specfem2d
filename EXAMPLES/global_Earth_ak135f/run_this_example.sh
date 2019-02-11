#!/bin/bash

mkdir -p OUTPUT_FILES
mkdir -p bin
rm -rf OUTPUT_FILES/*

# links executables
mkdir -p bin
cd bin/
rm -f xmeshfem2D xspecfem2D
ln -s ../../../bin/xmeshfem2D
ln -s ../../../bin/xspecfem2D
cd ../

# Build the mesh
./make_specific_mesher_for_2D_Earth.csh
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo

# Runs the mesher (to create the database)
./bin/xmeshfem2D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# Get the number of processors from Par_file
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

if [ "$NPROC" -eq 1 ]; then # If NPROC == 1 this is a serial simulation
  echo
  echo " Running solver..."
  echo
  ./bin/xspecfem2D
else # Else this is a MPI simulation
  echo
  echo " Running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi