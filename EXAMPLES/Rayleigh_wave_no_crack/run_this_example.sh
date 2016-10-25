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
cp ../Par_file_Rayleigh_2D Par_file
cp ../interfaces_Rayleigh_flat.dat .
cp ../SOURCE_Rayleigh_2D SOURCE
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

cd $currentdir

# links executables
rm -f xmeshfem2D xspecfem2D
ln -s ../../bin/xmeshfem2D
ln -s ../../bin/xspecfem2D

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# Get the number of processors from Par_file
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
echo
echo "  running mesher..."
echo
./xmeshfem2D

# runs simulation
if [ "$NPROC" -eq 1 ]; then # This is a serial simulation
  echo
  echo " Running solver..."
  echo
  ./xspecfem2D
else # This is a MPI simulation
  echo
  echo " Running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./xspecfem2D
fi

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
