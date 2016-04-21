#!/bin/bash
#
# This script runs Specfem2d mesher and solver
# with AXISYM=.true.
# This is then an axisymmetric simulation
#
# To clean this directory : rm -rf DATA/ OUTPUT_FILES/ x* *~

echo "Simple axisymmetric simulation : `date`"
currentdir=`pwd`

# Sets up directory structure in current example directoy
echo
echo " Setting up example..."
echo

mkdir -p OUTPUT_FILES
mkdir -p DATA

# Sets up local DATA/ directory
cd DATA/
ln -s ../Par_file_axisym Par_file
ln -s ../SOURCE_axisym SOURCE
cd ../

# Clean output files
rm -rf OUTPUT_FILES/*

#cd ../../
#make -j 16
cd $currentdir

# Link executables
rm -f xmeshfem2D xspecfem2D
ln -s ../../bin/xmeshfem2D
ln -s ../../bin/xspecfem2D

# Store setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# Get the number of processors
NPROC=`grep -i NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

echo $NPROC

# Run database generation
echo
echo " Running mesher..."
echo
./xmeshfem2D

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

# Store output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "See results in directory: OUTPUT_FILES/"
echo
echo "Done!"
date
