#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`


# sets up directory structure in current example directoy
echo
echo "setting up example..."
echo

mkdir -p OUTPUT_FILES
# cleans output files
rm -rf OUTPUT_FILES/*

# links executables
mkdir -p bin
cd bin/
rm -f xmeshfem2D xspecfem2D
ln -s ../../../bin/xmeshfem2D
ln -s ../../../bin/xspecfem2D
cd ../

# scripts
if [ ! -e createTomographyFile.py ]; then
  ln -s ../../utils/createTomographyFile.py
fi
if [ ! -e plot_model_from_gll_or_binary_file.py ]; then
  ln -s ../../utils/Visualization/plot_model_from_gll_or_binary_file.py
fi

# creates tomography file
echo
echo "setting up tomography model..."
echo
./createTomographyFile.py

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# links tomo file
cd DATA/
rm -f tomo_file.xyz
ln -s ../profile.xyz tomo_file.xyz
cd ../

# stores setup
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
  ./bin/xmeshfem2D
else
  # This is a MPI simulation
  echo
  echo "running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xmeshfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver..."
  echo
  ./bin/xspecfem2D
else
  # This is a MPI simulation
  echo
  echo "running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo

## plots model figures
echo
echo "plotting model figures..."
echo
./plot_model_from_gll_or_binary_file.py DATA/proc000000_rho.bin
./plot_model_from_gll_or_binary_file.py DATA/proc000000_vs.bin
./plot_model_from_gll_or_binary_file.py DATA/proc000000_vp.bin

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo `date`

