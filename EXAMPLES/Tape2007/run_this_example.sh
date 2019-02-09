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

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p OUTPUT_FILES

# cleans output files
rm -rf OUTPUT_FILES/*

# sets up local DATA/ directory
cd DATA/
rm -f Par_file
ln -s Par_file_Tape2007_132rec_checker Par_file
rm -f SOURCE
ln -s SOURCE_005 SOURCE
rm -f proc000000_model_velocity.dat_input
ln -s model_velocity.dat_checker proc000000_model_velocity.dat_input
rm -f STATIONS
ln -s STATIONS_checker STATIONS
cd ../

# compiles executables in root directory
#cd ../../
#make > tmp.log

cd $currentdir

# links executables
mkdir -p bin
cd bin/
rm -f xmeshfem2D xspecfem2D
ln -s ../../../bin/xmeshfem2D
ln -s ../../../bin/xspecfem2D
cd ../

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`
if [ "$NPROC" -ne 1 ]; then
  echo "Checkerboard example requires to run simulation with a single process only."
  echo "Please set NPROC to 1 in DATA/Par_file"
  exit 1
fi

# runs database generation
# This is a serial simulation
echo
echo "running mesher..."
echo
./bin/xmeshfem2D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
# This is a serial simulation
echo
echo "running solver..."
echo
./bin/xspecfem2D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# stores output
cp DATA/SOURCE OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/
#cp DATA/STATIONS_checker OUTPUT_FILES/

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date

