#!/bin/bash -eu
#
# script runs mesher and solver indirectly through the program 'interpolate'
#

# gets the default compiler
f90=`grep "^FC " ../../Makefile | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# selects compiler
#f90=gfortran
#f90=ifort
flags="-O3"  # debug: -Wall -g -fbounds-check


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

# backup
cp -v DATA/Par_file DATA/Par_file.bak


echo
echo "compiling xinterpolate..."
echo
$f90 $flags -o bin/xinterpolate interpolate.f90
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs database generation and forward simulations needed
echo
echo "running interpolation..."
echo
./bin/xinterpolate
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# restores original Par_file
cp -v DATA/Par_file.bak DATA/Par_file

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`

