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
mkdir -p SEM

# cleans output files (only if a command line argument is NOT given, because the second time when we run the adjoint run we should not erase anything)
noargument=0
if [ "$#" -eq "$noargument" ]
then
  echo "erasing all database files"
  rm -rf OUTPUT_FILES/*
else
  echo "not erasing any of the database files"
fi

cd $currentdir

# links executables
mkdir -p bin
cd bin/
rm -f xmeshfem2D xspecfem2D
ln -s ../../../bin/xmeshfem2D
ln -s ../../../bin/xspecfem2D
cd ../

# stores setup
cp -f DATA/Par_file OUTPUT_FILES
cp -f DATA/SOURCE OUTPUT_FILES

# runs database generation
echo
echo "  running mesher..."
echo
./bin/xmeshfem2D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
echo
echo "  running solver..."
echo
./bin/xspecfem2D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# stores output
cp -f DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
