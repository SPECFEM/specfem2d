#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#
##########################################

# master station #id
master_id=1

##########################################

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p OUTPUT_FILES
mkdir -p DATA
mkdir -p NOISE_TOMOGRAPHY

# sets up local DATA/ directory
cd DATA/
cp ../Par_file_noise_1 Par_file
cp ../SOURCE_noise SOURCE
cp ../STATIONS_noise STATIONS
cp ../uniform.dat ./
cd ../

# sets a master station
echo $master_id > NOISE_TOMOGRAPHY/irec_master_noise

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

# runs database generation
echo
echo "  running mesher..."
echo
./xmeshfem2D

# runs simulation
echo
echo "  running solver..."
echo
./xspecfem2D

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
