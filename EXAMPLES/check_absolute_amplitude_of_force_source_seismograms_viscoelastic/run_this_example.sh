#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p OUTPUT_FILES
mkdir -p DATA

# sets up local DATA/ directory
cd DATA/

# use this to test a source that is at a GLL point in the corner of an element, shared by several elements and assembled in the mass matrix
cp -f ./Par_file_no_attenuation_2D_at_the_corner_between_several_spectral_elements Par_file
cp -f ./SOURCE_no_attenuation_2D_at_the_corner_between_several_spectral_elements SOURCE

# or use this to test a source that is inside a given spectral element
# cp -f ../Par_file_no_attenuation_2D_inside_a_given_spectral_element Par_file
# cp -f ../SOURCE_no_attenuation_2D_inside_a_given_spectral_element SOURCE

#cp ../interfaces_attenuation_analytic.dat .
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

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
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
