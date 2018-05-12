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

# use this to test a source that is located at a corner point shared between four spectral elements
cp ../Par_file_attenuation_2D_source_in_point_shared_by_four_elements Par_file
cp ../SOURCE_attenuation_2D_source_in_point_shared_by_four_elements SOURCE

# or use this to test a source that is located inside a given element, which does NOT coincide with a GLL point
# cp ../Par_file_attenuation_2D_source_in_point_inside_an_element Par_file
# cp ../SOURCE_attenuation_2D_source_in_point_inside_an_element SOURCE

cp ../interfaces_attenuation_analytic.dat .
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

cd $currentdir

# links executables
rm -f xmeshfem2D xspecfem2D
ln -s ../../../bin/xmeshfem2D
ln -s ../../../bin/xspecfem2D

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
# use this to test a source that is located at a corner point shared between four spectral elements
cp ../Par_file_attenuation_2D_source_in_point_shared_by_four_elements Par_file
cp ../SOURCE_attenuation_2D_source_in_point_shared_by_four_elements SOURCE

cp ../interfaces_attenuation_analytic.dat .
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

cd $currentdir

# links executables
rm -f xmeshfem2D xspecfem2D
ln -s ../../../bin/xmeshfem2D
ln -s ../../../bin/xspecfem2D

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
