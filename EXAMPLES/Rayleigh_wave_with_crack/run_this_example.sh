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
# cleans output files
rm -rf OUTPUT_FILES/*

echo
echo "*********************************************************************************"
echo "This example requires the setting in the source code file setup/constants.h:"
echo "! select fast (Paul Fischer) or slow (topology only) global numbering algorithm"
echo "  logical, parameter :: FAST_NUMBERING = .false."
echo "  .."
echo "! add a small crack (discontinuity) in the medium manually"
echo "  logical, parameter :: ADD_A_SMALL_CRACK_IN_THE_MEDIUM = .true."
echo "*********************************************************************************"
echo
echo "will change settings and recompile the executables..."

rootdir=../../
cd $rootdir

# checks directory
if [ ! -e setup/constants.h ]; then
  echo "Please make sure to run this script within the EXAMPLES/Rayleigh_wave_with_crack/ directory..."
  exit 1
fi

# backup
cp -v setup/constants.h setup/constants.h.bak

# changes settings for example
sed -i "s:FAST_NUMBERING .*:FAST_NUMBERING = .false.:" setup/constants.h
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

sed -i "s:ADD_A_SMALL_CRACK_IN_THE_MEDIUM .*:ADD_A_SMALL_CRACK_IN_THE_MEDIUM = .true.:" setup/constants.h
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "set new constants.h"
echo "re-compile executables..."
echo

# recompiles
make clean  > $currentdir/tmp_compilation.log 2>&1
make -j4 all >> $currentdir/tmp_compilation.log 2>&1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# example directory
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

# Get the number of processors from Par_file
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
echo
echo "  running mesher..."
echo
./bin/xmeshfem2D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
if [ "$NPROC" -eq 1 ]; then # This is a serial simulation
  echo
  echo " Running solver..."
  echo
  ./bin/xspecfem2D
else # This is a MPI simulation
  echo
  echo " Running solver on $NPROC processors..."
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
echo "done"
date

echo
echo "will reset setup/constants.h to original setup (before running this script) and re-compile executables..."
echo
cd $rootdir

mv -v setup/constants.h.bak setup/constants.h
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "set original constants.h (reverted back)"
echo "re-compile executables..."
echo

# recompiles
make clean  > $currentdir/tmp_compilation.log 2>&1
make -j4 all >> $currentdir/tmp_compilation.log 2>&1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "original setup & state done"
echo
echo "done"
date

