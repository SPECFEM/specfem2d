#!/bin/bash -eu
#################################
# USER PARAMETER

NPROC=20

#################################

echo "setting up example: `pwd`"
echo

mkdir -p OUTPUT_FILES
# cleanup
rm -rf OUTPUT_FILES/*

# links executables
mkdir -p bin
cd bin/
rm -f xmeshfem2D xspecfem2D
ln -s ../../../bin/xmeshfem2D
ln -s ../../../bin/xspecfem2D
cd ../

echo
echo "running mesher..."
echo
./bin/xmeshfem2D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "running solver..."
echo
mpirun -np $NPROC ./bin/xspecfem2D
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo

