#!/bin/bash -eu

# selects compiler
f90=gfortran
#f90=ifort
flags="-lm -O4"


echo "setting up example..."

DIR_SEM="../.."

cp -v $DIR_SEM/bin/xspecfem2D $DIR_SEM/bin/xmeshfem2D ./

mkdir -p DATA
mkdir -p OUTPUT_FILES

# cleanup
rm -rf DATA/*
rm -rf OUTPUT_FILES/*

cp -v Par_file                      ./DATA/
cp -v SOURCE                        ./DATA/


echo
echo "compiling xinterpolate..."
echo
$f90 $flags -o xinterpolate interpolate.f90

echo
echo "running interpolation..."
echo
./xinterpolate

echo
echo "done"
echo

