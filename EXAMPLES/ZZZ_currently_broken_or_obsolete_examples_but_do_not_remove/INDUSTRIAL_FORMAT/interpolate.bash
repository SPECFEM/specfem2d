#!/bin/bash -eu

DIR_SEM="../../.."

cp $DIR_SEM/bin/xspecfem2D $DIR_SEM/bin/xmeshfem2D ./

rm -r DATA
rm -r OUTPUT_FILES

mkdir DATA
mkdir OUTPUT_FILES

cp Par_file                      ./DATA/
cp SOURCE                        ./DATA/

gfortran interpolate.f90 -o xinterpolate -lm -O4

./xinterpolate
