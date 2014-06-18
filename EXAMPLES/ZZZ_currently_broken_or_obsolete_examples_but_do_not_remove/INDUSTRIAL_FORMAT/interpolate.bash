#!/bin/bash -eu

DIR_SEM="../../"

cp Par_file                      $DIR_SEM/DATA/
cp interfaces_M2_UPPA_curved.dat $DIR_SEM/DATA/
cp SOURCE                        $DIR_SEM/DATA/
cp constants.h                   $DIR_SEM/setup/

ifort interpolate.f90 -o $DIR_SEM/xinterpolate
cd $DIR_SEM

make
cp bin/* ./
./xinterpolate
