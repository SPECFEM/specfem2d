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

cd $currentdir

rm -f change_simulation_type.pl
ln -s ../../utils/change_simulation_type.pl

##
## 1. forward simulation
##
echo
echo "running forward simulation (with saving forward wavefield)"
echo
./change_simulation_type.pl -F

./run_this_example.sh
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


##
## adjoint sources
##
echo
echo "setting up adjoint sources"
echo

# setup adjoint sources directory
mkdir -p SEM
rm -rf SEM/*

# check needs traces
# SH default case
if [ ! -e OUTPUT_FILES/AA.S0001.BXY.semd ]; then echo "no traces in OUTPUT_FILES/, please run forward simulation first..."; exit 1; fi

## execute this in case a local file is used instead of default in directory src/auxiliaries/
if [ "$1" == "1" ]; then
  my_file="adj_seismogram_Tape2007.f90"

  # fortran compiler (as specified in Makefile)
  FC=`grep '^FC .*' ../../Makefile | cut -d = -f 2 | sed "s/^[ \t]*//"`
  if [ "$FC" == "" ]; then echo "fortran compiler not found, exiting..."; exit 1; fi
  echo "using compiler: $FC"
  echo

  # creates adjoint sources
  rm -rf bin/xadj_seismogram
  $FC $my_file -o bin/xadj_seismogram

  echo
  echo "running adjoint source creation"
  echo
  ./bin/xadj_seismogram
else
  # links executable
  cd bin/
  rm -f xadj_seismogram
  ln -s ../../../bin/xadj_seismogram
  cd ../

  echo
  echo "running adjoint source creation"
  echo
  ./bin/xadj_seismogram 95.0 130.0  AA.S0001  2
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# stores forward output
rm -rf OUTPUT_FILES.for
cp -rp OUTPUT_FILES OUTPUT_FILES.for

##
## 2. kernel simulation
##
echo
echo "running kernel simulation"
echo
./change_simulation_type.pl -b

# runs database generation
# In principle we do not need rerun xmeshfem2D in the adjoint run.
# However, in the current structure of the code, the xspecfem2D program can not detect the
# the parameter change in Par_file. Since in the adjoint run we change the SIMULATION_TYPE and the save_forward
# Get the number of processors

./run_this_example.sh noclean
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

##
## plot the traveltime kernel
##

# GMT plotting
#rm -f plot_wavefield.pl
#ln -s ../../utils/Visualization/plot_wavefield.pl
#./plot_wavefield.pl 400/3000/400 400/2800/800 0/480/0/480 120/20/120/20 -8/-8/-8 1/1/1 -48.0/0.06 0/1/0/0/0 2.0/1/0 1/0/1/120 Tape2007_kernel onerec_homo

# matplotlib plotting
rm -f plot_kernel.py
ln -s ../../utils/Visualization/plot_kernel.py

cat OUTPUT_FILES/proc000*_rho_kappa_mu_kernel.dat > OUTPUT_FILES/rho_kappa_mu_kernel.dat

./plot_kernel.py OUTPUT_FILES/rho_kappa_mu_kernel.dat

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
