#!/bin/bash
#
# script runs mesher and solver
# using this example setup
#

echo "running example: `date`"
currentdir=`pwd`

echo
echo "(will take about 5 minutes)"
echo
cd $currentdir

rm -f change_simulation_type.pl
ln -s ../../utils/change_simulation_type.pl

# backup
cp DATA/Par_file DATA/Par_file.bak

# changes seismogram output to displacement (for adjoint source creation)
sed -i "s:^seismotype .*:seismotype                      = 1:g" DATA/Par_file

# first P-arrival only
#sed -i "s:^NSTEP .*:NSTEP                           = 1500:g" DATA/Par_file

# image output
#sed -i "s:^output_color_image .*:output_color_image              = .false.:g" DATA/Par_file

# ABC
#sed -i "s:^STACEY_ABSORBING_CONDITIONS .*:STACEY_ABSORBING_CONDITIONS     = .false.:g" DATA/Par_file


##
## 1. forward simulation
##
echo "############################################################"
echo
echo "running forward simulation (with saving forward wavefield)"
echo
echo "############################################################"
./change_simulation_type.pl -F

./run_this_example.sh
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


##
## adjoint sources
##
echo "############################################################"
echo
echo "setting up adjoint sources"
echo
echo "############################################################"

# station
station_A=AA.S0001

echo "using station: ${station_A}"

# check needs traces
if [ ! -e OUTPUT_FILES/${station_A}.BXZ.semd ]; then echo "no displacement traces in OUTPUT_FILES/, please run forward simulation first..."; exit 1; fi

# setup adjoint sources directory
mkdir -p SEM
rm -rf SEM/*

## generates adjoint sources
# links executable
cd bin/
rm -f xadj_seismogram
ln -s ../../../bin/xadj_seismogram
cd ../

echo "############################################################"
echo
echo "running adjoint source creation"
echo
echo "############################################################"
# adjoint source for X&Z-component
# first P-arrival
./bin/xadj_seismogram 0.015 0.095 ${station_A} 4

# or slow P-arrival:
#./bin/xadj_seismogram 0.12 0.18 ${station_A}  4

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# zero out other components (optional)
#awk '{print $1,0.0}' SEM/${station_A}.BXZ.adj > SEM/${station_A}.BXX.adj

# stores forward output
rm -rf OUTPUT_FILES.for
cp -rp OUTPUT_FILES OUTPUT_FILES.for

##
## 2. kernel simulation
##
echo "############################################################"
echo
echo "running kernel simulation"
echo
echo "############################################################"
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
## plots a kernel
##
echo "############################################################"
echo
echo "plotting kernels"
echo
echo "############################################################"

# GMT plotting
#rm -f plot_wavefield.pl
#ln -s ../../utils/Visualization/plot_wavefield.pl
#./plot_wavefield.pl 400/3000/400 400/2800/800 0/480/0/480 120/20/120/20 -8/-8/-8 1/1/1 -48.0/0.06 0/1/0/0/0 2.0/1/0 1/0/1/120 Tape2007_kernel onerec_homo

# matplotlib plotting
rm -f plot_kernel.py
ln -s ../../utils/Visualization/plot_kernel.py

# kernels
names=( cpI_cpII_cs_kernel mu_B_C_kernel mub_Bb_Cb_kernel M_rho_rhof_kernel Mb_rhob_rhofb_kernel)

for name in ${names[@]}; do
  # combines kernel outputs from all processes
  cat OUTPUT_FILES/proc000*_${name}.dat > OUTPUT_FILES/${name}.dat
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

  # plots kernel
  ./plot_kernel.py OUTPUT_FILES/${name}.dat
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
done

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
date
