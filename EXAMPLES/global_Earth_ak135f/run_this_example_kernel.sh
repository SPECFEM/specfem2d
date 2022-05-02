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

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

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
if [ ! -e OUTPUT_FILES/AA.S0001.BXZ.semd ]; then echo "no traces in OUTPUT_FILES/, please run forward simulation first..."; exit 1; fi

# waveforms as sources
#cp -v OUTPUT_FILES/AA.S*.semd SEM/
#cd SEM/
#rename 's/.semd/.adj/' *.semd
#cd ../

# traveltime sources
rm -f create_adjoint_source.py
ln -s ../../utils/adjoint_sources/create_adjoint_source.py

for file in OUTPUT_FILES/AA.*.semd
do
  ./create_adjoint_source.py 1 $file
done

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

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver..."
  echo
  ./bin/xspecfem2D
else
  # This is a MPI simulation
  echo
  echo "running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# stores output
cp -v DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

##
## plot the traveltime kernel
##
echo
echo "plotting..."
echo

# GMT plotting
#rm -f plot_wavefield.pl
#ln -s ../../utils/Visualization/plot_wavefield.pl
#./plot_wavefield.pl 400/3000/400 400/2800/800 0/480/0/480 120/20/120/20 -8/-8/-8 1/1/1 -48.0/0.06 0/1/0/0/0 2.0/1/0 1/0/1/120 Tape2007_kernel onerec_homo

# matplotlib plotting
rm -f plot_kernel.py
ln -s ../../utils/Visualization/plot_kernel.py

# elastic domain kernels
cat OUTPUT_FILES/proc000*_rho_kappa_mu_kernel.dat > OUTPUT_FILES/rho_kappa_mu_kernel.dat
./plot_kernel.py OUTPUT_FILES/rho_kappa_mu_kernel.dat 0 1.e-16

# acoustic domain kernels
cat OUTPUT_FILES/proc000*_rho_kappa_kernel.dat > OUTPUT_FILES/rho_kappa_kernel.dat
./plot_kernel.py OUTPUT_FILES/rho_kappa_kernel.dat 0 1.e-16

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done: `date`"
echo

