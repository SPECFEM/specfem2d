#!/bin/bash
#
# test script to run a kernel simulation with a acoustic/elastic coupled model
#

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directoy
echo
echo "setting up example..."
echo

mkdir -p OUTPUT_FILES

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

# helper scripts
ln -s ../../utils/change_simulation_type.pl
ln -s ../../utils/Visualization/plot_kernel.py

# backup
cp -v DATA/Par_file DATA/Par_file.bak

# simulation setup
# switches to Stacey condition to be able to run also on GPU (for testing purposes)
sed -i "s:^PML_BOUNDARY_CONDITIONS .*:PML_BOUNDARY_CONDITIONS = .false.:" DATA/Par_file
sed -i "s:^STACEY_ABSORBING_CONDITIONS .*:STACEY_ABSORBING_CONDITIONS = .true.:" DATA/Par_file

# reduce image outputs
sed -i "s:^output_postscript_snapshot .*:output_postscript_snapshot = .false.:" DATA/Par_file

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

####################################################################
##
## forward simulation
##
####################################################################

## forward run
./change_simulation_type.pl -F

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running mesher..."
  echo
  ./bin/xmeshfem2D
else
  # This is a MPI simulation
  echo
  echo "running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xmeshfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

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
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`

####################################################################
##
## creating adjoint sources
##
####################################################################

# for testing purpose: uses the full forward traces as adjoint traces

# adjoint setup
mkdir -p SEM/
rm -rf SEM/*

if [ -f OUTPUT_FILES/AA.S0001.PRE.semp ]; then
# pressure
end=.semp
else
# displacement
end=.semd
fi

cp -v OUTPUT_FILES/AA.*$end SEM/
cd SEM/
rename "s/$end/.adj/" AA*$end
cd ../

# backup seismograms
cd OUTPUT_FILES/
mkdir -p forward_run
rm -rf forward_run/*
mv -v AA.*$end forward_run/
mv -v output_* forward_run/
cd ../

####################################################################
##
## adjoint/kernel simulation
##
####################################################################

## kernel run
./change_simulation_type.pl -b

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

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


echo
echo


# plots kernel image
if [ -e OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat ]; then
  # elastic kernel
  ./plot_kernel.py OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat
fi
if [ -e OUTPUT_FILES/proc000000_rhop_c_kernel.dat.dat ]; then
  # acoustic kernel
  ./plot_kernel.py OUTPUT_FILES/proc000000_rhop_c_kernel.dat.dat
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo

# restore backup version
cp -v DATA/Par_file.bak DATA/Par_file

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`


