#!/bin/bash -eu

#############################################
## USER PARAMETERS

# input parameters
# secondary station id
SECONDARY=$1

# main station id
MAIN=$2

# cross-correlation branch
BRANCH=$3

# define your Fortran compiler here
FC=gfortran         # or: ifort

############################################

if [ $# -ne 3 ]; then
 echo "USAGE: ./run_kernel_branch.sh  secondary-station-id  main-station-id branch(0=negative/1=positive)"
 exit 1
fi

if [ "$BRANCH" != "0" ] && [ "$BRANCH" != "1" ]; then
  echo "Invalid branch number, must be 0 for negative or 1 for positive branch"
  exit 1
fi


echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo
echo "main: $MAIN"
echo "secondary : $SECONDARY"
if [ "$BRANCH" == "0" ]; then echo "using negative branch"; fi
if [ "$BRANCH" == "1" ]; then echo "using positive branch"; fi

# prepare directories
mkdir -p NOISE_TOMOGRAPHY
mkdir -p SEM
mkdir -p OUTPUT_FILES

rm -rf OUTPUT_ALL
mkdir -p OUTPUT_ALL

# sets up local DATA/ directory
echo "$MAIN" > NOISE_TOMOGRAPHY/irec_main_noise

# noise source
#if [ -f S_squared ]; then cp -v S_squared NOISE_TOMOGRAPHY/; fi
# velocity model
#if [ -f model_velocity.dat_input ]; then cp -v model_velocity.dat_input DATA/; fi

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

##
## simulation 1
##
echo "###############################################"
echo "##                                           ##"
echo "## noise simulation: step 1                  ##"
echo "##                                           ##"
echo "###############################################"
# Par_file settings
sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 1:" DATA/Par_file
sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 1:" DATA/Par_file
sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .false.:" DATA/Par_file
cp -v DATA/Par_file  DATA/Par_file_noise_1
echo

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running mesher..."
  echo
  ./bin/xmeshfem2D > OUTPUT_FILES/output_mesher.log
else
  # This is a MPI simulation
  echo
  echo "running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xmeshfem2D > OUTPUT_FILES/output_mesher.log
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver..."
  echo
  ./bin/xspecfem2D > OUTPUT_FILES/output_solver.log
else
  # This is a MPI simulation
  echo
  echo "running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem2D > OUTPUT_FILES/output_solver.log
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup
mkdir -p OUTPUT_ALL/step_1

echo
mv -v OUTPUT_FILES/*image*.jpg  OUTPUT_ALL/step_1
mv -v OUTPUT_FILES/*.semd       OUTPUT_ALL/step_1
mv -v OUTPUT_FILES/plot*        OUTPUT_ALL/step_1
mv -v OUTPUT_FILES/output*.log  OUTPUT_ALL/step_1
cp -v DATA/Par_file             OUTPUT_ALL/step_1

mv -v OUTPUT_FILES/mask*        OUTPUT_ALL/
mv -v OUTPUT_FILES/mesh_????    OUTPUT_ALL/
mv -v OUTPUT_FILES/model*       OUTPUT_ALL/
echo

##
## simulation 2
##
echo "###############################################"
echo "##                                           ##"
echo "## noise simulation: step 2                  ##"
echo "##                                           ##"
echo "###############################################"
# Par_file settings
sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 1:" DATA/Par_file
sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 2:" DATA/Par_file
sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .true.:" DATA/Par_file
cp -v DATA/Par_file  DATA/Par_file_noise_2
echo

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running mesher..."
  echo
  ./bin/xmeshfem2D > OUTPUT_FILES/output_mesher.log
else
  # This is a MPI simulation
  echo
  echo "running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xmeshfem2D > OUTPUT_FILES/output_mesher.log
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver..."
  echo
  ./bin/xspecfem2D > OUTPUT_FILES/output_solver.log
else
  # This is a MPI simulation
  echo
  echo "running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem2D > OUTPUT_FILES/output_solver.log
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup
mkdir -p OUTPUT_ALL/step_2

echo
mv -v OUTPUT_FILES/*image*.jpg  OUTPUT_ALL/step_2
cp -v OUTPUT_FILES/*.semd       OUTPUT_ALL/step_2
mv -v OUTPUT_FILES/output*.log  OUTPUT_ALL/step_2
cp -v DATA/Par_file             OUTPUT_ALL/step_2
echo

##
## adjoint source
##
echo "###############################################"
echo "##                                           ##"
echo "## creating adjoint source                   ##"
echo "##                                           ##"
echo "###############################################"

TRACE=`printf 'AA.S%04d.BXY.semd' $SECONDARY`
TRACE_ADJ=`printf 'AA.S%04d.BXY.adj' $SECONDARY`

echo "using trace: OUTPUT_FILES/$TRACE"
if [ ! -f OUTPUT_FILES/$TRACE ]; then echo "trace file OUTPUT_FILES/$TRACE is missing"; exit 1; fi

# clean
rm -f SEM/*

# write zero traces
awk '{printf(" %12.6f %12.6f\n",$1,0.0)}' < OUTPUT_FILES/$TRACE > SEM/zero
cd SEM/
for ((ii=1; ii<=3; ++ii))
do
  if [ "$ii" -ne "$SECONDARY" ]; then
    #cp zero `printf AA.S%04d.BXX.adj $ii`
    cp -v zero `printf AA.S%04d.BXY.adj $ii`
    #cp zero `printf AA.S%04d.BXZ.adj $ii`
  fi
done
cd ../

# compile and write main trace
ADJCC='adj_traveltime_filter.f90'
if [ "$BRANCH" == "0" ]; then
  echo "negative branch"
  sed -i'.bak' 's/branch_type = .*/branch_type = 0/' $ADJCC
elif [ "$BRANCH" == "1" ]; then
  echo "positive branch"
  sed -i'.bak' 's/branch_type = .*/branch_type = 1/' $ADJCC
fi
# no time reversal of signal
sed -i'.bak' 's/time_reverse = .[a-z]*./time_reverse = .false./' $ADJCC
echo ""

$FC $ADJCC -o xadj_run
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

cp OUTPUT_FILES/$TRACE SEM/

./xadj_run SEM/$TRACE
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
# note: using > rename .semd.adj .adj SEM/$TRACE.adj
#       might fails on different OS due to different rename command-line versions
cp -p -v SEM/$TRACE.adj SEM/$TRACE_ADJ
echo

##
## simulation 3
##
echo "###############################################"
echo "##                                           ##"
echo "## noise simulation: step 3                  ##"
echo "##                                           ##"
echo "###############################################"
# Par_file settings
sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 3:" DATA/Par_file
sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 3:" DATA/Par_file
sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .false.:" DATA/Par_file
cp -v DATA/Par_file  DATA/Par_file_noise_3
echo

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running mesher..."
  echo
  ./bin/xmeshfem2D > OUTPUT_FILES/output_mesher.log
else
  # This is a MPI simulation
  echo
  echo "running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xmeshfem2D > OUTPUT_FILES/output_mesher.log
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver..."
  echo
  ./bin/xspecfem2D > OUTPUT_FILES/output_solver.log
else
  # This is a MPI simulation
  echo
  echo "running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./bin/xspecfem2D > OUTPUT_FILES/output_solver.log
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup
mkdir -p OUTPUT_ALL/step_3

echo
cp -v OUTPUT_FILES/*image*.jpg  OUTPUT_ALL/step_3
cp -v OUTPUT_FILES/*.semd       OUTPUT_ALL/step_3
cp -v OUTPUT_FILES/output*.log  OUTPUT_ALL/step_3
cp -v SEM/*Y.adj                OUTPUT_ALL/step_3
cp -v DATA/Par_file             OUTPUT_ALL/step_3
# kernels
cp -v OUTPUT_FILES/proc*        OUTPUT_ALL/
echo

echo
echo "see results in directory: OUTPUT_ALL/"
echo "for plotting kernels, run script: ./plot_kernel.sh"
echo "done: `date`"
echo

