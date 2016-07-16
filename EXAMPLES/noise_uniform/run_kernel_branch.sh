#!/bin/bash -eu

#############################################
## USER PARAMETERS

# input parameters
# slave station id
SLAVE=$1

# master station id
MASTER=$2

# cross-correlation branch
BRANCH=$3

# define your Fortran compiler here
FC=gfortran         # or: ifort

############################################

if [ $# -ne 3 ]; then
 echo "USAGE: ./run_kernel_branch.sh  slave-station-id  master-station-id branch(0=negative/1=positive)"
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
echo "master: $MASTER"
echo "slave : $SLAVE"
if [ "$BRANCH" == "0" ]; then echo "using negative branch"; fi
if [ "$BRANCH" == "1" ]; then echo "using positive branch"; fi

# prepare directories
mkdir -p NOISE_TOMOGRAPHY
mkdir -p SEM
mkdir -p OUTPUT_FILES

rm -rf OUTPUT_ALL
mkdir -p OUTPUT_ALL

# sets up local DATA/ directory
echo "$MASTER" > NOISE_TOMOGRAPHY/irec_master_noise

# noise source
#if [ -f S_squared ]; then cp -v S_squared NOISE_TOMOGRAPHY/; fi
# velocity model
#if [ -f model_velocity.dat_input ]; then cp -v model_velocity.dat_input DATA/; fi

# cleans output files
rm -rf OUTPUT_FILES/*

cd $currentdir

# links executables
rm -f xmeshfem2D xspecfem2D
ln -s ../../bin/xmeshfem2D
ln -s ../../bin/xspecfem2D

##
## simulation 1
##
echo "###############################################"
echo "##                                           ##"
echo "## noise simulation: step 1                  ##"
echo "##                                           ##"
echo "###############################################"
cp -v DATA/Par_file_noise_1  DATA/Par_file
echo

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running mesher..."
  echo
  ./xmeshfem2D
else
  # This is a MPI simulation
  echo
  echo "running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC ./xmeshfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver..."
  echo
  ./xspecfem2D
else
  # This is a MPI simulation
  echo
  echo "running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./xspecfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup
mkdir -p OUTPUT_ALL/step_1
mv OUTPUT_FILES/image*       OUTPUT_ALL/step_1
mv OUTPUT_FILES/*.semd       OUTPUT_ALL/step_1
mv OUTPUT_FILES/plot*        OUTPUT_ALL/step_1
mv DATA/Par_file             OUTPUT_ALL/step_1

mv OUTPUT_FILES/mask*        OUTPUT_ALL/
mv OUTPUT_FILES/mesh_????    OUTPUT_ALL/
mv OUTPUT_FILES/model*       OUTPUT_ALL/

##
## simulation 2
##
echo "###############################################"
echo "##                                           ##"
echo "## noise simulation: step 2                  ##"
echo "##                                           ##"
echo "###############################################"
cp -v DATA/Par_file_noise_2  DATA/Par_file
echo

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running mesher..."
  echo
  ./xmeshfem2D
else
  # This is a MPI simulation
  echo
  echo "running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC ./xmeshfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver..."
  echo
  ./xspecfem2D
else
  # This is a MPI simulation
  echo
  echo "running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./xspecfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup
mkdir -p OUTPUT_ALL/step_2
mv OUTPUT_FILES/image*       OUTPUT_ALL/step_2
cp OUTPUT_FILES/*.semd       OUTPUT_ALL/step_2
mv DATA/Par_file             OUTPUT_ALL/step_2

##
## adjoint source
##
echo "###############################################"
echo "##                                           ##"
echo "## creating adjoint source                   ##"
echo "##                                           ##"
echo "###############################################"

TRACE=`printf 'AA.S%04d.BXY.semd' $SLAVE`
TRACE_ADJ=`printf 'AA.S%04d.BXY.adj' $SLAVE`

echo "using trace: OUTPUT_FILES/$TRACE"
if [ ! -f OUTPUT_FILES/$TRACE ]; then echo "trace file OUTPUT_FILES/$TRACE is missing"; exit 1; fi

# clean
rm -f SEM/*

# write zero traces
awk '{printf(" %12.6f %12.6f\n",$1,0.0)}' < OUTPUT_FILES/$TRACE > SEM/zero
cd SEM/
for ((ii=1; ii<=3; ++ii))
do
  if [ "$ii" -ne "$SLAVE" ]; then
    #cp zero `printf AA.S%04d.BXX.adj $ii`
    cp -v zero `printf AA.S%04d.BXY.adj $ii`
    #cp zero `printf AA.S%04d.BXZ.adj $ii`
  fi
done
cd ../

# compile and write master trace
ADJCC='adj_cc.f90'
if [ "$BRANCH" == "0" ]; then
  echo "negative branch"
  sed -i'.bak' 's/use_positive_branch = .[a-z]*./use_positive_branch = .false./' $ADJCC
  sed -i'.bak' 's/use_negative_branch = .[a-z]*./use_negative_branch = .true./' $ADJCC
elif [ "$BRANCH" == "1" ]; then
  echo "positive branch"
  sed -i'.bak' 's/use_positive_branch = .[a-z]*./use_positive_branch = .true./' $ADJCC
  sed -i'.bak' 's/use_negative_branch = .[a-z]*./use_negative_branch = .false./' $ADJCC
fi
sed -i'.bak' 's/use_custom_window = .[a-z]*./use_custom_window = .false./' $ADJCC
sed -i'.bak' 's/time_reverse = .[a-z]*./time_reverse = .false./' $ADJCC


$FC $ADJCC -o adj_run
cp OUTPUT_FILES/$TRACE SEM/

./adj_run SEM/$TRACE
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
# note: using > rename .semd.adj .adj SEM/$TRACE.adj
#       might fails on different OS due to different rename command-line versions
cp -p -v SEM/$TRACE.adj SEM/$TRACE_ADJ


##
## simulation 3
##
echo "###############################################"
echo "##                                           ##"
echo "## noise simulation: step 3                  ##"
echo "##                                           ##"
echo "###############################################"
cp -v DATA/Par_file_noise_3  DATA/Par_file
echo

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running mesher..."
  echo
  ./xmeshfem2D
else
  # This is a MPI simulation
  echo
  echo "running mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC ./xmeshfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver..."
  echo
  ./xspecfem2D
else
  # This is a MPI simulation
  echo
  echo "running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./xspecfem2D
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# backup
mkdir -p OUTPUT_ALL/step_3
cp OUTPUT_FILES/image*       OUTPUT_ALL/step_3
cp OUTPUT_FILES/*.semd       OUTPUT_ALL/step_3
cp SEM/*Y.adj                OUTPUT_ALL/step_3
cp DATA/Par_file             OUTPUT_ALL/step_3
# kernels
cp OUTPUT_FILES/proc*        OUTPUT_ALL/


echo
echo "see results in directory: OUTPUT_ALL/"
echo "for plotting kernels, run script: ./plot_kernel.sh"
echo "done: `date`"
echo

