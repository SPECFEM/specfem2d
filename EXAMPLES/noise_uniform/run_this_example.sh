#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#
##########################################

# select main station #id
#main_id=1

##########################################

echo "running example: `date`"
currentdir=`pwd`

# sets up directory structure in current example directoy
echo
echo "   setting up example..."
echo

mkdir -p OUTPUT_FILES
mkdir -p NOISE_TOMOGRAPHY

# sets up local DATA/ directory
if [ "$1" == "" ]; then
cp -v DATA/Par_file DATA/Par_file_noise_1
else
echo "using noise Par_file $1"
step=$1
case $step in
1) sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 1:" DATA/Par_file
   sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 1:" DATA/Par_file
   sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .false.:" DATA/Par_file
   ;;
2) sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 1:" DATA/Par_file
   sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 2:" DATA/Par_file
   sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .true.:" DATA/Par_file
   ;;
3) sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 3:" DATA/Par_file
   sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 3:" DATA/Par_file
   sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .false.:" DATA/Par_file
   ;;
*) echo "step not recognized: $step"; echo "please use as step number 1, 2 or 3"; exit 1 ;;
esac
cp -v DATA/Par_file DATA/Par_file_noise_${step}
echo
fi

# sets a main station (or takes the default one from the existing file in NOISE_TOMOGRAPHY/)
#echo "main id: $main_id"
#echo $main_id > NOISE_TOMOGRAPHY/irec_main_noise

# cleans output files
rm -rf OUTPUT_FILES/*

# links executables
mkdir -p bin
cd bin/
rm -f xmeshfem2D xspecfem2D
ln -s ../../../bin/xmeshfem2D
ln -s ../../../bin/xspecfem2D
cd ../

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

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
