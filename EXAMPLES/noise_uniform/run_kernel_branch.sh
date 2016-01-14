#!/bin/sh

if [ $# -ne 3 ]; then echo "USAGE: ./run_branch.sh  slave-station-id  master-station-id branch(0=negative/1=positive)"; exit; fi

#############################################

# input parameters
# slave station id
SLAVE=$1

# master station id
MASTER=$2

# cross-correlation branch
BRANCH=$3

############################################

if [ "$BRANCH" != "0" ] && [ "$BRANCH" != "1" ]; then echo "Invalid branch number, must be 0 for negative or 1 for positive branch"; exit 1; fi

# define your Fortran compiler here
FC="gfortran"

TRACE=`printf 'AA.S%04d.BXY.semd' $SLAVE`
ADJCC='adj_cc.f90'

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
#rm -rf   DATA SEM OUTPUT_FILES
mkdir -p DATA
mkdir -p NOISE_TOMOGRAPHY
mkdir -p SEM
mkdir -p OUTPUT_FILES

rm -rf OUTPUT_ALL
mkdir -p OUTPUT_ALL

# sets up local DATA/ directory
cp -v SOURCE_noise    DATA/SOURCE
cp -v STATIONS_noise  DATA/STATIONS
cp -v uniform.dat     DATA/
echo $MASTER > NOISE_TOMOGRAPHY/irec_master_noise
# noise source
if [ -f S_squared ]; then cp -v S_squared NOISE_TOMOGRAPHY/; fi
# velocity model
if [ -f model_velocity.dat_input ]; then cp -v model_velocity.dat_input DATA/; fi

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
echo
echo "noise simulation 1"
echo
cp -v Par_file_noise_1  DATA/Par_file
echo

./xmeshfem2D
./xspecfem2D

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
echo
echo "noise simulation 2"
echo
cp -v Par_file_noise_2  DATA/Par_file
echo

./xmeshfem2D
./xspecfem2D

# backup
mkdir -p OUTPUT_ALL/step_2
mv OUTPUT_FILES/image*       OUTPUT_ALL/step_2
cp OUTPUT_FILES/*.semd       OUTPUT_ALL/step_2
mv DATA/Par_file             OUTPUT_ALL/step_2

##
## adjoint source
##
echo
echo "creating adjoint source"
echo
if [ ! -f OUTPUT_FILES/$TRACE ]; then echo "trace file OUTPUT_FILES/$TRACE is missing"; exit 1; fi

# write zero traces
awk '{printf(" %12.6f %12.6f\n",$1,0.)}' < OUTPUT_FILES/$TRACE > SEM/zero
cd SEM/
for ((ii=1; ii<=3; ++ii))
do
  cp zero `printf AA.S%04d.BXX.adj $ii`
  cp zero `printf AA.S%04d.BXY.adj $ii`
  cp zero `printf AA.S%04d.BXZ.adj $ii`
done
cd ../

# compile and write master trace
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

echo
rename '.semd' '' SEM/$TRACE.adj


##
## simulation 3
##
echo
echo "noise simulation 3"
echo
cp -v Par_file_noise_3  DATA/Par_file
echo

./xmeshfem2D
./xspecfem2D


# kernel visualization
if [ -f OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat ]; then
gnuplot plot_kernel.gnu
mv -v image_rho_kappa_mu_kernels.png OUTPUT_FILES/
fi

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
echo
echo "done"
date

