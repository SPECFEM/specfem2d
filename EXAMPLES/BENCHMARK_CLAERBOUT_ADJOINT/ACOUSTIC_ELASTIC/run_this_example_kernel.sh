#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#
echo "running example: `date`"
currentdir=`pwd`

cd $currentdir

##############################################
## Elastic benchmark

# Simulation type 1 == acoustic / 2 == elastic P-SV / 3 == elastic SH / 4 == coupled acoustic-elastic
SIM_TYPE=4

# perturbation model parameter rho/vp/vs (e.g. "rho" or "vp" or "rhovp")
perturb_param="vp"

# perturbation (should be small enough for approximating S(m - m0) ~ S(m) - S(m0)
perturb_percent=-0.02

# number of stations along x/z lines
nlinesx=1
nlinesz=1

##############################################

echo
echo "setup:"
echo "  SIM_TYPE                : $SIM_TYPE     (1 == acoustic / 2 == elastic P-SV / 3 == elastic SH / 4 == coupled acoustic-elastic)"
echo "  perturbation parameter  : $perturb_param"
echo "  perturbation percent    : $perturb_percent"
echo "  number of stations/lines: $nlinesx $nlinesz"
echo

# SH-wave simulation setup
if [ "$SIM_TYPE" == "3" ]; then
  # source
  sed -i "s:^source_type .*:source_type  = 1:" DATA/SOURCE
  sed -i "s:^time_function_type .*:time_function_type  = 1:" DATA/SOURCE
  sed -i "s:^factor .*:factor  = 1.d10:" DATA/SOURCE
  # Par_file
  sed -i "s:^P_SV .*:P_SV  = .false.:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP  = 900:" DATA/Par_file
  # checks param selection for SH-wave simulations
  if [ "$perturb_param" == "vp" ]; then
    # switch to vs perturbations instead
    echo "SH-wave simulation: switching perturbation parameter vp to vs"
    echo
    perturb_param="vs"
  fi
fi

# gets Par_file parameters
# Get the number of processors
NPROC=`grep '^NPROC ' DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`
NSTEP=`grep '^NSTEP ' DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`
DT=`grep '^DT ' DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

echo "Par_file parameters:"
echo "  NPROC = $NPROC"
echo "  NSTEP = $NSTEP"
echo "  DT    = $DT"
echo

# create and compile all setup
do_setup=1

##############################################

if [ "$do_setup" == "1" ]; then
echo
echo "setting up example..."
echo

# cleans files
rm -rf DATA/*.bin

mkdir -p OUTPUT_FILES
rm -rf OUTPUT_FILES/*
mkdir -p OUTPUT_FILES/

mkdir -p SEM
rm -rf SEM/*
mkdir -p SEM/dat SEM/syn

mkdir -p MODELS
rm -rf MODELS/*

mkdir -p MODELS/initial_model MODELS/target_model

mkdir -p KERNELS
rm -rf KERNELS/*

rm -f adj_seismogram.py
ln -s ../adj_seismogram.py
rm -f model_add_Gaussian_perturbation.py
ln -s ../model_add_Gaussian_perturbation.py
rm -f model_update.py
ln -s ../model_update.py
rm -f kernel_evaluation_postprocessing.py
ln -s ../kernel_evaluation_postprocessing.py
rm -f helper_functions.py
ln -s ../helper_functions.py

rm -f change_simulation_type.pl
ln -s ../../../utils/change_simulation_type.pl

rm -f create_STATIONS_file.py
ln -s ../create_STATIONS_file.py

# creates STATIONS file
./create_STATIONS_file.py $nlinesx $nlinesz

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

fi

########################### data ########################################

##
## data simulation
##

## forward simulation
echo
echo "running data forward simulation"
echo
./change_simulation_type.pl -f

# saving model files
sed -i "s:^SAVE_MODEL .*=.*:SAVE_MODEL = gll:" ./DATA/Par_file

./run_this_example.sh > output.log

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
mv -v output.log OUTPUT_FILES/

# backup copy
rm -rf OUTPUT_FILES.dat.forward
cp -rp OUTPUT_FILES OUTPUT_FILES.dat.forward

cp -v OUTPUT_FILES/*.su SEM/dat/

# target model
cp -v DATA/*rho.bin MODELS/target_model/
cp -v DATA/*vp.bin  MODELS/target_model/
cp -v DATA/*vs.bin  MODELS/target_model/

cp -v DATA/*NSPEC_ibool.bin  MODELS/target_model/
cp -v DATA/*x.bin  MODELS/target_model/
cp -v DATA/*z.bin  MODELS/target_model/

########################### model perturbation ################################

echo
echo "setting up perturbed model..."
echo "> ./model_add_Gaussian_perturbation.py $perturb_param $perturb_percent $NPROC "
echo
./model_add_Gaussian_perturbation.py $perturb_param $perturb_percent $NPROC

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# replaces model files with perturbed ones
for (( iproc = 0; iproc < $NPROC; iproc++ )); do
  rank=`printf "%06i\n" $iproc`
  cp -v DATA/proc${rank}_rho_gaussian.bin DATA/proc${rank}_rho.bin
  cp -v DATA/proc${rank}_vp_gaussian.bin DATA/proc${rank}_vp.bin
  cp -v DATA/proc${rank}_vs_gaussian.bin DATA/proc${rank}_vs.bin
  if [[ $? -ne 0 ]]; then exit 1; fi
done

########################### synthetics ################################

##
## synthetics simulation
##

echo
echo "running synthetics forward simulation (with saving forward wavefield)"
echo
./change_simulation_type.pl -F

# Par_file using GLL model
sed -i "s:^MODEL .*=.*:MODEL = gll:" ./DATA/Par_file
sed -i "s:^SAVE_MODEL .*=.*:SAVE_MODEL = .false.:" ./DATA/Par_file

./run_this_example.sh > output.log

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
mv -v output.log OUTPUT_FILES/output.forward.log

# backup copy
rm -rf OUTPUT_FILES.syn.forward
cp -rp OUTPUT_FILES OUTPUT_FILES.syn.forward

# initial model
mv -v OUTPUT_FILES/*.su SEM/syn/
cp -v DATA/*rho.bin MODELS/initial_model/
cp -v DATA/*vp.bin  MODELS/initial_model/
cp -v DATA/*vs.bin  MODELS/initial_model/

cp -v DATA/*NSPEC_ibool.bin  MODELS/initial_model/
cp -v DATA/*x.bin  MODELS/initial_model/
cp -v DATA/*z.bin  MODELS/initial_model/

########################### adj sources ################################
## adjoint sources
echo
echo "creating adjoint sources..."

# x-component
if [ -e OUTPUT_FILES.syn.forward/Ux_file_single_d.su ]; then
  syn=OUTPUT_FILES.syn.forward/Ux_file_single_d.su
  dat=OUTPUT_FILES.dat.forward/Ux_file_single_d.su
  echo "> ./adj_seismogram.py $syn $dat"
  echo
  # adjoint source f^adj = (s - d)
  ./adj_seismogram.py $syn $dat
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
fi

# y-component
if [ -e OUTPUT_FILES.syn.forward/Uy_file_single_d.su ]; then
  syn=OUTPUT_FILES.syn.forward/Uy_file_single_d.su
  dat=OUTPUT_FILES.dat.forward/Uy_file_single_d.su
  echo "> ./adj_seismogram.py $syn $dat"
  echo
  # adjoint source f^adj = (s - d)
  ./adj_seismogram.py $syn $dat
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
fi

# z-component
if [ -e OUTPUT_FILES.syn.forward/Uz_file_single_d.su ]; then
  syn=OUTPUT_FILES.syn.forward/Uz_file_single_d.su
  dat=OUTPUT_FILES.dat.forward/Uz_file_single_d.su
  echo "> ./adj_seismogram.py $syn $dat"
  echo
  # adjoint source f^adj = (s - d)
  ./adj_seismogram.py $syn $dat
  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi
fi

########################### kernel ################################

## kernel simulation
echo
echo "running kernel simulation"
echo
./change_simulation_type.pl -b

# In principle we do not need rerun xmeshfem2D in the adjoint run.
# However, in the current structure of the code, the xspecfem2D program can not detect the
# the parameter change in Par_file. Since in the adjoint run we change the SIMULATION_TYPE and the save_forward

./run_this_example.sh noclean > output.log

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
mv -v output.log OUTPUT_FILES/output.kernel.log

# backup
rm -rf OUTPUT_FILES.syn.adjoint
cp -rp OUTPUT_FILES OUTPUT_FILES.syn.adjoint

# kernels
cp -v OUTPUT_FILES/output.kernel.log KERNELS/
cp -v OUTPUT_FILES/*_kernel.* KERNELS/


########################### model update ################################

echo
echo "model update"
echo

# takes absolute value of percent
update_percent=$( sed "s/-//" <<< $perturb_percent )

echo "> ./model_update.py $NPROC $SIM_TYPE $update_percent"
echo
./model_update.py $NPROC $SIM_TYPE $update_percent

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# replaces model files with perturbed ones
for (( iproc = 0; iproc < $NPROC; iproc++ )); do
  rank=`printf "%06i\n" $iproc`
  cp -v DATA/proc${rank}_rho_new.bin DATA/proc${rank}_rho.bin
  cp -v DATA/proc${rank}_vp_new.bin DATA/proc${rank}_vp.bin
  cp -v DATA/proc${rank}_vs_new.bin DATA/proc${rank}_vs.bin
  if [[ $? -ne 0 ]]; then exit 1; fi
done

########################### final forward ################################

## forward simulation
echo
echo "running forward simulation (updated model)"
echo
./change_simulation_type.pl -f

# In principle we do not need rerun xmeshfem2D in the adjoint run.
# However, in the current structure of the code, the xspecfem2D program can not detect the
# the parameter change in Par_file. Since in the adjoint run we change the SIMULATION_TYPE and the save_forward

./run_this_example.sh > output.log
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
mv -v output.log OUTPUT_FILES/output.log

# backup
rm -rf OUTPUT_FILES.syn.updated
cp -rp OUTPUT_FILES OUTPUT_FILES.syn.updated


########################### kernel ################################

echo
echo "postprocessing..."
echo "> ./kernel_evaluation_postprocessing.py $NSTEP $DT $NPROC $SIM_TYPE"
echo
./kernel_evaluation_postprocessing.py $NSTEP $DT $NPROC $SIM_TYPE

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done: `date`"
echo
