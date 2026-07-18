#!/bin/bash
#
# runs a test example case
#

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=`pwd`
dir=${TESTDIR}
TESTID=${TESTID:-}
TESTCOV=${TESTCOV:-}

# info
echo "work directory: $WORKDIR"
echo `date`
echo
echo "**********************************************************"
echo
echo "test directory: $dir"
echo
echo "**********************************************************"
echo

# bash function for checking seismogram output with reference solutions
my_test(){
  echo "######################################################################################################################"
  echo "testing seismograms:"
  ln -s $WORKDIR/utils/compare_seismogram_correlations.py
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/
  if [[ $? -ne 0 ]]; then exit 1; fi
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.999 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
  if [[ $? -ne 0 ]]; then exit 1; fi
  echo "######################################################################################################################"
}

# test example
cd $dir

# default setup
# limit number of time steps
#sed -i "s:^NSTEP .*:NSTEP = 200:" DATA/Par_file
# shortens output interval to avoid timeouts
#sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file

# specific example setups
if [ "${TESTDIR}" == "EXAMPLES/applications/moving_sources_acoustic" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 3000:" DATA/Par_file
fi
if [ "${TESTDIR}" == "EXAMPLES/real_world/Industrial_Format_SEP" ]; then
  sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
fi
if [ "${TESTDIR}" == "EXAMPLES/applications/axisymmetric_examples/axisymmetric_case_AXISYM_option" ]; then
  sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP    = 1400:" DATA/Par_file
fi
if [ "${TESTDIR}" == "EXAMPLES/reproducible_study/Komatitsch2000_fluid_solid/fluid_solid_external_mesh" ]; then
  sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
fi
if [ "${TESTDIR}" == "EXAMPLES/reproducible_study/Morency2008_poroelastic_semi_infinite_homogeneous" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 2000:" DATA/Par_file
fi
if [ "${TESTDIR}" == "EXAMPLES/applications/Rayleigh_wave_no_crack" ]; then
  sed -i "s:^NPROC .*:NPROC    = 4:" DATA/Par_file
fi
if [ "${TESTDIR}" == "EXAMPLES/applications/Rayleigh_wave_with_crack" ]; then
  sed -i "s:^NPROC .*:NPROC    = 4:" DATA/Par_file
fi


# debug
if [ "${DEBUG}" == "true" ]; then
  # limit for debugging
  sed -i "s:^NSTEP .*:NSTEP    = 5:" DATA/Par_file
fi

# coverage runs use short steps
if [ "$TESTCOV" == "true" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
fi

# selects kernel script for kernel benchmark examples
RUN_KERNEL=${RUN_KERNEL:-}
if [ "${TESTDIR}" == "EXAMPLES/benchmarks/BENCHMARK_CLAERBOUT_ADJOINT/ACOUSTIC" ]; then RUN_KERNEL=true; fi
if [ "${TESTDIR}" == "EXAMPLES/benchmarks/BENCHMARK_CLAERBOUT_ADJOINT/ELASTIC" ]; then RUN_KERNEL=true; fi
if [ "${TESTDIR}" == "EXAMPLES/benchmarks/BENCHMARK_CLAERBOUT_ADJOINT/ACOUSTIC_ELASTIC" ]; then RUN_KERNEL=true; fi

# setup elastic kernel for SH simulations
if [ "${TESTDIR}" == "EXAMPLES/benchmarks/BENCHMARK_CLAERBOUT_ADJOINT/ELASTIC" ] && [ "${TESTCASE}" == "SH" ]; then
  # sets simulation type for SH-waves
  sed -i "s:^SIM_TYPE.*:SIM_TYPE=3:" run_this_example_kernel.sh
fi

# save Par_file state
if [ -e DATA/Par_file ]; then
  cp -v DATA/Par_file DATA/Par_file.bak
fi

# runs simulation
if [ "${RUN_KERNEL}" == "true" ]; then
  # use kernel script
  ./run_this_example_kernel.sh | tee output.log
else
  # default script
  ./run_this_example.sh
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# simulation done
echo
echo "simulation done: `pwd`"
echo `date`
echo

# seismogram comparison
RUN_COMPARE=true
# turn off for non-default runs
if [ "${TESTCOV}" == "true" ]; then RUN_COMPARE=false; fi
if [ "${DEBUG}" == "true" ]; then RUN_COMPARE=false; fi
if [ "${RUN_KERNEL}" == "true" ]; then RUN_COMPARE=false; fi

if [ "${RUN_COMPARE}" == "true" ]; then
  my_test
else
  # no comparisons
  :     # do nothing
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# restore original Par_file
if [ -e DATA/Par_file.bak ]; then
  cp -v DATA/Par_file.bak DATA/Par_file
fi

# cleanup
rm -rf OUTPUT_FILES*
if [ -e DATABASES_MPI ]; then rm -rf DATABASES_MPI/; fi
if [ -e SEM ]; then rm -rf SEM/; fi

echo
echo "all good"
echo `date`
echo
