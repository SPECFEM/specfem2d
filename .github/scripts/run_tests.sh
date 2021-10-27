#!/bin/bash
#
# runs a test example case
#

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=`pwd`
dir=${TESTDIR}

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
  echo "testing seismograms:"
  ln -s $WORKDIR/utils/compare_seismogram_correlations.py
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/
  if [[ $? -ne 0 ]]; then exit 1; fi
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.999 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
  if [[ $? -ne 0 ]]; then exit 1; fi
}

# test example
cd $dir

# default setup
# limit number of time steps
#sed -i "s:^NSTEP .*:NSTEP = 200:" DATA/Par_file
# shortens output interval to avoid timeouts
#sed -i "s:^NTSTEP_BETWEEN_OUTPUT_INFO .*:NTSTEP_BETWEEN_OUTPUT_INFO    = 50:" DATA/Par_file

# specific example setups
if [ "${TESTDIR}" == "EXAMPLES/moving_sources_acoustic" ]; then
  sed -i "s:^NSTEP .*:NSTEP    = 3000:" DATA/Par_file
fi
if [ "${TESTDIR}" == "EXAMPLES/Industrial_Format_SEP" ]; then
  sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
fi
if [ "${TESTDIR}" == "EXAMPLES/axisymmetric_case_AXISYM_option" ]; then
  sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP    = 1400:" DATA/Par_file
fi


# debug
if [ "${DEBUG}" == "true" ]; then
  # limit for debugging
  sed -i "s:^NSTEP .*:NSTEP = 5:" DATA/Par_file
fi

# selects kernel script for kernel benchmark examples
do_kernel_script=0

if [ "${TESTDIR}" == "EXAMPLES/BENCHMARK_CLAERBOUT_ADJOINT/ACOUSTIC" ]; then do_kernel_script=1; fi
if [ "${TESTDIR}" == "EXAMPLES/BENCHMARK_CLAERBOUT_ADJOINT/ELASTIC" ]; then do_kernel_script=1; fi
if [ "${TESTDIR}" == "EXAMPLES/BENCHMARK_CLAERBOUT_ADJOINT/ACOUSTIC_ELASTIC" ]; then do_kernel_script=1; fi

# setup elastic kernel for SH simulations
if [ "${TESTDIR}" == "EXAMPLES/BENCHMARK_CLAERBOUT_ADJOINT/ELASTIC" ] && [ "${TESTCASE}" == "SH" ]; then
  # sets simulation type for SH-waves
  sed -i "s:^SIM_TYPE.*:SIM_TYPE=3:" run_this_example_kernel.sh
fi

# runs simulation
if [ "$do_kernel_script" == "1" ]; then
  # kernel test script
  ./run_this_example_kernel.sh
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
if [ "${DEBUG}" == "true" ] || [ "$do_kernel_script" == "1" ]; then
  # no comparisons
  :       # do nothing
else
  # tests seismograms
  my_test
fi

# cleanup
rm -rf OUTPUT_FILES*

echo
echo "all good"
echo `date`
echo
