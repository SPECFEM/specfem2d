#!/bin/bash

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

###########################################################
# setup
###########################################################
# chooses example directory
case "$TESTDIR" in
0) dir=./ ;;
1) dir=EXAMPLES/simple_topography_and_also_a_simple_fluid_layer/ ;;
2) dir=EXAMPLES/semi_infinite_homogeneous/ ;;
3) dir=EXAMPLES/Gmsh_example_Stacey_MPI/ ;;
4) dir=EXAMPLES/Tromp2005_kernel/ ;;
5) dir=EXAMPLES/poroelastic_acoustic/ ;;
6) dir=EXAMPLES/noise_uniform/ ;;
7) dir=EXAMPLES/axisymmetric_case_AXISYM_option/ ;;
8) dir=EXAMPLES/simple_topography_and_also_a_simple_fluid_layer/ ;;
9) dir=EXAMPLES/Industrial_Format_SEP/ ;;
10) dir=EXAMPLES/Rayleigh_wave_no_crack/ ;;
11) dir=EXAMPLES/Rayleigh_wave_with_crack/ ;;
12) dir=EXAMPLES/Tape2007/ ;;
13) dir=EXAMPLES/check_absolute_amplitude_of_pressure_source_seismograms_acoustic/ ;;
14) dir=EXAMPLES/check_absolute_amplitude_of_force_source_seismograms_elastic/ ;;
15) dir=EXAMPLES/check_absolute_amplitude_of_force_source_seismograms_viscoelastic/ ;;
16) dir=EXAMPLES/fluid_solid/fluid_solid_external_mesh/ ;;
17) dir=EXAMPLES/poroelastic_semi_infinite_homogeneous/ ;;
*) dir=EXAMPLES/simple_topography_and_also_a_simple_fluid_layer/ ;;
esac


# info
echo $TRAVIS_BUILD_DIR
echo $WORKDIR
echo
echo "**********************************************************"
echo
echo "configuration test: TESTID=${TESTID} TESTDIR=${TESTDIR} TESTCOV=${TESTCOV} TESTFLAGS=${TESTFLAGS}"
echo
echo "    test directory: $dir"
echo
echo "**********************************************************"
echo

# bash function for checking seismogram output with reference solutions
my_test(){
  echo "testing seismograms:"
  ln -s $WORKDIR/utils/compare_seismogram_correlations.py
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/
  if [[ $? -ne 0 ]]; then exit 1; fi
  # correlations
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.9 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
  if [[ $? -ne 0 ]]; then exit 1; fi
  # L2 errors
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/ | grep min/max | cut -d \| -f 4 | awk '{print "L2 error:",$1; if ($1 > 0.01 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
  if [[ $? -ne 0 ]]; then exit 1; fi
  # cleanup
  rm -rf OUTPUT_FILES/
}

# bash function for checking seismogram output with reference solutions
my_report(){
  # report example
  if [ -f OUTPUT_FILES/output_mesher.txt ]; then cat OUTPUT_FILES/output_mesher.txt; else exit 1; fi
  if [ -f OUTPUT_FILES/output_solver.txt ]; then cat OUTPUT_FILES/output_solver.txt; else exit 1; fi
}

# ASCII kernel files comparison
my_test_kernel(){
  echo "testing kernels:"
  ln -s $WORKDIR/utils/compare_kernel_correlations.py
  ./compare_kernel_correlations.py OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat REF_KERNEL/proc000000_rhop_alpha_beta_kernel.dat
  if [[ $? -ne 0 ]]; then exit 1; fi
  ./compare_kernel_correlations.py OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat REF_KERNEL/proc000000_rhop_alpha_beta_kernel.dat | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.9 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
  if [[ $? -ne 0 ]]; then exit 1; fi
  rm -rf OUTPUT_FILES/
}


###########################################################
# configuration & compilation
###########################################################
# configuration
echo 'Configure...' && echo -en 'travis_fold:start:configure\\r'
echo "configuration:"

if [ "$TESTCOV" == "1" ]; then
  echo "configuration: for coverage"
  ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TESTFLAGS} FLAGS_CHECK="-fprofile-arcs -ftest-coverage -O0" CFLAGS="-coverage -O0"
else
  if [ "$CUDA" == "true" ]; then
    if [ "$OPENCL" == "true" ]; then
      echo "configuration: for opencl" # uses libOpenCL provided from CUDA package
      ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TESTFLAGS} OCL_CPU_FLAGS="-g -Wall -std=c99 -DWITH_MPI" OCL_GPU_FLAGS="-Werror" OCL_INC="${CUDA_HOME}/include" OCL_LIB="${CUDA_HOME}/lib64" OCL_LIBS="-lOpenCL"
    else
      echo "configuration: for cuda"
      ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TESTFLAGS} CUDA_LIB="${CUDA_HOME}/lib64" CUDA_INC="${CUDA_HOME}/include" CUDA_FLAGS="-Xcompiler -Wall,-Wno-unused-function,-Wno-unused-const-variable,-Wfatal-errors -g -G"
    fi
  else
    echo "configuration: default"
    ./configure FC=${FC} MPIFC=${MPIFC} CC=${CC} ${TESTFLAGS}
  fi
fi

# (default already) we output to console
#sed -i "s:IMAIN .*:IMAIN = ISTANDARD_OUTPUT:" setup/constants.h

# Rayleigh wave and manual crack
if [ "$TESTID" == "17" ]; then
  sed -i "s:FAST_NUMBERING .*:FAST_NUMBERING = .false.:" setup/constants.h
  sed -i "s:ADD_A_SMALL_CRACK_IN_THE_MEDIUM .*:ADD_A_SMALL_CRACK_IN_THE_MEDIUM = .true.:" setup/constants.h
fi

echo -en 'travis_fold:end:configure\\r'

# compilation  (only cleaning)
echo 'Build...' && echo -en 'travis_fold:start:build\\r'
echo "compilation:"
make clean; make -j2 all
echo -en 'travis_fold:end:build\\r'


###########################################################
# test examples
###########################################################
# testing internal mesher example (short & quick for all configuration)
# runs test
echo "test directory: $dir"
echo
cd $dir
echo 'Tests...' && echo -en 'travis_fold:start:tests\\r'
if [ "$TESTID" == "2" ]; then
  # runs default tests
  make tests
else
  # specific example settings needed for testing
  if [ "$TESTDIR" == "5" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 1000:" DATA/Par_file
  fi
  if [ "$TESTDIR" == "7" ]; then
    sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
    sed -i "s:^NSTEP .*:NSTEP    = 1400:" DATA/Par_file
  fi
  if [ "$TESTDIR" == "9" ]; then
    sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
    sed -i "s:^NSTEP .*:NSTEP    = 3000:" DATA/Par_file
  fi
  if [ "$TESTDIR" == "17" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 2000:" DATA/Par_file
  fi

  # coverage run
  if [ "$TESTCOV" == "1" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 400:" DATA/Par_file
  fi

  # default
  if [ "$TESTDIR" == "4" ]; then
    # kernel script
    ./run_this_example_kernel.sh
    my_test_kernel
  else
    # default script
    ./run_this_example.sh
    my_test
  fi
  cd $WORKDIR
fi
echo -en 'travis_fold:end:tests\\r'


# code coverage: https://codecov.io/gh/geodynamics/specfem2d/
# additional runs for coverage
#
# note: log becomes too long, trying to fold each test output
# first coverage tester (without mpi)
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.noise\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing noise example
  ##
  cd EXAMPLES/noise_uniform/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file_noise_1
  ./run_this_example.sh
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.noise\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.tape\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing Tape2007 example
  ##
  cd EXAMPLES/Tape2007/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.tape\\r'

# second coverage tester (with mpi)
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.semi_infinite\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing example with pml (longer testing only together with mpi and code coverage)
  ##
  cd EXAMPLES/semi_infinite_homogeneous/
  sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  ./run_this_example.sh
  my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.semi_infinite\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.gmsh\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing external mesher example with mpi and stacey
  ##
  cd EXAMPLES/Gmsh_example_Stacey_MPI/
  sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  ./run_this_example.sh
  # tests mesh quality output
  awk '{if(NR==1){dy=sqrt(($2-13.3242693)^2);if(dy>1.e-5){print $0,"failed",dy;exit 1;}else{print $0,"good",dy;exit 0;}}}' OUTPUT_FILES/mesh_quality_histogram.txt
  if [[ $? -ne 0 ]]; then exit 1; fi
  my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.gmsh\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.tromp\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing kernel example
  ##
  cd EXAMPLES/Tromp2005_kernel/
  sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  ./run_this_example_kernel.sh
  # no kernel value testing: only execution failure
  #my_test_kernel
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.tromp\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.poroelastic_acoustic\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing poroelastic example
  ##
  cd EXAMPLES/poroelastic_acoustic/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.poroelastic_acoustic\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.axisym\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing axisymmetric example
  ##
  cd EXAMPLES/axisymmetric_case_AXISYM_option/
  sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.axisym\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.simple_topo\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing PML & MPI example
  ##
  cd EXAMPLES/simple_topography_and_also_a_simple_fluid_layer/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.simple_topo\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.industrial\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing MPI, SEP example
  ##
  cd EXAMPLES/Industrial_Format_SEP/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.industrial\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.rayleigh\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing plane wave example
  ##
  cd EXAMPLES/Rayleigh_wave_no_crack/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  # coarser resolution
  sed -i "s:60:15:g" DATA/Par_file
  sed -i "s:28:7:g" DATA/Par_file
  sed -i "s:28:7:g" DATA/interfaces_Rayleigh_flat.dat
  ./run_this_example.sh
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.rayleigh\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.fluid_solid\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing fluid solid w/ external mesh
  ##
  cd EXAMPLES/fluid_solid/fluid_solid_external_mesh/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.fluid_solid\\r'


# done
echo "done `pwd`"

