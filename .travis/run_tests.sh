#!/bin/bash

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

# checks if anything to do
echo "run checks: $RUN_CHECKS"
if [ "$RUN_CHECKS" == "0" ]; then
  echo "  no run checks required, exiting..."
  exit 0
else
  echo "  run checks required, start testing..."
fi
echo

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
18) dir=EXAMPLES/initial_mode_LDDRK/ ;;
19) dir=EXAMPLES/moving_sources_acoustic/ ;;
20) dir=EXAMPLES/anisotropic_isotropic_model/ ;;
21) dir=EXAMPLES/infinite_homogeneous_moment_tensor_vertical_dip_slip/ ;;
*) dir=EXAMPLES/simple_topography_and_also_a_simple_fluid_layer/ ;;
esac

# info
#echo $TRAVIS_BUILD_DIR
echo $WORKDIR
echo `date`
echo
echo "**********************************************************"
echo
echo "run test: TESTID=${TESTID} TESTDIR=${TESTDIR} TESTCOV=${TESTCOV}"
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
  ./compare_seismogram_correlations.py REF_SEIS/ OUTPUT_FILES/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.999 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}'
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
# test examples
###########################################################

# testing internal mesher example (short & quick for all configuration)
echo 'Tests...' && echo -en 'travis_fold:start:tests\\r'

# runs test
echo "test directory: $dir"
echo
cd $dir

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
    #sed -i "s:^NSTEP .*:NSTEP    = 3000:" DATA/Par_file
  fi
  # Rayleigh wave
  if [ "$TESTDIR" == "10" ]; then
    sed -i "s:^NPROC .*:NPROC    = 4:" DATA/Par_file
  fi
  if [ "$TESTDIR" == "11" ]; then
    sed -i "s:^NPROC .*:NPROC    = 4:" DATA/Par_file
  fi
  # fluid-solid
  if [ "$TESTDIR" == "16" ]; then
    sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
  fi
  # porous
  if [ "$TESTDIR" == "17" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 2000:" DATA/Par_file
  fi
  # moving sources
  if [ "$TESTDIR" == "19" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 3000:" DATA/Par_file
  fi

  # elastic kernel example Tromp2005_kernel/ w/ NO_BACKWARD_RECONSTRUCTION
  if [ "$TESTID" == "26" ]; then
    sed -i "s:^NO_BACKWARD_RECONSTRUCTION .*:NO_BACKWARD_RECONSTRUCTION = .true.:" DATA/Par_file
    sed -i "s:^NTSTEP_BETWEEN_COMPUTE_KERNELS .*:NTSTEP_BETWEEN_COMPUTE_KERNELS = 12:" DATA/Par_file
  fi

  # coverage run
  if [ "$TESTCOV" == "1" ]; then
    sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  fi

  # default run script
  if [ "$TESTDIR" == "4" ] || [ "$TESTID" == "24" ]; then
    # kernel script
    ./run_this_example_kernel.sh
    if [[ $? -ne 0 ]]; then exit 1; fi

    # simulation done
    echo
    echo "simulation done: `pwd`"
    echo `date`
    echo

    # kernel comparison
    if [ "$TESTDIR" == "4" ]; then
      my_test_kernel
    fi
  else
    # default script
    ./run_this_example.sh
    if [[ $? -ne 0 ]]; then exit 1; fi

    # simulation done
    echo
    echo "simulation done: `pwd`"
    echo `date`
    echo

    # seismogram comparison
    my_test
  fi
fi
if [[ $? -ne 0 ]]; then exit 1; fi

# simulation done
echo
echo "test done: `pwd`"
echo `date`
echo

echo -en 'travis_fold:end:tests\\r'
echo

# code coverage: https://codecov.io/gh/geodynamics/specfem2d/
# additional runs for coverage
#
# note: log becomes too long, trying to fold each test output
cd $WORKDIR

# first coverage tester (without mpi)
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.noise\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing noise example
  ##
  echo "##################################################################"
  echo "EXAMPLES/noise_uniform/"
  echo
  cd EXAMPLES/noise_uniform/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.noise\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.tape\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing Tape2007 example
  ##
  echo "##################################################################"
  echo "EXAMPLES/Tape2007/"
  echo
  cd EXAMPLES/Tape2007/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.tape\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.moment_tensor\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing infinite_homogeneous_moment_tensor_vertical_dip_slip example
  ##
  echo "##################################################################"
  echo "EXAMPLES/infinite_homogeneous_moment_tensor_vertical_dip_slip/"
  echo
  cd EXAMPLES/infinite_homogeneous_moment_tensor_vertical_dip_slip/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.moment_tensor\\r'

# second coverage tester (with mpi)
echo 'Coverage...' && echo -en 'travis_fold:start:coverage.semi_infinite\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing example with pml (longer testing only together with mpi and code coverage)
  ##
  echo "##################################################################"
  echo "EXAMPLES/semi_infinite_homogeneous/"
  echo
  cd EXAMPLES/semi_infinite_homogeneous/
  sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.semi_infinite\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.gmsh\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing external mesher example with mpi and stacey
  ##
  echo "##################################################################"
  echo "EXAMPLES/Gmsh_example_Stacey_MPI/"
  echo
  cd EXAMPLES/Gmsh_example_Stacey_MPI/
  sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
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
  echo "##################################################################"
  echo "EXAMPLES/Tromp2005_kernel/"
  echo
  cd EXAMPLES/Tromp2005_kernel/
  sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  ./run_this_example_kernel.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
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
  echo "##################################################################"
  echo "EXAMPLES/poroelastic_acoustic/"
  echo
  cd EXAMPLES/poroelastic_acoustic/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.poroelastic_acoustic\\r'


echo 'Coverage...' && echo -en 'travis_fold:start:coverage.rayleigh\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing plane wave example
  ##
  echo "##################################################################"
  echo "EXAMPLES/Rayleigh_wave_no_crack/"
  echo
  cd EXAMPLES/Rayleigh_wave_no_crack/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  # coarser resolution
  sed -i "s:60:15:g" DATA/Par_file
  sed -i "s:28:7:g" DATA/Par_file
  sed -i "s:28:7:g" DATA/interfaces_Rayleigh_flat.dat
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.rayleigh\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.fluid_solid\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing fluid solid w/ external mesh
  ##
  echo "##################################################################"
  echo "EXAMPLES/fluid_solid/fluid_solid_external_mesh/"
  echo
  cd EXAMPLES/fluid_solid/fluid_solid_external_mesh/
  sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.fluid_solid\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.initial_mode_LDDRK\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing fluid solid w/ external mesh
  ##
  echo "##################################################################"
  echo "EXAMPLES/initial_mode_LDDRK"
  echo
  cd EXAMPLES/initial_mode_LDDRK
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.initial_mode_LDDRK\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.no_backward\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## elastic kernel example Tromp2005_kernel/ w/ NO_BACKWARD_RECONSTRUCTION
  ##
  echo "##################################################################"
  echo "EXAMPLES/Tromp2005_kernel/"
  echo
  cd EXAMPLES/Tromp2005_kernel/
  sed -i "s:^NSTEP .*:NSTEP    = 500:" DATA/Par_file
  sed -i "s:^NO_BACKWARD_RECONSTRUCTION .*:NO_BACKWARD_RECONSTRUCTION = .true.:" DATA/Par_file
  sed -i "s:^NTSTEP_BETWEEN_COMPUTE_KERNELS .*:NTSTEP_BETWEEN_COMPUTE_KERNELS = 12:" DATA/Par_file
  ./run_this_example_kernel.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # no kernel value testing: only execution failure
  #my_test_kernel
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.no_backward\\r'



#################################################################
##
## tested by github actions
##
#################################################################

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.moving_sources\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## moving sources
  ##
  echo "##################################################################"
  echo "EXAMPLES/moving_sources_acoustic/"
  echo
  cd EXAMPLES/moving_sources_acoustic/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.moving_sources\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.anisotropic_isotropic_model\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## anisotropy
  ##
  echo "##################################################################"
  echo "EXAMPLES/anisotropic_isotropic_model/"
  echo
  cd EXAMPLES/anisotropic_isotropic_model/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.anisotropic_isotropic_model\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.industrial\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing MPI, SEP example
  ##
  echo "##################################################################"
  echo "EXAMPLES/Industrial_Format_SEP/"
  echo
  cd EXAMPLES/Industrial_Format_SEP/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.industrial\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.axisym\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "1" ]; then
  ##
  ## testing axisymmetric example
  ##
  echo "##################################################################"
  echo "EXAMPLES/axisymmetric_case_AXISYM_option/"
  echo
  cd EXAMPLES/axisymmetric_case_AXISYM_option/
  sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.axisym\\r'

## w/out MPI compilation

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.simple_topo\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing PML & MPI example
  ##
  echo "##################################################################"
  echo "EXAMPLES/simple_topography_and_also_a_simple_fluid_layer/"
  echo
  cd EXAMPLES/simple_topography_and_also_a_simple_fluid_layer/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.simple_topo\\r'

echo 'Coverage...' && echo -en 'travis_fold:start:coverage.tomo_ocean\\r'
if [ "$TESTCOV" == "1" ] && [ "$TESTID" == "0" ]; then
  ##
  ## testing tomographic_ocean_model example
  ##
  echo "##################################################################"
  echo "EXAMPLES/tomographic_ocean_model/"
  echo
  cd EXAMPLES/tomographic_ocean_model/
  sed -i "s:^NSTEP .*:NSTEP    = 10:" DATA/Par_file
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi
  # only for coverage, comparison would fail: my_test
  cd $WORKDIR
fi
echo -en 'travis_fold:end:coverage.tomo_ocean\\r'


# done
echo "all done"
echo `date`
echo

