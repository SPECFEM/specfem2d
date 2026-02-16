#!/bin/bash
#
# runs additional coverage examples
#

set -euo pipefail

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=$(pwd)
TESTID=${TESTID:-}
TESTCOV=${TESTCOV:-}

if [ "$TESTCOV" != "true" ]; then
  echo "TESTCOV=${TESTCOV} (not coverage), skipping run_coverage.sh"
  exit 0
fi

run_simple() {
  local rel_dir="$1"
  local nstep="$2"
  local model="${3:-}"

  echo "##################################################################"
  echo "${rel_dir}"
  echo

  cd "${WORKDIR}/${rel_dir}"

  # setup
  cp -v DATA/Par_file DATA/Par_file.org
  sed -i "s:^NSTEP .*:NSTEP    = ${nstep}:" DATA/Par_file

  # parallel
  if [ "${rel_dir}" == "EXAMPLES/benchmarks/semi_infinite_homogeneous/" ]; then
    sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
  fi
  if [ "${rel_dir}" == "EXAMPLES/reproducible_study/Komatitsch2000_fluid_solid/fluid_solid_external_mesh/" ]; then
    sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
  fi
  if [ "${rel_dir}" == "EXAMPLES/applications/axisymmetric_examples/axisymmetric_case_AXISYM_option/" ]; then
    sed -i "s:^NPROC .*:NPROC    = 2:" DATA/Par_file
  fi

  # coarser mesh
  if [ "${rel_dir}" == "EXAMPLES/applications/Rayleigh_wave_no_crack/" ]; then
    # coarser resolution
    sed -i "s:60:15:g" DATA/Par_file
    sed -i "s:28:7:g" DATA/Par_file
    sed -i "s:28:7:g" DATA/interfaces_Rayleigh_flat.dat
  fi

  # run
  ./run_this_example.sh
  if [[ $? -ne 0 ]]; then exit 1; fi

  # mesher check
  if [ "${rel_dir}" == "EXAMPLES/applications/meshing/Gmsh_example_Stacey_MPI/" ]; then
    # tests mesh quality output
    awk '{if(NR==1){dy=sqrt(($2-13.3242693)^2);if(dy>1.e-5){print $0,"failed",dy;exit 1;}else{print $0,"good",dy;exit 0;}}}' OUTPUT_FILES/mesh_quality_histogram.txt
    if [[ $? -ne 0 ]]; then exit 1; fi
  fi
  
  # cleanup
  mv -v DATA/Par_file.org DATA/Par_file
  rm -rf OUTPUT_FILES/*
  if [ -e DATABASES_MPI ]; then rm -rf DATABASES_MPI/*; fi
  cd "$WORKDIR"
}

run_kernel() {
  local rel_dir="$1"
  local nstep="$2"
  local backward="$3"

  echo "##################################################################"
  echo "${rel_dir} (kernel coverage)"
  echo

  cd "${WORKDIR}/${rel_dir}"

  # setup
  cp -v DATA/Par_file DATA/Par_file.org
  sed -i "s:^NSTEP .*:NSTEP    = ${nstep}:" DATA/Par_file

  if [ "${backward}" == "no_backward" ]; then
    sed -i "s:^NO_BACKWARD_RECONSTRUCTION .*:NO_BACKWARD_RECONSTRUCTION = .true.:" DATA/Par_file
    sed -i "s:^NTSTEP_BETWEEN_COMPUTE_KERNELS .*:NTSTEP_BETWEEN_COMPUTE_KERNELS = 12:" DATA/Par_file
  fi

  # run
  ./run_this_example_kernel.sh
  if [[ $? -ne 0 ]]; then exit 1; fi

  # cleanup
  mv -v DATA/Par_file.org DATA/Par_file
  rm -rf OUTPUT_FILES/*
  if [ -e DATABASES_MPI ]; then rm -rf DATABASES_MPI/*; fi
  if [ -e SEM ]; then rm -rf SEM/*; fi
  cd "$WORKDIR"
}

echo
echo "coverage run: TESTID=${TESTID}"
echo "work directory: ${WORKDIR}"
echo

# additional example tests (after base to avoid repeating code setup/configuration/compilation)
case "$TESTID" in
  0) # serial bunch
    run_simple "EXAMPLES/applications/ocean_acoustics/tomographic_ocean_model/" 10
    run_simple "EXAMPLES/benchmarks/infinite_homogeneous_moment_tensor_vertical_dip_slip/" 10
    run_simple "EXAMPLES/reproducible_study/Tape2007/" 10
    run_simple "EXAMPLES/reproducible_study/Tromp2010_noise_uniform/" 10
    ;;
  1) # parallel bunch 1
    run_simple "EXAMPLES/applications/meshing/Gmsh_example_Stacey_MPI/" 10
    run_simple "EXAMPLES/applications/moving_sources_acoustic/" 10
    run_simple "EXAMPLES/applications/anisotropy/anisotropic_isotropic_model/" 10
    run_simple "EXAMPLES/benchmarks/semi_infinite_homogeneous/" 10
    run_simple "EXAMPLES/real_world/Marmousi2" 10
    run_simple "EXAMPLES/real_world/Industrial_Format_SEP/" 10
    run_kernel "EXAMPLES/reproducible_study/Tromp2005_kernel/" 500
    ;;
  2) # parallel bunch 2 - vectorization
    run_simple "EXAMPLES/applications/Rayleigh_wave_no_crack/" 10
    run_simple "EXAMPLES/applications/initial_mode_LDDRK" 10
    run_simple "EXAMPLES/applications/axisymmetric_examples/axisymmetric_case_AXISYM_option/" 10
    run_simple "EXAMPLES/reproducible_study/Morency2008_poroelastic_acoustic/" 10
    run_simple "EXAMPLES/reproducible_study/Komatitsch2000_fluid_solid/fluid_solid_external_mesh/" 10
    run_kernel "EXAMPLES/reproducible_study/Tromp2005_kernel" 500 "no_backward"
    ;;
  3) # NGLL 6
    echo "TESTID=3: no additional coverage examples (base run already executed)"
    ;;
  *)
    echo "TESTID=${TESTID}: no additional coverage examples configured"
    ;;
esac

echo
echo "coverage examples done"
echo "$(date)"
echo
