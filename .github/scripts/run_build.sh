#!/bin/bash
#
# builds all executables
#

# getting updated environment (CUDA_HOME, PATH, ..)
if [ -f $HOME/.tmprc ]; then source $HOME/.tmprc; fi

WORKDIR=`pwd`
TESTCOV=${TESTCOV:-}

# info
echo "work directory: $WORKDIR"
echo `date`
echo
echo "**********************************************************"
echo
echo "configuration test: TESTFLAGS=${TESTFLAGS} TESTNGLL=${TESTNGLL} TESTCOV=${TESTCOV}"
echo "                    CUDA=${CUDA}"
echo
echo "**********************************************************"
echo

# compiler infos
echo "compiler versions:"
echo "gcc --version"
gcc --version
echo "gfortran --version"
gfortran --version
echo "mpif90 --version"
mpif90 --version
echo

## CUDA
if [ "${CUDA}" == "true" ]; then
  echo
  echo "enabling CUDA"
  echo
  cuda=(--with-cuda=cuda13 CUDA_LIB="${CUDA_HOME}/lib64" CUDA_INC="${CUDA_HOME}/include" \
        CUDA_FLAGS="-Xcompiler -Wall,-Wno-unused-function,-Wno-unused-const-variable,-Wfatal-errors -g -G")
else
  cuda=()
fi

## special testflags
if [ "${TESTFLAGS}" == "check-mcmodel-medium" ]; then
  # note: this is a work-around as using the 'env:' parameter in the workflow 'CI.yml' with TESTFLAGS: FLAGS_CHECK=".."
  #       won't work as the FLAGS_CHECK string will then get split up and ./configure .. complains about unknown parameters.
  #       here, we re-define TESTFLAGS with a single quote around FLAGS_CHECK=".." to avoid the splitting.
  # use FLAGS_CHECK
  flags=(FLAGS_CHECK="-O3 -mcmodel=medium -std=f2008 -Wall -Wno-do-subscript -Wno-conversion -Wno-maybe-uninitialized")
  TESTFLAGS=""  # reset
else
  flags=()
fi

# configuration
echo
echo "configuration:"
echo

# split TESTFLAGS into individual items
set -- ${TESTFLAGS}

###########################################################
# configuration & compilation
###########################################################
# configuration

if [ "${TESTCOV}" == "true" ]; then
  echo "configuration: for coverage"
  ./configure \
    "${adios[@]}" \
    "${hdf[@]}" \
    "${cuda[@]}" \
    "${hip[@]}" \
    "${flags[@]}" \
    FLAGS_CHECK="-fprofile-arcs -ftest-coverage -O0" CFLAGS="-coverage -O0" \
    FC=${FC} MPIFC=${MPIFC} CC=${CC} "$@"
else
  if [ "${CUDA}" == "true" ]; then
    echo "configuration: for cuda"
  else
    echo "configuration: default"
  fi
  ./configure \
    "${cuda[@]}" \
    "${flags[@]}" \
    FC=${FC} MPIFC=${MPIFC} CC=${CC} "$@"
fi

# checks
if [[ $? -ne 0 ]]; then echo "configuration failed:"; cat config.log; echo ""; echo "exiting..."; exit 1; fi

# w/ NGLL = 6
if [ "$TESTNGLL" == "6" ]; then
  sed -i "s:NGLLX =.*:NGLLX = 6:" setup/constants.h
fi

# we output to console
sed -i "s:IMAIN .*:IMAIN = ISTANDARD_OUTPUT:" setup/constants.h

# compilation
echo
echo "compilation:"
make clean; make -j2 all

# checks
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done "
echo `date`
echo
