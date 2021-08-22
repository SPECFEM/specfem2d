#!/bin/bash

# fortran/openMPI compiler
sudo apt-get install -y gfortran libgomp1 openmpi-bin libopenmpi-dev
echo
# python script needs numpy
# for distribution precise
#sudo apt-get install -qq python-numpy #python-scipy
# trusty:
# (Sep2017 update: adding flag --user : see https://github.com/travis-ci/travis-ci/issues/8382)
#pip install --user --upgrade pip setuptools wheel
#pip install --user --only-binary=numpy numpy
# bionic:
sudo apt-get install -qq python-numpy           # gets numpy version 1.13.3
#pip install --user --upgrade pip setuptools wheel
#pip install --user --only-binary=numpy numpy   # requirements misses ppc64 architecture
#pip install --user numpy                       # installation takes quite long (~3-4min), gets numpy version 1.19.5
echo

# version info
echo "Python on path: $(which python)"
python --version
echo
echo "pip on path   : $(which pip)"
pip --version
echo
echo "numpy version : "
python -c "import numpy; print(numpy.__version__)"
echo

# CPU architecture
echo
echo "CPU architecture: ${TRAVIS_CPU_ARCH}"
echo

# version infos
echo "compiler versions:" ${FC} ${MPIFC} ${CC}
${FC} --version
${MPIFC} --version
${CC} --version
echo

# installs the CUDA toolkit
if [ "$CUDA" == "true" ]; then
  # Linux environment
  # https://docs.travis-ci.com/user/reference/linux
  ## distribution precise: from ubuntu 12.04
  #UBUNTU_VERSION=ubuntu1204
  ## distribution trusty: from ubuntu 14.04
  #UBUNTU_VERSION=ubuntu1404
  ## distribution xenial: from ubuntu 16.04
  #UBUNTU_VERSION=ubuntu1604
  ## distribution bionic: from ubuntu 18.04
  UBUNTU_VERSION=ubuntu1804
  ## distribution focal: from ubuntu 20.04
  #UBUNTU_VERSION=ubuntu2004

  # CUDA_VERSION - specifies CUDA toolkit version
  # http://developer.download.nvidia.com/compute/cuda/repos/
  ## trusty
  #CUDA_VERSION=6.5-14
  ## xenial
  #CUDA_VERSION=9.2.148-1
  ## bionic
  CUDA_VERSION=10.2.89-1

  # default architecture amd64
  CUDA_OS=x86_64
  CUDA_ARCH=amd64
  if [ "${TRAVIS_CPU_ARCH}" == "ppc64le" ]; then
    CUDA_OS=ppc64el
    CUDA_ARCH=ppc64el # IBM Power
  elif [ "${TRAVIS_CPU_ARCH}" == "arm64" ]; then
    CUDA_OS=sbsa
    CUDA_ARCH=arm64 # ARM
  fi

  echo "Installing CUDA library"
  echo "CUDA version  : ${CUDA_VERSION}"
  echo "UBUNTU version: ${UBUNTU_VERSION}"
  echo "CUDA OS       : ${CUDA_OS}"
  echo "CUDA arch     : ${CUDA_ARCH}"

  # note: travis could stall and time out here
  #       one could try to add: travis_retry sudo dpgk -i ..
  #       https://docs.travis-ci.com/user/common-build-problems/#travis_retry
  #
  # remove old nvidia-cuda packages
  #sudo apt-get remove nvidia-cuda-* ;

  # package needs key
  sudo apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/${UBUNTU_VERSION}/${CUDA_OS}/7fa2af80.pub
  echo

  # gets packages
  INSTALLER=cuda-repo-${UBUNTU_VERSION}_${CUDA_VERSION}_${CUDA_ARCH}.deb
  wget http://developer.download.nvidia.com/compute/cuda/repos/${UBUNTU_VERSION}/${CUDA_OS}/${INSTALLER}
  sudo dpkg -i ${INSTALLER}
  echo

  # update
  echo "Updating libraries"
  sudo apt-get update -qq
  dpkg -l | grep cuda
  #export CUDA_APT=${CUDA_VERSION:0:3}  # version 9.2
  export CUDA_APT=${CUDA_VERSION:0:4}  # version 10.2
  export CUDA_APT=${CUDA_APT/./-}
  echo "CUDA: ${CUDA_APT}"  # apt version 10-2 -> package name: cuda-compiler-10-2

  # installs packages
  # CUDA_PACKAGES="cuda-drivers cuda-core-${CUDA_APT} cuda-cudart-dev-${CUDA_APT} cuda-cufft-dev-${CUDA_APT}";
  CUDA_PACKAGES="cuda-drivers cuda-compiler-${CUDA_APT} cuda-cudart-dev-${CUDA_APT}"
  echo "Installing ${CUDA_PACKAGES}"
  sudo apt-get install -y --no-install-recommends ${CUDA_PACKAGES}
  sudo apt-get clean
  #export CUDA_HOME=/usr/local/cuda-${CUDA_VERSION:0:3}   # version 9.2
  export CUDA_HOME=/usr/local/cuda-${CUDA_VERSION:0:4}    # version 10.2
  export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}
  export PATH=${CUDA_HOME}/bin:${PATH}
  echo ""
  nvcc --version
else
  export CUDA_HOME=""
fi

# storing updated environment parameters for following bash-script
echo "export PATH=${PATH}" > $HOME/.tmprc
echo "export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}" >> $HOME/.tmprc
echo "export CUDA_HOME=${CUDA_HOME}" >> $HOME/.tmprc

# to avoid mpi issues on travis
# see: https://github.com/open-mpi/ompi/issues/1393#issuecomment-187980473
#      https://github.com/open-mpi/ompi/issues/4948
echo "export export OMPI_MCA_btl_vader_single_copy_mechanism=none" >> $HOME/.tmprc
echo "export export OMPI_MCA_btl=^openib" >> $HOME/.tmprc

echo "exports:"
export
echo ""
