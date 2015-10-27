/*
  !========================================================================
  !
  !                   S P E C F E M 2 D  Version 7 . 0
  !                   --------------------------------
  !
  !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
  !                        Princeton University, USA
  !                and CNRS / University of Marseille, France
  !                 (there are currently many more authors!)
  ! (c) Princeton University and CNRS / University of Marseille, April 2014
  !
  ! This software is a computer program whose purpose is to solve
  ! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
  ! using a spectral-element method (SEM).
  !
  ! This software is governed by the CeCILL license under French law and
  ! abiding by the rules of distribution of free software. You can use,
  ! modify and/or redistribute the software under the terms of the CeCILL
  ! license as circulated by CEA, CNRS and Inria at the following URL
  ! "http://www.cecill.info".
  !
  ! As a counterpart to the access to the source code and rights to copy,
  ! modify and redistribute granted by the license, users are provided only
  ! with a limited warranty and the software's author, the holder of the
  ! economic rights, and the successive licensors have only limited
  ! liability.
  !
  ! In this respect, the user's attention is drawn to the risks associated
  ! with loading, using, modifying and/or developing or reproducing the
  ! software by the user in light of its specific status of free software,
  ! that may mean that it is complicated to manipulate, and that also
  ! therefore means that it is reserved for developers and experienced
  ! professionals having in-depth computer knowledge. Users are therefore
  ! encouraged to load and test the software's suitability as regards their
  ! requirements in conditions enabling the security of their systems and/or
  ! data to be ensured and, more generally, to use and operate it in the
  ! same conditions as regards security.
  !
  ! The full text of the license is available in file "LICENSE".
  !
  !========================================================================

 */

#include <stdio.h>
#include <stdlib.h>

#include "config.h"

typedef float realw;

void FC_FUNC_(initialize_cuda_device,
              INITIALIZE_CUDA_DEVICE)(int* myrank_f,int* ncuda_devices){}


void FC_FUNC_(prepare_gpu,
              PREPARE_GPU)(long * Container,
                          realw * xstore_me,
                          realw * zstore_me,
                          realw sigma_h2_inv,
                          realw sigma_v2_inv,
                          realw h_criterion,
                          realw v_criterion,
                          int nspec_me,
                          int nker,
                          realw gll){}

void FC_FUNC_(compute_smooth,
              COMPUTE_SMOOTH)(long * smooth_pointer,
                              realw * jacobian,
                              realw * xstore_other,
                              realw * zstore_other,
                              realw * data_other,
                              const int * nspec_other){}


void FC_FUNC_(get_smooth,
              GET_SMOOTH)(long * smooth_pointer,realw * data_smooth){}
