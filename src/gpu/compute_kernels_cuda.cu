/*
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================
*/

#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// ELASTIC SIMULATIONS

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_kernels_elastic_cuda,
              COMPUTE_KERNELS_ELASTIC_CUDA)(long* Mesh_pointer,realw * deltat) {

  TRACE("compute_kernels_elastic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL2; // NGLLX*NGLLZ

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_AB,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  compute_kernels_cudakernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ispec_is_elastic,
                                                                    mp->p_sv,
                                                                    mp->d_ibool,
                                                                    mp->d_displ,
                                                                    mp->d_accel,
                                                                    mp->d_b_displ,
                                                                    mp->d_rho_kl,
                                                                    mp->d_mu_kl,
                                                                    mp->d_kappa_kl,
                                                                    mp->NSPEC_AB,
                                                                    mp->d_hprime_xx,
                                                                    mp->d_xix,mp->d_xiz,
                                                                    mp->d_gammax,mp->d_gammaz);

  GPU_ERROR_CHECKING ("compute_kernels_elastic_cuda");
}


/* ----------------------------------------------------------------------------------------------- */

// ACOUSTIC SIMULATIONS

/* ----------------------------------------------------------------------------------------------- */



/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_kernels_acoustic_cuda,
              COMPUTE_KERNELS_ACOUSTIC_CUDA)(long* Mesh_pointer,realw * deltat) {

TRACE("compute_kernels_acoustic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL2; // NGLLX*NGLLZ

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_AB,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  compute_kernels_acoustic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ispec_is_acoustic,
                                                                         mp->d_ibool,
                                                                         mp->d_rhostore,
                                                                         mp->d_kappastore,
                                                                         mp->d_hprime_xx,
                                                                         mp->d_xix,mp->d_xiz,
                                                                         mp->d_gammax,mp->d_gammaz,
                                                                         mp->d_potential_acoustic,
                                                                         mp->d_potential_dot_dot_acoustic,
                                                                         mp->d_b_potential_acoustic,
                                                                         mp->d_rho_ac_kl,
                                                                         mp->d_kappa_ac_kl,
                                                                         mp->NSPEC_AB,
                                                                         *deltat);

  GPU_ERROR_CHECKING ("compute_kernels_acoustic_kernel");
}

/* ----------------------------------------------------------------------------------------------- */

// preconditioner (approximate Hessian kernel)

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_kernels_hess_cuda,
              COMPUTE_KERNELS_HESS_CUDA)(long* Mesh_pointer,
                                         int* ELASTIC_SIMULATION,
                                         int* ACOUSTIC_SIMULATION) {
  TRACE("compute_kernels_hess_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL2; // NGLLX*NGLLZ

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_AB,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if (*ELASTIC_SIMULATION) {
    compute_kernels_hess_el_cudakernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ispec_is_elastic,
                                                                              mp->d_ibool,
                                                                              mp->d_accel,
                                                                              mp->d_b_accel,
                                                                              mp->d_hess_el_kl,
                                                                              mp->NSPEC_AB);
  }

  if (*ACOUSTIC_SIMULATION) {
    compute_kernels_hess_ac_cudakernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ispec_is_acoustic,
                                                                              mp->d_ibool,
                                                                              mp->d_potential_dot_dot_acoustic,
                                                                              mp->d_b_potential_dot_dot_acoustic,
                                                                              mp->d_rhostore,
                                                                              mp->d_hprime_xx,
                                                                              mp->d_xix,mp->d_xiz,
                                                                              mp->d_gammax,mp->d_gammaz,
                                                                              mp->d_hess_ac_kl,
                                                                              mp->NSPEC_AB);
  }


  GPU_ERROR_CHECKING ("compute_kernels_hess_cuda");
}

