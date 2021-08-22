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

/* KERNEL for PML boundary acoustic */

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(pml_boundary_acoustic_cuda,
              PML_BOUNDARY_ACOUSTIC_CUDA)(long* Mesh_pointer,int* compute_wavefield_1,int* compute_wavefield_2) {

  TRACE("pml_boundary_acoustic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int size = mp->pml_nglob_abs_acoustic;

  // checks if anything to do
  if (size == 0) return;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // sets potentials to zero at free surface
  if (*compute_wavefield_1) {
    pml_boundary_acoustic_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_acoustic,
                                                                             mp->d_potential_dot_acoustic,
                                                                             mp->d_potential_dot_dot_acoustic,
                                                                             mp->pml_nglob_abs_acoustic,
                                                                             mp->d_pml_abs_points_acoustic);
  }

  // for backward/reconstructed potentials
  if (*compute_wavefield_2) {
    pml_boundary_acoustic_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_acoustic,
                                                                             mp->d_b_potential_dot_acoustic,
                                                                             mp->d_b_potential_dot_dot_acoustic,
                                                                             mp->pml_nglob_abs_acoustic,
                                                                             mp->d_pml_abs_points_acoustic);
  }

  GPU_ERROR_CHECKING ("pml_boundary_acoustic_cuda");
}

