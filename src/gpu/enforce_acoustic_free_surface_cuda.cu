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

/* KERNEL for enforce free surface */

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(acoustic_enforce_free_surf_cuda,
              ACOUSTIC_ENFORCE_FREE_SURF_CUDA)(long* Mesh_pointer,int* compute_wavefield_1,int* compute_wavefield_2) {

  TRACE("acoustic_enforce_free_surf_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // does not absorb free surface, thus we enforce the potential to be zero at surface

  // checks if anything to do
  if (mp->num_free_surface_faces == 0) return;

  // block sizes
  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->num_free_surface_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLLX,1,1);


  // sets potentials to zero at free surface
  if (*compute_wavefield_1) {
  enforce_free_surface_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_acoustic,
                                                                          mp->d_potential_dot_acoustic,
                                                                          mp->d_potential_dot_dot_acoustic,
                                                                          mp->num_free_surface_faces,
                                                                          mp->d_free_surface_ispec,
                                                                          mp->d_free_surface_ijk,
                                                                          mp->d_ibool,
                                                                          mp->d_ispec_is_acoustic);
  }
  // for backward/reconstructed potentials
  if (*compute_wavefield_2) {
    enforce_free_surface_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_acoustic,
                                                                            mp->d_b_potential_dot_acoustic,
                                                                            mp->d_b_potential_dot_dot_acoustic,
                                                                            mp->num_free_surface_faces,
                                                                            mp->d_free_surface_ispec,
                                                                            mp->d_free_surface_ijk,
                                                                            mp->d_ibool,
                                                                            mp->d_ispec_is_acoustic);
  }

  GPU_ERROR_CHECKING ("enforce_free_surface_cuda");
}

