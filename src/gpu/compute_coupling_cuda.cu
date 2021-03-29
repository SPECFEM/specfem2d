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

// ACOUSTIC - ELASTIC coupling

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_coupling_ac_el_cuda,
              COMPUTE_COUPLING_AC_EL_CUDA)(long* Mesh_pointer,
                                           int* iphasef,
                                           int* num_coupling_ac_el_facesf) {
  TRACE("compute_coupling_ac_el_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  int iphase            = *iphasef;

  // only adds this contribution for first pass
  if (iphase != 1) return;

  int num_coupling_ac_el_faces  = *num_coupling_ac_el_facesf;

  // way 1: exact blocksize to match NGLLSQUARE
  int blocksize = NGLLX;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(num_coupling_ac_el_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // launches GPU kernel
  compute_coupling_acoustic_el_kernel<<<grid,threads>>>(mp->d_displ,
                                                       mp->d_potential_dot_dot_acoustic,
                                                       num_coupling_ac_el_faces,
                                                       mp->d_coupling_ac_el_ispec,
                                                       mp->d_coupling_ac_el_ijk,
                                                       mp->d_coupling_ac_el_normal,
                                                       mp->d_coupling_ac_el_jacobian2Dw,
                                                       mp->d_ibool);

  //  adjoint simulations
  if (mp->simulation_type == 3) {
    compute_coupling_acoustic_el_kernel<<<grid,threads>>>(mp->d_b_displ,
                                                          mp->d_b_potential_dot_dot_acoustic,
                                                          num_coupling_ac_el_faces,
                                                          mp->d_coupling_ac_el_ispec,
                                                          mp->d_coupling_ac_el_ijk,
                                                          mp->d_coupling_ac_el_normal,
                                                          mp->d_coupling_ac_el_jacobian2Dw,
                                                          mp->d_ibool);
  }


  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);

  GPU_ERROR_CHECKING ("compute_coupling_acoustic_el_kernel");
}


/* ----------------------------------------------------------------------------------------------- */

// ELASTIC - ACOUSTIC coupling

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_coupling_el_ac_cuda,
              COMPUTE_COUPLING_EL_AC_CUDA)(long* Mesh_pointer,
                                           int* iphasef,
                                           int* num_coupling_ac_el_facesf) {
  TRACE("compute_coupling_el_ac_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  int iphase            = *iphasef;

  // only adds this contribution for first pass
  if (iphase != 1) return;

  int num_coupling_ac_el_faces  = *num_coupling_ac_el_facesf;

  // way 1: exact blocksize to match NGLLX
  int blocksize = NGLLX;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(num_coupling_ac_el_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // launches GPU kernel
  compute_coupling_elastic_ac_kernel<<<grid,threads>>>(mp->d_potential_dot_dot_acoustic,
                                                       mp->d_accel,
                                                       num_coupling_ac_el_faces,
                                                       mp->d_coupling_ac_el_ispec,
                                                       mp->d_coupling_ac_el_ijk,
                                                       mp->d_coupling_ac_el_normal,
                                                       mp->d_coupling_ac_el_jacobian2Dw,
                                                       mp->d_ibool);

  //  adjoint simulations
  if (mp->simulation_type == 3) {
    compute_coupling_elastic_ac_kernel<<<grid,threads>>>(mp->d_b_potential_dot_dot_acoustic,
                                                         mp->d_b_accel,
                                                         num_coupling_ac_el_faces,
                                                         mp->d_coupling_ac_el_ispec,
                                                         mp->d_coupling_ac_el_ijk,
                                                         mp->d_coupling_ac_el_normal,
                                                         mp->d_coupling_ac_el_jacobian2Dw,
                                                         mp->d_ibool);
  }

  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);

  GPU_ERROR_CHECKING ("compute_coupling_el_ac_cuda");
}
