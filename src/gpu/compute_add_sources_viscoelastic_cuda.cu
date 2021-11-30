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

// elastic domain sources

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_el_cuda,
              COMPUTE_ADD_SOURCES_EL_CUDA)(long* Mesh_pointer,
                                           int* iphasef,
                                           int* itf) {

  TRACE("compute_add_sources_el_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // check if anything to do
  if (mp->nsources_local == 0) return;

  int iphase = *iphasef;

  // only adds this contribution for first pass
  if (iphase != 1) return;

  int it = *itf - 1;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nsources_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLX,1);

  compute_add_sources_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,mp->d_ibool,
                                                                    mp->d_sourcearrays,
                                                                    mp->d_source_time_function,
                                                                    mp->myrank,
                                                                    mp->d_ispec_selected_source,
                                                                    mp->d_ispec_is_elastic,
                                                                    it,mp->nsources_local);


  GPU_ERROR_CHECKING ("compute_add_sources_el_cuda");
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_el_s3_cuda,
              COMPUTE_ADD_SOURCES_EL_S3_CUDA)(long* Mesh_pointer,
                                              int* iphasef,
                                              int* itf) {

  TRACE("compute_add_sources_el_s3_cuda");
  // EPIK_TRACER("compute_add_sources_el_s3_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // check if anything to do
  if (mp->nsources_local == 0) return;

  int iphase = *iphasef;

  // only adds this contribution for first pass
  if (iphase != 1) return;

  int it = *itf - 1;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nsources_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLX,1);

  compute_add_sources_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel,mp->d_ibool,
                                                                    mp->d_sourcearrays,
                                                                    mp->d_source_time_function,
                                                                    mp->myrank,
                                                                    mp->d_ispec_selected_source,
                                                                    mp->d_ispec_is_elastic,
                                                                    it,mp->nsources_local);

  GPU_ERROR_CHECKING ("compute_add_sources_el_s3_cuda");
}

/* ----------------------------------------------------------------------------------------------- */

// ADJOINT sources

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(add_sources_el_sim_type_2_or_3,
              ADD_SOURCES_EL_SIM_TYPE_2_OR_3)(long* Mesh_pointer,
                                               int* iphasef,
                                               int* itf,
                                               int* nadj_rec_local,
                                               int* NSTEP) {

  TRACE("add_sources_el_sim_type_2_or_3");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks
  if (*nadj_rec_local != mp->nadj_rec_local) exit_on_error("add_sources_el_sim_type_2_or_3: nadj_rec_local not equal\n");

  // checks if anything to do
  if (mp->nadj_rec_local == 0) return;

  int iphase = *iphasef;

  // only add this contribution for first pass
  if (iphase != 1) return;

  int it = *itf - 1; // C-arrays start at 0

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nadj_rec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLLX,NGLLX,1);

  // the irec_local variable needs to be precomputed (as
  // h_pre_comp..), because normally it is in the loop updating accel,
  // and due to how it's incremented, it cannot be parallelized

  add_sources_el_SIM_TYPE_2_OR_3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                               mp->d_source_adjoint,
                                                                               mp->d_xir_store_loc,
                                                                               mp->d_gammar_store_loc,
                                                                               mp->d_ibool,
                                                                               mp->d_ispec_is_elastic,
                                                                               mp->d_ispec_selected_rec_loc,
                                                                               it,
                                                                               mp->nadj_rec_local,
                                                                               *NSTEP);

  GPU_ERROR_CHECKING ("add_sources_SIM_TYPE_2_OR_3_kernel");
}

