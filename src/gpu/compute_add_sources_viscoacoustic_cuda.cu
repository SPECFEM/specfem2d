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

// acoustic sources

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_ac_cuda,
              COMPUTE_ADD_SOURCES_AC_CUDA)(long* Mesh_pointer,
                                           int* iphasef,
                                           int * itf) {

  TRACE("compute_add_sources_ac_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // check if anything to do
  if (mp->nsources_local == 0) return;

  int iphase = *iphasef;

  // only adds this contribution for first pass
  if (iphase != 1) return;

  int it = *itf - 1;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nsources_local,&num_blocks_x,&num_blocks_y);

  // One block per source
  dim3 grid(num_blocks_x,num_blocks_y);
  // One thread per GLL
  dim3 threads(NGLLX,NGLLX,1);

  compute_add_sources_acoustic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                              mp->d_ibool,
                                                                              mp->d_sourcearrays,
                                                                              mp->d_source_time_function,
                                                                              mp->myrank,
                                                                              mp->d_ispec_selected_source,
                                                                              mp->d_ispec_is_acoustic,
                                                                              mp->d_kappastore,
                                                                              it,mp->nsources_local);

  // print_CUDA_error_if_any(cudaStreamSynchronize(mp->compute_stream),37);

  GPU_ERROR_CHECKING ("compute_add_sources_ac_cuda");
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_ac_s3_cuda,
              COMPUTE_ADD_SOURCES_AC_s3_CUDA)(long* Mesh_pointer,
                                              int* iphasef,
                                              int* itf) {
  // Same function than previous one but applying to mp->d_b_potential_dot_dot_acoustic instead
  // of mp->d_potential_dot_dot_acoustic

  TRACE("compute_add_sources_ac_s3_cuda");

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

  compute_add_sources_acoustic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_dot_dot_acoustic,
                                                                              mp->d_ibool,
                                                                              mp->d_sourcearrays,
                                                                              mp->d_source_time_function,
                                                                              mp->myrank,
                                                                              mp->d_ispec_selected_source,
                                                                              mp->d_ispec_is_acoustic,
                                                                              mp->d_kappastore,
                                                                              it,mp->nsources_local);


    //  print_CUDA_error_if_any(cudaStreamSynchronize(mp->compute_stream),38);

  GPU_ERROR_CHECKING ("compute_add_sources_ac_s3_cuda");
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_moving_sources_ac_cuda,
              COMPUTE_ADD_MOVING_SOURCES_AC_CUDA)(long* Mesh_pointer,
                                                  int* iphase_f,
                                                  int* nsources_local_moving,
                                                  int* itf,
                                                  int* NSTEP_f,
                                                  int* nsources_f) {

  // Same function but for moving sources

  TRACE("compute_add_moving_sources_ac_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int nsources_local = *nsources_local_moving;
  // Beware! nsources_local_moving is a pointer to an integer, not to the array of size NSTEP

  // check if anything to do
  if (nsources_local == 0) return;

  int iphase = *iphase_f;

  // only adds this contribution for first pass
  if (iphase != 1) return;

  int it = *itf - 1;
  int NSTEP = *NSTEP_f;
  int nsources = *nsources_f;

  // Look up for the best way to distribute the sources in a grid (one block per source)
  int num_blocks_x, num_blocks_y;
  get_blocks_xy(nsources_local, &num_blocks_x, &num_blocks_y);

   // create the grid (one block per source)
  dim3 grid(num_blocks_x,num_blocks_y);
  // create the threads in the block (one thread per GLL), it has NGLLX,NGLLX,1 structure
  dim3 threads(NGLLX,NGLLX,1);

  // note: we will launch kernels only for local sources since sourcearrays_moving(..) and ispec_selected_souce_movie(..)
  //       are sorted along local sources only.
  //       this avoids the need to have an array like islice_selected_source(NSOURCES,NSTEP) on the GPU.

  // Useful to check memory
  // cudaMemoryTest(3);

  compute_add_moving_sources_acoustic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                                    mp->d_ibool,
                                                                                    mp->d_sourcearrays_moving,
                                                                                    mp->d_source_time_function_moving,
                                                                                    mp->myrank,
                                                                                    mp->d_ispec_selected_source_moving,
                                                                                    mp->d_ispec_is_acoustic,
                                                                                    mp->d_kappastore,
                                                                                    it,
                                                                                    nsources_local,
                                                                                    nsources,
                                                                                    NSTEP);
  // cudaMemoryTest(334); // Useful to check memory

  GPU_ERROR_CHECKING ("compute_add_moving_sources_ac_cuda");
}


/* ----------------------------------------------------------------------------------------------- */

// acoustic adjoint sources

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(add_sources_ac_sim_2_or_3_cuda,
              ADD_SOURCES_AC_SIM_2_OR_3_CUDA)(long* Mesh_pointer,
                                               int* iphasef,
                                               int* itf,
                                               int* nadj_rec_local,
                                               int* NSTEP) {

  TRACE("add_sources_ac_sim_2_or_3_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int iphase = *iphasef;

  // only adds this contributions for first pass
  if (iphase != 1) return;

  // checks
  if (*nadj_rec_local != mp->nadj_rec_local) exit_on_cuda_error("add_sources_ac_sim_type_2_or_3: nadj_rec_local not equal\n");

  // checks if anything to do
  if (mp->nadj_rec_local == 0) return;

  int it = *itf - 1; // C-arrays start at 0

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nadj_rec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLLX,NGLLX,1);

  // launches cuda kernel for acoustic adjoint sources
  add_sources_ac_SIM_TYPE_2_OR_3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                               mp->d_source_adjoint,
                                                                               mp->d_xir_store_loc,
                                                                               mp->d_gammar_store_loc,
                                                                               mp->d_ibool,
                                                                               mp->d_ispec_is_acoustic,
                                                                               mp->d_ispec_selected_rec_loc,
                                                                               it,
                                                                               mp->nadj_rec_local,
                                                                               //mp->d_kappastore,
                                                                               *NSTEP);

//  print_CUDA_error_if_any(cudaStreamSynchronize(mp->compute_stream),38);

  GPU_ERROR_CHECKING ("add_sources_acoustic_SIM_TYPE_2_OR_3_kernel");
}
