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

__global__ void compute_add_sources_acoustic_kernel(realw* potential_dot_dot_acoustic,
                                                    int* d_ibool,
                                                    realw* sourcearrays,
                                                    realw* source_time_function,
                                                    int myrank,
                                                    int* ispec_selected_source,
                                                    int* ispec_is_acoustic,
                                                    realw* kappastore,
                                                    int it,
                                                    int nsources_local) {

  int i = threadIdx.x;
  int j = threadIdx.y;

  int isource  = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int ispec,iglob;
  realw stf,kappal,accel;

  if (isource < nsources_local) {

      ispec = ispec_selected_source[isource] - 1;

      if (ispec_is_acoustic[ispec]) {
        iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

        kappal = kappastore[INDEX3(NGLLX,NGLLX,i,j,ispec)];
        stf = source_time_function[INDEX2(nsources_local,isource,it)]/kappal;
        accel = sourcearrays[INDEX4(NDIM,NGLLX,NGLLX, 0,i,j,isource)] * stf;

        atomicAdd(&potential_dot_dot_acoustic[iglob], accel);
    }
  }
}


/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_add_moving_sources_acoustic_kernel(realw* potential_dot_dot_acoustic,
                                                           int* d_ibool,
                                                           realw* sourcearrays_moving,
                                                           realw* source_time_function_moving,
                                                           int myrank,
                                                           int* ispec_selected_source_moving,
                                                           int* ispec_is_acoustic,
                                                           realw* kappastore,
                                                           int it,
                                                           int nsources_local,
                                                           int nsources,
                                                           int NSTEP) {

  // Same but for moving sources
  int i = threadIdx.x;  // correspond to GLLx id (each thread has its own x and y id)
  int j = threadIdx.y;  // correspond to GLLz id (each thread has its own x and y id)

  int isource = blockIdx.x + gridDim.x*blockIdx.y;

  int ispec, iglob;
  realw stf, kappal, accel;

  if (isource < nsources_local) {
    ispec = ispec_selected_source_moving[INDEX2(nsources, isource, it)] - 1;
    if (ispec_is_acoustic[ispec]) {
      iglob = d_ibool[INDEX3_PADDED(NGLLX, NGLLX, i, j, ispec)] - 1;

      kappal = kappastore[INDEX3(NGLLX, NGLLX, i, j, ispec)];
      stf = source_time_function_moving[INDEX2(nsources, isource, it)]/kappal;
      accel = sourcearrays_moving[INDEX5(NDIM, NGLLX, NGLLX, nsources, 0, i, j, isource, it)] * stf;

      atomicAdd(&potential_dot_dot_acoustic[iglob], accel);
    }
  }
}


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

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nsources_local,&num_blocks_x,&num_blocks_y);

  // One block per source
  dim3 grid(num_blocks_x,num_blocks_y);
  // One thread per GLL
  dim3 threads(NGLLX,NGLLX,1);

  int it = *itf - 1;

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


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_ac_cuda");
#endif
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

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nsources_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLX,1);

  int it = *itf - 1;

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

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_ac_s3_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_moving_sources_ac_cuda,
              COMPUTE_ADD_MOVING_SOURCES_AC_CUDA)(long* Mesh_pointer,
                                                  int* iphase_f,
                                                  int* nsources_local_moving,
                                                  int* it_f,
                                                  int* NSTEP_f,
                                                  int* nsources_f) {

  // Same function but for moving sources

  TRACE("compute_add_moving_sources_ac_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int nsources_local = *nsources_local_moving;
  // Beware! nsources_local_moving is a pointer to an integer, not to the array of size NSTEP

  // check if anything to do
  if (nsources_local == 0) return;

  int nsources = *nsources_f;
  int NSTEP = *NSTEP_f;

  int it = *it_f - 1;
  int iphase = *iphase_f;

  // only adds this contribution for first pass
  if (iphase != 1) return;

  int num_blocks_x, num_blocks_y;

  // note: we will launch kernels only for local sources since sourcearrays_moving(..) and ispec_selected_souce_movie(..)
  //       are sorted along local sources only.
  //       this avoids the need to have an array like islice_selected_source(NSOURCES,NSTEP) on the GPU.

  // Look up for the best way to distribute the sources in a grid (one block per source)
  get_blocks_xy(nsources_local, &num_blocks_x, &num_blocks_y);
   // create the grid (one block per source)
  dim3 grid(num_blocks_x,num_blocks_y);
  // create the threads in the block (one thread per GLL), it has NGLLX,NGLLX,1 structure
  dim3 threads(NGLLX,NGLLX,1);

  // cudaMemoryTest(3); // Useful to check memory
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

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_moving_sources_ac_cuda");
#endif

}


/* ----------------------------------------------------------------------------------------------- */

// acoustic adjoint sources

/* ----------------------------------------------------------------------------------------------- */

__global__ void add_sources_ac_SIM_TYPE_2_OR_3_kernel(realw* potential_dot_dot_acoustic,
                                                      realw* source_adjointe,
                                                      realw* xir_store,
                                                      realw* gammar_store,
                                                      int* d_ibool,
                                                      int* ispec_is_acoustic,
                                                      int* ispec_selected_rec_loc,
                                                      int it,
                                                      int nadj_rec_local,
                                                      realw* kappastore,
                                                      int NSTEP ) {

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  // because of grid shape, irec_local can be too big
  if (irec_local < nadj_rec_local) {

    int ispec = ispec_selected_rec_loc[irec_local] - 1;

    if (ispec_is_acoustic[ispec]) {
      int i = threadIdx.x;
      int j = threadIdx.y;

      int iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

      realw  kappal = kappastore[INDEX3(NGLLX,NGLLX,i,j,ispec)];
      realw  xir = xir_store[INDEX2(nadj_rec_local,irec_local,i)];
      realw  gammar = gammar_store[INDEX2(nadj_rec_local,irec_local,j)];
      realw  source_adj = source_adjointe[INDEX3(nadj_rec_local,NSTEP,irec_local,it,0)];


      // beware, for acoustic medium, a pressure source would be taking the negative
      // and divide by Kappa of the fluid;
      // this would have to be done when constructing the adjoint source.
      //
      //          the idea is to have e.g. a pressure source, where all 3 components would be the same
      realw stf = source_adj * gammar * xir / kappal ;

      atomicAdd(&potential_dot_dot_acoustic[iglob],-stf);
      // Alexis Bottero added a - sign for consistency with CPU version
      //atomicAdd(&potential_dot_dot_acoustic[iglob],stf);
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(add_sources_ac_sim_2_or_3_cuda,
              ADD_SOURCES_AC_SIM_2_OR_3_CUDA)(long* Mesh_pointer,
                                               int* iphasef,
                                               int* it,
                                               int* nadj_rec_local,
                                               int* NSTEP) {

  TRACE("add_sources_ac_sim_2_or_3_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int iphase = *iphasef;

  // only adds this contributions for first pass
  if (iphase != 1) return;

  // checks
  if (*nadj_rec_local != mp->nadj_rec_local) exit_on_cuda_error("add_sources_ac_sim_type_2_or_3: nadj_rec_local not equal\n");

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nadj_rec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLLX,NGLLX,1);

  int it_index = (*it) - 1;

  // launches cuda kernel for acoustic adjoint sources
  add_sources_ac_SIM_TYPE_2_OR_3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                                mp->d_source_adjointe,
                                                                                mp->d_xir_store_loc,
                                                                                mp->d_gammar_store_loc,
                                                                                mp->d_ibool,
                                                                                mp->d_ispec_is_acoustic,
                                                                                mp->d_ispec_selected_rec_loc,
                                                                                it_index,
                                                                                mp->nadj_rec_local,
                                                                                mp->d_kappastore,
                                                                                *NSTEP);

//  print_CUDA_error_if_any(cudaStreamSynchronize(mp->compute_stream),38);
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("add_sources_acoustic_SIM_TYPE_2_OR_3_kernel");
#endif
}
