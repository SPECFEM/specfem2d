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

__global__ void compute_add_sources_kernel(realw* accel,
                                           int* d_ibool,
                                           realw* sourcearrays,
                                           realw* d_source_time_function,
                                           int myrank,
                                           int* ispec_selected_source,
                                           int* ispec_is_elastic,
                                           int it,int nsources_local) {
  int i = threadIdx.x;
  int j = threadIdx.y;

  int isource  = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int ispec,iglob;
  realw stf;

  if (isource < nsources_local) { // when NSOURCES > 65535, but mod(nspec_top,2) > 0, we end up with an extra block.

      ispec = ispec_selected_source[isource]-1;

      if (ispec_is_elastic[ispec]) {


        stf = d_source_time_function[INDEX2(nsources_local,isource,it)];
        iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;


        atomicAdd(&accel[iglob*2],sourcearrays[INDEX4(nsources_local,NDIM,NGLLX,isource, 0,i,j)]*stf);
        atomicAdd(&accel[iglob*2+1],sourcearrays[INDEX4(nsources_local,NDIM,NGLLX,isource, 1,i,j)]*stf);



    }
  }

}


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


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_el_cuda");
#endif
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

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nsources_local,&num_blocks_x,&num_blocks_y);
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLX,1);
  int it = *itf -1;

  compute_add_sources_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel,mp->d_ibool,
                                                                    mp->d_sourcearrays,
                                                                    mp->d_source_time_function,
                                                                    mp->myrank,
                                                                    mp->d_ispec_selected_source,
                                                                    mp->d_ispec_is_elastic,
                                                                    it,mp->nsources_local);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_el_s3_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// ADJOINT sources

/* ----------------------------------------------------------------------------------------------- */

__global__ void add_sources_el_SIM_TYPE_2_OR_3_kernel(realw* accel,
                                                      realw* source_adjointe,
                                                      realw* xir_store,
                                                      realw* gammar_store,
                                                      int* d_ibool,
                                                      int* ispec_is_elastic,
                                                      int* ispec_selected_rec_loc,
                                                      int it,
                                                      int nadj_rec_local,
                                                      int NSTEP) {

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  if (irec_local < nadj_rec_local) { // when nrec > 65535, but mod(nspec_top,2) > 0, we end up with an extra block.

    int ispec = ispec_selected_rec_loc[irec_local]-1;

    if (ispec_is_elastic[ispec]) {
      int i = threadIdx.x;
      int j = threadIdx.y;
      int iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)]-1;

      realw  xir = xir_store[INDEX2(nadj_rec_local,irec_local,i)];

      realw  gammar = gammar_store[INDEX2(nadj_rec_local,irec_local,j)];

      realw  source_adjx = source_adjointe[INDEX3(nadj_rec_local,NSTEP,irec_local,it,0)];

      realw  source_adjz = source_adjointe[INDEX3(nadj_rec_local,NSTEP,irec_local,it,1)];


      // atomic operations are absolutely necessary for correctness!
      atomicAdd(&accel[2*iglob],source_adjx * gammar * xir);
      atomicAdd(&accel[1+2*iglob], source_adjz * gammar * xir);
    } // ispec_is_elastic
  }

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(add_sources_el_sim_type_2_or_3,
              ADD_SOURCES_EL_SIM_TYPE_2_OR_3)(long* Mesh_pointer,
                                               int* iphasef,
                                               int* it,
                                               int* nadj_rec_local,
                                               int* NSTEP) {

  TRACE("add_sources_el_sim_type_2_or_3");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks
  if (*nadj_rec_local != mp->nadj_rec_local) exit_on_error("add_sources_el_sim_type_2_or_3: nadj_rec_local not equal\n");

  int iphase = *iphasef;

  // only add this contribution for first pass
  if (iphase != 1) return;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nadj_rec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLLX,NGLLX,1);

  int it_index = (*it) - 1;


  // the irec_local variable needs to be precomputed (as
  // h_pre_comp..), because normally it is in the loop updating accel,
  // and due to how it's incremented, it cannot be parallelized

  add_sources_el_SIM_TYPE_2_OR_3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,
                                                                               mp->d_source_adjointe,
                                                                               mp->d_xir_store_loc,
                                                                               mp->d_gammar_store_loc,
                                                                               mp->d_ibool,
                                                                               mp->d_ispec_is_elastic,
                                                                               mp->d_ispec_selected_rec_loc,
                                                                               it_index,
                                                                               mp->nadj_rec_local,
                                                                               *NSTEP);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("add_sources_SIM_TYPE_2_OR_3_kernel");
#endif
}

