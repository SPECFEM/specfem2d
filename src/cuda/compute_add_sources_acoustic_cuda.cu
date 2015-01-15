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
#include <cuda.h>
#include <cublas.h>

#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// acoustic sources

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_add_sources_acoustic_kernel(realw* potential_dot_dot_acoustic,
                                                    int* d_ibool,
                                                    int* ispec_is_inner,
                                                    int phase_is_inner,
                                                    realw* sourcearrays,
                                                    realw* source_time_function,
                                                    int myrank,
                                                    int* ispec_selected_source,
                                                    int* ispec_is_acoustic,
                                                    realw* kappastore,
                                                    int it,int nsources_local) {
  int i = threadIdx.x;
  int j = threadIdx.y;

  int isource  = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int ispec,iglob;
  realw stf,kappal;

  if( isource < nsources_local ){

      ispec = ispec_selected_source[isource]-1;

      if(ispec_is_inner[ispec] == phase_is_inner && ispec_is_acoustic[ispec] ) {

        iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

        kappal = kappastore[INDEX3(NGLLX,NGLLX,i,j,ispec)];

        stf = source_time_function[INDEX2(nsources_local,isource,it)]/kappal;
        atomicAdd(&potential_dot_dot_acoustic[iglob],
                  -sourcearrays[INDEX4(nsources_local,NDIM,NGLLX,isource, 0,i,j)]*stf);



    }
  }
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_ac_cuda,
              COMPUTE_ADD_SOURCES_AC_CUDA)(long* Mesh_pointer,
                                           int* phase_is_innerf,
                                           int * itf) {

  TRACE("compute_add_sources_ac_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // check if anything to do
  if( mp->nsources_local == 0 ) return;

  int phase_is_inner = *phase_is_innerf;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nsources_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,1);

  int it = *itf -1;



  compute_add_sources_acoustic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                              mp->d_ibool,
                                                                              mp->d_ispec_is_inner,
                                                                              phase_is_inner,
                                                                              mp->d_sourcearrays,
                                                                              mp->d_source_time_function,
                                                                              mp->myrank,
                                                                              mp->d_ispec_selected_source,
                                                                              mp->d_ispec_is_acoustic,
                                                                              mp->d_kappastore,
                                                                              it,mp->nsources_local);

  print_CUDA_error_if_any(cudaStreamSynchronize(mp->compute_stream),mp->myrank);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_ac_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_ac_s3_cuda,
              COMPUTE_ADD_SOURCES_AC_s3_CUDA)(long* Mesh_pointer,
                                              int* phase_is_innerf,
                                              int* itf) {

  TRACE("compute_add_sources_ac_s3_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // check if anything to do
  if( mp->nsources_local == 0 ) return;

  int phase_is_inner = *phase_is_innerf;
  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nsources_local,&num_blocks_x,&num_blocks_y);
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,1);
  int it = *itf -1;

  compute_add_sources_acoustic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_dot_dot_acoustic,
                                                                              mp->d_ibool,
                                                                              mp->d_ispec_is_inner,
                                                                              phase_is_inner,
                                                                              mp->d_sourcearrays,
                                                                              mp->d_source_time_function,
                                                                              mp->myrank,
                                                                              mp->d_ispec_selected_source,
                                                                              mp->d_ispec_is_acoustic,
                                                                              mp->d_kappastore,
                                                                              it,mp->nsources_local);


      print_CUDA_error_if_any(cudaStreamSynchronize(mp->compute_stream),38);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_ac_s3_cuda");
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
                                                      int* ispec_is_inner,
                                                      int* ispec_is_acoustic,
                                                      int* ispec_selected_rec,
                                                      int phase_is_inner,
                                                      int it,
                                                      int* pre_computed_irec,
                                                      int nadj_rec_local,
                                                      realw* kappastore,
                                                      int NSTEP  ) {

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  // because of grid shape, irec_local can be too big
  if(irec_local < nadj_rec_local) {

    int irec = pre_computed_irec[irec_local];

    int ispec = ispec_selected_rec[irec]-1;
    if( ispec_is_acoustic[ispec] ){

      // checks if element is in phase_is_inner run
      if(ispec_is_inner[ispec] == phase_is_inner) {
        int i = threadIdx.x;
        int j = threadIdx.y;


        int iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)]-1;

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

        atomicAdd(&potential_dot_dot_acoustic[iglob],stf);


      }
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(add_sources_ac_sim_2_or_3_cuda,
              ADD_SOURCES_AC_SIM_2_OR_3_CUDA)(long* Mesh_pointer,
                                               int* phase_is_inner,
                                               int* it,
                                               int* nadj_rec_local,
                                               int* NSTEP) {

  TRACE("add_sources_ac_sim_2_or_3_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks
  if( *nadj_rec_local != mp->nadj_rec_local) exit_on_cuda_error("add_sources_ac_sim_type_2_or_3: nadj_rec_local not equal\n");

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nadj_rec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(5,5,1);

  int it_index = (*it) - 1;


  // launches cuda kernel for acoustic adjoint sources
  add_sources_ac_SIM_TYPE_2_OR_3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                                mp->d_source_adjointe,
                                                                                mp->d_xir_store_loc,
                                                                                mp->d_gammar_store_loc,
                                                                                mp->d_ibool,
                                                                                mp->d_ispec_is_inner,
                                                                                mp->d_ispec_is_acoustic,
                                                                                mp->d_ispec_selected_rec,
                                                                                *phase_is_inner,
                                                                                it_index,
                                                                                mp->d_pre_computed_irec,
                                                                                mp->nadj_rec_local,
                                                                                mp->d_kappastore,
                                                                                *NSTEP);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("add_sources_acoustic_SIM_TYPE_2_OR_3_kernel");
#endif
}
