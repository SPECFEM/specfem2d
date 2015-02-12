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

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"




/* ----------------------------------------------------------------------------------------------- */

// ACOUSTIC - ELASTIC coupling

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_coupling_acoustic_el_kernel(realw* displ,
                                                    realw* potential_dot_dot_acoustic,
                                                    int num_coupling_ac_el_faces,
                                                    int* coupling_ac_el_ispec,
                                                    int* coupling_ac_el_ij,
                                                    realw* coupling_ac_el_normal,
                                                    realw* coupling_ac_el_jacobian1Dw,
                                                    int* d_ibool,
                                                    int* ispec_is_inner,
                                                    int phase_is_inner) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,iglob,ispec;
  realw displ_x,displ_z,displ_n;
  realw nx,nz;
  realw jacobianw;

  if( iface < num_coupling_ac_el_faces){

    // don't compute points outside NGLLSQUARE==NGLL2==25
    // way 2: no further check needed since blocksize = 25
    //  if(igll<NGLL2) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = coupling_ac_el_ispec[iface] - 1;

    if(ispec_is_inner[ispec] == phase_is_inner ) {

      i = coupling_ac_el_ij[INDEX3(NDIM,NGLLX,0,igll,iface)] - 1;
      j = coupling_ac_el_ij[INDEX3(NDIM,NGLLX,1,igll,iface)] - 1;

      iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

      // elastic displacement on global point
      displ_x = displ[iglob*2] ; // (1,iglob)
      displ_z = displ[iglob*2+1] ; // (3,iglob)

      // gets associated normal on GLL point
      nx = coupling_ac_el_normal[INDEX3(NDIM,NGLLX,0,igll,iface)]; // (1,igll,iface)
      nz = coupling_ac_el_normal[INDEX3(NDIM,NGLLX,1,igll,iface)]; // (3,igll,iface)

      // calculates displacement component along normal
      // (normal points outwards of acoustic element)
      displ_n = displ_x*nx + displ_z*nz;

      // gets associated, weighted jacobian
      jacobianw = coupling_ac_el_jacobian1Dw[INDEX2(NGLLX,igll,iface)];

      // continuity of pressure and normal displacement on global point

      // note: Newmark time scheme together with definition of scalar potential:
      //          pressure = - chi_dot_dot
      //          requires that this coupling term uses the updated displacement at time step [t+delta_t],
      //          which is done at the very beginning of the time loop
      //          (see e.g. Chaljub & Vilotte, Nissen-Meyer thesis...)
      //          it also means you have to calculate and update this here first before
      //          calculating the coupling on the elastic side for the acceleration...
      atomicAdd(&potential_dot_dot_acoustic[iglob],+ jacobianw*displ_n);

    }
  //  }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_coupling_ac_el_cuda,
              COMPUTE_COUPLING_AC_EL_CUDA)(long* Mesh_pointer,
                                           int* phase_is_innerf,
                                           int* num_coupling_ac_el_facesf) {
  TRACE("compute_coupling_ac_el_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  int phase_is_inner            = *phase_is_innerf;
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
                                                       mp->d_ibool,
                                                       mp->d_ispec_is_inner,
                                                       phase_is_inner);

  //  adjoint simulations
  if (mp->simulation_type == 3 ){
    compute_coupling_acoustic_el_kernel<<<grid,threads>>>(mp->d_b_displ,
                                                          mp->d_b_potential_dot_dot_acoustic,
                                                          num_coupling_ac_el_faces,
                                                          mp->d_coupling_ac_el_ispec,
                                                          mp->d_coupling_ac_el_ijk,
                                                          mp->d_coupling_ac_el_normal,
                                                          mp->d_coupling_ac_el_jacobian2Dw,
                                                          mp->d_ibool,
                                                          mp->d_ispec_is_inner,
                                                          phase_is_inner);

  }


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_coupling_acoustic_el_kernel");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// ELASTIC - ACOUSTIC coupling

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_coupling_elastic_ac_kernel(realw* potential_dot_dot_acoustic,
                                                    realw* accel,
                                                    int num_coupling_ac_el_faces,
                                                    int* coupling_ac_el_ispec,
                                                    int* coupling_ac_el_ij,
                                                    realw* coupling_ac_el_normal,
                                                    realw* coupling_ac_el_jacobian1Dw,
                                                    int* d_ibool,
                                                    int* ispec_is_inner,
                                                    int phase_is_inner) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,iglob,ispec;
  realw pressure;
  realw nx,nz;
  realw jacobianw;


  if( iface < num_coupling_ac_el_faces){

    // don't compute points outside NGLLSQUARE==NGLL2==25
    // way 2: no further check needed since blocksize = 25
    //  if(igll<NGLL2) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = coupling_ac_el_ispec[iface] - 1;

    if(ispec_is_inner[ispec] == phase_is_inner ) {

      i = coupling_ac_el_ij[INDEX3(NDIM,NGLLX,0,igll,iface)] - 1;
      j = coupling_ac_el_ij[INDEX3(NDIM,NGLLX,1,igll,iface)] - 1;


      iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

      // gets associated normal on GLL point
      // note: normal points away from acoustic element
      nx = coupling_ac_el_normal[INDEX3(NDIM,NGLLX,0,igll,iface)]; // (1,igll,iface)
      nz = coupling_ac_el_normal[INDEX3(NDIM,NGLLX,1,igll,iface)]; // (3,igll,iface)

      // gets associated, weighted jacobian
      jacobianw = coupling_ac_el_jacobian1Dw[INDEX2(NGLLX,igll,iface)];


        pressure = - potential_dot_dot_acoustic[iglob];


      // continuity of displacement and pressure on global point
      //
      // note: Newmark time scheme together with definition of scalar potential:
      //          pressure = - chi_dot_dot
      //          requires that this coupling term uses the *UPDATED* pressure (chi_dot_dot), i.e.
      //          pressure at time step [t + delta_t]
      //          (see e.g. Chaljub & Vilotte, Nissen-Meyer thesis...)
      //          it means you have to calculate and update the acoustic pressure first before
      //          calculating this term...
      atomicAdd(&accel[iglob*2],+ jacobianw*nx*pressure);
      atomicAdd(&accel[iglob*2+1],+ jacobianw*nz*pressure);
    }
    //  }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_coupling_el_ac_cuda,
              COMPUTE_COUPLING_EL_AC_CUDA)(long* Mesh_pointer,
                                           int* phase_is_innerf,
                                           int* num_coupling_ac_el_facesf) {
  TRACE("compute_coupling_el_ac_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  int phase_is_inner            = *phase_is_innerf;
  int num_coupling_ac_el_faces  = *num_coupling_ac_el_facesf;

  // way 1: exact blocksize to match NGLLX
  int blocksize = 5;

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
                                                       mp->d_ibool,
                                                       mp->d_ispec_is_inner,
                                                       phase_is_inner);

  //  adjoint simulations
  if (mp->simulation_type == 3 ){
    compute_coupling_elastic_ac_kernel<<<grid,threads>>>(mp->d_b_potential_dot_dot_acoustic,
                                                         mp->d_b_accel,
                                                         num_coupling_ac_el_faces,
                                                         mp->d_coupling_ac_el_ispec,
                                                         mp->d_coupling_ac_el_ijk,
                                                         mp->d_coupling_ac_el_normal,
                                                         mp->d_coupling_ac_el_jacobian2Dw,
                                                         mp->d_ibool,
                                                         mp->d_ispec_is_inner,
                                                         phase_is_inner);

  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_coupling_el_ac_cuda");
#endif
}
