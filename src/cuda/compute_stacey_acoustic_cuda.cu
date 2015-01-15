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

__global__ void compute_stacey_acoustic_kernel(realw* potential_dot_acoustic,
                                               realw* potential_dot_dot_acoustic,
                                               int* abs_boundary_ispec,
                                               int* abs_boundary_ij,
                                               realw* abs_boundary_jacobian1Dw,
                                               int* d_ibool,
                                               realw* rhostore,
                                               realw* kappastore,
                                               int* ispec_is_inner,
                                               int* ispec_is_acoustic,
                                               int phase_is_inner,
                                               int SIMULATION_TYPE,
                                               int SAVE_FORWARD,
                                               int num_abs_boundary_faces,
                                               realw* b_potential_dot_acoustic,
                                               realw* b_potential_dot_dot_acoustic,
                                               realw* b_absorb_potential_left,
                                               realw* b_absorb_potential_right,
                                               realw* b_absorb_potential_top,
                                               realw* b_absorb_potential_bottom,
                                               int* ib_left,
                                               int* ib_right,
                                               int* ib_top,
                                               int* ib_bottom,
                                               int* cote_abs) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,iglob,ispec,num_local;
  realw rhol,kappal,cpl;
  realw jacobianw;
  realw vel;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if( iface < num_abs_boundary_faces){

  //  if(igll<NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if(ispec_is_inner[ispec] == phase_is_inner && ispec_is_acoustic[ispec] ) {

      i = abs_boundary_ij[INDEX3(NDIM,NGLLX,0,igll,iface)]-1;
      j = abs_boundary_ij[INDEX3(NDIM,NGLLX,1,igll,iface)]-1;


      iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)]-1;

      // determines bulk sound speed
      rhol = rhostore[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)];

      kappal = kappastore[INDEX3(NGLLX,NGLLX,i,j,ispec)];

      cpl = sqrt( kappal / rhol );


        // uses a potential definition of: s = 1/rho grad(chi)
        vel = potential_dot_acoustic[iglob] / rhol;


      // gets associated, weighted jacobian
      jacobianw = abs_boundary_jacobian1Dw[INDEX2(NGLLX,igll,iface)];

     // Sommerfeld condition
      atomicAdd(&potential_dot_dot_acoustic[iglob],-vel*jacobianw/cpl);

      // adjoint simulations
      if( SIMULATION_TYPE == 3 ){

if (cote_abs[iface] == 1){ num_local = ib_bottom[iface] -1 ; atomicAdd(&b_potential_dot_dot_acoustic[iglob],
                                     -b_absorb_potential_bottom[INDEX2(NGLLX,igll,num_local)]);}
else if (cote_abs[iface] == 2){ num_local = ib_right[iface] -1 ; atomicAdd(&b_potential_dot_dot_acoustic[iglob],
                                     -b_absorb_potential_right[INDEX2(NGLLX,igll,num_local)]);}
else if (cote_abs[iface] == 3){ num_local = ib_top[iface] -1 ; atomicAdd(&b_potential_dot_dot_acoustic[iglob],
                                     -b_absorb_potential_top[INDEX2(NGLLX,igll,num_local)]);}
else if (cote_abs[iface] == 4){ num_local = ib_left[iface] -1 ; atomicAdd(&b_potential_dot_dot_acoustic[iglob],
                                     -b_absorb_potential_left[INDEX2(NGLLX,igll,num_local)]);}

}else if( SIMULATION_TYPE == 1 && SAVE_FORWARD ){
        // saves boundary values
if (cote_abs[iface] == 1) { num_local = ib_bottom[iface] -1 ; b_absorb_potential_bottom[INDEX2(NGLLX,igll,num_local)] = vel*jacobianw/cpl;}
else if (cote_abs[iface] == 2) { num_local = ib_right[iface] -1 ; b_absorb_potential_right[INDEX2(NGLLX,igll,num_local)] = vel*jacobianw/cpl;}
else if (cote_abs[iface] == 3) { num_local = ib_top[iface] -1 ; b_absorb_potential_top[INDEX2(NGLLX,igll,num_local)] = vel*jacobianw/cpl;}
else if (cote_abs[iface] == 4) { num_local = ib_left[iface] -1 ; b_absorb_potential_left[INDEX2(NGLLX,igll,num_local)] = vel*jacobianw/cpl;}

      }
    }
//  }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_stacey_acoustic_cuda,
              COMPUTE_STACEY_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                            int* phase_is_innerf,
                                            realw* h_b_absorb_potential_left,
                                            realw* h_b_absorb_potential_right,
                                            realw* h_b_absorb_potential_top,
                                            realw* h_b_absorb_potential_bottom) {
TRACE("compute_stacey_acoustic_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( mp->d_num_abs_boundary_faces == 0 ) return;

  int phase_is_inner          = *phase_is_innerf;

  // way 1: Elapsed time: 4.385948e-03
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //  int blocksize = 32;

  // way 2: Elapsed time: 4.379034e-03
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLLX;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->d_num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);




  //  adjoint simulations: reads in absorbing boundary
  if (mp->simulation_type == 3 && phase_is_inner == 0){
    // copies array to GPU
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_potential_left,h_b_absorb_potential_left,
                                       mp->d_nspec_left*sizeof(realw)*NGLLX,cudaMemcpyHostToDevice),7700);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_potential_right,h_b_absorb_potential_right,
                                       mp->d_nspec_right*sizeof(realw)*NGLLX,cudaMemcpyHostToDevice),7700);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_potential_top,h_b_absorb_potential_top,
                                       mp->d_nspec_top*sizeof(realw)*NGLLX,cudaMemcpyHostToDevice),7700);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_potential_bottom,h_b_absorb_potential_bottom,
                                       mp->d_nspec_bottom*sizeof(realw)*NGLLX,cudaMemcpyHostToDevice),7700);
  }


  compute_stacey_acoustic_kernel<<<grid,threads>>>(mp->d_potential_dot_acoustic,
                                                   mp->d_potential_dot_dot_acoustic,
                                                   mp->d_abs_boundary_ispec,
                                                   mp->d_abs_boundary_ijk,
                                                   mp->d_abs_boundary_jacobian2Dw,
                                                   mp->d_ibool,
                                                   mp->d_rhostore,
                                                   mp->d_kappastore,
                                                   mp->d_ispec_is_inner,
                                                   mp->d_ispec_is_acoustic,
                                                   phase_is_inner,
                                                   mp->simulation_type,
                                                   mp->save_forward,
                                                   mp->d_num_abs_boundary_faces,
                                                   mp->d_b_potential_dot_acoustic,
                                                   mp->d_b_potential_dot_dot_acoustic,
                                                   mp->d_b_absorb_potential_left,
                                                   mp->d_b_absorb_potential_right,
                                                   mp->d_b_absorb_potential_top,
                                                   mp->d_b_absorb_potential_bottom,
                                                   mp->d_ib_left,
                                                   mp->d_ib_right,
                                                   mp->d_ib_top,
                                                   mp->d_ib_bottom,
                                                   mp->d_cote_abs);

  //  adjoint simulations: stores absorbed wavefield part
  if (mp->simulation_type == 1 && mp->save_forward && phase_is_inner==1 ){
    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    // copies array to CPU
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_potential_left,mp->d_b_absorb_potential_left,
                                       mp->d_nspec_left*sizeof(realw)*NGLLX,cudaMemcpyDeviceToHost),7701);
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_potential_right,mp->d_b_absorb_potential_right,
                                       mp->d_nspec_right*sizeof(realw)*NGLLX,cudaMemcpyDeviceToHost),7702);
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_potential_top,mp->d_b_absorb_potential_top,
                                       mp->d_nspec_top*sizeof(realw)*NGLLX,cudaMemcpyDeviceToHost),7703);
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_potential_bottom,mp->d_b_absorb_potential_bottom,
                                       mp->d_nspec_bottom*sizeof(realw)*NGLLX,cudaMemcpyDeviceToHost),7704);

  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_stacey_acoustic_kernel");
#endif
}

