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

__global__ void compute_stacey_elastic_kernel(realw* veloc,
                                              realw* accel,
                                              int* abs_boundary_ispec,
                                              int* abs_boundary_ij,
                                              realw* abs_boundary_normal,
                                              realw* abs_boundary_jacobian1Dw,
                                              int* d_ibool,
                                              realw* rho_vp,
                                              realw* rho_vs,
                                              int* ispec_is_elastic,
                                              int SIMULATION_TYPE,
                                              int SAVE_FORWARD,
                                              int num_abs_boundary_faces,
                                              realw* b_absorb_elastic_left,
                                              realw* b_absorb_elastic_right,
                                              realw* b_absorb_elastic_top,
                                              realw* b_absorb_elastic_bottom,
                                              int* ib_left,
                                              int* ib_right,
                                              int* ib_top,
                                              int* ib_bottom,
                                              int* cote_abs) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,iglob,ispec,num_local;
  realw vx,vz,vn;
  realw nx,nz;
  realw rho_vp_temp,rho_vs_temp;
  realw tx,tz;
  realw jacobianw;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if (iface < num_abs_boundary_faces){

  //if (igll < NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if (ispec_is_elastic[ispec]) {

      i = abs_boundary_ij[INDEX3(NDIM,NGLLX,0,igll,iface)]-1;
      j = abs_boundary_ij[INDEX3(NDIM,NGLLX,1,igll,iface)]-1;
      //check if the point must be computed
      if (i==NGLLX-1 || j==NGLLX-1) return;

      iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)]-1;

      // gets associated velocity

      vx = veloc[iglob*2];
      vz = veloc[iglob*2+1];

      // gets associated normal
      nx = abs_boundary_normal[INDEX3(NDIM,NGLLX,0,igll,iface)];
      nz = abs_boundary_normal[INDEX3(NDIM,NGLLX,1,igll,iface)];

      // // velocity component in normal direction (normal points out of element)
      vn = vx*nx + vz*nz;

      rho_vp_temp = rho_vp[INDEX3(NGLLX,NGLLX,i,j,ispec)];
      rho_vs_temp = rho_vs[INDEX3(NGLLX,NGLLX,i,j,ispec)];

      tx = rho_vp_temp*vn*nx + rho_vs_temp*(vx-vn*nx);
      tz = rho_vp_temp*vn*nz + rho_vs_temp*(vz-vn*nz);

      jacobianw = abs_boundary_jacobian1Dw[INDEX2(NGLLX,igll,iface)];

      atomicAdd(&accel[iglob*2],-tx*jacobianw);
      atomicAdd(&accel[iglob*2+1],-tz*jacobianw);

      if (SAVE_FORWARD && SIMULATION_TYPE == 1) {

      if (cote_abs[iface] == 1) {num_local = ib_bottom[iface]-1;
                                b_absorb_elastic_bottom[INDEX3(NDIM,NGLLX,0,igll,num_local)] = tx*jacobianw;
                                b_absorb_elastic_bottom[INDEX3(NDIM,NGLLX,1,igll,num_local)] = tz*jacobianw;}
      else if (cote_abs[iface] == 2)   {num_local = ib_right[iface]-1;
                                 b_absorb_elastic_right[INDEX3(NDIM,NGLLX,0,igll,num_local)] = tx*jacobianw;
                                b_absorb_elastic_right[INDEX3(NDIM,NGLLX,1,igll,num_local)] = tz*jacobianw;}
      else if (cote_abs[iface] == 3)  {num_local = ib_top[iface]-1;
                                 b_absorb_elastic_top[INDEX3(NDIM,NGLLX,0,igll,num_local)] = tx*jacobianw;
                                b_absorb_elastic_top[INDEX3(NDIM,NGLLX,1,igll,num_local)] = tz*jacobianw;}
      else if (cote_abs[iface] == 4){num_local = ib_left[iface]-1;
                                b_absorb_elastic_left[INDEX3(NDIM,NGLLX,0,igll,num_local)] = tx*jacobianw;
                                b_absorb_elastic_left[INDEX3(NDIM,NGLLX,1,igll,num_local)] = tz*jacobianw;}

      } // SIMULATION_TYPE
    }
  } // num_abs_boundary_faces

}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_stacey_elastic_sim3_kernel(int* abs_boundary_ispec,
                                                   int* abs_boundary_ijk,
                                                   int* d_ibool,
                                                   int* ispec_is_elastic,
                                                   int num_abs_boundary_faces,
                                                   realw* b_accel,
                                                   realw* b_absorb_elastic_left,
                                                   realw* b_absorb_elastic_right,
                                                   realw* b_absorb_elastic_top,
                                                   realw* b_absorb_elastic_bottom,
                                                   int* ib_left,
                                                   int* ib_right,
                                                   int* ib_top,
                                                   int* ib_bottom,
                                                   int* d_cote_abs) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,iglob,ispec,num_local;

  if (iface < num_abs_boundary_faces){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if (ispec_is_elastic[ispec]) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLLX,0,igll,iface)]-1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLLX,1,igll,iface)]-1;

      iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)]-1;

      if (d_cote_abs[iface] == 1){
        num_local= ib_bottom[iface]-1;
        atomicAdd(&b_accel[iglob*2 ], -b_absorb_elastic_bottom[INDEX3(NDIM,NGLLX,0,igll,num_local)]);
        atomicAdd(&b_accel[iglob*2+1 ], -b_absorb_elastic_bottom[INDEX3(NDIM,NGLLX,1,igll,num_local)]);

      } else if (d_cote_abs[iface] == 2){
        num_local= ib_right[iface]-1;
        atomicAdd(&b_accel[iglob*2 ], -b_absorb_elastic_right[INDEX3(NDIM,NGLLX,0,igll,num_local)]);
        atomicAdd(&b_accel[iglob*2+1 ], -b_absorb_elastic_right[INDEX3(NDIM,NGLLX,1,igll,num_local)]);

      } else if (d_cote_abs[iface] == 3){
        num_local= ib_top[iface]-1;
        atomicAdd(&b_accel[iglob*2 ], -b_absorb_elastic_top[INDEX3(NDIM,NGLLX,0,igll,num_local)]);
        atomicAdd(&b_accel[iglob*2+1 ], -b_absorb_elastic_top[INDEX3(NDIM,NGLLX,1,igll,num_local)]);

      } else if (d_cote_abs[iface] == 4){
        num_local= ib_left[iface]-1;
        atomicAdd(&b_accel[iglob*2 ], -b_absorb_elastic_left[INDEX3(NDIM,NGLLX,0,igll,num_local)]);
        atomicAdd(&b_accel[iglob*2+1 ], -b_absorb_elastic_left[INDEX3(NDIM,NGLLX,1,igll,num_local)]);
      }
    }
  } // num_abs_boundary_faces

}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_stacey_viscoelastic_cuda,
              COMPUTE_STACEY_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                           int* iphasef,
                                           realw* h_b_absorb_elastic_left,
                                           realw* h_b_absorb_elastic_right,
                                           realw* h_b_absorb_elastic_top,
                                           realw* h_b_absorb_elastic_bottom) {

  TRACE("compute_stacey_viscoelastic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if (mp->d_num_abs_boundary_faces == 0) return;

  int iphase    = *iphasef;

  // only add this contribution for first pass
  if (iphase != 1) return;

  int blocksize = NGLLX;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->d_num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if (mp->simulation_type == 3) {
    // reading is done in fortran routine
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_elastic_left,h_b_absorb_elastic_left,
                                       2*mp->d_nspec_left*sizeof(realw)*NGLLX,cudaMemcpyHostToDevice),7700);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_elastic_right,h_b_absorb_elastic_right,
                                       2*mp->d_nspec_right*sizeof(realw)*NGLLX,cudaMemcpyHostToDevice),7700);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_elastic_top,h_b_absorb_elastic_top,
                                       2*mp->d_nspec_top*sizeof(realw)*NGLLX,cudaMemcpyHostToDevice),7700);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_elastic_bottom,h_b_absorb_elastic_bottom,
                                       2*mp->d_nspec_bottom*sizeof(realw)*NGLLX,cudaMemcpyHostToDevice),7700);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("between cudamemcpy and compute_stacey_elastic_kernel");
#endif

  compute_stacey_elastic_kernel<<<grid,threads,0,mp->compute_stream>>>( mp->d_veloc,
                                                                        mp->d_accel,
                                                                        mp->d_abs_boundary_ispec,
                                                                        mp->d_abs_boundary_ijk,
                                                                        mp->d_abs_boundary_normal,
                                                                        mp->d_abs_boundary_jacobian2Dw,
                                                                        mp->d_ibool,
                                                                        mp->d_rho_vp,
                                                                        mp->d_rho_vs,
                                                                        mp->d_ispec_is_elastic,
                                                                        mp->simulation_type,
                                                                        mp->save_forward,
                                                                        mp->d_num_abs_boundary_faces,
                                                                        mp->d_b_absorb_elastic_left,
                                                                        mp->d_b_absorb_elastic_right,
                                                                        mp->d_b_absorb_elastic_top,
                                                                        mp->d_b_absorb_elastic_bottom,
                                                                        mp->d_ib_left,
                                                                        mp->d_ib_right,
                                                                        mp->d_ib_top,
                                                                        mp->d_ib_bottom,
                                                                        mp->d_cote_abs);

  // adjoint simulations
  if (mp->simulation_type == 3) {
    compute_stacey_elastic_sim3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_abs_boundary_ispec,
                                                                              mp->d_abs_boundary_ijk,
                                                                              mp->d_ibool,
                                                                              mp->d_ispec_is_elastic,
                                                                              mp->d_num_abs_boundary_faces,
                                                                              mp->d_b_accel,
                                                                              mp->d_b_absorb_elastic_left,
                                                                              mp->d_b_absorb_elastic_right,
                                                                              mp->d_b_absorb_elastic_top,
                                                                              mp->d_b_absorb_elastic_bottom,
                                                                              mp->d_ib_left,
                                                                              mp->d_ib_right,
                                                                              mp->d_ib_top,
                                                                              mp->d_ib_bottom,
                                                                              mp->d_cote_abs);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_stacey_elastic_kernel");
#endif

  if (mp->simulation_type == 1 && mp->save_forward) {
    // explicitly wait until compute stream is done
    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    cudaStreamSynchronize(mp->compute_stream);

    // writing is done in fortran routine
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_elastic_left,mp->d_b_absorb_elastic_left,
                                       2*mp->d_nspec_left*sizeof(realw)*NGLLX,cudaMemcpyDeviceToHost),7701);
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_elastic_right,mp->d_b_absorb_elastic_right,
                                       2*mp->d_nspec_right*sizeof(realw)*NGLLX,cudaMemcpyDeviceToHost),7702);
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_elastic_top,mp->d_b_absorb_elastic_top,
                                       2*mp->d_nspec_top*sizeof(realw)*NGLLX,cudaMemcpyDeviceToHost),7703);
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_elastic_bottom,mp->d_b_absorb_elastic_bottom,
                                       2*mp->d_nspec_bottom*sizeof(realw)*NGLLX,cudaMemcpyDeviceToHost),7704);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after compute_stacey_elastic after cudamemcpy");
#endif
}

