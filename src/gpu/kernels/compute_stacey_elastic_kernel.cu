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


__global__ void compute_stacey_elastic_kernel(realw* veloc,
                                              realw* accel,
                                              const int* abs_boundary_ispec,
                                              const int* abs_boundary_ijk,
                                              const realw* abs_boundary_normal,
                                              const realw* abs_boundary_jacobian1Dw,
                                              const int* d_ibool,
                                              const realw* rho_vp,
                                              const realw* rho_vs,
                                              const int* ispec_is_elastic,
                                              const int simulation_type,
                                              const int p_sv,
                                              const int SAVE_FORWARD,
                                              const int num_abs_boundary_faces,
                                              realw* b_absorb_elastic_left,
                                              realw* b_absorb_elastic_right,
                                              realw* b_absorb_elastic_top,
                                              realw* b_absorb_elastic_bottom,
                                              const int* ib_left,
                                              const int* ib_right,
                                              const int* ib_top,
                                              const int* ib_bottom,
                                              const int* edge_abs) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,iglob,ispec,num_local;
  realw vx,vz,vn;
  realw nx,nz;
  realw rho_vp_temp,rho_vs_temp;
  realw tx,tz;
  realw jacobianw;
  realw absorblx,absorblz;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if (iface < num_abs_boundary_faces){

  //if (igll < NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface] - 1;

    if (ispec_is_elastic[ispec]) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLLX,0,igll,iface)] - 1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLLX,1,igll,iface)] - 1;

      //daniel todo: check if we can simplify.
      //       in fortran routine, we set i == NGLLX+1 or j == NGLLX+1
      //       to indicate points which duplicate contributions and can be left out
      //
      //check if the point must be computed
      if (i==NGLLX || j==NGLLX) return;

      iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

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

      if (p_sv){
        // P_SV case
        tx = rho_vp_temp*vn*nx + rho_vs_temp*(vx-vn*nx);
        tz = rho_vp_temp*vn*nz + rho_vs_temp*(vz-vn*nz);
      }else{
        // SH-case
        tx = rho_vs_temp*vx;  // would be vy = veloc_elastic(1,iglob); ty = rho_vs*vy for CPU-case
        tz = 0.f;             // second component not needed
      }

      jacobianw = abs_boundary_jacobian1Dw[INDEX2(NGLLX,igll,iface)];

      absorblx = tx * jacobianw;
      absorblz = tz * jacobianw;

      atomicAdd(&accel[iglob*2],-absorblx);
      atomicAdd(&accel[iglob*2+1],-absorblz);

      if (SAVE_FORWARD && simulation_type == 1) {
        if (edge_abs[iface] == 1) {
          num_local = ib_bottom[iface] - 1;
          b_absorb_elastic_bottom[INDEX3(NDIM,NGLLX,0,igll,num_local)] = absorblx;
          b_absorb_elastic_bottom[INDEX3(NDIM,NGLLX,1,igll,num_local)] = absorblz;

        }else if (edge_abs[iface] == 2) {
          num_local = ib_right[iface] - 1;
          b_absorb_elastic_right[INDEX3(NDIM,NGLLX,0,igll,num_local)] = absorblx;
          b_absorb_elastic_right[INDEX3(NDIM,NGLLX,1,igll,num_local)] = absorblz;

        }else if (edge_abs[iface] == 3) {
          num_local = ib_top[iface] - 1;
          b_absorb_elastic_top[INDEX3(NDIM,NGLLX,0,igll,num_local)] = absorblx;
          b_absorb_elastic_top[INDEX3(NDIM,NGLLX,1,igll,num_local)] = absorblz;

        }else if (edge_abs[iface] == 4) {
          num_local = ib_left[iface] - 1;
          b_absorb_elastic_left[INDEX3(NDIM,NGLLX,0,igll,num_local)] = absorblx;
          b_absorb_elastic_left[INDEX3(NDIM,NGLLX,1,igll,num_local)] = absorblz;
        }
      }
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
                                                   int* d_edge_abs) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,iglob,ispec,num_local;

  if (iface < num_abs_boundary_faces){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface] - 1;

    if (ispec_is_elastic[ispec]) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLLX,0,igll,iface)] - 1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLLX,1,igll,iface)] - 1;

      //daniel todo: check if we can simplify.
      //       in fortran routine, we set i == NGLLX+1 or j == NGLLX+1
      //       to indicate points which duplicate contributions and can be left out
      //
      //check if the point must be computed
      if (i==NGLLX || j==NGLLX) return;

      iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

      if (d_edge_abs[iface] == 1){
        num_local= ib_bottom[iface]-1;
        atomicAdd(&b_accel[iglob*2 ], -b_absorb_elastic_bottom[INDEX3(NDIM,NGLLX,0,igll,num_local)]);
        atomicAdd(&b_accel[iglob*2+1 ], -b_absorb_elastic_bottom[INDEX3(NDIM,NGLLX,1,igll,num_local)]);

      } else if (d_edge_abs[iface] == 2){
        num_local= ib_right[iface]-1;
        atomicAdd(&b_accel[iglob*2 ], -b_absorb_elastic_right[INDEX3(NDIM,NGLLX,0,igll,num_local)]);
        atomicAdd(&b_accel[iglob*2+1 ], -b_absorb_elastic_right[INDEX3(NDIM,NGLLX,1,igll,num_local)]);

      } else if (d_edge_abs[iface] == 3){
        num_local= ib_top[iface]-1;
        atomicAdd(&b_accel[iglob*2 ], -b_absorb_elastic_top[INDEX3(NDIM,NGLLX,0,igll,num_local)]);
        atomicAdd(&b_accel[iglob*2+1 ], -b_absorb_elastic_top[INDEX3(NDIM,NGLLX,1,igll,num_local)]);

      } else if (d_edge_abs[iface] == 4){
        num_local= ib_left[iface]-1;
        atomicAdd(&b_accel[iglob*2 ], -b_absorb_elastic_left[INDEX3(NDIM,NGLLX,0,igll,num_local)]);
        atomicAdd(&b_accel[iglob*2+1 ], -b_absorb_elastic_left[INDEX3(NDIM,NGLLX,1,igll,num_local)]);
      }
    }
  } // num_abs_boundary_faces
}

