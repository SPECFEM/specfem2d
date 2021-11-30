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
  realw stf,accel_x,accel_z;

  // when NSOURCES > 65535, but mod(nspec_top,2) > 0, we end up with an extra block.

  if (isource < nsources_local) {
    ispec = ispec_selected_source[isource] - 1;

    if (ispec_is_elastic[ispec]) {
      iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

      stf = d_source_time_function[INDEX2(nsources_local,isource,it)];
      accel_x = sourcearrays[INDEX4(NDIM,NGLLX,NGLLX, 0,i,j,isource)] * stf;
      accel_z = sourcearrays[INDEX4(NDIM,NGLLX,NGLLX, 1,i,j,isource)] * stf;

      atomicAdd(&accel[iglob*2],accel_x);
      atomicAdd(&accel[iglob*2+1],accel_z);
    }
  }
}

