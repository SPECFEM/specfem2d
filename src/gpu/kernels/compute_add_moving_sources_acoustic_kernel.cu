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
      stf = source_time_function_moving[INDEX2(nsources, isource, it)] / kappal;

      accel = sourcearrays_moving[INDEX5(NDIM, NGLLX, NGLLX, nsources, 0, i, j, isource, it)] * stf;

      atomicAdd(&potential_dot_dot_acoustic[iglob], accel);
    }
  }
}

