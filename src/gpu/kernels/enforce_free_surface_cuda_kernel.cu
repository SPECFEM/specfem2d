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


__global__ void enforce_free_surface_cuda_kernel(realw_p potential_acoustic,
                                                 realw_p potential_dot_acoustic,
                                                 realw_p potential_dot_dot_acoustic,
                                                 const int num_free_surface_faces,
                                                 const int* free_surface_ispec,
                                                 const int* free_surface_ijk,
                                                 const int* d_ibool,
                                                 const int* ispec_is_acoustic) {
  // gets spectral element face id
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  // for all faces on free surface
  if (iface < num_free_surface_faces) {

    int ispec = free_surface_ispec[iface]-1;

    // checks if element is in acoustic domain
    if (ispec_is_acoustic[ispec]) {

      // gets global point index
      int igll = threadIdx.x + threadIdx.y*blockDim.x;

      int i = free_surface_ijk[INDEX3(NDIM,NGLLX,0,igll,iface)] - 1; // (1,igll,iface)
      int j = free_surface_ijk[INDEX3(NDIM,NGLLX,1,igll,iface)] - 1;

      int iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

      // sets potentials to zero at free surface
      potential_acoustic[iglob] = 0.f;
      potential_dot_acoustic[iglob] = 0.f;
      potential_dot_dot_acoustic[iglob] = 0.f;
    }
  }
}

