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


__global__ void pml_boundary_acoustic_cuda_kernel(realw_p potential_acoustic,
                                                  realw_p potential_dot_acoustic,
                                                  realw_p potential_dot_dot_acoustic,
                                                  int pml_nglob_abs_acoustic,
                                                  const int* d_pml_abs_points_acoustic) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // for all PML points on absorbing boundary
  if (id < pml_nglob_abs_acoustic) {

    int iglob = d_pml_abs_points_acoustic[id] - 1;

    // sets potentials to zero at absorbing boundary
    potential_acoustic[iglob] = 0.f;
    potential_dot_acoustic[iglob] = 0.f;
    potential_dot_dot_acoustic[iglob] = 0.f;
  }
}


