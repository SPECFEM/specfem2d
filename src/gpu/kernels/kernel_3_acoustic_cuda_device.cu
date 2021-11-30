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


__global__ void kernel_3_acoustic_cuda_device(realw* potential_dot_dot_acoustic,
                                              realw* b_potential_dot_dot_acoustic,
                                              realw* potential_dot_acoustic,
                                              realw* b_potential_dot_acoustic,
                                              int size,
                                              int compute_wavefield_1,
                                              int compute_wavefield_2,
                                              realw deltatover2,
                                              realw b_deltatover2,
                                              realw* rmass_acoustic) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  realw p_dot_dot;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if (id < size) {
    // multiplies pressure with the inverse of the mass matrix
    realw rmass = rmass_acoustic[id];

    if (compute_wavefield_1){
      p_dot_dot = potential_dot_dot_acoustic[id]*rmass;
      potential_dot_dot_acoustic[id] = p_dot_dot;
      // corrector:
      // updates the chi_dot term which requires chi_dot_dot(t+delta)
      potential_dot_acoustic[id] += deltatover2*p_dot_dot;
    }

    if (compute_wavefield_2){
      p_dot_dot = b_potential_dot_dot_acoustic[id]*rmass;
      b_potential_dot_dot_acoustic[id] = p_dot_dot;
      // corrector:
      // updates the chi_dot term which requires chi_dot_dot(t+delta)
      b_potential_dot_acoustic[id] += b_deltatover2*p_dot_dot;
    }

  } // id<size
}

