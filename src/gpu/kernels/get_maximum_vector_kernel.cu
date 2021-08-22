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


__global__ void get_maximum_vector_kernel(realw* array, int size, realw* d_max){

  // reduction example:
  __shared__ realw sdata[BLOCKSIZE_TRANSFER] ;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int bx = blockIdx.y*gridDim.x+blockIdx.x;
  unsigned int i = tid + bx*blockDim.x;

  // loads values into shared memory: assume array is a vector array
  sdata[tid] = (i < size) ? (array[i*2]*array[i*2] + array[i*2+1]*array[i*2+1]) : 0.0f ;

  __syncthreads();

  // do reduction in shared mem
  for(unsigned int s=blockDim.x/2; s>0; s>>=1)
  {
    if (tid < s){
      // summation:
      //sdata[tid] += sdata[tid + s];
      // maximum:
      if (sdata[tid] < sdata[tid + s] ) sdata[tid] = sdata[tid + s];
    }
    __syncthreads();
  }

  // write result for this block to global mem
  if (tid == 0) d_max[bx] = sdata[0];

}

