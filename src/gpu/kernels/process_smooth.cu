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


__global__ void process_smooth(realw_const_p xstore_me,
                               realw_const_p zstore_me,
                               realw_const_p xstore_other,
                               realw_const_p zstore_other,
                               realw_const_p data_other,
                               const realw sigma_h2_inv,
                               const realw sigma_v2_inv,
                               const int iker,
                               const int nspec_me,
                               const int nspec_other,
                               const realw v_criterion,
                               const realw h_criterion,
                               realw_const_p jacobian,
                               realw_p sum_data_smooth,
                               realw_p normalisation,
                               realw_const_p wgll_sq){

  int ispec = blockIdx.x + gridDim.x*blockIdx.y;
  int igll = threadIdx.x;
  int gll_other;
  realw x_me, z_me, x_other, z_other, coef, normalisation_slice;
  realw dat;
  __shared__ int sh_test[NGLL2];
  __shared__ realw sh_x_other[NGLL2];
  __shared__ realw sh_z_other[NGLL2];
  __shared__ realw sh_jacobian[NGLL2];
  __shared__ realw sh_wgll_sq[NGLL2];
  __shared__ realw sh_data[NGLL2];

  int n_loop = nspec_other/NGLL2 + 1;
  x_me = xstore_me[NGLL2*ispec + igll ];
  z_me = zstore_me[NGLL2*ispec + igll ];
  sh_wgll_sq[igll]=wgll_sq[igll];
  __syncthreads();

  dat=0.f;
  normalisation_slice=0.f;
  //We test 32 spectral elements at a time
  for (int i=0;i<n_loop;i++)
  {
     __syncthreads();
    if (NGLL2*i + threadIdx.x < nspec_other){
      x_other = xstore_other[i*NGLL2*NGLL2 + threadIdx.x*NGLL2 ];
      z_other = zstore_other[i*NGLL2*NGLL2 + threadIdx.x*NGLL2 ];
    }
    sh_test[threadIdx.x] = ( NGLL2*i + threadIdx.x >= nspec_other || (x_me-x_other)*(x_me-x_other) > h_criterion || (z_me-z_other)*(z_me-z_other) > v_criterion ) ? 1 : 0 ;
    __syncthreads();

    //loop over each spectral element tested
    for (int k=0;k<NGLL2;k++)
    {
      if (sh_test[k]) continue ;

      //Load data from other slice to shared memory
      sh_x_other[igll] = xstore_other[i*NGLL2*NGLL2 + k*NGLL2 + igll ];
      sh_z_other[igll] = zstore_other[i*NGLL2*NGLL2 + k*NGLL2 + igll ];
      sh_data[igll] = data_other[i*NGLL2*NGLL2 + k*NGLL2 + igll ];
      sh_jacobian[igll] = jacobian[i*NGLL2*NGLL2 + k*NGLL2 + igll ];
      __syncthreads();

      for (int j=0;j<NGLL2;j++){

        gll_other = (igll + j) % NGLL2;

        x_other = sh_x_other[gll_other];
        z_other = sh_z_other[gll_other];
        coef = expf(- sigma_h2_inv*(x_me-x_other)*(x_me-x_other) - sigma_v2_inv*(z_me-z_other)*(z_me-z_other))*sh_jacobian[gll_other]*sh_wgll_sq[gll_other];
        normalisation_slice = normalisation_slice + coef;
        dat += sh_data[gll_other]*coef;
      } //loop on each gll_other
    } //loop on each spec_other tested
  } //loop on each serie of 32 spec_other

  sum_data_smooth[NGLL2*nspec_me*iker+NGLL2*ispec + igll] += dat;
  normalisation[NGLL2*ispec + igll] += normalisation_slice;
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void normalize_data(realw_p data_smooth,
                               realw_const_p normalisation,
                               int nker,
                               int nspec_me){
  int ispec = blockIdx.x + gridDim.x*blockIdx.y;
  realw norm = normalisation[NGLL2*ispec + threadIdx.x];
  for (int j=0;j<nker;j++) data_smooth[NGLL2*nspec_me*j + NGLL2*ispec + threadIdx.x] /= norm/nker;
}

