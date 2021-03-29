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


__global__ void compute_kernels_cudakernel(const int* ispec_is_elastic,
                                           const int p_sv,
                                           const int* d_ibool,
                                           realw* accel,
                                           realw* b_displ,
                                           realw* rho_kl,
                                           realw* mu_kl,
                                           realw* kappa_kl,
                                           int NSPEC_AB,
                                           realw* dsxx,
                                           realw* dsxz,
                                           realw* dszz,
                                           realw* b_dsxx,
                                           realw* b_dsxz,
                                           realw* b_dszz) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ij = threadIdx.x;

  // handles case when there is 1 extra block (due to rectangular grid)
  if (ispec < NSPEC_AB) {

    // elastic elements only
    if (ispec_is_elastic[ispec]) {
      int iglob = d_ibool[ij + NGLL2_PADDED*ispec] - 1 ;
      int ij_ispec = ij + NGLL2*ispec;

      realw prod = (dsxx[iglob] + dszz[iglob]) * (b_dsxx[iglob] + b_dszz[iglob]);

      // isotropic kernels:
      // density kernel
      rho_kl[ij_ispec] += (accel[2*iglob]*b_displ[2*iglob] + accel[2*iglob+1]*b_displ[2*iglob+1]);

      if (p_sv){
        // P_SV case
        // shear modulus kernel
        mu_kl[ij_ispec] += dsxx[iglob] * b_dsxx[iglob] + dszz[iglob] * b_dszz[iglob] +
                            2.0f * dsxz[iglob] * b_dsxz[iglob] - prod/3.0f;

        // bulk modulus kernel
        kappa_kl[ij_ispec] += prod;
      }else{
        // SH-case (membrane) waves
        mu_kl[ij_ispec] += dsxx[iglob] * b_dsxx[iglob] + dsxz[iglob] * b_dsxz[iglob]; // dux_dxl * b_dux_dxl + dux_dzl * b_dux_dzl
        // zero kappa kernel
        //kappa_kl[ij_ispec] = 0.f;
      }
    }
  }
}


//  !!!!!!!!!!  BEWARE this kernel needs to be adapted, it only reflect the 3D case for now
/*
__global__ void compute_kernels_ani_cudakernel(int* ispec_is_elastic,
                                               int* d_ibool,
                                               realw* accel,
                                               realw* b_displ,
                                               realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                               realw* epsilondev_xz,realw* epsilondev_yz,
                                               realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,
                                               realw* b_epsilondev_xz,realw* b_epsilondev_yz,
                                               realw* rho_kl,
                                               realw deltat,
                                               realw* cijkl_kl,
                                               realw* epsilon_trace_over_3,
                                               realw* b_epsilon_trace_over_3,
                                               int NSPEC_AB) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;
  int ijk_ispec = ijk + NGLL2*ispec;
  int ijk21_ispec = ijk + 9*NGLL2*ispec;

  realw prod[21];
  realw eps[6];
  realw b_eps[6];
  realw epsdev[6];
  realw b_epsdev[6];
  realw eps_trace_over_3,b_eps_trace_over_3;
  int i,j;

  // handles case when there is 1 extra block (due to rectangular grid)
  if (ispec < NSPEC_AB) {

    // elastic elements only
    if (ispec_is_elastic[ispec]) {
      int iglob = d_ibool[ijk + NGLL2_PADDED*ispec] - 1;

      // anisotropic kernels:
      // density kernel
      rho_kl[ijk_ispec] += deltat * (accel[3*iglob]*b_displ[3*iglob]+
                                     accel[3*iglob+1]*b_displ[3*iglob+1]+
                                     accel[3*iglob+2]*b_displ[3*iglob+2]);


      // anisotropic kernel
      epsdev[0] = epsilondev_xx[ijk_ispec];
      epsdev[1] = epsilondev_yy[ijk_ispec];
      epsdev[2] = epsilondev_xy[ijk_ispec];
      epsdev[3] = epsilondev_xz[ijk_ispec];
      epsdev[4] = epsilondev_yz[ijk_ispec];

      b_epsdev[0] = b_epsilondev_xx[ijk_ispec];
      b_epsdev[1] = b_epsilondev_yy[ijk_ispec];
      b_epsdev[2] = b_epsilondev_xy[ijk_ispec];
      b_epsdev[3] = b_epsilondev_xz[ijk_ispec];
      b_epsdev[4] = b_epsilondev_yz[ijk_ispec];

      eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec];
      b_eps_trace_over_3 = b_epsilon_trace_over_3[ijk_ispec];

      //! Building of the local matrix of the strain tensor
      //! for the adjoint field and the regular backward field
      //!eps11 et eps22
      eps[0] = epsdev[0] + eps_trace_over_3;
      eps[1] = epsdev[1] + eps_trace_over_3;
      //!eps33
      eps[2] = -(eps[0]+eps[1])+3*eps_trace_over_3;
      //!eps23
      eps[3] = epsdev[4];
      //!eps13
      eps[4] = epsdev[3];
      //!eps12
      eps[5] = epsdev[2];

      // backward arrays
      b_eps[0] = b_epsdev[0] + b_eps_trace_over_3;
      b_eps[1] = b_epsdev[1] + b_eps_trace_over_3;
      b_eps[2] = -(b_eps[0]+b_eps[1])+3*b_eps_trace_over_3;
      b_eps[3] = b_epsdev[4];
      b_eps[4] = b_epsdev[3];
      b_eps[5] = b_epsdev[2];

      //! Computing the 21 strain products without assuming eps(i)*b_eps(j) = eps(j)*b_eps(i)
      int p = 0;
      for( i=0; i<6; i++){
        for( j=i; j<6; j++){
          prod[p] = eps[i] * b_eps[j];
          if (j > i) {
            prod[p] = prod[p] + eps[j]*b_eps[i];
            if (j > 2 && i < 3) { prod[p] = prod[p]*2; }
          }
          if (i > 2) { prod[p] = prod[p]*4; }
          p++;
        }
      }

      // all 21 anisotropic coefficients
      for( i=0; i<21; i++){
        cijkl_kl[i+ijk21_ispec] += deltat * prod[i];
      }

    }
  }
}
*/

