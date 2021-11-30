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
                                           realw* displ,
                                           realw* accel,
                                           realw* b_displ,
                                           realw* rho_kl,
                                           realw* mu_kl,
                                           realw* kappa_kl,
                                           int NSPEC_AB,
                                           realw* d_hprime_xx,
                                           realw* d_xix,realw* d_xiz,
                                           realw* d_gammax,realw* d_gammaz) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ij = threadIdx.x;

  // local and global indices
  int ij_ispec = ij + NGLL2*ispec;
  int ij_ispec_padded = ij + NGLL2_PADDED*ispec;
  int iglob;

  int offset1,offset3;
  realw hp1,hp3;

  realw xixl,xizl;
  realw gammaxl,gammazl;

  realw duxdxl,duxdzl;
  realw duzdxl,duzdzl;

  realw b_duxdxl,b_duxdzl;
  realw b_duzdxl,b_duzdzl;

  realw dsxx,dszz,dsxz;
  realw b_dsxx,b_dszz,b_dsxz;

  // shared memory between all threads within this block
  __shared__ realw field_displ[2*NGLL2];
  __shared__ realw field_b_displ[2*NGLL2];

  int active = 0;

  // computing the gradient of displ and accel

  // handles case when there is 1 extra block (due to rectangular grid)
  if (ispec < NSPEC_AB) {
    // elastic elements only
    if (ispec_is_elastic[ispec]) {
      active = 1;

      // global index
      iglob = d_ibool[ij_ispec_padded] - 1;

      // copy field values
      field_displ[2*ij] = displ[2*iglob];
      field_displ[2*ij+1] = displ[2*iglob+1];

      field_b_displ[2*ij] = b_displ[2*iglob];
      field_b_displ[2*ij+1] = b_displ[2*iglob+1];
    }
  }

  // synchronizes threads
  __syncthreads();

  if (active) {
    int J = (ij/NGLLX);
    int I = (ij-J*NGLLX);

    // derivative along xi
    realw tempx1l = 0.f;
    realw tempz1l = 0.f;

    realw b_tempx1l = 0.f;
    realw b_tempz1l = 0.f;

    for(int l=0; l<NGLLX;l++){
      hp1 = d_hprime_xx[l*NGLLX+I];
      offset1 = J*NGLLX+l;

      tempx1l += field_displ[2*offset1] * hp1;
      tempz1l += field_displ[2*offset1+1] * hp1;

      b_tempx1l += field_b_displ[2*offset1] * hp1;
      b_tempz1l += field_b_displ[2*offset1+1] * hp1;
    }

    // derivative along gamma
    realw tempx3l = 0.f;
    realw tempz3l = 0.f;

    realw b_tempx3l = 0.f;
    realw b_tempz3l = 0.f;

    for(int l=0; l<NGLLX;l++){
      // assumes hprime_xx == hprime_zz
      hp3 = d_hprime_xx[l*NGLLX+J];
      offset3 = l*NGLLX+I;

      tempx3l += field_displ[2*offset3] * hp3;
      tempz3l += field_displ[2*offset3+1] * hp3;

      b_tempx3l += field_b_displ[2*offset3] * hp3;
      b_tempz3l += field_b_displ[2*offset3+1] * hp3;
    }

    xixl = d_xix[ij_ispec_padded];
    xizl = d_xiz[ij_ispec_padded];
    gammaxl = d_gammax[ij_ispec_padded];
    gammazl = d_gammaz[ij_ispec_padded];

    // derivatives of ux and uz with respect to x and z
    duxdxl = xixl*tempx1l + gammaxl*tempx3l;
    duxdzl = xizl*tempx1l + gammazl*tempx3l;

    duzdxl = xixl*tempz1l + gammaxl*tempz3l;
    duzdzl = xizl*tempz1l + gammazl*tempz3l;

    b_duxdxl = xixl*b_tempx1l + gammaxl*b_tempx3l;
    b_duxdzl = xizl*b_tempx1l + gammazl*b_tempx3l;

    b_duzdxl = xixl*b_tempz1l + gammaxl*b_tempz3l;
    b_duzdzl = xizl*b_tempz1l + gammazl*b_tempz3l;

    // strain components
    dsxx = duxdxl;
    b_dsxx = b_duxdxl;

    dszz = duzdzl;
    b_dszz = b_duzdzl;

    if (p_sv){
      // P_SV case
      dsxz = 0.5f * (duzdxl + duxdzl);
      b_dsxz = 0.5f * (b_duzdxl + b_duxdzl);
    }else{
      // SH-case
      dsxz = duxdzl;
      b_dsxz = b_duxdzl;
    }

    // isotropic kernels:
    // density kernel
    rho_kl[ij_ispec] += (accel[2*iglob]*b_displ[2*iglob] + accel[2*iglob+1]*b_displ[2*iglob+1]);

    if (p_sv){
      // P_SV case
      realw prod = (dsxx + dszz) * (b_dsxx + b_dszz);

      // shear modulus kernel
      mu_kl[ij_ispec] += dsxx * b_dsxx + dszz * b_dszz +
                          2.0f * dsxz * b_dsxz - prod/3.0f;

      // bulk modulus kernel
      kappa_kl[ij_ispec] += prod;
    }else{
      // SH-case (membrane) waves
      mu_kl[ij_ispec] += 0.5f * (dsxx * b_dsxx + dsxz * b_dsxz);
      // zero kappa kernel
      //kappa_kl[ij_ispec] = 0.f;
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

