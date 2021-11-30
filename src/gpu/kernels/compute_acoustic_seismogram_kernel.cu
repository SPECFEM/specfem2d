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


__global__ void compute_acoustic_seismogram_kernel(int nrec_local,
                                                   realw* displ,
                                                   realw* pressure,
                                                   int* d_ibool,
                                                   realw* hxir, realw* hgammar,
                                                   realw* seismograms,
                                                   int* ispec_selected_rec_loc,
                                                   int* ispec_is_elastic,
                                                   int* ispec_is_acoustic,
                                                   realw* d_kappav,
                                                   realw* d_muv,
                                                   realw* d_hprime_xx,
                                                   realw* d_xix,realw* d_xiz,
                                                   realw* d_gammax,realw* d_gammaz,
                                                   int* ispec_is_anisotropic,
                                                   realw* d_c11store,realw* d_c12store,realw* d_c13store,
                                                   realw* d_c15store,realw* d_c23store,realw* d_c25store,
                                                   realw* d_c33store,realw* d_c35store,
                                                   int it,
                                                   int NSTEP){

  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  __shared__ realw sh_dxd[NGLL2_PADDED];
  __shared__ realw sh_tempx[NGLL2];
  __shared__ realw sh_tempz[NGLL2];

  if (irec_local < nrec_local) {
    // initializes
    sh_dxd[tx] = 0.0f;

    // receiver element
    int ispec = ispec_selected_rec_loc[irec_local]-1;

    // acoustic domains
    if (ispec_is_acoustic[ispec]) {
      if (tx < NGLL2) {
        int iglob = d_ibool[tx+NGLL2_PADDED*ispec]-1;

        realw hlagrange = hxir[irec_local + nrec_local*I]*hgammar[irec_local + nrec_local*J];

        // Signe moins car pression = -potential_dot_dot
        sh_dxd[tx] = - hlagrange * pressure[iglob];
      }
    }

    // elastic domains
    if (ispec_is_elastic[ispec]) {
      // see pressure comments in compute_pressure.f90
      // this routines limits the pressure computations to: non-Axisem, non-anisotropic, non-attenuation case
      // todo for the future...

      // loads displ into shared memory
      if (tx < NGLL2) {
        int iglob = d_ibool[tx+NGLL2_PADDED*ispec] - 1;
        sh_tempx[tx] = displ[iglob*2];
        sh_tempz[tx] = displ[iglob*2+1];
      }
      // synchronizes threads
      __syncthreads();

      // derivative along x (dux_xi,duz_xi)
      realw temp1x = 0.0f;
      realw temp1z = 0.0f;
      for(int l=0; l<NGLLX;l++){
        realw hp1 = d_hprime_xx[l*NGLLX+I];
        temp1x += sh_tempx[J*NGLLX+l] * hp1;
        temp1z += sh_tempz[J*NGLLX+l] * hp1;
      }

      // derivative along z (dux_gamma,duz_gamma)
      realw temp3x = 0.0f;
      realw temp3z = 0.0f;
      for(int l=0; l<NGLLX;l++){
        // assumes hprime_xx == hprime_zz
        realw hp3 = d_hprime_xx[l*NGLLX+J];
        temp3x += sh_tempx[l*NGLLX+I] * hp3;
        temp3z += sh_tempz[l*NGLLX+I] * hp3;
      }

      int ij_ispec_padded = tx + NGLL2_PADDED*ispec;

      realw xixl = d_xix[ij_ispec_padded];
      realw xizl = d_xiz[ij_ispec_padded];
      realw gammaxl = d_gammax[ij_ispec_padded];
      realw gammazl = d_gammaz[ij_ispec_padded];

      // derivatives of displ field on GLL points
      realw duxdxl = (temp1x*xixl + temp3x*gammaxl);
      realw duzdzl = (temp1z*xizl + temp3z*gammazl);

      // isotropic
      realw kappal = d_kappav[ij_ispec_padded];
      realw mul = d_muv[ij_ispec_padded];

      realw lambdal = kappal - mul;
      realw lambdalplus2mul = kappal + mul;

      // compute diagonal components of the stress tensor
      // non-anisotropic, non-attenuation
      // pressure for P_SV case only
      realw sigma_xx = lambdalplus2mul * duxdxl + lambdal * duzdzl;
      realw sigma_zz = lambdalplus2mul * duzdzl + lambdal * duxdxl;
      // sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
      realw sigma_yy = lambdal * (duxdxl + duzdzl);

      // sigma_xz = mul * (duzdxl+duxdzl); not needed for pressure

      // anisotropic case
      if (ispec_is_anisotropic[ispec]){
        // anisotropy
        realw c11 = d_c11store[ij_ispec_padded];
        realw c12 = d_c12store[ij_ispec_padded];
        realw c13 = d_c13store[ij_ispec_padded];
        realw c15 = d_c15store[ij_ispec_padded];
        realw c23 = d_c23store[ij_ispec_padded];
        realw c25 = d_c25store[ij_ispec_padded];
        realw c33 = d_c33store[ij_ispec_padded];
        realw c35 = d_c35store[ij_ispec_padded];

        realw duxdzl = (temp1x*xizl + temp3x*gammazl);
        realw duzdxl = (temp1z*xixl + temp3z*gammaxl);

        // compute the three components of the stress tensor sigma (full anisotropy)
        // implement anisotropy in 2D
        sigma_xx = c11 * duxdxl + c13 * duzdzl + c15 * (duzdxl + duxdzl);
        sigma_zz = c13 * duxdxl + c33 * duzdzl + c35 * (duzdxl + duxdzl);
        // sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
        sigma_yy = c12 * duxdxl + c23 * duzdzl + c25 * (duzdxl + duxdzl);

        if (c12 < 1.e-7 || c23 < 1.e-7){
          // cannot compute pressure for an anisotropic material if c12 or c23 are zero
          sigma_xx = 0.0f; sigma_zz = 0.0f; sigma_yy = 0.0f;
        }
      }

      realw hlagrange = hxir[irec_local + nrec_local*I] * hgammar[irec_local + nrec_local*J];

      // stores pressure
      sh_dxd[tx] = - hlagrange * (sigma_xx + sigma_yy + sigma_zz) / 3.0f;
    }

    // synchronizes threads
    __syncthreads();

    // reduction
    for (unsigned int s=1; s<NGLL2_PADDED ; s *= 2) {
      if (tx % (2*s) == 0) sh_dxd[tx] += sh_dxd[tx + s];
      __syncthreads();
    }

    if (tx == 0) {
      seismograms[irec_local*NSTEP + it ] = sh_dxd[0];
    }
    if (tx == 1) {
      seismograms[irec_local*NSTEP + it + nrec_local*NSTEP] = 0;
    }
  }
}

