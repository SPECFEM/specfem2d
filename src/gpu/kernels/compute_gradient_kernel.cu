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


__device__ void compute_gradient_kernel(int ij,
                                        int ispec,
                                        realw* scalar_field,
                                        realw* vector_field_element,
                                        realw* d_hprime_xx,
                                        realw* d_xix,realw* d_xiz,
                                        realw* d_gammax,realw* d_gammaz,
                                        realw rhol) {

  realw temp1l,temp3l;
  realw hp1,hp3;
  realw xixl,xizl,gammaxl,gammazl;
  realw rho_invl;
  int l,offset,offset1,offset3;

  const int NGLL2_ALIGN = NGLL2_PADDED;

  int J = (ij/NGLLX);
  int I = (ij-J*NGLLX);

  // derivative along x
  temp1l = 0.f;
  for( l=0; l<NGLLX;l++){
    hp1 = d_hprime_xx[l*NGLLX+I];
    offset1 = J*NGLLX+l;
    temp1l += scalar_field[offset1]*hp1;
  }

  // derivative along z
  temp3l = 0.f;
  for( l=0; l<NGLLX;l++){
    // assumes hprime_xx == hprime_zz
    hp3 = d_hprime_xx[l*NGLLX+J];
    offset3 = l*NGLLX+I;
    temp3l += scalar_field[offset3]*hp3;
  }

  offset = ispec*NGLL2_ALIGN + ij;

  xixl = d_xix[offset];
  xizl = d_xiz[offset];
  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  rho_invl = 1.0f / rhol;

  // derivatives of acoustic scalar potential field on GLL points
  vector_field_element[0] = (temp1l*xixl + temp3l*gammaxl) * rho_invl;
  vector_field_element[1] = (temp1l*xizl + temp3l*gammazl) * rho_invl;
}

