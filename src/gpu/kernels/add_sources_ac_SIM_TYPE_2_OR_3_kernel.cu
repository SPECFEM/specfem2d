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


__global__ void add_sources_ac_SIM_TYPE_2_OR_3_kernel(realw* potential_dot_dot_acoustic,
                                                      realw* source_adjoint,
                                                      realw* xir_store,
                                                      realw* gammar_store,
                                                      int* d_ibool,
                                                      int* ispec_is_acoustic,
                                                      int* ispec_selected_rec_loc,
                                                      int it,
                                                      int nadj_rec_local,
                                                      //realw* kappastore,
                                                      int NSTEP ) {

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  // because of grid shape, irec_local can be too big
  if (irec_local < nadj_rec_local) {

    int ispec = ispec_selected_rec_loc[irec_local] - 1;

    if (ispec_is_acoustic[ispec]) {
      int i = threadIdx.x;
      int j = threadIdx.y;

      int iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

      //realw  kappal = kappastore[INDEX3(NGLLX,NGLLX,i,j,ispec)];
      realw  xir = xir_store[INDEX2(nadj_rec_local,irec_local,i)];
      realw  gammar = gammar_store[INDEX2(nadj_rec_local,irec_local,j)];
      realw  source_adj = source_adjoint[INDEX3(nadj_rec_local,NSTEP,irec_local,it,0)];

      // adjoint source of Peter et al. (A8):
      //   f^adj = - sum_i \partial_t^2 (p^syn - p^obs)(T-t) \delta(x - x_i)
      // note that using the adjoint source derived from the optimization problem, there is no 1/kappa term applied
      // to the adjoint source. the negative sign also is part of the construction of the adjoint source.
      //
      // since we don't know which formulation of adjoint source is used for the input, we add the adjoint source as is,
      // without 1/kappa factor, and with a positive sign.
      realw stf = source_adj * gammar * xir;
      atomicAdd(&potential_dot_dot_acoustic[iglob],stf);
    }
  }
}

