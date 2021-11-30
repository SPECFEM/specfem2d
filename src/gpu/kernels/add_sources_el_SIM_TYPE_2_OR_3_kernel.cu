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


__global__ void add_sources_el_SIM_TYPE_2_OR_3_kernel(realw* accel,
                                                      realw* source_adjoint,
                                                      realw* xir_store,
                                                      realw* gammar_store,
                                                      int* d_ibool,
                                                      int* ispec_is_elastic,
                                                      int* ispec_selected_rec_loc,
                                                      int it,
                                                      int nadj_rec_local,
                                                      int NSTEP) {

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  // when nrec > 65535, but mod(nspec_top,2) > 0, we end up with an extra block.

  if (irec_local < nadj_rec_local) {

    int ispec = ispec_selected_rec_loc[irec_local] - 1;

    if (ispec_is_elastic[ispec]) {
      int i = threadIdx.x;
      int j = threadIdx.y;
      int iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

      realw  xir = xir_store[INDEX2(nadj_rec_local,irec_local,i)];
      realw  gammar = gammar_store[INDEX2(nadj_rec_local,irec_local,j)];

      realw  source_adjx = source_adjoint[INDEX3(nadj_rec_local,NSTEP,irec_local,it,0)];
      realw  source_adjz = source_adjoint[INDEX3(nadj_rec_local,NSTEP,irec_local,it,1)];

      // atomic operations are absolutely necessary for correctness!
      atomicAdd(&accel[2*iglob],source_adjx * gammar * xir);
      atomicAdd(&accel[2*iglob+1], source_adjz * gammar * xir);
    } // ispec_is_elastic
  }
}


