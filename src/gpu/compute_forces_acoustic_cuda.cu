/*
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
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

#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2 - acoustic compute forces kernel

/* ----------------------------------------------------------------------------------------------- */

void Kernel_2_acoustic(int nb_blocks_to_compute, Mesh* mp, int d_iphase,
                       int* d_ibool,
                       realw* d_xix,realw* d_xiz,
                       realw* d_gammax,realw* d_gammaz,
                       realw* d_rhostore,
                       int ATTENUATION_VISCOACOUSTIC,
                       int compute_wavefield_1,
                       int compute_wavefield_2) {

TRACE("Kernel_2_acoustic");

  // if the grid can handle the number of blocks, we let it be 1D
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(nb_blocks_to_compute,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  cudaEvent_t start, stop;
  if (CUDA_TIMING) {
    start_timing_cuda(&start,&stop);
  }

  // forward and/or backward fields
  int nb_field;
  if (compute_wavefield_1 && compute_wavefield_2){
    nb_field = 2;  // both fields
  }else{
    nb_field = 1;  // single field only (either forward wavefield1 or backward wavefield2)
  }

  if ( ! ATTENUATION_VISCOACOUSTIC){
    // no attenuation
    if (nb_field == 2){
      // forward wavefields -> FORWARD_OR_ADJOINT == 1
      Kernel_2_acoustic_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                       d_ibool,
                                                                       mp->d_phase_ispec_inner_acoustic,
                                                                       mp->num_phase_ispec_acoustic,
                                                                       d_iphase,
                                                                       mp->d_potential_acoustic,
                                                                       mp->d_potential_dot_dot_acoustic,
                                                                       mp->d_b_potential_acoustic,
                                                                       mp->d_b_potential_dot_dot_acoustic,
                                                                       nb_field,
                                                                       d_xix, d_xiz,
                                                                       d_gammax, d_gammaz,
                                                                       mp->d_hprime_xx,
                                                                       mp->d_hprimewgll_xx,
                                                                       mp->d_wxgll,
                                                                       d_rhostore,
                                                                       mp->pml_boundary_conditions,
                                                                       mp->d_spec_to_pml);
    }else{
      // nb_field == 1
      if (compute_wavefield_1){
        // forward wavefields -> FORWARD_OR_ADJOINT == 1
        Kernel_2_acoustic_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                         d_ibool,
                                                                         mp->d_phase_ispec_inner_acoustic,
                                                                         mp->num_phase_ispec_acoustic,
                                                                         d_iphase,
                                                                         mp->d_potential_acoustic,
                                                                         mp->d_potential_dot_dot_acoustic,
                                                                         mp->d_b_potential_acoustic,
                                                                         mp->d_b_potential_dot_dot_acoustic,
                                                                         nb_field,
                                                                         d_xix, d_xiz,
                                                                         d_gammax, d_gammaz,
                                                                         mp->d_hprime_xx,
                                                                         mp->d_hprimewgll_xx,
                                                                         mp->d_wxgll,
                                                                         d_rhostore,
                                                                         mp->pml_boundary_conditions,
                                                                         mp->d_spec_to_pml);

        if (mp->pml_boundary_conditions){
          Kernel_2_acoustic_PML_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                               d_ibool,
                                                                               mp->d_phase_ispec_inner_acoustic,
                                                                               mp->num_phase_ispec_acoustic,
                                                                               d_iphase,
                                                                               mp->d_potential_acoustic,
                                                                               mp->d_potential_dot_dot_acoustic,
                                                                               d_xix, d_xiz,
                                                                               d_gammax, d_gammaz,
                                                                               mp->d_hprime_xx,
                                                                               mp->d_hprimewgll_xx,
                                                                               mp->d_wxgll,
                                                                               d_rhostore,
                                                                               mp->d_spec_to_pml,
                                                                               mp->ALPHA_MAX_PML,
                                                                               mp->d0_max,
                                                                               mp->abscissa_norm,
                                                                               mp->nspec_pml_x,
                                                                               mp->nspec_pml_z,
                                                                               mp->deltat,
                                                                               mp->PML_dpotentialdxl_old,
                                                                               mp->PML_dpotentialdzl_old,
                                                                               mp->d_potential_old,
                                                                               mp->rmemory_acoustic_dux_dx,
                                                                               mp->rmemory_acoustic_dux_dz,
                                                                               mp->rmemory_acoustic_dux_dx2,
                                                                               mp->rmemory_acoustic_dux_dz2,
                                                                               mp->rmemory_pot_acoustic,
                                                                               mp->rmemory_pot_acoustic2,
                                                                               mp->d_potential_dot_acoustic,
                                                                               mp->d_kappastore,
                                                                               mp->alphax_store,
                                                                               mp->alphaz_store,
                                                                               mp->betax_store,
                                                                               mp->betaz_store);
        } //PML
      } // compute_wavefield1
      if (compute_wavefield_2){
        // this run only happens with UNDO_ATTENUATION_AND_OR_PML on
        // adjoint wavefields -> FORWARD_OR_ADJOINT == 3
        Kernel_2_acoustic_impl<3><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                         d_ibool,
                                                                         mp->d_phase_ispec_inner_acoustic,
                                                                         mp->num_phase_ispec_acoustic,
                                                                         d_iphase,
                                                                         mp->d_b_potential_acoustic,
                                                                         mp->d_b_potential_dot_dot_acoustic,
                                                                         mp->d_b_potential_acoustic,
                                                                         mp->d_b_potential_dot_dot_acoustic,
                                                                         nb_field,
                                                                         d_xix, d_xiz,
                                                                         d_gammax, d_gammaz,
                                                                         mp->d_hprime_xx,
                                                                         mp->d_hprimewgll_xx,
                                                                         mp->d_wxgll,
                                                                         d_rhostore,
                                                                         mp->pml_boundary_conditions,
                                                                         mp->d_spec_to_pml);
      } //compute_wavefield_1
    } //nb_field
  }else{
    // attenuation
    // ATTENUATION_VISCOACOUSTIC == .true. below
    if (compute_wavefield_1) {
      Kernel_2_viscoacoustic_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                            d_ibool,
                                                                            mp->d_phase_ispec_inner_acoustic,
                                                                            mp->num_phase_ispec_acoustic,
                                                                            d_iphase,
                                                                            mp->d_potential_acoustic,
                                                                            mp->d_potential_dot_dot_acoustic,
                                                                            d_xix, d_xiz,
                                                                            d_gammax, d_gammaz,
                                                                            mp->d_hprime_xx,
                                                                            mp->d_hprimewgll_xx,
                                                                            mp->d_wxgll,
                                                                            d_rhostore,
                                                                            mp->d_e1_acous,
                                                                            mp->d_A_newmark_acous,
                                                                            mp->d_B_newmark_acous,
                                                                            mp->d_sum_forces_old);
    }
    if (compute_wavefield_2) {
      Kernel_2_viscoacoustic_impl<3><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                            d_ibool,
                                                                            mp->d_phase_ispec_inner_acoustic,
                                                                            mp->num_phase_ispec_acoustic,
                                                                            d_iphase,
                                                                            mp->d_b_potential_acoustic,
                                                                            mp->d_b_potential_dot_dot_acoustic,
                                                                            d_xix, d_xiz,
                                                                            d_gammax, d_gammaz,
                                                                            mp->d_hprime_xx,
                                                                            mp->d_hprimewgll_xx,
                                                                            mp->d_wxgll,
                                                                            d_rhostore,
                                                                            mp->d_b_e1_acous,
                                                                            mp->d_A_newmark_acous,
                                                                            mp->d_B_newmark_acous,
                                                                            mp->d_b_sum_forces_old);
    }
  } // ATTENUATION_VISCOACOUSTIC

  // Cuda timing
  if (CUDA_TIMING) {
    realw flops,time;
    stop_timing_cuda(&start,&stop,"Kernel_2_acoustic_impl",&time);
    // time in seconds
    time = time / 1000.;
    flops = 15559 * nb_blocks_to_compute;
    printf("  performance: %f GFlop/s\n", flops/time * 1.e-9);
  }

  GPU_ERROR_CHECKING ("Kernel_2_acoustic");
}

/* ----------------------------------------------------------------------------------------------- */

// main compute_forces_acoustic CUDA routine

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_forces_acoustic_cuda,
              COMPUTE_FORCES_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                            int* iphase,
                                            int* nspec_outer_acoustic,
                                            int* nspec_inner_acoustic,
                                            int* ATTENUATION_VISCOACOUSTIC,
                                            int* compute_wavefield_1,
                                            int* compute_wavefield_2) {
  TRACE("compute_forces_acoustic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int num_elements;

  if (*iphase == 1)
    num_elements = *nspec_outer_acoustic;
  else
    num_elements = *nspec_inner_acoustic;

  // checks if anything to do
  if (num_elements == 0) return;

  // no mesh coloring: uses atomic updates
  Kernel_2_acoustic(num_elements, mp, *iphase,
                    mp->d_ibool,
                    mp->d_xix,mp->d_xiz,
                    mp->d_gammax,mp->d_gammaz,
                    mp->d_rhostore,
                    *ATTENUATION_VISCOACOUSTIC,
                    *compute_wavefield_1,
                    *compute_wavefield_2);

}

