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

// KERNEL 2 - visco-elastic forces main kernel

/* ----------------------------------------------------------------------------------------------- */

void Kernel_2(int nb_blocks_to_compute,Mesh* mp,int d_iphase,
              realw d_deltat,
              int ANISOTROPY,
              int ATTENUATION_VISCOELASTIC,
              int compute_wavefield_1,
              int compute_wavefield_2) {

TRACE("Kernel_2");

  // if the grid can handle the number of blocks, we let it be 1D

  int blocksize = NGLL2_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(nb_blocks_to_compute,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  cudaEvent_t start,stop;
  if (CUDA_TIMING) {
    start_timing_cuda(&start,&stop);
  }

  // safety check
  if (mp->simulation_type == 3) {
    if (ATTENUATION_VISCOELASTIC) exit_on_error("GPU_MODE not supported yet for adjoint simulations with attenuation viscoelastic");
  }

  // compute kernels without attenuation
  if (ANISOTROPY) {
    // full anisotropy
    if (ATTENUATION_VISCOELASTIC){
      // anisotropy, attenuation
      if (compute_wavefield_1){
        // forward wavefields -> FORWARD_OR_ADJOINT == 1
        TRACE("\tKernel_2_att_ani_impl 1");
        Kernel_2_att_ani_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                        mp->d_ibool,
                                                                        mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                        d_iphase,
                                                                        mp->d_displ,
                                                                        mp->d_accel,
                                                                        mp->d_xix, mp->d_xiz,
                                                                        mp->d_gammax, mp->d_gammaz,
                                                                        mp->d_hprime_xx,
                                                                        mp->d_hprimewgll_xx,
                                                                        mp->d_wxgll,
                                                                        mp->d_kappav,
                                                                        mp->d_muv,
                                                                        mp->simulation_type,
                                                                        mp->p_sv,
                                                                        mp->d_ispec_is_anisotropic,
                                                                        mp->d_c11store,mp->d_c12store,mp->d_c13store,
                                                                        mp->d_c15store,mp->d_c23store,mp->d_c25store,
                                                                        mp->d_c33store,mp->d_c35store,mp->d_c55store,
                                                                        mp->d_A_newmark_mu,
                                                                        mp->d_B_newmark_mu,
                                                                        mp->d_A_newmark_kappa,
                                                                        mp->d_B_newmark_kappa,
                                                                        mp->d_e1,
                                                                        mp->d_e11,
                                                                        mp->d_e13,
                                                                        mp->d_dux_dxl_old,
                                                                        mp->d_duz_dzl_old,
                                                                        mp->d_dux_dzl_plus_duz_dxl_old);
      }
      if (compute_wavefield_2){
        // this run only happens with UNDO_ATTENUATION_AND_OR_PML on
        // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
        TRACE("\tKernel_2_att_ani_impl 3");
        Kernel_2_att_ani_impl<3><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                        mp->d_ibool,
                                                                        mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                        d_iphase,
                                                                        mp->d_b_displ,
                                                                        mp->d_b_accel,
                                                                        mp->d_xix, mp->d_xiz,
                                                                        mp->d_gammax, mp->d_gammaz,
                                                                        mp->d_hprime_xx,
                                                                        mp->d_hprimewgll_xx,
                                                                        mp->d_wxgll,
                                                                        mp->d_kappav,
                                                                        mp->d_muv,
                                                                        mp->simulation_type,
                                                                        mp->p_sv,
                                                                        mp->d_ispec_is_anisotropic,
                                                                        mp->d_c11store,mp->d_c12store,mp->d_c13store,
                                                                        mp->d_c15store,mp->d_c23store,mp->d_c25store,
                                                                        mp->d_c33store,mp->d_c35store,mp->d_c55store,
                                                                        mp->d_A_newmark_mu,
                                                                        mp->d_B_newmark_mu,
                                                                        mp->d_A_newmark_kappa,
                                                                        mp->d_B_newmark_kappa,
                                                                        mp->d_b_e1,
                                                                        mp->d_b_e11,
                                                                        mp->d_b_e13,
                                                                        mp->d_b_dux_dxl_old,
                                                                        mp->d_b_duz_dzl_old,
                                                                        mp->d_b_dux_dzl_plus_duz_dxl_old);
      }
    }else{
      // anisotropy, no attenuation
      if (compute_wavefield_1){
        // forward wavefields -> FORWARD_OR_ADJOINT == 1
        TRACE("\tKernel_2_noatt_ani_impl 1");
        Kernel_2_noatt_ani_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                          mp->d_ibool,
                                                                          mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                          d_iphase,
                                                                          mp->d_displ,
                                                                          mp->d_accel,
                                                                          mp->d_xix, mp->d_xiz,
                                                                          mp->d_gammax, mp->d_gammaz,
                                                                          mp->d_hprime_xx,
                                                                          mp->d_hprimewgll_xx,
                                                                          mp->d_wxgll,
                                                                          mp->d_kappav,
                                                                          mp->d_muv,
                                                                          mp->simulation_type,
                                                                          mp->p_sv,
                                                                          mp->d_ispec_is_anisotropic,
                                                                          mp->d_c11store,mp->d_c12store,mp->d_c13store,
                                                                          mp->d_c15store,mp->d_c23store,mp->d_c25store,
                                                                          mp->d_c33store,mp->d_c35store,mp->d_c55store);
      }
      if (compute_wavefield_2){
        // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
        TRACE("\tKernel_2_noatt_ani_impl 3");
        Kernel_2_noatt_ani_impl<3><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                          mp->d_ibool,
                                                                          mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                          d_iphase,
                                                                          mp->d_b_displ,
                                                                          mp->d_b_accel,
                                                                          mp->d_xix, mp->d_xiz,
                                                                          mp->d_gammax, mp->d_gammaz,
                                                                          mp->d_hprime_xx,
                                                                          mp->d_hprimewgll_xx,
                                                                          mp->d_wxgll,
                                                                          mp->d_kappav,
                                                                          mp->d_muv,
                                                                          mp->simulation_type,
                                                                          mp->p_sv,
                                                                          mp->d_ispec_is_anisotropic,
                                                                          mp->d_c11store,mp->d_c12store,mp->d_c13store,
                                                                          mp->d_c15store,mp->d_c23store,mp->d_c25store,
                                                                          mp->d_c33store,mp->d_c35store,mp->d_c55store);
      }
    }
  }else{
    // isotropic
    if (ATTENUATION_VISCOELASTIC){
      // isotropic, attenuation
      if (compute_wavefield_1){
        // forward wavefields -> FORWARD_OR_ADJOINT == 1
        TRACE("\tKernel_2_att_iso_impl 1");
        Kernel_2_att_iso_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                        mp->d_ibool,
                                                                        mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                        d_iphase,
                                                                        mp->d_displ,
                                                                        mp->d_accel,
                                                                        mp->d_xix, mp->d_xiz,
                                                                        mp->d_gammax, mp->d_gammaz,
                                                                        mp->d_hprime_xx,
                                                                        mp->d_hprimewgll_xx,
                                                                        mp->d_wxgll,
                                                                        mp->d_kappav,
                                                                        mp->d_muv,
                                                                        mp->simulation_type,
                                                                        mp->p_sv,
                                                                        mp->d_A_newmark_mu,
                                                                        mp->d_B_newmark_mu,
                                                                        mp->d_A_newmark_kappa,
                                                                        mp->d_B_newmark_kappa,
                                                                        mp->d_e1,
                                                                        mp->d_e11,
                                                                        mp->d_e13,
                                                                        mp->d_dux_dxl_old,
                                                                        mp->d_duz_dzl_old,
                                                                        mp->d_dux_dzl_plus_duz_dxl_old);
      }
      if (compute_wavefield_2){
        // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
        TRACE("\tKernel_2_att_iso_impl 3");
        Kernel_2_att_iso_impl<3><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                        mp->d_ibool,
                                                                        mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                        d_iphase,
                                                                        mp->d_b_displ,
                                                                        mp->d_b_accel,
                                                                        mp->d_xix, mp->d_xiz,
                                                                        mp->d_gammax, mp->d_gammaz,
                                                                        mp->d_hprime_xx,
                                                                        mp->d_hprimewgll_xx,
                                                                        mp->d_wxgll,
                                                                        mp->d_kappav,
                                                                        mp->d_muv,
                                                                        mp->simulation_type,
                                                                        mp->p_sv,
                                                                        mp->d_A_newmark_mu,
                                                                        mp->d_B_newmark_mu,
                                                                        mp->d_A_newmark_kappa,
                                                                        mp->d_B_newmark_kappa,
                                                                        mp->d_b_e1,
                                                                        mp->d_b_e11,
                                                                        mp->d_b_e13,
                                                                        mp->d_b_dux_dxl_old,
                                                                        mp->d_b_duz_dzl_old,
                                                                        mp->d_b_dux_dzl_plus_duz_dxl_old);
      }
    } else {
      // isotropic, no attenuation
      // without storing strains
      if (compute_wavefield_1){
        // forward wavefields -> FORWARD_OR_ADJOINT == 1
        TRACE("\tKernel_2_noatt_iso_impl 1");
        Kernel_2_noatt_iso_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                          mp->d_ibool,
                                                                          mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                          d_iphase,
                                                                          mp->d_displ,
                                                                          mp->d_accel,
                                                                          mp->d_xix, mp->d_xiz,
                                                                          mp->d_gammax, mp->d_gammaz,
                                                                          mp->d_hprime_xx,
                                                                          mp->d_hprimewgll_xx,
                                                                          mp->d_wxgll,
                                                                          mp->d_kappav,
                                                                          mp->d_muv,
                                                                          mp->simulation_type,
                                                                          mp->p_sv);
      }
      // backward/reconstructed wavefield
      if (compute_wavefield_2){
        // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
        TRACE("\tKernel_2_noatt_iso_impl 3");
        Kernel_2_noatt_iso_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                           mp->d_ibool,
                                                                           mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                           d_iphase,
                                                                           mp->d_b_displ,
                                                                           mp->d_b_accel,
                                                                           mp->d_xix, mp->d_xiz,
                                                                           mp->d_gammax,mp->d_gammaz,
                                                                           mp->d_hprime_xx,
                                                                           mp->d_hprimewgll_xx,
                                                                           mp->d_wxgll,
                                                                           mp->d_kappav,
                                                                           mp->d_muv,
                                                                           mp->simulation_type,
                                                                           mp->p_sv);
      }
    }
  } // ANISOTROPY

  // Cuda timing
  if (CUDA_TIMING) {
    if (ANISOTROPY) {
      stop_timing_cuda(&start,&stop,"Kernel_2_noatt_ani_impl");
    }else{
      realw time;
      stop_timing_cuda(&start,&stop,"Kernel_2_noatt_iso_impl",&time);
      // time in seconds
      time = time / 1000.;
      // performance
      // see with: nvprof --metrics flops_sp ./xspecfem3D -> using 883146240 FLOPS (Single) floating-point operations
      // hand-counts: 89344 * number-of-blocks
      realw flops = 89344 * nb_blocks_to_compute;
      printf("  performance: %f GFlops/s\n", flops/time *(1./1000./1000./1000.));
    }
  }

  GPU_ERROR_CHECKING ("Kernel_2_impl");
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_forces_viscoelastic_cuda,
              COMPUTE_FORCES_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                                int* iphase,
                                                realw* deltat,
                                                int* nspec_outer_elastic,
                                                int* nspec_inner_elastic,
                                                int* ANISOTROPY,
                                                int* ATTENUATION_VISCOELASTIC,
                                                int* compute_wavefield_1,
                                                int* compute_wavefield_2) {

  TRACE("compute_forces_viscoelastic_cuda");
  // EPIK_TRACER("compute_forces_viscoelastic_cuda");

  //printf("Running compute_forces\n");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int num_elements;

  if (*iphase == 1){
    num_elements = *nspec_outer_elastic;
  }else{
    num_elements = *nspec_inner_elastic;
  }

  // checks if anything to do
  if (num_elements == 0) return;

  // no mesh coloring: uses atomic updates
  Kernel_2(num_elements,mp,*iphase,
           *deltat,
           *ANISOTROPY,
           *ATTENUATION_VISCOELASTIC,
           *compute_wavefield_1,
           *compute_wavefield_2);

}
