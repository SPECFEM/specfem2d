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


extern "C"
void FC_FUNC_(compute_stacey_viscoelastic_cuda,
              COMPUTE_STACEY_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                                int* iphasef,
                                                realw* h_b_absorb_elastic_left,
                                                realw* h_b_absorb_elastic_right,
                                                realw* h_b_absorb_elastic_top,
                                                realw* h_b_absorb_elastic_bottom,
                                                int* compute_wavefield_1,
                                                int* compute_wavefield_2,
                                                int *UNDO_ATTENUATION_AND_OR_PML) {

  TRACE("compute_stacey_viscoelastic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if (mp->d_num_abs_boundary_faces == 0) return;

  int iphase = *iphasef;

  // only add this contribution for first pass
  if (iphase != 1) return;

  int blocksize = NGLLX;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->d_num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // We have to distinguish between a UNDO_ATTENUATION_AND_OR_PML run or not to know if read/write operations are necessary
  int read_abs = (mp->simulation_type == 3 && (! *UNDO_ATTENUATION_AND_OR_PML)) ? 1 : 0;
  int write_abs = (mp->simulation_type == 1 && mp->save_forward && (! *UNDO_ATTENUATION_AND_OR_PML)) ? 1 : 0;

  if (read_abs) {
    // reading is done in fortran routine
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_elastic_left,h_b_absorb_elastic_left,
                                       2*mp->d_nspec_left*sizeof(realw)*NGLLX,cudaMemcpyHostToDevice),7700);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_elastic_right,h_b_absorb_elastic_right,
                                       2*mp->d_nspec_right*sizeof(realw)*NGLLX,cudaMemcpyHostToDevice),7700);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_elastic_top,h_b_absorb_elastic_top,
                                       2*mp->d_nspec_top*sizeof(realw)*NGLLX,cudaMemcpyHostToDevice),7700);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_elastic_bottom,h_b_absorb_elastic_bottom,
                                       2*mp->d_nspec_bottom*sizeof(realw)*NGLLX,cudaMemcpyHostToDevice),7700);
  }

  GPU_ERROR_CHECKING ("between cudamemcpy and compute_stacey_elastic_kernel");

  if (*compute_wavefield_1){
    // forward/adjoint wavefield
    compute_stacey_elastic_kernel<<<grid,threads,0,mp->compute_stream>>>( mp->d_veloc,
                                                                          mp->d_accel,
                                                                          mp->d_abs_boundary_ispec,
                                                                          mp->d_abs_boundary_ijk,
                                                                          mp->d_abs_boundary_normal,
                                                                          mp->d_abs_boundary_jacobian2Dw,
                                                                          mp->d_ibool,
                                                                          mp->d_rho_vp,
                                                                          mp->d_rho_vs,
                                                                          mp->d_ispec_is_elastic,
                                                                          mp->simulation_type,
                                                                          mp->p_sv,
                                                                          mp->save_forward,
                                                                          mp->d_num_abs_boundary_faces,
                                                                          mp->d_b_absorb_elastic_left,
                                                                          mp->d_b_absorb_elastic_right,
                                                                          mp->d_b_absorb_elastic_top,
                                                                          mp->d_b_absorb_elastic_bottom,
                                                                          mp->d_ib_left,
                                                                          mp->d_ib_right,
                                                                          mp->d_ib_top,
                                                                          mp->d_ib_bottom,
                                                                          mp->d_edge_abs);
  }

  // kernel simulations
  if (*compute_wavefield_2){
    // backward/reconstructed wavefield
    compute_stacey_elastic_sim3_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_abs_boundary_ispec,
                                                                              mp->d_abs_boundary_ijk,
                                                                              mp->d_ibool,
                                                                              mp->d_ispec_is_elastic,
                                                                              mp->d_num_abs_boundary_faces,
                                                                              mp->d_b_accel,
                                                                              mp->d_b_absorb_elastic_left,
                                                                              mp->d_b_absorb_elastic_right,
                                                                              mp->d_b_absorb_elastic_top,
                                                                              mp->d_b_absorb_elastic_bottom,
                                                                              mp->d_ib_left,
                                                                              mp->d_ib_right,
                                                                              mp->d_ib_top,
                                                                              mp->d_ib_bottom,
                                                                              mp->d_edge_abs);
  }

  GPU_ERROR_CHECKING ("compute_stacey_elastic_kernel");

  if (write_abs) {
    // explicitly wait until compute stream is done
    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    cudaStreamSynchronize(mp->compute_stream);

    // writing is done in fortran routine
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_elastic_left,mp->d_b_absorb_elastic_left,
                                       2*mp->d_nspec_left*sizeof(realw)*NGLLX,cudaMemcpyDeviceToHost),7801);
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_elastic_right,mp->d_b_absorb_elastic_right,
                                       2*mp->d_nspec_right*sizeof(realw)*NGLLX,cudaMemcpyDeviceToHost),7802);
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_elastic_top,mp->d_b_absorb_elastic_top,
                                       2*mp->d_nspec_top*sizeof(realw)*NGLLX,cudaMemcpyDeviceToHost),7803);
    print_CUDA_error_if_any(cudaMemcpy(h_b_absorb_elastic_bottom,mp->d_b_absorb_elastic_bottom,
                                       2*mp->d_nspec_bottom*sizeof(realw)*NGLLX,cudaMemcpyDeviceToHost),7804);
  }

  GPU_ERROR_CHECKING ("after compute_stacey_elastic after cudamemcpy");
}

