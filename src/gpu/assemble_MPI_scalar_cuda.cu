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

// ASSEMBLY - mpi data transfer between CPU-GPU

/* ----------------------------------------------------------------------------------------------- */

// prepares and transfers the inter-element edge-nodes to the host to be MPI'd

extern "C"
void FC_FUNC_(transfer_boun_pot_from_device,
              TRANSFER_BOUN_POT_FROM_DEVICE)(long* Mesh_pointer,
                                             realw* h_send_potential_dot_dot_buffer,
                                             const int* FORWARD_OR_ADJOINT){

  TRACE("transfer_boun_pot_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if (mp->size_mpi_buffer_potential == 0) return;
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) return;

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)(mp->max_nibool_interfaces_ext_mesh))/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // sets gpu arrays
  realw* potential;
  realw* send_buffer;
  if (*FORWARD_OR_ADJOINT == 1) {
    // forward
    potential = mp->d_potential_dot_dot_acoustic;
    send_buffer = mp->d_send_potential_dot_dot_buffer;
  }else if (*FORWARD_OR_ADJOINT == 3) {
    // backward/reconstructed
    potential = mp->d_b_potential_dot_dot_acoustic;
    send_buffer = mp->d_b_send_potential_dot_dot_buffer;
  }

  prepare_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(potential,
                                                                              send_buffer,
                                                                              mp->ninterface_acoustic,
                                                                              mp->max_nibool_interfaces_ext_mesh,
                                                                              mp->d_nibool_interfaces_ext_mesh,
                                                                              mp->d_ibool_interfaces_ext_mesh,
                                                                              mp->d_inum_interfaces_acoustic);

  // synchronizes
  // explicitly waits until previous compute stream finishes
  // (cudaMemcpy implicitly synchronizes all other cuda operations)
  cudaStreamSynchronize(mp->compute_stream);

  print_CUDA_error_if_any(cudaMemcpy(h_send_potential_dot_dot_buffer,send_buffer,
                                     mp->size_mpi_buffer_potential*sizeof(realw),cudaMemcpyDeviceToHost),98000);

  // finish timing of kernel+memcpy
  // cudaEventRecord( stop, 0 );
  // cudaEventSynchronize( stop );
  // cudaEventElapsedTime( &time, start, stop );
  // cudaEventDestroy( start );
  // cudaEventDestroy( stop );
  // printf("boundary xfer d->h Time: %f ms\n",time);

  GPU_ERROR_CHECKING ("transfer_boun_pot_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

// Assembly

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(transfer_asmbl_pot_to_device,
              TRANSFER_ASMBL_POT_TO_DEVICE)(long* Mesh_pointer,
                                            realw* h_buffer_recv_scalar_gpu,
                                            const int* FORWARD_OR_ADJOINT) {

  TRACE("transfer_asmbl_pot_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // Cuda timing
  //cudaEvent_t start, stop;
  //start_timing_cuda(&start,&stop);

  // checks if anything to do
  if (mp->size_mpi_buffer_potential == 0) return;
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) return;

  // assembles on GPU
  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // sets gpu arrays
  realw* potential;
  realw* send_buffer;
  if (*FORWARD_OR_ADJOINT == 1) {
    //assemble forward field
    potential = mp->d_potential_dot_dot_acoustic;
    send_buffer = mp->d_send_potential_dot_dot_buffer;
  }else if (*FORWARD_OR_ADJOINT == 3) {
    //assemble reconstructed/backward field
    potential = mp->d_b_potential_dot_dot_acoustic;
    send_buffer = mp->d_b_send_potential_dot_dot_buffer;
  }

  // copies buffer onto GPU
  print_CUDA_error_if_any(cudaMemcpy(send_buffer, h_buffer_recv_scalar_gpu,
                                     mp->size_mpi_buffer_potential*sizeof(realw), cudaMemcpyHostToDevice),98010);

  //assembles field
  assemble_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(potential,
                                                                               send_buffer,
                                                                               mp->ninterface_acoustic,
                                                                               mp->max_nibool_interfaces_ext_mesh,
                                                                               mp->d_nibool_interfaces_ext_mesh,
                                                                               mp->d_ibool_interfaces_ext_mesh,
                                                                               mp->d_inum_interfaces_acoustic);

  // Cuda timing
  //stop_timing_cuda(&start,&stop,"assemble_boundary_potential_on_device");

  GPU_ERROR_CHECKING ("transfer_asmbl_pot_to_device");
}

