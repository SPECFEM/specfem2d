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
// (elements on boundary)

extern "C"
void FC_FUNC_(transfer_boun_accel_from_device,
              TRANSFER_BOUN_ACCEL_FROM_DEVICE)(long* Mesh_pointer,
                                               realw* h_send_accel_buffer,
                                               const int* FORWARD_OR_ADJOINT){
TRACE("\ttransfer_boun_accel_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if (mp->size_mpi_buffer == 0) return;
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) return;

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // sets gpu arrays
  realw* accel;
  realw* send_buffer;
  if (*FORWARD_OR_ADJOINT == 1) {
    // forward
    accel = mp->d_accel;
    send_buffer = mp->d_send_accel_buffer;
  }else if (*FORWARD_OR_ADJOINT == 3) {
    // backward/reconstructed
    accel = mp->d_b_accel;
    send_buffer = mp->d_b_send_accel_buffer;
  }

  // Cuda timing
  //cudaEvent_t start, stop;
  //start_timing_cuda(&start,&stop);

  prepare_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(accel,
                                                                          send_buffer,
                                                                          mp->ninterface_elastic,
                                                                          mp->max_nibool_interfaces_ext_mesh,
                                                                          mp->d_nibool_interfaces_ext_mesh,
                                                                          mp->d_ibool_interfaces_ext_mesh,
                                                                          mp->d_inum_interfaces_elastic);


  // synchronizes
  // explicitly waits until previous compute stream finishes
  // (cudaMemcpy implicitly synchronizes all other cuda operations)
  cudaStreamSynchronize(mp->compute_stream);

  // copies buffer from GPU to CPU host
  print_CUDA_error_if_any(cudaMemcpy(h_send_accel_buffer,send_buffer,
                          mp->size_mpi_buffer*sizeof(realw),cudaMemcpyDeviceToHost),97001);

  // Cuda timing
  // finish timing of kernel+memcpy
  //stop_timing_cuda(&start,&stop,"prepare_boundary_accel_on_device");

  GPU_ERROR_CHECKING ("transfer_boun_accel_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_boundary_from_device_a,
              TRANSFER_BOUNDARY_FROM_DEVICE_A)(long* Mesh_pointer) {

// asynchronous transfer from device to host

  TRACE("\ttransfer_boundary_from_device_a");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  // checks if anything to do
  if (mp->size_mpi_buffer == 0) return;

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  prepare_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,mp->d_send_accel_buffer,
                                                                          mp->ninterface_elastic,
                                                                          mp->max_nibool_interfaces_ext_mesh,
                                                                          mp->d_nibool_interfaces_ext_mesh,
                                                                          mp->d_ibool_interfaces_ext_mesh,
                                                                          mp->d_inum_interfaces_elastic);
  // waits until kernel is finished before starting async memcpy
  //synchronize_cuda();
  // waits until previous compute stream finishes
  cudaStreamSynchronize(mp->compute_stream);

  cudaMemcpyAsync(mp->h_send_accel_buffer,mp->d_send_accel_buffer,
                  mp->size_mpi_buffer*sizeof(realw),cudaMemcpyDeviceToHost,mp->copy_stream);

  GPU_ERROR_CHECKING ("transfer_boundary_from_device_a");
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(prepare_boundary_on_device,
              PREPARE_BOUNDARY_ON_DEVICE)(long* Mesh_pointer) {

// asynchronous transfer from device to host

  TRACE("\tprepare_boundary_on_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  // checks if anything to do
  if (mp->size_mpi_buffer == 0) return;

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  prepare_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,mp->d_send_accel_buffer,
                                                                          mp->ninterface_elastic,
                                                                          mp->max_nibool_interfaces_ext_mesh,
                                                                          mp->d_nibool_interfaces_ext_mesh,
                                                                          mp->d_ibool_interfaces_ext_mesh,
                                                                          mp->d_inum_interfaces_elastic);
  // waits until kernel is finished before starting async memcpy
  //synchronize_cuda();
  // waits until previous compute stream finishes
  cudaStreamSynchronize(mp->compute_stream);

  GPU_ERROR_CHECKING ("prepare_boundary_on_device");
}



/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_boundary_to_device_a,
              TRANSFER_BOUNDARY_TO_DEVICE_A)(long* Mesh_pointer,
                                             realw* buffer_recv_vector_gpu,
                                             const int* max_nibool_interfaces_ext_mesh) {

// asynchronous transfer from host to device

  TRACE("transfer_boundary_to_device_a");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  if (mp->size_mpi_buffer > 0) {
    // copy on host memory
    memcpy(mp->h_recv_accel_buffer,buffer_recv_vector_gpu,mp->size_mpi_buffer*sizeof(realw));

    // asynchronous copy to GPU using copy_stream
    cudaMemcpyAsync(mp->d_send_accel_buffer,mp->h_recv_accel_buffer,
                    mp->size_mpi_buffer*sizeof(realw),cudaMemcpyHostToDevice,mp->copy_stream);
  }

  GPU_ERROR_CHECKING ("transfer_boundary_to_device_a");
}


/* ----------------------------------------------------------------------------------------------- */

// Assembly

/* ----------------------------------------------------------------------------------------------- */


// FORWARD_OR_ADJOINT == 1 for accel, and == 3 for b_accel
extern "C"
void FC_FUNC_(transfer_asmbl_accel_to_device,
              TRANSFER_ASMBL_ACCEL_TO_DEVICE)(long* Mesh_pointer,
                                              realw* buffer_recv_vector_gpu,
                                              const int* max_nibool_interfaces_ext_mesh,
                                              const int* nibool_interfaces_ext_mesh,
                                              const int* ibool_interfaces_ext_mesh,
                                              const int* FORWARD_OR_ADJOINT) {
  TRACE("\ttransfer_asmbl_accel_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if (mp->size_mpi_buffer == 0) return;
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) return;

  //daniel: todo - check if this copy is only needed for adjoint simulation, otherwise it is called asynchronously?
  if (*FORWARD_OR_ADJOINT == 1) {
    // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
    cudaStreamSynchronize(mp->copy_stream);
  }
  else if (*FORWARD_OR_ADJOINT == 3) {
    // explicitly synchronizes
    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    synchronize_cuda();

    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_send_accel_buffer, buffer_recv_vector_gpu,
                            mp->size_mpi_buffer*sizeof(realw),cudaMemcpyHostToDevice),97001);
  }

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // sets gpu arrays
  realw* accel;
  realw* send_buffer;
  if (*FORWARD_OR_ADJOINT == 1) {
    // forward
    accel = mp->d_accel;
    send_buffer = mp->d_send_accel_buffer;
  }else if (*FORWARD_OR_ADJOINT == 3) {
    // backward/reconstructed
    accel = mp->d_b_accel;
    send_buffer = mp->d_b_send_accel_buffer;
  }

  //double start_time = get_time();
  // cudaEvent_t start, stop;
  // realw time;
  // cudaEventCreate(&start);
  // cudaEventCreate(&stop);
  // cudaEventRecord( start, 0 );

  //assembles
  assemble_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(accel,
                                                                           send_buffer,
                                                                           mp->ninterface_elastic,
                                                                           mp->max_nibool_interfaces_ext_mesh,
                                                                           mp->d_nibool_interfaces_ext_mesh,
                                                                           mp->d_ibool_interfaces_ext_mesh,
                                                                           mp->d_inum_interfaces_elastic);

  // cudaEventRecord( stop, 0 );
  // cudaEventSynchronize( stop );
  // cudaEventElapsedTime( &time, start, stop );
  // cudaEventDestroy( start );
  // cudaEventDestroy( stop );
  // printf("Boundary Assemble Kernel Execution Time: %f ms\n",time);

  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);

  GPU_ERROR_CHECKING ("transfer_asmbl_accel_to_device");
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(sync_copy_from_device,
              SYNC_copy_FROM_DEVICE)(long* Mesh_pointer,
                                     int* iphase,
                                     realw* send_buffer) {

  TRACE("sync_copy_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  // Wait until async-memcpy of outer elements is finished and start MPI.
  if (*iphase != 2) { exit_on_cuda_error("sync_copy_from_device must be called for iphase == 2"); }

  if (mp->size_mpi_buffer > 0) {
    // waits for asynchronous copy to finish
    cudaStreamSynchronize(mp->copy_stream);

    // There have been problems using the pinned-memory with MPI, so
    // we copy the buffer into a non-pinned region.
    memcpy(send_buffer,mp->h_send_accel_buffer,mp->size_mpi_buffer*sizeof(float));
  }
  // memory copy is now finished, so non-blocking MPI send can proceed

  GPU_ERROR_CHECKING ("sync_copy_from_device");
}

