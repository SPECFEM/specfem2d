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

/*
//#define CUBLAS_ERROR(s,n)  if (s != CUBLAS_STATUS_SUCCESS) {  \
//fprintf (stderr, "CUBLAS Memory Write Error @ %d\n",n); \
//exit(EXIT_FAILURE); }
*/

/* ----------------------------------------------------------------------------------------------- */

// elastic wavefield

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(update_displacement_cuda,
              UPDATE_DISPLACMENT_CUDA)(long* Mesh_pointer,
                                          realw* deltat_F,
                                          realw* deltatsqover2_F,
                                          realw* deltatover2_F,
                                          realw* b_deltat_F,
                                          realw* b_deltatsqover2_F,
                                          realw* b_deltatover2_F) {

  TRACE("\tupdate_displacement_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  realw deltat = *deltat_F;
  realw deltatsqover2 = *deltatsqover2_F;
  realw deltatover2 = *deltatover2_F;

  int size = NDIM * mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  cudaEvent_t start,stop;
  if (CUDA_TIMING_UPDATE) {
    start_timing_cuda(&start,&stop);
  }

  // debug
  //realw max_d,max_v,max_a;
  //max_d = get_device_array_maximum_value(mp->d_displ, size);
  //max_v = get_device_array_maximum_value(mp->d_veloc, size);
  //max_a = get_device_array_maximum_value(mp->d_accel, size);
  //printf("rank %d - max displ: %f veloc: %f accel: %f\n",mp->myrank,max_d,max_v,max_a);

  //launch kernel
  UpdateDispVeloc_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_displ,mp->d_veloc,mp->d_accel,
                                                                size,deltat,deltatsqover2,deltatover2);

  // kernel for backward fields
  if (mp->simulation_type == 3) {
    realw b_deltat = *b_deltat_F;
    realw b_deltatsqover2 = *b_deltatsqover2_F;
    realw b_deltatover2 = *b_deltatover2_F;

    UpdateDispVeloc_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_displ,mp->d_b_veloc,mp->d_b_accel,
                                                                  size,b_deltat,b_deltatsqover2,b_deltatover2);
  }

  // Cuda timing
  if (CUDA_TIMING_UPDATE) {
    realw flops,time;
    stop_timing_cuda(&start,&stop,"UpdateDispVeloc_kernel",&time);
    // time in seconds
    time = time / 1000.;
    // performance: 6 FLOPS per thread
    flops = 6.0 * size;
    //printf("  performance: %f GFlop/s num_blocks x/y: %d %d threads: %d\n", flops/time * 1.e-9,num_blocks_x,num_blocks_y,size);
    printf("  performance: %f GFlop/s\n", flops/time * 1.e-9);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("update_displacement_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// acoustic wavefield

// KERNEL 1
/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(update_displacement_ac_cuda,
              UPDATE_DISPLACEMENT_AC_CUDA)(long* Mesh_pointer,
                                           realw* deltat_F,
                                           realw* deltatsqover2_F,
                                           realw* deltatover2_F,
                                           realw* b_deltat_F,
                                           realw* b_deltatsqover2_F,
                                           realw* b_deltatover2_F,
                                           int* compute_b_wavefield,
                                           int* UNDO_ATTENUATION_AND_OR_PML) {
  TRACE("\tupdate_displacement_ac_cuda");
  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  //launch kernel
  // forward wavefields
  realw deltat = *deltat_F;
  realw deltatsqover2 = *deltatsqover2_F;
  realw deltatover2 = *deltatover2_F;

  //cudaEventRecord(mp->end_of_iteration,mp->compute_stream);
  //cudaEventSynchronize(mp->end_of_iteration);

  // Cuda timing
  cudaEvent_t start,stop;
  if (CUDA_TIMING_UPDATE) {
    start_timing_cuda(&start,&stop);
  }

  if (!(*UNDO_ATTENUATION_AND_OR_PML && *compute_b_wavefield)) {
    UpdatePotential_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_acoustic,
                                                                  mp->d_potential_dot_acoustic,
                                                                  mp->d_potential_dot_dot_acoustic,
                                                                  size,deltat,deltatsqover2,deltatover2);
  }

  // backward/reconstructed wavefields
  if (mp->simulation_type == 3 && *compute_b_wavefield) {
    realw b_deltat = *b_deltat_F;
    realw b_deltatsqover2 = *b_deltatsqover2_F;
    realw b_deltatover2 = *b_deltatover2_F;
    UpdatePotential_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_acoustic,
                                                                  mp->d_b_potential_dot_acoustic,
                                                                  mp->d_b_potential_dot_dot_acoustic,
                                                                  size,b_deltat,b_deltatsqover2,b_deltatover2);
  }

  // Cuda timing
  if (CUDA_TIMING_UPDATE) {
    realw flops,time;
    stop_timing_cuda(&start,&stop,"UpdatePotential_kernel",&time);
    // time in seconds
    time = time / 1000.;
    // performance
    // see with: nvprof --metrics flops_sp ./xspecfem3D
    //           -> using 8199750 FLOPS (Single) floating-point operations for 1366625 threads
    //              = 6 FLOPS per thread
    flops = 6.0 * size;
    //printf("  performance: %f GFlop/s num_blocks x/y: %d %d threads: %d\n", flops/time * 1.e-9,num_blocks_x,num_blocks_y,size);
    printf("  performance: %f GFlop/s\n", flops/time * 1.e-9);
  }

  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);

  GPU_ERROR_CHECKING ("update_displacement_ac_cuda");
}


/* ----------------------------------------------------------------------------------------------- */

// elastic domains

// KERNEL 3

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(kernel_3_a_cuda,
              KERNEL_3_A_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               realw* b_deltatover2_F,
                               int* compute_wavefield_1,
                               int* compute_wavefield_2) {

  TRACE("\tkernel_3_a_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  realw deltatover2 = *deltatover2_F;

  // updates both, accel and veloc
  if (*compute_wavefield_1){
    kernel_3_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_veloc,
                                                                  mp->d_accel,
                                                                  size, deltatover2,
                                                                  mp->d_rmassx,mp->d_rmassz);
  }
  // kernel for backward fields
  if (*compute_wavefield_2) {
    realw b_deltatover2 = *b_deltatover2_F;
    kernel_3_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_b_veloc,
                                                                  mp->d_b_accel,
                                                                  size, b_deltatover2,
                                                                  mp->d_rmassx,mp->d_rmassz);
  }

  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);

  GPU_ERROR_CHECKING ("kernel_3_a_cuda");
}

/* ----------------------------------------------------------------------------------------------- */

/* unused so far...

// updates velocity only

extern "C"
void FC_FUNC_(kernel_3_b_cuda,
              KERNEL_3_B_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               realw* b_deltatover2_F) {
  TRACE("\tkernel_3_b_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL3;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  realw deltatover2 = *deltatover2_F;

  // updates only veloc at this point
  kernel_3_veloc_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_veloc,
                                                                      mp->d_accel,
                                                                      size,deltatover2);
  // kernel for backward fields
  if (mp->simulation_type == 3) {
    realw b_deltatover2 = *b_deltatover2_F;
    kernel_3_veloc_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_b_veloc,
                                                                        mp->d_b_accel,
                                                                        size,b_deltatover2);
  }

  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);

  GPU_ERROR_CHECKING ("kernel_3_b_cuda");
}
*/

/* ----------------------------------------------------------------------------------------------- */

// acoustic domains

// KERNEL 3

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(kernel_3_acoustic_cuda,
              KERNEL_3_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                      realw* deltatover2,
                                      realw* b_deltatover2,
                                      int* compute_wavefield_1,
                                      int* compute_wavefield_2) {

  TRACE("\tkernel_3_acoustic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_AB;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  kernel_3_acoustic_cuda_device<<< grid, threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                         mp->d_b_potential_dot_dot_acoustic,
                                                                         mp->d_potential_dot_acoustic,
                                                                         mp->d_b_potential_dot_acoustic,
                                                                         size,
                                                                         *compute_wavefield_1,
                                                                         *compute_wavefield_2,
                                                                         *deltatover2,
                                                                         *b_deltatover2,
                                                                          mp->d_rmass_acoustic);

  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);

 GPU_ERROR_CHECKING ("kernel_3_acoustic_cuda");
}


