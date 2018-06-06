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

/* trivia

- for most working arrays we use now "realw" instead of "float" type declarations to make it easier to switch
  between a real or double precision simulation
  (matching CUSTOM_REAL == 4 or 8 in fortran routines).

- instead of boolean "logical" declared in fortran routines, in C (or Cuda-C) we have to use "int" variables.
  ifort / gfortran caveat:
    to check whether it is true or false, do not check for == 1 to test for true values since ifort just uses
    non-zero values for true (e.g. can be -1 for true). however, false will be always == 0.
  thus, rather use: if (var) {...}  for testing if true instead of if (var == 1){...} (alternative: one could use if (var != 0) {...}

*/

#ifndef MESH_CONSTANTS_CUDA_H
#define MESH_CONSTANTS_CUDA_H

#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <cuda.h>
//#include <cuda_runtime.h>
//#include <cublas.h>

#include "config.h"


/* ----------------------------------------------------------------------------------------------- */

// for debugging and benchmarking

/* ----------------------------------------------------------------------------------------------- */

#define DEBUG 0
#if DEBUG == 1
#define TRACE(x) printf("%s\n",x);
#else
#define TRACE(x) // printf("%s\n",x);
#endif

#define MAXDEBUG 0
#if MAXDEBUG == 1
#define LOG(x) printf("%s\n",x)
#define PRINT5(var,offset) for(;print_count<5;print_count++) printf("var(%d)=%2.20f\n",print_count,var[offset+print_count]);
#define PRINT10(var) if (print_count<10) { printf("var=%1.20e\n",var); print_count++; }
#define PRINT10i(var) if (print_count<10) { printf("var=%d\n",var); print_count++; }
#else
#define LOG(x) // printf("%s\n",x);
#define PRINT5(var,offset) // for(i=0;i<10;i++) printf("var(%d)=%f\n",i,var[offset+i]);
#endif

// performance timers
#define CUDA_TIMING 0
#define CUDA_TIMING_UPDATE 0

// error checking after cuda function calls
// (note: this synchronizes many calls, thus e.g. no asynchronuous memcpy possible)
//#define ENABLE_VERY_SLOW_ERROR_CHECKING

// maximum function
#define MAX(x,y)                    (((x) < (y)) ? (y) : (x))

/* ----------------------------------------------------------------------------------------------- */

// cuda constant arrays

/* ----------------------------------------------------------------------------------------------- */

// dimensions
#define NDIM 2

// Gauss-Lobatto-Legendre
#define NGLLX 5
#define NGLL2 25

#define NGLL2_PADDED 32

// number of standard linear solids
#define N_SLS 3

/* ----------------------------------------------------------------------------------------------- */

//For Cuda aware MPI
#define ENV_LOCAL_RANK "OMPI_COMM_WORLD_LOCAL_RANK"

/* ----------------------------------------------------------------------------------------------- */

// Output paths, see setup/constants.h
#define OUTPUT_FILES_PATH "./OUTPUT_FILES/"

/* ----------------------------------------------------------------------------------------------- */

// macros for version output
#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)
#define VAR_NAME_VALUE(var) #var " = "  VALUE(var)

#pragma message ("Compiling with: " VAR_NAME_VALUE(CUDA_VERSION) "\n")
#if defined(__CUDA_ARCH__)
#pragma message ("Compiling with: " VAR_NAME_VALUE(__CUDA_ARCH__) "\n")
#endif

// Texture memory usage:
// requires CUDA version >= 4.0, see check below
// Use textures for d_displ and d_accel -- 10% performance boost
//#define USE_TEXTURES_FIELDS

// Using texture memory for the hprime-style constants is slower on
// Fermi generation hardware, but *may* be faster on Kepler
// generation.
// Use textures for hprime_xx
//#define USE_TEXTURES_CONSTANTS  // might not working properly yet, please test on your card...

// CUDA version >= 4.0 needed for cudaTextureType1D and cudaDeviceSynchronize()
#if CUDA_VERSION < 4000
#undef USE_TEXTURES_FIELDS
#undef USE_TEXTURES_CONSTANTS
#endif

#ifdef USE_TEXTURES_FIELDS
#pragma message ("Compiling with: USE_TEXTURES_FIELDS enabled\n")
#endif
#ifdef USE_TEXTURES_CONSTANTS
#pragma message ("Compiling with: USE_TEXTURES_CONSTANTS enabled\n")
#endif

// CUDA compiler specifications
// (optional) use launch_bounds specification to increase compiler optimization
//
// Kepler architecture
#ifdef GPU_DEVICE_K20
// (depending on GPU type, register spilling might slow down the performance)
// (uncomment if not desired)
#define USE_LAUNCH_BOUNDS
// elastic kernel
// note: main kernel is Kernel_2_***_impl() which is limited by shared memory usage to 8 active blocks
//       while register usage might use up to 9 blocks
//
// performance statistics: kernel Kernel_2_noatt_impl():
//       shared memory per block = 1700    for Kepler: total = 49152 -> limits active blocks to 16
//       registers per thread    = 48
//       registers per block     = 6144                total = 65536 -> limits active blocks to 10
//
// performance statistics: kernel Kernel_2_att_impl():
//       shared memory per block = 6100    for Kepler: total = 49152 -> limits active blocks to 8
//       registers per thread    = 59
//       registers per block     = 8192                total = 65536 -> limits active blocks to 8
#define LAUNCH_MIN_BLOCKS 10
// acoustic kernel
// performance statistics: kernel Kernel_2_acoustic_impl():
//       shared memory per block = 2200    for Kepler: -> limits active blocks to 16 (maximum possible)
//       registers per thread    = 40
//       registers per block     = 5120                -> limits active blocks to 12
// note: for K20x, using a minimum of 16 blocks leads to register spilling.
//       this slows down the kernel by ~ 4%
#define LAUNCH_MIN_BLOCKS_ACOUSTIC 16
#endif
//
// Pascal architecture
#ifdef GPU_DEVICE_Pascal
// Pascal P100: Pascal: total of 65536 register size
//              careful, as using launch bounds to increase the number of blocks might lead to register spilling.
#undef USE_LAUNCH_BOUNDS
#define LAUNCH_MIN_BLOCKS 10
#define LAUNCH_MIN_BLOCKS_ACOUSTIC 16
#endif

#ifdef USE_LAUNCH_BOUNDS
#pragma message ("\nCompiling with: USE_LAUNCH_BOUNDS enabled\n")
#endif

/* ----------------------------------------------------------------------------------------------- */

// cuda kernel block size for updating displacements/potential (newmark time scheme)
// current hardware: 128 is slightly faster than 256 ( ~ 4%)
#define BLOCKSIZE_KERNEL1 32
#define BLOCKSIZE_KERNEL3 128
#define BLOCKSIZE_TRANSFER 256

// maximum grid dimension in one direction of GPU
#define MAXIMUM_GRID_DIM 65535

/* ----------------------------------------------------------------------------------------------- */

// indexing
#define INDEX2(xsize,x,y) x + (y)*xsize
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))
#define INDEX6(xsize,ysize,zsize,isize,jsize,x,y,z,i,j,k) x + xsize*(y + ysize*(z + zsize*(i + isize*(j + jsize*k))))


#define INDEX3_PADDED(xsize,ysize,x,y,i) x + (y)*xsize + (i)*NGLL2_PADDED
#define INDEX4_PADDED(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*z) + (i)*NGLL3_PADDED

/* ----------------------------------------------------------------------------------------------- */

// custom type declarations

/* ----------------------------------------------------------------------------------------------- */

// type of "working" variables: see also CUSTOM_REAL
// double precision temporary variables leads to 10% performance decrease
// in Kernel_2_impl (not very much..)
typedef float realw;

// textures
typedef texture<float, cudaTextureType1D, cudaReadModeElementType> realw_texture;

// pointer declarations
// restricted pointers: can improve performance on Kepler ~ 10%
//   see: http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#restrict
//   however, compiler tends to use texture loads for restricted memory arrays, which might slow down performance
//
// non-restricted (default)
//typedef const realw* realw_const_p;
// restricted
typedef const realw* __restrict__ realw_const_p;
//
// non-restricted (default)
//typedef realw* realw_p;
// restricted
typedef realw* __restrict__ realw_p;

// wrapper for global memory load function
// usage:  val = get_global_cr( &A[index] );
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 350)
// Device has ldg
__device__ __forceinline__ realw get_global_cr(realw_const_p ptr) { return __ldg(ptr); }
#else
//Device does not, fall back.
__device__ __forceinline__ realw get_global_cr(realw_const_p ptr) { return (*ptr); }
#endif

/* ----------------------------------------------------------------------------------------------- */

// utility functions: defined in check_fields_cuda.cu

/* ----------------------------------------------------------------------------------------------- */

double get_time();
void get_free_memory(double* free_db, double* used_db, double* total_db);
void print_CUDA_error_if_any(cudaError_t err, int num);
void pause_for_debugger(int pause);
void exit_on_cuda_error(const char* kernel_name);
void exit_on_error(const char* info);
void synchronize_cuda();
void synchronize_mpi();
void start_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop);
void stop_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop,const char* info_str);
void stop_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop,const char* info_str,realw* t);
void get_blocks_xy(int num_blocks,int* num_blocks_x,int* num_blocks_y);
realw get_device_array_maximum_value(realw* array,int size);


/* ----------------------------------------------------------------------------------------------- */

// mesh pointer wrapper structure

/* ----------------------------------------------------------------------------------------------- */

typedef struct mesh_ {

  // mesh resolution
  int NSPEC_AB;
  int NGLOB_AB;

  // mpi process
  int myrank;

  // constants
  int simulation_type;
  int save_forward;
  int absorbing_conditions;

  // ------------------------------------------------------------------ //
  // GLL points & weights
  // ------------------------------------------------------------------ //

  // interpolators
  realw* d_xix; realw* d_xiz;
  realw* d_gammax; realw* d_gammaz;

  // model parameters
  realw* d_kappav; realw* d_muv;

  // global indexing
  int* d_ibool;

  // inner / outer elements
  int* d_ispec_is_inner;

  // pointers to constant memory arrays
  realw* d_hprime_xx;
  realw* d_hprimewgll_xx;
  realw* d_wxgll;

  // A buffer for mpi-send/recv, which is duplicated in fortran but is
  // allocated with pinned memory to facilitate asynchronus device <->
  // host memory transfers
  float* h_send_accel_buffer;
  float* h_send_b_accel_buffer;

  float* send_buffer;
  float* h_recv_accel_buffer;
  float* h_recv_b_accel_buffer;
  //float* recv_buffer;

  int size_mpi_buffer;
  int size_mpi_buffer_potential;

  // mpi interfaces
  int num_interfaces_ext_mesh;
  int max_nibool_interfaces_ext_mesh;
  int* d_nibool_interfaces_ext_mesh;
  int* d_ibool_interfaces_ext_mesh;

  // overlapped memcpy streams
  cudaStream_t compute_stream;
  cudaStream_t copy_stream;
  // stream to copy wavefield over several iterations
  cudaStream_t copy_stream_no_backward;
  //cudaStream_t b_copy_stream;

  // sources
  int nsources_local;
  realw* d_sourcearrays;
  int* d_ispec_selected_source;
  realw* d_source_time_function;

  // receivers
  int nrec_local;
  int* d_ispec_selected_rec_loc;
  realw* h_seismograms;
  realw* d_seismograms;
  realw* d_cosrot;
  realw* d_sinrot;

  // adjoint receivers/sources
  int nadj_rec_local;
  realw* d_adj_sourcearrays;
  realw* h_adj_sourcearrays_slice;
  realw* d_source_adjointe;
  realw* d_xir_store_loc;
  realw* d_gammar_store_loc;

  // ------------------------------------------------------------------ //
  // elastic wavefield parameters
  // ------------------------------------------------------------------ //

  // displacement, velocity, acceleration
  realw* d_displ; realw* d_veloc; realw* d_accel;
  // backward/reconstructed elastic wavefield
  realw* d_b_displ; realw* d_b_veloc; realw* d_b_accel;

  // elastic elements
  int nspec_elastic;
  int* d_ispec_is_elastic;

  // elastic domain parameters
  int* d_phase_ispec_inner_elastic;
  int num_phase_ispec_elastic;
  int ninterface_elastic;
  int * d_inum_interfaces_elastic;

  realw* d_rmassx;
  realw* d_rmassz;

  //attenuation
  realw* d_e1;
  realw* d_e11;
  realw* d_e13;
  realw* d_A_newmark_mu;
  realw* d_B_newmark_mu;
  realw* d_A_newmark_kappa;
  realw* d_B_newmark_kappa;
  realw* d_dux_dxl_old;
  realw* d_duz_dzl_old;
  realw* d_dux_dzl_plus_duz_dxl_old;

  // mpi buffer
  realw* d_send_accel_buffer;
  realw* d_b_send_accel_buffer;
  realw* d_recv_accel_buffer;
  realw* d_b_recv_accel_buffer;

  //used for absorbing stacey boundaries
  int d_num_abs_boundary_faces;
  int* d_abs_boundary_ispec;
  int* d_abs_boundary_ijk;
  realw* d_abs_boundary_normal;
  realw* d_abs_boundary_jacobian2Dw;
  int* d_cote_abs;
  int* d_ib_left;
  int* d_ib_right;
  int* d_ib_top;
  int* d_ib_bottom;
  realw* d_b_absorb_potential_bottom;
  realw* d_b_absorb_potential_left;
  realw* d_b_absorb_potential_right;
  realw* d_b_absorb_potential_top;
  int d_nspec_bottom;
  int d_nspec_left;
  int d_nspec_right;
  int d_nspec_top;
  realw* d_b_absorb_elastic_bottom;
  realw* d_b_absorb_elastic_left;
  realw* d_b_absorb_elastic_right;
  realw* d_b_absorb_elastic_top;

  realw* d_rho_vp;
  realw* d_rho_vs;

  // surface elements (to save for noise tomography and acoustic simulations)
  int* d_free_surface_ispec;
  int* d_free_surface_ijk;
  int num_free_surface_faces;

  // anisotropy
  realw* d_c11store;
  realw* d_c12store;
  realw* d_c13store;
  realw* d_c15store;
  realw* d_c23store;
  realw* d_c25store;
  realw* d_c33store;
  realw* d_c35store;
  realw* d_c55store;

  // sensitivity kernels
  realw* d_rho_kl;
  realw* d_mu_kl;
  realw* d_kappa_kl;
  realw* d_hess_el_kl;
  realw* d_dsxx;
  realw* d_dsxz;
  realw* d_dszz;
  realw* d_b_dsxx;
  realw* d_b_dsxz;
  realw* d_b_dszz;

  // JC JC here we will need to add GPU support for the new C-PML routines

  // ------------------------------------------------------------------ //
  // acoustic wavefield
  // ------------------------------------------------------------------ //
  // potential and first and second time derivative
  realw* d_potential_acoustic; realw* d_potential_dot_acoustic; realw* d_potential_dot_dot_acoustic;
  // backward/reconstructed wavefield
  realw* d_b_potential_acoustic; realw* d_b_potential_dot_acoustic; realw* d_b_potential_dot_dot_acoustic;
  // buffer for NO_BACKWARD_RECONSTRUCTION
  realw* d_potential_acoustic_buffer;

  cudaEvent_t transfer_is_complete1;
  cudaEvent_t transfer_is_complete2;

  // acoustic domain parameters
  int nspec_acoustic;
  int* d_ispec_is_acoustic;

  int* d_phase_ispec_inner_acoustic;
  int num_phase_ispec_acoustic;
  int ninterface_acoustic;
  int * d_inum_interfaces_acoustic;

  realw* d_rhostore;
  realw* d_kappastore;
  realw* d_rmass_acoustic;

  // attenuation
  realw* d_A_newmark_acous;
  realw* d_B_newmark_acous;
  realw* d_e1_acous;
  realw* d_sum_forces_old;
  realw* d_b_e1_acous;
  realw* d_b_sum_forces_old;

  // mpi buffer
  realw* d_send_potential_dot_dot_buffer;
  realw* d_b_send_potential_dot_dot_buffer;

  // sensitivity kernels
  realw* d_rho_ac_kl;
  realw* d_kappa_ac_kl;

  // approximative hessian for preconditioning kernels
  realw* d_hess_ac_kl;

  // coupling acoustic-elastic
  int* d_coupling_ac_el_ispec;
  int* d_coupling_ac_el_ijk;
  realw* d_coupling_ac_el_normal;
  realw* d_coupling_ac_el_jacobian2Dw;

} Mesh;


#endif
