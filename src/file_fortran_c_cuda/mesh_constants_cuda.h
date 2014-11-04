/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 1
 !               ---------------------------------------
 !
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                        Princeton University, USA
 !                and CNRS / University of Marseille, France
 !                 (there are currently many more authors!)
 ! (c) Princeton University and CNRS / University of Marseille, July 2012
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !=====================================================================
 */

/* trivia

- for most working arrays we use now "realw" instead of "float" type declarations to make it easier to switch
  between a real or double precision simulation
  (matching CUSTOM_REAL == 4 or 8 in fortran routines).

- instead of boolean "logical" declared in fortran routines, in C (or Cuda-C) we have to use "int" variables.
  ifort / gfortran caveat:
    to check whether it is true or false, do not check for == 1 to test for true values since ifort just uses
    non-zero values for true (e.g. can be -1 for true). however, false will be always == 0.
  thus, rather use: if( var ) {...}  for testing if true instead of if( var == 1){...} (alternative: one could use if( var != 0 ){...}

*/

#ifndef GPU_MESH_
#define GPU_MESH_

#include <sys/types.h>
#include <unistd.h>

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
#define PRINT10(var) if(print_count<10) { printf("var=%1.20e\n",var); print_count++; }
#define PRINT10i(var) if(print_count<10) { printf("var=%d\n",var); print_count++; }
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


//For Cuda aware MPI
#define ENV_LOCAL_RANK "OMPI_COMM_WORLD_LOCAL_RANK"
/* ----------------------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------------------------- */

// mesh pointer wrapper structure

/* ----------------------------------------------------------------------------------------------- */

typedef float realw;


typedef struct mesh_ {

  // mesh resolution
  int NSPEC_AB;
  int NGLOB_AB;

  // mpi process
  int myrank;

  // constants
  int simulation_type;
  int save_forward;
  int use_mesh_coloring_gpu;
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
  float* recv_buffer;

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
  //cudaStream_t b_copy_stream;

  // ------------------------------------------------------------------ //
  // elastic wavefield parameters
  // ------------------------------------------------------------------ //

  // displacement, velocity, acceleration
  realw* d_displ; realw* d_veloc; realw* d_accel;
  // backward/reconstructed elastic wavefield
  realw* d_b_displ; realw* d_b_veloc; realw* d_b_accel;

  // elastic elements
  int* d_ispec_is_elastic;

  // elastic domain parameters
  int* d_phase_ispec_inner_elastic;
  int num_phase_ispec_elastic;
  int ninterface_elastic;
  int * d_inum_interfaces_elastic;


  // mesh coloring
  int* h_num_elem_colors_elastic;
  int num_colors_outer_elastic,num_colors_inner_elastic;
  int nspec_elastic;

  realw* d_rmassx;
  realw* d_rmassz;

  // mpi buffer
  realw* d_send_accel_buffer;
  realw* d_b_send_accel_buffer;
  realw* d_recv_accel_buffer;
  realw* d_b_recv_accel_buffer;
  int* h_d_send_buffer;
  int* d_recv_buffer_ptr;

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

  // sources
  int nsources_local;
  realw* d_sourcearrays;
  int* d_islice_selected_source;
  int* d_ispec_selected_source;
  realw* d_source_time_function;
  int* d_num_src_loc;

  // receivers
  int* d_number_receiver_global;
  int* d_ispec_selected_rec;
  int nrec_local;

  realw* h_seismograms;
  realw* d_seismograms;
  realw* d_cosrot;
  realw* d_sinrot;

  // adjoint receivers/sources
  int nadj_rec_local;
  realw* d_adj_sourcearrays;
  realw* h_adj_sourcearrays_slice;
  int* d_pre_computed_irec;
  realw* d_source_adjointe;
  realw* d_xir_store_loc;
  realw* d_gammar_store_loc;

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

  // acoustic domain parameters
  int* d_ispec_is_acoustic;

  int* d_phase_ispec_inner_acoustic;
  int num_phase_ispec_acoustic;
  int ninterface_acoustic;
  int * d_inum_interfaces_acoustic;


  // mesh coloring
  int* h_num_elem_colors_acoustic;
  int num_colors_outer_acoustic,num_colors_inner_acoustic;
  int nspec_acoustic;

  realw* d_rhostore;
  realw* d_kappastore;
  realw* d_rmass_acoustic;

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
