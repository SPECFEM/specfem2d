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


#ifdef USE_TEXTURES_FIELDS
realw_texture d_potential_tex;
realw_texture d_potential_dot_dot_tex;
//backward/reconstructed
realw_texture d_b_potential_tex;
realw_texture d_b_potential_dot_dot_tex;

//note: texture variables are implicitly static, and cannot be passed as arguments to cuda kernels;
//      thus, 1) we thus use if-statements (FORWARD_OR_ADJOINT) to determine from which texture to fetch from
//            2) we use templates
//      since if-statements are a bit slower as the variable is only known at runtime, we use option 2)

// templates definitions
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_potential(int x);
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_potential_dot_dot(int x);

// templates for texture fetching
// FORWARD_OR_ADJOINT == 1 <- forward arrays
template<> __device__ float texfetch_potential<1>(int x) { return tex1Dfetch(d_potential_tex, x); }
template<> __device__ float texfetch_potential_dot_dot<1>(int x) { return tex1Dfetch(d_potential_dot_dot_tex, x); }
// FORWARD_OR_ADJOINT == 3 <- backward/reconstructed arrays
template<> __device__ float texfetch_potential<3>(int x) { return tex1Dfetch(d_b_potential_tex, x); }
template<> __device__ float texfetch_potential_dot_dot<3>(int x) { return tex1Dfetch(d_b_potential_dot_dot_tex, x); }

#endif

#ifdef USE_TEXTURES_CONSTANTS
// already defined in compute_forces_viscoelastic_cuda.cu
extern realw_texture d_hprime_xx_tex;
//extern realw_texture d_hprimewgll_xx_tex;
extern realw_texture d_wxgll_xx_tex;
#endif


// note on performance optimizations:
//
//   performance tests done:
//   - registers: we were trying to reduce the number of registers, as this is the main limiter for the
//                occupancy of the kernel. however, there is only little difference in register pressure for one "general" kernel
//                or multiple "spezialized" kernels. reducing registers is mainly achieved through the launch_bonds() directive.
//   - branching: we were trying to reduce code branches, such as the if-active check in earlier code versions.
//                reducing the branching helps the compiler to better optimize the executable.
//   - memory accesses: the global memory accesses are avoiding texture reads for coalescent arrays, as this is
//                still faster. thus we were using no __ldg() loads or __restricted__ pointer usage,
//                as those implicitly lead the compiler to use texture reads.
//   - arithmetic intensity: ratio of floating-point operations vs. memory accesses is still low for our kernels.
//                tests with using a loop over elements to re-use the constant arrays (like hprime, wgllwgll,..) and thus
//                increasing the arithmetic intensity failed because the number of registers increased as well.
//                this increased register pressure reduced the occupancy and slowed down the kernel performance.
//   - hiding memory latency: to minimize waiting times to retrieve a memory value from global memory, we put
//                some more calculations into the same code block before calling syncthreads(). this should help the
//                compiler to move independent calculations to wherever it can overlap it with memory access operations.
//                note, especially the if (gravity )-block locations are very sensitive
//                for optimal register usage and compiler optimizations
//

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2 - acoustic compute forces kernel

/* ----------------------------------------------------------------------------------------------- */

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS_ACOUSTIC)
#endif
Kernel_2_acoustic_impl(const int nb_blocks_to_compute,
                       const int* d_ibool,
                       const int* d_phase_ispec_inner_acoustic,
                       const int num_phase_ispec_acoustic,
                       const int d_iphase,
                       realw_const_p d_potential_acoustic,
                       realw_p d_potential_dot_dot_acoustic,
                       realw_const_p d_b_potential_acoustic,
                       realw_p d_b_potential_dot_dot_acoustic,
                       const int nb_field,
                       const realw* d_xix, const realw* d_xiz,
                       const realw* d_gammax,const realw* d_gammaz,
                       realw_const_p d_hprime_xx,
                       realw_const_p d_hprimewgll_xx,
                       realw_const_p d_wxgll,
                       const realw* d_rhostore){

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // thread-id == GLL node id
  // note: use only NGLL^2 = 25 active threads, plus 7 inactive/ghost threads,
  //       because we used memory padding from NGLL^2 = 25 to 32 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;

  int I,J;
  int iglob,offset;

  realw temp1l,temp3l;
  realw xixl,xizl,gammaxl,gammazl;

  realw dpotentialdxl,dpotentialdzl;
  realw rho_invl_times_jacobianl;

  realw sum_terms;

  __shared__ realw s_dummy_loc[2*NGLL2];

  __shared__ realw s_temp1[NGLL2];
  __shared__ realw s_temp3[NGLL2];

  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];
  __shared__ realw sh_wxgll[NGLLX];



// arithmetic intensity: ratio of number-of-arithmetic-operations / number-of-bytes-accessed-on-DRAM
//
// hand-counts on floating-point operations: counts addition/subtraction/multiplication/division
//                                           no counts for operations on indices in for-loops (compiler will likely unrool loops)
//
//                                           counts accesses to global memory, but no shared memory or register loads/stores
//                                           float has 4 bytes

//         counts floating-point operations (FLOP) per thread
//         counts global memory accesses in bytes (BYTES) per block
// 2 FLOP
//
// 0 BYTES

  // checks if anything to do
  if (bx >= nb_blocks_to_compute ) return;

// counts:
// + 1 FLOP
//
// + 0 BYTE

  // local padded index
  offset = (d_phase_ispec_inner_acoustic[bx + num_phase_ispec_acoustic*(d_iphase-1)]-1)*NGLL2_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1;


// counts:
// + 7 FLOP
//
// + 2 float * 32 threads = 256 BYTE

#ifdef USE_TEXTURES_FIELDS
  s_dummy_loc[tx] = texfetch_potential<FORWARD_OR_ADJOINT>(iglob);
  if (nb_field==2) s_dummy_loc[NGLL2+tx]=texfetch_potential<3>(iglob);
#else
  // changing iglob indexing to match fortran row changes fast style
  s_dummy_loc[tx] = d_potential_acoustic[iglob];
  if (nb_field==2) s_dummy_loc[NGLL2+tx]=d_b_potential_acoustic[iglob];
#endif


// counts:
// + 0 FLOP
//
// + 1 float * 25 threads = 100 BYTE

  // local index
  J = (tx/NGLLX);
  I = (tx-J*NGLLX);

// counts:
// + 3 FLOP
//
// + 0 BYTES

  // note: loads mesh values here to give compiler possibility to overlap memory fetches with some computations;
  //       arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads (arrays accesses are coalescent, thus no need for texture reads)
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] );
  xizl = d_xiz[offset];
  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  rho_invl_times_jacobianl = 1.f /(d_rhostore[offset] * (xixl*gammazl-gammaxl*xizl));

// counts:
// + 5 FLOP
//
// + 5 float * 32 threads = 160 BYTE

  // loads hprime into shared memory

#ifdef USE_TEXTURES_CONSTANTS
  sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_tex,tx);
#else
  sh_hprime_xx[tx] = d_hprime_xx[tx];
#endif
  // loads hprimewgll into shared memory
  sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];

  if (threadIdx.x < NGLLX){
#ifdef USE_TEXTURES_CONSTANTS
    sh_wxgll[tx] = tex1Dfetch(d_wxgll_xx_tex,tx);
#else
    // changing iglob indexing to match fortran row changes fast style
    sh_wxgll[tx] = d_wxgll[tx];
#endif
  }


// counts:
// + 0 FLOP
//
// + 2 * 1 float * 25 threads = 200 BYTE

  for (int k=0 ; k < nb_field ; k++) {

    // synchronize all the threads (one thread for each of the NGLL grid points of the
    // current spectral element) because we need the whole element to be ready in order
    // to be able to compute the matrix products along cut planes of the 3D element below
    __syncthreads();

    // computes first matrix product
    temp1l = 0.f;
    temp3l = 0.f;

    for (int l=0;l<NGLLX;l++) {

      //assumes that hprime_xx = hprime_yy = hprime_zz
      // 1. cut-plane along xi-direction
      temp1l += s_dummy_loc[NGLL2*k+J*NGLLX+l] * sh_hprime_xx[l*NGLLX+I];
      // 3. cut-plane along gamma-direction
      temp3l += s_dummy_loc[NGLL2*k+l*NGLLX+I] * sh_hprime_xx[l*NGLLX+J];
    }

// counts:
// + NGLLX * 2 * 6 FLOP = 60 FLOP
//
// + 0 BYTE

    // compute derivatives of ux, uy and uz with respect to x, y and z
    // derivatives of potential
    dpotentialdxl = xixl*temp1l +  gammaxl*temp3l;
    dpotentialdzl = xizl*temp1l +  gammazl*temp3l;

// counts:
// + 2 * 3 FLOP = 6 FLOP
//
// + 0 BYTE

    // form the dot product with the test vector
    s_temp1[tx] = sh_wxgll[J]*rho_invl_times_jacobianl  * (dpotentialdxl*xixl  + dpotentialdzl*xizl)  ;
    s_temp3[tx] = sh_wxgll[I]*rho_invl_times_jacobianl  * (dpotentialdxl*gammaxl + dpotentialdzl*gammazl)  ;

// counts:
// + 2 * 6 FLOP = 12 FLOP
//
// + 2 BYTE

    // synchronize all the threads (one thread for each of the NGLL grid points of the
    // current spectral element) because we need the whole element to be ready in order
    // to be able to compute the matrix products along cut planes of the 3D element below
    __syncthreads();

    sum_terms = 0.f;
    for (int l=0;l<NGLLX;l++) {
      //assumes hprimewgll_xx = hprimewgll_zz
      sum_terms -= s_temp1[J*NGLLX+l] * sh_hprimewgll_xx[I*NGLLX+l] + s_temp3[l*NGLLX+I] * sh_hprimewgll_xx[J*NGLLX+l];
    }

// counts:
// + NGLLX * 11 FLOP = 55 FLOP
//
// + 0 BYTE

    // assembles potential array
    if (k==0) {
      atomicAdd(&d_potential_dot_dot_acoustic[iglob],sum_terms);
    } else {
      atomicAdd(&d_b_potential_dot_dot_acoustic[iglob],sum_terms);
    }
// counts:
// + 1 FLOP
//
// + 1 float * 25 threads = 100 BYTE

// -----------------
// total of: 149 FLOP per thread
//           ~ 32 * 149 = 4768 FLOP per block
//
//           818 BYTE DRAM accesses per block
//
//           -> arithmetic intensity: 4768 FLOP / 818 BYTES ~ 5.83 FLOP/BYTE (hand-count)
  }
}

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2 - viscoacoustic compute forces kernel

/* ----------------------------------------------------------------------------------------------- */

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS_ACOUSTIC)
#endif
Kernel_2_viscoacoustic_impl(const int nb_blocks_to_compute,
                            const int* d_ibool,
                            const int* d_phase_ispec_inner_acoustic,
                            const int num_phase_ispec_acoustic,
                            const int d_iphase,
                            realw_const_p d_potential_acoustic,
                            realw_p d_potential_dot_dot_acoustic,
                            const realw* d_xix, const realw* d_xiz,
                            const realw* d_gammax,const realw* d_gammaz,
                            realw_const_p d_hprime_xx,
                            realw_const_p d_hprimewgll_xx,
                            realw_const_p d_wxgll,
                            const realw* d_rhostore,
                            realw_p d_e1_acous,
                            const realw* d_A_newmark,
                            const realw* d_B_newmark,
                            realw_p d_sum_forces_old){

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  int tx = threadIdx.x;
  int I,J;
  int iglob,offset,offset_align,i_sls;

  realw temp1l,temp3l;
  realw xixl,xizl,gammaxl,gammazl;
  realw dpotentialdxl,dpotentialdzl;
  realw rho_invl_times_jacobianl;
  realw sum_terms;
  realw sum_forces_old,forces_attenuation,a_newmark;
  realw e1_acous_load[N_SLS];

  __shared__ realw s_dummy_loc[NGLL2];
  __shared__ realw s_temp1[NGLL2];
  __shared__ realw s_temp3[NGLL2];
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];
  __shared__ realw sh_wxgll[NGLLX];

  if (bx >= nb_blocks_to_compute ) return;

  I =d_phase_ispec_inner_acoustic[bx + num_phase_ispec_acoustic*(d_iphase-1)]-1;
  offset = I*NGLL2_PADDED + tx;
  offset_align = I*NGLL2 + tx;
  iglob = d_ibool[offset] - 1;

#ifdef USE_TEXTURES_FIELDS
  s_dummy_loc[tx] = texfetch_potential<FORWARD_OR_ADJOINT>(iglob);
#else
  s_dummy_loc[tx] = d_potential_acoustic[iglob];
#endif

  // local index
  J = (tx/NGLLX);
  I = (tx-J*NGLLX);

  xixl = get_global_cr( &d_xix[offset] );
  xizl = d_xiz[offset];
  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  rho_invl_times_jacobianl = 1.f /(d_rhostore[offset] * (xixl*gammazl-gammaxl*xizl));

  for (i_sls=0;i_sls<N_SLS;i_sls++)  e1_acous_load[i_sls] = d_e1_acous[N_SLS*offset_align+i_sls];

#ifdef USE_TEXTURES_CONSTANTS
  sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_tex,tx);
#else
  sh_hprime_xx[tx] = d_hprime_xx[tx];
#endif
  // loads hprimewgll into shared memory
  sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];

  if (threadIdx.x < NGLLX){
#ifdef USE_TEXTURES_CONSTANTS
    sh_wxgll[tx] = tex1Dfetch(d_wxgll_xx_tex,tx);
#else
    sh_wxgll[tx] = d_wxgll[tx];
#endif
  }

  __syncthreads();

  // computes first matrix product
  temp1l = 0.f;
  temp3l = 0.f;

  for (int l=0;l<NGLLX;l++) {
    //assumes that hprime_xx = hprime_yy = hprime_zz
    // 1. cut-plane along xi-direction
    temp1l += s_dummy_loc[J*NGLLX+l] * sh_hprime_xx[l*NGLLX+I];
    // 3. cut-plane along gamma-direction
    temp3l += s_dummy_loc[l*NGLLX+I] * sh_hprime_xx[l*NGLLX+J];
  }

  dpotentialdxl = xixl*temp1l +  gammaxl*temp3l;
  dpotentialdzl = xizl*temp1l +  gammazl*temp3l;
  s_temp1[tx] = sh_wxgll[J]*rho_invl_times_jacobianl  * (dpotentialdxl*xixl  + dpotentialdzl*xizl)  ;
  s_temp3[tx] = sh_wxgll[I]*rho_invl_times_jacobianl  * (dpotentialdxl*gammaxl + dpotentialdzl*gammazl)  ;

  __syncthreads();

  sum_terms = 0.f;
  for (int l=0;l<NGLLX;l++) {
    //assumes hprimewgll_xx = hprimewgll_zz
    sum_terms -= s_temp1[J*NGLLX+l] * sh_hprimewgll_xx[I*NGLLX+l] + s_temp3[l*NGLLX+I] * sh_hprimewgll_xx[J*NGLLX+l];
  }

  sum_forces_old = d_sum_forces_old[offset_align];
  forces_attenuation = 0.f;

  for (i_sls=0;i_sls<N_SLS;i_sls++){
    a_newmark = d_A_newmark[N_SLS * offset_align + i_sls];
    e1_acous_load[i_sls] = a_newmark * a_newmark * e1_acous_load[i_sls] + d_B_newmark[N_SLS * offset_align + i_sls] * (sum_terms + a_newmark * sum_forces_old);
    forces_attenuation += e1_acous_load[i_sls];
    d_e1_acous[N_SLS*offset_align+i_sls] = e1_acous_load[i_sls];
  }

  d_sum_forces_old[offset_align] = sum_terms;
  sum_terms += forces_attenuation;

  atomicAdd(&d_potential_dot_dot_acoustic[iglob],sum_terms);
}




/* ----------------------------------------------------------------------------------------------- */

void Kernel_2_acoustic(int nb_blocks_to_compute, Mesh* mp, int d_iphase,
                       int* d_ibool,
                       realw* d_xix,realw* d_xiz,
                       realw* d_gammax,realw* d_gammaz,
                       realw* d_rhostore,
                       int ATTENUATION_VISCOACOUSTIC,
                       int compute_wavefield_1,
                       int compute_wavefield_2) {

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before acoustic kernel Kernel 2");
#endif

  // if the grid can handle the number of blocks, we let it be 1D
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y, nb_field;
  get_blocks_xy(nb_blocks_to_compute,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  cudaEvent_t start, stop;
  if (CUDA_TIMING) {
    start_timing_cuda(&start,&stop);
  }

  if (compute_wavefield_1 && compute_wavefield_2){
    nb_field=2;
  }else{
    nb_field=1;
  }
  if ( ! ATTENUATION_VISCOACOUSTIC){
    if (nb_field==2){
      // forward wavefields -> FORWARD_OR_ADJOINT == 1
      Kernel_2_acoustic_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                       d_ibool,
                                                                       mp->d_phase_ispec_inner_acoustic,
                                                                       mp->num_phase_ispec_acoustic,
                                                                       d_iphase,
                                                                       mp->d_potential_acoustic, mp->d_potential_dot_dot_acoustic,
                                                                       mp->d_b_potential_acoustic,mp->d_b_potential_dot_dot_acoustic,
                                                                       nb_field,
                                                                       d_xix, d_xiz,
                                                                       d_gammax, d_gammaz,
                                                                       mp->d_hprime_xx,
                                                                       mp->d_hprimewgll_xx,
                                                                       mp->d_wxgll,
                                                                       d_rhostore);
    }else{ // nb_field==1
      if (compute_wavefield_1){
        // forward wavefields -> FORWARD_OR_ADJOINT == 1
        Kernel_2_acoustic_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                         d_ibool,
                                                                         mp->d_phase_ispec_inner_acoustic,
                                                                         mp->num_phase_ispec_acoustic,
                                                                         d_iphase,
                                                                         mp->d_potential_acoustic, mp->d_potential_dot_dot_acoustic,
                                                                         mp->d_b_potential_acoustic,mp->d_b_potential_dot_dot_acoustic,
                                                                         nb_field,
                                                                         d_xix, d_xiz,
                                                                         d_gammax, d_gammaz,
                                                                         mp->d_hprime_xx,
                                                                         mp->d_hprimewgll_xx,
                                                                         mp->d_wxgll,
                                                                         d_rhostore);
      }
      if (compute_wavefield_2){
        // this run only happens with UNDO_ATTENUATION_AND_OR_PML on
        // adjoint wavefields -> FORWARD_OR_ADJOINT == 3
        Kernel_2_acoustic_impl<3><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                         d_ibool,
                                                                         mp->d_phase_ispec_inner_acoustic,
                                                                         mp->num_phase_ispec_acoustic,
                                                                         d_iphase,
                                                                         mp->d_b_potential_acoustic, mp->d_b_potential_dot_dot_acoustic,
                                                                         mp->d_b_potential_acoustic,mp->d_b_potential_dot_dot_acoustic,
                                                                         nb_field,
                                                                         d_xix, d_xiz,
                                                                         d_gammax, d_gammaz,
                                                                         mp->d_hprime_xx,
                                                                         mp->d_hprimewgll_xx,
                                                                         mp->d_wxgll,
                                                                         d_rhostore);
      } //compute_wavefield_1
    } //nb_field
  }else{ // ATTENUATION_VISCOACOUSTIC== .true. below
    if (compute_wavefield_1) {
      Kernel_2_viscoacoustic_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                            d_ibool,
                                                                            mp->d_phase_ispec_inner_acoustic,
                                                                            mp->num_phase_ispec_acoustic,
                                                                            d_iphase,
                                                                            mp->d_potential_acoustic, mp->d_potential_dot_dot_acoustic,
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
                                                                            mp->d_b_potential_acoustic, mp->d_b_potential_dot_dot_acoustic,
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

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("kernel Kernel_2");
#endif
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



/* ----------------------------------------------------------------------------------------------- */

/* KERNEL for enforce free surface */

/* ----------------------------------------------------------------------------------------------- */


__global__ void enforce_free_surface_cuda_kernel(realw_p potential_acoustic,
                                                 realw_p potential_dot_acoustic,
                                                 realw_p potential_dot_dot_acoustic,
                                                 const int num_free_surface_faces,
                                                 const int* free_surface_ispec,
                                                 const int* free_surface_ij,
                                                 const int* d_ibool,
                                                 const int* ispec_is_acoustic) {
  // gets spectral element face id
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  // for all faces on free surface
  if (iface < num_free_surface_faces) {

    int ispec = free_surface_ispec[iface]-1;

    // checks if element is in acoustic domain
    if (ispec_is_acoustic[ispec]) {

      // gets global point index
      int igll = threadIdx.x + threadIdx.y*blockDim.x;

      int i = free_surface_ij[INDEX3(NDIM,NGLLX,0,igll,iface)] - 1; // (1,igll,iface)
      int j = free_surface_ij[INDEX3(NDIM,NGLLX,1,igll,iface)] - 1;

      int iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

      // sets potentials to zero at free surface
      potential_acoustic[iglob] = 0.f;
      potential_dot_acoustic[iglob] = 0.f;
      potential_dot_dot_acoustic[iglob] = 0.f;
    }
  }
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(acoustic_enforce_free_surf_cuda,
              ACOUSTIC_ENFORCE_FREE_SURF_CUDA)(long* Mesh_pointer,int* compute_wavefield_1,int* compute_wavefield_2) {

  TRACE("acoustic_enforce_free_surf_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // does not absorb free surface, thus we enforce the potential to be zero at surface

  // checks if anything to do
  if (mp->num_free_surface_faces == 0) return;

  // block sizes
  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->num_free_surface_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLLX,1,1);


  // sets potentials to zero at free surface
  if (*compute_wavefield_1) {
  enforce_free_surface_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_acoustic,
                                                                          mp->d_potential_dot_acoustic,
                                                                          mp->d_potential_dot_dot_acoustic,
                                                                          mp->num_free_surface_faces,
                                                                          mp->d_free_surface_ispec,
                                                                          mp->d_free_surface_ijk,
                                                                          mp->d_ibool,
                                                                          mp->d_ispec_is_acoustic);
  }
  // for backward/reconstructed potentials
  if (*compute_wavefield_2) {
    enforce_free_surface_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_acoustic,
                                                                            mp->d_b_potential_dot_acoustic,
                                                                            mp->d_b_potential_dot_dot_acoustic,
                                                                            mp->num_free_surface_faces,
                                                                            mp->d_free_surface_ispec,
                                                                            mp->d_free_surface_ijk,
                                                                            mp->d_ibool,
                                                                            mp->d_ispec_is_acoustic);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("enforce_free_surface_cuda");
#endif
}


