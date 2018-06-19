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

realw_texture d_displ_tex;
realw_texture d_accel_tex;
// backward/reconstructed
realw_texture d_b_displ_tex;
realw_texture d_b_accel_tex;

//note: texture variables are implicitly static, and cannot be passed as arguments to cuda kernels;
//      thus, 1) we thus use if-statements (FORWARD_OR_ADJOINT) to determine from which texture to fetch from
//            2) we use templates
//      since if-statements are a bit slower as the variable is only known at runtime, we use option 2)

// templates definitions
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_displ(int x);
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_accel(int x);


// templates for texture fetching
// FORWARD_OR_ADJOINT == 1 <- forward arrays
template<> __device__ float texfetch_displ<1>(int x) { return tex1Dfetch(d_displ_tex, x); }
template<> __device__ float texfetch_accel<1>(int x) { return tex1Dfetch(d_accel_tex, x); }
// FORWARD_OR_ADJOINT == 3 <- backward/reconstructed arrays
template<> __device__ float texfetch_displ<3>(int x) { return tex1Dfetch(d_b_displ_tex, x); }
template<> __device__ float texfetch_accel<3>(int x) { return tex1Dfetch(d_b_accel_tex, x); }

#endif

#ifdef USE_TEXTURES_CONSTANTS
realw_texture d_hprime_xx_tex;
__constant__ size_t d_hprime_xx_tex_offset;
realw_texture d_wxgll_xx_tex;
__constant__ size_t d_wxgll_xx_tex_offset;
#endif


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2

/* ----------------------------------------------------------------------------------------------- */


// loads displacement into shared memory for element

template<int FORWARD_OR_ADJOINT>
__device__  __forceinline__ void load_shared_memory_displ(const int* tx, const int* iglob,
                                                          realw_const_p d_displ,
                                                          realw* sh_tempx,
                                                          realw* sh_tempz){

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^2 = 25 data points
#ifdef USE_TEXTURES_FIELDS
  sh_tempx[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*2);
  sh_tempz[(*tx)] = texfetch_displ<FORWARD_OR_ADJOINT>((*iglob)*2 + 1);
#else
  // changing iglob indexing to match fortran row changes fast style
  sh_tempx[(*tx)] = d_displ[(*iglob)*2];
  sh_tempz[(*tx)] = d_displ[(*iglob)*2 + 1];
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// loads hprime into shared memory for element

__device__  __forceinline__ void load_shared_memory_hprime(const int* tx,
                                                           realw_const_p d_hprime_xx,
                                                           realw* sh_hprime_xx){

  // each thread reads its corresponding value
  // (might be faster sometimes...)
#ifdef USE_TEXTURES_CONSTANTS
  // hprime
  sh_hprime_xx[(*tx)] = tex1Dfetch(d_hprime_xx_tex,(*tx) + d_hprime_xx_tex_offset);
#else
  // hprime
  sh_hprime_xx[(*tx)] = d_hprime_xx[(*tx)];
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// loads hprime into shared memory for element

__device__  __forceinline__ void load_shared_memory_wxgll(const int* tx,
                                                           realw_const_p d_wxgll,
                                                           realw* sh_wxgll){

  // each thread reads its corresponding value
  // (might be faster sometimes...)
#ifdef USE_TEXTURES_CONSTANTS
  // hprime
  sh_wxgll[(*tx)] = tex1Dfetch(d_wxgll_xx_tex,(*tx) + d_wxgll_xx_tex_offset);
#else
  // hprime
  sh_wxgll[(*tx)] = d_wxgll[(*tx)];
#endif
}




/* ----------------------------------------------------------------------------------------------- */

// loads hprimewgll into shared memory for element

__device__  __forceinline__ void load_shared_memory_hprimewgll(const int* tx,
                                                               realw_const_p d_hprimewgll_xx,
                                                               realw* sh_hprimewgll_xx) {

  // each thread reads its corresponding value
  // weighted hprime
//#ifdef USE_TEXTURES_CONSTANTS
  // hprime
//  sh_hprimewgll_xx[(*tx)] = tex1Dfetch(d_hprimewgll_xx_tex,(*tx));
//#else
  sh_hprimewgll_xx[(*tx)] = d_hprimewgll_xx[(*tx)];
//#endif
}

/* ----------------------------------------------------------------------------------------------- */



__device__  __forceinline__ void sum_hprime_xi(int I, int J,
                                              realw* tempxl,realw* tempzl,
                                              realw* sh_tempx,realw* sh_tempz, realw* sh_hprime) {

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumz = 0.f;

  // 1. cut-plane along xi-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprime[l*NGLLX+I];

    sumx += sh_tempx[J*NGLLX+l] * fac;
    sumz += sh_tempz[J*NGLLX+l] * fac;
  }

// counts:
// + NGLLX * ( 2 + 3*6 ) FLOP = 100 FLOP
//
// + 0 BYTE

  *tempxl = sumx;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */


__device__  __forceinline__ void sum_hprime_gamma(int I, int J,
                                                 realw* tempxl,realw* tempzl,
                                                 realw* sh_tempx,realw* sh_tempz, realw* sh_hprime) {

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumz = 0.f;

  // 3. cut-plane along gamma-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprime[l*NGLLX+J];

    sumx += sh_tempx[l*NGLLX+I] * fac;
    sumz += sh_tempz[l*NGLLX+I] * fac;
  }

  *tempxl = sumx;
  *tempzl = sumz;
}

/* ----------------------------------------------------------------------------------------------- */



__device__  __forceinline__ void sum_hprimewgll_xi(int I, int J,
                                                   realw* tempxl,realw* tempzl,
                                                   realw* sh_tempx,realw* sh_tempz, realw* sh_hprimewgll) {

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumz = 0.f;

  // 1. cut-plane along xi-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprimewgll[I*NGLLX+l]; //  d_hprimewgll_xx[I*NGLLX+l];

    sumx += sh_tempx[J*NGLLX+l] * fac;
    sumz += sh_tempz[J*NGLLX+l] * fac;
  }

  *tempxl = sumx;
  *tempzl = sumz;
}


/* ----------------------------------------------------------------------------------------------- */

// computes a 3D matrix-vector product along a 2D cut-plane

__device__  __forceinline__ void sum_hprimewgll_gamma(int I, int J,
                                                 realw* tempxl,realw* tempzl,
                                                 realw* sh_tempx,realw* sh_tempz, realw* sh_hprimewgll) {

  realw fac;

  // initializes
  realw sumx = 0.f;
  realw sumz = 0.f;

  // 3. cut-plane along gamma-direction
  #pragma unroll
  for (int l=0;l<NGLLX;l++) {
    fac = sh_hprimewgll[J*NGLLX+l]; // d_hprimewgll_xx[K*NGLLX+l];

    sumx += sh_tempx[l*NGLLX+I] * fac;
    sumz += sh_tempz[l*NGLLX+I] * fac;
  }

  *tempxl = sumx;
  *tempzl = sumz;
}


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2
//
// for elastic domains

/* ----------------------------------------------------------------------------------------------- */

// note:
// kernel_2 is split into 2 kernels:
//  - a kernel without attenuation and for isotropic media: Kernel_2_noatt_iso_impl()
//  - a kernel without attenuation and for anisotropic media: Kernel_2_noatt_ani_impl()
//
// this should help with performance:
// the high number of registers needed for our kernels limits the occupancy; separation tries to reduce this.


// kernel without attenuation
//
// we use templates to distinguish between calls with forward or adjoint texture fields

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_iso_impl(const int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        realw_const_p d_displ,
                        realw_p d_accel,
                        realw* d_xix,realw* d_xiz,
                        realw* d_gammax,realw* d_gammaz,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p wxgll,
                        realw* d_kappav,
                        realw* d_muv,
                        int simulation_type,
                        realw* dsxx,realw* dsxz,realw* dszz){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .false.
//  COMPUTE_AND_STORE_STRAIN  = .false.

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;

  int iglob,offset;
  int working_element;

  realw tempx1l,tempx3l,tempz1l,tempz3l;
  realw xixl,xizl,gammaxl,gammazl,jacobianl;
  realw duxdxl,duxdzl,duzdxl,duzdzl;
  realw duzdxl_plus_duxdzl;

  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_zz,sigma_xz;
  realw sum_terms1,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL2];
  __shared__ realw sh_tempz[NGLL2];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
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

// counts:
// 2 FLOP

  // checks if anything to do
  if (bx >= nb_blocks_to_compute ) return;


  // limits thread ids to range [0,25-1]
  if (tx >= NGLL2 ) tx = tx - NGLL2 ;

// counts:
// + 1 FLOP
//
// + 0 BYTE

  // loads hprime's into shared memory
  if (threadIdx.x < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);

    // copy hprimewgll from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }
  else if (threadIdx.x < NGLL2 + NGLLX ) load_shared_memory_wxgll(&tx,wxgll,sh_wxgll);

// counts:
// + 0 FLOP
//
// 2 * 1 float * 25 threads = 200 BYTE

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;

  // local padded index
  offset = working_element*NGLL2_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

// counts:
// + 7 FLOP
//
// + 2 float * 128 threads = 1024 BYTE

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^2 = 25 data points
  if (threadIdx.x < NGLL2) {
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempz);
  }

// counts:
// + 5 FLOP
//
// + 3 float * 125 threads = 1500 BYTE

  kappal = d_kappav[offset];
  mul = d_muv[offset];

// counts:
// + 0 FLOP
//
// + 2 * 1 float * 128 threads = 1024 BYTE

  // loads mesh values here to give compiler possibility to overlap memory fetches with some computations
  // note: arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads
  //       we only use the first loads explicitly by texture loads, all subsequent without. this should lead/trick
  //       the compiler to use global memory loads for all the subsequent accesses.
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xizl = get_global_cr( &d_xiz[offset] ); // first array with texture load

//  xixl = d_xix[offset]; // first array with texture load
//  xiyl = d_xiy[offset]; // all subsequent without to avoid over-use of texture for coalescent access
//  xizl = d_xiz[offset];

  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*gammazl-gammaxl*xizl);

// counts:
// + 15 FLOP
//
// + 9 float * 128 threads = 4608 BYTE

  // local index
  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

// counts:
// + 8 FLOP
//
// + 0 BYTE

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

 // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + gammaxl*tempx3l;
  duxdzl = xizl*tempx1l + gammazl*tempx3l;

  duzdxl = xixl*tempz1l + gammaxl*tempz3l;
  duzdzl = xizl*tempz1l + gammazl*tempz3l;

  // precompute some sums to save CPU time
  duzdxl_plus_duxdzl = duzdxl + duxdzl;

  // stress calculations

  // isotropic case
  // compute elements with an elastic isotropic rheology

  // note:
  // here, kappal and mul are taken from arrays kappastore and mustore,
  // while the CPU-routine takes values lambda and mu from poroelastcoef array
  //
  // conversion from kappa/mu to lambda/mu
  // AXISYM    : kappal = lambdal + TWO_THIRDS * mul
  // non-AXISYM: kappal = lambdal + mul

  // original
  //lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  //lambdal = lambdalplus2mul - 2.0f * mul;

  // new
  lambdal = kappal - mul;
  lambdalplus2mul = kappal + mul;

  // compute the three components of the stress tensor sigma

  sigma_xx = lambdalplus2mul*duxdxl + lambdal*duzdzl;
  sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl;
  sigma_xz = mul*duzdxl_plus_duxdzl;

// counts:
// + 22 FLOP
//
// + 0 BYTE

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  if (threadIdx.x < NGLL2) {
    sh_tempx[tx] = sh_wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_xz*xizl); // sh_tempx1
    sh_tempz[tx] = sh_wxgll[J] *jacobianl * (sigma_xz*xixl + sigma_zz*xizl); // sh_tempz1
  }
  __syncthreads();

  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);
  __syncthreads();

  if (threadIdx.x < NGLL2) {
    sh_tempx[tx] = sh_wxgll[I] * jacobianl * (sigma_xx*gammaxl +  sigma_xz*gammazl); // sh_tempx3
    sh_tempz[tx] = sh_wxgll[I] * jacobianl * (sigma_xz*gammaxl +  sigma_zz*gammazl); // sh_tempz3
  }
  __syncthreads();

  // 3. cut-plane gamma
  sum_hprimewgll_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprimewgll_xx);
  __syncthreads();

  sum_terms1= -tempx1l - tempx3l;
  sum_terms3= -tempz1l - tempz3l;

  // assembles acceleration array
  if (threadIdx.x < NGLL2) {
    atomicAdd(&d_accel[iglob*2], sum_terms1);
    atomicAdd(&d_accel[iglob*2+1], sum_terms3);
  }

// counts:
// + 8 FLOP
//
// + 3 float * 125 threads = 1500 BYTE

// Servira pour calcul futur des noyaux
  if (simulation_type == 3){
    dsxx[iglob] = duxdxl;
    dszz[iglob] = duzdzl;
    dsxz[iglob] = duzdxl_plus_duxdzl;
  }

// counts:
// -----------------
// total of: 790 FLOP per thread
//           ~ 128 * 790 = 101120 FLOP per block
//
//           11392 BYTE DRAM accesses per block
//
// arithmetic intensity: 101120 FLOP / 11392 BYTES ~ 8.9 FLOP/BYTE
// -----------------
//
// nvprof: nvprof --metrics flops_sp ./xspecfem3D
//          -> 883146240 FLOPS (Single) floating-point operations for 20736 elements
//          -> 42590 FLOP per block
// arithmetic intensity: 42590 FLOP / 11392 BYTES ~ 3.74 FLOP/BYTE
//
// roofline model: Kepler K20x
// ---------------------------
//   for a Kepler K20x card, the peak single-precision performance is about 3.95 TFlop/s.
//   global memory access has a bandwidth of ~ 250 GB/s.
//
//   memory bandwidth: 250 GB/s
//   single-precision peak performance: 3.95 TFlop/s -> corner arithmetic intensity = 3950./250. ~ 15.8 flop/byte
//
//   elastic kernel has an arithmetic intensity of: hand-counts   ~ 8.9 flop/byte
//                                                  nvprof-counts ~ 42590./11392. flop/byte = 3.74 flop/byte
//
//   -> we can only achieve about: (hand-counts)   56% of the peak performance
//                                 (nvprof-counts) 24% of the peak performance -> 935.0 GFlop/s
//
// roofline model: Tesla K20c (Kepler architecture: http://www.nvidia.com/content/tesla/pdf/Tesla-KSeries-Overview-LR.pdf)
// ---------------------------
//   memory bandwidth: 208 GB/s
//   single-precision peak performance: 3.52 TFlop/s -> corner arithmetic intensity = 3520 / 208 ~ 16.9 flop/byte
//
//   we can only achieve about: (hand-counts)   52% of the peak performance
//                              (nvprof-counts) 22% of the peak performance -> 779.0 GFlop/s - measured: 647.3 GFlop/s


} // kernel_2_noatt_iso_impl()

/* ----------------------------------------------------------------------------------------------- */


template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_noatt_ani_impl(int nb_blocks_to_compute,
                        const int* d_ibool,
                        const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                        const int d_iphase,
                        realw_const_p d_displ,
                        realw_p d_accel,
                        realw* d_xix,realw* d_xiz,
                        realw* d_gammax,realw* d_gammaz,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p wxgll,
                        realw_const_p d_kappav,
                        realw_const_p d_muv,
                        const int SIMULATION_TYPE,
                        const int ANISOTROPY,
                        realw* d_c11store,realw* d_c12store,realw* d_c13store,
                        realw* d_c15store,
                        realw* d_c23store,
                        realw* d_c25store,realw* d_c33store,
                        realw* d_c35store,
                        realw* d_c55store) {

// elastic compute kernel without attenuation for anisotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .true.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if (bx >= nb_blocks_to_compute ) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if (tx >= NGLL2 ) tx = NGLL2-1;

  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  int iglob,offset;
  int working_element;

  realw tempx1l,tempx3l,tempz1l,tempz3l;
  realw xixl,xizl,gammaxl,gammazl,jacobianl;
  realw duxdxl,duxdzl,duzdxl,duzdzl;
  realw duzdxl_plus_duxdzl;

  realw lambdal,mul,lambdalplus2mul,kappal;

  realw sigma_xx,sigma_zz,sigma_xz,sigma_zx;

  realw c11,c13,c15,c33,c35,c55;
  realw sum_terms1,sum_terms3;

  // shared memory
  __shared__ realw sh_tempx[NGLL2];
  __shared__ realw sh_tempz[NGLL2];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime's into shared memory
  if (tx < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;

  // local padded index
  offset = working_element*NGLL2_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (threadIdx.x < NGLL2) {
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempz);
  }

  // loads mesh values here to give compiler possibility to overlap memory fetches with some computations
  // note: arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads
  //       we only use the first loads explicitly by texture loads, all subsequent without. this should lead/trick
  //       the compiler to use global memory loads for all the subsequent accesses.
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
                              // all subsequent without to avoid over-use of texture for coalescent access
  xizl = d_xiz[offset];

  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*gammazl-gammaxl*xizl);

  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + gammaxl*tempx3l;
  duxdzl = xizl*tempx1l + gammazl*tempx3l;

  duzdxl = xixl*tempz1l + gammaxl*tempz3l;
  duzdzl = xizl*tempz1l + gammazl*tempz3l;

  // precompute some sums to save CPU time
  duzdxl_plus_duxdzl = duzdxl + duxdzl;

  // full anisotropic case, stress calculations
  if (ANISOTROPY){
    c11 = d_c11store[offset];
    c13 = d_c13store[offset];
    c15 = d_c15store[offset];
    c33 = d_c33store[offset];
    c35 = d_c35store[offset];
    c55 = d_c55store[offset];

    sigma_xx = c11*duxdxl + c15*duzdxl_plus_duxdzl + c13*duzdzl;
    sigma_zz = c13*duxdxl + c35*duzdxl_plus_duxdzl + c33*duzdzl;
    sigma_xz = c15*duxdxl + c55*duzdxl_plus_duxdzl + c35*duzdzl;
    sigma_zx = sigma_xz;
  }else{
    // isotropic case

    // compute elements with an elastic isotropic rheology
    kappal = d_kappav[offset];
    mul = d_muv[offset];

    lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
    lambdal = lambdalplus2mul - 2.0f * mul;

    // compute the three components of the stress tensor sigma
    sigma_xx = lambdalplus2mul*duxdxl + lambdal*duzdzl;
    sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl;
    sigma_xz = mul*duzdxl_plus_duxdzl;
    sigma_zx = sigma_xz;
  }

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  if (threadIdx.x < NGLL2) {
    sh_tempx[tx] = wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_zx*xizl); // sh_tempx1
    sh_tempz[tx] = wxgll[J] *jacobianl * (sigma_xz*xixl + sigma_zz*xizl); // sh_tempz1
  }
  __syncthreads();

  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  if (threadIdx.x < NGLL2) {
    sh_tempx[tx] = wxgll[I] * jacobianl * (sigma_xx*gammaxl +  sigma_zx*gammazl); // sh_tempx3
    sh_tempz[tx] = wxgll[I] * jacobianl * (sigma_xz*gammaxl +  sigma_zz*gammazl); // sh_tempz3
  }
  __syncthreads();

  // 3. cut-plane gamma
  sum_hprimewgll_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprimewgll_xx);
  __syncthreads();

  sum_terms1 = - tempx1l - tempx3l;
  sum_terms3 = - tempz1l - tempz3l;

  // assembles acceleration array
  if (threadIdx.x < NGLL2) {
    atomicAdd(&d_accel[iglob*2], sum_terms1);
    atomicAdd(&d_accel[iglob*2+1], sum_terms3);
  } // threadIdx.x

} // kernel_2_noatt_ani_impl()



/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2
//
// for viscoelastic domains

/* ----------------------------------------------------------------------------------------------- */

template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_att_iso_impl(const int nb_blocks_to_compute,
                      const int* d_ibool,
                      const int* d_phase_ispec_inner_elastic,const int num_phase_ispec_elastic,
                      const int d_iphase,
                      realw_const_p d_displ,
                      realw_p d_accel,
                      realw* d_xix,realw* d_xiz,
                      realw* d_gammax,realw* d_gammaz,
                      realw_const_p d_hprime_xx,
                      realw_const_p d_hprimewgll_xx,
                      realw_const_p wxgll,
                      realw* d_kappav,
                      realw* d_muv,
                      int simulation_type,
                      realw* dsxx,realw* dsxz,realw* dszz,
                      realw_const_p A_newmark_mu,realw_const_p B_newmark_mu,
                      realw_const_p A_newmark_kappa,realw_const_p B_newmark_kappa,
                      realw_p e1,realw_p e11,realw_p e13,
                      realw_p dux_dxl_old,realw_p duz_dzl_old,realw_p dux_dzl_plus_duz_dxl_old){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .true.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .false.
//  COMPUTE_AND_STORE_STRAIN  = .false.

  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  int tx = threadIdx.x;

  int iglob,offset,offset_align,i_sls;
  int working_element;

  realw tempx1l,tempx3l,tempz1l,tempz3l;
  realw xixl,xizl,gammaxl,gammazl,jacobianl;
  realw duxdxl,duxdzl,duzdxl,duzdzl;
  realw duzdxl_plus_duxdzl,duxdxl_plus_duzdzl;
  realw duxdxl_old,duzdzl_old,duxdzl_plus_duzdxl_old,duxdxl_plus_duzdzl_old;

  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_zz,sigma_xz;
  realw sum_terms1,sum_terms3;
  realw e1_load[N_SLS],e11_load[N_SLS],e13_load[N_SLS];
  realw e1_sum,e11_sum,e13_sum,a_newmark,b_newmark;

  // shared memory
  __shared__ realw sh_tempx[NGLL2];
  __shared__ realw sh_tempz[NGLL2];
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];
  __shared__ realw sh_wxgll[NGLLX];

  // checks if anything to do
  if (bx >= nb_blocks_to_compute ) return;


  // limits thread ids to range [0,25-1]
  if (tx >= NGLL2 ) tx = tx - NGLL2 ;

  // loads hprime's into shared memory
  if (threadIdx.x < NGLL2) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);

    // copy hprimewgll from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }
  else if (threadIdx.x < NGLL2 + NGLLX ) load_shared_memory_wxgll(&tx,wxgll,sh_wxgll);

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
  working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;

  // local padded index
  offset = working_element*NGLL2_PADDED + tx;
  offset_align = working_element*NGLL2 + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^2 = 25 data points
  if (threadIdx.x < NGLL2) {
    // copy displacement from global memory to shared memory
    load_shared_memory_displ<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,sh_tempx,sh_tempz);
  }

  kappal = d_kappav[offset];
  mul = d_muv[offset];

  for (i_sls=0;i_sls<N_SLS;i_sls++){
    e1_load[i_sls] = e1[N_SLS*offset_align+i_sls];
    e11_load[i_sls] = e11[N_SLS*offset_align+i_sls];
    e13_load[i_sls] = e13[N_SLS*offset_align+i_sls];
  }

  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xizl = get_global_cr( &d_xiz[offset] ); // first array with texture load
  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*gammazl-gammaxl*xizl);

  // local index
  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  __syncthreads();

 // computes first matrix products
  // 1. cut-plane
  sum_hprime_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprime_xx);
  // 3. cut-plane
  sum_hprime_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprime_xx);

  // compute derivatives of ux, uy and uz with respect to x, y and z
  duxdxl = xixl*tempx1l + gammaxl*tempx3l;
  duxdzl = xizl*tempx1l + gammazl*tempx3l;

  duzdxl = xixl*tempz1l + gammaxl*tempz3l;
  duzdzl = xizl*tempz1l + gammazl*tempz3l;

  // precompute some sums to save CPU time
  duzdxl_plus_duxdzl = duzdxl + duxdzl;

  // new
  lambdal = kappal - mul;
  lambdalplus2mul = kappal + mul;

  // compute the three components of the stress tensor sigma
  sigma_xx = lambdalplus2mul*duxdxl + lambdal*duzdzl;
  sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl;
  sigma_xz = mul*duzdxl_plus_duxdzl;

  //get the contribution of attenuation and update the memory variables
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duxdxl_old = dux_dxl_old[offset_align];
  duzdzl_old = duz_dzl_old[offset_align];
  duxdxl_plus_duzdzl_old = duxdxl_old + duzdzl_old;
  duxdzl_plus_duzdxl_old = dux_dzl_plus_duz_dxl_old[offset_align];

  e1_sum = 0.f;
  e11_sum = 0.f;
  e13_sum = 0.f;

  for (i_sls=0;i_sls<N_SLS;i_sls++){
    a_newmark = A_newmark_kappa[N_SLS * offset_align + i_sls];
    b_newmark = B_newmark_kappa[N_SLS * offset_align + i_sls];

    e1_load[i_sls] = a_newmark * a_newmark * e1_load[i_sls] + b_newmark * (duxdxl_plus_duzdzl + a_newmark * (duxdxl_plus_duzdzl_old));
    e1_sum += e1_load[i_sls];
    e1[N_SLS*offset_align+i_sls] = e1_load[i_sls];

    a_newmark = A_newmark_mu[N_SLS * offset_align + i_sls];
    b_newmark = B_newmark_mu[N_SLS * offset_align + i_sls];

    e11_load[i_sls] = a_newmark * a_newmark * e11_load[i_sls] + b_newmark * (duxdxl - 0.5f*duxdxl_plus_duzdzl + a_newmark * (duxdxl_old-0.5f*duxdxl_plus_duzdzl_old));
    e11_sum += e11_load[i_sls];
    e11[N_SLS*offset_align+i_sls] = e11_load[i_sls];

    e13_load[i_sls] = a_newmark * a_newmark * e13_load[i_sls] + b_newmark * (duzdxl_plus_duxdzl + a_newmark * duxdzl_plus_duzdxl_old);
    e13_sum += e13_load[i_sls];
    e13[N_SLS*offset_align+i_sls] = e13_load[i_sls];
  }

  // add the contribution of the attenuation
  sigma_xx += (lambdalplus2mul-mul) * e1_sum + 2.0f * mul * e11_sum;
  sigma_xz += mul * e13_sum;
  sigma_zz += (lambdalplus2mul-mul) * e1_sum - 2.0f * mul * e11_sum;

  // saves the grad(displ) to use at the next iteration
  dux_dxl_old[offset_align] = duxdxl;
  duz_dzl_old[offset_align] = duzdzl;
  dux_dzl_plus_duz_dxl_old[offset_align] = duzdxl_plus_duxdzl;

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  if (threadIdx.x < NGLL2) {
    sh_tempx[tx] = sh_wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_xz*xizl); // sh_tempx1
    sh_tempz[tx] = sh_wxgll[J] *jacobianl * (sigma_xz*xixl + sigma_zz*xizl); // sh_tempz1
  }
  __syncthreads();

  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);
  __syncthreads();

  if (threadIdx.x < NGLL2) {
    sh_tempx[tx] = sh_wxgll[I] * jacobianl * (sigma_xx*gammaxl +  sigma_xz*gammazl); // sh_tempx3
    sh_tempz[tx] = sh_wxgll[I] * jacobianl * (sigma_xz*gammaxl +  sigma_zz*gammazl); // sh_tempz3
  }
  __syncthreads();

  // 3. cut-plane gamma
  sum_hprimewgll_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprimewgll_xx);
  __syncthreads();

  sum_terms1= -tempx1l - tempx3l;
  sum_terms3= -tempz1l - tempz3l;

  // assembles acceleration array
  if (threadIdx.x < NGLL2) {
    atomicAdd(&d_accel[iglob*2], sum_terms1);
    atomicAdd(&d_accel[iglob*2+1], sum_terms3);
  }

// Servira pour calcul futur des noyaux
  if (simulation_type == 3){
    dsxx[iglob] = duxdxl;
    dszz[iglob] = duzdzl;
    dsxz[iglob] = duzdxl_plus_duxdzl;
  }

} // kernel_2_att_iso_impl()

/* ----------------------------------------------------------------------------------------------- */






/* ----------------------------------------------------------------------------------------------- */


void Kernel_2(int nb_blocks_to_compute,Mesh* mp,int d_iphase,realw d_deltat,
              int ANISOTROPY,int ATTENUATION_VISCOELASTIC,
              int* d_ibool,
              realw* d_xix,realw* d_xiz,
              realw* d_gammax,realw* d_gammaz,
              realw* d_kappav,
              realw* d_muv,
              realw* d_c11store,realw* d_c12store,realw* d_c13store,
              realw* d_c15store,realw* d_c23store,realw* d_c25store,
              realw* d_c33store,realw* d_c35store,realw* d_c55store) {

  TRACE("\tKernel_2");

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before kernel Kernel 2");
#endif

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

  // compute kernels without attenuation
  if (ANISOTROPY) {
    // full anisotropy
    // forward wavefields -> FORWARD_OR_ADJOINT == 1
    Kernel_2_noatt_ani_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                      d_ibool,
                                                                      mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                      d_iphase,
                                                                      mp->d_displ,
                                                                      mp->d_accel,
                                                                      d_xix, d_xiz,
                                                                      d_gammax, d_gammaz,
                                                                      mp->d_hprime_xx,
                                                                      mp->d_hprimewgll_xx,
                                                                      mp->d_wxgll,
                                                                      d_kappav,
                                                                      d_muv,
                                                                      mp->simulation_type,
                                                                      ANISOTROPY,
                                                                      d_c11store,d_c12store,d_c13store,
                                                                      d_c15store,
                                                                      d_c23store,
                                                                      d_c25store,d_c33store,
                                                                      d_c35store,
                                                                      d_c55store);

    // backward/reconstructed wavefield
    if (mp->simulation_type == 3) {
      // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
      Kernel_2_noatt_ani_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                      d_ibool,
                                                                      mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                      d_iphase,
                                                                      mp->d_b_displ,
                                                                      mp->d_b_accel,
                                                                      d_xix, d_xiz,
                                                                      d_gammax, d_gammaz,
                                                                      mp->d_hprime_xx,
                                                                      mp->d_hprimewgll_xx,
                                                                      mp->d_wxgll,
                                                                      d_kappav,
                                                                      d_muv,
                                                                      mp->simulation_type,
                                                                      ANISOTROPY,
                                                                      d_c11store,d_c12store,d_c13store,
                                                                      d_c15store,
                                                                      d_c23store,
                                                                      d_c25store,d_c33store,
                                                                      d_c35store,
                                                                      d_c55store);
    }
  }else{
    if (ATTENUATION_VISCOELASTIC){
      Kernel_2_att_iso_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                      d_ibool,
                                                                      mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                      d_iphase,
                                                                      mp->d_displ,
                                                                      mp->d_accel,
                                                                      d_xix, d_xiz,
                                                                      d_gammax, d_gammaz,
                                                                      mp->d_hprime_xx,
                                                                      mp->d_hprimewgll_xx,
                                                                      mp->d_wxgll,
                                                                      d_kappav,
                                                                      d_muv,
                                                                      mp->simulation_type,
                                                                      mp->d_dsxx,
                                                                      mp->d_dsxz,
                                                                      mp->d_dszz,
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
    else{
    // without storing strains
    // forward wavefields -> FORWARD_OR_ADJOINT == 1
    Kernel_2_noatt_iso_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                      d_ibool,
                                                                      mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                      d_iphase,
                                                                      mp->d_displ,
                                                                      mp->d_accel,
                                                                      d_xix, d_xiz,
                                                                      d_gammax, d_gammaz,
                                                                      mp->d_hprime_xx,
                                                                      mp->d_hprimewgll_xx,
                                                                      mp->d_wxgll,
                                                                      d_kappav,
                                                                      d_muv,
                                                                      mp->simulation_type,
                                                                      mp->d_dsxx,
                                                                      mp->d_dsxz,
                                                                      mp->d_dszz);

    }
    // backward/reconstructed wavefield
    if (mp->simulation_type == 3) {
      // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
      Kernel_2_noatt_iso_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                         d_ibool,
                                                                         mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                         d_iphase,
                                                                         mp->d_b_displ,
                                                                         mp->d_b_accel,
                                                                         d_xix, d_xiz,
                                                                         d_gammax,d_gammaz,
                                                                         mp->d_hprime_xx,
                                                                         mp->d_hprimewgll_xx,
                                                                         mp->d_wxgll,
                                                                         d_kappav,
                                                                         d_muv,
                                                                         mp->simulation_type,
                                                                         mp->d_b_dsxx,
                                                                         mp->d_b_dsxz,
                                                                         mp->d_b_dszz);
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

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("Kernel_2_impl");
#endif
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
                                                int* ATTENUATION_VISCOELASTIC) {

  TRACE("compute_forces_viscoelastic_cuda");
  // EPIK_TRACER("compute_forces_viscoelastic_cuda");

  //printf("Running compute_forces\n");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int num_elements;

  if (*iphase == 1)
    num_elements = *nspec_outer_elastic;
  else
    num_elements = *nspec_inner_elastic;

  // checks if anything to do
  if (num_elements == 0) return;

  // no mesh coloring: uses atomic updates
  Kernel_2(num_elements,mp,*iphase,*deltat,*ANISOTROPY,*ATTENUATION_VISCOELASTIC,
           mp->d_ibool,
           mp->d_xix,mp->d_xiz,
           mp->d_gammax,mp->d_gammaz,
           mp->d_kappav,
           mp->d_muv,
           mp->d_c11store,mp->d_c12store,mp->d_c13store,
           mp->d_c15store,mp->d_c23store,mp->d_c25store,
           mp->d_c33store,mp->d_c35store,mp->d_c55store);
}
