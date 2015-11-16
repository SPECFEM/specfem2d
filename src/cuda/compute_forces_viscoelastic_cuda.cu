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
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

*/
#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
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
//realw_texture d_hprimewgll_xx_tex;
realw_texture d_wxgll_xx_tex;
#endif


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2

/* ----------------------------------------------------------------------------------------------- */


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
  sh_hprime_xx[(*tx)] = tex1Dfetch(d_hprime_xx_tex,(*tx));
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
  sh_wxgll[(*tx)] = tex1Dfetch(d_wxgll_xx_tex,(*tx));
#else
  // hprime
  sh_wxgll[(*tx)] = d_wxgll[(*tx)];
#endif
}




/* ----------------------------------------------------------------------------------------------- */

// loads hprimewgll into shared memory for element

__device__  __forceinline__ void load_shared_memory_hprimewgll(const int* tx,
                                                               realw_const_p d_hprimewgll_xx,
                                                               realw* sh_hprimewgll_xx ){

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
                                              realw* sh_tempx,realw* sh_tempz, realw* sh_hprime ){

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
                                                 realw* sh_tempx,realw* sh_tempz, realw* sh_hprime ){

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
                                                   realw* sh_tempx,realw* sh_tempz, realw* sh_hprimewgll ){

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
                                                 realw* sh_tempx,realw* sh_tempz, realw* sh_hprimewgll ){

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
                        realw* d_kappav,realw* d_muv, int simulation_type,
                        realw* dsxx,realw* dsxz,realw* dszz){

// elastic compute kernel without attenuation for isotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .false.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)
//  gravity                   = .false.
//  use_mesh_coloring_gpu     = .false.
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
  if( bx >= nb_blocks_to_compute ) return;


  // limits thread ids to range [0,25-1]
  if( tx >= NGLL2 ) tx = tx - NGLL2 ;


// counts:
// + 1 FLOP
//
// + 0 BYTE

  // loads hprime's into shared memory
  if( threadIdx.x < NGLL2 ) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);

    // copy hprimewgll from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }
  else if ( threadIdx.x < NGLL2 + NGLLX ) load_shared_memory_wxgll(&tx,wxgll,sh_wxgll);




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
  if( threadIdx.x < NGLL2 ){
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

  lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  lambdal = lambdalplus2mul - 2.0f * mul;

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
  if( threadIdx.x < NGLL2 ) {
    sh_tempx[tx] = sh_wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_xz*xizl); // sh_tempx1
    sh_tempz[tx] = sh_wxgll[J] *jacobianl * (sigma_xz*xixl + sigma_zz*xizl); // sh_tempz1
  }
  __syncthreads();
  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);


  __syncthreads();
  if( threadIdx.x < NGLL2 ) {
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
  if(threadIdx.x < NGLL2) {
    atomicAdd(&d_accel[iglob*2], sum_terms1);
    atomicAdd(&d_accel[iglob*2+1], sum_terms3);
  }


// counts:
// + 8 FLOP
//
// + 3 float * 125 threads = 1500 BYTE



// Servira pour calcul futur des noyaux
if (simulation_type ==3){
dsxx[iglob]=duxdxl;
dszz[iglob]=duzdzl;
dsxz[iglob]=duzdxl_plus_duxdzl;
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
                        const int use_mesh_coloring_gpu,
                        realw_const_p d_displ,
                        realw_p d_accel,
                        realw* d_xix,realw* d_xiz,
                        realw* d_gammax,realw* d_gammaz,
                        realw_const_p d_hprime_xx,
                        realw_const_p d_hprimewgll_xx,
                        realw_const_p wxgll,
                        realw_const_p d_kappav,realw_const_p d_muv,
                        const int SIMULATION_TYPE,
                        const int ANISOTROPY,
                        realw* d_c11store,realw* d_c12store,realw* d_c13store,
                        realw* d_c15store,
                        realw* d_c23store,
                        realw* d_c25store,realw* d_c33store,
                        realw* d_c35store,
                        realw* d_c55store ){

// elastic compute kernel without attenuation for anisotropic elements
//
// holds for:
//  ATTENUATION               = .false.
//  ANISOTROPY                = .true.
//  COMPUTE_AND_STORE_STRAIN  = .true. or .false. (true for kernel simulations)


  // block-id == number of local element id in phase_ispec array
  int bx = blockIdx.y*gridDim.x+blockIdx.x;

  // checks if anything to do
  if( bx >= nb_blocks_to_compute ) return;

  // thread-id == GLL node id
  // note: use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  //       because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses;
  //       to avoid execution branching and the need of registers to store an active state variable,
  //       the thread ids are put in valid range
  int tx = threadIdx.x;
  if( tx >= NGLL2 ) tx = NGLL2-1;

  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  int iglob,offset;
  int working_element;

  realw tempx1l,tempx3l,tempz1l,tempz3l;
  realw xixl,xizl,gammaxl,gammazl,jacobianl;
  realw duxdxl,duxdzl,duzdxl,duzdzl;
  realw duzdxl_plus_duxdzl;

  realw lambdal,mul,lambdalplus2mul,kappal;

  realw sigma_xx,sigma_zz,sigma_xz;

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
  if( tx < NGLL2 ) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

  // spectral-element id
  // iphase-1 and working_element-1 for Fortran->C array conventions
#ifdef USE_MESH_COLORING_GPU
  working_element = bx;
#else
  //mesh coloring
  if( use_mesh_coloring_gpu ){
    working_element = bx;
  }else{
    // iphase-1 and working_element-1 for Fortran->C array conventions
    working_element = d_phase_ispec_inner_elastic[bx + num_phase_ispec_elastic*(d_iphase-1)] - 1;
  }
#endif
  // local padded index
  offset = working_element*NGLL2_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1 ;

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if( threadIdx.x < NGLL2 ){
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
  if(ANISOTROPY){
    c11 = d_c11store[offset];
    c13 = d_c13store[offset];
    c15 = d_c15store[offset];
    c33 = d_c33store[offset];
    c35 = d_c35store[offset];
    c55 = d_c55store[offset];


    sigma_xx = c11*duxdxl + c15*duzdxl_plus_duxdzl + c13*duzdzl;
    sigma_zz = c13*duxdxl + c35*duzdxl_plus_duxdzl + c33*duzdzl;
    sigma_xz = c15*duxdxl + c55*duzdxl_plus_duxdzl + c35*duzdzl;

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

  }


  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  if( threadIdx.x < NGLL2 ) {
    sh_tempx[tx] = wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_xz*xizl); // sh_tempx1
    sh_tempz[tx] = wxgll[J] *jacobianl * (sigma_xz*xixl + sigma_zz*xizl); // sh_tempz1
  }
  __syncthreads();
  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  if( threadIdx.x < NGLL2 ) {
    sh_tempx[tx] =wxgll[I] * jacobianl * (sigma_xx*gammaxl +  sigma_xz*gammazl); // sh_tempx3
    sh_tempz[tx] =wxgll[I] * jacobianl * (sigma_xz*gammaxl +  sigma_zz*gammazl); // sh_tempz3
  }
  __syncthreads();
  // 3. cut-plane gamma
  sum_hprimewgll_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprimewgll_xx);
  __syncthreads();


  sum_terms1= - tempx1l - tempx3l;
  sum_terms3= - tempz1l - tempz3l;


  // assembles acceleration array
  if(threadIdx.x < NGLL2) {

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*2]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*2) + sum_terms1;
    d_accel[iglob*2 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*2 + 1) + sum_terms3;
#else
    d_accel[iglob*2]     += sum_terms1;
    d_accel[iglob*2 + 1] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

#else // MESH_COLORING

    //mesh coloring
    if( use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*2]     = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*2) + sum_terms1;
      d_accel[iglob*2 + 1] = texfetch_accel<FORWARD_OR_ADJOINT>(iglob*2 + 1) + sum_terms3;
#else
      d_accel[iglob*2]     += sum_terms1;
      d_accel[iglob*2 + 1] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }else {
      atomicAdd(&d_accel[iglob*2], sum_terms1);
      atomicAdd(&d_accel[iglob*2+1], sum_terms3);
    } // if(use_mesh_coloring_gpu)

#endif // MESH_COLORING

  } // threadIdx.x

} // kernel_2_noatt_ani_impl()


/* ----------------------------------------------------------------------------------------------- */


void Kernel_2(int nb_blocks_to_compute,Mesh* mp,int d_iphase,realw d_deltat,
              int ANISOTROPY,
              int* d_ibool,
              realw* d_xix,realw* d_xiz,
              realw* d_gammax,realw* d_gammaz,
              realw* d_kappav,
              realw* d_muv,
              realw* d_c11store,realw* d_c12store,realw* d_c13store,
              realw* d_c15store,realw* d_c23store,realw* d_c25store,
              realw* d_c33store,realw* d_c35store,realw* d_c55store
              ){

  TRACE("\tKernel_2");

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before kernel Kernel 2");
#endif

  // if the grid can handle the number of blocks, we let it be 1D
  // grid_2_x = nb_elem_color;
  // nb_elem_color is just how many blocks we are computing now

  int blocksize = NGLL2_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(nb_blocks_to_compute,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  cudaEvent_t start,stop;
  if( CUDA_TIMING ){
    start_timing_cuda(&start,&stop);
  }



    // compute kernels without attenuation
    if( ANISOTROPY ){
      // full anisotropy
      // forward wavefields -> FORWARD_OR_ADJOINT == 1
      Kernel_2_noatt_ani_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                        d_ibool,
                                                                        mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                        d_iphase,
                                                                        mp->use_mesh_coloring_gpu,
                                                                        mp->d_displ,
                                                                        mp->d_accel,
                                                                        d_xix, d_xiz,
                                                                        d_gammax, d_gammaz,
                                                                        mp->d_hprime_xx,
                                                                        mp->d_hprimewgll_xx,
                                                                        mp->d_wxgll,
                                                                        d_kappav, d_muv,
                                                                        mp->simulation_type,
                                                                        ANISOTROPY,
                                                                        d_c11store,d_c12store,d_c13store,
                                                                        d_c15store,
                                                                        d_c23store,
                                                                        d_c25store,d_c33store,
                                                                        d_c35store,
                                                                        d_c55store);

      // backward/reconstructed wavefield
      if(mp->simulation_type == 3) {
        // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
        Kernel_2_noatt_ani_impl<3><<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                        d_ibool,
                                                                        mp->d_phase_ispec_inner_elastic,mp->num_phase_ispec_elastic,
                                                                        d_iphase,
                                                                        mp->use_mesh_coloring_gpu,
                                                                        mp->d_b_displ,
                                                                        mp->d_b_accel,
                                                                        d_xix, d_xiz,
                                                                        d_gammax, d_gammaz,
                                                                        mp->d_hprime_xx,
                                                                        mp->d_hprimewgll_xx,
                                                                        mp->d_wxgll,
                                                                        d_kappav, d_muv,
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
                                                                              d_kappav, d_muv,
                                                                              mp->simulation_type,
                                                                              mp->d_dsxx,
                                                                              mp->d_dsxz,
                                                                              mp->d_dszz);


            // backward/reconstructed wavefield
            if(mp->simulation_type == 3) {
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
                                                                                 d_kappav, d_muv,
                                                                                 mp->simulation_type,
                                                                                 mp->d_b_dsxx,
                                                                                 mp->d_b_dsxz,
                                                                                 mp->d_b_dszz);
            }

    } // ANISOTROPY


  // Cuda timing
  if( CUDA_TIMING ){

      if( ANISOTROPY ){
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
                                                int* ANISOTROPY) {

  TRACE("\tcompute_forces_viscoelastic_cuda");
  // EPIK_TRACER("compute_forces_viscoelastic_cuda");
  //printf("Running compute_forces\n");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
  int num_elements;

  if( *iphase == 1 )
    num_elements = *nspec_outer_elastic;
  else
    num_elements = *nspec_inner_elastic;

  // checks if anything to do
  if( num_elements == 0 ) return;

  // mesh coloring
  if( mp->use_mesh_coloring_gpu ){
    // note: array offsets require sorted arrays, such that e.g. ibool starts with elastic elements
    //         and followed by acoustic ones.
    //         elastic elements also start with outer than inner element ordering
    int nb_colors,nb_blocks_to_compute;
    int istart;
    int offset,offset_nonpadded;

    // sets up color loop
    if( *iphase == 1 ){
      // outer elements
      nb_colors = mp->num_colors_outer_elastic;
      istart = 0;

      // array offsets
      offset = 0;
      offset_nonpadded = 0;
    }else{
      // inner elements (start after outer elements)
      nb_colors = mp->num_colors_outer_elastic + mp->num_colors_inner_elastic;
      istart = mp->num_colors_outer_elastic;

      // array offsets
      offset = (*nspec_outer_elastic) * NGLL2_PADDED;
      offset_nonpadded = (*nspec_outer_elastic) * NGLL2;
    }

    // loops over colors
    for(int icolor = istart; icolor < nb_colors; icolor++){

      nb_blocks_to_compute = mp->h_num_elem_colors_elastic[icolor];

      // checks
      //if( nb_blocks_to_compute <= 0 ){
      //  printf("error number of elastic color blocks: %d -- color = %d \n",nb_blocks_to_compute,icolor);
      //  exit(EXIT_FAILURE);
      //}

      Kernel_2(nb_blocks_to_compute,mp,*iphase,*deltat,*ANISOTROPY,
               mp->d_ibool + offset,
               mp->d_xix + offset,mp->d_xiz + offset,
               mp->d_gammax + offset,mp->d_gammaz + offset,
               mp->d_kappav + offset,
               mp->d_muv + offset,
               mp->d_c11store + offset,mp->d_c12store + offset,mp->d_c13store + offset,
               mp->d_c15store + offset,mp->d_c23store + offset,mp->d_c25store + offset,
               mp->d_c33store + offset,mp->d_c35store + offset,mp->d_c55store + offset
               );

      // for padded and aligned arrays
      offset += nb_blocks_to_compute * NGLL2_PADDED;
      // for no-aligned arrays
      offset_nonpadded += nb_blocks_to_compute * NGLL2;


      //note: we use the same stream, so kernels are executed one after the other
      //      thus, there should be no need to synchronize in case we run on only 1 process to avoid race-conditions

    }

  }else{
    // no mesh coloring: uses atomic updates
    Kernel_2(num_elements,mp,*iphase,*deltat,*ANISOTROPY,
               mp->d_ibool,
               mp->d_xix,mp->d_xiz,
               mp->d_gammax,mp->d_gammaz,
               mp->d_kappav,
               mp->d_muv,
               mp->d_c11store,mp->d_c12store,mp->d_c13store,
               mp->d_c15store,mp->d_c23store,mp->d_c25store,
               mp->d_c33store,mp->d_c35store,mp->d_c55store
               );

  }
}
