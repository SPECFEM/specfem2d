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
                        const int simulation_type,
                        const int p_sv){

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

  // checks if anything to do
  if (bx >= nb_blocks_to_compute ) return;

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
  if (p_sv){
    // P_SV case
    sigma_xx = lambdalplus2mul*duxdxl + lambdal*duzdzl;
    sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl;
    sigma_xz = mul*duzdxl_plus_duxdzl;
  }else{
    // SH-case
    sigma_xx = mul * duxdxl;  // would be sigma_xy in CPU-version
    sigma_xz = mul * duxdzl;  // sigma_zy
  }

// counts:
// + 22 FLOP
//
// + 0 BYTE

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  if (threadIdx.x < NGLL2) {
    if (p_sv){
      // P_SV case
      sh_tempx[tx] = sh_wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_xz*xizl); // sh_tempx1
      sh_tempz[tx] = sh_wxgll[J] *jacobianl * (sigma_xz*xixl + sigma_zz*xizl); // sh_tempz1
    }else{
      // SH-case
      sh_tempx[tx] = sh_wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_xz*xizl); // sh_tempx1
      sh_tempz[tx] = 0.f;
    }
  }
  __syncthreads();

  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);
  __syncthreads();

  if (threadIdx.x < NGLL2) {
    if (p_sv){
      // P_SV case
      sh_tempx[tx] = sh_wxgll[I] * jacobianl * (sigma_xx*gammaxl +  sigma_xz*gammazl); // sh_tempx3
      sh_tempz[tx] = sh_wxgll[I] * jacobianl * (sigma_xz*gammaxl +  sigma_zz*gammazl); // sh_tempz3
    }else{
      // SH-case
      sh_tempx[tx] = sh_wxgll[I] * jacobianl * (sigma_xx*gammaxl +  sigma_xz*gammazl); // sh_tempx3
      sh_tempz[tx] = 0.f; // sh_tempz3
    }
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
// + 2 float * 25 threads = 50 BYTE


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
                        const int simulation_type,
                        const int p_sv,
                        const int* ispec_is_anisotropic,
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
  if (threadIdx.x < NGLL2) {
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

  // stress calculations
  if (ispec_is_anisotropic[working_element]){
    // full anisotropic case
    c11 = d_c11store[offset];
    c13 = d_c13store[offset];
    c15 = d_c15store[offset];
    c33 = d_c33store[offset];
    c35 = d_c35store[offset];
    c55 = d_c55store[offset];

    // compute the three components of the stress tensor sigma (full anisotropy)
    if (p_sv){
      // P_SV case
      sigma_xx = c11*duxdxl + c13*duzdzl + c15*duzdxl_plus_duxdzl;
      sigma_zz = c13*duxdxl + c33*duzdzl + c35*duzdxl_plus_duxdzl;
      sigma_xz = c15*duxdxl + c35*duzdzl + c55*duzdxl_plus_duxdzl;
      sigma_zx = sigma_xz;
    }else{
      // SH-case
      sigma_xx = c55 * duxdxl;  // assumes c55 == mu, and still isotropic in both directions - no anisotropy implemented yet...
      sigma_xz = c55 * duxdzl;
    }
  }else{
    // isotropic case

    // compute elements with an elastic isotropic rheology
    kappal = d_kappav[offset];
    mul = d_muv[offset];

    // original
    //lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
    //lambdal = lambdalplus2mul - 2.0f * mul;

    // new
    lambdal = kappal - mul;
    lambdalplus2mul = kappal + mul;

    // compute the three components of the stress tensor sigma
    if (p_sv){
      // P_SV case
      sigma_xx = lambdalplus2mul*duxdxl + lambdal*duzdzl;
      sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl;
      sigma_xz = mul*duzdxl_plus_duxdzl;
      sigma_zx = sigma_xz;
    }else{
      // SH-case
      sigma_xx = mul * duxdxl;  // would be sigma_xy in CPU-version
      sigma_xz = mul * duxdzl;  // sigma_zy
    }
  }

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  if (threadIdx.x < NGLL2) {
    if (p_sv){
      // P_SV case
      sh_tempx[tx] = wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_zx*xizl); // sh_tempx1
      sh_tempz[tx] = wxgll[J] *jacobianl * (sigma_xz*xixl + sigma_zz*xizl); // sh_tempz1
    }else{
      // SH-case
      sh_tempx[tx] = wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_xz*xizl); // sh_tempx1
      sh_tempz[tx] = 0.f;
    }
  }
  __syncthreads();

  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  if (threadIdx.x < NGLL2) {
    if (p_sv){
      // P_SV case
      sh_tempx[tx] = wxgll[I] * jacobianl * (sigma_xx*gammaxl + sigma_zx*gammazl); // sh_tempx3
      sh_tempz[tx] = wxgll[I] * jacobianl * (sigma_xz*gammaxl + sigma_zz*gammazl); // sh_tempz3
    }else{
      // SH-case
      sh_tempx[tx] = wxgll[I] * jacobianl * (sigma_xx*gammaxl + sigma_xz*gammazl); // sh_tempx3
      sh_tempz[tx] = 0.f; // sh_tempz3
    }
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
                      const int simulation_type,
                      const int p_sv,
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

  // attenuation
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

  // attenuation
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
  if (p_sv){
    // P_SV case
    sigma_xx = lambdalplus2mul*duxdxl + lambdal*duzdzl;
    sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl;
    sigma_xz = mul*duzdxl_plus_duxdzl;
  }else{
    // SH-case
    sigma_xx = mul * duxdxl;  // would be sigma_xy in CPU-version
    sigma_xz = mul * duxdzl;  // sigma_zy
  }

  // attenuation
  // get the contribution of attenuation and update the memory variables
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duxdxl_old = dux_dxl_old[offset_align];
  duzdzl_old = duz_dzl_old[offset_align];
  duxdxl_plus_duzdzl_old = duxdxl_old + duzdzl_old;
  duxdzl_plus_duzdxl_old = dux_dzl_plus_duz_dxl_old[offset_align];

  e1_sum = 0.f;
  e11_sum = 0.f;
  e13_sum = 0.f;
  for (i_sls=0;i_sls<N_SLS;i_sls++){
    // bulk attenuation
    a_newmark = A_newmark_kappa[N_SLS * offset_align + i_sls];
    b_newmark = B_newmark_kappa[N_SLS * offset_align + i_sls];

    e1_load[i_sls] = a_newmark * a_newmark * e1_load[i_sls] + b_newmark * (duxdxl_plus_duzdzl + a_newmark * (duxdxl_plus_duzdzl_old));
    e1_sum += e1_load[i_sls];
    e1[N_SLS*offset_align+i_sls] = e1_load[i_sls];

    // shear attenuation
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
  if (p_sv){
    // P_SV case
    sigma_xx += (lambdalplus2mul-mul) * e1_sum + 2.0f * mul * e11_sum;
    sigma_xz += mul * e13_sum;
    sigma_zz += (lambdalplus2mul-mul) * e1_sum - 2.0f * mul * e11_sum;
  }else{
    // SH-case
    sigma_xx += 0.f;  // attenuation not implemented yet for SH
    sigma_xz += 0.f;
  }

  // saves the grad(displ) to use at the next iteration
  dux_dxl_old[offset_align] = duxdxl;
  duz_dzl_old[offset_align] = duzdzl;
  dux_dzl_plus_duz_dxl_old[offset_align] = duzdxl_plus_duxdzl;

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  if (threadIdx.x < NGLL2) {
    if (p_sv){
      // P_SV case
      sh_tempx[tx] = sh_wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_xz*xizl); // sh_tempx1
      sh_tempz[tx] = sh_wxgll[J] *jacobianl * (sigma_xz*xixl + sigma_zz*xizl); // sh_tempz1
    }else{
      // SH-case
      sh_tempx[tx] = sh_wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_xz*xizl); // sh_tempx1
      sh_tempz[tx] = 0.f;
    }
  }
  __syncthreads();

  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);
  __syncthreads();

  if (threadIdx.x < NGLL2) {
    if (p_sv){
      // P_SV case
      sh_tempx[tx] = sh_wxgll[I] * jacobianl * (sigma_xx*gammaxl +  sigma_xz*gammazl); // sh_tempx3
      sh_tempz[tx] = sh_wxgll[I] * jacobianl * (sigma_xz*gammaxl +  sigma_zz*gammazl); // sh_tempz3
    }else{
      // SH-case
      sh_tempx[tx] = sh_wxgll[I] * jacobianl * (sigma_xx*gammaxl +  sigma_xz*gammazl); // sh_tempx3
      sh_tempz[tx] = 0.f; // sh_tempz3
    }
  }
  __syncthreads();

  // 3. cut-plane gamma
  sum_hprimewgll_gamma(I,J,&tempx3l,&tempz3l,sh_tempx,sh_tempz,sh_hprimewgll_xx);
  __syncthreads();

  sum_terms1 = -tempx1l - tempx3l;
  sum_terms3 = -tempz1l - tempz3l;

  // assembles acceleration array
  if (threadIdx.x < NGLL2) {
    atomicAdd(&d_accel[iglob*2], sum_terms1);
    atomicAdd(&d_accel[iglob*2+1], sum_terms3);
  }
} // kernel_2_att_iso_impl()

/* ----------------------------------------------------------------------------------------------- */


template<int FORWARD_OR_ADJOINT> __global__ void
#ifdef USE_LAUNCH_BOUNDS
// adds compiler specification
__launch_bounds__(NGLL2_PADDED,LAUNCH_MIN_BLOCKS)
#endif
// main kernel
Kernel_2_att_ani_impl(int nb_blocks_to_compute,
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
                      const int simulation_type,
                      const int p_sv,
                      const int* ispec_is_anisotropic,
                      realw* d_c11store,realw* d_c12store,realw* d_c13store,
                      realw* d_c15store,
                      realw* d_c23store,
                      realw* d_c25store,realw* d_c33store,
                      realw* d_c35store,
                      realw* d_c55store,
                      realw_const_p A_newmark_mu,realw_const_p B_newmark_mu,
                      realw_const_p A_newmark_kappa,realw_const_p B_newmark_kappa,
                      realw_p e1,realw_p e11,realw_p e13,
                      realw_p dux_dxl_old,realw_p duz_dzl_old,realw_p dux_dzl_plus_duz_dxl_old) {

// elastic compute kernel without attenuation for anisotropic elements
//
// holds for:
//  ATTENUATION               = .true.
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
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw sigma_xx,sigma_zz,sigma_xz,sigma_zx;
  realw c11,c13,c15,c33,c35,c55;
  realw sum_terms1,sum_terms3;

  // attenuation
  int offset_align;
  realw duzdxl_plus_duxdzl,duxdxl_plus_duzdzl;
  realw duxdxl_old,duzdzl_old;
  realw duxdzl_plus_duzdxl_old,duxdxl_plus_duzdzl_old;
  realw e1_load[N_SLS],e11_load[N_SLS],e13_load[N_SLS];
  realw e1_sum,e11_sum,e13_sum,a_newmark,b_newmark;

  // shared memory
  __shared__ realw sh_tempx[NGLL2];
  __shared__ realw sh_tempz[NGLL2];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // loads hprime's into shared memory
  if (threadIdx.x < NGLL2) {
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

  // attenuation
  offset_align = working_element*NGLL2 + tx;
  for (int i_sls=0;i_sls<N_SLS;i_sls++){
    e1_load[i_sls] = e1[N_SLS*offset_align+i_sls];
    e11_load[i_sls] = e11[N_SLS*offset_align+i_sls];
    e13_load[i_sls] = e13[N_SLS*offset_align+i_sls];
  }

  // loads mesh values here to give compiler possibility to overlap memory fetches with some computations
  // note: arguments defined as realw* instead of const realw* __restrict__ to avoid that the compiler
  //       loads all memory by texture loads
  //       we only use the first loads explicitly by texture loads, all subsequent without. this should lead/trick
  //       the compiler to use global memory loads for all the subsequent accesses.
  //
  // calculates laplacian
  xixl = get_global_cr( &d_xix[offset] ); // first array with texture load
  xizl = d_xiz[offset]; // all subsequent without to avoid over-use of texture for coalescent access

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

  // compute elements with an elastic isotropic rheology
  // note: also needed for anisotropy with attenuation case
  kappal = d_kappav[offset];
  mul = d_muv[offset];

  // original
  //lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  //lambdal = lambdalplus2mul - 2.0f * mul;
  // new
  lambdal = kappal - mul;
  lambdalplus2mul = kappal + mul;

  // stress calculations
  if (ispec_is_anisotropic[working_element]){
    // full anisotropic case
    c11 = d_c11store[offset];
    c13 = d_c13store[offset];
    c15 = d_c15store[offset];
    c33 = d_c33store[offset];
    c35 = d_c35store[offset];
    c55 = d_c55store[offset];

    // compute the three components of the stress tensor sigma (full anisotropy)
    if (p_sv){
      // P_SV case
      sigma_xx = c11*duxdxl + c13*duzdzl + c15*duzdxl_plus_duxdzl;
      sigma_zz = c13*duxdxl + c33*duzdzl + c35*duzdxl_plus_duxdzl;
      sigma_xz = c15*duxdxl + c35*duzdzl + c55*duzdxl_plus_duxdzl;
      sigma_zx = sigma_xz;
    }else{
      // SH-case
      sigma_xx = c55 * duxdxl;  // assumes c55 == mu, and still isotropic in both directions - no anisotropy implemented yet...
      sigma_xz = c55 * duxdzl;
    }
  }else{
    // isotropic case
    // compute the three components of the stress tensor sigma
    if (p_sv){
      // P_SV case
      sigma_xx = lambdalplus2mul*duxdxl + lambdal*duzdzl;
      sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl;
      sigma_xz = mul*duzdxl_plus_duxdzl;
      sigma_zx = sigma_xz;
    }else{
      // SH-case
      sigma_xx = mul * duxdxl;  // would be sigma_xy in CPU-version
      sigma_xz = mul * duxdzl;  // sigma_zy
    }
  }

  // attenuation
  // get the contribution of attenuation and update the memory variables
  duxdxl_plus_duzdzl = duxdxl + duzdzl;
  duxdxl_old = dux_dxl_old[offset_align];
  duzdzl_old = duz_dzl_old[offset_align];
  duxdxl_plus_duzdzl_old = duxdxl_old + duzdzl_old;
  duxdzl_plus_duzdxl_old = dux_dzl_plus_duz_dxl_old[offset_align];

  e1_sum = 0.f;
  e11_sum = 0.f;
  e13_sum = 0.f;
  for (int i_sls=0;i_sls<N_SLS;i_sls++){
    // bulk attenuation
    a_newmark = A_newmark_kappa[N_SLS * offset_align + i_sls];
    b_newmark = B_newmark_kappa[N_SLS * offset_align + i_sls];

    e1_load[i_sls] = a_newmark * a_newmark * e1_load[i_sls] + b_newmark * (duxdxl_plus_duzdzl + a_newmark * (duxdxl_plus_duzdzl_old));
    e1_sum += e1_load[i_sls];
    e1[N_SLS*offset_align+i_sls] = e1_load[i_sls];

    // shear attenuation
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
  if (p_sv){
    // P_SV case
    sigma_xx += (lambdalplus2mul-mul) * e1_sum + 2.0f * mul * e11_sum;
    sigma_zz += (lambdalplus2mul-mul) * e1_sum - 2.0f * mul * e11_sum;
    sigma_xz += mul * e13_sum;
    sigma_zx = sigma_xz;
  }else{
    // SH-case
    sigma_xx += 0.f;  // attenuation not implemented yet for SH
    sigma_xz += 0.f;
  }

  // saves the grad(displ) to use at the next iteration
  dux_dxl_old[offset_align] = duxdxl;
  duz_dzl_old[offset_align] = duzdzl;
  dux_dzl_plus_duz_dxl_old[offset_align] = duzdxl_plus_duxdzl;

  // form dot product with test vector, non-symmetric form
  // 1. cut-plane xi
  __syncthreads();
  if (threadIdx.x < NGLL2) {
    if (p_sv){
      // P_SV case
      sh_tempx[tx] = wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_zx*xizl); // sh_tempx1
      sh_tempz[tx] = wxgll[J] *jacobianl * (sigma_xz*xixl + sigma_zz*xizl); // sh_tempz1
    }else{
      // SH-case
      sh_tempx[tx] = wxgll[J] *jacobianl * (sigma_xx*xixl + sigma_xz*xizl); // sh_tempx1
      sh_tempz[tx] = 0.f;
    }
  }
  __syncthreads();

  // 1. cut-plane xi
  sum_hprimewgll_xi(I,J,&tempx1l,&tempz1l,sh_tempx,sh_tempz,sh_hprimewgll_xx);

  // 3. cut-plane gamma
  __syncthreads();
  if (threadIdx.x < NGLL2) {
    if (p_sv){
      // P_SV case
      sh_tempx[tx] = wxgll[I] * jacobianl * (sigma_xx*gammaxl + sigma_zx*gammazl); // sh_tempx3
      sh_tempz[tx] = wxgll[I] * jacobianl * (sigma_xz*gammaxl + sigma_zz*gammazl); // sh_tempz3
    }else{
      // SH-case
      sh_tempx[tx] = wxgll[I] * jacobianl * (sigma_xx*gammaxl + sigma_xz*gammazl); // sh_tempx3
      sh_tempz[tx] = 0.f; // sh_tempz3
    }
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
} // kernel_2_att_ani_impl()


/* ----------------------------------------------------------------------------------------------- */

// note: we used templating to be able to call the same kernel_2 twice for both,
//       forward and backward wavefields. that is, calling it by
//          Kernel_2_noatt_iso_impl<1>
//       and
//          Kernel_2_noatt_iso_impl<3>
//       the templating helped to use textures for forward/backward fields.
//
//       most of this has become obsolete, textures are hardly needed for speedup anymore
//       and the Kernel_2 has become more and more specialized for different cases to
//       reduce register pressure and increase occupancy for better performance.
//       thus, in future we might re-evaluate and remove this template-feature.
//
// "forced" template instantiation
// see: https://isocpp.org/wiki/faq/templates#separate-template-fn-defn-from-decl
//      https://stackoverflow.com/questions/31705764/cuda-c-using-a-template-function-which-calls-a-template-kernel
//
// for compute_forces_viscoelastic_cuda.cu:
// Kernel_2_noatt_iso_impl<1> needs an explicit instantiation here to be able to link against it from a different .cu file, ..

template __global__ void Kernel_2_noatt_iso_impl<1>(const int,const int*,const int*,const int,const int,
                                                    realw_const_p,realw_p,
                                                    realw*,realw*,realw*,realw*,
                                                    realw_const_p,realw_const_p,realw_const_p,
                                                    realw*,realw*,const int,const int);

template __global__ void Kernel_2_noatt_iso_impl<3>(const int,const int*,const int*,const int,const int,
                                                    realw_const_p,realw_p,
                                                    realw*,realw*,realw*,realw*,
                                                    realw_const_p,realw_const_p,realw_const_p,
                                                    realw*,realw*,const int,const int);

template __global__ void Kernel_2_noatt_ani_impl<1>(int,const int*,const int*,const int,const int,
                                                    realw_const_p,realw_p,
                                                    realw*,realw*,realw*,realw*,
                                                    realw_const_p,realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                    const int,const int,const int*,
                                                    realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*);

template __global__ void Kernel_2_noatt_ani_impl<3>(int,const int*,const int*,const int,const int,
                                                    realw_const_p,realw_p,
                                                    realw*,realw*,realw*,realw*,
                                                    realw_const_p,realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                    const int,const int,const int*,
                                                    realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*,realw*);

template __global__ void Kernel_2_att_iso_impl<1>(const int,const int*,const int*,const int,const int,
                                                  realw_const_p,realw_p,
                                                  realw*,realw*,realw*,realw*,
                                                  realw_const_p,realw_const_p,realw_const_p,
                                                  realw*,realw*,const int,const int,
                                                  realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                  realw_p,realw_p,realw_p,realw_p,realw_p,realw_p);

template __global__ void Kernel_2_att_iso_impl<3>(const int,const int*,const int*,const int,const int,
                                                  realw_const_p,realw_p,
                                                  realw*,realw*,realw*,realw*,
                                                  realw_const_p,realw_const_p,realw_const_p,
                                                  realw*,realw*,const int,const int,
                                                  realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                  realw_p,realw_p,realw_p,realw_p,realw_p,realw_p);

template __global__ void Kernel_2_att_ani_impl<1>(int,const int*,const int*,const int,const int,
                                                  realw_const_p,realw_p,
                                                  realw*,realw*,realw*,realw*,
                                                  realw_const_p,realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                  const int,const int,const int*,
                                                  realw*,realw*,realw*,realw*,realw*,realw*,realw*,
                                                  realw*,realw*,
                                                  realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                  realw_p,realw_p,realw_p,realw_p,realw_p,realw_p);

template __global__ void Kernel_2_att_ani_impl<3>(int,const int*,const int*,const int,const int,
                                                  realw_const_p,realw_p,
                                                  realw*,realw*,realw*,realw*,
                                                  realw_const_p,realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                  const int,const int,const int*,
                                                  realw*,realw*,realw*,realw*,realw*,realw*,realw*,
                                                  realw*,realw*,
                                                  realw_const_p,realw_const_p,realw_const_p,realw_const_p,
                                                  realw_p,realw_p,realw_p,realw_p,realw_p,realw_p);
