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
realw_texture d_hprime_xx_tex;
//realw_texture d_hprimewgll_xx_tex;
realw_texture d_wxgll_xx_tex;
#endif


// note on performance optimizations:
//
//   instead of providing spezialized kernel routines (without mesh coloring, without gravity, etc.),
//   we only provide one "general" kernel to handle all cases. this reduces code redundancy and improves code readability.
//   as tradeoff, we take a little performance hit of around ~ 3%
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
//                note, especially the if( gravity )-block locations are very sensitive
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
                       realw* d_xix,realw* d_xiz,
                       realw* d_gammax,realw* d_gammaz,
                       realw_const_p d_hprime_xx,
                       realw_const_p hprimewgll_xx,
                       realw_const_p wxgll,
                       realw* d_rhostore,
                       const int use_mesh_coloring_gpu,
                       realw* d_kappastore){



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
  int working_element;

  realw temp1l,temp3l;
  realw xixl,xizl,gammaxl,gammazl;
  realw jacobianl;

  realw dpotentialdxl,dpotentialdzl;
  realw rho_invl;

  realw sum_terms;

  __shared__ realw s_dummy_loc[NGLL2];

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

// counts: for simulations without gravity, without mesh_coloring
//         counts floating-point operations (FLOP) per thread
//         counts global memory accesses in bytes (BYTES) per block
// 2 FLOP
//
// 0 BYTES

  // checks if anything to do
  if( bx >= nb_blocks_to_compute ) return;

  // limits thread ids to range [0,25-1]
  if( tx >= NGLL2 ) tx = tx - NGLL2;

// counts:
// + 1 FLOP
//
// + 0 BYTE

  // spectral-element id
#ifdef USE_MESH_COLORING_GPU
  working_element = bx;
#else
  //mesh coloring
  if( use_mesh_coloring_gpu ){
    working_element = bx;
  }else{
    // iphase-1 and working_element-1 for Fortran->C array conventions
    working_element = d_phase_ispec_inner_acoustic[bx + num_phase_ispec_acoustic*(d_iphase-1)]-1;
  }
#endif

  // local padded index
  offset = working_element*NGLL2_PADDED + tx;

  // global index
  iglob = d_ibool[offset] - 1;

// counts:
// + 7 FLOP
//
// + 2 float * 32 threads = 256 BYTE

  // loads potential values into shared memory
  if(threadIdx.x < NGLL2) {
#ifdef USE_TEXTURES_FIELDS
    s_dummy_loc[tx] = texfetch_potential<FORWARD_OR_ADJOINT>(iglob);
#else
    // changing iglob indexing to match fortran row changes fast style
    s_dummy_loc[tx] = d_potential_acoustic[iglob];
#endif
  }


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
  xixl = d_xix[offset];
  xizl = d_xiz[offset];
  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];

  jacobianl = 1.f / (xixl*gammazl-gammaxl*xizl);


  // density (reciproc)
  rho_invl = 1.f / d_rhostore[offset];

// counts:
// + 5 FLOP
//
// + 5 float * 32 threads = 160 BYTE

  // loads hprime into shared memory
  if (threadIdx.x < NGLL2) {
#ifdef USE_TEXTURES_CONSTANTS
    sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_tex,tx);
#else
    sh_hprime_xx[tx] = d_hprime_xx[tx];
#endif
    // loads hprimewgll into shared memory
    sh_hprimewgll_xx[tx] = hprimewgll_xx[tx];
  }
  else if (threadIdx.x < NGLL2 + NGLLX) {
#ifdef USE_TEXTURES_CONSTANTS
    sh_wxgll[tx] = tex1Dfetch(d_wxgll_xx_tex,tx);
#else
    // changing iglob indexing to match fortran row changes fast style
    sh_wxgll[tx] = wxgll[tx];
#endif
     }


// counts:
// + 0 FLOP
//
// + 2 * 1 float * 25 threads = 200 BYTE

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
    temp1l += s_dummy_loc[J*NGLLX+l] * sh_hprime_xx[l*NGLLX+I];
    // 3. cut-plane along gamma-direction
    temp3l += s_dummy_loc[l*NGLLX+I] * sh_hprime_xx[l*NGLLX+J];
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
  if( threadIdx.x < NGLL2 ) {
    s_temp1[tx] = sh_wxgll[J] * jacobianl * rho_invl * (dpotentialdxl*xixl  + dpotentialdzl*xizl)  ;
    s_temp3[tx] = sh_wxgll[I] * jacobianl * rho_invl * (dpotentialdxl*gammaxl + dpotentialdzl*gammazl)  ;
  }

// counts:
// + 2 * 6 FLOP = 12 FLOP
//
// + 2 BYTE


  // synchronize all the threads (one thread for each of the NGLL grid points of the
  // current spectral element) because we need the whole element to be ready in order
  // to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  sum_terms=0;

  for (int l=0;l<NGLLX;l++) {
    //assumes hprimewgll_xx = hprimewgll_zz

  sum_terms -= s_temp1[J*NGLLX+l] * sh_hprimewgll_xx[I*NGLLX+l] + s_temp3[l*NGLLX+I] * sh_hprimewgll_xx[J*NGLLX+l];

  }


// counts:
// + NGLLX * 11 FLOP = 55 FLOP
//
// + 0 BYTE



  // assembles potential array
  if(threadIdx.x < NGLL2) {
#ifdef USE_MESH_COLORING_GPU
  // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
    d_potential_dot_dot_acoustic[iglob] = texfetch_potential_dot_dot<FORWARD_OR_ADJOINT>(iglob) + sum_terms;
#else
    d_potential_dot_dot_acoustic[iglob] += sum_terms;
#endif // USE_TEXTURES_FIELDS
#else  // MESH_COLORING
    //mesh coloring
    if( use_mesh_coloring_gpu ){
      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_potential_dot_dot_acoustic[iglob] = texfetch_potential_dot_dot<FORWARD_OR_ADJOINT>(iglob) + sum_terms;
#else
      d_potential_dot_dot_acoustic[iglob] += sum_terms;
#endif // USE_TEXTURES_FIELDS
    }else{
      atomicAdd(&d_potential_dot_dot_acoustic[iglob],sum_terms);
    }
#endif // MESH_COLORING
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


/* ----------------------------------------------------------------------------------------------- */

void Kernel_2_acoustic(int nb_blocks_to_compute, Mesh* mp, int d_iphase,
                       int* d_ibool,
                       realw* d_xix,realw* d_xiz,
                       realw* d_gammax,realw* d_gammaz,
                       realw* d_rhostore,
                       realw* d_kappastore){

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before acoustic kernel Kernel 2");
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
  cudaEvent_t start, stop;
  if( CUDA_TIMING ){
    start_timing_cuda(&start,&stop);
  }

  // forward wavefields -> FORWARD_OR_ADJOINT == 1
  Kernel_2_acoustic_impl<1><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
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
                                                                    mp->use_mesh_coloring_gpu,
                                                                    d_kappastore);

  if(mp->simulation_type == 3) {
    // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
    Kernel_2_acoustic_impl<3><<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
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
                                                                      mp->use_mesh_coloring_gpu,
                                                                      d_kappastore);
  }

  // Cuda timing
  if( CUDA_TIMING ){
    realw flops,time;
    stop_timing_cuda(&start,&stop,"Kernel_2_acoustic_impl",&time);
    // time in seconds
    time = time / 1000.;

      if( ! mp->use_mesh_coloring_gpu ){
        // see with: nvprof --metrics flops_sp ./xspecfem3D
        //           -> using 322631424 FLOPS (Single) floating-point operations for 20736 elements
        //              = 15559 FLOPS per block
        flops = 15559 * nb_blocks_to_compute;
      }else{
        // coloring
        flops = 15559 * nb_blocks_to_compute;
      }

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
                                            int* nspec_inner_acoustic) {

  TRACE("compute_forces_acoustic_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper


if(mp->simulation_type==3){
  #ifdef USE_TEXTURES_FIELDS
  {
      int size = mp->NGLOB_AB;
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      const textureReference* d_b_potential_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_b_potential_tex_ref_ptr, "d_b_potential_tex"), 3001);
      print_CUDA_error_if_any(cudaBindTexture(0, d_b_potential_tex_ref_ptr, mp->d_b_potential_acoustic, &channelDesc, sizeof(realw)*size), 3001);

      const textureReference* d_b_potential_dot_dot_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_b_potential_dot_dot_tex_ref_ptr, "d_b_potential_dot_dot_tex"),3003);
      print_CUDA_error_if_any(cudaBindTexture(0, d_b_potential_dot_dot_tex_ref_ptr, mp->d_b_potential_dot_dot_acoustic, &channelDesc, sizeof(realw)*size), 3003);
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_b_potential_tex, mp->d_b_potential_acoustic, &channelDesc, sizeof(realw)*size), 3001);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_b_potential_dot_dot_tex, mp->d_b_potential_dot_dot_acoustic, &channelDesc, sizeof(realw)*size), 3003);
    #endif
  }
  #endif

}

  #ifdef USE_TEXTURES_FIELDS
  {
      int size = mp->NGLOB_AB;
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      const textureReference* d_potential_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_potential_tex_ref_ptr, "d_potential_tex"), 2001);
      print_CUDA_error_if_any(cudaBindTexture(0, d_potential_tex_ref_ptr, mp->d_potential_acoustic, &channelDesc, sizeof(realw)*size), 2001);

      const textureReference* d_potential_dot_dot_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_potential_dot_dot_tex_ref_ptr, "d_potential_dot_dot_tex"), 2003);
      print_CUDA_error_if_any(cudaBindTexture(0, d_potential_dot_dot_tex_ref_ptr, mp->d_potential_dot_dot_acoustic, &channelDesc, sizeof(realw)*size), 2003);
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_potential_tex, mp->d_potential_acoustic, &channelDesc, sizeof(realw)*size), 2001);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_potential_dot_dot_tex, mp->d_potential_dot_dot_acoustic, &channelDesc, sizeof(realw)*size), 2003);
    #endif
  }
  #endif



  #ifdef USE_TEXTURES_CONSTANTS
  {
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      const textureReference* d_hprime_xx_tex_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_hprime_xx_tex_ptr, "d_hprime_xx_tex"), 4101);
      print_CUDA_error_if_any(cudaBindTexture(0, d_hprime_xx_tex_ptr, mp->d_hprime_xx, &channelDesc, sizeof(realw)*(NGLL2)), 4001);
      const textureReference* d_wxgll_xx_tex_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_wxgll_xx_tex_ptr, "d_wxgll_xx_tex"), 4101);
      print_CUDA_error_if_any(cudaBindTexture(0, d_wxgll_xx_tex_ptr, mp->d_wxgll, &channelDesc, sizeof(realw)*(NGLL2)), 4001);
   #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
//      print_CUDA_error_if_any(cudaBindTexture(0, &d_hprimewgll_xx_tex, mp->d_hprimewgll_xx, &channelDesc, sizeof(realw)*(NGLL2)), 40010);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_hprime_xx_tex, mp->d_hprime_xx, &channelDesc, sizeof(realw)*(NGLL2)), 4002);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_wxgll_xx_tex, mp->d_wxgll, &channelDesc, sizeof(realw)*(NGLLX)), 40013);
   #endif
  }
  #endif


  int num_elements;

  if( *iphase == 1 )
    num_elements = *nspec_outer_acoustic;
  else
    num_elements = *nspec_inner_acoustic;

  if( num_elements == 0 ) return;

  // mesh coloring
  if( mp->use_mesh_coloring_gpu ){

    // note: array offsets require sorted arrays, such that e.g. ibool starts with elastic elements
    //         and followed by acoustic ones.
    //         acoustic elements also start with outer than inner element ordering

    int nb_colors,nb_blocks_to_compute;
    int istart;
    int offset,offset_nonpadded;

    // sets up color loop
    if( *iphase == 1 ){
      // outer elements
      nb_colors = mp->num_colors_outer_acoustic;
      istart = 0;

      // array offsets (acoustic elements start after elastic ones)
      offset = mp->nspec_elastic * NGLL2_PADDED;
      offset_nonpadded = mp->nspec_elastic * NGLL2;
    }else{
      // inner element colors (start after outer elements)
      nb_colors = mp->num_colors_outer_acoustic + mp->num_colors_inner_acoustic;
      istart = mp->num_colors_outer_acoustic;

      // array offsets (inner elements start after outer ones)
      offset = ( mp->nspec_elastic + (*nspec_outer_acoustic) ) * NGLL2_PADDED;
      offset_nonpadded = ( mp->nspec_elastic + (*nspec_outer_acoustic) ) * NGLL2;
    }

    // loops over colors
    for(int icolor = istart; icolor < nb_colors; icolor++){

      nb_blocks_to_compute = mp->h_num_elem_colors_acoustic[icolor];

      Kernel_2_acoustic(nb_blocks_to_compute,mp,*iphase,
                         mp->d_ibool + offset_nonpadded,
                         mp->d_xix + offset,mp->d_xiz + offset,
                         mp->d_gammax + offset,mp->d_gammaz + offset,
                         mp->d_rhostore + offset,
                         mp->d_kappastore + offset_nonpadded);

      // for padded and aligned arrays
      offset += nb_blocks_to_compute * NGLL2_PADDED;
      // for no-aligned arrays
      offset_nonpadded += nb_blocks_to_compute * NGLL2;
    }

  }else{

    // no mesh coloring: uses atomic updates
    Kernel_2_acoustic(num_elements, mp, *iphase,
                      mp->d_ibool,
                      mp->d_xix,mp->d_xiz,
                      mp->d_gammax,mp->d_gammaz,
                      mp->d_rhostore,
                      mp->d_kappastore);

  }
}



/* ----------------------------------------------------------------------------------------------- */

/* KERNEL for enforce free surface */

/* ----------------------------------------------------------------------------------------------- */


__global__ void enforce_free_surface_cuda_kernel(
                                       realw_p potential_acoustic,
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
  if( iface < num_free_surface_faces ){

    int ispec = free_surface_ispec[iface]-1;

    // checks if element is in acoustic domain
    if( ispec_is_acoustic[ispec] ){

      // gets global point index
      int igll = threadIdx.x + threadIdx.y*blockDim.x;

      int i = free_surface_ij[INDEX3(NDIM,NGLLX,0,igll,iface)] - 1; // (1,igll,iface)
      int j = free_surface_ij[INDEX3(NDIM,NGLLX,1,igll,iface)] - 1;

      int iglob = d_ibool[INDEX3_PADDED(NGLLX,NGLLX,i,j,ispec)] - 1;

      // sets potentials to zero at free surface
      potential_acoustic[iglob] = 0;
      potential_dot_acoustic[iglob] = 0;
      potential_dot_dot_acoustic[iglob] = 0;
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(acoustic_enforce_free_surf_cuda,
              ACOUSTIC_ENFORCE_FREE_SURF_CUDA)(long* Mesh_pointer) {

TRACE("acoustic_enforce_free_surf_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container


    // does not absorb free surface, thus we enforce the potential to be zero at surface

    // block sizes
    int num_blocks_x, num_blocks_y;
    get_blocks_xy(mp->num_free_surface_faces,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y,1);
    dim3 threads(NGLLX,1,1);

    // sets potentials to zero at free surface
    enforce_free_surface_cuda_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_acoustic,
                                                                             mp->d_potential_dot_acoustic,
                                                                             mp->d_potential_dot_dot_acoustic,
                                                                             mp->num_free_surface_faces,
                                                                             mp->d_free_surface_ispec,
                                                                             mp->d_free_surface_ijk,
                                                                             mp->d_ibool,
                                                                             mp->d_ispec_is_acoustic);
    // for backward/reconstructed potentials
    if(mp->simulation_type == 3) {
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

