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

#ifndef CUDA_HEADER_H
#define CUDA_HEADER_H





typedef float realw;  // type of "working" variables

// CUDA version >= 5.0 needed for new symbol addressing and texture binding
#if CUDA_VERSION < 5000
  #ifndef USE_OLDER_CUDA4_GPU
    #define USE_OLDER_CUDA4_GPU
  #endif
#else
  #undef USE_OLDER_CUDA4_GPU
#endif

#ifdef USE_OLDER_CUDA4_GPU
#pragma message ("\nCompiling with: USE_OLDER_CUDA4_GPU enabled\n")
#endif






/* ----------------------------------------------------------------------------------------------- */

// CONSTANT arrays setup

/* ----------------------------------------------------------------------------------------------- */

/* note:
 constant arrays when used in compute_forces_acoustic_cuda.cu routines stay zero,
 constant declaration and cudaMemcpyToSymbol would have to be in the same file...

 extern keyword doesn't work for __constant__ declarations.

 also:
 cudaMemcpyToSymbol("deviceCaseParams", caseParams, sizeof(CaseParams));
 ..
 and compile with -arch=sm_20

 see also: http://stackoverflow.com/questions/4008031/how-to-use-cuda-constant-memory-in-a-programmer-pleasant-way
 doesn't seem to work.

 we could keep arrays separated for acoustic and elastic routines...

 for now, we store pointers with cudaGetSymbolAddress() function calls.

 */

// cuda constant arrays
//
// note: we use definition __device__ to use global memory rather than constant memory registers
//          to avoid over-loading registers; this should help increasing the occupancy on the GPU

__device__ realw d_hprime_xx[NGLL2];
//__device__ realw d_hprime_yy[NGLL2]; // only needed if NGLLX != NGLLY != NGLLZ
//__device__ realw d_hprime_zz[NGLL2]; // only needed if NGLLX != NGLLY != NGLLZ

__device__ realw d_hprimewgll_xx[NGLL2];

__device__ realw d_wxgll[NGLLX];

//__device__ realw d_hprimewgll_yy[NGLL2]; // only needed if NGLLX != NGLLY != NGLLZ
//__device__ realw d_hprimewgll_zz[NGLL2]; // only needed if NGLLX != NGLLY != NGLLZ

void setConst_wxgll(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_wxgll, array, NGLLX*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_wxgll: %s\n", cudaGetErrorString(err));
    exit(1);
  }
#ifdef USE_OLDER_CUDA4_GPU
  err = cudaGetSymbolAddress((void**)&(mp->d_wxgll),"d_wxgll");
#else
  err = cudaGetSymbolAddress((void**)&(mp->d_wxgll),d_wxgll);
#endif
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_wxgll: %s\n", cudaGetErrorString(err));
    exit(1);
  }

}

void setConst_hprime_xx(realw* array,Mesh* mp)
{

  cudaError_t err = cudaMemcpyToSymbol(d_hprime_xx, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprime_xx: %s\n", cudaGetErrorString(err));
    fprintf(stderr, "The problem is maybe -arch sm_13 instead of -arch sm_11 in the Makefile, please doublecheck\n");
    exit(1);
  }

#ifdef USE_OLDER_CUDA4_GPU
  err = cudaGetSymbolAddress((void**)&(mp->d_hprime_xx),"d_hprime_xx");
#else
  err = cudaGetSymbolAddress((void**)&(mp->d_hprime_xx),d_hprime_xx);
#endif
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprime_xx: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}


// void setConst_hprime_zz(realw* array,Mesh* mp)
// {

//   cudaError_t err = cudaMemcpyToSymbol(d_hprime_zz, array, NGLL2*sizeof(realw));
//   if (err != cudaSuccess)
//   {
//     fprintf(stderr, "Error in setConst_hprime_zz: %s\n", cudaGetErrorString(err));
//     fprintf(stderr, "The problem is maybe -arch sm_13 instead of -arch sm_11 in the Makefile, please doublecheck\n");
//     exit(1);
//   }

//   err = cudaGetSymbolAddress((void**)&(mp->d_hprime_zz),"d_hprime_zz");
//   if(err != cudaSuccess) {
//     fprintf(stderr, "Error with d_hprime_zz: %s\n", cudaGetErrorString(err));
//     exit(1);
//   }
// }


void setConst_hprimewgll_xx(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_hprimewgll_xx, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprimewgll_xx: %s\n", cudaGetErrorString(err));
    exit(1);
  }

#ifdef USE_OLDER_CUDA4_GPU
  err = cudaGetSymbolAddress((void**)&(mp->d_hprimewgll_xx),"d_hprimewgll_xx");
#else
  err = cudaGetSymbolAddress((void**)&(mp->d_hprimewgll_xx),d_hprimewgll_xx);
#endif
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprimewgll_xx: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}


/*
// only needed if NGLLX != NGLLY != NGLLZ
void setConst_hprimewgll_zz(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_hprimewgll_zz, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprimewgll_zz: %s\n", cudaGetErrorString(err));
    exit(1);
  }

  err = cudaGetSymbolAddress((void**)&(mp->d_hprimewgll_zz),"d_hprimewgll_zz");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprimewgll_zz: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}
*/



#endif //CUDA_HEADER_H
