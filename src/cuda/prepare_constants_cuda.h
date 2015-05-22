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

__device__ realw d_hprimewgll_xx[NGLL2];

__device__ realw d_wxgll[NGLLX];


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
