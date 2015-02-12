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

#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// ELASTIC SIMULATIONS

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_kernels_ani_cudakernel(int* ispec_is_elastic,
                                           int* d_ibool,
                                           realw* accel,
                                           realw* b_displ,
                                           realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                           realw* epsilondev_xz,realw* epsilondev_yz,
                                           realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,
                                           realw* b_epsilondev_xz,realw* b_epsilondev_yz,
                                           realw* rho_kl,
                                           realw deltat,
                                           realw* cijkl_kl,
                                           realw* epsilon_trace_over_3,
                                           realw* b_epsilon_trace_over_3,
                                           int NSPEC_AB) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk = threadIdx.x;
  int ijk_ispec = ijk + NGLL3*ispec;
  int ijk21_ispec = ijk + 21*NGLL3*ispec;

  realw prod[21];
  realw eps[6];
  realw b_eps[6];
  realw epsdev[6];
  realw b_epsdev[6];
  realw eps_trace_over_3,b_eps_trace_over_3;
  int i,j;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC_AB) {

    // elastic elements only
    if( ispec_is_elastic[ispec] ) {
      int iglob = d_ibool[ijk + NGLL3_PADDED*ispec] - 1;

      // anisotropic kernels:
      // density kernel
      rho_kl[ijk_ispec] += deltat * (accel[3*iglob]*b_displ[3*iglob]+
                                     accel[3*iglob+1]*b_displ[3*iglob+1]+
                                     accel[3*iglob+2]*b_displ[3*iglob+2]);


      // anisotropic kernel
      epsdev[0] = epsilondev_xx[ijk_ispec];
      epsdev[1] = epsilondev_yy[ijk_ispec];
      epsdev[2] = epsilondev_xy[ijk_ispec];
      epsdev[3] = epsilondev_xz[ijk_ispec];
      epsdev[4] = epsilondev_yz[ijk_ispec];

      b_epsdev[0] = b_epsilondev_xx[ijk_ispec];
      b_epsdev[1] = b_epsilondev_yy[ijk_ispec];
      b_epsdev[2] = b_epsilondev_xy[ijk_ispec];
      b_epsdev[3] = b_epsilondev_xz[ijk_ispec];
      b_epsdev[4] = b_epsilondev_yz[ijk_ispec];

      eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec];
      b_eps_trace_over_3 = b_epsilon_trace_over_3[ijk_ispec];

      //! Building of the local matrix of the strain tensor
      //! for the adjoint field and the regular backward field
      //!eps11 et eps22
      eps[0] = epsdev[0] + eps_trace_over_3;
      eps[1] = epsdev[1] + eps_trace_over_3;
      //!eps33
      eps[2] = -(eps[0]+eps[1])+3*eps_trace_over_3;
      //!eps23
      eps[3] = epsdev[4];
      //!eps13
      eps[4] = epsdev[3];
      //!eps12
      eps[5] = epsdev[2];

      // backward arrays
      b_eps[0] = b_epsdev[0] + b_eps_trace_over_3;
      b_eps[1] = b_epsdev[1] + b_eps_trace_over_3;
      b_eps[2] = -(b_eps[0]+b_eps[1])+3*b_eps_trace_over_3;
      b_eps[3] = b_epsdev[4];
      b_eps[4] = b_epsdev[3];
      b_eps[5] = b_epsdev[2];

      //! Computing the 21 strain products without assuming eps(i)*b_eps(j) = eps(j)*b_eps(i)
      int p = 0;
      for( i=0; i<6; i++){
        for( j=i; j<6; j++){
          prod[p] = eps[i] * b_eps[j];
          if( j > i ){
            prod[p] = prod[p] + eps[j]*b_eps[i];
            if( j > 2 && i < 3 ){ prod[p] = prod[p]*2; }
          }
          if(i > 2 ){ prod[p] = prod[p]*4; }
          p++;
        }
      }

      // all 21 anisotropic coefficients
      for( i=0; i<21; i++){
        cijkl_kl[i+ijk21_ispec] += deltat * prod[i];
      }

    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_kernels_cudakernel(int* ispec_is_elastic,
                                           int* d_ibool,
                                           realw* accel,
                                           realw* b_displ,
                                           realw* rho_kl,
                                           realw* mu_kl,
                                           realw* kappa_kl,
                                           int NSPEC_AB,
                                           realw* dsxx,
                                           realw* dsxz,
                                           realw* dszz,
                                           realw* b_dsxx,
                                           realw* b_dsxz,
                                           realw* b_dszz,
                                           realw deltat) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ij = threadIdx.x;
  int ij_ispec = ij + NGLL2*ispec;


  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC_AB) {

    // elastic elements only
    if( ispec_is_elastic[ispec] ) {
      int iglob = d_ibool[ij + NGLL2_PADDED*ispec] - 1 ;

      // isotropic kernels:
      // density kernel
      rho_kl[ij_ispec] +=  (accel[2*iglob]*b_displ[2*iglob]+
                                     accel[2*iglob+1]*b_displ[2*iglob+1]);


      realw prod = (dsxx[iglob] + dszz[iglob]) *  (b_dsxx[iglob] + b_dszz[iglob]);

      // shear modulus kernel
      mu_kl[ij_ispec] += (dsxx[iglob] * b_dsxx[iglob] + dszz[iglob] * b_dszz[iglob] +
                           0.5* dsxz[iglob] * b_dsxz[iglob] -
                           prod/3) * deltat;

      // bulk modulus kernel
      kappa_kl[ij_ispec] += prod * deltat;

    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_kernels_elastic_cuda,
              COMPUTE_KERNELS_ELASTIC_CUDA)(long* Mesh_pointer,realw * deltat) {

  TRACE("compute_kernels_elastic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL2; // NGLLX*NGLLZ


  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_AB,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);


    compute_kernels_cudakernel<<<grid,threads>>>(mp->d_ispec_is_elastic,mp->d_ibool,
                                                 mp->d_accel, mp->d_b_displ,
                                                 mp->d_rho_kl,
                                                 mp->d_mu_kl,
                                                 mp->d_kappa_kl,
                                                 mp->NSPEC_AB,
                                                 mp->d_dsxx,
                                                 mp->d_dsxz,
                                                 mp->d_dszz,
                                                 mp->d_b_dsxx,
                                                 mp->d_b_dsxz,
                                                 mp->d_b_dszz,
                                                 *deltat);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_kernels_elastic_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// ACOUSTIC SIMULATIONS

/* ----------------------------------------------------------------------------------------------- */


__device__ void compute_gradient_kernel(int ij,
                                        int ispec,
                                        realw* scalar_field,
                                        realw* vector_field_element,
                                        realw* d_hprime_xx,
                                        realw* d_xix,realw* d_xiz,
                                        realw* d_gammax,realw* d_gammaz,
                                        realw rhol) {

  realw temp1l,temp3l;
  realw hp1,hp3;
  realw xixl,xizl,gammaxl,gammazl;
  realw rho_invl;
  int l,offset,offset1,offset3;

  const int NGLL2_ALIGN = NGLL2_PADDED;


  int J = (ij/NGLLX);
  int I = (ij-J*NGLLX);

  // derivative along x
  temp1l = 0.f;
  for( l=0; l<NGLLX;l++){
    hp1 = d_hprime_xx[l*NGLLX+I];
    offset1 = J*NGLLX+l;
    temp1l += scalar_field[offset1]*hp1;
  }

  // derivative along z
  temp3l = 0.f;
  for( l=0; l<NGLLX;l++){
    // assumes hprime_xx == hprime_zz
    hp3 = d_hprime_xx[l*NGLLX+J];
    offset3 = l*NGLLX+I;
    temp3l += scalar_field[offset3]*hp3;
  }

  offset = ispec*NGLL2_ALIGN + ij;

  xixl = d_xix[offset];
  xizl = d_xiz[offset];
  gammaxl = d_gammax[offset];
  gammazl = d_gammaz[offset];


    rho_invl = 1.0f / rhol;

  // derivatives of acoustic scalar potential field on GLL points
  vector_field_element[0] = (temp1l*xixl + temp3l*gammaxl) * rho_invl;
  vector_field_element[1] = (temp1l*xizl + temp3l*gammazl) * rho_invl;

}

/* ----------------------------------------------------------------------------------------------- */


__global__ void compute_kernels_acoustic_kernel(int* ispec_is_acoustic,
                                                int* d_ibool,
                                                realw* rhostore,
                                                realw* kappastore,
                                                realw* d_hprime_xx,
                                                realw* d_xix,realw* d_xiz,
                                                realw* d_gammax,realw* d_gammaz,
                                                realw* potential_dot_dot_acoustic,
                                                realw* b_potential_acoustic,
                                                realw* b_potential_dot_dot_acoustic,
                                                realw* rho_ac_kl,
                                                realw* kappa_ac_kl,
                                                int NSPEC_AB,
                                                realw deltat) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ij = threadIdx.x;

  // local and global indices
  int ij_ispec = ij + NGLL2*ispec;
  int ij_ispec_padded = ij + NGLL2_PADDED*ispec;
  int iglob;

  // shared memory between all threads within this block
  __shared__ realw scalar_field_displ[NGLL2];
  __shared__ realw scalar_field_accel[NGLL2];

  int active = 0;

  // handles case when there is 1 extra block (due to rectangular grid)
  if( ispec < NSPEC_AB ){
    // acoustic elements only
    if( ispec_is_acoustic[ispec] ){
      active = 1;

      // copy field values
      iglob = d_ibool[ij_ispec_padded] - 1;
      scalar_field_displ[ij] = b_potential_acoustic[iglob];
      scalar_field_accel[ij] = potential_dot_dot_acoustic[iglob];
    }
  }

  // synchronizes threads
  __syncthreads();

  if( active ){
    realw accel_elm[2];
    realw b_displ_elm[2];
    realw rhol,kappal;

    // gets material parameter
    rhol = rhostore[ij_ispec_padded];

    // displacement vector from backward field
    compute_gradient_kernel(ij,ispec,scalar_field_displ,b_displ_elm,
                            d_hprime_xx,
                            d_xix,d_xiz,d_gammax,d_gammaz,
                            rhol);

    // acceleration vector
    compute_gradient_kernel(ij,ispec,scalar_field_accel,accel_elm,
                            d_hprime_xx,
                            d_xix,d_xiz,d_gammax,d_gammaz,
                            rhol);

    // density kernel
    rho_ac_kl[ij_ispec] -=  rhol * (accel_elm[0]*b_displ_elm[0] +
                                             accel_elm[1]*b_displ_elm[1] ) * deltat;

    // bulk modulus kernel
    kappal = kappastore[ij_ispec];
    kappa_ac_kl[ij_ispec] -= potential_dot_dot_acoustic[iglob] * b_potential_acoustic[iglob] * deltat / kappal ;
  } // active
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_kernels_acoustic_cuda,
              COMPUTE_KERNELS_ACOUSTIC_CUDA)(long* Mesh_pointer,realw * deltat) {

TRACE("compute_kernels_acoustic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL2; // NGLLX*NGLLZ

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_AB,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  compute_kernels_acoustic_kernel<<<grid,threads>>>(mp->d_ispec_is_acoustic,
                                                    mp->d_ibool,
                                                    mp->d_rhostore,
                                                    mp->d_kappastore,
                                                    mp->d_hprime_xx,
                                                    mp->d_xix,mp->d_xiz,
                                                    mp->d_gammax,mp->d_gammaz,
                                                    mp->d_potential_dot_dot_acoustic,
                                                    mp->d_b_potential_acoustic,
                                                    mp->d_b_potential_dot_dot_acoustic,
                                                    mp->d_rho_ac_kl,
                                                    mp->d_kappa_ac_kl,
                                                    mp->NSPEC_AB,
                                                    *deltat);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_kernels_acoustic_kernel");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// preconditioner (approximate Hessian kernel)

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_kernels_hess_el_cudakernel(int* ispec_is_elastic,
                                                   int* d_ibool,
                                                   realw* accel,
                                                   realw* b_accel,
                                                   realw* hess_kl,
                                                   int NSPEC_AB) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ij = threadIdx.x;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC_AB) {

    // elastic elements only
    if( ispec_is_elastic[ispec] ) {
      int iglob = d_ibool[ij + NGLL2_PADDED*ispec] - 1;

      // approximate hessian
      hess_kl[ij + NGLL2*ispec] +=  (accel[2*iglob]*b_accel[2*iglob]+
                                              accel[2*iglob+1]*b_accel[2*iglob+1]);
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_kernels_hess_ac_cudakernel(int* ispec_is_acoustic,
                                                   int* d_ibool,
                                                   realw* potential_dot_dot_acoustic,
                                                   realw* b_potential_dot_dot_acoustic,
                                                   realw* rhostore,
                                                   realw* d_hprime_xx,
                                                   realw* d_xix,realw* d_xiz,
                                                   realw* d_gammax,realw* d_gammaz,
                                                   realw* hess_kl,
                                                   int NSPEC_AB) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ij = threadIdx.x;
  int ij_ispec_padded = ij + NGLL2_PADDED*ispec;
  int iglob;

  // shared memory between all threads within this block
  __shared__ realw scalar_field_accel[NGLL2];
  __shared__ realw scalar_field_b_accel[NGLL2];

  int active = 0;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC_AB) {

    // acoustic elements only
    if( ispec_is_acoustic[ispec] ){
      active = 1;

      // global indices
      iglob = d_ibool[ij_ispec_padded] - 1;

      // copy field values
      scalar_field_accel[ij] = potential_dot_dot_acoustic[iglob];
      scalar_field_b_accel[ij] = b_potential_dot_dot_acoustic[iglob];
    }
  }

  // synchronizes threads
  __syncthreads();

  if( active ){
    realw accel_elm[2];
    realw b_accel_elm[2];
    realw rhol;

    // gets material parameter
    rhol = rhostore[ij_ispec_padded];

    // acceleration vector
    compute_gradient_kernel(ij,ispec,
                            scalar_field_accel,accel_elm,
                            d_hprime_xx,
                            d_xix,d_xiz,d_gammax,d_gammaz,
                            rhol);

    // acceleration vector from backward field
    compute_gradient_kernel(ij,ispec,
                            scalar_field_b_accel,b_accel_elm,
                            d_hprime_xx,
                            d_xix,d_xiz,d_gammax,d_gammaz,
                            rhol);
    // approximates hessian
    hess_kl[ij + NGLL2*ispec] +=  (accel_elm[0]*b_accel_elm[0] +
                                            accel_elm[1]*b_accel_elm[1]);

  } // active
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_kernels_hess_cuda,
              COMPUTE_KERNELS_HESS_CUDA)(long* Mesh_pointer,
                                         int* ELASTIC_SIMULATION,
                                         int* ACOUSTIC_SIMULATION) {
  TRACE("compute_kernels_hess_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL2; // NGLLX*NGLLZ

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_AB,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if( *ELASTIC_SIMULATION ) {
    compute_kernels_hess_el_cudakernel<<<grid,threads>>>(mp->d_ispec_is_elastic,
                                                         mp->d_ibool,
                                                         mp->d_accel,
                                                         mp->d_b_accel,
                                                         mp->d_hess_el_kl,
                                                         mp->NSPEC_AB);
  }

  if( *ACOUSTIC_SIMULATION ) {
    compute_kernels_hess_ac_cudakernel<<<grid,threads>>>(mp->d_ispec_is_acoustic,
                                                         mp->d_ibool,
                                                         mp->d_potential_dot_dot_acoustic,
                                                         mp->d_b_potential_dot_dot_acoustic,
                                                         mp->d_rhostore,
                                                         mp->d_hprime_xx,
                                                         mp->d_xix,mp->d_xiz,
                                                         mp->d_gammax,mp->d_gammaz,
                                                         mp->d_hess_ac_kl,
                                                         mp->NSPEC_AB);
  }


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_kernels_hess_cuda");
#endif
}

