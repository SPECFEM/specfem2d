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


/* ----------------------------------------------------------------------------------------------- */

// ELASTIC simulations

/* ----------------------------------------------------------------------------------------------- */

//fortran code snippet...
/*
  ! gets global number of that receiver
  irec = number_receiver_global(irec_local)

  ! gets local receiver interpolators
  ! (1-D Lagrange interpolators)
  hxir(:) = hxir_store(irec_local,:)
  hetar(:) = hetar_store(irec_local,:)
  hgammar(:) = hgammar_store(irec_local,:)
*/

/* ----------------------------------------------------------------------------------------------- */

// unused...
/*
__device__ double my_atomicAdd(double* address, double val) {

    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do{
      assumed = old;
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
*/

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_elastic_seismogram_kernel(int nrec_local,
                                                  realw* field,
                                                  int* d_ibool,
                                                  realw* hxir, realw* hgammar,
                                                  realw* seismograms,
                                                  realw* cosrot,
                                                  realw* sinrot,
                                                  int* ispec_selected_rec_loc,
                                                  int it,
                                                  int NSTEP)
{


  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  __shared__ realw sh_dxd[NGLL2_PADDED];
  __shared__ realw sh_dzd[NGLL2_PADDED];


  if (irec_local < nrec_local) {

    int ispec = ispec_selected_rec_loc[irec_local] - 1;

    sh_dxd[tx] = 0;
    sh_dzd[tx] = 0;

    if (tx < NGLL2) {
      realw hlagrange = hxir[irec_local + nrec_local*I] * hgammar[irec_local + nrec_local*J];
      int iglob = d_ibool[tx+NGLL2_PADDED*ispec] - 1;

      sh_dxd[tx] = hlagrange * field[0 + 2*iglob];
      sh_dzd[tx] = hlagrange * field[1 + 2*iglob];

      //debug
      //if (tx == 0) printf("thread %d %d %d - %f %f %f\n",ispec,iglob,irec_local,hlagrange,field[0 + 2*iglob],field[1 + 2*iglob]);
    }
    __syncthreads();

    // reduction
    for (unsigned int s=1; s<NGLL2_PADDED ; s *= 2) {
      if (tx % (2*s) == 0){ sh_dxd[tx] += sh_dxd[tx + s];sh_dzd[tx] += sh_dzd[tx + s];}
      __syncthreads();
    }

    if (tx == 0) {seismograms[irec_local*NSTEP + it]                    = cosrot[irec_local]*sh_dxd[0]  + sinrot[irec_local]*sh_dzd[0];}
    if (tx == 1) {seismograms[irec_local*NSTEP + it + nrec_local*NSTEP] = cosrot[irec_local]*sh_dzd[0]  - sinrot[irec_local]*sh_dxd[0];}

    /*
    // simple, single-thread reduction
    if (tx == 0) {
      // a loop in thread 0 is faster than atomic operations
      for(int s=1;s<NGLL2;s++) {
        sh_dxd[0] += sh_dxd[s];
        sh_dzd[0] += sh_dzd[s];
      }

      // rotate seismogram components
      seismograms[irec_local]            =    cosrot[irec_local]*sh_dxd[0] + sinrot[irec_local]*sh_dzd[0];
      seismograms[irec_local+nrec_local] =  - sinrot[irec_local]*sh_dxd[0] + cosrot[irec_local]*sh_dzd[0];
    }
    */
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_acoustic_seismogram_kernel(int nrec_local,
                                                   realw* pressure,
                                                   int* d_ibool,
                                                   realw* hxir, realw* hgammar,
                                                   realw* seismograms,
                                                   int* ispec_selected_rec_loc,
                                                   int it,
                                                   int NSTEP){
  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  __shared__ realw sh_dxd[NGLL2_PADDED];

  if (irec_local < nrec_local) {

    int ispec = ispec_selected_rec_loc[irec_local]-1;

    sh_dxd[tx] = 0;

    if (tx < NGLL2) {

      int iglob = d_ibool[tx+NGLL2_PADDED*ispec]-1;

      realw hlagrange = hxir[irec_local + nrec_local*I]*hgammar[irec_local + nrec_local*J];
      sh_dxd[tx] = hlagrange*pressure[iglob];
    }
    __syncthreads();

    for (unsigned int s=1; s<NGLL2_PADDED ; s *= 2) {
      if (tx % (2*s) == 0) sh_dxd[tx] += sh_dxd[tx + s];
      __syncthreads();
    }

    // Signe moins car pression = -potential_dot_dot
   if (tx == 0) {seismograms[irec_local*NSTEP + it ]                   = -sh_dxd[0];}
   if (tx == 1) {seismograms[irec_local*NSTEP + it + nrec_local*NSTEP] = 0;}
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_seismograms_cuda,
              COMPUTE_SEISMOGRAMS_CUDA)(long* Mesh_pointer_f,
                                        int* i_sigf,
                                        double* sisux, double* sisuz,
                                        int* seismo_currentf,
                                        int* nlength_seismogramf,
                                        int* ELASTIC_SIMULATION,
                                        int* ACOUSTIC_SIMULATION,
                                        int* USE_TRICK_FOR_BETTER_PRESSURE,
                                        int* itf,
                                        int* it_endf) {

// compute_seismograms
  TRACE("compute_seismograms_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  // synchronizes to wait for cuda computations to be completed before copying wavefield values
  synchronize_cuda();

  //checks if anything to do
  if (mp->nrec_local == 0) return;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nrec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLL2_PADDED,1,1);

  int i_sig = *i_sigf - 1;
  int seismo_current = *seismo_currentf - 1 ;
  int nlength_seismogram = *nlength_seismogramf;

  int seismotype = mp->h_seismotypeVec[i_sig];

  // warnings
  if (seismo_current == 0){
    switch (seismotype){
    case 1 :
      //Deplacement
      // warnings
      if (! *ELASTIC_SIMULATION)
        printf("\nWarning: Wrong type of seismogram for a pure fluid simulation, use pressure in seismotype\n");
      if (*ELASTIC_SIMULATION && *ACOUSTIC_SIMULATION)
        printf("\nWarning: Coupled elastic/fluid simulation has only valid displacement seismograms in elastic domain for GPU simulation\n\n");
      break;
    case 2 :
      //Vitesse
      if (! *ELASTIC_SIMULATION)
        printf("\nWarning: Wrong type of seismogram for a pure fluid simulation, use pressure in seismotype\n");
      if (*ELASTIC_SIMULATION && *ACOUSTIC_SIMULATION)
        printf("\nWarning: Coupled elastic/fluid simulation has only valid velocity seismograms in elastic domain for GPU simulation\n\n");
      break;
    case 3 :
      //Acceleration
      if (! *ELASTIC_SIMULATION)
        printf("\nWarning: Wrong type of seismogram for a pure fluid simulation, use pressure in seismotype\n");
      if (*ELASTIC_SIMULATION && *ACOUSTIC_SIMULATION)
        printf("\nWarning: Coupled elastic/fluid simulation has only valid acceleration seismograms in elastic domain for GPU simulation\n\n");
      break;
    case 4 :
      //Pression
      if (! *ACOUSTIC_SIMULATION)
        printf("\nWarning: Wrong type of seismogram for a pure elastic simulation, use displ veloc or accel in seismotype\n");
      if (*ELASTIC_SIMULATION && *ACOUSTIC_SIMULATION)
        printf("\nWarning: Coupled elastic/fluid simulation has only valid pressure seismograms in fluid domain for GPU simulation\n\n");
      break;
    }
  }

  // setup for copy
  realw* h_seismo = mp->h_seismograms[i_sig];
  realw* d_seismo = mp->d_seismograms[i_sig];

  // todo: for coupled simulations, one should check in which domain the receiver lies to output displacement
  //       similar to what routine compute_vector_one_element(..) is doing

  // computes current seismograms value
  switch (seismotype){
    case 1 :
      //Deplacement
      compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                               mp->d_displ,
                                                                               mp->d_ibool,
                                                                               mp->d_xir_store_loc, mp->d_gammar_store_loc,
                                                                               d_seismo,
                                                                               mp->d_cosrot,
                                                                               mp->d_sinrot,
                                                                               mp->d_ispec_selected_rec_loc,
                                                                               seismo_current,
                                                                               nlength_seismogram);
      break;

    case 2 :
      //Vitesse
      compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                               mp->d_veloc,
                                                                               mp->d_ibool,
                                                                               mp->d_xir_store_loc, mp->d_gammar_store_loc,
                                                                               d_seismo,
                                                                               mp->d_cosrot,
                                                                               mp->d_sinrot,
                                                                               mp->d_ispec_selected_rec_loc,
                                                                               seismo_current,
                                                                               nlength_seismogram);
      break;

    case 3 :
      //Acceleration
      compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                               mp->d_accel,
                                                                               mp->d_ibool,
                                                                               mp->d_xir_store_loc, mp->d_gammar_store_loc,
                                                                               d_seismo,
                                                                               mp->d_cosrot,
                                                                               mp->d_sinrot,
                                                                               mp->d_ispec_selected_rec_loc,
                                                                               seismo_current,
                                                                               nlength_seismogram);
      break;

    case 4 :
      //Pression
      if (*USE_TRICK_FOR_BETTER_PRESSURE){
        compute_acoustic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                  mp->d_potential_acoustic,
                                                                                  mp->d_ibool,
                                                                                  mp->d_xir_store_loc, mp->d_gammar_store_loc,
                                                                                  d_seismo,
                                                                                  mp->d_ispec_selected_rec_loc,
                                                                                  seismo_current,
                                                                                  nlength_seismogram);

      }else{
        compute_acoustic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                  mp->d_potential_dot_dot_acoustic,
                                                                                  mp->d_ibool,
                                                                                  mp->d_xir_store_loc, mp->d_gammar_store_loc,
                                                                                  d_seismo,
                                                                                  mp->d_ispec_selected_rec_loc,
                                                                                  seismo_current,
                                                                                  nlength_seismogram);

      }
      break;
  }//switch

  // note: due to subsampling, the last time step it == it_end might not be reached,
  //       but computing seismogram entries might end before.
  //       thus, both checks
  //         it%NTSTEP_BETWEEN_OUTPUT_SEISMOS == 0 || it == it_end
  //       might not be reached. instead we test if the seismogram array is full by
  //         seismo_current == nlength_seismogram - 1
  //       and copy it back whenever.

  // copies array to CPU host
  if (seismo_current == nlength_seismogram - 1 || it == it_end){
    // seismogram buffers are 1D and components appended; size for one single component record
    int size = mp->nrec_local * nlength_seismogram;

    // copies from GPU to CPU (note: could use async mem copy in future...)
    print_CUDA_error_if_any(cudaMemcpy(h_seismo, d_seismo, sizeof(realw) * 2 * size, cudaMemcpyDeviceToHost),72001);

    // copies values into host array
    for (int irec=0; irec < mp->nrec_local; irec++){
      for (int j=0; j < nlength_seismogram; j++){
        sisux[j + nlength_seismogram * irec] = (double) h_seismo[j + nlength_seismogram * irec];
        sisuz[j + nlength_seismogram * irec] = (double) h_seismo[j + nlength_seismogram * irec + size];
      }
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after compute_seismograms_cuda");
#endif
}
