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
! the Free Software Foundation; either version 2 of the License, or
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

#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#include <sys/types.h>
#include <unistd.h>

#include "config.h"
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
                                                         int* number_receiver_global,
                                                         int* ispec_selected_rec) {


  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;
  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  __shared__ realw sh_dxd[NGLL2_PADDED];
  __shared__ realw sh_dzd[NGLL2_PADDED];


  if (irec_local < nrec_local) {

    int irec = number_receiver_global[irec_local]-1;
    int ispec = ispec_selected_rec[irec]-1;

   sh_dxd[tx] = 0;
   sh_dzd[tx] = 0;


  if (tx < NGLL2) {

    int iglob = d_ibool[tx+NGLL2_PADDED*ispec]-1;

    realw hlagrange = hxir[irec_local + nrec_local*I]*hgammar[irec_local + nrec_local*J];
    sh_dxd[tx] = hlagrange*field[0+2*iglob];
    sh_dzd[tx] = hlagrange*field[1+2*iglob];
    __syncthreads();}

for (unsigned int s=1; s<NGLL2_PADDED ; s *= 2) {
  if (tx % (2*s) == 0){ sh_dxd[tx] += sh_dxd[tx + s];sh_dzd[tx] += sh_dzd[tx + s];}
  __syncthreads();
}

  if (tx == 0) {seismograms[irec_local] = cosrot[irec_local]*sh_dxd[0]  + sinrot[irec_local]*sh_dzd[0];}
  if (tx == 1) {seismograms[irec_local+nrec_local] = cosrot[irec_local]*sh_dzd[0]  - sinrot[irec_local]*sh_dxd[0];}
}

}
/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_acoustic_seismogram_kernel(int nrec_local,
                                                         realw* pressure,
                                                         int* d_ibool,
                                                         realw* hxir, realw* hgammar,
                                                         realw* seismograms,
                                                         int* number_receiver_global,
                                                         int* ispec_selected_rec) {
  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;
  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  __shared__ realw sh_dxd[NGLL2_PADDED];



  if (irec_local < nrec_local) {

    int irec = number_receiver_global[irec_local]-1;
    int ispec = ispec_selected_rec[irec]-1;

   sh_dxd[tx] = 0;
realw hlagrange;
int iglob;
  if (tx < NGLL2) {

    iglob = d_ibool[tx+NGLL2_PADDED*ispec]-1;

    hlagrange = hxir[irec_local + nrec_local*I]*hgammar[irec_local + nrec_local*J];
    sh_dxd[tx] = hlagrange*pressure[iglob];
    __syncthreads();}

for (unsigned int s=1; s<NGLL2_PADDED ; s *= 2) {
  if (tx % (2*s) == 0) sh_dxd[tx] += sh_dxd[tx + s];
  __syncthreads();}


// Signe moins car pression = -minus_pressure
  if (tx == 0) {seismograms[irec_local] = -sh_dxd[0];}
  if (tx == 1) {seismograms[irec_local+nrec_local] = 0;}

    }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_seismograms_cuda,
              COMPUTE_SEISMOGRAMS_CUDA)(long* Mesh_pointer_f,
                                        int* seismotypef,
                                        double* sisux, double* sisuz,
                                        int* seismo_currentf,
                                        int* NSTEP_BETWEEN_OUTPUT_SEISMOSf,
                                        int * ELASTIC_SIMULATION,int * ACOUSTIC_SIMULATION,
                                        int* USE_TRICK_FOR_BETTER_PRESSURE) {

// compute_seismograms
  TRACE("\tcompute_seismograms");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper
  int seismotype = *seismotypef;
  int seismo_current = *seismo_currentf - 1 ;
  int NSTEP_BETWEEN_OUTPUT_SEISMOS = *NSTEP_BETWEEN_OUTPUT_SEISMOSf;
  int i;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nrec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLL2_PADDED,1,1);

  switch (seismotype){

  case 1 :  //Deplacement
  if (! *ELASTIC_SIMULATION) printf("\nWrong type of seismogram for a pure fluid simulation, use pressure in seismotype\n");

  compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(       mp->nrec_local,
                                                                                  mp->d_displ,
                                                                                  mp->d_ibool,
                                                                                  mp->d_xir_store_loc, mp->d_gammar_store_loc,
                                                                                  mp->d_seismograms,
                                                                                  mp->d_cosrot,
                                                                                  mp->d_sinrot,
                                                                                  mp->d_number_receiver_global,
                                                                                  mp->d_ispec_selected_rec
                                                                                  );

  break;

  case 2 :  //Vitesse
  if (! *ELASTIC_SIMULATION) printf("\nWrong type of seismogram for a pure fluid simulation, use pressure in seismotype\n");

  compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(       mp->nrec_local,
                                                                                  mp->d_veloc,
                                                                                  mp->d_ibool,
                                                                                  mp->d_xir_store_loc, mp->d_gammar_store_loc,
                                                                                  mp->d_seismograms,
                                                                                  mp->d_cosrot,
                                                                                  mp->d_sinrot,
                                                                                  mp->d_number_receiver_global,
                                                                                  mp->d_ispec_selected_rec
                                                                                  );
  break;

  case 3 :  //Acceleration
  if (! *ELASTIC_SIMULATION) printf("\nWrong type of seismogram for a pure fluid simulation, use pressure in seismotype\n");

  compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(       mp->nrec_local,
                                                                                  mp->d_accel,
                                                                                  mp->d_ibool,
                                                                                  mp->d_xir_store_loc, mp->d_gammar_store_loc,
                                                                                  mp->d_seismograms,
                                                                                  mp->d_cosrot,
                                                                                  mp->d_sinrot,
                                                                                  mp->d_number_receiver_global,
                                                                                  mp->d_ispec_selected_rec
                                                                                  );
  break;

  case 4 :  //Pression
  if (! *ACOUSTIC_SIMULATION) printf("\nWrong type of seismogram for a pure elasticsimulation, use displ veloc or accel in seismotype\n");

  if (*USE_TRICK_FOR_BETTER_PRESSURE)
  compute_acoustic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(      mp->nrec_local,
                                                                                  mp->d_minus_int_int_pressure_acoustic,
                                                                                  mp->d_ibool,
                                                                                  mp->d_xir_store_loc, mp->d_gammar_store_loc,
                                                                                  mp->d_seismograms,
                                                                                  mp->d_number_receiver_global,
                                                                                  mp->d_ispec_selected_rec
                                                                                  );
  else
  compute_acoustic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(      mp->nrec_local,
                                                                                  mp->d_minus_pressure_acoustic,
                                                                                  mp->d_ibool,
                                                                                  mp->d_xir_store_loc, mp->d_gammar_store_loc,
                                                                                  mp->d_seismograms,
                                                                                  mp->d_number_receiver_global,
                                                                                  mp->d_ispec_selected_rec
                                                                                  );

  break;
  }

  int size = mp->nrec_local;

  // (cudaMemcpy implicitly synchronizes all other cuda operations)
  print_CUDA_error_if_any(cudaMemcpy(mp->h_seismograms,mp->d_seismograms,
                                    sizeof(realw)*2* size,cudaMemcpyDeviceToHost),72001);

  for (i=0;i<size;i++)
   { sisux[seismo_current + NSTEP_BETWEEN_OUTPUT_SEISMOS * i ] = (double)*(mp->h_seismograms+i);
     sisuz[seismo_current + NSTEP_BETWEEN_OUTPUT_SEISMOS * i ] = (double)*(mp->h_seismograms+i+size);
   }

}
