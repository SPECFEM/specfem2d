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
                                                  realw* field_acoustic,
                                                  int* d_ibool,
                                                  realw* hxir, realw* hgammar,
                                                  realw* seismograms,
                                                  realw* cosrot,
                                                  realw* sinrot,
                                                  int* ispec_selected_rec_loc,
                                                  int* ispec_is_elastic,
                                                  int* ispec_is_acoustic,
                                                  realw* rhostore,
                                                  realw* d_hprime_xx,
                                                  realw* d_xix,realw* d_xiz,
                                                  realw* d_gammax,realw* d_gammaz,
                                                  int it,
                                                  int NSTEP){


  int irec_local = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  int J = (tx/NGLLX);
  int I = (tx-J*NGLLX);

  __shared__ realw sh_dxd[NGLL2_PADDED];
  __shared__ realw sh_dzd[NGLL2_PADDED];
  __shared__ realw scalar_field[NGLL2];

  if (irec_local < nrec_local) {
    // initializes
    sh_dxd[tx] = 0.0f;
    sh_dzd[tx] = 0.0f;

    // receiver element
    int ispec = ispec_selected_rec_loc[irec_local] - 1;

    // elastic domains
    if (ispec_is_elastic[ispec]) {
      if (tx < NGLL2) {
        realw hlagrange = hxir[irec_local + nrec_local*I] * hgammar[irec_local + nrec_local*J];
        int iglob = d_ibool[tx+NGLL2_PADDED*ispec] - 1;

        sh_dxd[tx] = hlagrange * field[2*iglob];
        sh_dzd[tx] = hlagrange * field[2*iglob+1];

        //debug
        //if (tx == 0) printf("thread %d %d %d - %f %f %f\n",ispec,iglob,irec_local,hlagrange,field[0 + 2*iglob],field[1 + 2*iglob]);
      }
    }

    // acoustic domains
    if (ispec_is_acoustic[ispec]) {
      // loads scalar into shared memory
      if (tx < NGLL2) {
        int iglob = d_ibool[tx+NGLL2_PADDED*ispec] - 1;
        scalar_field[tx] = field_acoustic[iglob];
      }
      // synchronizes threads
      __syncthreads();

      if (tx < NGLL2) {
        // compute gradient of potential to calculate vector if acoustic element
        // we then need to divide by density because the potential is a potential of (density * displacement)

        // gets material parameter
        int ij_ispec_padded = tx + NGLL2_PADDED*ispec;
        realw rhol = rhostore[ij_ispec_padded];

        // vector from scalar field
        realw vec_elem[2];

        compute_gradient_kernel(tx,ispec,scalar_field,vec_elem,
                                d_hprime_xx,
                                d_xix,d_xiz,d_gammax,d_gammaz,
                                rhol);

        realw hlagrange = hxir[irec_local + nrec_local*I] * hgammar[irec_local + nrec_local*J];

        sh_dxd[tx] = hlagrange * vec_elem[0];
        sh_dzd[tx] = hlagrange * vec_elem[1];

        //debug
        //if (tx == 0) printf("thread %d %d %d - %f %f %f\n",ispec,iglob,irec_local,hlagrange,field[0 + 2*iglob],field[1 + 2*iglob]);
      }
    }

    // synchronizes threads
    __syncthreads();

    // reduction
    for (unsigned int s=1; s<NGLL2_PADDED ; s *= 2) {
      if (tx % (2*s) == 0){
        sh_dxd[tx] += sh_dxd[tx + s];
        sh_dzd[tx] += sh_dzd[tx + s];
      }
      __syncthreads();
    }

    // component rotation
    if (tx == 0) {
      seismograms[irec_local*NSTEP + it] = cosrot[irec_local]*sh_dxd[0]  + sinrot[irec_local]*sh_dzd[0];
    }
    if (tx == 1) {
      seismograms[irec_local*NSTEP + it + nrec_local*NSTEP] = cosrot[irec_local]*sh_dzd[0] - sinrot[irec_local]*sh_dxd[0];
    }

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

