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
                                        int* ATTENUATION_VISCOELASTIC,
                                        int* itf,
                                        int* it_endf) {

  // compute_seismograms
  TRACE("compute_seismograms_cuda");

  // flag to indicate that traces for kernel runs are taken from backward/reconstructed wavefields instead of adjoint wavefields;
  // useful for debugging.
  // default (0) is to output adjoint wavefield
  const int OUTPUT_BACKWARD_WAVEFIELD = 0;

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

  // warnings (deprecated)
  //if (seismo_current == 0){
  if (1 == 0){
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
  // warnings
  if (seismotype == 4 && seismo_current == 0){
    if (*ELASTIC_SIMULATION && ! mp->p_sv)
      printf("\nWarning: Pressure seismograms in elastic domain for GPU simulation only valid for P_SV cases\n\n");
    if (*ELASTIC_SIMULATION && *ATTENUATION_VISCOELASTIC)
      printf("\nWarning: Pressure seismograms in elastic domain for GPU simulation only valid for non-attenuation P_SV cases\n\n");
  }

  // setup for copy
  realw* h_seismo = mp->h_seismograms[i_sig];
  realw* d_seismo = mp->d_seismograms[i_sig];

  // selects field
  realw *displ;
  realw *potential;

  if (seismotype == 1){
    // deplacement
    if (OUTPUT_BACKWARD_WAVEFIELD && mp->simulation_type == 3){
      displ = mp->d_b_displ;
      potential = mp->d_b_potential_acoustic;
    }else{
      displ = mp->d_displ;
      potential = mp->d_potential_acoustic;
    }
  }else if (seismotype == 2){
    // vitesse
    if (OUTPUT_BACKWARD_WAVEFIELD && mp->simulation_type == 3){
      displ = mp->d_b_veloc;
      potential = mp->d_b_potential_dot_acoustic;
    }else{
      displ = mp->d_veloc;
      potential = mp->d_potential_dot_acoustic;
    }
  }else if (seismotype == 3){
    // acceleration
    if (OUTPUT_BACKWARD_WAVEFIELD && mp->simulation_type == 3){
      displ = mp->d_b_accel;
      potential = mp->d_b_potential_dot_dot_acoustic;
    }else{
      displ = mp->d_accel;
      potential = mp->d_potential_dot_dot_acoustic;
    }
  }else if (seismotype == 4){
    // pression
    if (OUTPUT_BACKWARD_WAVEFIELD && mp->simulation_type == 3){
      displ = mp->d_b_displ;
      if (*USE_TRICK_FOR_BETTER_PRESSURE){
        potential = mp->d_b_potential_acoustic;
      }else{
        potential = mp->d_b_potential_dot_dot_acoustic;
      }
    }else{
      displ = mp->d_displ;
      if (*USE_TRICK_FOR_BETTER_PRESSURE){
        potential = mp->d_potential_acoustic;
      }else{
        potential = mp->d_potential_dot_dot_acoustic;
      }
    }
  }

  // computes current seismograms value
  switch (seismotype){
    case 1 :
    case 2 :
    case 3 :
      //Displ/Veloc/Accel
      compute_elastic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                               displ,
                                                                               potential,
                                                                               mp->d_ibool,
                                                                               mp->d_xir_store_loc, mp->d_gammar_store_loc,
                                                                               d_seismo,
                                                                               mp->d_cosrot,
                                                                               mp->d_sinrot,
                                                                               mp->d_ispec_selected_rec_loc,
                                                                               mp->d_ispec_is_elastic,
                                                                               mp->d_ispec_is_acoustic,
                                                                               mp->d_rhostore,
                                                                               mp->d_hprime_xx,
                                                                               mp->d_xix,mp->d_xiz,
                                                                               mp->d_gammax,mp->d_gammaz,
                                                                               seismo_current,
                                                                               nlength_seismogram);
      break;

    case 4 :
      //Pression
      compute_acoustic_seismogram_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                                displ,
                                                                                potential,
                                                                                mp->d_ibool,
                                                                                mp->d_xir_store_loc, mp->d_gammar_store_loc,
                                                                                d_seismo,
                                                                                mp->d_ispec_selected_rec_loc,
                                                                                mp->d_ispec_is_elastic,
                                                                                mp->d_ispec_is_acoustic,
                                                                                mp->d_kappav,mp->d_muv,
                                                                                mp->d_hprime_xx,
                                                                                mp->d_xix,mp->d_xiz,
                                                                                mp->d_gammax,mp->d_gammaz,
                                                                                mp->d_ispec_is_anisotropic,
                                                                                mp->d_c11store,mp->d_c12store,mp->d_c13store,
                                                                                mp->d_c15store,mp->d_c23store,mp->d_c25store,
                                                                                mp->d_c33store,mp->d_c35store,
                                                                                seismo_current,
                                                                                nlength_seismogram);
      break;
  }//switch

  // note: due to subsampling, the last time step it == it_end might not be reached,
  //       but computing seismogram entries might end before.
  //       thus, both checks
  //         it%NTSTEP_BETWEEN_OUTPUT_SEISMOS == 0 || it == it_end
  //       might not be reached. instead we test if the seismogram array is full by
  //         seismo_current == nlength_seismogram - 1
  //       and copy it back whenever.

  int it = *itf;
  int it_end = *it_endf;

  // copies array to CPU host
  if (seismo_current == nlength_seismogram - 1 || it == it_end){
    // waits until previous compute stream finishes
    cudaStreamSynchronize(mp->compute_stream);

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

 GPU_ERROR_CHECKING ("compute_seismograms_cuda");
}
