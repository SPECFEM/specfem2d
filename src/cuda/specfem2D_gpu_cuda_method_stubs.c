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
#include <stdlib.h>

#include "config.h"

typedef float realw;



//
// src/cuda/assemble_MPI_scalar_cuda.cu
//

void FC_FUNC_(transfer_boun_pot_from_device,
              TRANSFER_BOUN_POT_FROM_DEVICE)(long* Mesh_pointer,
                                             realw* send_potential_dot_dot_buffer,
                                             const int* FORWARD_OR_ADJOINT){}

void FC_FUNC_(transfer_asmbl_pot_to_device,
              TRANSFER_ASMBL_POT_TO_DEVICE)(long* Mesh_pointer,
                                            realw* buffer_recv_scalar_ext_mesh,
                                            const int* FORWARD_OR_ADJOINT) {}


//
// src/cuda/assemble_MPI_vector_cuda.cu
//

void FC_FUNC_(transfer_boun_accel_from_device,
              TRANSFER_BOUN_ACCEL_FROM_DEVICE)(long* Mesh_pointer,
                                               realw* send_accel_buffer,
                                               const int* FORWARD_OR_ADJOINT){}

void FC_FUNC_(transfer_boundary_from_device_a,
              TRANSFER_BOUNDARY_FROM_DEVICE_A)(long* Mesh_pointer) {}

void FC_FUNC_(prepare_boundary_on_device,
              PREPARE_BOUNDARY_ON_DEVICE)(long* Mesh_pointer) {}

void FC_FUNC_(transfer_boundary_to_device_a,
              TRANSFER_BOUNDARY_TO_DEVICE_A)(long* Mesh_pointer,
                                             realw* buffer_recv_vector_ext_mesh,
                                             const int* max_nibool_interfaces_ext_mesh) {}

void FC_FUNC_(transfer_asmbl_accel_to_device,
              TRANSFER_ASMBL_ACCEL_TO_DEVICE)(long* Mesh_pointer,
                                              realw* buffer_recv_vector_ext_mesh,
                                              const int* max_nibool_interfaces_ext_mesh,
                                              const int* nibool_interfaces_ext_mesh,
                                              const int* ibool_interfaces_ext_mesh,
                                              const int* FORWARD_OR_ADJOINT,const int *CUDA_AWARE_MPI) {}

//void FC_FUNC_(assemble_accel_on_device,
//              ASSEMBLE_ACCEL_on_DEVICE)(long* Mesh_pointer, realw* accel,
//                                              realw* buffer_recv_vector_ext_mesh,
//                                              int* num_interfaces_ext_mesh,
//                                              int* max_nibool_interfaces_ext_mesh,
//                                              int* nibool_interfaces_ext_mesh,
//                                              int* ibool_interfaces_ext_mesh,
//                                              int* FORWARD_OR_ADJOINT) {}


void FC_FUNC_(sync_copy_from_device,
              SYNC_copy_FROM_DEVICE)(long* Mesh_pointer,
                                     int* iphase,
                                     realw* send_buffer) {}


//
// src/cuda/check_fields_cuda.cu
//

void FC_FUNC_(pause_for_debug,PAUSE_FOR_DEBUG)() {}

void FC_FUNC_(output_free_device_memory,
              OUTPUT_FREE_DEVICE_MEMORY)(int* myrank_f) {}

void FC_FUNC_(get_free_device_memory,
              get_FREE_DEVICE_MEMORY)(realw* free, realw* used, realw* total ) {}

void FC_FUNC_(get_norm_acoustic_from_device,
              GET_NORM_ACOUSTIC_FROM_DEVICE)(realw* norm,long* Mesh_pointer,int* sim_type) {}

void FC_FUNC_(get_norm_elastic_from_device,
              GET_NORM_ELASTIC_FROM_DEVICE)(realw* norm,
                                            long* Mesh_pointer,
                                            int* type,int *it) {}

void FC_FUNC_(get_max_accel,
              GET_MAX_ACCEL)(int* itf,int* sizef,long* Mesh_pointer) {}

void FC_FUNC_(check_max_norm_displ_gpu,
              CHECK_MAX_NORM_DISPL_GPU)(int* size, realw* displ,long* Mesh_pointer,int* announceID) {}

void FC_FUNC_(check_max_norm_vector,
              CHECK_MAX_NORM_VECTOR)(int* size, realw* vector1, int* announceID) {}

void FC_FUNC_(check_max_norm_displ,
              CHECK_MAX_NORM_DISPL)(int* size, realw* displ, int* announceID) {}

void FC_FUNC_(check_max_norm_b_displ_gpu,
              CHECK_MAX_NORM_B_DISPL_GPU)(int* size, realw* b_displ,long* Mesh_pointer,int* announceID) {}

void FC_FUNC_(check_max_norm_b_accel_gpu,
              CHECK_MAX_NORM_B_ACCEL_GPU)(int* size, realw* b_accel,long* Mesh_pointer,int* announceID) {}

void FC_FUNC_(check_max_norm_b_veloc_gpu,
              CHECK_MAX_NORM_B_VELOC_GPU)(int* size, realw* b_veloc,long* Mesh_pointer,int* announceID) {}

void FC_FUNC_(check_max_norm_b_displ,
              CHECK_MAX_NORM_B_DISPL)(int* size, realw* b_displ,int* announceID) {}

void FC_FUNC_(check_max_norm_b_accel,
              CHECK_MAX_NORM_B_ACCEL)(int* size, realw* b_accel,int* announceID) {}

void FC_FUNC_(check_error_vectors,
              CHECK_ERROR_VECTORS)(int* sizef, realw* vector1,realw* vector2) {}


//
// src/cuda/compute_add_sources_acoustic_cuda.cu
//

void FC_FUNC_(compute_add_sources_ac_cuda,
              COMPUTE_ADD_SOURCES_AC_CUDA)(long* Mesh_pointer,
                                           int* phase_is_innerf,
                                           int* NSOURCESf,
                                           int * itf) {}

void FC_FUNC_(compute_add_sources_ac_s3_cuda,
              COMPUTE_ADD_SOURCES_AC_s3_CUDA)(long* Mesh_pointer,
                                              int* phase_is_innerf,
                                              int* NSOURCESf,
                                              double* h_stf_pre_compute) {}

void FC_FUNC_(add_sources_ac_sim_2_or_3_cuda,
              ADD_SOURCES_AC_SIM_2_OR_3_CUDA)(long* Mesh_pointer,
                                               realw* h_adj_sourcearrays,
                                               int* phase_is_inner,
                                               int* h_ispec_is_inner,
                                               int* h_ispec_is_acoustic,
                                               int* h_ispec_selected_rec,
                                               int* nrec,
                                               int* time_index,
                                               int* h_islice_selected_rec,
                                               int* nadj_rec_local,
                                               int* NTSTEP_BETWEEN_READ_ADJSRC) {}


//
// src/cuda/compute_add_sources_viscoelastic_cuda.cu
//

void FC_FUNC_(compute_add_sources_el_cuda,
              COMPUTE_ADD_SOURCES_EL_CUDA)(long* Mesh_pointer,
                                           double* h_stf_pre_compute,
                                           int* h_NSOURCES,
                                           int* h_phase_is_inner) {}

void FC_FUNC_(compute_add_sources_el_s3_cuda,
              COMPUTE_ADD_SOURCES_EL_S3_CUDA)(long* Mesh_pointer,
                                              int* phase_is_innerf,
                                              int* NSOURCESf,
                                              int* itf) {}

void FC_FUNC_(add_sources_el_sim_type_2_or_3,
              ADD_SOURCES_EL_SIM_TYPE_2_OR_3)(realw* potential_dot_dot_acoustic,
                                                      realw* source_adjointe,
                                                      realw* xir_store,
                                                      realw* gammar_store,
                                                      int* d_ibool,
                                                      int* ispec_is_inner,
                                                      int* ispec_is_acoustic,
                                                      int* ispec_selected_rec,
                                                      int phase_is_inner,
                                                      int it,
                                                      int* pre_computed_irec,
                                                      int nadj_rec_local,
                                                      realw* kappastore,
                                                      int NSTEP  ) {}


//
// src/cuda/compute_coupling_cuda.cu
//

void FC_FUNC_(compute_coupling_ac_el_cuda,
              COMPUTE_COUPLING_AC_EL_CUDA)(long* Mesh_pointer,
                                           int* phase_is_innerf,
                                           int* num_coupling_ac_el_facesf) {}

void FC_FUNC_(compute_coupling_el_ac_cuda,
              COMPUTE_COUPLING_EL_AC_CUDA)(long* Mesh_pointer,
                                           int* phase_is_innerf,
                                           int* num_coupling_ac_el_facesf) {}

//
// src/cuda/compute_forces_acoustic_cuda.cu
//

void FC_FUNC_(compute_forces_acoustic_cuda,
              COMPUTE_FORCES_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                            int* iphase,
                                            int* nspec_outer_acoustic,
                                            int* nspec_inner_acoustic) {}

void FC_FUNC_(acoustic_enforce_free_surf_cuda,
              ACOUSTIC_ENFORCE_FREE_SURF_CUDA)(long* Mesh_pointer) {}


//
// src/cuda/compute_forces_viscoelastic_cuda.cu
//

void FC_FUNC_(compute_forces_viscoelastic_cuda,
              COMPUTE_FORCES_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                                int* iphase,
                                                realw* deltat,
                                                int* nspec_outer_elastic,
                                                int* nspec_inner_elastic,
                                                int* ANISOTROPY) {}




//
// src/cuda/compute_kernels_cuda.cu
//

void FC_FUNC_(compute_kernels_elastic_cuda,
              COMPUTE_KERNELS_ELASTIC_CUDA)(long* Mesh_pointer) {}


void FC_FUNC_(compute_kernels_acoustic_cuda,
              COMPUTE_KERNELS_ACOUSTIC_CUDA)(long* Mesh_pointer) {}

void FC_FUNC_(compute_kernels_hess_cuda,
              COMPUTE_KERNELS_HESS_CUDA)(long* Mesh_pointer,
                                         int* ELASTIC_SIMULATION,
                                         int* ACOUSTIC_SIMULATION) {}


//
// src/cuda/compute_stacey_acoustic_cuda.cu
//

void FC_FUNC_(compute_stacey_acoustic_cuda,
              COMPUTE_STACEY_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                            int* phase_is_innerf,
                                            realw* h_b_absorb_potential_left,
                                            realw* h_b_absorb_potential_right,
                                            realw* h_b_absorb_potential_top,
                                            realw* h_b_absorb_potential_bottom) {}


//
// src/cuda/compute_stacey_viscoelastic_cuda.cu
//

void FC_FUNC_(compute_stacey_viscoelastic_cuda,
              COMPUTE_STACEY_VISCOELASTIC_CUDA)(long* Mesh_pointer,
                                           int* phase_is_innerf,
                                           realw* h_b_absorb_elastic_left,
                                           realw* h_b_absorb_elastic_right,
                                           realw* h_b_absorb_elastic_top,
                                           realw* h_b_absorb_elastic_bottom) {}


//
// src/cuda/initialize_cuda.cu
//

void FC_FUNC_(initialize_cuda_device,
              INITIALIZE_CUDA_DEVICE)(int* myrank_f,int* ncuda_devices) {
 fprintf(stderr,"ERROR: GPU_MODE enabled without GPU/CUDA Support. To enable GPU support, reconfigure with --with-cuda flag.\n");
 exit(1);
}

//
// src/cuda/prepare_mesh_constants_cuda.cu
//

void FC_FUNC_(prepare_constants_device,
              PREPARE_CONSTANTS_DEVICE)(long* Mesh_pointer,
                                        int* h_NGLLX, int* NSPEC_AB, int* NGLOB_AB,
                                        realw* h_xix, realw* h_xiz,
                                        realw* h_gammax, realw* h_gammaz,
                                        realw* h_kappav, realw* h_muv,
                                        int* h_ibool,
                                        int* num_interfaces_ext_mesh, int* max_nibool_interfaces_ext_mesh,
                                        int* h_nibool_interfaces_ext_mesh, int* h_ibool_interfaces_ext_mesh,
                                        realw* h_hprime_xx, realw* h_hprimewgll_xx,
                                        realw* h_wxgll,
                                        int* ABSORBING_CONDITIONS,
                                        int* h_nspec_bottom,
                                        int* h_nspec_left,
                                        int* h_nspec_right,
                                        int* h_nspec_top,
                                        int* h_abs_boundary_ispec, int* h_abs_boundary_ij,
                                        realw* h_abs_boundary_normal,
                                        realw* h_abs_boundary_jacobian1Dw,
                                        int* h_num_abs_boundary_faces,
                                        int* h_cote_abs,
                                        int* h_ib_bottom,
                                        int* h_ib_left,
                                        int* h_ib_right,
                                        int* h_ib_top,
                                        int* h_ispec_is_inner,
                                        int* NSOURCES, int* nsources_local_f, int* h_num_src_loc,
                                        realw* h_sourcearrays, realw * h_source_time_function,
                                        int* NSTEP,
                                        int* h_islice_selected_source, int* h_ispec_selected_source,
                                        int* h_number_receiver_global, int* h_ispec_selected_rec,
                                        int* nrec,int* nrec_local,
                                        realw * h_cosrot,realw * h_sinrot,
                                        int* SIMULATION_TYPE,
                                        int* USE_MESH_COLORING_GPU_f,
                                        int* nspec_acoustic,int* nspec_elastic,
                                        int* h_myrank,
                                        int* SAVE_FORWARD,
                                        realw* h_xir_store, realw* h_gammar_store ) {}

void FC_FUNC_(prepare_fields_acoustic_device,
              PREPARE_FIELDS_ACOUSTIC_DEVICE)(long* Mesh_pointer,
                                              realw* rmass_acoustic, realw* rhostore, realw* kappastore,
                                              int* num_phase_ispec_acoustic, int* phase_ispec_inner_acoustic,
                                              int* ispec_is_acoustic,
                                              int* num_free_surface_faces,
                                              int* free_surface_ispec,
                                              int* free_surface_ijk,
                                              int* ELASTIC_SIMULATION,
                                              int* num_coupling_ac_el_faces,
                                              int* coupling_ac_el_ispec,
                                              int* coupling_ac_el_ijk,
                                              realw* coupling_ac_el_normal,
                                              realw* coupling_ac_el_jacobian2Dw,
                                              int * h_ninterface_acoustic,int * h_inum_interfaces_acoustic,
                                              int* num_colors_outer_acoustic,
                                              int* num_colors_inner_acoustic,
                                              int* num_elem_colors_acoustic) {}

void FC_FUNC_(prepare_fields_acoustic_adj_dev,
              PREPARE_FIELDS_ACOUSTIC_ADJ_DEV)(long* Mesh_pointer,
                                              int* APPROXIMATE_HESS_KL) {}

void FC_FUNC_(prepare_fields_elastic_device,
              PREPARE_FIELDS_ELASTIC_DEVICE)(long* Mesh_pointer,
                                             realw* rmassx, realw* rmassz,
                                             realw* rho_vp, realw* rho_vs,
                                             int* num_phase_ispec_elastic,
                                             int* phase_ispec_inner_elastic,
                                             int* ispec_is_elastic,
                                             int* h_nspec_left,
                                             int* h_nspec_right,
                                             int* h_nspec_top,
                                             int* h_nspec_bottom,
                                             int* ACOUSTIC_SIMULATION,
                                             int* num_colors_outer_elastic,
                                             int* num_colors_inner_elastic,
                                             int* num_elem_colors_elastic,
                                             int* ANISOTROPY,
                                             realw *c11store,realw *c12store,realw *c13store,
                                             realw *c15store,
                                             realw *c23store,
                                             realw *c25store,realw *c33store,
                                             realw *c35store,
                                             realw *c55store,int* h_ninterface_elastic,int * h_inum_interfaces_elastic ){}

void FC_FUNC_(prepare_fields_elastic_adj_dev,
              PREPARE_FIELDS_ELASTIC_ADJ_DEV)(long* Mesh_pointer,
                                             int* size_f,
                                             int* APPROXIMATE_HESS_KL){}

void FC_FUNC_(prepare_sim2_or_3_const_device,
              PREPARE_SIM2_OR_3_CONST_DEVICE)(long* Mesh_pointer,
                                              int* islice_selected_rec,
                                              int* nadj_rec_local,
                                              int* nrec,realw* h_source_adjointe,int* NSTEP) {}

void FC_FUNC_(prepare_cleanup_device,
              PREPARE_CLEANUP_DEVICE)(long* Mesh_pointer,
                                      int* ACOUSTIC_SIMULATION,
                                      int* ELASTIC_SIMULATION,
                                      int* ABSORBING_CONDITIONS,
                                      int* ANISOTROPY,
                                      int* APPROXIMATE_HESS_KL) {}


//
// src/cuda/transfer_fields_cuda.cu
//

void FC_FUNC_(transfer_fields_el_to_device,
              TRANSFER_FIELDS_EL_TO_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer) {}

void FC_FUNC_(transfer_fields_el_from_device,
              TRANSFER_FIELDS_EL_FROM_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_fields_to_device,
              TRANSFER_B_FIELDS_TO_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                           long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_fields_from_device,
              TRANSFER_B_FIELDS_FROM_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,long* Mesh_pointer) {}

void FC_FUNC_(transfer_accel_to_device,
              TRNASFER_ACCEL_TO_DEVICE)(int* size, realw* accel,long* Mesh_pointer) {}

void FC_FUNC_(transfer_accel_from_device,
              TRANSFER_ACCEL_FROM_DEVICE)(int* size, realw* accel,long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_accel_from_device,
              TRNASFER_B_ACCEL_FROM_DEVICE)(int* size, realw* b_accel,long* Mesh_pointer) {}

void FC_FUNC_(transfer_sigma_from_device,
              TRANSFER_SIGMA_FROM_DEVICE)(int* size, realw* sigma_kl,long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_displ_from_device,
              TRANSFER_B_DISPL_FROM_DEVICE)(int* size, realw* displ,long* Mesh_pointer) {}

void FC_FUNC_(transfer_displ_from_device,
              TRANSFER_DISPL_FROM_DEVICE)(int* size, realw* displ,long* Mesh_pointer) {}


void FC_FUNC_(transfer_kernels_el_to_host,
              TRANSFER_KERNELS_EL_TO_HOST)(long* Mesh_pointer,
                                            realw* h_rho_kl,
                                            realw* h_mu_kl,
                                            realw* h_kappa_kl,
                                            int* NSPEC_AB) {}

void FC_FUNC_(transfer_fields_ac_to_device,
              TRANSFER_FIELDS_AC_TO_DEVICE)(int* size,
                                            realw* potential_acoustic,
                                            realw* potential_dot_acoustic,
                                            realw* potential_dot_dot_acoustic,
                                            long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_fields_ac_to_device,
              TRANSFER_B_FIELDS_AC_TO_DEVICE)(int* size,
                                              realw* b_potential_acoustic,
                                              realw* b_potential_dot_acoustic,
                                              realw* b_potential_dot_dot_acoustic,
                                              long* Mesh_pointer) {}

void FC_FUNC_(transfer_fields_ac_from_device,
              TRANSFER_FIELDS_AC_FROM_DEVICE)(int* size,
                                              realw* potential_acoustic,
                                              realw* potential_dot_acoustic,
                                              realw* potential_dot_dot_acoustic,
                                              long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_fields_ac_from_device,
              TRANSFER_B_FIELDS_AC_FROM_DEVICE)(int* size,
                                                realw* b_potential_acoustic,
                                                realw* b_potential_dot_acoustic,
                                                realw* b_potential_dot_dot_acoustic,
                                                long* Mesh_pointer) {}

void FC_FUNC_(transfer_dot_dot_from_device,
              TRNASFER_DOT_DOT_FROM_DEVICE)(int* size, realw* potential_dot_dot_acoustic,long* Mesh_pointer) {}

void FC_FUNC_(transfer_b_dot_dot_from_device,
              TRNASFER_B_DOT_DOT_FROM_DEVICE)(int* size, realw* b_potential_dot_dot_acoustic,long* Mesh_pointer) {}

void FC_FUNC_(transfer_kernels_ac_to_host,
              TRANSFER_KERNELS_AC_TO_HOST)(long* Mesh_pointer,realw* h_rho_ac_kl,realw* h_kappa_ac_kl,int* NSPEC_AB) {}

void FC_FUNC_(transfer_kernels_hess_el_tohost,
              TRANSFER_KERNELS_HESS_EL_TOHOST)(long* Mesh_pointer,realw* h_hess_kl,int* NSPEC_AB) {}

void FC_FUNC_(transfer_kernels_hess_ac_tohost,
              TRANSFER_KERNELS_HESS_AC_TOHOST)(long* Mesh_pointer,realw* h_hess_ac_kl,int* NSPEC_AB) {}

//
// src/cuda/update_displacement_cuda.cu
//

void FC_FUNC_(update_displacement_cuda,
              UPDATE_DISPLACMENT_CUDA)(long* Mesh_pointer,
                                          realw* deltat_F,
                                          realw* deltatsqover2_F,
                                          realw* deltatover2_F,
                                          realw* b_deltat_F,
                                          realw* b_deltatsqover2_F,
                                          realw* b_deltatover2_F) {}

void FC_FUNC_(it_update_displacement_ac_cuda,
              it_update_displacement_ac_cuda)(long* Mesh_pointer,
                                               realw* deltat_F,
                                               realw* deltatsqover2_F,
                                               realw* deltatover2_F,
                                               realw* b_deltat_F,
                                               realw* b_deltatsqover2_F,
                                               realw* b_deltatover2_F) {}

void FC_FUNC_(kernel_3_a_cuda,
              KERNEL_3_A_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               realw* b_deltatover2_F) {}

void FC_FUNC_(kernel_3_b_cuda,
              KERNEL_3_B_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               realw* b_deltatover2_F) {}

void FC_FUNC_(kernel_3_a_acoustic_cuda,
              KERNEL_3_ACOUSTIC_CUDA)(long* Mesh_pointer ) {}

void FC_FUNC_(kernel_3_b_acoustic_cuda,
              KERNEL_3_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                      realw* deltatover2_F,
                                      realw* b_deltatover2_F) {}


//
// src/cuda/write_seismograms_cuda.cu
//

void FC_FUNC_(compute_seismograms_cuda,
              COMPUTE_SEISMOGRAMS_CUDA)(long* Mesh_pointer_f,int* seismotypef,double* sisux, double* sisuz,int* seismo_currentf,
                                   int* NSTEP_BETWEEN_OUTPUT_SEISMOSf,int * any_elastic_glob,int * any_acoustic_glob) {}



