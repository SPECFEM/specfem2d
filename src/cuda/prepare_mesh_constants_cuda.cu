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
#include "prepare_constants_cuda.h"


#ifdef USE_OLDER_CUDA4_GPU
#else
  #ifdef USE_TEXTURES_FIELDS
    // elastic
    extern realw_texture d_displ_tex;
    extern realw_texture d_accel_tex;
    // backward/reconstructed
    extern realw_texture d_b_displ_tex;
    extern realw_texture d_b_accel_tex;
    // acoustic
    extern realw_texture d_potential_tex;
    extern realw_texture d_potential_dot_dot_tex;
    // backward/reconstructed
    extern realw_texture d_b_potential_tex;
    extern realw_texture d_b_potential_dot_dot_tex;
  #endif
  #ifdef USE_TEXTURES_CONSTANTS
    extern realw_texture d_hprime_xx_tex;
    extern size_t d_hprime_xx_tex_offset;
    extern realw_texture d_wxgll_xx_tex;
    extern size_t d_wxgll_xx_tex_offset;
  #endif
#endif


/* ----------------------------------------------------------------------------------------------- */

// helper functions

/* ----------------------------------------------------------------------------------------------- */


// copies integer array from CPU host to GPU device
void copy_todevice_int(void** d_array_addr_ptr,int* h_array,int size){
  TRACE("  copy_todevice_int");

  // allocates memory on GPU
  //
  // note: cudaMalloc uses a double-pointer, such that it can return an error code in case it fails
  //          we thus pass the address to the pointer above (as void double-pointer) to have it
  //          pointing to the correct pointer of the array here
  print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(int)),
                          12001);

  // copies values onto GPU
  //
  // note: cudaMemcpy uses the pointer to the array, we thus re-cast the value of
  //          the double-pointer above to have the correct pointer to the array
  print_CUDA_error_if_any(cudaMemcpy((int*) *d_array_addr_ptr,h_array,size*sizeof(int),cudaMemcpyHostToDevice),
                          12002);
}

/* ----------------------------------------------------------------------------------------------- */

// copies integer array from CPU host to GPU device
void copy_todevice_realw(void** d_array_addr_ptr,realw* h_array,int size){
  TRACE("  copy_todevice_realw");

  // allocates memory on GPU
  print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(realw)),
                          22001);
  // copies values onto GPU
  print_CUDA_error_if_any(cudaMemcpy((realw*) *d_array_addr_ptr,h_array,size*sizeof(realw),cudaMemcpyHostToDevice),
                          22002);
}


/*
__global__ void check_field(int* ibool,int* nibool, int max_nibool,int num_interfaces_ext_mesh)
{


  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int ientry,iglob;

  for(int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    if (id<nibool[iinterface]) {

      // entry in interface array
      ientry = id + max_nibool*iinterface;
      // global index in wavefield
      iglob = ibool[ientry] - 1;

cuPrintf("valeurs de iglob %d, de l'indice d'entree %d,du nombre de points dans l'interface %d : %d\n", iglob,ientry,iinterface,nibool[iinterface]);
  }}
}*/
/* ----------------------------------------------------------------------------------------------- */

// GPU preparation

/* ----------------------------------------------------------------------------------------------- */

extern "C"
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
                                        int* nsources_local_f,
                                        realw* h_sourcearrays, realw * h_source_time_function,
                                        int* NSTEP,
                                        int* h_ispec_selected_source,
                                        int* h_ispec_selected_rec_loc,
                                        int* nrec_local,
                                        realw * h_cosrot,realw * h_sinrot,
                                        int* SIMULATION_TYPE,
                                        int* nspec_acoustic,int* nspec_elastic,
                                        int* h_myrank,
                                        int* SAVE_FORWARD,
                                        realw* h_xir_store, realw* h_gammar_store) {

  TRACE("prepare_constants_device");

  // allocates mesh parameter structure
  Mesh* mp = (Mesh*) malloc( sizeof(Mesh) );
  if (mp == NULL) exit_on_error("error allocating mesh pointer");
  *Mesh_pointer = (long)mp;

  if (*h_NGLLX != NGLLX) {
    exit_on_error("NGLLX defined in constants.h must match the NGLLX defined in src/cuda/mesh_constants_cuda.h");
  }

  // sets processes mpi rank
  mp->myrank = *h_myrank;

  // sets global parameters
  mp->NSPEC_AB = *NSPEC_AB;
  mp->NGLOB_AB = *NGLOB_AB;

  // constants
  mp->simulation_type = *SIMULATION_TYPE;
  mp->absorbing_conditions = *ABSORBING_CONDITIONS;
  mp->save_forward = *SAVE_FORWARD;

  // sets constant arrays
  setConst_hprime_xx(h_hprime_xx,mp);
  // setConst_hprime_zz(h_hprime_zz,mp); // only needed if NGLLX != NGLLY != NGLLZ

  setConst_hprimewgll_xx(h_hprimewgll_xx,mp);
  //setConst_hprimewgll_zz(h_hprimewgll_zz,mp); // only needed if NGLLX != NGLLY != NGLLZ

  setConst_wxgll(h_wxgll,mp);

  // Using texture memory for the hprime-style constants is slower on
  // Fermi generation hardware, but *may* be faster on Kepler
  // generation. We will reevaluate this again, so might as well leave
  // in the code with with #USE_TEXTURES_FIELDS not-defined.
  #ifdef USE_TEXTURES_CONSTANTS
  {
    // checks that realw is a float
    if (sizeof(realw) != sizeof(float)) exit_on_error("TEXTURES only work with realw selected as float");

    // note: device memory returned by cudaMalloc guarantees that the offset is 0,
    //       however here we use the global memory array d_hprime_xx and need to provide an offset variable for the function call
    // binds texture
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<realw>();

      const textureReference* d_hprime_xx_tex_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_hprime_xx_tex_ptr, "d_hprime_xx_tex"), 1101);
      print_CUDA_error_if_any(cudaBindTexture(&d_hprime_xx_tex_offset, d_hprime_xx_tex_ptr, mp->d_hprime_xx, &channelDesc, sizeof(realw)*(NGLL2)), 1102);

      const textureReference* d_wxgll_xx_tex_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_wxgll_xx_tex_ptr, "d_wxgll_xx_tex"), 1103);
      print_CUDA_error_if_any(cudaBindTexture(&d_wxgll_xx_tex_offset, d_wxgll_xx_tex_ptr, mp->d_wxgll, &channelDesc, sizeof(realw)*(NGLL2)), 1104);

   #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<realw>();

      print_CUDA_error_if_any(cudaBindTexture(&d_hprime_xx_tex_offset, &d_hprime_xx_tex, mp->d_hprime_xx,
                                              &channelDesc, sizeof(realw)*(NGLL2)), 1105);
      //printf("Bind texture hprime_xx offset = %lu\n",d_hprime_xx_tex_offset);

      print_CUDA_error_if_any(cudaBindTexture(&d_wxgll_xx_tex_offset, &d_wxgll_xx_tex, mp->d_wxgll,
                                              &channelDesc, sizeof(realw)*(NGLLX)), 1106);
      //printf("Bind texture wxgll_xx offset = %lu\n",d_wxgll_xx_tex_offset);
   #endif
  }
  #endif

  // mesh
  // Assuming NGLLX=5. Padded is then 32 (5^2+3)
  int size_padded = NGLL2_PADDED * (mp->NSPEC_AB);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xix, size_padded*sizeof(realw)),1001);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiz, size_padded*sizeof(realw)),1002);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammax, size_padded*sizeof(realw)),1003);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammaz, size_padded*sizeof(realw)),1004);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_kappav, size_padded*sizeof(realw)),1005);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_muv, size_padded*sizeof(realw)),1006);

  print_CUDA_error_if_any(cudaMemcpy2D(mp->d_xix, NGLL2_PADDED*sizeof(realw),
                                       h_xix, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                       mp->NSPEC_AB, cudaMemcpyHostToDevice),1501);
  print_CUDA_error_if_any(cudaMemcpy2D(mp->d_xiz, NGLL2_PADDED*sizeof(realw),
                                       h_xiz, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                       mp->NSPEC_AB, cudaMemcpyHostToDevice),1503);
  print_CUDA_error_if_any(cudaMemcpy2D(mp->d_gammax, NGLL2_PADDED*sizeof(realw),
                                       h_gammax, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                       mp->NSPEC_AB, cudaMemcpyHostToDevice),1507);
  print_CUDA_error_if_any(cudaMemcpy2D(mp->d_gammaz, NGLL2_PADDED*sizeof(realw),
                                       h_gammaz, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                       mp->NSPEC_AB, cudaMemcpyHostToDevice),1509);
  print_CUDA_error_if_any(cudaMemcpy2D(mp->d_kappav, NGLL2_PADDED*sizeof(realw),
                                       h_kappav, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                       mp->NSPEC_AB, cudaMemcpyHostToDevice),1510);
  print_CUDA_error_if_any(cudaMemcpy2D(mp->d_muv, NGLL2_PADDED*sizeof(realw),
                                       h_muv, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                       mp->NSPEC_AB, cudaMemcpyHostToDevice),1511);

  // global indexing (padded)
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibool, size_padded*sizeof(int)),1600);
  print_CUDA_error_if_any(cudaMemcpy2D(mp->d_ibool, NGLL2_PADDED*sizeof(int),
                                       h_ibool, NGLL2*sizeof(int), NGLL2*sizeof(int),
                                       mp->NSPEC_AB, cudaMemcpyHostToDevice),1601);

  // prepare interprocess-edge exchange information
  mp->num_interfaces_ext_mesh = *num_interfaces_ext_mesh;
  mp->max_nibool_interfaces_ext_mesh = *max_nibool_interfaces_ext_mesh;
  if (mp->num_interfaces_ext_mesh > 0) {
    copy_todevice_int((void**)&mp->d_nibool_interfaces_ext_mesh,h_nibool_interfaces_ext_mesh,
                      mp->num_interfaces_ext_mesh);
    copy_todevice_int((void**)&mp->d_ibool_interfaces_ext_mesh,h_ibool_interfaces_ext_mesh,
                      (mp->num_interfaces_ext_mesh)*(mp->max_nibool_interfaces_ext_mesh));
    //int blocksize = BLOCKSIZE_TRANSFER;
    //int size_padded = ((int)ceil(((double)(mp->max_nibool_interfaces_ext_mesh))/((double)blocksize)))*blocksize;
  }
  mp->size_mpi_buffer = 0;
  mp->size_mpi_buffer_potential = 0;

  // streams
  cudaStreamCreate(&mp->compute_stream);
  // copy stream (needed to transfer mpi buffers)
  if (mp->num_interfaces_ext_mesh * mp->max_nibool_interfaces_ext_mesh > 0) {
  cudaStreamCreate(&mp->copy_stream);
  }

  // inner elements
  copy_todevice_int((void**)&mp->d_ispec_is_inner,h_ispec_is_inner,mp->NSPEC_AB);

  // absorbing boundaries
  mp->d_num_abs_boundary_faces = *h_num_abs_boundary_faces;
  if (mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0) {
    copy_todevice_int((void**)&mp->d_abs_boundary_ispec,h_abs_boundary_ispec,mp->d_num_abs_boundary_faces);
    copy_todevice_int((void**)&mp->d_abs_boundary_ijk,h_abs_boundary_ij,
                      2*NGLLX*(mp->d_num_abs_boundary_faces));
    copy_todevice_realw((void**)&mp->d_abs_boundary_normal,h_abs_boundary_normal,
                        NDIM*NGLLX*(mp->d_num_abs_boundary_faces));
    copy_todevice_realw((void**)&mp->d_abs_boundary_jacobian2Dw,h_abs_boundary_jacobian1Dw,
                        NGLLX*(mp->d_num_abs_boundary_faces));
    copy_todevice_int((void**)&mp->d_cote_abs,h_cote_abs,(mp->d_num_abs_boundary_faces));
    copy_todevice_int((void**)&mp->d_ib_left,h_ib_left,(mp->d_num_abs_boundary_faces));
    copy_todevice_int((void**)&mp->d_ib_right,h_ib_right,(mp->d_num_abs_boundary_faces));
    copy_todevice_int((void**)&mp->d_ib_top,h_ib_top,(mp->d_num_abs_boundary_faces));
    copy_todevice_int((void**)&mp->d_ib_bottom,h_ib_bottom,(mp->d_num_abs_boundary_faces));

    mp->d_nspec_bottom = *h_nspec_bottom;
    mp->d_nspec_left = *h_nspec_left;
    mp->d_nspec_right = *h_nspec_right;
    mp->d_nspec_top = *h_nspec_top;
  }

  // sources
  mp->nsources_local = *nsources_local_f;

  if (mp->nsources_local > 0){
    copy_todevice_realw((void**)&mp->d_source_time_function,h_source_time_function,(*NSTEP)*(mp->nsources_local));
    copy_todevice_realw((void**)&mp->d_sourcearrays,h_sourcearrays,mp->nsources_local*NDIM*NGLL2);
    copy_todevice_int((void**)&mp->d_ispec_selected_source,h_ispec_selected_source,mp->nsources_local);
  }

  // receiver stations
  mp->nrec_local = *nrec_local; // number of receiver located in this partition
  // note that: size of size(ispec_selected_rec_loc) = nrec_local
  if (mp->nrec_local > 0) {
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_seismograms,2*(*NSTEP)*sizeof(realw)*(mp->nrec_local)*2),1303);
    // pinned memory
    print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_seismograms),2*(*NSTEP)*sizeof(realw)*(mp->nrec_local)*2),8004);
    // host memory
    //mp->h_seismograms = (float*)malloc((mp->nrec_local)*2*sizeof(float));
    //if (mp->h_seismograms == NULL) exit_on_error("h_seismograms not allocated \n");

    copy_todevice_realw((void**)&mp->d_cosrot,h_cosrot,mp->nrec_local);
    copy_todevice_realw((void**)&mp->d_sinrot,h_sinrot,mp->nrec_local);

    copy_todevice_realw((void**)&mp->d_xir_store_loc,h_xir_store,(mp->nrec_local)*NGLLX);
    copy_todevice_realw((void**)&mp->d_gammar_store_loc,h_gammar_store,(mp->nrec_local)*NGLLX);

    copy_todevice_int((void**)&mp->d_ispec_selected_rec_loc,h_ispec_selected_rec_loc,mp->nrec_local);
  }

  // number of elements per domain
  mp->nspec_acoustic = *nspec_acoustic;
  mp->nspec_elastic  = *nspec_elastic;

  // JC JC here we will need to add GPU support for the new C-PML routines

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_constants_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// for ACOUSTIC simulations

/* ----------------------------------------------------------------------------------------------- */


extern "C"
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
                                              int * h_ninterface_acoustic,int * h_inum_interfaces_acoustic,int* ATTENUATION_VISCOACOUSTIC,
                                              realw* h_A_newmark,realw* h_B_newmark,int* NO_BACKWARD_RECONSTRUCTION,realw* h_no_backward_acoustic_buffer) {

  TRACE("prepare_fields_acoustic_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // allocates arrays on device (GPU)
  int size = mp->NGLOB_AB;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_potential_acoustic),sizeof(realw)*size),2001);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_potential_dot_acoustic),sizeof(realw)*size),2002);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_potential_dot_dot_acoustic),sizeof(realw)*size),2003);
  // initializes values to zero
  print_CUDA_error_if_any(cudaMemset(mp->d_potential_acoustic,0,sizeof(realw)*size),2007);
  print_CUDA_error_if_any(cudaMemset(mp->d_potential_dot_acoustic,0,sizeof(realw)*size),2007);
  print_CUDA_error_if_any(cudaMemset(mp->d_potential_dot_dot_acoustic,0,sizeof(realw)*size),2007);

  #ifdef USE_TEXTURES_FIELDS
  {
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

  // mpi buffer
  mp->size_mpi_buffer_potential = (mp->num_interfaces_ext_mesh) * (mp->max_nibool_interfaces_ext_mesh);
  if (mp->size_mpi_buffer_potential > 0) {
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_send_potential_dot_dot_buffer),mp->size_mpi_buffer_potential *sizeof(realw)),2004);
  }

  // mass matrix
  copy_todevice_realw((void**)&mp->d_rmass_acoustic,rmass_acoustic,mp->NGLOB_AB);

  // density
  // padded array
  // Assuming NGLLX==5. Padded is then 32 (5^2+3)
  int size_padded = NGLL2_PADDED * mp->NSPEC_AB;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rhostore),size_padded*sizeof(realw)),2006);
  // transfer constant element data with padding
  print_CUDA_error_if_any(cudaMemcpy2D(mp->d_rhostore, NGLL2_PADDED*sizeof(realw),
                                       rhostore, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                       mp->NSPEC_AB, cudaMemcpyHostToDevice),2106);

  // non-padded array
  copy_todevice_realw((void**)&mp->d_kappastore,kappastore,NGLL2*mp->NSPEC_AB);

  // phase elements
  mp->num_phase_ispec_acoustic = *num_phase_ispec_acoustic;
  copy_todevice_int((void**)&mp->d_phase_ispec_inner_acoustic,phase_ispec_inner_acoustic,
                    2*mp->num_phase_ispec_acoustic);
  copy_todevice_int((void**)&mp->d_ispec_is_acoustic,ispec_is_acoustic,mp->NSPEC_AB);

  // allocate surface arrays
  mp->num_free_surface_faces = *num_free_surface_faces;
  if (mp->num_free_surface_faces > 0) {
    copy_todevice_int((void**)&mp->d_free_surface_ispec,free_surface_ispec,mp->num_free_surface_faces);
    copy_todevice_int((void**)&mp->d_free_surface_ijk,free_surface_ijk,2*NGLLX*mp->num_free_surface_faces);
  }

  // absorbing boundaries
  if (mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0) {
    // absorb_field array used for file i/o
    if (mp->simulation_type == 3 || ( mp->simulation_type == 1 && mp->save_forward )){
      print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_b_absorb_potential_left,mp->d_nspec_left*sizeof(realw)*NGLLX),2201);
      print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_b_absorb_potential_right,mp->d_nspec_right*sizeof(realw)*NGLLX),2201);
      print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_b_absorb_potential_top,mp->d_nspec_top*sizeof(realw)*NGLLX),2201);
      print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_b_absorb_potential_bottom,mp->d_nspec_bottom*sizeof(realw)*NGLLX),2201);

    }
  }

  // coupling with elastic parts
  if (*ELASTIC_SIMULATION && *num_coupling_ac_el_faces > 0) {
    copy_todevice_int((void**)&mp->d_coupling_ac_el_ispec,coupling_ac_el_ispec,(*num_coupling_ac_el_faces));
    copy_todevice_int((void**)&mp->d_coupling_ac_el_ijk,coupling_ac_el_ijk,2*NGLLX*(*num_coupling_ac_el_faces));
    copy_todevice_realw((void**)&mp->d_coupling_ac_el_normal,coupling_ac_el_normal,
                        2*NGLLX*(*num_coupling_ac_el_faces));
    copy_todevice_realw((void**)&mp->d_coupling_ac_el_jacobian2Dw,coupling_ac_el_jacobian2Dw,
                        NGLLX*(*num_coupling_ac_el_faces));
  }

  mp->ninterface_acoustic = *h_ninterface_acoustic;
  copy_todevice_int((void**)&mp->d_inum_interfaces_acoustic,h_inum_interfaces_acoustic,mp->num_interfaces_ext_mesh);

  if (*ATTENUATION_VISCOACOUSTIC) {
    copy_todevice_realw((void**)&mp->d_A_newmark_acous,h_A_newmark,NGLL2*mp->NSPEC_AB*N_SLS);
    copy_todevice_realw((void**)&mp->d_B_newmark_acous,h_B_newmark,NGLL2*mp->NSPEC_AB*N_SLS);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_e1_acous,mp->NSPEC_AB*sizeof(realw)*NGLL2*N_SLS),2202);
    print_CUDA_error_if_any(cudaMemset(mp->d_e1_acous,0,mp->NSPEC_AB*sizeof(realw)*NGLL2*N_SLS),2203);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_sum_forces_old,mp->NSPEC_AB*sizeof(realw)*NGLL2),2204);
    print_CUDA_error_if_any(cudaMemset(mp->d_sum_forces_old,0,mp->NSPEC_AB*sizeof(realw)*NGLL2),2205);
  }

  if (*NO_BACKWARD_RECONSTRUCTION){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_potential_acoustic_buffer),mp->NGLOB_AB*sizeof(realw)),2206);
    cudaStreamCreateWithFlags(&mp->copy_stream_no_backward,cudaStreamNonBlocking);
    cudaHostRegister(h_no_backward_acoustic_buffer,3*mp->NGLOB_AB*sizeof(realw),0);
    cudaEventCreate(&mp->transfer_is_complete1);
    cudaEventCreate(&mp->transfer_is_complete2);
  }


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_fields_acoustic_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_acoustic_adj_dev,
              PREPARE_FIELDS_ACOUSTIC_ADJ_DEV)(long* Mesh_pointer,
                                               int* APPROXIMATE_HESS_KL,
                                               int* ATTENUATION_VISCOACOUSTIC,
                                               int* NO_BACKWARD_RECONSTRUCTION) {

  TRACE("prepare_fields_acoustic_adj_dev");

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // kernel simulations
  if (mp->simulation_type != 3 ) return;

  // allocates backward/reconstructed arrays on device (GPU)
  int size = mp->NGLOB_AB * sizeof(realw);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_potential_acoustic),size),3014);
  if (! *NO_BACKWARD_RECONSTRUCTION){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_potential_dot_acoustic),size),3015);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_potential_dot_dot_acoustic),size),3016);
  }
  // initializes values to zero
  //print_CUDA_error_if_any(cudaMemset(mp->d_b_potential_acoustic,0,size),3007);
  //print_CUDA_error_if_any(cudaMemset(mp->d_b_potential_dot_acoustic,0,size),3007);
  //print_CUDA_error_if_any(cudaMemset(mp->d_b_potential_dot_dot_acoustic,0,size),3007);

  #ifdef USE_TEXTURES_FIELDS
  {
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      const textureReference* d_b_potential_tex_ref_ptr;

      print_CUDA_error_if_any(cudaGetTextureReference(&d_b_potential_tex_ref_ptr, "d_b_potential_tex"), 3001);
      print_CUDA_error_if_any(cudaBindTexture(0, d_b_potential_tex_ref_ptr, mp->d_b_potential_acoustic, &channelDesc, size), 3001);

      if (! *NO_BACKWARD_RECONSTRUCTION){
        const textureReference* d_b_potential_dot_dot_tex_ref_ptr;
        print_CUDA_error_if_any(cudaGetTextureReference(&d_b_potential_dot_dot_tex_ref_ptr, "d_b_potential_dot_dot_tex"),3003);
        print_CUDA_error_if_any(cudaBindTexture(0, d_b_potential_dot_dot_tex_ref_ptr, mp->d_b_potential_dot_dot_acoustic, &channelDesc, size), 3003);
      }
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_b_potential_tex, mp->d_b_potential_acoustic, &channelDesc, size), 3001);
      if (! *NO_BACKWARD_RECONSTRUCTION) print_CUDA_error_if_any(cudaBindTexture(0, &d_b_potential_dot_dot_tex, mp->d_b_potential_dot_dot_acoustic, &channelDesc, size), 3003);
    #endif
  }
  #endif

  // allocates kernels
  size = NGLL2 * mp->NSPEC_AB * sizeof(realw);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_ac_kl),size),3017);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_kappa_ac_kl),size),3018);
  // initializes kernel values to zero
  print_CUDA_error_if_any(cudaMemset(mp->d_rho_ac_kl,0,size),3019);
  print_CUDA_error_if_any(cudaMemset(mp->d_kappa_ac_kl,0,size),3020);

  // preconditioner
  if (*APPROXIMATE_HESS_KL) {
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_hess_ac_kl),size),3030);
    // initializes with zeros
    print_CUDA_error_if_any(cudaMemset(mp->d_hess_ac_kl,0,size),3031);
  }

  if (*ATTENUATION_VISCOACOUSTIC && (! *NO_BACKWARD_RECONSTRUCTION) ) {
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_sum_forces_old),size),3040);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_e1_acous),size*N_SLS),3041);
  }

  // mpi buffer
  if (mp->size_mpi_buffer_potential > 0 && (! *NO_BACKWARD_RECONSTRUCTION)) {
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_send_potential_dot_dot_buffer),mp->size_mpi_buffer_potential*sizeof(realw)),3014);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_fields_acoustic_adj_dev");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// for ELASTIC simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
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
                                             int* ANISOTROPY,
                                             realw *c11store,realw *c12store,realw *c13store,
                                             realw *c15store,
                                             realw *c23store,
                                             realw *c25store,realw *c33store,
                                             realw *c35store,
                                             realw *c55store,int* h_ninterface_elastic,int * h_inum_interfaces_elastic,int* ATTENUATION_VISCOELASTIC,
                                             realw* h_A_newmark_mu,realw* h_B_newmark_mu,realw* h_A_newmark_kappa,realw* h_B_newmark_kappa) {



  TRACE("prepare_fields_elastic_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  int size,size_padded;

  // debug
  //printf("prepare_fields_elastic_device: rank %d - wavefield setup\n",mp->myrank);
  //synchronize_mpi();

  // elastic wavefields
  size = NDIM * mp->NGLOB_AB;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_displ),sizeof(realw)*size),4001);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_veloc),sizeof(realw)*size),4002);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_accel),sizeof(realw)*size),4003);
  // initializes values to zero
  //print_CUDA_error_if_any(cudaMemset(mp->d_displ,0,sizeof(realw)*size),4007);
  //print_CUDA_error_if_any(cudaMemset(mp->d_veloc,0,sizeof(realw)*size),4007);
  //print_CUDA_error_if_any(cudaMemset(mp->d_accel,0,sizeof(realw)*size),4007);

 #ifdef USE_TEXTURES_FIELDS
  {
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      const textureReference* d_displ_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_displ_tex_ref_ptr, "d_displ_tex"), 4004);
      print_CUDA_error_if_any(cudaBindTexture(0, d_displ_tex_ref_ptr, mp->d_displ, &channelDesc, sizeof(realw)*size), 4005);
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_displ_tex, mp->d_displ, &channelDesc, sizeof(realw)*size), 4008);
    #endif
  }
  #endif


  // debug
  //synchronize_mpi();

  // MPI buffer
  mp->size_mpi_buffer = NDIM * (mp->num_interfaces_ext_mesh) * (mp->max_nibool_interfaces_ext_mesh);
  if (mp->size_mpi_buffer > 0) {
    // note: Allocate pinned mpi-buffers.
    //       MPI buffers use pinned memory allocated by cudaMallocHost, which
    //       enables the use of asynchronous memory copies from host <-> device
    // send buffer
    print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_send_accel_buffer),sizeof(realw)*(mp->size_mpi_buffer)),8004);
    //mp->send_buffer = (float*)malloc((mp->size_mpi_buffer)*sizeof(float));
    // adjoint
    //print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_send_b_accel_buffer),sizeof(float)*(mp->size_mpi_buffer)),8004);
    // mp->b_send_buffer = (float*)malloc((size_mpi_buffer)*sizeof(float));
    // receive buffer
    print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_recv_accel_buffer),sizeof(realw)*(mp->size_mpi_buffer)),8004);
    //mp->recv_buffer = (float*)malloc((mp->size_mpi_buffer)*sizeof(float));

    // non-pinned buffer
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_send_accel_buffer),mp->size_mpi_buffer*sizeof(realw)),4004);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_recv_accel_buffer),mp->size_mpi_buffer*sizeof(realw)),4004);

    // adjoint
    if (mp->simulation_type == 3) {
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_send_accel_buffer),mp->size_mpi_buffer*sizeof(realw)),4004);
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_recv_accel_buffer),mp->size_mpi_buffer*sizeof(realw)),4004);
    }
  }

  // debug
  //printf("prepare_fields_elastic_device: rank %d - mass matrix\n",mp->myrank);
  //synchronize_mpi();

  // mass matrix
  copy_todevice_realw((void**)&mp->d_rmassx,rmassx,mp->NGLOB_AB);
  copy_todevice_realw((void**)&mp->d_rmassz,rmassz,mp->NGLOB_AB);

  // element indices
  copy_todevice_int((void**)&mp->d_ispec_is_elastic,ispec_is_elastic,mp->NSPEC_AB);

  // phase elements
  mp->num_phase_ispec_elastic = *num_phase_ispec_elastic;

  copy_todevice_int((void**)&mp->d_phase_ispec_inner_elastic,phase_ispec_inner_elastic,2*mp->num_phase_ispec_elastic);

  // debug
  //synchronize_mpi();

  // absorbing conditions
  if (mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0){

    // debug
    //printf("prepare_fields_elastic_device: rank %d - absorbing boundary setup\n",mp->myrank);

    // non-padded arrays
    // rho_vp, rho_vs non-padded; they are needed for stacey boundary condition
    copy_todevice_realw((void**)&mp->d_rho_vp,rho_vp,NGLL2*mp->NSPEC_AB);
    copy_todevice_realw((void**)&mp->d_rho_vs,rho_vs,NGLL2*mp->NSPEC_AB);

    // absorb_field array used for file i/o
    if (mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0) {
      // absorb_field array used for file i/o
      if (mp->simulation_type == 3 || ( mp->simulation_type == 1 && mp->save_forward )){
        mp->d_nspec_left = *h_nspec_left;
        print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_b_absorb_elastic_left,2*mp->d_nspec_left*sizeof(realw)*NGLLX),2201);

        mp->d_nspec_right = *h_nspec_right;
        print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_b_absorb_elastic_right,2*mp->d_nspec_right*sizeof(realw)*NGLLX),2201);

        mp->d_nspec_top = *h_nspec_top;
        print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_b_absorb_elastic_top,2*mp->d_nspec_top*sizeof(realw)*NGLLX),2201);

        mp->d_nspec_bottom = *h_nspec_bottom;
        print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_b_absorb_elastic_bottom,2*mp->d_nspec_bottom*sizeof(realw)*NGLLX),2201);
      }
    }
  }

  // debug
  //synchronize_mpi();

  // anisotropy
  if (*ANISOTROPY) {
    // debug
    //printf("prepare_fields_elastic_device: rank %d - attenuation setup\n",mp->myrank);
    //synchronize_mpi();

    // Assuming NGLLX==5. Padded is then 32 (5^2+3)
    size_padded = NGLL2_PADDED * (mp->NSPEC_AB);

    // allocates memory on GPU
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c11store),size_padded*sizeof(realw)),4700);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c12store),size_padded*sizeof(realw)),4701);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c13store),size_padded*sizeof(realw)),4702);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c15store),size_padded*sizeof(realw)),4704);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c23store),size_padded*sizeof(realw)),4707);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c25store),size_padded*sizeof(realw)),4709);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c33store),size_padded*sizeof(realw)),4711);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c35store),size_padded*sizeof(realw)),4711);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c55store),size_padded*sizeof(realw)),4718);

    print_CUDA_error_if_any(cudaMemcpy2D(mp->d_c11store, NGLL2_PADDED*sizeof(realw),
                                         c11store, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                         mp->NSPEC_AB, cudaMemcpyHostToDevice),4800);
    print_CUDA_error_if_any(cudaMemcpy2D(mp->d_c12store, NGLL2_PADDED*sizeof(realw),
                                         c12store, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                         mp->NSPEC_AB, cudaMemcpyHostToDevice),4800);
    print_CUDA_error_if_any(cudaMemcpy2D(mp->d_c13store, NGLL2_PADDED*sizeof(realw),
                                         c13store, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                         mp->NSPEC_AB, cudaMemcpyHostToDevice),4800);
    print_CUDA_error_if_any(cudaMemcpy2D(mp->d_c15store, NGLL2_PADDED*sizeof(realw),
                                         c15store, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                         mp->NSPEC_AB, cudaMemcpyHostToDevice),4800);
    print_CUDA_error_if_any(cudaMemcpy2D(mp->d_c23store, NGLL2_PADDED*sizeof(realw),
                                         c23store, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                         mp->NSPEC_AB, cudaMemcpyHostToDevice),4800);
    print_CUDA_error_if_any(cudaMemcpy2D(mp->d_c25store, NGLL2_PADDED*sizeof(realw),
                                         c25store, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                         mp->NSPEC_AB, cudaMemcpyHostToDevice),4800);
    print_CUDA_error_if_any(cudaMemcpy2D(mp->d_c33store, NGLL2_PADDED*sizeof(realw),
                                         c33store, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                         mp->NSPEC_AB, cudaMemcpyHostToDevice),4800);
    print_CUDA_error_if_any(cudaMemcpy2D(mp->d_c35store, NGLL2_PADDED*sizeof(realw),
                                         c35store, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                         mp->NSPEC_AB, cudaMemcpyHostToDevice),4800);
    print_CUDA_error_if_any(cudaMemcpy2D(mp->d_c55store, NGLL2_PADDED*sizeof(realw),
                                         c55store, NGLL2*sizeof(realw), NGLL2*sizeof(realw),
                                         mp->NSPEC_AB, cudaMemcpyHostToDevice),4800);
  }

  mp->ninterface_elastic = *h_ninterface_elastic;
  copy_todevice_int((void**)&mp->d_inum_interfaces_elastic,h_inum_interfaces_elastic,mp->num_interfaces_ext_mesh);


  if (*ATTENUATION_VISCOELASTIC) {
    copy_todevice_realw((void**)&mp->d_A_newmark_mu,h_A_newmark_mu,NGLL2*mp->NSPEC_AB*N_SLS);
    copy_todevice_realw((void**)&mp->d_B_newmark_mu,h_B_newmark_mu,NGLL2*mp->NSPEC_AB*N_SLS);
    copy_todevice_realw((void**)&mp->d_A_newmark_kappa,h_A_newmark_kappa,NGLL2*mp->NSPEC_AB*N_SLS);
    copy_todevice_realw((void**)&mp->d_B_newmark_kappa,h_B_newmark_kappa,NGLL2*mp->NSPEC_AB*N_SLS);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_e1,mp->NSPEC_AB*sizeof(realw)*NGLL2*N_SLS),4801);
    print_CUDA_error_if_any(cudaMemset(mp->d_e1,0,mp->NSPEC_AB*sizeof(realw)*NGLL2*N_SLS),4802);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_e11,mp->NSPEC_AB*sizeof(realw)*NGLL2*N_SLS),4803);
    print_CUDA_error_if_any(cudaMemset(mp->d_e11,0,mp->NSPEC_AB*sizeof(realw)*NGLL2*N_SLS),4804);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_e13,mp->NSPEC_AB*sizeof(realw)*NGLL2*N_SLS),4805);
    print_CUDA_error_if_any(cudaMemset(mp->d_e13,0,mp->NSPEC_AB*sizeof(realw)*NGLL2*N_SLS),4806);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_dux_dxl_old,mp->NSPEC_AB*sizeof(realw)*NGLL2),4807);
    print_CUDA_error_if_any(cudaMemset(mp->d_dux_dxl_old,0,mp->NSPEC_AB*sizeof(realw)*NGLL2),4808);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_duz_dzl_old,mp->NSPEC_AB*sizeof(realw)*NGLL2),4809);
    print_CUDA_error_if_any(cudaMemset(mp->d_duz_dzl_old,0,mp->NSPEC_AB*sizeof(realw)*NGLL2),4810);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_dux_dzl_plus_duz_dxl_old,mp->NSPEC_AB*sizeof(realw)*NGLL2),4811);
    print_CUDA_error_if_any(cudaMemset(mp->d_dux_dzl_plus_duz_dxl_old,0,mp->NSPEC_AB*sizeof(realw)*NGLL2),4812);
  }



  // JC JC here we will need to add GPU support for the new C-PML routines

  // debug
  //printf("prepare_fields_elastic_device: rank %d - done\n",mp->myrank);
  //synchronize_mpi();

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_fields_elastic_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_elastic_adj_dev,
              PREPARE_FIELDS_ELASTIC_ADJ_DEV)(long* Mesh_pointer,
                                              int* size_f,
                                              int* APPROXIMATE_HESS_KL){

  TRACE("prepare_fields_elastic_adj_dev");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  int size;

  // checks if kernel simulation
  if (mp->simulation_type != 3 ) return;

  // kernel simulations
  // debug
  //printf("prepare_fields_elastic_adj_dev: rank %d - kernel setup\n",mp->myrank);
  //synchronize_mpi();

  // backward/reconstructed wavefields
  // allocates backward/reconstructed arrays on device (GPU)
  size = *size_f;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_displ),sizeof(realw)*size),5201);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_veloc),sizeof(realw)*size),5202);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_accel),sizeof(realw)*size),5203);
  // initializes values to zero
  //print_CUDA_error_if_any(cudaMemset(mp->d_b_displ,0,sizeof(realw)*size),5207);
  //print_CUDA_error_if_any(cudaMemset(mp->d_b_veloc,0,sizeof(realw)*size),5207);
  //print_CUDA_error_if_any(cudaMemset(mp->d_b_accel,0,sizeof(realw)*size),5207);

 #ifdef USE_TEXTURES_FIELDS
  {
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      const textureReference* d_b_displ_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_b_displ_tex_ref_ptr, "d_b_displ_tex"), 5204);
      print_CUDA_error_if_any(cudaBindTexture(0, d_b_displ_tex_ref_ptr, mp->d_b_displ, &channelDesc, sizeof(realw)*size), 5205);
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_b_displ_tex, mp->d_b_displ, &channelDesc, sizeof(realw)*size), 5208);
    #endif
  }
  #endif

  // anisotropic/isotropic kernels
  // debug
  //printf("prepare_fields_elastic_adj_dev: rank %d -  anisotropic/isotropic kernels\n",mp->myrank);
  //synchronize_mpi();

  // allocates kernels
  size = NGLL2 * mp->NSPEC_AB; // note: non-aligned; if align, check memcpy below and indexing
  // density kernel
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_kl),size*sizeof(realw)),5211);
  // initializes kernel values to zero
  print_CUDA_error_if_any(cudaMemset(mp->d_rho_kl,0,size*sizeof(realw)),5212);


    // isotropic kernels
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_mu_kl),size*sizeof(realw)),5213);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_kappa_kl),size*sizeof(realw)),5214);
    print_CUDA_error_if_any(cudaMemset(mp->d_mu_kl,0,size*sizeof(realw)),5216);
    print_CUDA_error_if_any(cudaMemset(mp->d_kappa_kl,0,size*sizeof(realw)),5217);

    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_dsxx),size*sizeof(realw)),5218);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_dsxz),size*sizeof(realw)),5219);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_dszz),size*sizeof(realw)),5220);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_dsxx),size*sizeof(realw)),5221);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_dsxz),size*sizeof(realw)),5222);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_dszz),size*sizeof(realw)),5223);

  // approximate hessian kernel
  if (*APPROXIMATE_HESS_KL) {
    // debug
    //printf("prepare_fields_elastic_adj_dev: rank %d - hessian kernel\n",mp->myrank);
    //synchronize_mpi();

    size = NGLL2 * mp->NSPEC_AB; // note: non-aligned; if align, check memcpy below and indexing
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_hess_el_kl),size*sizeof(realw)),5450);
    print_CUDA_error_if_any(cudaMemset(mp->d_hess_el_kl,0,size*sizeof(realw)),5451);
  }

  // debug
  //printf("prepare_fields_elastic_adj_dev: rank %d - done\n",mp->myrank);
  //synchronize_mpi();

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_fields_elastic_adj_dev");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// purely adjoint & kernel simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_sim2_or_3_const_device,
              PREPARE_SIM2_OR_3_CONST_DEVICE)(long* Mesh_pointer,
                                              int* nadj_rec_local,
                                              realw* h_source_adjointe,
                                              int* NSTEP) {

  TRACE("prepare_sim2_or_3_const_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // adjoint source arrays
  mp->nadj_rec_local = *nadj_rec_local;
  if (mp->nadj_rec_local > 0) {
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_adj_sourcearrays,(mp->nadj_rec_local)*2*NGLL2*sizeof(realw)),6003);

    copy_todevice_realw((void**)&mp->d_source_adjointe,h_source_adjointe,(*NSTEP)*(*nadj_rec_local)*NDIM);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_sim2_or_3_const_device");
#endif
}



/* ----------------------------------------------------------------------------------------------- */

// cleanup

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_cleanup_device,
              PREPARE_CLEANUP_DEVICE)(long* Mesh_pointer,
                                      int* ACOUSTIC_SIMULATION,
                                      int* ELASTIC_SIMULATION,
                                      int* ABSORBING_CONDITIONS,
                                      int* ANISOTROPY,
                                      int* APPROXIMATE_HESS_KL,
                                      int* ATTENUATION_VISCOACOUSTIC,
                                      int* ATTENUATION_VISCOELASTIC,
                                      int* NO_BACKWARD_RECONSTRUCTION,
                                      realw * h_no_backward_acoustic_buffer) {


TRACE("prepare_cleanup_device");

  // frees allocated memory arrays
  Mesh* mp = (Mesh*)(*Mesh_pointer);

#ifdef USE_TEXTURES_CONSTANTS
  cudaUnbindTexture(d_hprime_xx_tex);
  cudaUnbindTexture(d_wxgll_xx_tex);
#endif

  // frees memory on GPU
  // mesh
  cudaFree(mp->d_xix);
  cudaFree(mp->d_xiz);
  cudaFree(mp->d_gammax);
  cudaFree(mp->d_gammaz);
  cudaFree(mp->d_muv);
  cudaFree(mp->d_kappav);

  // absorbing boundaries
  if (*ABSORBING_CONDITIONS && mp->d_num_abs_boundary_faces > 0) {
    cudaFree(mp->d_abs_boundary_ispec);
    cudaFree(mp->d_abs_boundary_ijk);
    cudaFree(mp->d_abs_boundary_normal);
    cudaFree(mp->d_abs_boundary_jacobian2Dw);
    cudaFree(mp->d_cote_abs);
    cudaFree(mp->d_ib_left);
    cudaFree(mp->d_ib_right);
    cudaFree(mp->d_ib_top);
    cudaFree(mp->d_ib_bottom);
  }

  // interfaces
  if (mp->num_interfaces_ext_mesh > 0) {
    cudaFree(mp->d_nibool_interfaces_ext_mesh);
    cudaFree(mp->d_ibool_interfaces_ext_mesh);
  }

  // global indexing
  cudaFree(mp->d_ispec_is_inner);
  cudaFree(mp->d_ibool);

  // sources
  if (mp->nsources_local > 0){
    cudaFree(mp->d_sourcearrays);
    cudaFree(mp->d_source_time_function);
    cudaFree(mp->d_ispec_selected_source);
  }

  // receivers
  if (mp->nrec_local > 0) {
    cudaFree(mp->d_seismograms);
    cudaFree(mp->d_cosrot),cudaFree(mp->d_sinrot);
    cudaFree(mp->d_gammar_store_loc);
    cudaFree(mp->d_xir_store_loc);
    cudaFree(mp->d_ispec_selected_rec_loc);
    cudaFreeHost(mp->h_seismograms);
  }

  // ACOUSTIC arrays
  if (*ACOUSTIC_SIMULATION) {
    cudaFree(mp->d_potential_acoustic);
    cudaFree(mp->d_potential_dot_acoustic);
    cudaFree(mp->d_potential_dot_dot_acoustic);
    if (mp->size_mpi_buffer_potential > 0 ) cudaFree(mp->d_send_potential_dot_dot_buffer);
    cudaFree(mp->d_rmass_acoustic);
    cudaFree(mp->d_rhostore);
    cudaFree(mp->d_kappastore);
    cudaFree(mp->d_phase_ispec_inner_acoustic);
    cudaFree(mp->d_ispec_is_acoustic);
    cudaFree(mp->d_inum_interfaces_acoustic);

    if (*NO_BACKWARD_RECONSTRUCTION){
      cudaFree(mp->d_potential_acoustic_buffer);
      cudaHostUnregister(h_no_backward_acoustic_buffer);
      cudaEventDestroy(mp->transfer_is_complete1);
      cudaEventDestroy(mp->transfer_is_complete2);

    }
    if (mp->simulation_type == 3) {
      cudaFree(mp->d_b_potential_acoustic);
      if (! *NO_BACKWARD_RECONSTRUCTION){
        cudaFree(mp->d_b_potential_dot_acoustic);
        cudaFree(mp->d_b_potential_dot_dot_acoustic);
      }
      cudaFree(mp->d_rho_ac_kl);
      cudaFree(mp->d_kappa_ac_kl);
      if (*APPROXIMATE_HESS_KL) cudaFree(mp->d_hess_ac_kl);
      if (mp->size_mpi_buffer_potential > 0 && ! *NO_BACKWARD_RECONSTRUCTION) cudaFree(mp->d_b_send_potential_dot_dot_buffer);
      if (*ATTENUATION_VISCOACOUSTIC && ! *NO_BACKWARD_RECONSTRUCTION) {
        cudaFree(mp->d_b_sum_forces_old);
        cudaFree(mp->d_b_e1_acous);
      }
    }

    if (*ABSORBING_CONDITIONS && mp->d_num_abs_boundary_faces > 0){
      if (mp->simulation_type == 3 || ( mp->simulation_type == 1 && mp->save_forward )){
        cudaFree(mp->d_b_absorb_potential_bottom);
        cudaFree(mp->d_b_absorb_potential_left);
        cudaFree(mp->d_b_absorb_potential_right);
        cudaFree(mp->d_b_absorb_potential_top);
       }
    }
    if (*ATTENUATION_VISCOACOUSTIC){
      cudaFree(mp->d_e1_acous);
      cudaFree(mp->d_A_newmark_acous);
      cudaFree(mp->d_B_newmark_acous);
      cudaFree(mp->d_sum_forces_old);
    }

  } // ACOUSTIC_SIMULATION

  // ELASTIC arrays
  if (*ELASTIC_SIMULATION) {
    cudaFree(mp->d_displ);
    cudaFree(mp->d_veloc);
    cudaFree(mp->d_accel);

    if (mp->size_mpi_buffer > 0){
      cudaFree(mp->d_send_accel_buffer);
      cudaFree(mp->d_recv_accel_buffer);
      cudaFreeHost(mp->h_send_accel_buffer);
      cudaFreeHost(mp->h_recv_accel_buffer);
      if (mp->simulation_type == 3){
        cudaFree(mp->d_b_send_accel_buffer);
        cudaFree(mp->d_b_recv_accel_buffer);
      }
    }

    cudaFree(mp->d_rmassx);
    cudaFree(mp->d_rmassz);

    cudaFree(mp->d_phase_ispec_inner_elastic);
    cudaFree(mp->d_ispec_is_elastic);
    cudaFree(mp->d_inum_interfaces_elastic);

    if (*ABSORBING_CONDITIONS && mp->d_num_abs_boundary_faces > 0){
      cudaFree(mp->d_rho_vp);
      cudaFree(mp->d_rho_vs);
      cudaFree(mp->d_b_absorb_elastic_bottom);
      cudaFree(mp->d_b_absorb_elastic_left);
      cudaFree(mp->d_b_absorb_elastic_right);
      cudaFree(mp->d_b_absorb_elastic_top);
    }

    if (mp->simulation_type == 3) {
      cudaFree(mp->d_b_displ);
      cudaFree(mp->d_b_veloc);
      cudaFree(mp->d_b_accel);
      cudaFree(mp->d_rho_kl);
      cudaFree(mp->d_mu_kl);
      cudaFree(mp->d_kappa_kl);
      cudaFree(mp->d_b_dsxx);
      cudaFree(mp->d_b_dsxz);
      cudaFree(mp->d_b_dszz);
      cudaFree(mp->d_dsxx);
      cudaFree(mp->d_dsxz);
      cudaFree(mp->d_dszz);
      if (*APPROXIMATE_HESS_KL ) cudaFree(mp->d_hess_el_kl);
    }

    if (*ANISOTROPY) {
      cudaFree(mp->d_c11store);
      cudaFree(mp->d_c12store);
      cudaFree(mp->d_c13store);
      cudaFree(mp->d_c15store);
      cudaFree(mp->d_c23store);
      cudaFree(mp->d_c25store);
      cudaFree(mp->d_c33store);
      cudaFree(mp->d_c35store);
      cudaFree(mp->d_c55store);
    }

  if (*ATTENUATION_VISCOELASTIC) {
    cudaFree(mp->d_A_newmark_mu);
    cudaFree(mp->d_B_newmark_mu);
    cudaFree(mp->d_A_newmark_kappa);
    cudaFree(mp->d_B_newmark_kappa);
    cudaFree(mp->d_e1);
    cudaFree(mp->d_e11);
    cudaFree(mp->d_e13);
    cudaFree(mp->d_dux_dxl_old);
    cudaFree(mp->d_duz_dzl_old);
    cudaFree(mp->d_dux_dzl_plus_duz_dxl_old);
  }



  } // ELASTIC_SIMULATION

  // purely adjoint & kernel array
  if (mp->simulation_type == 3) {
    if (mp->nadj_rec_local > 0) {
      cudaFree(mp->d_adj_sourcearrays);
    }
  }

  // mesh pointer - not needed anymore
  free(mp);
}
