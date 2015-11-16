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

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
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
    //extern realw_texture d_hprimewgll_xx_tex;
    extern realw_texture d_wxgll_xx_tex;
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
    if(id<nibool[iinterface]) {

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
                                        int* h_number_receiver_global, int* h_ispec_selected_rec,
                                        int* nrec,int* nrec_local,
                                        realw * h_cosrot,realw * h_sinrot,
                                        int* SIMULATION_TYPE,
                                        int* USE_MESH_COLORING_GPU_f,
                                        int* nspec_acoustic,int* nspec_elastic,
                                        int* h_myrank,
                                        int* SAVE_FORWARD,
                                        realw* h_xir_store, realw* h_gammar_store ) {

  TRACE("prepare_constants_device");

  // allocates mesh parameter structure
  Mesh* mp = (Mesh*) malloc( sizeof(Mesh) );
  if (mp == NULL) exit_on_error("error allocating mesh pointer");
  *Mesh_pointer = (long)mp;

  // checks if NGLLX == 5
  if( *h_NGLLX != NGLLX ){
    exit_on_error("NGLLX must be 5 for CUDA devices");
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

  setConst_wxgll(h_wxgll,mp);


  // sets constant arrays
  setConst_hprime_xx(h_hprime_xx,mp);

  // setConst_hprime_zz(h_hprime_zz,mp); // only needed if NGLLX != NGLLY != NGLLZ

  setConst_hprimewgll_xx(h_hprimewgll_xx,mp);

  //setConst_hprimewgll_zz(h_hprimewgll_zz,mp); // only needed if NGLLX != NGLLY != NGLLZ

  // Using texture memory for the hprime-style constants is slower on
  // Fermi generation hardware, but *may* be faster on Kepler
  // generation. We will reevaluate this again, so might as well leave
  // in the code with with #USE_TEXTURES_FIELDS not-defined.
  #ifdef USE_TEXTURES_CONSTANTS
  {
    #ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

      const textureReference* d_hprime_xx_tex_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_hprime_xx_tex_ptr, "d_hprime_xx_tex"), 4101);
      print_CUDA_error_if_any(cudaBindTexture(0, d_hprime_xx_tex_ptr, mp->d_hprime_xx, &channelDesc, sizeof(realw)*(NGLL2)), 4001);

      const textureReference* d_wxgll_xx_tex_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_wxgll_xx_tex_ptr, "d_wxgll_xx_tex"), 4102);
      print_CUDA_error_if_any(cudaBindTexture(0, d_wxgll_xx_tex_ptr, mp->d_wxgll, &channelDesc, sizeof(realw)*(NGLL2)), 4002);

   #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

      print_CUDA_error_if_any(cudaBindTexture(0, &d_hprime_xx_tex, mp->d_hprime_xx, &channelDesc, sizeof(realw)*(NGLL2)), 4001);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_wxgll_xx_tex, mp->d_wxgll, &channelDesc, sizeof(realw)*(NGLLX)), 40013);
      //print_CUDA_error_if_any(cudaBindTexture(0, &d_hprimewgll_xx_tex, mp->d_hprimewgll_xx, &channelDesc, sizeof(realw)*(NGLL2)), 40010);
   #endif
  }
  #endif


  // mesh
  // Assuming NGLLX=5. Padded is then 32 (5^2+3)
  int size_padded = NGLL2_PADDED * (mp->NSPEC_AB);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xix, size_padded*sizeof(realw)),1001);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiz, size_padded*sizeof(realw)),1003);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammax, size_padded*sizeof(realw)),1007);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammaz, size_padded*sizeof(realw)),1009);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_kappav, size_padded*sizeof(realw)),1010);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_muv, size_padded*sizeof(realw)),1011);


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
  if( mp->num_interfaces_ext_mesh > 0 ){
    copy_todevice_int((void**)&mp->d_nibool_interfaces_ext_mesh,h_nibool_interfaces_ext_mesh,
                      mp->num_interfaces_ext_mesh);
    copy_todevice_int((void**)&mp->d_ibool_interfaces_ext_mesh,h_ibool_interfaces_ext_mesh,
                      (mp->num_interfaces_ext_mesh)*(mp->max_nibool_interfaces_ext_mesh));


  int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)(mp->max_nibool_interfaces_ext_mesh))/((double)blocksize)))*blocksize;


  }

  cudaStreamCreate(&mp->compute_stream);
  // copy stream (needed to transfer mpi buffers)
  if( mp->num_interfaces_ext_mesh * mp->max_nibool_interfaces_ext_mesh > 0 ){
    cudaStreamCreate(&mp->copy_stream);
  }


  // inner elements
  copy_todevice_int((void**)&mp->d_ispec_is_inner,h_ispec_is_inner,mp->NSPEC_AB);

  // absorbing boundaries
  mp->d_num_abs_boundary_faces = *h_num_abs_boundary_faces;
  if( mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0 ){
    copy_todevice_int((void**)&mp->d_abs_boundary_ispec,h_abs_boundary_ispec,mp->d_num_abs_boundary_faces);
    copy_todevice_int((void**)&mp->d_abs_boundary_ijk,h_abs_boundary_ij,
                      2*NGLL*(mp->d_num_abs_boundary_faces));
    copy_todevice_realw((void**)&mp->d_abs_boundary_normal,h_abs_boundary_normal,
                        NDIM*NGLL*(mp->d_num_abs_boundary_faces));
    copy_todevice_realw((void**)&mp->d_abs_boundary_jacobian2Dw,h_abs_boundary_jacobian1Dw,
                        NGLL*(mp->d_num_abs_boundary_faces));
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



  if( mp->nsources_local > 0){
    copy_todevice_realw((void**)&mp->d_source_time_function,h_source_time_function,(*NSTEP)*(mp->nsources_local));
    copy_todevice_realw((void**)&mp->d_sourcearrays,h_sourcearrays,mp->nsources_local*NDIM*NGLL2);
    copy_todevice_int((void**)&mp->d_ispec_selected_source,h_ispec_selected_source,mp->nsources_local);
    }



  // receiver stations
  mp->nrec_local = *nrec_local; // number of receiver located in this partition
  // note that:
  // size(number_receiver_global) = nrec_local
  // size(ispec_selected_rec) = nrec
  if( mp->nrec_local > 0 ){
    copy_todevice_int((void**)&mp->d_number_receiver_global,h_number_receiver_global,mp->nrec_local);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_seismograms,(mp->nrec_local)*sizeof(realw)*2),1303);
    print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_seismograms),sizeof(float)*(mp->nrec_local)*2),8004);
    copy_todevice_realw((void**)&mp->d_cosrot,h_cosrot,mp->nrec_local);
    copy_todevice_realw((void**)&mp->d_sinrot,h_sinrot,mp->nrec_local);
  }
  copy_todevice_int((void**)&mp->d_ispec_selected_rec,h_ispec_selected_rec,(*nrec));


#ifdef USE_MESH_COLORING_GPU
  mp->use_mesh_coloring_gpu = 1;
  if( ! *USE_MESH_COLORING_GPU_f ) exit_on_error("error with USE_MESH_COLORING_GPU constant; please re-compile\n");
#else
  // mesh coloring
  // note: this here passes the coloring as an option to the kernel routines
  //          the performance seems to be the same if one uses the pre-processing directives above or not
  mp->use_mesh_coloring_gpu = *USE_MESH_COLORING_GPU_f;
#endif

  // number of elements per domain
  mp->nspec_acoustic = *nspec_acoustic;
  mp->nspec_elastic = *nspec_elastic;

  copy_todevice_realw((void**)&mp->d_xir_store_loc,h_xir_store,(*nrec_local)*NGLLX);

  copy_todevice_realw((void**)&mp->d_gammar_store_loc,h_gammar_store,(*nrec_local)*NGLLX);



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
                                              int * h_ninterface_acoustic,int * h_inum_interfaces_acoustic,
                                              int* num_colors_outer_acoustic,
                                              int* num_colors_inner_acoustic,
                                              int* num_elem_colors_acoustic) {

  TRACE("prepare_fields_acoustic_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // allocates arrays on device (GPU)
  int size = mp->NGLOB_AB;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_potential_acoustic),sizeof(realw)*size),2001);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_potential_dot_acoustic),sizeof(realw)*size),2002);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_potential_dot_dot_acoustic),sizeof(realw)*size),2003);
  // initializes values to zero
  //print_CUDA_error_if_any(cudaMemset(mp->d_potential_acoustic,0,sizeof(realw)*size),2007);
  //print_CUDA_error_if_any(cudaMemset(mp->d_potential_dot_acoustic,0,sizeof(realw)*size),2007);
  //print_CUDA_error_if_any(cudaMemset(mp->d_potential_dot_dot_acoustic,0,sizeof(realw)*size),2007);

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
  if( mp->size_mpi_buffer_potential > 0 ){
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
    if( mp->num_free_surface_faces > 0 ){
      copy_todevice_int((void**)&mp->d_free_surface_ispec,free_surface_ispec,mp->num_free_surface_faces);
      copy_todevice_int((void**)&mp->d_free_surface_ijk,free_surface_ijk,
                        2*NGLLX*mp->num_free_surface_faces);

  }

  // absorbing boundaries
  if( mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0 ){
    // absorb_field array used for file i/o
    if(mp->simulation_type == 3 || ( mp->simulation_type == 1 && mp->save_forward )){
      print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_b_absorb_potential_left,mp->d_nspec_left*sizeof(realw)*NGLLX),2201);
      print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_b_absorb_potential_right,mp->d_nspec_right*sizeof(realw)*NGLLX),2201);
      print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_b_absorb_potential_top,mp->d_nspec_top*sizeof(realw)*NGLLX),2201);
      print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_b_absorb_potential_bottom,mp->d_nspec_bottom*sizeof(realw)*NGLLX),2201);

    }
  }


  // coupling with elastic parts
  if( *ELASTIC_SIMULATION && *num_coupling_ac_el_faces > 0 ){
    copy_todevice_int((void**)&mp->d_coupling_ac_el_ispec,coupling_ac_el_ispec,(*num_coupling_ac_el_faces));
    copy_todevice_int((void**)&mp->d_coupling_ac_el_ijk,coupling_ac_el_ijk,2*NGLL*(*num_coupling_ac_el_faces));
    copy_todevice_realw((void**)&mp->d_coupling_ac_el_normal,coupling_ac_el_normal,
                        2*NGLL*(*num_coupling_ac_el_faces));
    copy_todevice_realw((void**)&mp->d_coupling_ac_el_jacobian2Dw,coupling_ac_el_jacobian2Dw,
                        NGLL*(*num_coupling_ac_el_faces));
  }

  // mesh coloring
  if( mp->use_mesh_coloring_gpu ){
    mp->num_colors_outer_acoustic = *num_colors_outer_acoustic;
    mp->num_colors_inner_acoustic = *num_colors_inner_acoustic;
    mp->h_num_elem_colors_acoustic = (int*) num_elem_colors_acoustic;
  }

  mp->ninterface_acoustic = *h_ninterface_acoustic;
  copy_todevice_int((void**)&mp->d_inum_interfaces_acoustic,h_inum_interfaces_acoustic,mp->num_interfaces_ext_mesh);



#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_fields_acoustic_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_acoustic_adj_dev,
              PREPARE_FIELDS_ACOUSTIC_ADJ_DEV)(long* Mesh_pointer,
                                              int* APPROXIMATE_HESS_KL) {

  TRACE("prepare_fields_acoustic_adj_dev");

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // kernel simulations
  if( mp->simulation_type != 3 ) return;

  // allocates backward/reconstructed arrays on device (GPU)
  int size = mp->NGLOB_AB;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_potential_acoustic),sizeof(realw)*size),3014);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_potential_dot_acoustic),sizeof(realw)*size),3015);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_potential_dot_dot_acoustic),sizeof(realw)*size),3016);
  // initializes values to zero
  print_CUDA_error_if_any(cudaMemset(mp->d_b_potential_acoustic,0,sizeof(realw)*size),3007);
  print_CUDA_error_if_any(cudaMemset(mp->d_b_potential_dot_acoustic,0,sizeof(realw)*size),3007);
  print_CUDA_error_if_any(cudaMemset(mp->d_b_potential_dot_dot_acoustic,0,sizeof(realw)*size),3007);


  #ifdef USE_TEXTURES_FIELDS
  {
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

  // allocates kernels
  size = NGLL2*mp->NSPEC_AB;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_ac_kl),size*sizeof(realw)),3017);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_kappa_ac_kl),size*sizeof(realw)),3018);
  // initializes kernel values to zero
  print_CUDA_error_if_any(cudaMemset(mp->d_rho_ac_kl,0,size*sizeof(realw)),3019);
  print_CUDA_error_if_any(cudaMemset(mp->d_kappa_ac_kl,0,size*sizeof(realw)),3020);

  // preconditioner
  if( *APPROXIMATE_HESS_KL ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_hess_ac_kl),size*sizeof(realw)),3030);
    // initializes with zeros
    print_CUDA_error_if_any(cudaMemset(mp->d_hess_ac_kl,0,size*sizeof(realw)),3031);
  }

  // mpi buffer
  if( mp->size_mpi_buffer_potential > 0 ){
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
                                             realw *c55store,int* h_ninterface_elastic,int * h_inum_interfaces_elastic ){



  TRACE("prepare_fields_elastic_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  int size;


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
      print_CUDA_error_if_any(cudaGetTextureReference(&d_displ_tex_ref_ptr, "d_displ_tex"), 4001);
      print_CUDA_error_if_any(cudaBindTexture(0, d_displ_tex_ref_ptr, mp->d_displ, &channelDesc, sizeof(realw)*size), 4001);
      if( mp->use_mesh_coloring_gpu ){
        const textureReference* d_accel_tex_ref_ptr;
        print_CUDA_error_if_any(cudaGetTextureReference(&d_accel_tex_ref_ptr, "d_accel_tex"), 4003);
        print_CUDA_error_if_any(cudaBindTexture(0, d_accel_tex_ref_ptr, mp->d_accel, &channelDesc, sizeof(realw)*size), 4003);
      }
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_displ_tex, mp->d_displ, &channelDesc, sizeof(realw)*size), 4001);
      if( mp->use_mesh_coloring_gpu ) print_CUDA_error_if_any(cudaBindTexture(0, &d_accel_tex, mp->d_accel, &channelDesc, sizeof(realw)*size), 4003);
    #endif
  }
  #endif


  // debug
  //synchronize_mpi();

  // MPI buffer
  mp->size_mpi_buffer = NDIM * (mp->num_interfaces_ext_mesh) * (mp->max_nibool_interfaces_ext_mesh);
  if( mp->size_mpi_buffer > 0 ){
    // note: Allocate pinned mpi-buffers.
    //       MPI buffers use pinned memory allocated by cudaMallocHost, which
    //       enables the use of asynchronous memory copies from host <-> device
    // send buffer
    print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_send_accel_buffer),sizeof(float)*(mp->size_mpi_buffer)),8004);
    //mp->send_buffer = (float*)malloc((mp->size_mpi_buffer)*sizeof(float));
    // adjoint
    //print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_send_b_accel_buffer),sizeof(float)*(mp->size_mpi_buffer)),8004);
    // mp->b_send_buffer = (float*)malloc((size_mpi_buffer)*sizeof(float));
    // receive buffer
    print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_recv_accel_buffer),sizeof(float)*(mp->size_mpi_buffer)),8004);
    mp->recv_buffer = (float*)malloc((mp->size_mpi_buffer)*sizeof(float));

    // non-pinned buffer
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_send_accel_buffer),mp->size_mpi_buffer*sizeof(realw)),4004);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_recv_accel_buffer),mp->size_mpi_buffer*sizeof(realw)),4004);


    // adjoint
    if( mp->simulation_type == 3 ){
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



  // debug
  //synchronize_mpi();

  // absorbing conditions
  if( mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0){

    // debug
    //printf("prepare_fields_elastic_device: rank %d - absorbing boundary setup\n",mp->myrank);

    // non-padded arrays
    // rho_vp, rho_vs non-padded; they are needed for stacey boundary condition
    copy_todevice_realw((void**)&mp->d_rho_vp,rho_vp,NGLL2*mp->NSPEC_AB);
    copy_todevice_realw((void**)&mp->d_rho_vs,rho_vs,NGLL2*mp->NSPEC_AB);

    // absorb_field array used for file i/o

  if( mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0 ){
    // absorb_field array used for file i/o
    if(mp->simulation_type == 3 || ( mp->simulation_type == 1 && mp->save_forward )){
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
  if( *ANISOTROPY ){
    // debug
    //printf("prepare_fields_elastic_device: rank %d - attenuation setup\n",mp->myrank);
    //synchronize_mpi();

    // Assuming NGLLX==5. Padded is then 32 (5^2+3)
    int size_padded = NGLL2_PADDED * (mp->NSPEC_AB);

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


  // mesh coloring
  if( mp->use_mesh_coloring_gpu ){
    mp->num_colors_outer_elastic = *num_colors_outer_elastic;
    mp->num_colors_inner_elastic = *num_colors_inner_elastic;
    mp->h_num_elem_colors_elastic = (int*) num_elem_colors_elastic;
  }


  mp->ninterface_elastic = *h_ninterface_elastic;
  copy_todevice_int((void**)&mp->d_inum_interfaces_elastic,h_inum_interfaces_elastic,mp->num_interfaces_ext_mesh);

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
  if( mp->simulation_type != 3 ) return;

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
      print_CUDA_error_if_any(cudaGetTextureReference(&d_b_displ_tex_ref_ptr, "d_b_displ_tex"), 4001);
      print_CUDA_error_if_any(cudaBindTexture(0, d_b_displ_tex_ref_ptr, mp->d_b_displ, &channelDesc, sizeof(realw)*size), 4001);
      if( mp->use_mesh_coloring_gpu ){
        const textureReference* d_b_accel_tex_ref_ptr;
        print_CUDA_error_if_any(cudaGetTextureReference(&d_b_accel_tex_ref_ptr, "d_b_accel_tex"), 4003);
        print_CUDA_error_if_any(cudaBindTexture(0, d_b_accel_tex_ref_ptr, mp->d_b_accel, &channelDesc, sizeof(realw)*size), 4003);
      }
    #else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_b_displ_tex, mp->d_b_displ, &channelDesc, sizeof(realw)*size), 4001);
      if( mp->use_mesh_coloring_gpu ) print_CUDA_error_if_any(cudaBindTexture(0, &d_b_accel_tex, mp->d_b_accel, &channelDesc, sizeof(realw)*size), 4003);
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
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_kl),size*sizeof(realw)),5204);
  // initializes kernel values to zero
  print_CUDA_error_if_any(cudaMemset(mp->d_rho_kl,0,size*sizeof(realw)),5214);


    // isotropic kernels
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_mu_kl),size*sizeof(realw)),5206);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_kappa_kl),size*sizeof(realw)),5207);
    print_CUDA_error_if_any(cudaMemset(mp->d_mu_kl,0,size*sizeof(realw)),5216);
    print_CUDA_error_if_any(cudaMemset(mp->d_kappa_kl,0,size*sizeof(realw)),5217);

    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_dsxx),size*sizeof(realw)),5207);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_dsxz),size*sizeof(realw)),5207);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_dszz),size*sizeof(realw)),5207);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_dsxx),size*sizeof(realw)),5207);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_dsxz),size*sizeof(realw)),5207);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_dszz),size*sizeof(realw)),5207);

  // approximate hessian kernel
  if( *APPROXIMATE_HESS_KL ){
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
                                              int* islice_selected_rec,
                                              int* nadj_rec_local,
                                              int* nrec,realw* h_source_adjointe,int* NSTEP) {

  TRACE("prepare_sim2_or_3_const_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // adjoint source arrays
  mp->nadj_rec_local = *nadj_rec_local;
  if( mp->nadj_rec_local > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_adj_sourcearrays,
                                       (mp->nadj_rec_local)*2*NGLL2*sizeof(realw)),6003);

    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_pre_computed_irec,
                                       (mp->nadj_rec_local)*sizeof(int)),6004);

    // prepares local irec array:
    // the irec_local variable needs to be precomputed (as
    // h_pre_comp..), because normally it is in the loop updating accel,
    // and due to how it's incremented, it cannot be parallelized
    int* h_pre_computed_irec = (int*) malloc( (mp->nadj_rec_local)*sizeof(int) );
    if( h_pre_computed_irec == NULL ) exit_on_error("prepare_sim2_or_3_const_device: h_pre_computed_irec not allocated\n");

    int irec_local = 0;
    for(int irec = 0; irec < *nrec; irec++) {
      if(mp->myrank == islice_selected_rec[irec]) {
        irec_local++;
        h_pre_computed_irec[irec_local-1] = irec;
      }
    }
    // checks if all local receivers have been found
    if( irec_local != mp->nadj_rec_local ) exit_on_error("prepare_sim2_or_3_const_device: irec_local not equal\n");

    // copies values onto GPU
    print_CUDA_error_if_any(cudaMemcpy(mp->d_pre_computed_irec,h_pre_computed_irec,
                                       (mp->nadj_rec_local)*sizeof(int),cudaMemcpyHostToDevice),6010);
    free(h_pre_computed_irec);

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
                                      int* APPROXIMATE_HESS_KL) {


TRACE("prepare_cleanup_device");

  // frees allocated memory arrays
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  // frees memory on GPU
  // mesh
  cudaFree(mp->d_xix);
  cudaFree(mp->d_xiz);
  cudaFree(mp->d_gammax);
  cudaFree(mp->d_gammaz);
  cudaFree(mp->d_muv);

  // absorbing boundaries
  if( *ABSORBING_CONDITIONS && mp->d_num_abs_boundary_faces > 0 ){
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
  if( mp->num_interfaces_ext_mesh > 0 ){
    cudaFree(mp->d_nibool_interfaces_ext_mesh);
    cudaFree(mp->d_ibool_interfaces_ext_mesh);
  }

  // global indexing
  cudaFree(mp->d_ispec_is_inner);
  cudaFree(mp->d_ibool);

  // sources
  if( mp->nsources_local > 0){
    cudaFree(mp->d_sourcearrays);
    cudaFree(mp->d_source_time_function);
  cudaFree(mp->d_ispec_selected_source);
  }



  // receivers
  if( mp->nrec_local > 0 ){
  cudaFree(mp->d_number_receiver_global);cudaFree(mp->d_seismograms);
  cudaFree(mp->d_cosrot),cudaFree(mp->d_sinrot);
  }

  cudaFree(mp->d_ispec_selected_rec);

  cudaFree(mp->d_gammar_store_loc);
  cudaFree(mp->d_xir_store_loc);

  // ACOUSTIC arrays
  if( *ACOUSTIC_SIMULATION ){
    cudaFree(mp->d_potential_acoustic);
    cudaFree(mp->d_potential_dot_acoustic);
    cudaFree(mp->d_potential_dot_dot_acoustic);
    cudaFree(mp->d_send_potential_dot_dot_buffer);
    cudaFree(mp->d_rmass_acoustic);
    cudaFree(mp->d_rhostore);
    cudaFree(mp->d_kappastore);
    cudaFree(mp->d_phase_ispec_inner_acoustic);
    cudaFree(mp->d_ispec_is_acoustic);
    cudaFree(mp->d_inum_interfaces_acoustic);

    if( mp->simulation_type == 3 ) {
      cudaFree(mp->d_b_potential_acoustic);
      cudaFree(mp->d_b_potential_dot_acoustic);
      cudaFree(mp->d_b_potential_dot_dot_acoustic);
      cudaFree(mp->d_rho_ac_kl);
      cudaFree(mp->d_kappa_ac_kl);
      if( *APPROXIMATE_HESS_KL) cudaFree(mp->d_hess_ac_kl);
      if( mp->size_mpi_buffer_potential > 0 ) cudaFree(mp->d_b_send_potential_dot_dot_buffer);
  }


 if( *ABSORBING_CONDITIONS && mp->d_num_abs_boundary_faces > 0){
 if(mp->simulation_type == 3 || ( mp->simulation_type == 1 && mp->save_forward )){
      cudaFree(mp->d_b_absorb_potential_bottom);
      cudaFree(mp->d_b_absorb_potential_left);
      cudaFree(mp->d_b_absorb_potential_right);
      cudaFree(mp->d_b_absorb_potential_top); }}

  } // ACOUSTIC_SIMULATION

  // ELASTIC arrays
  if( *ELASTIC_SIMULATION ){
    cudaFree(mp->d_displ);
    cudaFree(mp->d_veloc);
    cudaFree(mp->d_accel);

    cudaFree(mp->d_send_accel_buffer);
    if( mp->simulation_type == 3) cudaFree(mp->d_b_send_accel_buffer);

    cudaFree(mp->d_rmassx);
    cudaFree(mp->d_rmassz);

    cudaFree(mp->d_phase_ispec_inner_elastic);
    cudaFree(mp->d_ispec_is_elastic);
    cudaFree(mp->d_inum_interfaces_elastic);



    if( *ABSORBING_CONDITIONS && mp->d_num_abs_boundary_faces > 0){
      cudaFree(mp->d_rho_vp);
      cudaFree(mp->d_rho_vs);
      cudaFree(mp->d_b_absorb_elastic_bottom);
      cudaFree(mp->d_b_absorb_elastic_left);
      cudaFree(mp->d_b_absorb_elastic_right);
      cudaFree(mp->d_b_absorb_elastic_top);

    }

    if( mp->simulation_type == 3 ) {
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

      if( *APPROXIMATE_HESS_KL ) cudaFree(mp->d_hess_el_kl);
    }

    if( *ANISOTROPY ){
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


  } // ELASTIC_SIMULATION

  // purely adjoint & kernel array
  if(mp->simulation_type == 3 ){
    if(mp->nadj_rec_local > 0 ){
      cudaFree(mp->d_adj_sourcearrays);
      cudaFree(mp->d_pre_computed_irec);
    }
  }


  // mesh pointer - not needed anymore
  free(mp);
}
