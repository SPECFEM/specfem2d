!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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
!
! United States and French Government Sponsorship Acknowledged.

module constants

  include "constants.h"

end module constants

!=====================================================================

module specfem_par

! main parameter module for specfem simulations

  use constants

  implicit none

! parameters deduced from parameters read from file
  integer :: NPROC,nproc_read_from_database
  integer :: NSPEC_AB, NGLOB_AB

! mesh parameters
  integer, dimension(:,:,:), allocatable :: ibool
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore,zstore

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: xix,xiz,gammax,gammaz,jacobian

! material properties
  ! isotropic
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: kappastore,mustore
  real(kind=CUSTOM_REAL) :: mul_unrelaxed_elastic, kappal

! density
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhostore

! CUDA mesh pointer<->integer wrapper
  integer(kind=8) :: Mesh_pointer

! Global GPU toggle. Set in Par_file
  logical :: GPU_MODE

! use integer array to store topography values
  integer :: NX_TOPO,NY_TOPO
  integer, dimension(:,:), allocatable :: itopo_bathy

! absorbing boundary arrays (for all boundaries) - keeps all infos, allowing for irregular surfaces
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: abs_boundary_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: abs_boundary_jacobian1Dw
  integer, dimension(:,:,:), allocatable :: abs_boundary_ij
  integer, dimension(:), allocatable :: abs_boundary_ispec
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,NSPEC2D_TOP
  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top

! free surface arrays
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: free_surface_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: free_surface_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: free_surface_ij
  integer, dimension(:), allocatable :: free_ac_ispec
  integer :: num_free_surface_faces

! VM for new method
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: Veloc_dsm_boundary,Tract_dsm_boundary

! attenuation
  integer :: NSPEC_ATTENUATION_AB,NSPEC_ATTENUATION_AB_kappa
  character(len=256) prname_Q



! time scheme
  double precision deltat,deltatover2,deltatsquareover2
  real(kind=CUSTOM_REAL) deltatf,deltatover2f,deltatsquareover2f

  double precision  xixl,xizl,gammaxl,gammazl

! time loop step
  integer :: it

! VM for new  method
  integer :: it_dsm

! parameters for the source
  integer, dimension(:), allocatable :: islice_selected_source,ispec_selected_source
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sourcearray
  double precision, dimension(:), allocatable :: Mxx,Mzz,Mxz
  double precision, dimension(:), allocatable :: xi_source,gamma_source
  double precision, dimension(:), allocatable :: tshift_src,hdur,hdur_gaussian,hdur_tiny
  double precision, dimension(:), allocatable :: utm_x_source,utm_y_source

  double precision :: t0
  real(kind=CUSTOM_REAL) :: stf_used_total
  integer :: NSOURCES,nsources_local
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: source_time_function
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable ::  source_time_function_loc
  integer, dimension(:), allocatable :: num_src_loc

! source encoding
  ! for acoustic sources: takes +/- 1 sign, depending on sign(Mxx)[ = sign(Myy) = sign(Mzz)
  ! since they have to equal in the acoustic setting]
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: pm1_source_encoding

! receiver information
  character(len=256) :: rec_filename,filtered_rec_filename,dummystring
  integer :: nrec,nrec_tot_found
  integer :: nrec_simulation
  integer, dimension(:), allocatable :: islice_selected_rec,ispec_selected_rec
  double precision, dimension(:), allocatable :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(:,:), allocatable :: hpxir_store,hpetar_store,hpgammar_store
  integer, dimension(:,:), allocatable  :: gather_ispec_selected_rec

! timing information for the stations
  double precision, allocatable, dimension(:,:,:) :: nu
  character(len=MAX_LENGTH_STATION_NAME), allocatable, dimension(:) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), allocatable, dimension(:) :: network_name



! Gauss-Lobatto-Legendre points of integration and weights
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wzgll
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll


! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz


! Lagrange interpolators at receivers
  double precision, dimension(:), allocatable :: hxir,hpxir,hgammar,hpgammar
  double precision, dimension(:,:), allocatable :: hxir_store,hgammar_store

! proc numbers for MPI
  integer :: myrank, sizeprocs

! timer MPI
  double precision, external :: wtime
  double precision :: time_start

! parameters for a force source located exactly at a grid point
  logical :: USE_FORCE_POINT_SOURCE
  double precision, dimension(:), allocatable :: factor_force_source
  double precision, dimension(:), allocatable :: comp_dir_vect_source_E
  double precision, dimension(:), allocatable :: comp_dir_vect_source_N
  double precision, dimension(:), allocatable :: comp_dir_vect_source_Z_UP

! parameters
  integer :: SIMULATION_TYPE
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,UTM_PROJECTION_ZONE
  integer :: IMODEL,NGNOD,NGNOD2D

  double precision :: DT,OLSEN_ATTENUATION_RATIO,f0_FOR_PML

  logical :: ATTENUATION,USE_OLSEN_ATTENUATION, &
            TOPOGRAPHY, &
            STACEY_INSTEAD_OF_FREE_SURFACE

  logical :: FULL_ATTENUATION_SOLID,PML_BOUNDARY_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE

  logical :: GRAVITY

  logical :: SAVE_FORWARD,SAVE_MESH_FILES

  logical :: USE_RICKER_TIME_FUNCTION,PRINT_SOURCE_TIME_FUNCTION

  logical :: SUPPRESS_UTM_PROJECTION

  integer :: NTSTEP_BETWEEN_OUTPUT_INFO

! parameters read from mesh parameter file
  integer :: NPROC_XI,NPROC_ETA
  double precision :: LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

  character(len=256) OUTPUT_FILES,LOCAL_PATH,TOMOGRAPHY_PATH,prname,dsmname,TRAC_PATH

  logical :: ADIOS_ENABLED
  logical :: ADIOS_FOR_DATABASES, ADIOS_FOR_MESH, ADIOS_FOR_FORWARD_ARRAYS, &
             ADIOS_FOR_KERNELS

! names of the data files for all the processors in MPI
  character(len=256) outputname, outputname2

! for assembling in case of external mesh
!  integer :: num_interfaces_ext_mesh
  integer :: max_nibool_interfaces_ext_mesh
  integer, dimension(:), allocatable :: nibool_interfaces_ext_mesh
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_init, ibool_interfaces_ext_mesh
  integer, dimension(:), allocatable  :: my_neighbours
  integer, dimension(:), allocatable  :: tab_requests_send_recv_scalar
  integer, dimension(:), allocatable  :: b_tab_requests_send_recv_scalar
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_scalar_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_recv_scalar_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_buffer_recv_scalar_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_recv_vector_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_recv_vector_ext_mesh
  integer, dimension(:), allocatable  :: tab_requests_send_recv_vector
  integer, dimension(:), allocatable  :: b_tab_requests_send_recv_vector
  double precision, dimension(:), allocatable  :: tab_requests_send_recv_vector_c



! for detecting surface receivers and source in case of external mesh
  logical, dimension(:), allocatable :: iglob_is_surface_external_mesh
  logical, dimension(:), allocatable :: ispec_is_surface_external_mesh

! MPI partition surfaces
  logical, dimension(:), allocatable :: ispec_is_inner

! maximum speed in velocity model
  real(kind=CUSTOM_REAL):: model_speed_max


! ADJOINT parameters

  ! time scheme
  double precision b_deltat, b_deltatover2, b_deltatsquareover2
  real(kind=CUSTOM_REAL) b_deltatf,b_deltatover2f,b_deltatsquareover2f


  ! Moho mesh
  real(CUSTOM_REAL), dimension(:,:,:),allocatable :: normal_moho_top
  real(CUSTOM_REAL), dimension(:,:,:),allocatable :: normal_moho_bot
  integer,dimension(:,:,:),allocatable :: ijk_moho_top, ijk_moho_bot
  integer,dimension(:),allocatable :: ibelm_moho_top, ibelm_moho_bot
  integer :: NSPEC_BOUN,NSPEC2D_MOHO
  logical, dimension(:),allocatable :: is_moho_top, is_moho_bot

  ! adjoint sources
  character(len=256) adj_source_file
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: adj_sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: adj_sourcearrays
  integer :: nadj_rec_local
  real(kind=CUSTOM_REAL),  dimension(:,:,:), allocatable :: source_adjointe
  real(kind=CUSTOM_REAL),  dimension(:,:), allocatable :: xir_store_loc, gammar_store_loc


  ! adjoint source frechet derivatives
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: Mxx_der,Myy_der,&
    Mzz_der,Mxy_der,Mxz_der,Myz_der
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sloc_der
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: seismograms_eps

  ! adjoint elements
  integer :: NSPEC_ADJOINT, NGLOB_ADJOINT

  ! length of reading blocks
  integer :: NTSTEP_BETWEEN_READ_ADJSRC

  ! parameter module for noise simulations
  integer :: irec_master_noise, NOISE_TOMOGRAPHY
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sigma_kl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: noise_sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: noise_surface_movie
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
             normal_x_noise,normal_y_noise,normal_z_noise, mask_noise

  
  !Simulation acoustique
  integer :: num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic
  integer, dimension(:,:), allocatable  :: acoustic_surface
  integer :: nelem_acoustic_surface
  integer, dimension(:,:), allocatable :: phase_ispec_inner_acoustic
  integer, dimension(:), allocatable :: cote_abs
  integer, dimension(:), allocatable :: ib_left,ib_right,ib_top,ib_bottom
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable ::b_absorb_acoustic_right,b_absorb_acoustic_left&
                                                         ,b_absorb_acoustic_top,b_absorb_acoustic_bottom
  integer :: nspec_left,nspec_right,nspec_bottom,nspec_top
  


  ! Simulation elastique, couplage elastique/acoustique

  logical :: any_elastic, any_acoustic, any_poroelastic
  integer :: num_fluid_solid_edges
  integer, dimension(:), allocatable :: coupling_ac_el_ispec
  integer, dimension(:,:,:), allocatable :: coupling_ac_el_ij
  integer, dimension(NGLLX,NEDGES) :: ivalue,jvalue,ivalue_inverse,jvalue_inverse
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: coupling_ac_el_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: coupling_ac_el_jacobian1Dw
  integer :: num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic
  integer, dimension(:,:), allocatable :: phase_ispec_inner_elastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: displ_2D,veloc_2D,accel_2D
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_displ_2D,b_veloc_2D,b_accel_2D

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_absorb_elastic_left,b_absorb_elastic_right,b_absorb_elastic_top,b_absorb_elastic_bottom


  !ANISOTROPY

  logical :: ANY_ANISOTROPY
 ! anisotropic
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
            c11store,c12store,c13store,c15store,c23store,c25store,c33store,c35store,c55store

  !stacey
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rho_vp,rho_vs

! mesh coloring
  integer :: num_colors_outer_acoustic,num_colors_inner_acoustic
  integer, dimension(:), allocatable :: num_elem_colors_acoustic
  integer :: num_colors_outer_elastic,num_colors_inner_elastic
  integer, dimension(:), allocatable :: num_elem_colors_elastic

   ! Source

   integer, dimension(:), allocatable :: source_type


   ! Sismogrammes

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: cosrot_irecf, sinrot_irecf




end module specfem_par
