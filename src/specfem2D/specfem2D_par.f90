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

module specfem_par

! main parameter module for specfem simulations

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,MAX_LENGTH_NETWORK_NAME,MAX_LENGTH_STATION_NAME, &
    NEDGES,NGLLX,NGLLZ,NGLJ,NDIM

  use shared_parameters

  implicit none
  !=====================================================================
  ! input for simulation (its beginning)
  !=====================================================================
  !---------------------------------------------------------------------
  ! for material information
  !---------------------------------------------------------------------
  integer :: numat
  logical :: assign_external_model

  ! poroelastic and elastic coefficients
  double precision, dimension(:,:,:), allocatable :: poroelastcoef
  logical, dimension(:), allocatable :: already_shifted_velocity
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: vpext,vsext,rhoext
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: QKappa_attenuationext,Qmu_attenuationext

  ! anisotropy parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,c22ext
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: c11_k, c13_k, c15_k,c33_k, c35_k, c55_k
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: c11_kl, c13_kl, c15_kl,c33_kl, c35_kl, c55_kl

  logical :: all_anisotropic
  logical, dimension(:), allocatable :: ispec_is_anisotropic

  double precision, dimension(:,:), allocatable :: anisotropy

  ! for attenuation
  double precision, dimension(:), allocatable  :: QKappa_attenuation
  double precision, dimension(:), allocatable  :: Qmu_attenuation

  integer :: nspec_ATT_el,nspec_ATT_ac,nglob_att
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2
  real(kind=CUSTOM_REAL), dimension(:,:,:) , allocatable :: Mu_nu1,Mu_nu2

  ! material
  ! density
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhostore
  ! isotropic moduli
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: kappastore,mustore
  ! anisotropic
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    c11store,c12store,c13store,c15store,c23store,c25store,c33store,c35store,c55store

  ! for absorbing boundaries
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rho_vp,rho_vs

  !---------------------------------------------------------------------
  ! for boundary condition (physical BC or artificial BC)
  !---------------------------------------------------------------------
  logical :: anyabs_glob

  ! PML
  logical, dimension(:), allocatable :: ispec_is_PML
  integer :: nspec_PML

  integer, dimension(:), allocatable :: region_CPML
  integer, dimension(:), allocatable :: spec_to_PML

  logical, dimension(:,:), allocatable :: which_PML_elem
  logical, dimension(:), allocatable  :: mask_ibool_PML

  double precision, dimension(:,:,:), allocatable :: &
    K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store
  double precision :: min_distance_between_CPML_parameter

  ! Stacey BC
  logical, dimension(:,:), allocatable  :: codeabs
  integer, dimension(:), allocatable  :: typeabs
  ! for detection of corner element on absorbing boundary
  logical, dimension(:,:), allocatable  :: codeabs_corner

  ! horizontal periodicity distance for periodic conditions
  logical, dimension(:), allocatable :: this_ibool_is_a_periodic_edge

  ! edge detection
  integer, dimension(NEDGES) :: i_begin,j_begin,i_end,j_end
  integer, dimension(NGLLX,NEDGES) :: ivalue,jvalue,ivalue_inverse,jvalue_inverse

  ! fluid/solid interface
  integer :: num_fluid_solid_edges
  logical :: coupled_acoustic_elastic,any_fluid_solid_edges
  integer, dimension(:), allocatable :: fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
                                        fluid_solid_elastic_ispec,fluid_solid_elastic_iedge

  ! fluid/porous interface
  integer :: num_fluid_poro_edges
  logical :: coupled_acoustic_poro,any_fluid_poro_edges
  integer, dimension(:), allocatable :: fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge, &
                                        fluid_poro_poroelastic_ispec,fluid_poro_poroelastic_iedge

  ! solid/porous interface
  integer :: num_solid_poro_edges
  logical :: coupled_elastic_poro, any_solid_poro_edges
  integer, dimension(:), allocatable :: solid_poro_elastic_ispec,solid_poro_elastic_iedge, &
                                        solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge

  ! for absorbing boundary mass matrix correction
  integer, dimension(:), allocatable :: ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro, &
            iend_edge3_poro,ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro

  !---------------------------------------------------------------------
  ! for source-receiver information
  !---------------------------------------------------------------------
  ! source description
  integer, dimension(:), allocatable :: source_type,time_function_type
  character(len=MAX_STRING_LEN), dimension(:), allocatable :: name_of_source_file
  double precision, dimension(:), allocatable :: burst_band_width

  ! source locations
  double precision, dimension(:), allocatable :: x_source,z_source
  double precision, dimension(:), allocatable :: xi_source,gamma_source

  double precision, dimension(:), allocatable :: Mxx,Mzz,Mxz
  double precision, dimension(:), allocatable :: f0_source,tshift_src,factor,anglesource

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: sourcearrays
  double precision :: t0

  integer, dimension(:), allocatable :: ispec_selected_source
  integer, dimension(:), allocatable :: islice_selected_source
  integer, dimension(:), allocatable :: iglob_source

  ! source time function
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: source_time_function

  ! for absorbing and acoustic free surface conditions
  !acoustic free surface
  integer :: nelem_acoustic_surface
  integer, dimension(:,:), allocatable :: acoustic_surface
  integer, dimension(:,:), allocatable :: acoustic_edges
  logical :: any_acoustic_edges

  ! perform a forcing of an acoustic medium with a rigid boundary
  integer :: nelem_acforcing
  logical, dimension(:,:), allocatable  :: codeacforcing
  integer, dimension(:), allocatable  :: typeacforcing
  integer, dimension(:), allocatable :: numacforcing, &
     ibegin_edge1_acforcing,iend_edge1_acforcing,ibegin_edge3_acforcing,iend_edge3_acforcing, &
     ibegin_edge4_acforcing,iend_edge4_acforcing,ibegin_edge2_acforcing,iend_edge2_acforcing
  integer :: nspec_left_acforcing,nspec_right_acforcing,nspec_bottom_acforcing,nspec_top_acforcing
  integer, dimension(:), allocatable :: ib_left_acforcing,ib_right_acforcing,ib_bottom_acforcing,ib_top_acforcing

  ! Variables for forcing
  logical, dimension(:), allocatable :: iglob_is_forced,acoustic_iglob_is_forced,elastic_iglob_is_forced
  integer :: nLines ! Number of lines in input file
  double precision, dimension(:), allocatable  :: zmode,realMode,imagMode,modeAmplitude

  ! for plane wave incidence
  double precision :: anglesource_refl, c_inc, c_refl
  double precision, dimension(2) :: A_plane, B_plane, C_plane
  double precision :: time_offset
  ! beyond critical angle
  integer , dimension(:), allocatable :: left_bound,right_bound,bot_bound
  double precision , dimension(:,:), allocatable :: v0x_left,v0z_left,v0x_right,v0z_right,v0x_bot,v0z_bot
  double precision , dimension(:,:), allocatable :: t0x_left,t0z_left,t0x_right,t0z_right,t0x_bot,t0z_bot
  integer count_left,count_right,count_bottom
  logical :: over_critical_angle

  ! receivers description
  ! timing information for the stations
  character(len=MAX_LENGTH_STATION_NAME), allocatable, dimension(:) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), allocatable, dimension(:) :: network_name

  ! total number of receiver stations
  integer :: nrec
  ! local receivers in this slice
  integer :: nrecloc
  double precision :: xirec,gammarec
  integer, dimension(:), allocatable :: recloc
  integer, dimension(:), allocatable :: islice_selected_rec
  integer, dimension(:), allocatable :: ispec_selected_rec,ispec_selected_rec_loc
  double precision, dimension(:), allocatable :: xi_receiver,gamma_receiver,st_xval,st_zval

  ! tangential detection for source or receivers
  double precision, dimension(:), allocatable :: anglerec_irec
  double precision, dimension(:), allocatable :: cosrot_irec, sinrot_irec
  double precision, dimension(:), allocatable :: x_final_receiver, z_final_receiver

  integer, dimension(:), allocatable :: source_courbe_eros

  integer  :: nnodes_tangential_curve
  double precision, dimension(:,:), allocatable  :: nodes_tangential_curve
  logical  :: any_tangential_curve

  integer  :: n1_tangential_detection_curve
  integer, dimension(4)  :: n_tangential_detection_curve
  integer, dimension(:), allocatable  :: rec_tangential_detection_curve
  double precision :: distmin, dist_current, anglesource_recv
  double precision, dimension(:), allocatable :: dist_tangential_detection_curve

  !---------------------------------------------------------------------
  ! for SEM discretization of the model
  !---------------------------------------------------------------------

  ! Gauss-Lobatto-Legendre points and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll

  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  real(kind=CUSTOM_REAL), dimension(NGLLZ) :: wzgll

  ! derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

  double precision, dimension(:,:,:), allocatable :: shape2D
  double precision, dimension(:,:,:,:), allocatable :: dershape2D

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: xix,xiz,gammax,gammaz,jacobian
  integer, dimension(:,:,:), allocatable :: ibool

  integer, dimension(:,:), allocatable  :: knods
  integer, dimension(:), allocatable :: kmato

  integer, dimension(:), allocatable :: numabs, &
     ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3,ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2

  ! Lagrange interpolators at receivers
  double precision, dimension(:), allocatable :: hxir,hgammar,hpxir,hpgammar

  ! Lagrange interpolators at sources
  double precision, dimension(:), allocatable :: hxis,hgammas,hpxis,hpgammas
  double precision, dimension(:,:), allocatable :: hxis_store,hgammas_store

  !---------------------------------------------------------------------
  ! AXISYM parameters
  !---------------------------------------------------------------------
  ! Number of elements on the symmetry axis
  integer :: nelem_on_the_axis
  ! Flag to know if an element is on the axis
  logical, dimension(:), allocatable :: is_on_the_axis
  integer, dimension(:), allocatable :: ispec_of_axial_elements

  ! Gauss-Lobatto-Jacobi points and weights
  double precision, dimension(NGLJ) :: xiglj
  real(kind=CUSTOM_REAL), dimension(NGLJ) :: wxglj

  ! derivatives of GLJ polynomials
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLJ) :: hprimeBar_xx,hprimeBarwglj_xx

  ! Shape functions (and their derivatives) evaluated at the GLJ points
  double precision, dimension(:,:), allocatable :: flagrange_GLJ


  !---------------------------------------------------------------------
  ! for time discretization
  !---------------------------------------------------------------------
  ! for LDDRK46
  integer :: i_stage,stage_time_scheme

  ! coefficients of the explicit Newmark time scheme
  double precision :: deltat,deltatover2,deltatsquareover2

  ! for backward simulation in adjoint inversion
  double precision :: b_deltatover2,b_deltatsquareover2,b_deltat ! coefficients of the explicit Newmark time scheme

  ! current time
  double precision :: timeval

  ! UNDO_ATTENUATION_AND_OR_PML
  integer :: NSUBSET_ITERATIONS
  integer :: iteration_on_subset,it_of_this_subset
  integer :: it_subset_end

  ! to determine date and time at which the run will finish
  double precision :: timestamp_seconds_start

  !---------------------------------------------------------------------
  ! global variable shared by acoustic/elastic/poroelastic simulation
  !---------------------------------------------------------------------
  double precision, dimension(:,:), allocatable :: coord
  double precision, dimension(:,:), allocatable :: coorg

  !---------------------------------------------------------------------
  ! for acoustic simulation
  !---------------------------------------------------------------------
  ! for acoustic medium
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic

  ! PML
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic_old

  ! RK time schemes
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_dot_acoustic_LDDRK, potential_acoustic_LDDRK
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_dot_acoustic_temp
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic_init_rk, potential_dot_acoustic_init_rk
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: potential_dot_dot_acoustic_rk, potential_dot_acoustic_rk

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic_adj_coupling

  ! for acoustic detection
  integer :: nglob_acoustic

  ! number of purely acoustic elements in this slice
  integer :: nspec_acoustic
  integer :: nspec_acoustic_b

  ! local flag to determine if any acoustic elements in this slice
  logical :: any_acoustic

  ! global flag for acoustic simulations
  logical :: ACOUSTIC_SIMULATION

  logical, dimension(:), allocatable :: ispec_is_acoustic

  ! MPI inner/outer
  integer :: nspec_inner_acoustic,nspec_outer_acoustic
  integer :: num_phase_ispec_acoustic
  integer, dimension(:,:), allocatable :: phase_ispec_inner_acoustic

  ! inverse mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_acoustic

  ! the variable for PML
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_potential_acoustic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                          rmemory_acoustic_dux_dx,rmemory_acoustic_dux_dz

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_potential_acoustic_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                          rmemory_acoustic_dux_dx_LDDRK,rmemory_acoustic_dux_dz_LDDRK

  ! for backward simulation in adjoint inversion
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic

  ! PML
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: b_potential_acoustic_old

  ! store potential, potential_dot, potential_dot_dot along interior interface of PML, shared by interior compuational domain
  integer :: nglob_interface !can be optimized
  integer, dimension(:), allocatable :: point_interface
  logical, dimension(:,:), allocatable :: PML_interior_interface
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: pml_interface_history_potential
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: pml_interface_history_potential_dot
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: pml_interface_history_potential_dot_dot

  ! buffer for I/O
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: no_backward_acoustic_buffer

  !---------------------------------------------------------------------
  ! for elastic simulation
  !---------------------------------------------------------------------

  ! number of node associated with elastic medium
  integer :: nglob_elastic

  ! number of purely elastic elements in this slice
  integer :: nspec_elastic
  integer :: nspec_elastic_b

  ! local flag if any elastic element is in this slice
  logical :: any_elastic

  ! global flag to determine if any elastic elements are present in all slices
  logical :: ELASTIC_SIMULATION

  logical, dimension(:), allocatable :: ispec_is_elastic

  ! inverse mass matrices
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rmass_inverse_elastic

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_elastic,veloc_elastic,displ_elastic

  ! PML/attenuation
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: displ_elastic_old

  ! LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: veloc_elastic_LDDRK,displ_elastic_LDDRK, &
                                                         veloc_elastic_LDDRK_temp

  ! RK
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: accel_elastic_rk,veloc_elastic_rk
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: veloc_elastic_initial_rk,displ_elastic_initial_rk

  ! variable for viscoelastic medium (also shared by solid in poroelastic-simulation)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: e1,e11,e13
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_e1,b_e11,b_e13 !for undo_attenuation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: e1_LDDRK,e11_LDDRK,e13_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: e1_initial_rk,e11_initial_rk,e13_initial_rk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: e1_force_rk,e11_force_rk,e13_force_rk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: A_newmark_nu1,B_newmark_nu1,A_newmark_nu2,B_newmark_nu2
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: dux_dxl_old,duz_dzl_old,dux_dzl_plus_duz_dxl_old
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_dux_dxl_old,b_duz_dzl_old,b_dux_dzl_plus_duz_dxl_old


  ! inverse mass matrix for viscoacoustic simulations
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rmass_inverse_e1
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: e1_acous
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: e1_LDDRK_acous
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: e1_acous_temp, dot_e1, dot_e1_old,A_newmark_e1, B_newmark_e1
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: A_newmark_e1_sf, B_newmark_e1_sf,e1_acous_sf,b_e1_acous_sf
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: e1_initial_rk_acous
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: e1_force_rk_acous,sum_forces_old,b_sum_forces_old

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_elastic_adj_coupling

  ! the variable for CPML in elastic simulation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                          rmemory_dux_dx,rmemory_duz_dx,rmemory_dux_dz,rmemory_duz_dz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                          rmemory_dux_dx_prime,rmemory_duz_dx_prime,rmemory_dux_dz_prime,rmemory_duz_dz_prime
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                          rmemory_dux_dx_LDDRK,rmemory_duz_dx_LDDRK,rmemory_dux_dz_LDDRK,rmemory_duz_dz_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_displ_elastic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_displ_elastic_LDDRK

  ! additional variables needed for CPML in viscoelastic simulation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tau_epsilon_nu1,tau_epsilon_nu2
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: kaPML_rmemory_dux_dxl,kaPML_rmemory_duz_dzl, &
                                                             muPML_rmemory_dux_dxl,muPML_rmemory_duz_dzl, &
                                                             muPML_rmemory_dux_dzl,muPML_rmemory_duz_dxl

  !for backward simulation in adjoint inversion
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    b_accel_elastic,b_veloc_elastic,b_displ_elastic,b_displ_elastic_old

  ! store potential, potential_dot, potential_dot_dot along interior interface of PML, shared by interior compuational domain
  ! for backward simulation in adjoint inversion
  ! integer :: nglob_interface !can be optimized
  ! integer, dimension(:), allocatable :: point_interface
  ! logical, dimension(:,:), allocatable :: PML_interior_interface
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: pml_interface_history_displ
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: pml_interface_history_veloc
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: pml_interface_history_accel

  ! MPI inner/outer
  integer :: nspec_inner_elastic,nspec_outer_elastic
  integer :: num_phase_ispec_elastic
  integer, dimension(:,:), allocatable :: phase_ispec_inner_elastic

  ! buffer for I/O
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: no_backward_displ_buffer,no_backward_accel_buffer
  integer :: no_backward_nframes,no_backward_iframe

  !---------------------------------------------------------------------
  ! for poroelastic simulation
  !---------------------------------------------------------------------
  integer :: nglob_poroelastic

  ! number of purely poroelastic elements in this slice
  integer :: nspec_poroelastic
  integer :: nspec_poroelastic_b

  ! local flag to determine if this slice has poroelastic elements
  logical :: any_poroelastic

  ! global flag for poroelastic simulations
  logical :: POROELASTIC_SIMULATION

  logical, dimension(:), allocatable :: ispec_is_poroelastic

  ! inverse mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic

  ! material properties of the poroelastic medium (solid phase:s and fluid phase [defined as w=phi(u_f-u_s)]: w)
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    accels_poroelastic,velocs_poroelastic,displs_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    accelw_poroelastic,velocw_poroelastic,displw_poroelastic

  ! PML/attenuation
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    displs_poroelastic_old

  ! LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    velocs_poroelastic_LDDRK,displs_poroelastic_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    velocw_poroelastic_LDDRK,displw_poroelastic_LDDRK

  ! RK
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    velocs_poroelastic_initial_rk,displs_poroelastic_initial_rk
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    velocw_poroelastic_initial_rk,displw_poroelastic_initial_rk

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    accels_poroelastic_rk,velocs_poroelastic_rk
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    accelw_poroelastic_rk,velocw_poroelastic_rk

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    accels_poroelastic_adj_coupling, accelw_poroelastic_adj_coupling

  double precision, dimension(:), allocatable :: porosity,tortuosity
  double precision, dimension(:,:), allocatable :: density,permeability

  ! for viscous attenuation in poroelastic_acoustic
  double precision :: theta_e,theta_s
  double precision :: alphaval,betaval,gammaval,thetainv

  double precision, dimension(:,:,:), allocatable :: rx_viscous,rz_viscous,viscox,viscoz
  double precision, dimension(:,:,:), allocatable :: rx_viscous_LDDRK,rz_viscous_LDDRK
  double precision, dimension(:,:,:), allocatable :: rx_viscous_initial_rk,rz_viscous_initial_rk
  double precision, dimension(:,:,:,:), allocatable :: rx_viscous_force_RK,rz_viscous_force_RK

  !for backward simulation in adjoint inversion
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    b_accels_poroelastic,b_velocs_poroelastic,b_displs_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    b_accelw_poroelastic,b_velocw_poroelastic,b_displw_poroelastic

  ! viscous damping
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_viscodampx,b_viscodampz

  ! strain
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: epsilondev_s,epsilondev_w
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_epsilondev_s,b_epsilondev_w

  ! MPI inner/outer
  integer :: nspec_inner_poroelastic,nspec_outer_poroelastic
  integer :: num_phase_ispec_poroelastic
  integer, dimension(:,:), allocatable :: phase_ispec_inner_poroelastic

  !---------------------------------------------------------------------

  ! PML parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_fsb_displ_elastic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_fsb_displ_elastic_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_sfb_potential_ddot_acoustic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_sfb_potential_ddot_acoustic_LDDRK

  ! for kernel computation
  integer :: tomo_material

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_ac,b_displ_ac,b_accel_ac

  ! elastic domain kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rho_kl, mu_kl, kappa_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: mu_k, kappa_k,rho_k
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhop_kl, beta_kl, alpha_kl
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: bulk_c_kl, bulk_beta_kl

  ! acoustic domain kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rho_ac_kl, kappa_ac_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhol_ac_global, kappal_ac_global
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhop_ac_kl, alpha_ac_kl

  double precision, dimension(:,:,:),allocatable:: rho_local,vp_local,vs_local

  ! approximate Hessians
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhorho_el_Hessian_final1, rhorho_el_Hessian_final2
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhorho_ac_Hessian_final1, rhorho_ac_Hessian_final2

  ! poro-elastic kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhot_kl, rhof_kl, sm_kl, eta_kl, mufr_kl, B_kl, &
    C_kl, M_kl, rhob_kl, rhofb_kl, phi_kl, mufrb_kl, rhobb_kl, rhofbb_kl, phib_kl, ratio_kl
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: cpI_kl, cpII_kl, cs_kl
  ! on global nodes
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhot_k, rhof_k, sm_k, eta_k, mufr_k, B_k, &
    C_k, M_k

  ! adjoint sources
  integer :: nadj_rec_local
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: adj_sourcearrays
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: source_adjoint
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: xir_store_loc,gammar_store_loc

  ! absorbing boundary
  logical :: anyabs
  ! number of absorbing elements
  integer :: nelemabs
  ! elastic contributions
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_absorb_elastic_left,b_absorb_elastic_right,b_absorb_elastic_bottom,b_absorb_elastic_top
  ! acoustic
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable ::  &
    b_absorb_acoustic_left,b_absorb_acoustic_right,b_absorb_acoustic_bottom, b_absorb_acoustic_top
  ! poroelastic solid
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_absorb_poro_s_left,b_absorb_poro_s_right,b_absorb_poro_s_bottom,b_absorb_poro_s_top
  ! poroelastic fluid
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_absorb_poro_w_left,b_absorb_poro_w_right,b_absorb_poro_w_bottom,b_absorb_poro_w_top

  integer :: nspec_left,nspec_right,nspec_bottom,nspec_top
  integer, dimension(:), allocatable :: ib_left,ib_right,ib_bottom,ib_top

  ! current time step
  integer :: it

  ! global points
  integer :: nglob,npgeo

  ! spectral elements
  integer :: nspec

  ! parameter read from parameter file
  integer :: nproc_read_from_database

  ! MPI interfaces
  integer :: ninterface
  integer :: max_interface_size
  integer, dimension(:), allocatable  :: my_neighbors
  integer, dimension(:), allocatable  :: my_nelmnts_neighbors
  integer, dimension(:,:,:), allocatable  :: my_interfaces

  integer, dimension(:), allocatable :: nibool_interfaces_ext_mesh
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_init
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh
  integer :: max_nibool_interfaces_ext_mesh

  integer, dimension(:), allocatable  :: nibool_interfaces_acoustic,nibool_interfaces_elastic,nibool_interfaces_poroelastic
  integer, dimension(:,:), allocatable  :: ibool_interfaces_acoustic,ibool_interfaces_elastic,ibool_interfaces_poroelastic

  integer :: ninterface_acoustic, ninterface_elastic,ninterface_poroelastic
  integer, dimension(:), allocatable :: inum_interfaces_acoustic, inum_interfaces_elastic, inum_interfaces_poroelastic
  integer :: max_ibool_interfaces_size_ac, max_ibool_interfaces_size_el, max_ibool_interfaces_size_po

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_send_faces_vector_ac
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_recv_faces_vector_ac
  integer, dimension(:), allocatable  :: request_send_recv_acoustic

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_send_faces_vector_el
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_recv_faces_vector_el
  integer, dimension(:), allocatable :: request_send_recv_elastic

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_send_faces_vector_pos,buffer_send_faces_vector_pow
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_recv_faces_vector_pos,buffer_recv_faces_vector_pow
  integer, dimension(:), allocatable  :: request_send_recv_poro

  ! for overlapping MPI communications with computation
  integer  :: nspec_outer, nspec_inner
  logical, dimension(:), allocatable :: ispec_is_inner

  ! to create a sorted version of the indirect addressing to reduce cache misses
  integer, dimension(:,:,:), allocatable :: copy_ibool_ori
  integer, dimension(:), allocatable :: integer_mask_ibool


  !---------------------------------------------------------------------
  ! for energy output
  !---------------------------------------------------------------------
  real(kind=CUSTOM_REAL) :: kinetic_energy,potential_energy

  ! Integrated energy field output int_0^t v^2 dt
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: integrated_kinetic_energy_field,max_kinetic_energy_field
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: integrated_potential_energy_field,max_potential_energy_field
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: kinetic_effective_duration_field,potential_effective_duration_field
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: total_integrated_energy_field,max_total_energy_field
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: total_effective_duration_field

  !---------------------------------------------------------------------
  ! for seismograms
  !---------------------------------------------------------------------
  integer :: seismo_offset, seismo_current

  ! for seismograms
  double precision, dimension(:,:), allocatable :: sisux,sisuz,siscurl

end module specfem_par

!=====================================================================

module specfem_par_noise

! parameter module for noise simulations

  use constants, only: CUSTOM_REAL

  use shared_parameters

  implicit none

  ! noise simulations:
  ! master station
  integer :: ispec_noise
  double precision :: xi_noise, gamma_noise

  ! for P_SV case, 0 for vertical along z-direction
  real(kind=CUSTOM_REAL) :: angle_noise

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: time_function_noise
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: source_array_noise
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: mask_noise

  ! The following array is used to hold snapshots of the generating
  ! wavefield or of the ensemble forward wavefield, depending on the type of
  ! noise simulation specified. In some cases, the entire generating wavefield
  ! or ensemble forward wavefield needs to be saved for all times steps. Since
  ! the disk space required to do this is usually quite large
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: surface_movie_y_or_z_noise

  ! For writing noise wavefields
  integer :: noise_output_ncol
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: noise_output_array
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: noise_output_rhokl

end module specfem_par_noise

!=====================================================================

module specfem_par_gpu

! parameter module for gpu simulations

  use constants, only: CUSTOM_REAL

  use shared_parameters

  implicit none

  ! CUDA mesh pointer to integer wrapper
  integer(kind=8) :: Mesh_pointer

  ! wavefield transfers
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: tmp_displ_2D,tmp_veloc_2D,tmp_accel_2D

  ! time steps
  real(kind=CUSTOM_REAL) :: deltatf,deltatover2f,deltatsquareover2f
  real(kind=CUSTOM_REAL) :: b_deltatf,b_deltatover2f,b_deltatsquareover2f

  ! mesh dimension
  integer :: NGLOB_AB, NSPEC_AB

  logical :: ANY_ANISOTROPY

  ! absorbing boundary
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: abs_boundary_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: abs_boundary_jacobian1Dw
  integer, dimension(:,:,:), allocatable :: abs_boundary_ij

  ! free surface
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: free_surface_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: free_surface_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: free_surface_ij
  integer, dimension(:), allocatable :: free_ac_ispec
  integer :: num_free_surface_faces

  integer, dimension(:), allocatable :: cote_abs

  ! sources
  integer :: nsources_local
  integer, dimension(:), allocatable :: ispec_selected_source_loc
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sourcearray_loc
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable ::  source_time_function_loc

  ! coupling
  integer, dimension(:), allocatable :: coupling_ac_el_ispec
  integer, dimension(:,:,:), allocatable :: coupling_ac_el_ij
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: coupling_ac_el_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: coupling_ac_el_jacobian1Dw

  ! buffers
  integer, dimension(:), allocatable  :: request_send_recv_scalar_gpu
  integer, dimension(:), allocatable  :: b_request_send_recv_scalar_gpu
  integer, dimension(:), allocatable  :: request_send_recv_vector_gpu
  integer, dimension(:), allocatable  :: b_request_send_recv_vector_gpu

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_scalar_gpu
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_buffer_send_scalar_gpu
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_recv_scalar_gpu
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_buffer_recv_scalar_gpu
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_gpu
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_gpu
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_recv_vector_gpu
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_recv_vector_gpu

end module specfem_par_gpu

!=====================================================================

module specfem_par_movie

! parameter module for noise simulations

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN

  use shared_parameters

  implicit none

  double precision, dimension(:,:), allocatable :: flagrange,xinterp,zinterp,Uxinterp,Uzinterp
  double precision, dimension(:,:), allocatable :: vector_field_display

  double precision, dimension(:,:,:), allocatable :: shape2D_display
  double precision, dimension(:,:,:,:), allocatable :: dershape2D_display

  !---------------------------------------------------------------------
  ! for color image
  !---------------------------------------------------------------------
  integer :: isnapshot_number
  integer :: nb_pixel_loc
  integer, dimension(:), allocatable :: ix_image_color_source,iy_image_color_source
  integer, dimension(:), allocatable :: ix_image_color_receiver,iy_image_color_receiver
  integer, dimension(:), allocatable :: nb_pixel_per_proc
  integer, dimension(:), allocatable :: num_pixel_loc
  integer, dimension(:,:), allocatable :: num_pixel_recv
  double precision, dimension(:), allocatable :: data_pixel_recv
  double precision, dimension(:), allocatable :: data_pixel_send

  integer :: NX_IMAGE_color,NZ_IMAGE_color
  double precision :: xmin_color_image,xmax_color_image, &
                      zmin_color_image,zmax_color_image
  integer, dimension(:,:), allocatable :: iglob_image_color,copy_iglob_image_color
  double precision, dimension(:,:), allocatable :: image_color_data
  double precision, dimension(:,:), allocatable :: image_color_vp_display

  !---------------------------------------------------------------------
  ! for plot_post
  !---------------------------------------------------------------------
  double precision :: vpImin,vpImax,vpIImin,vpIImax

  ! title of the plot
  character(len=MAX_STRING_LEN) simulation_title

  ! US letter paper or European A4
  double precision, dimension(:,:), allocatable  :: coorg_send_ps_velocity_model
  double precision, dimension(:,:), allocatable  :: coorg_recv_ps_velocity_model
  double precision, dimension(:,:), allocatable  :: RGB_send_ps_velocity_model
  double precision, dimension(:,:), allocatable  :: RGB_recv_ps_velocity_model

  double precision, dimension(:,:), allocatable  :: coorg_send_ps_element_mesh
  double precision, dimension(:,:), allocatable  :: coorg_recv_ps_element_mesh
  integer, dimension(:), allocatable  :: color_send_ps_element_mesh
  integer, dimension(:), allocatable  :: color_recv_ps_element_mesh

  double precision, dimension(:,:), allocatable  :: coorg_send_ps_abs
  double precision, dimension(:,:), allocatable  :: coorg_recv_ps_abs

  double precision, dimension(:,:), allocatable  :: coorg_send_ps_free_surface
  double precision, dimension(:,:), allocatable  :: coorg_recv_ps_free_surface

  double precision, dimension(:,:), allocatable  :: coorg_send_ps_vector_field
  double precision, dimension(:,:), allocatable  :: coorg_recv_ps_vector_field

  !---------------------------------------------------------------------
  ! for wavefield dump
  !---------------------------------------------------------------------
  logical :: this_is_the_first_time_we_dump
  logical, dimension(:), allocatable  :: mask_ibool
  logical, dimension(:), allocatable  :: mask_duplicate ! mask array for identifying duplicates between partitions
  integer, dimension(:), allocatable :: dump_recv_counts
  double precision, dimension(:,:), allocatable :: dump_recv, dump_send, dump_write, dump_gather
  logical, dimension(:), allocatable :: dump_duplicate_send, dump_duplicate_recv, dump_duplicate_gather

end module specfem_par_movie


