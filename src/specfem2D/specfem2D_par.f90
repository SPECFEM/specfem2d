
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

module constants

  include "constants.h"

end module constants

!=====================================================================

module specfem_par

! main parameter module for specfem simulations

  use constants

  implicit none

  !=====================================================================
  !input for simulation (its beginning)
  !=====================================================================

  !---------------------------------------------------------------------
  !for model description
  !---------------------------------------------------------------------
  character(len=100) :: MODEL 
  integer :: SIMULATION_TYPE  ! 1 = forward wavefield, 3 = backward and adjoint wavefields and kernels
  logical :: p_sv   ! for P-SV or SH (membrane) waves calculation
  logical :: SAVE_FORWARD ! whether or not the last frame is saved to reconstruct the forward field

  ! add a small crack (discontinuity) in the medium manually
  logical, parameter :: ADD_A_SMALL_CRACK_IN_THE_MEDIUM = .false.
  !! must be set equal to the number of spectral elements on one vertical side of the crack
  integer :: NB_POINTS_TO_ADD_TO_NPGEO = 3
  integer :: check_nb_points_to_add_to_npgeo,current_last_point,npgeo_ori,original_value,ignod
  logical :: already_found_a_crack_element

  logical :: read_external_mesh

  !---------------------------------------------------------------------
  ! 1.2 for material information
  !---------------------------------------------------------------------
  integer :: numat
  logical :: assign_external_model

  ! poroelastic and elastic coefficients
  double precision, dimension(:,:,:), allocatable :: poroelastcoef
  logical, dimension(:), allocatable :: already_shifted_velocity
  double precision, dimension(:,:,:), allocatable :: vpext,vsext,rhoext,gravityext,Nsqext
  double precision, dimension(:,:,:), allocatable :: QKappa_attenuationext,Qmu_attenuationext

  ! anisotropy parameters
  logical :: all_anisotropic
  double precision, dimension(:,:,:), allocatable :: c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext
  double precision ::  c11,c13,c15,c33,c35,c55,c12,c23,c25
  logical, dimension(:), allocatable :: anisotropic
  double precision, dimension(:,:), allocatable :: anisotropy

  ! for attenuation
  logical ATTENUATION_VISCOELASTIC_SOLID, ATTENUATION_PORO_FLUID_PART
  integer :: N_SLS
  double precision, dimension(:), allocatable  :: QKappa_attenuation
  double precision, dimension(:), allocatable  :: Qmu_attenuation
  double precision  :: f0_attenuation
  logical :: READ_VELOCITIES_AT_f0
  integer nspec_allocate
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2
  real(kind=CUSTOM_REAL), dimension(:,:,:) , allocatable :: Mu_nu1,Mu_nu2

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: tau_epsilon_nu1,tau_epsilon_nu2, &
                                                       inv_tau_sigma_nu1_sent,inv_tau_sigma_nu2_sent,&
                                                       phi_nu1_sent,phi_nu2_sent
  real(kind=CUSTOM_REAL) :: Mu_nu1_sent,Mu_nu2_sent

  !---------------------------------------------------------------------
  !for boundary condition (physical BC or artificial BC)
  !---------------------------------------------------------------------
  logical :: anyabs_glob

  !PML
  logical :: PML_BOUNDARY_CONDITIONS,ROTATE_PML_ACTIVATE
  double precision :: ROTATE_PML_ANGLE
  integer :: nspec_PML,NELEM_PML_THICKNESS
  logical, dimension(:), allocatable :: is_PML
  integer, dimension(:), allocatable :: region_CPML
  integer, dimension(:), allocatable :: spec_to_PML
  logical, dimension(:,:), allocatable :: which_PML_elem

  double precision, dimension(:,:,:), allocatable :: &
                    K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  ! stacey BC
  logical :: STACEY_BOUNDARY_CONDITIONS
  logical, dimension(:,:), allocatable  :: codeabs
  integer, dimension(:), allocatable  :: typeabs
  ! for detection of corner element on absorbing boundary
  logical, dimension(:,:), allocatable  :: codeabs_corner
  ! add spring to Stacey absorbing boundary condition
  logical :: ADD_SPRING_TO_STACEY
  double precision :: x_center_spring,z_center_spring
  double precision :: xmin,xmax,zmin,zmax
  double precision :: xmin_local,xmax_local,zmin_local,zmax_local

  ! for horizontal periodic conditions
  logical :: ADD_PERIODIC_CONDITIONS
  ! horizontal periodicity distance for periodic conditions
  double precision :: PERIODIC_HORIZ_DIST
  logical, dimension(:), allocatable :: this_ibool_is_a_periodic_edge
  double precision :: xmaxval,xminval,ymaxval,yminval,xtol,xtypdist
  integer :: counter

  ! fluid/solid interface
  integer :: num_fluid_solid_edges
  logical :: coupled_acoustic_elastic,any_fluid_solid_edges
  integer, dimension(NEDGES) :: i_begin,j_begin,i_end,j_end
  integer, dimension(NGLLX,NEDGES) :: ivalue,jvalue,ivalue_inverse,jvalue_inverse
  integer, dimension(:), allocatable :: fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
                                        fluid_solid_elastic_ispec,fluid_solid_elastic_iedge

  ! fluid/porous interface
  integer :: num_fluid_poro_edges
  logical :: coupled_acoustic_poro,any_fluid_poro_edges
  integer, dimension(:), allocatable :: fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge, &
                                        fluid_poro_poroelastic_ispec,fluid_poro_poroelastic_iedge

  ! solid/porous interface
  logical :: coupled_elastic_poro, any_solid_poro_edges
  integer, dimension(:), allocatable :: solid_poro_elastic_ispec,solid_poro_elastic_iedge, &
                                        solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge

  ! solid/porous interface
  integer :: num_solid_poro_edges
  integer, dimension(:), allocatable :: icount
  integer, dimension(:), allocatable :: ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro,&
            iend_edge3_poro,ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro

  !---------------------------------------------------------------------
  !for source-receiver information
  !---------------------------------------------------------------------
  ! source description
  integer NSOURCES
  integer, dimension(:), allocatable :: source_type,time_function_type
  double precision, dimension(:), allocatable :: x_source,z_source,xi_source,gamma_source,&
                                                 Mxx,Mzz,Mxz,f0,tshift_src,factor,anglesource
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: sourcearray
  double precision :: t0
  integer, dimension(:), allocatable :: ispec_selected_source,iglob_source,&
                                        is_proc_source,nb_proc_source
  double precision, dimension(:), allocatable :: aval
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: source_time_function
  double precision, external :: netlib_specfun_erf
  ! use this t0 as earliest starting time rather than the automatically calculated one
  ! (must be positive and bigger than the automatically one to be effective;
  !  simulation will start at t = - t0)
  double precision :: USER_T0

  ! for absorbing and acoustic free surface conditions
  integer :: ispec_acoustic_surface,inum
  real(kind=CUSTOM_REAL) :: nx,nz,weight,xxi,zgamma
  !acoustic free surface
  integer, dimension(:,:), allocatable :: acoustic_surface
  integer, dimension(:,:), allocatable :: acoustic_edges
  logical :: any_acoustic_edges
  integer  :: ixmin, ixmax, izmin, izmax

  ! perform a forcing of an acoustic medium with a rigid boundary
  logical :: ACOUSTIC_FORCING
  integer :: nelem_acforcing,nelem_acoustic_surface
  logical, dimension(:,:), allocatable  :: codeacforcing
  integer, dimension(:), allocatable  :: typeacforcing
  integer, dimension(:), allocatable :: numacforcing, &
     ibegin_edge1_acforcing,iend_edge1_acforcing,ibegin_edge3_acforcing,iend_edge3_acforcing, &
     ibegin_edge4_acforcing,iend_edge4_acforcing,ibegin_edge2_acforcing,iend_edge2_acforcing
  integer :: nspec_left_acforcing,nspec_right_acforcing,nspec_bottom_acforcing,nspec_top_acforcing
  integer, dimension(:), allocatable :: ib_left_acforcing,ib_right_acforcing,ib_bottom_acforcing,ib_top_acforcing

  ! for plane wave incidence
  ! to compute analytical initial plane wave field
  logical initialfield,add_Bielak_conditions
  double precision :: anglesource_refl, c_inc, c_refl, cploc, csloc
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

  integer  :: nrec,nrecloc
  double precision :: anglerec,xirec,gammarec
  integer, dimension(:), allocatable :: recloc, which_proc_receiver
  integer, dimension(:), allocatable :: ispec_selected_rec
  double precision, dimension(:), allocatable :: xi_receiver,gamma_receiver,st_xval,st_zval

  ! tangential detection for source or receivers
  double precision, dimension(:), allocatable :: anglerec_irec
  double precision, dimension(:), allocatable :: cosrot_irec, sinrot_irec
  double precision, dimension(:), allocatable :: x_final_receiver, z_final_receiver
  logical :: force_normal_to_surface,rec_normal_to_surface

  integer, dimension(:), allocatable :: source_courbe_eros

  integer  :: nnodes_tangential_curve
  double precision, dimension(:,:), allocatable  :: nodes_tangential_curve
  logical  :: any_tangential_curve

  integer  :: n1_tangential_detection_curve
  integer, dimension(4)  :: n_tangential_detection_curve
  integer, dimension(:), allocatable  :: rec_tangential_detection_curve
  double precision :: distmin, dist_current, anglesource_recv
  double precision, dimension(:), allocatable :: dist_tangential_detection_curve
  double precision :: x_final_receiver_dummy, z_final_receiver_dummy

  !---------------------------------------------------------------------
  ! for SEM discretization of the model
  !---------------------------------------------------------------------
  ! for Lagrange interpolants
  double precision, external :: hgll
  ! Gauss-Lobatto-Legendre points and weights
  double precision, dimension(NGLLX) :: xigll
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLZ) :: zigll
  real(kind=CUSTOM_REAL), dimension(NGLLZ) :: wzgll
  ! derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

  double precision, dimension(:,:,:), allocatable :: shape2D,shape2D_display
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: xix,xiz,gammax,gammaz,jacobian

  double precision, dimension(:,:,:,:), allocatable :: dershape2D,dershape2D_display

  integer, dimension(:,:,:), allocatable :: ibool,ibool_outer,ibool_inner
  integer, dimension(:,:), allocatable  :: knods
  integer, dimension(:), allocatable :: kmato,numabs, &
     ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3,ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2

  ! Lagrange interpolators at receivers
  double precision, dimension(:), allocatable :: hxir,hgammar,hpxir,hpgammar
  double precision, dimension(:,:), allocatable :: hxir_store,hgammar_store

  ! Lagrange interpolators at sources
  double precision, dimension(:), allocatable :: hxis,hgammas,hpxis,hpgammas
  double precision, dimension(:,:), allocatable :: hxis_store,hgammas_store

  !---------------------------------------------------------------------
  ! AXISYM parameters
  !---------------------------------------------------------------------
  logical :: AXISYM ! .true. if we are performing a 2.5D simulation
  ! Number of elements on the symmetry axis
  integer :: nelem_on_the_axis,nelem_on_the_axis_total
  ! Flag to know if an element is on the axis
  logical, dimension(:), allocatable :: is_on_the_axis
  integer, dimension(:), allocatable  ::ispec_of_axial_elements
  ! Gauss-Lobatto-Jacobi points and weights
  double precision, dimension(NGLJ) :: xiglj
  real(kind=CUSTOM_REAL), dimension(NGLJ) :: wxglj
  ! derivatives of GLJ polynomials
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLJ) :: hprimeBar_xx,hprimeBarwglj_xx
  ! Shape functions (and their derivatives) evaluated at the GLJ points
  double precision, dimension(:,:), allocatable :: flagrange_GLJ

  !---------------------------------------------------------------------
  ! for the check of mesh
  !---------------------------------------------------------------------
  integer :: UPPER_LIMIT_DISPLAY

  !---------------------------------------------------------------------
  ! for parallel simulation
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! for time discretization
  !---------------------------------------------------------------------
  ! value of time_stepping_scheme to decide which time scheme will be used
  ! 1 = Newmark (2nd order), 2 = LDDRK4-6 (4th-order 6-stage low storage Runge-Kutta)
  ! 3 = classical 4th-order 4-stage Runge-Kutta
  integer :: time_stepping_scheme

  ! coefficients of the explicit Newmark time scheme
  integer NSTEP
  double precision :: deltatover2,deltatsquareover2,timeval
  double precision :: deltat

  ! for backward simulation in adjoint inversion
  double precision :: b_deltatover2,b_deltatsquareover2,b_deltat ! coefficients of the explicit Newmark time scheme

  ! for LDDRK46
  integer :: i_stage,stage_time_scheme
  real(kind=CUSTOM_REAL), dimension(Nstages):: alpha_LDDRK,beta_LDDRK,c_LDDRK

  ! parameters used in LDDRK scheme, from equation (2) of
  ! Berland, J., Bogey, C., & Bailly, C.
  ! Low-dissipation and low-dispersion fourth-order Runge-Kutta algorithm, Computers & Fluids, 35(10), 1459-1463.
  Data alpha_LDDRK /0.0_CUSTOM_REAL,-0.737101392796_CUSTOM_REAL, &
                    -1.634740794341_CUSTOM_REAL,-0.744739003780_CUSTOM_REAL, &
                    -1.469897351522_CUSTOM_REAL,-2.813971388035_CUSTOM_REAL/

  Data beta_LDDRK /0.032918605146_CUSTOM_REAL,0.823256998200_CUSTOM_REAL,&
                   0.381530948900_CUSTOM_REAL,0.200092213184_CUSTOM_REAL,&
                   1.718581042715_CUSTOM_REAL,0.27_CUSTOM_REAL/

  Data c_LDDRK /0.0_CUSTOM_REAL,0.032918605146_CUSTOM_REAL,&
                0.249351723343_CUSTOM_REAL,0.466911705055_CUSTOM_REAL,&
                0.582030414044_CUSTOM_REAL,0.847252983783_CUSTOM_REAL/

  ! for rk44
  double precision :: weight_rk
  !=====================================================================
  ! input for simulation (its end)
  !=====================================================================


  !=====================================================================
  ! for simulation (its beginning)
  !=====================================================================
  ! to help locate elements with a negative Jacobian using OpenDX
  logical :: found_a_negative_jacobian

  ! to count the number of degrees of freedom
  integer :: count_nspec_acoustic,count_nspec_acoustic_total,nspec_total,nglob_total,nb_acoustic_DOFs,nb_elastic_DOFs
  double precision :: ratio_1DOF,ratio_2DOFs

  ! to determine date and time at which the run will finish
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  integer :: year,mon,day,hr,minutes,timestamp
  double precision :: timestamp_seconds_start

  !---------------------------------------------------------------------
  !global varable shared by acoustic/elastic/poroelastic simulation
  !---------------------------------------------------------------------
  double precision, dimension(:,:), allocatable :: &
    coord, flagrange,xinterp,zinterp,Uxinterp,Uzinterp,vector_field_display

  double precision, dimension(:,:), allocatable :: coorg

  !---------------------------------------------------------------------
  !for acoustic simulation
  !---------------------------------------------------------------------
  ! for acoustic medium
  logical :: any_acoustic,any_acoustic_glob
  integer :: nglob_acoustic
  integer :: nglob_gravitoacoustic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic,potential_acoustic_old
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_dot_acoustic_LDDRK, potential_acoustic_LDDRK
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_dot_acoustic_temp
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic_init_rk, potential_dot_acoustic_init_rk
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: potential_dot_dot_acoustic_rk, potential_dot_acoustic_rk
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic_adj_coupling

  ! for gravitoacoustic medium
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic,potential_gravitoacoustic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    potential_dot_dot_gravito,potential_dot_gravito,potential_gravito

  ! for acoustic and gravitoacoustic detection
  logical, dimension(:), allocatable :: acoustic,gravitoacoustic
  logical :: any_gravitoacoustic,any_gravitoacoustic_glob

  ! inverse mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_acoustic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_gravitoacoustic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_gravito

  ! the variable for PML
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_potential_acoustic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                          rmemory_acoustic_dux_dx,rmemory_acoustic_dux_dz

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_potential_acoustic_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                          rmemory_acoustic_dux_dx_LDDRK,rmemory_acoustic_dux_dz_LDDRK

  ! for backward simulation in adjoint inversion
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic,b_potential_acoustic_old

  ! store potential, potential_dot, potential_dot_dot along interior interface of PML, shared by interior compuational domain
  integer :: nglob_interface !can be optimized 
  integer, dimension(:), allocatable :: point_interface
  logical, dimension(:,:), allocatable :: PML_interior_interface
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: pml_interface_history_potential
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: pml_interface_history_potential_dot
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: pml_interface_history_potential_dot_dot

  !---------------------------------------------------------------------
  !for by elastic simulation
  !---------------------------------------------------------------------
  ! number of node associated with elastic medium
  logical :: any_elastic,any_elastic_glob
  integer :: nglob_elastic
  logical, dimension(:), allocatable :: elastic

  ! inverse mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_elastic_one
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_elastic_three

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: veloc_elastic_LDDRK,displ_elastic_LDDRK,&
                                                         veloc_elastic_LDDRK_temp

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: accel_elastic_rk,veloc_elastic_rk
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: veloc_elastic_initial_rk,displ_elastic_initial_rk

  ! variable for viscoelastic medium (also shared by solid in poroelastic-simulation)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: e1,e11,e13
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: e1_LDDRK,e11_LDDRK,e13_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: e1_initial_rk,e11_initial_rk,e13_initial_rk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: e1_force_rk,e11_force_rk,e13_force_rk

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_elastic_adj_coupling,accel_elastic_adj_coupling2

  ! the variable for PML
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                          rmemory_dux_dx,rmemory_duz_dx,rmemory_dux_dz,rmemory_duz_dz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                          rmemory_dux_dx_prime,rmemory_duz_dx_prime,rmemory_dux_dz_prime,rmemory_duz_dz_prime
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                          rmemory_dux_dx_LDDRK,rmemory_duz_dx_LDDRK,rmemory_dux_dz_LDDRK,rmemory_duz_dz_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_displ_elastic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_displ_elastic_LDDRK

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
  !---------------------------------------------------------------------
  !for by poroelastic simulation
  !---------------------------------------------------------------------
  logical :: any_poroelastic,any_poroelastic_glob
  integer :: nglob_poroelastic
  logical, dimension(:), allocatable :: poroelastic

  ! inverse mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic

  ! material properties of the poroelastic medium (solid phase:s and fluid phase [defined as w=phi(u_f-u_s)]: w)
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    accels_poroelastic,velocs_poroelastic,displs_poroelastic, displs_poroelastic_old
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    velocs_poroelastic_LDDRK,displs_poroelastic_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    velocs_poroelastic_initial_rk,displs_poroelastic_initial_rk
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    accels_poroelastic_rk,velocs_poroelastic_rk
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    accelw_poroelastic,velocw_poroelastic,displw_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    velocw_poroelastic_LDDRK,displw_poroelastic_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    velocw_poroelastic_initial_rk,displw_poroelastic_initial_rk
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    accelw_poroelastic_rk,velocw_poroelastic_rk

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    accels_poroelastic_adj_coupling, accelw_poroelastic_adj_coupling
  double precision, dimension(:), allocatable :: porosity,tortuosity
  double precision, dimension(:,:), allocatable :: density,permeability

  ! for viscous attenuation in poroelastic_acoustic
  double precision :: theta_e,theta_s
  double precision :: Q0,freq0
  double precision :: alphaval,betaval,gammaval,thetainv

  double precision, dimension(NGLLX,NGLLZ) :: viscox_loc,viscoz_loc
  double precision :: Sn,Snp1,etal_f
  double precision, dimension(3):: bl_unrelaxed_elastic
  double precision :: permlxx,permlxz,permlzz,invpermlxx,invpermlxz,invpermlzz,detk
  ! for shifting of velocities if needed in the case of viscoelasticity
  double precision :: vp,vs,rho,mu,lambda

  double precision, dimension(:,:,:), allocatable :: rx_viscous,rz_viscous,viscox,viscoz
  double precision, dimension(:,:,:), allocatable :: rx_viscous_LDDRK,rz_viscous_LDDRK
  double precision, dimension(:,:,:), allocatable :: rx_viscous_initial_rk,rz_viscous_initial_rk
  double precision, dimension(:,:,:,:), allocatable :: rx_viscous_force_RK,rz_viscous_force_RK

  !for backward simulation in adjoint inversion
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    b_accels_poroelastic,b_velocs_poroelastic,b_displs_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    b_accelw_poroelastic,b_velocw_poroelastic,b_displw_poroelastic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: b_viscodampx,b_viscodampz


  !---------------------------------------------------------------------
  !for fluid/solid coupling 
  !---------------------------------------------------------------------
  integer :: ispec_acoustic,ispec_elastic,iedge_acoustic,iedge_elastic,ipoin1D,iglob2
  real(kind=CUSTOM_REAL) :: displ_x,displ_z,displ_n,displw_x,displw_z,zxi,xgamma,jacobian1D,pressure
  ! PML parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_fsb_displ_elastic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_fsb_displ_elastic_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_sfb_potential_ddot_acoustic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_sfb_potential_ddot_acoustic_LDDRK

  !for adjoint
  real(kind=CUSTOM_REAL) :: b_displ_x,b_displ_z,b_displw_x,b_displw_z,b_pressure

  !---------------------------------------------------------------------
  !for fluid/porous coupling 
  !---------------------------------------------------------------------
  integer :: iedge_poroelastic
  double precision :: mul_G,lambdal_G,lambdalplus2mul_G
  double precision :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  double precision :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  double precision :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  double precision :: dwx_dxl,dwz_dxl,dwx_dzl,dwz_dzl

  !for adjoint
  double precision :: b_dux_dxi,b_dux_dgamma,b_duz_dxi,b_duz_dgamma
  double precision :: b_dwx_dxi,b_dwx_dgamma,b_dwz_dxi,b_dwz_dgamma
  double precision :: b_dux_dxl,b_duz_dxl,b_dux_dzl,b_duz_dzl
  double precision :: b_dwx_dxl,b_dwz_dxl,b_dwx_dzl,b_dwz_dzl

  !---------------------------------------------------------------------
  !for solid/porous coupling 
  !---------------------------------------------------------------------
  integer :: ispec_poroelastic,ii2,jj2
  double precision :: sigma_xx,sigma_xz,sigma_zz,sigmap

  !for adjoint
  double precision :: b_sigma_xx,b_sigma_xz,b_sigma_zz,b_sigmap


  ! for kernel compuation
  character(len=100) TOMOGRAPHY_FILE
  integer :: tomo_material
  logical :: save_ASCII_kernels
  integer reclen
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_ac,b_displ_ac,b_accel_ac
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rho_kl, mu_kl, kappa_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhol_global, mul_global, kappal_global
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: mu_k, kappa_k,rho_k
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhop_kl, beta_kl, alpha_kl
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rho_ac_kl, kappa_ac_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhol_ac_global, kappal_ac_global
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhop_ac_kl, alpha_ac_kl

  double precision, dimension(:,:,:),allocatable:: rho_local,vp_local,vs_local 
 
  !!!! hessian
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhorho_el_hessian_final1, rhorho_el_hessian_final2
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhorho_el_hessian_temp1, rhorho_el_hessian_temp2
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhorho_ac_hessian_final1, rhorho_ac_hessian_final2

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhot_kl, rhof_kl, sm_kl, eta_kl, mufr_kl, B_kl, &
    C_kl, M_kl, rhob_kl, rhofb_kl, phi_kl, Bb_kl, Cb_kl, Mb_kl, mufrb_kl, &
    rhobb_kl, rhofbb_kl, phib_kl, cpI_kl, cpII_kl, cs_kl, ratio_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhot_k, rhof_k, sm_k, eta_k, mufr_k, B_k, &
    C_k, M_k
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: phil_global,etal_f_global, &
    rhol_s_global,rhol_f_global,rhol_bar_global, &
    tortl_global,mulfr_global
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: permlxx_global,permlxz_global,permlzz_global
  character(len=150) :: adj_source_file
  integer :: irec_local,nadj_rec_local
  double precision :: xx,zz,rholb,tempx1l,tempx2l,b_tempx1l,b_tempx2l,bb_tempx1l,bb_tempx2l
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: adj_sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: adj_sourcearrays
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_absorb_elastic_left,b_absorb_poro_s_left,b_absorb_poro_w_left
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_absorb_elastic_right,b_absorb_poro_s_right,b_absorb_poro_w_right
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_absorb_elastic_bottom,b_absorb_poro_s_bottom,b_absorb_poro_w_bottom
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_absorb_elastic_top,b_absorb_poro_s_top,b_absorb_poro_w_top
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable ::  &
    b_absorb_acoustic_left,b_absorb_acoustic_right,&
    b_absorb_acoustic_bottom, b_absorb_acoustic_top
  integer :: nspec_left,nspec_right,nspec_bottom,nspec_top
  integer, dimension(:), allocatable :: ib_left,ib_right,ib_bottom,ib_top
  real(kind=CUSTOM_REAL),  dimension(:,:,:), allocatable :: source_adjointe
  real(kind=CUSTOM_REAL),  dimension(:,:), allocatable :: xir_store_loc, gammar_store_loc

  !---------------------------------------------------------------------
  !local varable used in the unclean part of code in iterate_time.F90
  !prepare_timerun_body.F90 et al
  !---------------------------------------------------------------------
  logical :: anyabs
  double precision :: dxd,dyd,dzd,dcurld,valux,valuy,valuz,valcurl,hlagrange,rhol,xi,gamma,x,z
  double precision :: gravityl,Nsql,hp1,hp2

  real(kind=CUSTOM_REAL) :: kinetic_energy,potential_energy,kinetic_energy_total,potential_energy_total
  double precision :: vpImin,vpImax,vpIImin,vpIImax
  integer :: iglobzero,ios
  integer :: it,id,n,nglob,npgeo
  character(len=150) dummystring
  ! material properties of the elastic medium
  double precision :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic
  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl

  ! to evaluate cpI, cpII, and cs, and rI (poroelastic medium)
  double precision :: rhol_s,rhol_f,rhol_bar,phil,tortl
  double precision :: mul_s,kappal_s
  double precision :: kappal_f
  double precision :: mul_fr,kappal_fr
  double precision :: D_biot,H_biot,C_biot,M_biot,B_biot,cpIsquare,cpIIsquare,cssquare
  real(kind=CUSTOM_REAL) :: ratio,dd1

  integer :: ngnod,nspec,pointsdisp, nelemabs

  ! for MPI and partitioning
  integer  :: ier
  integer  :: myrank,nproc,nproc_read_from_database
  character(len=150) :: inputname,outputname,outputname2
  integer  :: ninterface
  integer  :: max_interface_size
  integer, dimension(:), allocatable  :: my_neighbours
  integer, dimension(:), allocatable  :: my_nelmnts_neighbours
  integer, dimension(:,:,:), allocatable  :: my_interfaces
  integer, dimension(:,:), allocatable  :: ibool_interfaces_acoustic,ibool_interfaces_elastic,ibool_interfaces_poroelastic
  integer, dimension(:), allocatable  :: nibool_interfaces_acoustic,nibool_interfaces_elastic,nibool_interfaces_poroelastic
  integer, dimension(:), allocatable :: nibool_interfaces_ext_mesh
  integer, dimension(:,:), allocatable :: ibool_interfaces_ext_mesh_init, ibool_interfaces_ext_mesh
  integer :: max_nibool_interfaces_ext_mesh

  integer  :: ninterface_acoustic, ninterface_elastic,ninterface_poroelastic
  integer, dimension(:), allocatable  :: inum_interfaces_acoustic, inum_interfaces_elastic, inum_interfaces_poroelastic

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_send_faces_vector_ac
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_recv_faces_vector_ac
  integer, dimension(:), allocatable  :: tab_requests_send_recv_acoustic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_send_faces_vector_el
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_recv_faces_vector_el
  integer, dimension(:), allocatable  :: tab_requests_send_recv_elastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_send_faces_vector_pos,buffer_send_faces_vector_pow
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_recv_faces_vector_pos,buffer_recv_faces_vector_pow
  integer, dimension(:), allocatable  :: tab_requests_send_recv_poro
  integer :: max_ibool_interfaces_size_ac, max_ibool_interfaces_size_el, max_ibool_interfaces_size_po
  integer :: iproc

  ! for overlapping MPI communications with computation
  integer  :: nspec_outer, nspec_inner, num_ispec_outer, num_ispec_inner
  integer, dimension(:), allocatable  :: ispec_outer_to_glob, ispec_inner_to_glob
  logical, dimension(:), allocatable  :: mask_ispec_inner_outer

  ! inner/outer elements in the case of an MPI simulation
  integer :: nglob_outer,nglob_inner 

  ! to create a sorted version of the indirect addressing to reduce cache misses
  integer, dimension(:,:,:), allocatable :: copy_ibool_ori
  integer, dimension(:), allocatable :: integer_mask_ibool

  !=====================================================================
  ! for simulation (its end)
  !=====================================================================


  !=====================================================================
  ! output for simulation (the beginning)
  !=====================================================================
  !---------------------------------------------------------------------
  ! for information of the stability behavior during the simulation
  !---------------------------------------------------------------------
  integer :: NSTEP_BETWEEN_OUTPUT_INFO

  !---------------------------------------------------------------------
  ! for energy output
  !---------------------------------------------------------------------
  logical :: output_energy

  !---------------------------------------------------------------------
  ! for seismograms
  !---------------------------------------------------------------------
  integer :: seismotype,NSTEP_BETWEEN_OUTPUT_SEISMOS
  integer :: seismo_offset, seismo_current, subsamp_seismos

  logical :: USE_TRICK_FOR_BETTER_PRESSURE
  logical :: save_ASCII_seismograms,save_binary_seismograms_single,save_binary_seismograms_double

  ! output seismograms in Seismic Unix format (adjoint traces will be read in the same format)
  logical :: SU_FORMAT
  !<SU_FORMAT
  integer(kind=4) :: r4head(60)
  character(len=512) :: filename
  real(kind=4),dimension(:,:),allocatable :: adj_src_s
  integer(kind=2) :: header2(2)
  !>SU_FORMAT

  ! for seismograms
  double precision, dimension(:,:), allocatable :: sisux,sisuz,siscurl
  ! vector field in an element
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLX) :: vector_field_element
  ! pressure in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: pressure_element
  ! curl in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: curl_element
  !---------------------------------------------------------------------
  ! for color image
  !---------------------------------------------------------------------
  integer :: colors !also used in plotpost
  double precision :: cutsnaps !also used in plotpost
  logical :: output_color_image,DRAW_SOURCES_AND_RECEIVERS
  integer :: NSTEP_BETWEEN_OUTPUT_IMAGES,imagetype_JPEG
  integer :: isnapshot_number = 0  !remember which image are going to produce
  integer  :: nb_pixel_loc
  integer, dimension(:), allocatable :: ix_image_color_source,iy_image_color_source
  integer, dimension(:), allocatable :: ix_image_color_receiver,iy_image_color_receiver
  integer, dimension(:), allocatable :: nb_pixel_per_proc
  integer, dimension(:), allocatable :: num_pixel_loc
  integer, dimension(:,:), allocatable :: num_pixel_recv
  double precision, dimension(:), allocatable :: data_pixel_recv
  double precision, dimension(:), allocatable :: data_pixel_send

  ! factor to subsample color images output by the code (useful for very large models)
  double precision :: factor_subsample_image
  ! by default the code normalizes each image independently to its maximum; use this option to use the global maximum below instead
  logical :: USE_CONSTANT_MAX_AMPLITUDE
  ! constant maximum amplitude to use for all color images if the USE_CONSTANT_MAX_AMPLITUDE option is true
  double precision :: CONSTANT_MAX_AMPLITUDE_TO_USE
  ! use snapshot number in the file name of JPG color snapshots instead of the time step
  logical :: USE_SNAPSHOT_NUMBER_IN_FILENAME
  ! display acoustic layers as constant blue, because they likely correspond to water in the case of ocean acoustics
  ! or in the case of offshore oil industry experiments.
  ! (if off, display them as greyscale, as for elastic or poroelastic elements,
  !  for instance for acoustic-only oil industry models of solid media)
  logical :: DRAW_WATER_IN_BLUE
  ! non linear display to enhance small amplitudes in color images
  double precision :: POWER_DISPLAY_COLOR

  integer :: NX_IMAGE_color,NZ_IMAGE_color
  double precision :: xmin_color_image,xmax_color_image, &
                      zmin_color_image,zmax_color_image
  integer, dimension(:,:), allocatable :: iglob_image_color,copy_iglob_image_color
  double precision, dimension(:,:), allocatable :: image_color_data
  double precision, dimension(:,:), allocatable :: image_color_vp_display

  !---------------------------------------------------------------------
  ! for plotpost
  !---------------------------------------------------------------------
  integer :: subsamp_postscript,imagetype_postscript
  integer :: numbers
  double precision :: sizemax_arrows
  logical :: output_postscript_snapshot,US_LETTER,plot_lowerleft_corner_only
  logical :: interpol,meshvect,modelvect,boundvect
  ! title of the plot
  character(len=60) simulation_title
  ! US letter paper or European A4
  integer :: d1_coorg_send_ps_velocity_model,d2_coorg_send_ps_velocity_model, &
             d1_coorg_recv_ps_velocity_model,d2_coorg_recv_ps_velocity_model, &
             d1_RGB_send_ps_velocity_model,d2_RGB_send_ps_velocity_model, &
             d1_RGB_recv_ps_velocity_model,d2_RGB_recv_ps_velocity_model
  double precision, dimension(:,:), allocatable  :: coorg_send_ps_velocity_model
  double precision, dimension(:,:), allocatable  :: coorg_recv_ps_velocity_model
  double precision, dimension(:,:), allocatable  :: RGB_send_ps_velocity_model
  double precision, dimension(:,:), allocatable  :: RGB_recv_ps_velocity_model
  integer :: d1_coorg_send_ps_element_mesh,d2_coorg_send_ps_element_mesh, &
             d1_coorg_recv_ps_element_mesh,d2_coorg_recv_ps_element_mesh, &
             d1_color_send_ps_element_mesh, &
             d1_color_recv_ps_element_mesh
  double precision, dimension(:,:), allocatable  :: coorg_send_ps_element_mesh
  double precision, dimension(:,:), allocatable  :: coorg_recv_ps_element_mesh
  integer, dimension(:), allocatable  :: color_send_ps_element_mesh
  integer, dimension(:), allocatable  :: color_recv_ps_element_mesh
  integer :: d1_coorg_send_ps_abs, d2_coorg_send_ps_abs, &
             d1_coorg_recv_ps_abs, d2_coorg_recv_ps_abs
  double precision, dimension(:,:), allocatable  :: coorg_send_ps_abs
  double precision, dimension(:,:), allocatable  :: coorg_recv_ps_abs
  integer :: d1_coorg_send_ps_free_surface, d2_coorg_send_ps_free_surface, &
             d1_coorg_recv_ps_free_surface, d2_coorg_recv_ps_free_surface
  double precision, dimension(:,:), allocatable  :: coorg_send_ps_free_surface
  double precision, dimension(:,:), allocatable  :: coorg_recv_ps_free_surface
  integer :: d1_coorg_send_ps_vector_field, d2_coorg_send_ps_vector_field, &
             d1_coorg_recv_ps_vector_field, d2_coorg_recv_ps_vector_field
  double precision, dimension(:,:), allocatable  :: coorg_send_ps_vector_field
  double precision, dimension(:,:), allocatable  :: coorg_recv_ps_vector_field

  !---------------------------------------------------------------------
  ! for wavefield damp
  !---------------------------------------------------------------------
  integer :: NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS
  integer :: imagetype_wavefield_dumps
  ! name of wavefield snapshot file
  character(len=150) :: wavefield_file
  logical :: output_wavefield_dumps,use_binary_for_wavefield_dumps
  integer :: icounter,nb_of_values_to_save ! icounter is local variable iterate_time.F90
  logical :: this_is_the_first_time_we_dump
  logical, dimension(:), allocatable  :: mask_ibool,mask_ibool_pml


  !---------------------------------------------------------------------
  ! for wavefield snapshot file
  !---------------------------------------------------------------------
  logical :: output_grid_ASCII,output_grid_Gnuplot
  !=====================================================================
  ! output for simulation (end)
  !=====================================================================


!<NOISE_TOMOGRAPHY
  ! NOISE_TOMOGRAPHY = 0 - turn noise tomography subroutines off; setting
  ! NOISE_TOMOGRAPHY equal to 0, in other words, results in an earthquake
  ! simulation rather than a noise simulation

  ! NOISE_TOMOGRAPHY = 1 - compute "generating" wavefield and store the result;
  ! this stored wavefield is then used to compute the "ensemble forward"
  ! wavefield in the next noise simulation

  ! NOISE_TOMOGRAPHY = 2 - compute "ensemble forward" wavefield and store the
  ! result; if an adjoint simulation is planned, users need to store
  ! seismograms (actually, cross-correlograms) for later processing

  ! NOISE_TOMOGRAPHY = 3 - carry out adjoint simulation; users need to supply
  ! adjoint sources constructed from cross-correlograms computed during the
  ! "ensemble forward" step


  ! For an explanation of terms and concepts in noise tomography, see "Tromp et
  ! al., 2011, Noise Cross-Correlation Sensitivity Kernels, Geophysical Journal
  ! International"

  integer :: NOISE_TOMOGRAPHY
  integer :: irec_master, ispec_noise
  double precision :: xi_noise, gamma_noise, angle_noise
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: time_function_noise
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: source_array_noise
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: mask_noise

  ! The following three arrays are used to hold snapshots of the generating
  ! wavefield or of the ensemble forward wavefield, depending on the type of
  ! noise simulation specified. In some cases, the entire generating wavefield
  ! or ensemble forward wavefield needs to be saved for all times steps. Since
  ! the disk space required to do this is usually quite large, separate arrays
  ! are used for x,y,z to avoid having empty dimensions (one empty dimension in
  ! the case of P-SV or two empty dimensions in the case of SH).
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    surface_movie_x_noise, surface_movie_y_noise, surface_movie_z_noise

  ! For writing noise wavefields
  logical :: output_wavefields_noise = .true. ! this is output only in the case of noise tomography
  logical :: ex, od
  integer :: noise_output_ncol
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: noise_output_array
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: noise_output_rhokl
  character(len=512) :: noise_output_file

  ! For noise tomography only - specify whether to reconstruct ensemble forward
  ! wavefield by saving everywhere or by saving only at the boundaries (the
  ! latter usually much faster but prone to artefacts)
  logical :: save_everywhere = .false.
!<NOISE_TOMOGRAPHY

  !---------------------------------------------------------------------
  !global varable particular for computation with GPU
  !---------------------------------------------------------------------
  ! Global GPU toggle. Set in Par_file
  logical :: GPU_MODE

  ! CUDA mesh pointer<->integer wrapper
  integer(kind=8) :: Mesh_pointer
  integer :: ncuda_devices,ncuda_devices_min,ncuda_devices_max
  logical, dimension(:), allocatable :: ispec_is_inner
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: displ_2D,veloc_2D,accel_2D
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_displ_2D,b_veloc_2D,b_accel_2D
  integer :: NGLOB_AB, NSPEC_AB
  real(kind=CUSTOM_REAL) deltatf,deltatover2f,deltatsquareover2f
  logical :: ANY_ANISOTROPY
  integer, dimension(:,:), allocatable  :: gather_ispec_selected_rec
  real(kind=CUSTOM_REAL) b_deltatf,b_deltatover2f,b_deltatsquareover2f
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: kappastore,mustore, rhostore
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rho_vp,rho_vs
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: abs_boundary_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: abs_boundary_jacobian1Dw
  integer, dimension(:,:,:), allocatable :: abs_boundary_ij
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable ::  source_time_function_loc
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: free_surface_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: free_surface_jacobian2Dw
  integer, dimension(:,:,:), allocatable :: free_surface_ij
  integer, dimension(:), allocatable :: free_ac_ispec
  integer :: num_free_surface_faces
  integer, dimension(:), allocatable :: cote_abs
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
            c11store,c12store,c13store,c15store,c23store,c25store,c33store,c35store,c55store
  integer :: num_colors_outer_acoustic,num_colors_inner_acoustic
  integer, dimension(:), allocatable :: num_elem_colors_acoustic
  integer :: num_colors_outer_elastic,num_colors_inner_elastic
  integer, dimension(:), allocatable :: num_elem_colors_elastic
  integer :: nsources_local
  integer :: num_phase_ispec_elastic,nspec_inner_elastic,nspec_outer_elastic
  integer, dimension(:,:), allocatable :: phase_ispec_inner_elastic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: cosrot_irecf, sinrot_irecf
  integer, dimension(:), allocatable :: coupling_ac_el_ispec
  integer, dimension(:,:,:), allocatable :: coupling_ac_el_ij
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: coupling_ac_el_normal
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: coupling_ac_el_jacobian1Dw
  integer, dimension(:), allocatable :: ispec_selected_source_loc
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sourcearray_loc
  integer, dimension(:,:), allocatable :: phase_ispec_inner_acoustic
  integer, dimension(:), allocatable  :: tab_requests_send_recv_scalar
  integer, dimension(:), allocatable  :: b_tab_requests_send_recv_scalar
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_scalar_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_buffer_send_scalar_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_recv_scalar_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_buffer_recv_scalar_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_recv_vector_ext_mesh
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_recv_vector_ext_mesh
  integer, dimension(:), allocatable  :: tab_requests_send_recv_vector
  integer, dimension(:), allocatable  :: b_tab_requests_send_recv_vector
  integer :: num_phase_ispec_acoustic,nspec_inner_acoustic,nspec_outer_acoustic


end module specfem_par
