
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
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

!====================================================================================
!
!   An explicit 2D parallel MPI spectral element solver
!   for the anelastic anisotropic wave equation
!
!====================================================================================

! If you use this code for your own research, please cite:
!
! @ARTICLE{KoTr99,
! author={D. Komatitsch and J. Tromp},
! title={Introduction to the spectral-element method for 3-{D} seismic wave propagation},
! journal={Geophys. J. Int.},
! year=1999,
! volume=139,
! number=3,
! pages={806-822},
! doi={10.1046/j.1365-246x.1999.00967.x}}
!
! @ARTICLE{KoVi98,
! author={D. Komatitsch and J. P. Vilotte},
! title={The spectral-element method: an efficient tool to simulate the seismic response of 2{D} and 3{D} geological structures},
! journal={Bull. Seismol. Soc. Am.},
! year=1998,
! volume=88,
! number=2,
! pages={368-392}}
!
! If you use the METIS / SCOTCH / CUBIT non-structured version, please also cite:
!
! @INPROCEEDINGS{MaKoBlLe08,
! author = {R. Martin and D. Komatitsch and C. Blitz and N. {Le Goff}},
! title = {Simulation of seismic wave propagation in an asteroid based upon
! an unstructured {MPI} spectral-element method: blocking and non-blocking communication strategies}
! booktitle = {Proceedings of the VECPAR'2008 8th International Meeting
! on High Performance Computing for Computational Science},
! year = {2008},
! pages = {999998-999999},
! address = {Toulouse, France},
! note = {24-27 June 2008},
! url = {http://vecpar.fe.up.pt/2008}}

!
! version 5.2, Dimitri Komatitsch, Nicolas Le Goff and Roland Martin, February 2008:
!               - MPI implementation of the code based on domain decomposition
!                 with METIS or SCOTCH
!               - general fluid/solid implementation with any number, shape and orientation of
!                 matching edges
!               - fluid potential of density * displacement instead of displacement
!               - absorbing edges with any normal vector
!               - general numbering of absorbing and acoustic free surface edges
!               - cleaned implementation of attenuation as in Carcione (1993)
!               - merged loops in the solver for efficiency
!               - simplified input of external model
!               - added CPU time information
!               - translated many comments from French to English
!
! version 5.1, Dimitri Komatitsch, January 2005:
!               - more general mesher with any number of curved layers
!               - Dirac and Gaussian time sources and corresponding convolution routine
!               - option for acoustic medium instead of elastic
!               - receivers at any location, not only grid points
!               - moment-tensor source at any location, not only a grid point
!               - color snapshots
!               - more flexible DATA/Par_file with any number of comment lines
!               - Xsu scripts for seismograms
!               - subtract t0 from seismograms
!               - seismograms and snapshots in pressure in addition to vector field
!
! version 5.0, Dimitri Komatitsch, May 2004:
!               - got rid of useless routines, suppressed commons etc.
!               - weak formulation based explicitly on stress tensor
!               - implementation of full anisotropy
!               - implementation of attenuation based on memory variables
!
! based on SPECFEM2D version 4.2, June 1998
! (c) by Dimitri Komatitsch, Harvard University, USA
! and Jean-Pierre Vilotte, Institut de Physique du Globe de Paris, France
!
! itself based on SPECFEM2D version 1.0, 1995
! (c) by Dimitri Komatitsch and Jean-Pierre Vilotte,
! Institut de Physique du Globe de Paris, France
!

! in case of an acoustic medium, a potential Chi of (density * displacement) is used as in Chaljub and Valette,
! Geophysical Journal International, vol. 158, p. 131-141 (2004) and *NOT* a velocity potential
! as in Komatitsch and Tromp, Geophysical Journal International, vol. 150, p. 303-318 (2002).
! This permits acoustic-elastic coupling based on a non-iterative time scheme.
! Displacement is then: u = grad(Chi) / rho
! Velocity is then: v = grad(Chi_dot) / rho (Chi_dot being the time derivative of Chi)
! and pressure is: p = - Chi_dot_dot  (Chi_dot_dot being the time second derivative of Chi).
! The source in an acoustic element is a pressure source.
! First-order acoustic-acoustic discontinuities are also handled automatically
! because pressure is continuous at such an interface, therefore Chi_dot_dot
! is continuous, therefore Chi is also continuous, which is consistent with
! the spectral-element basis functions and with the assembling process.
! This is the reason why a simple displacement potential u = grad(Chi) would
! not work because it would be discontinuous at such an interface and would
! therefore not be consistent with the basis functions.

  program specfem2D

  implicit none

  include "constants.h"
#ifdef USE_MPI
  include "mpif.h"
#endif

  character(len=80) datlin

  integer NSOURCE,i_source
  integer, dimension(:), allocatable :: source_type,time_function_type
  double precision, dimension(:), allocatable :: x_source,z_source,xi_source,gamma_source,&
                  Mxx,Mzz,Mxz,f0,t0,factor,angleforce,hdur,hdur_gauss
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: sourcearray

  double precision, dimension(:,:), allocatable :: coorg
  double precision, dimension(:), allocatable :: coorgread

! receiver information
  integer :: nrec,ios
  integer, dimension(:), allocatable :: ispec_selected_rec
  double precision, dimension(:), allocatable :: xi_receiver,gamma_receiver,st_xval,st_zval
  character(len=150) dummystring

! for seismograms
  double precision, dimension(:,:), allocatable :: sisux,sisuz,siscurl
  integer :: seismo_offset, seismo_current

! vector field in an element
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLX) :: vector_field_element

! pressure in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: pressure_element

! curl in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: curl_element

  integer :: i,j,k,l,it,irec,ipoin,ip,id,n,ispec,npoin,npgeo,iglob
  logical :: anyabs
  double precision :: dxd,dzd,dcurld,valux,valuz,valcurl,hlagrange,rhol,cosrot,sinrot,xi,gamma,x,z

! coefficients of the explicit Newmark time scheme
  integer NSTEP
  double precision deltatover2,deltatsquareover2,time,deltat

! Gauss-Lobatto-Legendre points and weights
  double precision, dimension(NGLLX) :: xigll
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLZ) :: zigll
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wzgll

! derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl

! material properties of the elastic medium
  double precision :: mul_relaxed,lambdal_relaxed,lambdalplus2mul_relaxed,kappal

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_elastic,veloc_elastic,displ_elastic
  double precision, dimension(:,:), allocatable :: coord, flagrange,xinterp,zinterp,Uxinterp,Uzinterp,vector_field_display

! material properties of the poroelastic medium (solid phase:s and fluid phase [defined as w=phi(u_f-u_s)]: w)
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accels_poroelastic,velocs_poroelastic,displs_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accelw_poroelastic,velocw_poroelastic,displw_poroelastic
  double precision, dimension(:), allocatable :: porosity,tortuosity
  double precision, dimension(:,:), allocatable :: density,permeability

! poroelastic and elastic coefficients 
  double precision, dimension(:,:,:), allocatable :: poroelastcoef

! for acoustic medium
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic

! inverse mass matrices 
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_elastic,rmass_inverse_acoustic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic

! to evaluate cpI, cpII, and cs, and rI (poroelastic medium)
  real(kind=CUSTOM_REAL) :: rhol_s,rhol_f,rhol_bar,phil,tortl
  real(kind=CUSTOM_REAL) :: mul_s,kappal_s
  real(kind=CUSTOM_REAL) :: kappal_f
!  double precision :: etal_f
  real(kind=CUSTOM_REAL) :: mul_fr,kappal_fr
!  double precision :: permlxx,permlxz,permlzz
  real(kind=CUSTOM_REAL) :: afactor,bfactor,cfactor,D_biot,H_biot,C_biot,M_biot,B_biot,cpIsquare,cpIIsquare,cssquare
  real(kind=CUSTOM_REAL) :: gamma1,gamma2,gamma3,gamma4,ratio,dd1

  double precision, dimension(:), allocatable :: vp_display

  double precision, dimension(:,:,:), allocatable :: vpext,vsext,rhoext
  double precision :: previous_vsext,rho_at_source_location

  double precision, dimension(:,:,:), allocatable :: shape2D,shape2D_display
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: xix,xiz,gammax,gammaz,jacobian

  double precision, dimension(:,:,:,:), allocatable :: dershape2D,dershape2D_display

  integer, dimension(:,:,:), allocatable :: ibool,ibool_outer,ibool_inner
  integer, dimension(:,:), allocatable  :: knods
  integer, dimension(:), allocatable :: kmato,numabs, &
     ibegin_bottom,iend_bottom,ibegin_top,iend_top,jbegin_left,jend_left,jbegin_right,jend_right

  integer, dimension(:), allocatable :: ispec_selected_source,iglob_source,ix_source,iz_source,&
                                        is_proc_source,nb_proc_source
  double precision displnorm_all,displnorm_all_glob
  double precision, dimension(:), allocatable :: aval
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: source_time_function
  double precision, external :: netlib_specfun_erf

  double precision :: vpImin,vpImax,vpIImin,vpIImax

  integer :: colors,numbers,subsamp,imagetype,NTSTEP_BETWEEN_OUTPUT_INFO,NTSTEP_BETWEEN_OUTPUT_SEISMO,seismotype
  integer :: numat,ngnod,nspec,pointsdisp,nelemabs,nelem_acoustic_surface,ispecabs,UPPER_LIMIT_DISPLAY

  logical interpol,meshvect,modelvect,boundvect,assign_external_model,initialfield, &
    outputgrid,gnuplot,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON,output_postscript_snapshot,output_color_image, &
    plot_lowerleft_corner_only,add_Bielak_conditions

  double precision :: cutsnaps,sizemax_arrows,anglerec,xirec,gammarec

! for absorbing and acoustic free surface conditions
  integer :: ispec_acoustic_surface,inum,numabsread
  logical :: codeabsread(4)
  real(kind=CUSTOM_REAL) :: nx,nz,weight,xxi,zgamma

  logical, dimension(:,:), allocatable  :: codeabs

! for attenuation
  integer  :: N_SLS
  double precision, dimension(:), allocatable  :: Qp_attenuation
  double precision, dimension(:), allocatable  :: Qs_attenuation
  double precision  :: f0_attenuation
  integer nspec_allocate
  double precision :: deltatsquare,deltatcube,deltatfourth,twelvedeltat,fourdeltatsquare

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: e1,e11,e13
  double precision, dimension(:,:,:,:), allocatable :: inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2
  double precision, dimension(:), allocatable :: inv_tau_sigma_nu1_sent,phi_nu1_sent,inv_tau_sigma_nu2_sent,phi_nu2_sent
  double precision, dimension(:,:,:) , allocatable :: Mu_nu1,Mu_nu2
  double precision :: Mu_nu1_sent,Mu_nu2_sent

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n,dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1

! for viscous attenuation
  double precision, dimension(:,:,:), allocatable :: &
    rx_viscous,rz_viscous,viscox,viscoz
  double precision :: theta_e,theta_s,Q0,freq0
  double precision :: alphaval,betaval,gammaval,thetainv
  logical :: TURN_VISCATTENUATION_ON
  double precision, dimension(NGLLX,NGLLZ) :: viscox_loc,viscoz_loc
  double precision :: Sn,Snp1,etal_f
  double precision, dimension(3):: bl_relaxed
  double precision :: permlxx,permlxz,permlzz,invpermlxx,invpermlxz,invpermlzz,detk
! adjoint
  double precision, dimension(:), allocatable :: b_viscodampx,b_viscodampz
  integer reclen,reclen1,reclen2

! for fluid/solid coupling and edge detection
  logical, dimension(:), allocatable :: elastic
  integer, dimension(NEDGES) :: i_begin,j_begin,i_end,j_end
  integer, dimension(NGLLX,NEDGES) :: ivalue,jvalue,ivalue_inverse,jvalue_inverse
  integer, dimension(:), allocatable :: fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
                                        fluid_solid_elastic_ispec,fluid_solid_elastic_iedge
  integer :: fluid_solid_acoustic_ispec_read, fluid_solid_elastic_ispec_read
  integer :: num_fluid_solid_edges,ispec_acoustic,ispec_elastic, &
             iedge_acoustic,iedge_elastic,ipoin1D,iglob2
  logical :: any_acoustic,any_acoustic_glob,any_elastic,any_elastic_glob,coupled_acoustic_elastic
  real(kind=CUSTOM_REAL) :: displ_x,displ_z,displ_n,displw_x,displw_z,zxi,xgamma,jacobian1D,pressure
  real(kind=CUSTOM_REAL) :: b_displ_x,b_displ_z,b_displw_x,b_displw_z,b_pressure

! for fluid/porous medium coupling and edge detection
  logical, dimension(:), allocatable :: poroelastic
  logical :: any_poroelastic,any_poroelastic_glob
  integer, dimension(:), allocatable :: fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge, &
                                        fluid_poro_poroelastic_ispec,fluid_poro_poroelastic_iedge
  integer :: fluid_poro_acoustic_ispec_read, fluid_poro_poroelastic_ispec_read
  integer :: num_fluid_poro_edges,num_fluid_poro_edges_alloc,iedge_poroelastic
  logical :: coupled_acoustic_poroelastic
  double precision :: mul_G,lambdal_G,lambdalplus2mul_G
  double precision :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  double precision :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  double precision :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  double precision :: dwx_dxl,dwz_dxl,dwx_dzl,dwz_dzl
  double precision :: b_dux_dxi,b_dux_dgamma,b_duz_dxi,b_duz_dgamma
  double precision :: b_dwx_dxi,b_dwx_dgamma,b_dwz_dxi,b_dwz_dgamma
  double precision :: b_dux_dxl,b_duz_dxl,b_dux_dzl,b_duz_dzl
  double precision :: b_dwx_dxl,b_dwz_dxl,b_dwx_dzl,b_dwz_dzl

! for solid/porous medium coupling and edge detection
  integer, dimension(:), allocatable :: solid_poro_elastic_ispec,solid_poro_elastic_iedge, &
                                        solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge
  integer :: solid_poro_elastic_ispec_read, solid_poro_poroelastic_ispec_read
  integer :: num_solid_poro_edges,num_solid_poro_edges_alloc,ispec_poroelastic,ii2,jj2
  logical :: coupled_elastic_poroelastic
  double precision, dimension(:,:), allocatable :: displ,veloc
  double precision :: sigma_xx,sigma_xz,sigma_zz
  double precision :: b_sigma_xx,b_sigma_xz,b_sigma_zz
  integer, dimension(:), allocatable :: ibegin_bottom_poro,iend_bottom_poro,ibegin_top_poro,&
            iend_top_poro,jbegin_left_poro,jend_left_poro,jbegin_right_poro,jend_right_poro

! for adjoint method
  logical :: save_forward ! whether or not the last frame is saved to reconstruct the forward field
  integer :: isolver      ! 1 = forward wavefield, 2 = backward and adjoint wavefields and kernels
  double precision :: b_deltatover2,b_deltatsquareover2,b_deltat ! coefficients of the explicit Newmark time scheme
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_accels_poroelastic,b_velocs_poroelastic,b_displs_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_accelw_poroelastic,b_velocw_poroelastic,b_displw_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_accel_elastic,b_veloc_elastic,b_displ_elastic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rho_kl, mu_kl, kappa_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhol_global, mul_global, kappal_global
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: mu_k, kappa_k,rho_k
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhop_kl, beta_kl, alpha_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rho_ac_kl, kappa_ac_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhol_ac_global, kappal_ac_global
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: kappa_ac_k,rho_ac_k
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhop_ac_kl, alpha_ac_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhot_kl, rhof_kl, sm_kl, eta_kl, mufr_kl, B_kl, &
    C_kl, M_kl, rhob_kl, rhofb_kl, phi_kl, Bb_kl, Cb_kl, Mb_kl, mufrb_kl, &
    rhobb_kl, rhofbb_kl, phib_kl, cpI_kl, cpII_kl, cs_kl, ratio_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhot_k, rhof_k, sm_k, eta_k, mufr_k, B_k, &
    C_k, M_k
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: phil_global,etal_f_global,rhol_s_global,rhol_f_global,rhol_bar_global, &
    tortl_global,mulfr_global
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: permlxx_global,permlxz_global,permlzz_global
  character(len=150) :: adj_source_file,filename,filename2,filename3
  integer :: irec_local,nadj_rec_local
  double precision :: xx,zz,rholb
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: adj_sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: adj_sourcearrays
  double precision :: rhopmin,rhopmax,alphamin,alphamax,betamin,betamax
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_absorb_elastic_left,b_absorb_poro_s_left,b_absorb_poro_w_left
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_absorb_elastic_right,b_absorb_poro_s_right,b_absorb_poro_w_right
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_absorb_elastic_bottom,b_absorb_poro_s_bottom,b_absorb_poro_w_bottom
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_absorb_elastic_top,b_absorb_poro_s_top,b_absorb_poro_w_top
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable ::  b_absorb_acoustic_left,b_absorb_acoustic_right,&
                      b_absorb_acoustic_bottom, b_absorb_acoustic_top
  integer :: nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax
  integer, dimension(:), allocatable :: ib_xmin,ib_xmax,ib_zmin,ib_zmax

! for color images
  integer :: NX_IMAGE_color,NZ_IMAGE_color
  integer  :: npgeo_glob
  double precision :: xmin_color_image,xmax_color_image, &
    zmin_color_image,zmax_color_image,size_pixel_horizontal,size_pixel_vertical
  integer, dimension(:,:), allocatable :: iglob_image_color,copy_iglob_image_color
  double precision, dimension(:,:), allocatable :: image_color_data
  double precision, dimension(:,:), allocatable :: image_color_vp_display

  double precision  :: xmin_color_image_loc, xmax_color_image_loc, zmin_color_image_loc, &
       zmax_color_image_loc
  integer  :: min_i, min_j, max_i, max_j
  integer  :: nb_pixel_loc
  integer, dimension(:), allocatable  :: nb_pixel_per_proc
  double precision  :: i_coord, j_coord
  double precision, dimension(2,4)  :: elmnt_coords
  integer, dimension(:), allocatable  :: num_pixel_loc
  integer, dimension(:,:), allocatable  :: num_pixel_recv
  double precision, dimension(:), allocatable  :: data_pixel_recv
  double precision, dimension(:), allocatable  :: data_pixel_send
  logical  :: pixel_is_in
  double precision  :: dist_pixel, dist_min_pixel
#ifdef USE_MPI
  integer, dimension(MPI_STATUS_SIZE)  :: request_mpi_status
#endif

! timing information for the stations
  character(len=MAX_LENGTH_STATION_NAME), allocatable, dimension(:) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), allocatable, dimension(:) :: network_name

! title of the plot
  character(len=60) simulation_title

! Lagrange interpolators at receivers
  double precision, dimension(:), allocatable :: hxir,hgammar,hpxir,hpgammar
  double precision, dimension(:,:), allocatable :: hxir_store,hgammar_store

! for Lagrange interpolants
  double precision, external :: hgll

! timer to count elapsed time
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  integer ihours,iminutes,iseconds,int_tCPU
  double precision :: time_start,time_end,tCPU

! for MPI and partitioning
  integer  :: ier
  integer  :: nproc
  integer  :: myrank
  integer  :: iproc
  character(len=256)  :: prname

  integer  :: ninterface
  integer  :: max_interface_size
  integer, dimension(:), allocatable  :: my_neighbours
  integer, dimension(:), allocatable  :: my_nelmnts_neighbours
  integer, dimension(:,:,:), allocatable  :: my_interfaces
  integer, dimension(:,:), allocatable  :: ibool_interfaces_acoustic,ibool_interfaces_elastic,ibool_interfaces_poroelastic
  integer, dimension(:), allocatable  :: nibool_interfaces_acoustic,nibool_interfaces_elastic,nibool_interfaces_poroelastic

  integer  :: ninterface_acoustic, ninterface_elastic,ninterface_poroelastic
  integer, dimension(:), allocatable  :: inum_interfaces_acoustic, inum_interfaces_elastic, inum_interfaces_poroelastic

#ifdef USE_MPI
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_send_faces_vector_ac
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_recv_faces_vector_ac
  integer, dimension(:), allocatable  :: tab_requests_send_recv_acoustic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_send_faces_vector_el
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_recv_faces_vector_el
  integer, dimension(:), allocatable  :: tab_requests_send_recv_elastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_send_faces_vector_pos,buffer_send_faces_vector_pow
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_recv_faces_vector_pos,buffer_recv_faces_vector_pow
  integer, dimension(:), allocatable  :: tab_requests_send_recv_poroelastic
  integer  :: max_ibool_interfaces_size_ac, max_ibool_interfaces_size_el, max_ibool_interfaces_size_po
#endif

! for overlapping MPI communications with computation
  integer  :: nspec_outer, nspec_inner, num_ispec_outer, num_ispec_inner
  integer, dimension(:), allocatable  :: ispec_outer_to_glob, ispec_inner_to_glob
  logical, dimension(:), allocatable  :: mask_ispec_inner_outer

  integer, dimension(:,:), allocatable  :: acoustic_surface
  integer :: acoustic_edges_read
  integer, dimension(:,:), allocatable  :: acoustic_edges

  integer  :: ixmin, ixmax, izmin, izmax

  integer  :: ie, num_interface

  integer  :: nrecloc, irecloc
  integer, dimension(:), allocatable :: recloc, which_proc_receiver

! mask to sort ibool
  integer, dimension(:), allocatable :: mask_ibool
  integer, dimension(:,:,:), allocatable :: copy_ibool_ori
  integer :: inumber

! to compute analytical initial plane wave field
  double precision :: t
  double precision, external :: ricker_Bielak_displ,ricker_Bielak_veloc,ricker_Bielak_accel
  double precision :: angleforce_refl, c_inc, c_refl, cploc, csloc, denst, lambdaplus2mu, mu, p
  double precision, dimension(2) :: A_plane, B_plane, C_plane
  double precision :: PP, PS, SP, SS, z0_source, x0_source, xmax, xmin, zmax, zmin, time_offset
#ifdef USE_MPI
  double precision :: xmax_glob, xmin_glob, zmax_glob, zmin_glob
#endif

! beyond critical angle
  integer , dimension(:), allocatable :: left_bound,right_bound,bot_bound
  double precision , dimension(:,:), allocatable :: v0x_left,v0z_left,v0x_right,v0z_right,v0x_bot,v0z_bot
  double precision , dimension(:,:), allocatable :: t0x_left,t0z_left,t0x_right,t0z_right,t0x_bot,t0z_bot
  integer count_left,count_right,count_bot,ibegin,iend
  logical :: over_critical_angle

! further reduce cache misses inner/outer in two passes in the case of an MPI simulation
  integer :: ipass,ispec_inner,ispec_outer,NUMBER_OF_PASSES,kmato_read,my_interfaces_read
  integer :: npoin_outer,npoin_inner
  integer, dimension(:), allocatable :: knods_read
  integer, dimension(:), allocatable :: perm,antecedent_list,check_perm

! arrays for plotpost
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

! tangential detection
  double precision, dimension(:), allocatable :: anglerec_irec
  double precision, dimension(:), allocatable :: cosrot_irec, sinrot_irec
  double precision, dimension(:), allocatable :: x_final_receiver, z_final_receiver
  logical :: force_normal_to_surface,rec_normal_to_surface
  integer  :: nnodes_tangential_curve
  integer, dimension(:), allocatable :: source_courbe_eros
  double precision, dimension(:,:), allocatable  :: nodes_tangential_curve
  integer  :: n1_tangential_detection_curve
  integer, dimension(4)  :: n_tangential_detection_curve
  integer, dimension(:), allocatable  :: rec_tangential_detection_curve
  double precision :: distmin, dist_current, angleforce_recv
  double precision, dimension(:), allocatable :: dist_tangential_detection_curve
  double precision :: x_final_receiver_dummy, z_final_receiver_dummy

!***********************************************************************
!
!             i n i t i a l i z a t i o n    p h a s e
!
!***********************************************************************

#ifdef USE_MPI
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

! this is only used in the case of MPI because it distinguishes between inner and outer element
! in the MPI partitions, which is meaningless in the serial case
  if(FURTHER_REDUCE_CACHE_MISSES) then
    NUMBER_OF_PASSES = 2
  else
    NUMBER_OF_PASSES = 1
  endif

#else
  nproc = 1
  myrank = 0
  ier = 0
  ninterface_acoustic = 0
  ninterface_elastic = 0
  ninterface_poroelastic = 0
  iproc = 0
  ispec_inner = 0
  ispec_outer = 0

  if(PERFORM_CUTHILL_MCKEE) then
    NUMBER_OF_PASSES = 2
  else
    NUMBER_OF_PASSES = 1
  endif
#endif

! determine if we write to file instead of standard output
  if(IOUT /= ISTANDARD_OUTPUT) then
#ifdef USE_MPI
    write(prname,240) myrank
 240 format('simulation_results',i5.5,'.txt')
#else
    prname = 'simulation_results.txt'
#endif
    open(IOUT,file=prname,status='unknown',action='write')
  endif

! reduction of cache misses inner/outer in two passes
  do ipass = 1,NUMBER_OF_PASSES

  write(prname,230) myrank
 230 format('./OUTPUT_FILES/Database',i5.5)
  open(unit=IIN,file=prname,status='old',action='read')

!
!---  read job title and skip remaining titles of the input file
!
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a50)") simulation_title

!
!---- print the date, time and start-up banner
!
  if (myrank == 0 .and. ipass == 1) call datim(simulation_title)

  if (myrank == 0 .and. ipass == 1) then
    write(IOUT,*)
    write(IOUT,*)
    write(IOUT,*) '*********************'
    write(IOUT,*) '****             ****'
    write(IOUT,*) '****  SPECFEM2D  ****'
    write(IOUT,*) '****             ****'
    write(IOUT,*) '*********************'
  endif

!
!---- read parameters from input file
!

  read(IIN,"(a80)") datlin
  read(IIN,*) npgeo

  read(IIN,"(a80)") datlin
  read(IIN,*) gnuplot,interpol

  read(IIN,"(a80)") datlin
  read(IIN,*) NTSTEP_BETWEEN_OUTPUT_INFO

  read(IIN,"(a80)") datlin
  read(IIN,*) output_postscript_snapshot,output_color_image,colors,numbers

  read(IIN,"(a80)") datlin
  read(IIN,*) meshvect,modelvect,boundvect,cutsnaps,subsamp,sizemax_arrows
  cutsnaps = cutsnaps / 100.d0

  read(IIN,"(a80)") datlin
  read(IIN,*) anglerec

  read(IIN,"(a80)") datlin
  read(IIN,*) initialfield,add_Bielak_conditions
  if(add_Bielak_conditions .and. .not. initialfield) stop 'need to have an initial field to add Bielak plane wave conditions'

  read(IIN,"(a80)") datlin
  read(IIN,*) seismotype,imagetype,save_forward
  if(seismotype < 1 .or. seismotype > 6) call exit_MPI('Wrong type for seismogram output')
  if(imagetype < 1 .or. imagetype > 4) call exit_MPI('Wrong type for snapshots')

  if(save_forward .and. (seismotype /= 1 .and. seismotype /= 6)) then
  print*, '***** WARNING *****'
  print*, 'seismotype =',seismotype
  print*, 'Save forward wavefield => seismogram must be in displacement for (poro)elastic or potential for acoustic'
  print*, 'Seismotype must be changed to 1 (elastic/poroelastic adjoint source) or 6 (acoustic adjoint source)'
  stop
  endif

  read(IIN,"(a80)") datlin
  read(IIN,*) assign_external_model,outputgrid,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON

  read(IIN,"(a80)") datlin
  read(IIN,*) TURN_VISCATTENUATION_ON,Q0,freq0

!---- check parameters read
  if (myrank == 0 .and. ipass == 1) then
    write(IOUT,200) npgeo,NDIM
    write(IOUT,600) NTSTEP_BETWEEN_OUTPUT_INFO,colors,numbers
    write(IOUT,700) seismotype,anglerec
    write(IOUT,750) initialfield,add_Bielak_conditions,assign_external_model,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON,outputgrid
    write(IOUT,800) imagetype,100.d0*cutsnaps,subsamp
  endif

!---- read time step
  read(IIN,"(a80)") datlin
  read(IIN,*) NSTEP,deltat,isolver
  if (myrank == 0 .and. ipass == 1) write(IOUT,703) NSTEP,deltat,NSTEP*deltat

  if(isolver == 1 .and. save_forward .and. (TURN_ANISOTROPY_ON .or. TURN_ATTENUATION_ON .or. TURN_VISCATTENUATION_ON)) then
  print*, '*************** WARNING ***************'
  print*, 'Anisotropy & Attenuation & Viscous damping are not presently implemented for adjoint calculations'
  stop 
  endif

  NTSTEP_BETWEEN_OUTPUT_SEISMO = min(NSTEP,NTSTEP_BETWEEN_OUTPUT_INFO)

!
!----  read source information
!
  read(IIN,"(a80)") datlin
  read(IIN,*) NSOURCE
  if(ipass == 1) then
  allocate( source_type(NSOURCE) )
  allocate( time_function_type(NSOURCE) )
  allocate( x_source(NSOURCE) )
  allocate( z_source(NSOURCE) )
  allocate( f0(NSOURCE) )
  allocate( t0(NSOURCE) )
  allocate( factor(NSOURCE) )
  allocate( angleforce(NSOURCE) )
  allocate( hdur(NSOURCE) )
  allocate( hdur_gauss(NSOURCE) )
  allocate( Mxx(NSOURCE) )
  allocate( Mxz(NSOURCE) )
  allocate( Mzz(NSOURCE) )
  allocate( aval(NSOURCE) )
  allocate( ispec_selected_source(NSOURCE) )
  allocate( iglob_source(NSOURCE) )
  allocate( source_courbe_eros(NSOURCE) )
  allocate( ix_source(NSOURCE) )
  allocate( iz_source(NSOURCE) )
  allocate( xi_source(NSOURCE) )
  allocate( gamma_source(NSOURCE) )
  allocate( is_proc_source(NSOURCE) )
  allocate( nb_proc_source(NSOURCE) )
  allocate( sourcearray(NSOURCE,NDIM,NGLLX,NGLLZ) )
  endif

  do i_source=1,NSOURCE
     read(IIN,"(a80)") datlin
     read(IIN,*) source_type(i_source),time_function_type(i_source),x_source(i_source),z_source(i_source), &
                 f0(i_source),t0(i_source), &
                 factor(i_source),angleforce(i_source),Mxx(i_source),Mzz(i_source),Mxz(i_source)
  enddo

!
!----  read attenuation information
!
  read(IIN,"(a80)") datlin
  read(IIN,*) N_SLS, f0_attenuation

!
!-----  check the input
!
 do i_source=1,NSOURCE

 if(.not. initialfield) then
   if (source_type(i_source) == 1) then
     if ( myrank == 0 ) then
     write(IOUT,212) x_source(i_source),z_source(i_source),f0(i_source),t0(i_source), &
                     factor(i_source),angleforce(i_source)
     endif
   else if(source_type(i_source) == 2) then
     if ( myrank == 0 ) then
     write(IOUT,222) x_source(i_source),z_source(i_source),f0(i_source),t0(i_source), &
                     factor(i_source),Mxx(i_source),Mzz(i_source),Mxz(i_source)
     endif
   else
     call exit_MPI('Unknown source type number !')
   endif
 endif
! if Dirac source time function, use a very thin Gaussian instead
! if Heaviside source time function, use a very thin error function instead
! time delay of the source in seconds, use a 20 % security margin (use 2 / f0 if error function)
  if(time_function_type(i_source) == 4 .or. time_function_type(i_source) == 5) then
    f0(i_source) = 1.d0 / (10.d0 * deltat)
    if(time_function_type(i_source) == 5) then
      t0(i_source) = 2.0d0 / f0(i_source)
    else
      t0(i_source) = 1.20d0 / f0(i_source)
    endif
  endif

! for the source time function
  aval(i_source) = pi*pi*f0(i_source)*f0(i_source)

!-----  convert angle from degrees to radians
  angleforce(i_source) = angleforce(i_source) * pi / 180.d0

 enddo ! do i_source=1,NSOURCE

!
!---- read the spectral macrobloc nodal coordinates
!
  if(ipass == 1) allocate(coorg(NDIM,npgeo))

  ipoin = 0
  read(IIN,"(a80)") datlin
  allocate(coorgread(NDIM))
  do ip = 1,npgeo
   read(IIN,*) ipoin,(coorgread(id),id =1,NDIM)
   if(ipoin<1 .or. ipoin>npgeo) call exit_MPI('Wrong control point number')
   coorg(:,ipoin) = coorgread
  enddo
  deallocate(coorgread)

!
!---- read the basic properties of the spectral elements
!
  read(IIN,"(a80)") datlin
  read(IIN,*) numat,ngnod,nspec,pointsdisp,plot_lowerleft_corner_only
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,*) nelemabs,nelem_acoustic_surface,num_fluid_solid_edges,num_fluid_poro_edges,&
              num_solid_poro_edges,nnodes_tangential_curve

!
!---- allocate arrays
!
if(ipass == 1) then
  allocate(shape2D(ngnod,NGLLX,NGLLZ))
  allocate(dershape2D(NDIM,ngnod,NGLLX,NGLLZ))
  allocate(shape2D_display(ngnod,pointsdisp,pointsdisp))
  allocate(dershape2D_display(NDIM,ngnod,pointsdisp,pointsdisp))
  allocate(xix(NGLLX,NGLLZ,nspec))
  allocate(xiz(NGLLX,NGLLZ,nspec))
  allocate(gammax(NGLLX,NGLLZ,nspec))
  allocate(gammaz(NGLLX,NGLLZ,nspec))
  allocate(jacobian(NGLLX,NGLLZ,nspec))
  allocate(flagrange(NGLLX,pointsdisp))
  allocate(xinterp(pointsdisp,pointsdisp))
  allocate(zinterp(pointsdisp,pointsdisp))
  allocate(Uxinterp(pointsdisp,pointsdisp))
  allocate(Uzinterp(pointsdisp,pointsdisp))
  allocate(density(2,numat))
  allocate(porosity(numat))
  allocate(tortuosity(numat))
  allocate(permeability(3,numat))
  allocate(poroelastcoef(4,3,numat))
  allocate(Qp_attenuation(numat))
  allocate(Qs_attenuation(numat))
  allocate(kmato(nspec))
  allocate(knods(ngnod,nspec))
  allocate(ibool(NGLLX,NGLLZ,nspec))
  allocate(elastic(nspec))
  allocate(poroelastic(nspec))
  allocate(inv_tau_sigma_nu1(NGLLX,NGLLZ,nspec,N_SLS))
  allocate(inv_tau_sigma_nu2(NGLLX,NGLLZ,nspec,N_SLS))
  allocate(phi_nu1(NGLLX,NGLLZ,nspec,N_SLS))
  allocate(phi_nu2(NGLLX,NGLLZ,nspec,N_SLS))
  allocate(inv_tau_sigma_nu1_sent(N_SLS))
  allocate(inv_tau_sigma_nu2_sent(N_SLS))
  allocate(phi_nu1_sent(N_SLS))
  allocate(phi_nu2_sent(N_SLS))
endif

! --- allocate arrays for absorbing boundary conditions
  if(nelemabs <= 0) then
    nelemabs = 1
    anyabs = .false.
  else
    anyabs = .true.
  endif
if(ipass == 1) then
  allocate(numabs(nelemabs))
  allocate(codeabs(4,nelemabs))

  allocate(ibegin_bottom(nelemabs))
  allocate(iend_bottom(nelemabs))
  allocate(ibegin_top(nelemabs))
  allocate(iend_top(nelemabs))

  allocate(jbegin_left(nelemabs))
  allocate(jend_left(nelemabs))
  allocate(jbegin_right(nelemabs))
  allocate(jend_right(nelemabs))

  allocate(ibegin_bottom_poro(nelemabs))
  allocate(iend_bottom_poro(nelemabs))
  allocate(ibegin_top_poro(nelemabs))
  allocate(iend_top_poro(nelemabs))

  allocate(jbegin_left_poro(nelemabs))
  allocate(jend_left_poro(nelemabs))
  allocate(jbegin_right_poro(nelemabs))
  allocate(jend_right_poro(nelemabs))
endif

!
!---- print element group main parameters
!
  if (myrank == 0 .and. ipass == 1) then
    write(IOUT,107)
    write(IOUT,207) nspec,ngnod,NGLLX,NGLLZ,NGLLX*NGLLZ,pointsdisp,numat,nelemabs
  endif

! set up Gauss-Lobatto-Legendre derivation matrices
  call define_derivation_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,hprimewgll_xx,hprimewgll_zz)

!
!---- read the material properties
!
  call gmat01(density,porosity,tortuosity,permeability,poroelastcoef,numat,&
              myrank,ipass,Qp_attenuation,Qs_attenuation)

!
!----  read spectral macrobloc data
!
  n = 0
  read(IIN,"(a80)") datlin
  allocate(knods_read(ngnod))
  do ispec = 1,nspec
    read(IIN,*) n,kmato_read,(knods_read(k), k=1,ngnod)
if(ipass == 1) then
  kmato(n) = kmato_read
  knods(:,n)= knods_read(:)
else if(ipass == 2) then
  kmato(perm(antecedent_list(n))) = kmato_read
  knods(:,perm(antecedent_list(n)))= knods_read(:)
else
  stop 'error: maximum is 2 passes'
endif
  enddo
  deallocate(knods_read)

!-------------------------------------------------------------------------------
!----  determine if each spectral element is elastic, poroelastic, or acoustic
!-------------------------------------------------------------------------------
  any_acoustic = .false.
  any_elastic = .false.
  any_poroelastic = .false.
  do ispec = 1,nspec

    if(porosity(kmato(ispec)) == 1.d0) then  ! acoustic domain
      elastic(ispec) = .false.
      poroelastic(ispec) = .false.
      any_acoustic = .true.
    elseif(porosity(kmato(ispec)) < TINYVAL) then  ! elastic domain
      elastic(ispec) = .true.
      poroelastic(ispec) = .false.
      any_elastic = .true.
    else                                       ! poroelastic domain
      elastic(ispec) = .false.
      poroelastic(ispec) = .true.
      any_poroelastic = .true.
    endif

  enddo !do ispec = 1,nspec

  if(TURN_ATTENUATION_ON) then
    nspec_allocate = nspec
  else
    nspec_allocate = 1
  endif

! allocate memory variables for attenuation
if(ipass == 1) then
  allocate(e1(NGLLX,NGLLZ,nspec_allocate,N_SLS))
  allocate(e11(NGLLX,NGLLZ,nspec_allocate,N_SLS))
  allocate(e13(NGLLX,NGLLZ,nspec_allocate,N_SLS))
  e1(:,:,:,:) = 0._CUSTOM_REAL
  e11(:,:,:,:) = 0._CUSTOM_REAL
  e13(:,:,:,:) = 0._CUSTOM_REAL

  allocate(dux_dxl_n(NGLLX,NGLLZ,nspec_allocate))
  allocate(duz_dzl_n(NGLLX,NGLLZ,nspec_allocate))
  allocate(duz_dxl_n(NGLLX,NGLLZ,nspec_allocate))
  allocate(dux_dzl_n(NGLLX,NGLLZ,nspec_allocate))
  allocate(dux_dxl_np1(NGLLX,NGLLZ,nspec_allocate))
  allocate(duz_dzl_np1(NGLLX,NGLLZ,nspec_allocate))
  allocate(duz_dxl_np1(NGLLX,NGLLZ,nspec_allocate))
  allocate(dux_dzl_np1(NGLLX,NGLLZ,nspec_allocate))
  allocate(Mu_nu1(NGLLX,NGLLZ,nspec))
  allocate(Mu_nu2(NGLLX,NGLLZ,nspec))
endif

! define the attenuation quality factors.
! they can be different for each element.
!! DK DK if needed in the future, here the quality factor could be different for each point
  do ispec = 1,nspec
    call attenuation_model(N_SLS,Qp_attenuation(kmato(ispec)),Qs_attenuation(kmato(ispec)), &
            f0_attenuation,inv_tau_sigma_nu1_sent,phi_nu1_sent,inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu1_sent,Mu_nu2_sent)
    do j = 1,NGLLZ
      do i = 1,NGLLX
        inv_tau_sigma_nu1(i,j,ispec,:) = inv_tau_sigma_nu1_sent(:)
        phi_nu1(i,j,ispec,:) = phi_nu1_sent(:)
        inv_tau_sigma_nu2(i,j,ispec,:) = inv_tau_sigma_nu2_sent(:)
        phi_nu2(i,j,ispec,:) = phi_nu2_sent(:)
        Mu_nu1(i,j,ispec) = Mu_nu1_sent
        Mu_nu2(i,j,ispec) = Mu_nu2_sent
      enddo
    enddo
  enddo

! allocate memory variables for viscous attenuation (poroelastic media)
  if(ipass == 1) then
  if(TURN_VISCATTENUATION_ON) then
  allocate(rx_viscous(NGLLX,NGLLZ,nspec))
  allocate(rz_viscous(NGLLX,NGLLZ,nspec))
  allocate(viscox(NGLLX,NGLLZ,nspec))
  allocate(viscoz(NGLLX,NGLLZ,nspec))
  else
  allocate(rx_viscous(NGLLX,NGLLZ,1))
  allocate(rz_viscous(NGLLX,NGLLZ,1))
  allocate(viscox(NGLLX,NGLLZ,1))
  allocate(viscoz(NGLLX,NGLLZ,1))
  endif
  endif

!
!----  read interfaces data
!

  read(IIN,"(a80)") datlin
  read(IIN,*) ninterface, max_interface_size
  if ( ninterface > 0 ) then
if(ipass == 1) then
     allocate(my_neighbours(ninterface))
     allocate(my_nelmnts_neighbours(ninterface))
     allocate(my_interfaces(4,max_interface_size,ninterface))
     allocate(ibool_interfaces_acoustic(NGLLX*max_interface_size,ninterface))
     allocate(ibool_interfaces_elastic(NGLLX*max_interface_size,ninterface))
     allocate(ibool_interfaces_poroelastic(NGLLX*max_interface_size,ninterface))
     allocate(nibool_interfaces_acoustic(ninterface))
     allocate(nibool_interfaces_elastic(ninterface))
     allocate(nibool_interfaces_poroelastic(ninterface))
     allocate(inum_interfaces_acoustic(ninterface))
     allocate(inum_interfaces_elastic(ninterface))
     allocate(inum_interfaces_poroelastic(ninterface))
endif

     do num_interface = 1, ninterface
        read(IIN,*) my_neighbours(num_interface), my_nelmnts_neighbours(num_interface)
        do ie = 1, my_nelmnts_neighbours(num_interface)
           read(IIN,*) my_interfaces_read, my_interfaces(2,ie,num_interface), &
                my_interfaces(3,ie,num_interface), my_interfaces(4,ie,num_interface)

           if(ipass == 1) then
             my_interfaces(1,ie,num_interface) = my_interfaces_read
           else if(ipass == 2) then
             my_interfaces(1,ie,num_interface) = perm(antecedent_list(my_interfaces_read))
           else
             stop 'error: maximum number of passes is 2'
           endif

        enddo
     enddo
  endif

!
!----  read absorbing boundary data
!
  read(IIN,"(a80)") datlin
  if(anyabs) then
     do inum = 1,nelemabs
      read(IIN,*) numabsread,codeabsread(1),codeabsread(2),codeabsread(3),codeabsread(4), ibegin_bottom(inum), iend_bottom(inum), &
           jbegin_right(inum), jend_right(inum), ibegin_top(inum), iend_top(inum), jbegin_left(inum), jend_left(inum)
      if(numabsread < 1 .or. numabsread > nspec) call exit_MPI('Wrong absorbing element number')
      if(ipass == 1) then
        numabs(inum) = numabsread
      else if(ipass == 2) then
        numabs(inum) = perm(antecedent_list(numabsread))
      else
        call exit_MPI('error: maximum number of passes is 2')
      endif
      codeabs(IBOTTOM,inum) = codeabsread(1)
      codeabs(IRIGHT,inum) = codeabsread(2)
      codeabs(ITOP,inum) = codeabsread(3)
      codeabs(ILEFT,inum) = codeabsread(4)
    enddo
    if (myrank == 0 .and. ipass == 1) then
      write(IOUT,*)
      write(IOUT,*) 'Number of absorbing elements: ',nelemabs
    endif

    nspec_xmin = ZERO
    nspec_xmax = ZERO
    nspec_zmin = ZERO
    nspec_zmax = ZERO
    if(ipass == 1) then
    allocate(ib_xmin(nelemabs))
    allocate(ib_xmax(nelemabs))
    allocate(ib_zmin(nelemabs))
    allocate(ib_zmax(nelemabs))
    endif
    do inum = 1,nelemabs
       if (codeabs(IBOTTOM,inum)) then
         nspec_zmin = nspec_zmin + 1
         ib_zmin(inum) =  nspec_zmin
       endif
       if (codeabs(IRIGHT,inum)) then
         nspec_xmax = nspec_xmax + 1
         ib_xmax(inum) =  nspec_xmax
       endif
       if (codeabs(ITOP,inum)) then
         nspec_zmax = nspec_zmax + 1
         ib_zmax(inum) = nspec_zmax
       endif
       if (codeabs(ILEFT,inum)) then
         nspec_xmin = nspec_xmin + 1
         ib_xmin(inum) =  nspec_xmin
       endif
    enddo
! Files to save absorbed waves needed to reconstruct backward wavefield for adjoint method
   if(ipass == 1) then
     if(any_elastic .and. (save_forward .or. isolver == 2)) then
   allocate(b_absorb_elastic_left(NDIM,NGLLZ,nspec_xmin,NSTEP))
   allocate(b_absorb_elastic_right(NDIM,NGLLZ,nspec_xmax,NSTEP))
   allocate(b_absorb_elastic_bottom(NDIM,NGLLX,nspec_zmin,NSTEP))
   allocate(b_absorb_elastic_top(NDIM,NGLLX,nspec_zmax,NSTEP))
     endif
     if(any_poroelastic .and. (save_forward .or. isolver == 2)) then
   allocate(b_absorb_poro_s_left(NDIM,NGLLZ,nspec_xmin,NSTEP))
   allocate(b_absorb_poro_s_right(NDIM,NGLLZ,nspec_xmax,NSTEP))
   allocate(b_absorb_poro_s_bottom(NDIM,NGLLX,nspec_zmin,NSTEP))
   allocate(b_absorb_poro_s_top(NDIM,NGLLX,nspec_zmax,NSTEP))
   allocate(b_absorb_poro_w_left(NDIM,NGLLZ,nspec_xmin,NSTEP))
   allocate(b_absorb_poro_w_right(NDIM,NGLLZ,nspec_xmax,NSTEP))
   allocate(b_absorb_poro_w_bottom(NDIM,NGLLX,nspec_zmin,NSTEP))
   allocate(b_absorb_poro_w_top(NDIM,NGLLX,nspec_zmax,NSTEP))
     endif
     if(any_acoustic .and. (save_forward .or. isolver == 2)) then
   allocate(b_absorb_acoustic_left(NGLLZ,nspec_xmin,NSTEP))
   allocate(b_absorb_acoustic_right(NGLLZ,nspec_xmax,NSTEP))
   allocate(b_absorb_acoustic_bottom(NGLLX,nspec_zmin,NSTEP))
   allocate(b_absorb_acoustic_top(NGLLX,nspec_zmax,NSTEP))
     endif
   endif

    write(IOUT,*)
    write(IOUT,*) 'nspec_xmin = ',nspec_xmin
    write(IOUT,*) 'nspec_xmax = ',nspec_xmax
    write(IOUT,*) 'nspec_zmin = ',nspec_zmin
    write(IOUT,*) 'nspec_zmax = ',nspec_zmax

  endif
!
!----  read acoustic free surface data
!
  read(IIN,"(a80)") datlin
  if(nelem_acoustic_surface > 0) then
     if(ipass == 1) allocate(acoustic_edges(4,nelem_acoustic_surface))
      do inum = 1,nelem_acoustic_surface
        read(IIN,*) acoustic_edges_read, acoustic_edges(2,inum), acoustic_edges(3,inum), &
             acoustic_edges(4,inum)
        if(ipass == 1) then
          acoustic_edges(1,inum) = acoustic_edges_read
        else if(ipass == 2) then
          acoustic_edges(1,inum) = perm(antecedent_list(acoustic_edges_read))
        else
          call exit_MPI('error: maximum number of passes is 2')
        endif

     enddo
     if(ipass == 1) allocate(acoustic_surface(5,nelem_acoustic_surface))
     call construct_acoustic_surface ( nspec, ngnod, knods, nelem_acoustic_surface, &
          acoustic_edges, acoustic_surface)
    if (myrank == 0 .and. ipass == 1) then
      write(IOUT,*)
      write(IOUT,*) 'Number of free surface elements: ',nelem_acoustic_surface
    endif
  else
    if(ipass == 1) then
      allocate(acoustic_edges(4,1))
      allocate(acoustic_surface(5,1))
    endif
  endif

!
!---- read acoustic elastic coupled edges
!
  read(IIN,"(a80)") datlin
  if ( num_fluid_solid_edges > 0 ) then
if(ipass == 1) then
     allocate(fluid_solid_acoustic_ispec(num_fluid_solid_edges))
     allocate(fluid_solid_acoustic_iedge(num_fluid_solid_edges))
     allocate(fluid_solid_elastic_ispec(num_fluid_solid_edges))
     allocate(fluid_solid_elastic_iedge(num_fluid_solid_edges))
endif
     do inum = 1, num_fluid_solid_edges
        read(IIN,*) fluid_solid_acoustic_ispec_read,fluid_solid_elastic_ispec_read
        if(ipass == 1) then
          fluid_solid_acoustic_ispec(inum) = fluid_solid_acoustic_ispec_read
          fluid_solid_elastic_ispec(inum) = fluid_solid_elastic_ispec_read
        else if(ipass == 2) then
          fluid_solid_acoustic_ispec(inum) = perm(antecedent_list(fluid_solid_acoustic_ispec_read))
          fluid_solid_elastic_ispec(inum) = perm(antecedent_list(fluid_solid_elastic_ispec_read))
        else
          call exit_MPI('error: maximum number of passes is 2')
        endif
     enddo
  else
if(ipass == 1) then
     allocate(fluid_solid_acoustic_ispec(1))
     allocate(fluid_solid_acoustic_iedge(1))
     allocate(fluid_solid_elastic_ispec(1))
     allocate(fluid_solid_elastic_iedge(1))
endif
  endif

!
!---- read acoustic poroelastic coupled edges
!
  read(IIN,"(a80)") datlin
  if ( num_fluid_poro_edges > 0 ) then
if(ipass == 1) then
     allocate(fluid_poro_acoustic_ispec(num_fluid_poro_edges))
     allocate(fluid_poro_acoustic_iedge(num_fluid_poro_edges))
     allocate(fluid_poro_poroelastic_ispec(num_fluid_poro_edges))
     allocate(fluid_poro_poroelastic_iedge(num_fluid_poro_edges))
endif
     do inum = 1, num_fluid_poro_edges
        read(IIN,*) fluid_poro_acoustic_ispec_read,fluid_poro_poroelastic_ispec_read
        if(ipass == 1) then
          fluid_poro_acoustic_ispec(inum) = fluid_poro_acoustic_ispec_read
          fluid_poro_poroelastic_ispec(inum) = fluid_poro_poroelastic_ispec_read
        else if(ipass == 2) then
          fluid_poro_acoustic_ispec(inum) = perm(antecedent_list(fluid_poro_acoustic_ispec_read))
          fluid_poro_poroelastic_ispec(inum) = perm(antecedent_list(fluid_poro_poroelastic_ispec_read))
        else
          call exit_MPI('error: maximum number of passes is 2')
        endif
     enddo
  else
if(ipass == 1) then
     allocate(fluid_poro_acoustic_ispec(1))
     allocate(fluid_poro_acoustic_iedge(1))
     allocate(fluid_poro_poroelastic_ispec(1))
     allocate(fluid_poro_poroelastic_iedge(1))
endif
  endif

!
!---- read poroelastic elastic coupled edges
!
  read(IIN,"(a80)") datlin
  if ( num_solid_poro_edges > 0 ) then
if(ipass == 1) then
     allocate(solid_poro_elastic_ispec(num_solid_poro_edges))
     allocate(solid_poro_elastic_iedge(num_solid_poro_edges))
     allocate(solid_poro_poroelastic_ispec(num_solid_poro_edges))
     allocate(solid_poro_poroelastic_iedge(num_solid_poro_edges))
endif
     do inum = 1, num_solid_poro_edges
        read(IIN,*) solid_poro_poroelastic_ispec_read,solid_poro_elastic_ispec_read
        if(ipass == 1) then
          solid_poro_elastic_ispec(inum) = solid_poro_elastic_ispec_read
          solid_poro_poroelastic_ispec(inum) = solid_poro_poroelastic_ispec_read
        else if(ipass == 2) then
          solid_poro_elastic_ispec(inum) = perm(antecedent_list(solid_poro_elastic_ispec_read))
          solid_poro_poroelastic_ispec(inum) = perm(antecedent_list(solid_poro_poroelastic_ispec_read))
        else
          call exit_MPI('error: maximum number of passes is 2')
        endif
     enddo
  else
if(ipass == 1) then
     allocate(solid_poro_elastic_ispec(1))
     allocate(solid_poro_elastic_iedge(1))
     allocate(solid_poro_poroelastic_ispec(1))
     allocate(solid_poro_poroelastic_iedge(1))
endif
  endif

!
!---- read tangential detection curve
!
  read(IIN,"(a80)") datlin
  read(IIN,*) force_normal_to_surface,rec_normal_to_surface
  if (nnodes_tangential_curve > 0) then
if (ipass == 1) then
    allocate(nodes_tangential_curve(2,nnodes_tangential_curve))
    allocate(dist_tangential_detection_curve(nnodes_tangential_curve))
endif
    do i = 1, nnodes_tangential_curve
      read(IIN,*) nodes_tangential_curve(1,i),nodes_tangential_curve(2,i)
    enddo
  else
    force_normal_to_surface = .false.
    rec_normal_to_surface = .false.
    nnodes_tangential_curve = 0
if (ipass == 1) then
    allocate(nodes_tangential_curve(2,1))
    allocate(dist_tangential_detection_curve(1))
endif
  endif

!
!---- close input file
!
  close(IIN)

!
!---- compute shape functions and their derivatives for SEM grid
!
  do j = 1,NGLLZ
    do i = 1,NGLLX
      call define_shape_functions(shape2D(:,i,j),dershape2D(:,:,i,j),xigll(i),zigll(j),ngnod)
    enddo
  enddo

!
!---- generate the global numbering
!

! "slow and clean" or "quick and dirty" version
  if(FAST_NUMBERING) then
    call createnum_fast(knods,ibool,shape2D,coorg,npoin,npgeo,nspec,ngnod,myrank,ipass)
  else
    call createnum_slow(knods,ibool,npoin,nspec,ngnod,myrank,ipass)
  endif

! create a new indirect addressing array to reduce cache misses in memory access in the solver
  if(ipass == 2) then

  deallocate(perm)

  allocate(copy_ibool_ori(NGLLX,NGLLZ,nspec))
  allocate(mask_ibool(npoin))

  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:) = ibool(:,:,:)

  inumber = 0

  if(.not. ACTUALLY_IMPLEMENT_PERM_WHOLE) then

! first reduce cache misses in outer elements, since they are taken first
! loop over spectral elements
  do ispec = 1,nspec_outer
    do j=1,NGLLZ
      do i=1,NGLLX
        if(mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
! create a new point
          inumber = inumber + 1
          ibool(i,j,ispec) = inumber
          mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
        else
! use an existing point created previously
          ibool(i,j,ispec) = mask_ibool(copy_ibool_ori(i,j,ispec))
        endif
      enddo
    enddo
  enddo

! then reduce cache misses in inner elements, since they are taken second
! loop over spectral elements
  do ispec = nspec_outer+1,nspec
    do j=1,NGLLZ
      do i=1,NGLLX
        if(mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
! create a new point
          inumber = inumber + 1
          ibool(i,j,ispec) = inumber
          mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
        else
! use an existing point created previously
          ibool(i,j,ispec) = mask_ibool(copy_ibool_ori(i,j,ispec))
        endif
      enddo
    enddo
  enddo

  else ! if ACTUALLY_IMPLEMENT_PERM_WHOLE

! reduce cache misses in all the elements
! loop over spectral elements
  do ispec = 1,nspec
    do j=1,NGLLZ
      do i=1,NGLLX
        if(mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
! create a new point
          inumber = inumber + 1
          ibool(i,j,ispec) = inumber
          mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
        else
! use an existing point created previously
          ibool(i,j,ispec) = mask_ibool(copy_ibool_ori(i,j,ispec))
        endif
      enddo
    enddo
  enddo

  endif

  deallocate(copy_ibool_ori)
  deallocate(mask_ibool)

  else if(ipass /= 1) then
    stop 'incorrect pass number for reduction of cache misses'
  endif

!---- compute shape functions and their derivatives for regular interpolated display grid
  do j = 1,pointsdisp
    do i = 1,pointsdisp
      xirec  = 2.d0*dble(i-1)/dble(pointsdisp-1) - 1.d0
      gammarec  = 2.d0*dble(j-1)/dble(pointsdisp-1) - 1.d0
      call define_shape_functions(shape2D_display(:,i,j),dershape2D_display(:,:,i,j),xirec,gammarec,ngnod)
    enddo
  enddo

!---- compute Lagrange interpolants on a regular interpolated grid in (xi,gamma)
!---- for display (assumes NGLLX = NGLLZ)
  do j=1,NGLLX
    do i=1,pointsdisp
      xirec  = 2.d0*dble(i-1)/dble(pointsdisp-1) - 1.d0
      flagrange(j,i) = hgll(j-1,xirec,xigll,NGLLX)
    enddo
  enddo

! get number of stations from receiver file
  open(unit=IIN,file='DATA/STATIONS',iostat=ios,status='old',action='read')
  nrec = 0
  do while(ios == 0)
    read(IIN,"(a)",iostat=ios) dummystring
    if(ios == 0) nrec = nrec + 1
  enddo
  close(IIN)

  if (myrank == 0 .and. ipass == 1) then
    write(IOUT,*)
    write(IOUT,*) 'Total number of receivers = ',nrec
    write(IOUT,*)
  endif

  if(nrec < 1) call exit_MPI('need at least one receiver')

! receiver information
  if(ipass == 1) then

  allocate(ispec_selected_rec(nrec))
  allocate(st_xval(nrec))
  allocate(st_zval(nrec))
  allocate(xi_receiver(nrec))
  allocate(gamma_receiver(nrec))
  allocate(station_name(nrec))
  allocate(network_name(nrec))
  allocate(recloc(nrec))
  allocate(which_proc_receiver(nrec))
  allocate(x_final_receiver(nrec))
  allocate(z_final_receiver(nrec))

! allocate 1-D Lagrange interpolators and derivatives
  allocate(hxir(NGLLX))
  allocate(hpxir(NGLLX))
  allocate(hgammar(NGLLZ))
  allocate(hpgammar(NGLLZ))

! allocate Lagrange interpolators for receivers
  allocate(hxir_store(nrec,NGLLX))
  allocate(hgammar_store(nrec,NGLLZ))

! allocate other global arrays
  allocate(coord(NDIM,npoin))

! to display acoustic elements
  allocate(vector_field_display(NDIM,npoin))

  if(assign_external_model) then
    allocate(vpext(NGLLX,NGLLZ,nspec))
    allocate(vsext(NGLLX,NGLLZ,nspec))
    allocate(rhoext(NGLLX,NGLLZ,nspec))
  else
    allocate(vpext(1,1,1))
    allocate(vsext(1,1,1))
    allocate(rhoext(1,1,1))
  endif

  endif

!
!----  set the coordinates of the points of the global grid
!
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX

        xi = xigll(i)
        gamma = zigll(j)

        call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl,jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo)

        coord(1,ibool(i,j,ispec)) = x
        coord(2,ibool(i,j,ispec)) = z

        xix(i,j,ispec) = xixl
        xiz(i,j,ispec) = xizl
        gammax(i,j,ispec) = gammaxl
        gammaz(i,j,ispec) = gammazl
        jacobian(i,j,ispec) = jacobianl

      enddo
    enddo
  enddo

!
!--- save the grid of points in a file
!
  if(outputgrid .and. myrank == 0 .and. ipass == 1) then
     write(IOUT,*)
     write(IOUT,*) 'Saving the grid in a text file...'
     write(IOUT,*)
     open(unit=55,file='OUTPUT_FILES/grid_points_and_model.txt',status='unknown')
     zmax=maxval(coord(2,:))
     write(55,*) npoin
     do n = 1,npoin
        write(55,*) (coord(i,n), i=1,NDIM)
     enddo
     close(55)
  endif

!
!-----   plot the GLL mesh in a Gnuplot file
!
  if(gnuplot .and. myrank == 0 .and. ipass == 1) call plotgll(knods,ibool,coorg,coord,npoin,npgeo,ngnod,nspec)

!
!----  assign external velocity and density model if needed
!
  if(assign_external_model) then
    if (myrank == 0 .and. ipass == 1) then
      write(IOUT,*)
      write(IOUT,*) 'Assigning external velocity and density model...'
      write(IOUT,*)
    endif
    if(TURN_ANISOTROPY_ON .or. TURN_ATTENUATION_ON) &
         call exit_MPI('cannot have anisotropy nor attenuation if external model in current version')
    any_acoustic = .false.
    any_elastic = .false.
    do ispec = 1,nspec
      previous_vsext = -1.d0
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          call define_external_model(coord(1,iglob),coord(2,iglob),kmato(ispec),myrank,  &
                                         rhoext(i,j,ispec),vpext(i,j,ispec),vsext(i,j,ispec),myrank)
! stop if the same element is assigned both acoustic and elastic points in external model
          if(.not. (i == 1 .and. j == 1) .and. &
            ((vsext(i,j,ispec) >= TINYVAL .and. previous_vsext < TINYVAL) .or. &
             (vsext(i,j,ispec) < TINYVAL .and. previous_vsext >= TINYVAL)))  &
                call exit_MPI('external velocity model cannot be both fluid and solid inside the same spectral element')
          if(vsext(i,j,ispec) < TINYVAL) then
            elastic(ispec) = .false.
            poroelastic(ispec) = .false.
            any_acoustic = .true.
          else
            elastic(ispec) = .true.
            poroelastic(ispec) = .false.
            any_elastic = .true.
          endif
          previous_vsext = vsext(i,j,ispec)
        enddo
      enddo
    enddo
  endif

!
!----  perform basic checks on parameters read
!
  any_elastic_glob = any_elastic
#ifdef USE_MPI
  call MPI_ALLREDUCE(any_elastic, any_elastic_glob, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ier)
#endif

  any_poroelastic_glob = any_poroelastic
#ifdef USE_MPI
  call MPI_ALLREDUCE(any_poroelastic, any_poroelastic_glob, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ier)
#endif

  any_acoustic_glob = any_acoustic
#ifdef USE_MPI
  call MPI_ALLREDUCE(any_acoustic, any_acoustic_glob, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ier)
#endif

! for acoustic
  if(TURN_ANISOTROPY_ON .and. .not. any_elastic_glob) &
    call exit_MPI('cannot have anisotropy if acoustic/poroelastic simulation only')

  if(TURN_ATTENUATION_ON .and. .not. any_elastic_glob) &
    call exit_MPI('currently cannot have attenuation if acoustic/poroelastic simulation only')

! for attenuation
  if(TURN_ANISOTROPY_ON .and. TURN_ATTENUATION_ON) &
    call exit_MPI('cannot have anisotropy and attenuation both turned on in current version')

!
!----   define coefficients of the Newmark time scheme
!
  deltatover2 = HALF*deltat
  deltatsquareover2 = HALF*deltat*deltat

  if(isolver == 2) then
!  define coefficients of the Newmark time scheme for the backward wavefield
  b_deltat = - deltat
  b_deltatover2 = HALF*b_deltat
  b_deltatsquareover2 = HALF*b_deltat*b_deltat
  endif

!---- define actual location of source and receivers
  do i_source=1,NSOURCE

  if(source_type(i_source) == 1) then

! collocated force source
    call locate_source_force(coord,ibool,npoin,nspec,x_source(i_source),z_source(i_source), &
      ix_source(i_source),iz_source(i_source),ispec_selected_source(i_source),iglob_source(i_source), &
      is_proc_source(i_source),nb_proc_source(i_source),ipass)

! get density at the source in order to implement collocated force with the right
! amplitude later
    if(is_proc_source(i_source) == 1) then
      rho_at_source_location  = density(1,kmato(ispec_selected_source(i_source)))
! external velocity model
      if(assign_external_model) rho_at_source_location = &
          rhoext(ix_source(i_source),iz_source(i_source),ispec_selected_source(i_source))
    endif

! check that acoustic source is not exactly on the free surface because pressure is zero there
    if(is_proc_source(i_source) == 1) then
       do ispec_acoustic_surface = 1,nelem_acoustic_surface
          ispec = acoustic_surface(1,ispec_acoustic_surface)
          if( .not. elastic(ispec) .and. .not. poroelastic(ispec) .and. ispec == ispec_selected_source(i_source) ) then
             do j = acoustic_surface(4,ispec_acoustic_surface), acoustic_surface(5,ispec_acoustic_surface)
                do i = acoustic_surface(2,ispec_acoustic_surface), acoustic_surface(3,ispec_acoustic_surface)
                   iglob = ibool(i,j,ispec)
                   if ( iglob_source(i_source) == iglob ) then
 call exit_MPI('an acoustic source cannot be located exactly on the free surface because pressure is zero there')
                   endif
                enddo
             enddo
          endif
       enddo
    endif

  else if(source_type(i_source) == 2) then
! moment-tensor source
     call locate_source_moment_tensor(ibool,coord,nspec,npoin,xigll,zigll,x_source(i_source),z_source(i_source), &
          ispec_selected_source(i_source),is_proc_source(i_source),nb_proc_source(i_source),&
          nproc,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo,ipass)

! compute source array for moment-tensor source
    call compute_arrays_source(ispec_selected_source(i_source),xi_source(i_source),gamma_source(i_source),&
         sourcearray(i_source,:,:,:), &
         Mxx(i_source),Mzz(i_source),Mxz(i_source),xix,xiz,gammax,gammaz,xigll,zigll,nspec)

  else if(.not.initialfield) then
    call exit_MPI('incorrect source type')
  endif


! locate receivers in the mesh
  call locate_receivers(ibool,coord,nspec,npoin,xigll,zigll,nrec,nrecloc,recloc,which_proc_receiver,nproc,myrank,&
       st_xval,st_zval,ispec_selected_rec, &
       xi_receiver,gamma_receiver,station_name,network_name,x_source(i_source),z_source(i_source),coorg,knods,ngnod,npgeo,ipass, &
       x_final_receiver, z_final_receiver)

  enddo ! do i_source=1,NSOURCE

! compute source array for adjoint source
  if(isolver == 2) then  ! adjoint calculation
    nadj_rec_local = 0
    do irec = 1,nrec
      if(myrank == which_proc_receiver(irec))then
!   check that the source proc number is okay
        if(which_proc_receiver(irec) < 0 .or. which_proc_receiver(irec) > NPROC-1) &
              call exit_MPI(myrank,'something is wrong with the source proc number in adjoint simulation')
        nadj_rec_local = nadj_rec_local + 1
      endif
    enddo
  if(ipass == 1) allocate(adj_sourcearray(NSTEP,NDIM,NGLLX,NGLLZ))
  if (nadj_rec_local > 0 .and. ipass == 1)  allocate(adj_sourcearrays(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLZ))
    irec_local = 0
    do irec = 1, nrec
!   compute only adjoint source arrays in the local proc
      if(myrank == which_proc_receiver(irec))then
        irec_local = irec_local + 1
  adj_source_file = trim(station_name(irec))//'.'//trim(network_name(irec))
  call compute_arrays_adj_source(myrank,adj_source_file, &
              xi_receiver(irec), gamma_receiver(irec), &
              adj_sourcearray, xigll,zigll,NSTEP)
        adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
      endif
    enddo
  endif


if (ipass == 1) then
  if (nrecloc > 0) then
    allocate(anglerec_irec(nrecloc))
    allocate(cosrot_irec(nrecloc))
    allocate(sinrot_irec(nrecloc))
    allocate(rec_tangential_detection_curve(nrecloc))
  else
    allocate(anglerec_irec(1))
    allocate(cosrot_irec(1))
    allocate(sinrot_irec(1))
    allocate(rec_tangential_detection_curve(1))
  endif
  anglerec_irec(:) = anglerec * pi / 180.d0
  cosrot_irec(:) = cos(anglerec)
  sinrot_irec(:) = sin(anglerec)
endif

!
!--- tangential computation
!
if (ipass == NUMBER_OF_PASSES) then

! for receivers
  if (rec_normal_to_surface) then
    irecloc = 0
    do irec = 1, nrec
      if (which_proc_receiver(irec) == myrank) then
        irecloc = irecloc + 1
        distmin = HUGEVAL
        do i = 1, nnodes_tangential_curve
          dist_current = sqrt((x_final_receiver(irec)-nodes_tangential_curve(1,i))**2 + &
             (z_final_receiver(irec)-nodes_tangential_curve(2,i))**2)
          if ( dist_current < distmin ) then
            n1_tangential_detection_curve = i
            distmin = dist_current
          endif
       enddo

       rec_tangential_detection_curve(irecloc) = n1_tangential_detection_curve
       call tri_quad(n_tangential_detection_curve, n1_tangential_detection_curve, nnodes_tangential_curve)

       call calcul_normale( anglerec_irec(irecloc), nodes_tangential_curve(1,n_tangential_detection_curve(1)), &
         nodes_tangential_curve(1,n_tangential_detection_curve(2)), &
         nodes_tangential_curve(1,n_tangential_detection_curve(3)), &
         nodes_tangential_curve(1,n_tangential_detection_curve(4)), &
         nodes_tangential_curve(2,n_tangential_detection_curve(1)), &
         nodes_tangential_curve(2,n_tangential_detection_curve(2)), &
         nodes_tangential_curve(2,n_tangential_detection_curve(3)), &
         nodes_tangential_curve(2,n_tangential_detection_curve(4)) )
     endif

    enddo
    cosrot_irec(:) = cos(anglerec_irec(:))
    sinrot_irec(:) = sin(anglerec_irec(:))
  endif
! for the source
  if (force_normal_to_surface) then

    do i_source=1,NSOURCE
     if (is_proc_source(i_source) == 1) then
    distmin = HUGEVAL
    do i = 1, nnodes_tangential_curve
      dist_current = sqrt((coord(1,iglob_source(i_source))-nodes_tangential_curve(1,i))**2 + &
          (coord(2,iglob_source(i_source))-nodes_tangential_curve(2,i))**2)
      if ( dist_current < distmin ) then
        n1_tangential_detection_curve = i
        distmin = dist_current

      endif
    enddo

    call tri_quad(n_tangential_detection_curve, n1_tangential_detection_curve, nnodes_tangential_curve)

    call calcul_normale( angleforce(i_source), nodes_tangential_curve(1,n_tangential_detection_curve(1)), &
      nodes_tangential_curve(1,n_tangential_detection_curve(2)), &
      nodes_tangential_curve(1,n_tangential_detection_curve(3)), &
      nodes_tangential_curve(1,n_tangential_detection_curve(4)), &
      nodes_tangential_curve(2,n_tangential_detection_curve(1)), &
      nodes_tangential_curve(2,n_tangential_detection_curve(2)), &
      nodes_tangential_curve(2,n_tangential_detection_curve(3)), &
      nodes_tangential_curve(2,n_tangential_detection_curve(4)) )
   
    source_courbe_eros(i_source) = n1_tangential_detection_curve
    if ( myrank == 0 .and. is_proc_source(i_source) == 1 .and. nb_proc_source(i_source) == 1 ) then
      source_courbe_eros(i_source) = n1_tangential_detection_curve
      angleforce_recv = angleforce(i_source)
#ifdef USE_MPI
    else if ( myrank == 0 ) then
      do i = 1, nb_proc_source(i_source) - is_proc_source(i_source)
        call MPI_recv(source_courbe_eros(i_source),1,MPI_INTEGER,MPI_ANY_SOURCE,42,MPI_COMM_WORLD,request_mpi_status,ier)
        call MPI_recv(angleforce_recv,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,43,MPI_COMM_WORLD,request_mpi_status,ier)
      enddo
    else if ( is_proc_source(i_source) == 1 ) then
      call MPI_send(n1_tangential_detection_curve,1,MPI_INTEGER,0,42,MPI_COMM_WORLD,ier)
      call MPI_send(angleforce(i_source),1,MPI_DOUBLE_PRECISION,0,43,MPI_COMM_WORLD,ier)
#endif
    endif

#ifdef USE_MPI
    call MPI_bcast(angleforce_recv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    angleforce(i_source) = angleforce_recv
#endif
      endif !  if (is_proc_source(i_source) == 1)
     enddo ! do i_source=1,NSOURCE
  endif !  if (force_normal_to_surface)

! CHRIS --- how to deal with multiple source. Use first source now. ---
! compute distance from source to receivers following the curve
  if (force_normal_to_surface .and. rec_normal_to_surface) then
    dist_tangential_detection_curve(source_courbe_eros(1)) = 0
    do i = source_courbe_eros(1)+1, nnodes_tangential_curve
      dist_tangential_detection_curve(i) = dist_tangential_detection_curve(i-1) + &
          sqrt((nodes_tangential_curve(1,i)-nodes_tangential_curve(1,i-1))**2 + &
          (nodes_tangential_curve(2,i)-nodes_tangential_curve(2,i-1))**2)
    enddo
    dist_tangential_detection_curve(1) = dist_tangential_detection_curve(nnodes_tangential_curve) + &
         sqrt((nodes_tangential_curve(1,1)-nodes_tangential_curve(1,nnodes_tangential_curve))**2 + &
         (nodes_tangential_curve(2,1)-nodes_tangential_curve(2,nnodes_tangential_curve))**2)
    do i = 2, source_courbe_eros(1)-1
      dist_tangential_detection_curve(i) = dist_tangential_detection_curve(i-1) + &
          sqrt((nodes_tangential_curve(1,i)-nodes_tangential_curve(1,i-1))**2 + &
          (nodes_tangential_curve(2,i)-nodes_tangential_curve(2,i-1))**2)
    enddo
    do i = source_courbe_eros(1)-1, 1, -1
      dist_current = dist_tangential_detection_curve(i+1) + &
          sqrt((nodes_tangential_curve(1,i)-nodes_tangential_curve(1,i+1))**2 + &
          (nodes_tangential_curve(2,i)-nodes_tangential_curve(2,i+1))**2)
      if ( dist_current < dist_tangential_detection_curve(i) ) then
        dist_tangential_detection_curve(i) = dist_current
      endif
    enddo
    dist_current = dist_tangential_detection_curve(1) + &
       sqrt((nodes_tangential_curve(1,1)-nodes_tangential_curve(1,nnodes_tangential_curve))**2 + &
       (nodes_tangential_curve(2,1)-nodes_tangential_curve(2,nnodes_tangential_curve))**2)
    if ( dist_current < dist_tangential_detection_curve(nnodes_tangential_curve) ) then
      dist_tangential_detection_curve(nnodes_tangential_curve) = dist_current
    endif
    do i = nnodes_tangential_curve-1, source_courbe_eros(1)+1, -1
      dist_current = dist_tangential_detection_curve(i+1) + &
          sqrt((nodes_tangential_curve(1,i)-nodes_tangential_curve(1,i+1))**2 + &
          (nodes_tangential_curve(2,i)-nodes_tangential_curve(2,i+1))**2)
      if ( dist_current < dist_tangential_detection_curve(i) ) then
        dist_tangential_detection_curve(i) = dist_current
      endif
    enddo

    if ( myrank == 0 ) then
      open(unit=11,file='OUTPUT_FILES/dist_rec_tangential_detection_curve', form='formatted', status='unknown')
    endif
      irecloc = 0
    do irec = 1,nrec

      if ( myrank == 0 ) then
        if ( which_proc_receiver(irec) == myrank ) then
          irecloc = irecloc + 1
          n1_tangential_detection_curve = rec_tangential_detection_curve(irecloc)
          x_final_receiver_dummy = x_final_receiver(irec)
          z_final_receiver_dummy = z_final_receiver(irec)
#ifdef USE_MPI
        else

          call MPI_RECV(n1_tangential_detection_curve,1,MPI_INTEGER,&
             which_proc_receiver(irec),irec,MPI_COMM_WORLD,request_mpi_status,ier)
          call MPI_RECV(x_final_receiver_dummy,1,MPI_DOUBLE_PRECISION,&
             which_proc_receiver(irec),irec,MPI_COMM_WORLD,request_mpi_status,ier)
          call MPI_RECV(z_final_receiver_dummy,1,MPI_DOUBLE_PRECISION,&
             which_proc_receiver(irec),irec,MPI_COMM_WORLD,request_mpi_status,ier)

#endif
        endif

#ifdef USE_MPI
      else
        if ( which_proc_receiver(irec) == myrank ) then
          irecloc = irecloc + 1
          call MPI_SEND(rec_tangential_detection_curve(irecloc),1,MPI_INTEGER,0,irec,MPI_COMM_WORLD,ier)
          call MPI_SEND(x_final_receiver(irec),1,MPI_DOUBLE_PRECISION,0,irec,MPI_COMM_WORLD,ier)
          call MPI_SEND(z_final_receiver(irec),1,MPI_DOUBLE_PRECISION,0,irec,MPI_COMM_WORLD,ier)
        endif
#endif

      endif
      if ( myrank == 0 ) then
        write(11,*) dist_tangential_detection_curve(n1_tangential_detection_curve)
        write(12,*) x_final_receiver_dummy
        write(13,*) z_final_receiver_dummy
      endif
    enddo

    if ( myrank == 0 ) then
      close(11)
      close(12)
      close(13)
    endif

  endif
endif

!
!---
!

! allocate seismogram arrays
  if(ipass == 1) then
    allocate(sisux(NTSTEP_BETWEEN_OUTPUT_SEISMO,nrecloc))
    allocate(sisuz(NTSTEP_BETWEEN_OUTPUT_SEISMO,nrecloc))
    allocate(siscurl(NTSTEP_BETWEEN_OUTPUT_SEISMO,nrecloc))
  endif

! check if acoustic receiver is exactly on the free surface because pressure is zero there
  do ispec_acoustic_surface = 1,nelem_acoustic_surface
     ispec = acoustic_surface(1,ispec_acoustic_surface)
     ixmin = acoustic_surface(2,ispec_acoustic_surface)
     ixmax = acoustic_surface(3,ispec_acoustic_surface)
     izmin = acoustic_surface(4,ispec_acoustic_surface)
     izmax = acoustic_surface(5,ispec_acoustic_surface)
     do irecloc = 1,nrecloc
        irec = recloc(irecloc)
        if(.not. elastic(ispec) .and. .not. poroelastic(ispec) .and. ispec == ispec_selected_rec(irec)) then
           if ( (izmin==1 .and. izmax==1 .and. ixmin==1 .and. ixmax==NGLLX .and. &
                gamma_receiver(irec) < -0.99d0) .or.&
                (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==NGLLX .and. &
                gamma_receiver(irec) > 0.99d0) .or.&
                (izmin==1 .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==1 .and. &
                xi_receiver(irec) < -0.99d0) .or.&
                (izmin==1 .and. izmax==NGLLZ .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
                xi_receiver(irec) > 0.99d0) .or.&
                (izmin==1 .and. izmax==1 .and. ixmin==1 .and. ixmax==1 .and. &
                gamma_receiver(irec) < -0.99d0 .and. xi_receiver(irec) < -0.99d0) .or.&
                (izmin==1 .and. izmax==1 .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
                gamma_receiver(irec) < -0.99d0 .and. xi_receiver(irec) > 0.99d0) .or.&
                (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==1 .and. &
                gamma_receiver(irec) > 0.99d0 .and. xi_receiver(irec) < -0.99d0) .or.&
                (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
                gamma_receiver(irec) > 0.99d0 .and. xi_receiver(irec) > 0.99d0) ) then
              if(seismotype == 4) then
call exit_MPI('an acoustic pressure receiver cannot be located exactly on the free surface because pressure is zero there')
              else
                 print *, '**********************************************************************'
                 print *, '*** Warning: acoustic receiver located exactly on the free surface ***'
                 print *, '*** Warning: tangential component will be zero there               ***'
                 print *, '**********************************************************************'
                 print *
              endif
           endif
        endif
     enddo
  enddo

! define and store Lagrange interpolators at all the receivers
  do irec = 1,nrec
    call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
    call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
    hxir_store(irec,:) = hxir(:)
    hgammar_store(irec,:) = hgammar(:)
  enddo

! displacement, velocity, acceleration and inverse of the mass matrix for elastic elements
  if(ipass == 1) then

  if(any_elastic) then
    allocate(displ_elastic(NDIM,npoin))
    allocate(veloc_elastic(NDIM,npoin))
    allocate(accel_elastic(NDIM,npoin))
    allocate(rmass_inverse_elastic(npoin))
  else
! allocate unused arrays with fictitious size
    allocate(displ_elastic(1,1))
    allocate(veloc_elastic(1,1))
    allocate(accel_elastic(1,1))
    allocate(rmass_inverse_elastic(1))
  endif
! extra array if adjoint and kernels calculation
  if(isolver == 2 .and. any_elastic) then
    allocate(b_displ_elastic(NDIM,npoin))
    allocate(b_veloc_elastic(NDIM,npoin))
    allocate(b_accel_elastic(NDIM,npoin))
    allocate(rho_kl(npoin))
    allocate(rho_k(npoin))
    allocate(rhol_global(npoin))
    allocate(mu_kl(npoin))
    allocate(mu_k(npoin))
    allocate(mul_global(npoin))
    allocate(kappa_kl(npoin))
    allocate(kappa_k(npoin))
    allocate(kappal_global(npoin))
    allocate(rhop_kl(npoin))
    allocate(alpha_kl(npoin))
    allocate(beta_kl(npoin))
  else
    allocate(b_displ_elastic(1,1))
    allocate(b_veloc_elastic(1,1))
    allocate(b_accel_elastic(1,1))
    allocate(rho_kl(1))
    allocate(rho_k(1))
    allocate(rhol_global(1))
    allocate(mu_kl(1))
    allocate(mu_k(1))
    allocate(mul_global(1))
    allocate(kappa_kl(1))
    allocate(kappa_k(1))
    allocate(kappal_global(1))
    allocate(rhop_kl(1))
    allocate(alpha_kl(1))
    allocate(beta_kl(1))
  endif

  if(any_poroelastic) then
    allocate(displs_poroelastic(NDIM,npoin))
    allocate(velocs_poroelastic(NDIM,npoin))
    allocate(accels_poroelastic(NDIM,npoin))
    allocate(rmass_s_inverse_poroelastic(npoin))
    allocate(displw_poroelastic(NDIM,npoin))
    allocate(velocw_poroelastic(NDIM,npoin))
    allocate(accelw_poroelastic(NDIM,npoin))
    allocate(rmass_w_inverse_poroelastic(npoin))
  else
! allocate unused arrays with fictitious size
    allocate(displs_poroelastic(1,1))
    allocate(velocs_poroelastic(1,1))
    allocate(accels_poroelastic(1,1))
    allocate(rmass_s_inverse_poroelastic(1))
    allocate(displw_poroelastic(1,1))
    allocate(velocw_poroelastic(1,1))
    allocate(accelw_poroelastic(1,1))
    allocate(rmass_w_inverse_poroelastic(1))
  endif
! extra array if adjoint and kernels calculation
  if(isolver == 2 .and. any_poroelastic) then
    allocate(b_displs_poroelastic(NDIM,npoin))
    allocate(b_velocs_poroelastic(NDIM,npoin))
    allocate(b_accels_poroelastic(NDIM,npoin))
    allocate(b_displw_poroelastic(NDIM,npoin))
    allocate(b_velocw_poroelastic(NDIM,npoin))
    allocate(b_accelw_poroelastic(NDIM,npoin))
    allocate(rhot_kl(npoin))
    allocate(rhot_k(npoin))
    allocate(rhof_kl(npoin))
    allocate(rhof_k(npoin))
    allocate(sm_kl(npoin))
    allocate(sm_k(npoin))
    allocate(eta_kl(npoin))
    allocate(eta_k(npoin))
    allocate(mufr_kl(npoin))
    allocate(mufr_k(npoin))
    allocate(B_kl(npoin))
    allocate(B_k(npoin))
    allocate(C_kl(npoin))
    allocate(C_k(npoin))
    allocate(M_kl(npoin))
    allocate(M_k(npoin))
    allocate(rhob_kl(npoin))
    allocate(rhofb_kl(npoin))
    allocate(phi_kl(npoin))
    allocate(Bb_kl(npoin))
    allocate(Cb_kl(npoin))
    allocate(Mb_kl(npoin))
    allocate(mufrb_kl(npoin))
    allocate(rhobb_kl(npoin))
    allocate(rhofbb_kl(npoin))
    allocate(phib_kl(npoin))
    allocate(cpI_kl(npoin))
    allocate(cpII_kl(npoin))
    allocate(cs_kl(npoin))
    allocate(ratio_kl(npoin))
    allocate(phil_global(npoin))
    allocate(mulfr_global(npoin))
    allocate(etal_f_global(npoin))
    allocate(rhol_s_global(npoin))
    allocate(rhol_f_global(npoin))
    allocate(rhol_bar_global(npoin))
    allocate(tortl_global(npoin))
    allocate(permlxx_global(npoin))
    allocate(permlxz_global(npoin))
    allocate(permlzz_global(npoin))
  else
    allocate(b_displs_poroelastic(1,1))
    allocate(b_velocs_poroelastic(1,1))
    allocate(b_accels_poroelastic(1,1))
    allocate(b_displw_poroelastic(1,1))
    allocate(b_velocw_poroelastic(1,1))
    allocate(b_accelw_poroelastic(1,1))
    allocate(rhot_kl(1))
    allocate(rhot_k(1))
    allocate(rhof_kl(1))
    allocate(rhof_k(1))
    allocate(sm_kl(1))
    allocate(sm_k(1))
    allocate(eta_kl(1))
    allocate(eta_k(1))
    allocate(mufr_kl(1))
    allocate(mufr_k(1))
    allocate(B_kl(1))
    allocate(B_k(1))
    allocate(C_kl(1))
    allocate(C_k(1))
    allocate(M_kl(1))
    allocate(M_k(1))
    allocate(rhob_kl(1))
    allocate(rhofb_kl(1))
    allocate(phi_kl(1))
    allocate(Bb_kl(1))
    allocate(Cb_kl(1))
    allocate(Mb_kl(1))
    allocate(mufrb_kl(1))
    allocate(rhobb_kl(1))
    allocate(rhofbb_kl(1))
    allocate(phib_kl(1))
    allocate(cpI_kl(1))
    allocate(cpII_kl(1))
    allocate(cs_kl(1))
    allocate(ratio_kl(1))
    allocate(phil_global(1))
    allocate(mulfr_global(1))
    allocate(etal_f_global(1))
    allocate(rhol_s_global(1))
    allocate(rhol_f_global(1))
    allocate(rhol_bar_global(1))
    allocate(tortl_global(1))
    allocate(permlxx_global(1))
    allocate(permlxz_global(1))
    allocate(permlzz_global(1))
  endif

  if(any_poroelastic .and. any_elastic) then
    allocate(displ(2,npoin))
    allocate(veloc(2,npoin))
  else
    allocate(displ(2,1))
    allocate(veloc(2,1))
  endif

! potential, its first and second derivative, and inverse of the mass matrix for acoustic elements
  if(any_acoustic) then
    allocate(potential_acoustic(npoin))
    allocate(potential_dot_acoustic(npoin))
    allocate(potential_dot_dot_acoustic(npoin))
    allocate(rmass_inverse_acoustic(npoin))
  else
! allocate unused arrays with fictitious size
    allocate(potential_acoustic(1))
    allocate(potential_dot_acoustic(1))
    allocate(potential_dot_dot_acoustic(1))
    allocate(rmass_inverse_acoustic(1))
  endif
  if(isolver == 2 .and. any_acoustic) then
    allocate(b_potential_acoustic(npoin))
    allocate(b_potential_dot_acoustic(npoin))
    allocate(b_potential_dot_dot_acoustic(npoin))
    allocate(rho_ac_kl(npoin))
    allocate(rho_ac_k(npoin))
    allocate(rhol_ac_global(npoin))
    allocate(kappa_ac_kl(npoin))
    allocate(kappa_ac_k(npoin))
    allocate(kappal_ac_global(npoin))
    allocate(rhop_ac_kl(npoin))
    allocate(alpha_ac_kl(npoin))
  else
! allocate unused arrays with fictitious size
    allocate(b_potential_acoustic(1))
    allocate(b_potential_dot_acoustic(1))
    allocate(b_potential_dot_dot_acoustic(1))
    allocate(rho_ac_kl(1))
    allocate(rho_ac_k(1))
    allocate(rhol_ac_global(1))
    allocate(kappa_ac_kl(1))
    allocate(kappa_ac_k(1))
    allocate(kappal_ac_global(1))
    allocate(rhop_ac_kl(1))
    allocate(alpha_ac_kl(1))
  endif

  endif

!
!---- build the global mass matrix and invert it once and for all
!
  if(any_elastic) rmass_inverse_elastic(:) = ZERO
  if(any_poroelastic) rmass_s_inverse_poroelastic(:) = ZERO
  if(any_poroelastic) rmass_w_inverse_poroelastic(:) = ZERO
  if(any_acoustic) rmass_inverse_acoustic(:) = ZERO
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)

! if external density model (elastic or acoustic)
        if(assign_external_model) then
          rhol = rhoext(i,j,ispec)
          kappal = rhol * vpext(i,j,ispec)**2
        else
          rhol = density(1,kmato(ispec))
          lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
          mul_relaxed = poroelastcoef(2,1,kmato(ispec))
          kappal = lambdal_relaxed + 2.d0*mul_relaxed
        endif

        if(poroelastic(ispec)) then     ! material is poroelastic
          rhol_s = density(1,kmato(ispec))
          rhol_f = density(2,kmato(ispec))
          phil = porosity(kmato(ispec))
          tortl = tortuosity(kmato(ispec))
          rhol_bar = (1._CUSTOM_REAL-phil)*rhol_s + phil*rhol_f
! for the solid mass matrix
             rmass_s_inverse_poroelastic(iglob) = rmass_s_inverse_poroelastic(iglob) + &
       wxgll(i)*wzgll(j)*jacobian(i,j,ispec)*(rhol_bar - phil*rhol_f/tortl)
! for the fluid mass matrix
             rmass_w_inverse_poroelastic(iglob) = rmass_w_inverse_poroelastic(iglob) + &
      wxgll(i)*wzgll(j)*jacobian(i,j,ispec)*(rhol_bar*rhol_f*tortl - &
       phil*rhol_f*rhol_f)/(rhol_bar*phil)
        elseif(elastic(ispec)) then    ! for elastic medium
          rmass_inverse_elastic(iglob) = rmass_inverse_elastic(iglob) + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)
        else                           ! for acoustic medium
          rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) + wxgll(i)*wzgll(j)*jacobian(i,j,ispec) / kappal
        endif

      enddo
    enddo
  enddo ! do ispec = 1,nspec

#ifdef USE_MPI
  if ( nproc > 1 ) then
! preparing for MPI communications
    if(ipass == 1) allocate(mask_ispec_inner_outer(nspec))
    mask_ispec_inner_outer(:) = .false.

    call prepare_assemble_MPI (nspec,ibool,knods,ngnod,npoin,elastic,poroelastic, &
          ninterface, max_interface_size,my_nelmnts_neighbours, my_interfaces, &
          ibool_interfaces_acoustic, ibool_interfaces_elastic, ibool_interfaces_poroelastic, &
          nibool_interfaces_acoustic, nibool_interfaces_elastic, nibool_interfaces_poroelastic, &
          inum_interfaces_acoustic, inum_interfaces_elastic, inum_interfaces_poroelastic, &
          ninterface_acoustic, ninterface_elastic, ninterface_poroelastic,mask_ispec_inner_outer)

    nspec_outer = count(mask_ispec_inner_outer)
    nspec_inner = nspec - nspec_outer

    if(ipass == 1) then
      allocate(ispec_outer_to_glob(nspec_outer))
      allocate(ispec_inner_to_glob(nspec_inner))
    endif

! building of corresponding arrays between inner/outer elements and their global number
if(ipass == 1) then
    num_ispec_outer = 0
    num_ispec_inner = 0
    do ispec = 1, nspec
      if ( mask_ispec_inner_outer(ispec) ) then
        num_ispec_outer = num_ispec_outer + 1
        ispec_outer_to_glob(num_ispec_outer) = ispec
      else
        num_ispec_inner = num_ispec_inner + 1
        ispec_inner_to_glob(num_ispec_inner) = ispec
      endif
    enddo
endif

  max_ibool_interfaces_size_ac = maxval(nibool_interfaces_acoustic(:))
  max_ibool_interfaces_size_el = NDIM*maxval(nibool_interfaces_elastic(:))
  max_ibool_interfaces_size_po = NDIM*maxval(nibool_interfaces_poroelastic(:))
  if(ipass == 1) then
    allocate(tab_requests_send_recv_acoustic(ninterface_acoustic*2))
    allocate(buffer_send_faces_vector_ac(max_ibool_interfaces_size_ac,ninterface_acoustic))
    allocate(buffer_recv_faces_vector_ac(max_ibool_interfaces_size_ac,ninterface_acoustic))
    allocate(tab_requests_send_recv_elastic(ninterface_elastic*2))
    allocate(buffer_send_faces_vector_el(max_ibool_interfaces_size_el,ninterface_elastic))
    allocate(buffer_recv_faces_vector_el(max_ibool_interfaces_size_el,ninterface_elastic))
    allocate(tab_requests_send_recv_poroelastic(ninterface_poroelastic*4))
    allocate(buffer_send_faces_vector_pos(max_ibool_interfaces_size_po,ninterface_poroelastic))
    allocate(buffer_recv_faces_vector_pos(max_ibool_interfaces_size_po,ninterface_poroelastic))
    allocate(buffer_send_faces_vector_pow(max_ibool_interfaces_size_po,ninterface_poroelastic))
    allocate(buffer_recv_faces_vector_pow(max_ibool_interfaces_size_po,ninterface_poroelastic))
  endif

! assembling the mass matrix
  call assemble_MPI_scalar(rmass_inverse_acoustic,rmass_inverse_elastic,rmass_s_inverse_poroelastic, &
     rmass_w_inverse_poroelastic,npoin, &
     ninterface, max_interface_size, max_ibool_interfaces_size_ac, max_ibool_interfaces_size_el, &
     max_ibool_interfaces_size_po, &
     ibool_interfaces_acoustic,ibool_interfaces_elastic,ibool_interfaces_poroelastic, &
     nibool_interfaces_acoustic,nibool_interfaces_elastic,nibool_interfaces_poroelastic,my_neighbours)

  else
    ninterface_acoustic = 0
    ninterface_elastic = 0
    ninterface_poroelastic = 0

    num_ispec_outer = 0
    num_ispec_inner = 0
    if(ipass == 1) allocate(mask_ispec_inner_outer(1))

    nspec_outer = 0
    nspec_inner = nspec

    if(ipass == 1) allocate(ispec_inner_to_glob(nspec_inner))
    do ispec = 1, nspec
      ispec_inner_to_glob(ispec) = ispec
    enddo

  endif ! end of test on wether there is more than one process (nproc > 1)

#else
  num_ispec_outer = 0
  num_ispec_inner = 0
  if(ipass == 1) allocate(mask_ispec_inner_outer(1))

  nspec_outer = 0
  nspec_inner = nspec

  if(ipass == 1) then
    allocate(ispec_outer_to_glob(1))
    allocate(ispec_inner_to_glob(nspec_inner))
  endif
  do ispec = 1, nspec
     ispec_inner_to_glob(ispec) = ispec
  enddo

#endif

if(ipass == 1) then

  allocate(antecedent_list(nspec))

! loop over spectral elements
  do ispec_outer = 1,nspec_outer
! get global numbering for inner or outer elements
    ispec = ispec_outer_to_glob(ispec_outer)
    antecedent_list(ispec) = ispec_outer
  enddo

! loop over spectral elements
  do ispec_inner = 1,nspec_inner
! get global numbering for inner or outer elements
    ispec = ispec_inner_to_glob(ispec_inner)
    antecedent_list(ispec) = nspec_outer + ispec_inner
  enddo

  allocate(ibool_outer(NGLLX,NGLLZ,nspec_outer))
  allocate(ibool_inner(NGLLX,NGLLZ,nspec_inner))

! loop over spectral elements
  do ispec_outer = 1,nspec_outer
! get global numbering for inner or outer elements
    ispec = ispec_outer_to_glob(ispec_outer)
    ibool_outer(:,:,ispec_outer) = ibool(:,:,ispec)
  enddo

! loop over spectral elements
  do ispec_inner = 1,nspec_inner
! get global numbering for inner or outer elements
    ispec = ispec_inner_to_glob(ispec_inner)
    ibool_inner(:,:,ispec_inner) = ibool(:,:,ispec)
  enddo

  allocate(copy_ibool_ori(NGLLX,NGLLZ,nspec_outer))
  allocate(mask_ibool(npoin))

  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:) = ibool_outer(:,:,:)

  inumber = 0

  do ispec = 1,nspec_outer
    do j=1,NGLLZ
      do i=1,NGLLX
        if(mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
! create a new point
          inumber = inumber + 1
          ibool_outer(i,j,ispec) = inumber
          mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
        else
! use an existing point created previously
          ibool_outer(i,j,ispec) = mask_ibool(copy_ibool_ori(i,j,ispec))
        endif
      enddo
    enddo
  enddo

  deallocate(copy_ibool_ori)
  deallocate(mask_ibool)

! the total number of points without multiples in this region is now known
  npoin_outer = maxval(ibool_outer)

  allocate(copy_ibool_ori(NGLLX,NGLLZ,nspec_inner))
  allocate(mask_ibool(npoin))

  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:) = ibool_inner(:,:,:)

  inumber = 0

  do ispec = 1,nspec_inner
    do j=1,NGLLZ
      do i=1,NGLLX
        if(mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
! create a new point
          inumber = inumber + 1
          ibool_inner(i,j,ispec) = inumber
          mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
        else
! use an existing point created previously
          ibool_inner(i,j,ispec) = mask_ibool(copy_ibool_ori(i,j,ispec))
        endif
      enddo
    enddo
  enddo

  deallocate(copy_ibool_ori)
  deallocate(mask_ibool)

! the total number of points without multiples in this region is now known
  npoin_inner = maxval(ibool_inner)

  allocate(perm(nspec))

! use identity permutation by default
  do ispec = 1,nspec
    perm(ispec) = ispec
  enddo

  if(ACTUALLY_IMPLEMENT_PERM_WHOLE) then

    allocate(check_perm(nspec))
    call get_perm(ibool,perm,LIMIT_MULTI_CUTHILL,nspec,npoin)
! check that the permutation obtained is bijective
    check_perm(:) = -1
    do ispec = 1,nspec
      check_perm(perm(ispec)) = ispec
    enddo
    if(minval(check_perm) /= 1) stop 'minval check_perm is incorrect for whole'
    if(maxval(check_perm) /= nspec) stop 'maxval check_perm is incorrect for whole'
    deallocate(check_perm)
  else

  if(ACTUALLY_IMPLEMENT_PERM_OUT) then
    allocate(check_perm(nspec_outer))
    call get_perm(ibool_outer,perm(1:nspec_outer),LIMIT_MULTI_CUTHILL,nspec_outer,npoin_outer)
! check that the permutation obtained is bijective
    check_perm(:) = -1
    do ispec = 1,nspec_outer
      check_perm(perm(ispec)) = ispec
    enddo
    if(minval(check_perm) /= 1) stop 'minval check_perm is incorrect for outer'
    if(maxval(check_perm) /= nspec_outer) stop 'maxval check_perm is incorrect for outer'
    deallocate(check_perm)
    deallocate(ibool_outer)
  endif

  if(ACTUALLY_IMPLEMENT_PERM_INN) then
    allocate(check_perm(nspec_inner))
    call get_perm(ibool_inner,perm(nspec_outer+1:nspec),LIMIT_MULTI_CUTHILL,nspec_inner,npoin_inner)
! check that the permutation obtained is bijective
    check_perm(:) = -1
    do ispec = 1,nspec_inner
      check_perm(perm(nspec_outer+ispec)) = ispec
    enddo
    if(minval(check_perm) /= 1) stop 'minval check_perm is incorrect for inner'
    if(maxval(check_perm) /= nspec_inner) stop 'maxval check_perm is incorrect for inner'
    deallocate(check_perm)
! add the right offset
    perm(nspec_outer+1:nspec) = perm(nspec_outer+1:nspec) + nspec_outer
    deallocate(ibool_inner)
  endif

  endif

endif

  enddo ! end of further reduction of cache misses inner/outer in two passes

!---
!---  end of section performed in two passes
!---

! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
  if(any_elastic) where(rmass_inverse_elastic <= 0._CUSTOM_REAL) rmass_inverse_elastic = 1._CUSTOM_REAL
  if(any_poroelastic) where(rmass_s_inverse_poroelastic <= 0._CUSTOM_REAL) rmass_s_inverse_poroelastic = 1._CUSTOM_REAL
  if(any_poroelastic) where(rmass_w_inverse_poroelastic <= 0._CUSTOM_REAL) rmass_w_inverse_poroelastic = 1._CUSTOM_REAL
  if(any_acoustic) where(rmass_inverse_acoustic <= 0._CUSTOM_REAL) rmass_inverse_acoustic = 1._CUSTOM_REAL

! compute the inverse of the mass matrix
  if(any_elastic) rmass_inverse_elastic(:) = 1._CUSTOM_REAL / rmass_inverse_elastic(:)
  if(any_poroelastic) rmass_s_inverse_poroelastic(:) = 1._CUSTOM_REAL / rmass_s_inverse_poroelastic(:)
  if(any_poroelastic) rmass_w_inverse_poroelastic(:) = 1._CUSTOM_REAL / rmass_w_inverse_poroelastic(:)
  if(any_acoustic) rmass_inverse_acoustic(:) = 1._CUSTOM_REAL / rmass_inverse_acoustic(:)

! check the mesh, stability and number of points per wavelength
  if(DISPLAY_SUBSET_OPTION == 1) then
    UPPER_LIMIT_DISPLAY = nspec
  else if(DISPLAY_SUBSET_OPTION == 2) then
    UPPER_LIMIT_DISPLAY = nspec_inner
  else if(DISPLAY_SUBSET_OPTION == 3) then
    UPPER_LIMIT_DISPLAY = nspec_outer
  else if(DISPLAY_SUBSET_OPTION == 4) then
    UPPER_LIMIT_DISPLAY = NSPEC_DISPLAY_SUBSET
  else
    stop 'incorrect value of DISPLAY_SUBSET_OPTION'
  endif
  call checkgrid(vpext,vsext,rhoext,density,poroelastcoef,porosity,tortuosity,ibool,kmato, &
                 coord,npoin,vpImin,vpImax,vpIImin,vpIImax, &
                 assign_external_model,nspec,UPPER_LIMIT_DISPLAY,numat,deltat,f0,t0,initialfield, &
                 time_function_type,coorg,xinterp,zinterp,shape2D_display,knods,simulation_title, &
                 npgeo,pointsdisp,ngnod,any_elastic,any_poroelastic,myrank,nproc,NSOURCE,poroelastic)

! convert receiver angle to radians
  anglerec = anglerec * pi / 180.d0

!
!---- for color images
!

  if(output_color_image) then

! horizontal size of the image
  xmin_color_image_loc = minval(coord(1,:))
  xmax_color_image_loc = maxval(coord(1,:))

! vertical size of the image, slightly increase it to go beyond maximum topography
  zmin_color_image_loc = minval(coord(2,:))
  zmax_color_image_loc = maxval(coord(2,:))

! global values
  xmin_color_image = xmin_color_image_loc
  xmax_color_image = xmax_color_image_loc
  zmin_color_image = zmin_color_image_loc
  zmax_color_image = zmax_color_image_loc
  npgeo_glob = npgeo

#ifdef USE_MPI
  call MPI_ALLREDUCE(xmin_color_image_loc, xmin_color_image, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE(xmax_color_image_loc, xmax_color_image, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE(zmin_color_image_loc, zmin_color_image, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE(zmax_color_image_loc, zmax_color_image, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE(npgeo, npgeo_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ier)

#endif

  zmax_color_image = zmin_color_image + 1.05d0 * (zmax_color_image - zmin_color_image)

! compute number of pixels in the horizontal direction based on typical number
! of spectral elements in a given direction (may give bad results for very elongated models)
  NX_IMAGE_color = nint(sqrt(dble(npgeo_glob))) * (NGLLX-1) + 1

! compute number of pixels in the vertical direction based on ratio of sizes
  NZ_IMAGE_color = nint(NX_IMAGE_color * (zmax_color_image - zmin_color_image) / (xmax_color_image - xmin_color_image))

! convert pixel sizes to even numbers because easier to reduce size, create MPEG movies in postprocessing
  NX_IMAGE_color = 2 * (NX_IMAGE_color / 2)
  NZ_IMAGE_color = 2 * (NZ_IMAGE_color / 2)

! check that image size is not too big
  if (NX_IMAGE_color > 99999) call exit_MPI('output image too big : NX_IMAGE_color > 99999.')
  if (NZ_IMAGE_color > 99999) call exit_MPI('output image too big : NZ_IMAGE_color > 99999.')

! allocate an array for image data
  allocate(image_color_data(NX_IMAGE_color,NZ_IMAGE_color))
  allocate(image_color_vp_display(NX_IMAGE_color,NZ_IMAGE_color))

! allocate an array for the grid point that corresponds to a given image data point
  allocate(iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color))
  allocate(copy_iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color))

! create all the pixels
  if (myrank == 0) then
    write(IOUT,*)
    write(IOUT,*) 'locating all the pixels of color images'
  endif

  size_pixel_horizontal = (xmax_color_image - xmin_color_image) / dble(NX_IMAGE_color-1)
  size_pixel_vertical = (zmax_color_image - zmin_color_image) / dble(NZ_IMAGE_color-1)

  iglob_image_color(:,:) = -1

! cheking which pixels are inside each elements
  nb_pixel_loc = 0
  do ispec = 1, nspec

     do k = 1, 4
        elmnt_coords(1,k) = coorg(1,knods(k,ispec))
        elmnt_coords(2,k) = coorg(2,knods(k,ispec))

     enddo

! avoid working on the whole pixel grid
     min_i = floor(minval((elmnt_coords(1,:) - xmin_color_image))/size_pixel_horizontal) + 1
     max_i = ceiling(maxval((elmnt_coords(1,:) - xmin_color_image))/size_pixel_horizontal) + 1
     min_j = floor(minval((elmnt_coords(2,:) - zmin_color_image))/size_pixel_vertical) + 1
     max_j = ceiling(maxval((elmnt_coords(2,:) - zmin_color_image))/size_pixel_vertical) + 1

     do j = min_j, max_j
        do i = min_i, max_i
           i_coord = (i-1)*size_pixel_horizontal + xmin_color_image
           j_coord = (j-1)*size_pixel_vertical + zmin_color_image

! checking if the pixel is inside the element (must be a convex quadrilateral)
           call is_in_convex_quadrilateral( elmnt_coords, i_coord, j_coord, pixel_is_in)

! if inside, getting the nearest point inside the element
           if ( pixel_is_in ) then
              dist_min_pixel = HUGEVAL
              do k = 1, NGLLX
                 do l = 1, NGLLZ
                    iglob = ibool(k,l,ispec)
                    dist_pixel = (coord(1,iglob)-i_coord)**2 + (coord(2,iglob)-j_coord)**2
                    if (dist_pixel < dist_min_pixel) then
                       dist_min_pixel = dist_pixel
                       iglob_image_color(i,j) = iglob

                    endif

                 enddo
              enddo
              if ( dist_min_pixel >= HUGEVAL ) then
                 call exit_MPI('Error in detecting pixel for color image')

              endif
              nb_pixel_loc = nb_pixel_loc + 1

           endif

        enddo
     enddo
  enddo

! creating and filling array num_pixel_loc with the positions of each colored
! pixel owned by the local process (useful for parallel jobs)
  allocate(num_pixel_loc(nb_pixel_loc))

  nb_pixel_loc = 0
  do i = 1, NX_IMAGE_color
     do j = 1, NZ_IMAGE_color
        if ( iglob_image_color(i,j) /= -1 ) then
           nb_pixel_loc = nb_pixel_loc + 1
           num_pixel_loc(nb_pixel_loc) = (j-1)*NX_IMAGE_color + i

        endif

     enddo
  enddo

! filling array iglob_image_color, containing info on which process owns which pixels.
#ifdef USE_MPI
  allocate(nb_pixel_per_proc(nproc))

  call MPI_GATHER( nb_pixel_loc, 1, MPI_INTEGER, nb_pixel_per_proc(1), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

  if ( myrank == 0 ) then
     allocate(num_pixel_recv(maxval(nb_pixel_per_proc(:)),nproc))
     allocate(data_pixel_recv(maxval(nb_pixel_per_proc(:))))
  endif

  allocate(data_pixel_send(nb_pixel_loc))
  if (nproc > 1) then
     if (myrank == 0) then

        do iproc = 1, nproc-1

           call MPI_RECV(num_pixel_recv(1,iproc+1),nb_pixel_per_proc(iproc+1), MPI_INTEGER, &
                iproc, 42, MPI_COMM_WORLD, request_mpi_status, ier)
           do k = 1, nb_pixel_per_proc(iproc+1)
              j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
              i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color
              iglob_image_color(i,j) = iproc

           enddo
        enddo

     else
        call MPI_SEND(num_pixel_loc(1),nb_pixel_loc,MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)

     endif

  endif
#else
   allocate(nb_pixel_per_proc(1))
   deallocate(nb_pixel_per_proc)
   allocate(num_pixel_recv(1,1))
   deallocate(num_pixel_recv)
   allocate(data_pixel_recv(1))
   deallocate(data_pixel_recv)
   allocate(data_pixel_send(1))
   deallocate(data_pixel_send)
#endif

  if (myrank == 0) write(IOUT,*) 'done locating all the pixels of color images'

  endif

!
!---- initialize seismograms
!
  sisux = ZERO
  sisuz = ZERO

  cosrot = cos(anglerec)
  sinrot = sin(anglerec)

! initialize arrays to zero
  displ_elastic = ZERO
  veloc_elastic = ZERO
  accel_elastic = ZERO

  displs_poroelastic = ZERO
  velocs_poroelastic = ZERO
  accels_poroelastic = ZERO
  displw_poroelastic = ZERO
  velocw_poroelastic = ZERO
  accelw_poroelastic = ZERO

  potential_acoustic = ZERO
  potential_dot_acoustic = ZERO
  potential_dot_dot_acoustic = ZERO

!
!----- Files where absorbing signal are saved during forward wavefield calculation
!

  if( ((save_forward .and. isolver ==1) .or. isolver == 2) .and. anyabs ) then

     if(any_elastic) then

!--- left absorbing boundary
      if( nspec_xmin >0 ) then
        if(isolver == 2) then
          open(unit=35,file='OUTPUT_FILES/absorb_elastic_left.bin',status='old',&
                form='unformatted')
        else
          open(unit=35,file='OUTPUT_FILES/absorb_elastic_left.bin',status='unknown',&
                form='unformatted')
        endif

      endif  !  end of left absorbing boundary

!--- right absorbing boundary
      if( nspec_xmax >0 ) then
        if(isolver == 2) then
          open(unit=36,file='OUTPUT_FILES/absorb_elastic_right.bin',status='old',&
                form='unformatted')
        else
          open(unit=36,file='OUTPUT_FILES/absorb_elastic_right.bin',status='unknown',&
                form='unformatted')
        endif

      endif  !  end of right absorbing boundary

!--- bottom absorbing boundary
      if( nspec_zmin >0 ) then
        if(isolver == 2) then
          open(unit=37,file='OUTPUT_FILES/absorb_elastic_bottom.bin',status='old',&
                form='unformatted')
        else
          open(unit=37,file='OUTPUT_FILES/absorb_elastic_bottom.bin',status='unknown',&
                form='unformatted')
        endif

      endif  !  end of bottom absorbing boundary

!--- top absorbing boundary
      if( nspec_zmax >0 ) then
        if(isolver == 2) then
          open(unit=38,file='OUTPUT_FILES/absorb_elastic_top.bin',status='old',&
                form='unformatted')
        else
          open(unit=38,file='OUTPUT_FILES/absorb_elastic_top.bin',status='unknown',&
                form='unformatted')
        endif

      endif ! end of top absorbing boundary

     endif

     if(any_poroelastic) then

!--- left absorbing boundary
      if( nspec_xmin >0 ) then
        if(isolver == 2) then
          open(unit=45,file='OUTPUT_FILES/absorb_poro_s_left.bin',status='old',&
                form='unformatted')
          open(unit=25,file='OUTPUT_FILES/absorb_poro_w_left.bin',status='old',&
                form='unformatted')
        else
          open(unit=45,file='OUTPUT_FILES/absorb_poro_s_left.bin',status='unknown',&
                form='unformatted')
          open(unit=25,file='OUTPUT_FILES/absorb_poro_w_left.bin',status='unknown',&
                form='unformatted')
        endif

      endif  !  end of left absorbing boundary

!--- right absorbing boundary
      if( nspec_xmax >0 ) then
        if(isolver == 2) then
          open(unit=46,file='OUTPUT_FILES/absorb_poro_s_right.bin',status='old',&
                form='unformatted')
          open(unit=26,file='OUTPUT_FILES/absorb_poro_w_right.bin',status='old',&
                form='unformatted')
        else
          open(unit=46,file='OUTPUT_FILES/absorb_poro_s_right.bin',status='unknown',&
                form='unformatted')
          open(unit=26,file='OUTPUT_FILES/absorb_poro_w_right.bin',status='unknown',&
                form='unformatted')
        endif

      endif  !  end of right absorbing boundary

!--- bottom absorbing boundary
      if( nspec_zmin >0 ) then
        if(isolver == 2) then
          open(unit=47,file='OUTPUT_FILES/absorb_poro_s_bottom.bin',status='old',&
                form='unformatted')
          open(unit=29,file='OUTPUT_FILES/absorb_poro_w_bottom.bin',status='old',&
                form='unformatted')
        else
          open(unit=47,file='OUTPUT_FILES/absorb_poro_s_bottom.bin',status='unknown',&
                form='unformatted')
          open(unit=29,file='OUTPUT_FILES/absorb_poro_w_bottom.bin',status='unknown',&
                form='unformatted')
        endif

      endif  !  end of bottom absorbing boundary

!--- top absorbing boundary
      if( nspec_zmax >0 ) then
        if(isolver == 2) then
          open(unit=48,file='OUTPUT_FILES/absorb_poro_s_top.bin',status='old',&
                form='unformatted')
          open(unit=28,file='OUTPUT_FILES/absorb_poro_w_top.bin',status='old',&
                form='unformatted')
        else
          open(unit=48,file='OUTPUT_FILES/absorb_poro_s_top.bin',status='unknown',&
                form='unformatted')
          open(unit=28,file='OUTPUT_FILES/absorb_poro_w_top.bin',status='unknown',&
                form='unformatted')
        endif

      endif ! end of top absorbing boundary

     endif

     if(any_acoustic) then

!--- left absorbing boundary
      if( nspec_xmin >0 ) then
        print*,'Opening absorb_acoustic_left.bin....'
        if(isolver == 2) then
          open(unit=65,file='OUTPUT_FILES/absorb_acoustic_left.bin',status='old',&
                form='unformatted')
        else
          open(unit=65,file='OUTPUT_FILES/absorb_acoustic_left.bin',status='unknown',&
                form='unformatted')
        endif

      endif  !  end of left absorbing boundary

!--- right absorbing boundary
      if( nspec_xmax >0 ) then
        print*,'Opening absorb_acoustic_right.bin....'
        if(isolver == 2) then
          open(unit=66,file='OUTPUT_FILES/absorb_acoustic_right.bin',status='old',&
                form='unformatted')
        else
          open(unit=66,file='OUTPUT_FILES/absorb_acoustic_right.bin',status='unknown',&
                form='unformatted')
        endif

      endif  !  end of right absorbing boundary

!--- bottom absorbing boundary
      if( nspec_zmin >0 ) then
        print*,'Opening absorb_acoustic_bottom.bin....'
        if(isolver == 2) then
          open(unit=67,file='OUTPUT_FILES/absorb_acoustic_bottom.bin',status='old',&
                form='unformatted')
        else
          open(unit=67,file='OUTPUT_FILES/absorb_acoustic_bottom.bin',status='unknown',&
                form='unformatted')
        endif

      endif  !  end of bottom absorbing boundary

!--- top absorbing boundary
      if( nspec_zmax >0 ) then
        print*,'Opening absorb_acoustic_top.bin....'
        if(isolver == 2) then
          open(unit=68,file='OUTPUT_FILES/absorb_acoustic_top.bin',status='old',&
                form='unformatted')
        else
          open(unit=68,file='OUTPUT_FILES/absorb_acoustic_top.bin',status='unknown',&
                form='unformatted')
        endif

      endif ! end of top absorbing boundary

     endif

    endif !if( ((save_forward .and. isolver ==1) .or. isolver == 2) .and. anyabs )


    if(anyabs .and. isolver == 2) then

      if(any_elastic) then

     do it =1, NSTEP

!--- left absorbing boundary
      if(nspec_xmin >0) then
      do ispec = 1,nspec_xmin
       do id =1,2
         do i=1,NGLLZ
     read(35) b_absorb_elastic_left(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

!--- right absorbing boundary
      if(nspec_xmax >0) then
      do ispec = 1,nspec_xmax
       do id =1,2
         do i=1,NGLLZ
     read(36) b_absorb_elastic_right(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

!--- bottom absorbing boundary
      if(nspec_zmin >0) then
      do ispec = 1,nspec_zmin
       do id =1,2
         do i=1,NGLLX
     read(37) b_absorb_elastic_bottom(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

!--- top absorbing boundary
      if(nspec_zmax >0) then
      do ispec = 1,nspec_zmax
       do id =1,2
         do i=1,NGLLX
     read(38) b_absorb_elastic_top(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif


   enddo

      endif ! if(any_elastic)

      if(any_poroelastic) then

     do it =1, NSTEP

!--- left absorbing boundary
      if(nspec_xmin >0) then
      do ispec = 1,nspec_xmin
       do id =1,2
         do i=1,NGLLZ
     read(45) b_absorb_poro_s_left(id,i,ispec,it)
     read(25) b_absorb_poro_w_left(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

!--- right absorbing boundary
      if(nspec_xmax >0) then
      do ispec = 1,nspec_xmax
       do id =1,2
         do i=1,NGLLZ
     read(46) b_absorb_poro_s_right(id,i,ispec,it)
     read(26) b_absorb_poro_w_right(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

!--- bottom absorbing boundary
      if(nspec_zmin >0) then
      do ispec = 1,nspec_zmin
       do id =1,2
         do i=1,NGLLX
     read(47) b_absorb_poro_s_bottom(id,i,ispec,it)
     read(29) b_absorb_poro_w_bottom(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

!--- top absorbing boundary
      if(nspec_zmax >0) then
      do ispec = 1,nspec_zmax
       do id =1,2
         do i=1,NGLLX
     read(48) b_absorb_poro_s_top(id,i,ispec,it)
     read(28) b_absorb_poro_w_top(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif


   enddo

      endif ! if(any_poroelastic)

      if(any_acoustic) then

     do it =1, NSTEP

!--- left absorbing boundary
      if(nspec_xmin >0) then
      do ispec = 1,nspec_xmin
         do i=1,NGLLZ
     read(65) b_absorb_acoustic_left(i,ispec,it)
         enddo
      enddo
      endif

!--- right absorbing boundary
      if(nspec_xmax >0) then
      do ispec = 1,nspec_xmax
         do i=1,NGLLZ
     read(66) b_absorb_acoustic_right(i,ispec,it)
         enddo
      enddo
      endif

!--- bottom absorbing boundary
      if(nspec_zmin >0) then
      do ispec = 1,nspec_zmin
         do i=1,NGLLX
     read(67) b_absorb_acoustic_bottom(i,ispec,it)
         enddo
      enddo
      endif

!--- top absorbing boundary
      if(nspec_zmax >0) then
      do ispec = 1,nspec_zmax
         do i=1,NGLLX
     read(68) b_absorb_acoustic_top(i,ispec,it)
         enddo
      enddo
      endif


   enddo

      endif ! if(any_acoustic)


    endif ! if(anyabs .and. isolver == 2)



!
!----- Read last frame for backward wavefield calculation
!

  if(isolver == 2) then

   if(any_elastic) then
    open(unit=55,file='OUTPUT_FILES/lastframe_elastic.bin',status='old',action='read',form='unformatted')
       do j=1,npoin
      read(55) (b_displ_elastic(i,j), i=1,NDIM), &
                  (b_veloc_elastic(i,j), i=1,NDIM), &
                  (b_accel_elastic(i,j), i=1,NDIM)
       enddo
    close(55)

  rho_kl(:) = ZERO
  mu_kl(:) = ZERO
  kappa_kl(:) = ZERO
!
  rhop_kl(:) = ZERO
  beta_kl(:) = ZERO
  alpha_kl(:) = ZERO
   endif

   if(any_poroelastic) then
    open(unit=55,file='OUTPUT_FILES/lastframe_poroelastic_s.bin',status='old',action='read',form='unformatted')
    open(unit=56,file='OUTPUT_FILES/lastframe_poroelastic_w.bin',status='old',action='read',form='unformatted')
       do j=1,npoin
      read(55) (b_displs_poroelastic(i,j), i=1,NDIM), &
                  (b_velocs_poroelastic(i,j), i=1,NDIM), &
                  (b_accels_poroelastic(i,j), i=1,NDIM)
      read(56) (b_displw_poroelastic(i,j), i=1,NDIM), &
                  (b_velocw_poroelastic(i,j), i=1,NDIM), &
                  (b_accelw_poroelastic(i,j), i=1,NDIM)
       enddo
    close(55)
    close(56)

  rhot_kl(:) = ZERO
  rhof_kl(:) = ZERO
  eta_kl(:) = ZERO
  sm_kl(:) = ZERO
  mufr_kl(:) = ZERO
  B_kl(:) = ZERO
  C_kl(:) = ZERO
  M_kl(:) = ZERO
!
  rhob_kl(:) = ZERO
  rhofb_kl(:) = ZERO
  phi_kl(:) = ZERO
  mufrb_kl(:) = ZERO
  Bb_kl(:) = ZERO
  Cb_kl(:) = ZERO
  Mb_kl(:) = ZERO
!
  rhobb_kl(:) = ZERO
  rhofbb_kl(:) = ZERO
  phib_kl(:) = ZERO
  cs_kl(:) = ZERO
  cpI_kl(:) = ZERO
  cpII_kl(:) = ZERO
  ratio_kl(:) = ZERO
   endif

   if(any_acoustic) then
    open(unit=55,file='OUTPUT_FILES/lastframe_acoustic.bin',status='old',action='read',form='unformatted')
       do j=1,npoin
      read(55) b_potential_acoustic(j),&
               b_potential_dot_acoustic(j),&
               b_potential_dot_dot_acoustic(j)
       enddo
    close(55)

  rho_ac_kl(:) = ZERO
  kappa_ac_kl(:) = ZERO
!
  rhop_ac_kl(:) = ZERO
  alpha_ac_kl(:) = ZERO
   endif

  endif ! if(isover == 2)

!
!----  read initial fields from external file if needed
!

! if we are looking a plane wave beyond critical angle we use other method
  over_critical_angle = .false.

  if(initialfield) then
      if (myrank == 0) then
         write(IOUT,*)
!! DK DK reading of an initial field from an external file has been suppressed
!! DK DK and replaced with the implementation of an analytical plane wave
!! DK DK     write(IOUT,*) 'Reading initial fields from external file...'
         write(IOUT,*) 'Implementing an analytical initial plane wave...'
         write(IOUT,*)
      endif
      if(any_acoustic .or. any_poroelastic) call exit_MPI('initial field currently implemented for purely elastic simulation only')

      !=======================================================================
      !
      !     Calculation of the initial field for a plane wave
      !
      !=======================================================================

      if (myrank == 0) then
         write(IOUT,*) 'Number of grid points: ',npoin
         write(IOUT,*)
         write(IOUT,*) '*** calculation of the initial plane wave ***'
         write(IOUT,*)
         write(IOUT,*)  'To change the initial plane wave, change source_type in DATA/Par_file'
         write(IOUT,*)  'and use 1 for a plane P wave, 2 for a plane SV wave, 3 for a Rayleigh wave'
         write(IOUT,*)

! only implemented for one source
         if(NSOURCE > 1) call exit_MPI('calculation of the initial wave is only implemented for one source')
         if (source_type(1) == 1) then
            write(IOUT,*) 'initial P wave of', angleforce(1)*180.d0/pi, 'degrees introduced.'
         else if (source_type(1) == 2) then
            write(IOUT,*) 'initial SV wave of', angleforce(1)*180.d0/pi, ' degrees introduced.'

         else if (source_type(1) == 3) then
            write(IOUT,*) 'Rayleigh wave introduced.'
         else
            call exit_MPI('Unrecognized source_type: should be 1 for plane P waves, 2 for plane SV waves, 3 for Rayleigh wave')
         endif

         if ((angleforce(1) < 0.0d0 .or. angleforce(1) >= pi/2.d0) .and. source_type(1) /= 3) then
            call exit_MPI("incorrect angleforce: must have 0 <= angleforce < 90")
         endif
      endif
      ! only implemented for homogeneous media therefore only 1 material supported
      if (numat==1) then

         mu = poroelastcoef(2,1,numat)
         lambdaplus2mu  = poroelastcoef(3,1,numat)
         denst = density(1,numat)

         cploc = sqrt(lambdaplus2mu/denst)
         csloc = sqrt(mu/denst)
      
         ! P wave case
         if (source_type(1) == 1) then

            p=sin(angleforce(1))/cploc
            c_inc  = cploc
            c_refl = csloc

            angleforce_refl = asin(p*c_refl)

            ! from formulas (5.26) and (5.27) p 140 in Aki & Richards (1980)
            PP = (- cos(2.d0*angleforce_refl)**2/csloc**3 + 4.d0*p**2*cos(angleforce(1))*cos(angleforce_refl)/cploc) / &
                 (  cos(2.d0*angleforce_refl)**2/csloc**3 + 4.d0*p**2*cos(angleforce(1))*cos(angleforce_refl)/cploc)

            PS = 4.d0*p*cos(angleforce(1))*cos(2.d0*angleforce_refl) / &
                 (csloc**2*(cos(2.d0*angleforce_refl)**2/csloc**3 &
                 +4.d0*p**2*cos(angleforce(1))*cos(angleforce_refl)/cploc))

             if (myrank == 0) then
                write(IOUT,*) 'reflected convert plane wave angle: ', angleforce_refl*180.d0/pi
             endif

            ! from Table 5.1 p141 in Aki & Richards (1980)
            ! we put the opposite sign on z coefficients because z axis is oriented from bottom to top
            A_plane(1) = sin(angleforce(1));           A_plane(2) = cos(angleforce(1))
            B_plane(1) = PP * sin(angleforce(1));      B_plane(2) = - PP * cos(angleforce(1))
            C_plane(1) = PS * cos(angleforce_refl);    C_plane(2) = PS * sin(angleforce_refl)

         ! SV wave case
         else if (source_type(1) == 2) then

            p=sin(angleforce(1))/csloc
            c_inc  = csloc
            c_refl = cploc

            ! if this coefficient is greater than 1, we are beyond the critical SV wave angle and there cannot be a converted P wave
            if (p*c_refl<=1.d0) then
               angleforce_refl = asin(p*c_refl)

               ! from formulas (5.30) and (5.31) p 140 in Aki & Richards (1980)
               SS = (cos(2.d0*angleforce(1))**2/csloc**3 - 4.d0*p**2*cos(angleforce(1))*cos(angleforce_refl)/cploc) / &
                    (cos(2.d0*angleforce(1))**2/csloc**3 + 4.d0*p**2*cos(angleforce(1))*cos(angleforce_refl)/cploc)
               SP = 4.d0*p*cos(angleforce(1))*cos(2*angleforce(1)) / &
                    (cploc*csloc*(cos(2.d0*angleforce(1))**2/csloc**3&
                    +4.d0*p**2*cos(angleforce_refl)*cos(angleforce(1))/cploc))

               if (myrank == 0) then
                  write(IOUT,*) 'reflected convert plane wave angle: ', angleforce_refl*180.d0/pi
               endif

            ! SV45 degree incident plane wave is a particular case
            else if (angleforce(1)>pi/4.d0-1.0d-11 .and. angleforce(1)<pi/4.d0+1.0d-11) then
               angleforce_refl = 0.d0
               SS = -1.0d0
               SP = 0.d0
            else
               over_critical_angle=.true.
               angleforce_refl = 0.d0
               SS = 0.0d0
               SP = 0.d0
            endif

            ! from Table 5.1 p141 in Aki & Richards (1980)
            ! we put the opposite sign on z coefficients because z axis is oriented from bottom to top
            A_plane(1) = cos(angleforce(1));           A_plane(2) = - sin(angleforce(1))
            B_plane(1) = SS * cos(angleforce(1));      B_plane(2) = SS * sin(angleforce(1))
            C_plane(1) = SP * sin(angleforce_refl);    C_plane(2) = - SP * cos(angleforce_refl)

         ! Rayleigh case
         else if (source_type(1) == 3) then
            over_critical_angle=.true.
            A_plane(1)=0.d0; A_plane(2)=0.d0
            B_plane(1)=0.d0; B_plane(2)=0.d0
            C_plane(1)=0.d0; C_plane(2)=0.d0
         endif
      else
         call exit_MPI('not possible to have several materials with a plane wave')
      endif

      ! get minimum and maximum values of mesh coordinates
      xmin = minval(coord(1,:))
      zmin = minval(coord(2,:))
      xmax = maxval(coord(1,:))
      zmax = maxval(coord(2,:))

#ifdef USE_MPI
      call MPI_ALLREDUCE (xmin, xmin_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
      call MPI_ALLREDUCE (zmin, zmin_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
      call MPI_ALLREDUCE (xmax, xmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
      call MPI_ALLREDUCE (zmax, zmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
      xmin = xmin_glob
      zmin = zmin_glob
      xmax = xmax_glob
      zmax = zmax_glob
#endif

      ! initialize the time offset to put the plane wave not too close to the irregularity on the free surface
      if (abs(angleforce(1))<1.d0*pi/180.d0 .and. source_type(1)/=3) then
         time_offset=-1.d0*(zmax-zmin)/2.d0/c_inc
      else
         time_offset=0.d0
      endif

      ! to correctly center the initial plane wave in the mesh
      x0_source=x_source(1)
      z0_source=z_source(1)

      if (myrank == 0) then
         write(IOUT,*)
         write(IOUT,*) 'You can modify the location of the initial plane wave by changing xs and zs in DATA/Par_File.'
         write(IOUT,*) '   for instance: xs=',x_source(1),'   zs=',z_source(1), ' (zs must be the height of the free surface)'
         write(IOUT,*)
      endif

      if (.not. over_critical_angle) then

         do i = 1,npoin

            x = coord(1,i)
            z = coord(2,i)

            ! z is from bottom to top therefore we take -z to make parallel with Aki & Richards
            z = z0_source - z
            x = x - x0_source

            t = 0.d0 + time_offset

            ! formulas for the initial displacement for a plane wave from Aki & Richards (1980)
       displ_elastic(1,i) = A_plane(1) * ricker_Bielak_displ(t - sin(angleforce(1))*x/c_inc + cos(angleforce(1))*z/c_inc,f0(1)) &
                 + B_plane(1) * ricker_Bielak_displ(t - sin(angleforce(1))*x/c_inc - cos(angleforce(1))*z/c_inc,f0(1)) &
                 + C_plane(1) * ricker_Bielak_displ(t - sin(angleforce_refl)*x/c_refl - cos(angleforce_refl)*z/c_refl,f0(1))
       displ_elastic(2,i) = A_plane(2) * ricker_Bielak_displ(t - sin(angleforce(1))*x/c_inc + cos(angleforce(1))*z/c_inc,f0(1)) &
                 + B_plane(2) * ricker_Bielak_displ(t - sin(angleforce(1))*x/c_inc - cos(angleforce(1))*z/c_inc,f0(1)) &
                 + C_plane(2) * ricker_Bielak_displ(t - sin(angleforce_refl)*x/c_refl - cos(angleforce_refl)*z/c_refl,f0(1))

            ! formulas for the initial velocity for a plane wave (first derivative in time of the displacement)
       veloc_elastic(1,i) = A_plane(1) * ricker_Bielak_veloc(t - sin(angleforce(1))*x/c_inc + cos(angleforce(1))*z/c_inc,f0(1)) &
                 + B_plane(1) * ricker_Bielak_veloc(t - sin(angleforce(1))*x/c_inc - cos(angleforce(1))*z/c_inc,f0(1)) &
                 + C_plane(1) * ricker_Bielak_veloc(t - sin(angleforce_refl)*x/c_refl - cos(angleforce_refl)*z/c_refl,f0(1))
       veloc_elastic(2,i) = A_plane(2) * ricker_Bielak_veloc(t - sin(angleforce(1))*x/c_inc + cos(angleforce(1))*z/c_inc,f0(1)) &
                 + B_plane(2) * ricker_Bielak_veloc(t - sin(angleforce(1))*x/c_inc - cos(angleforce(1))*z/c_inc,f0(1)) &
                 + C_plane(2) * ricker_Bielak_veloc(t - sin(angleforce_refl)*x/c_refl - cos(angleforce_refl)*z/c_refl,f0(1))

            ! formulas for the initial acceleration for a plane wave (second derivative in time of the displacement)
       accel_elastic(1,i) = A_plane(1) * ricker_Bielak_accel(t - sin(angleforce(1))*x/c_inc + cos(angleforce(1))*z/c_inc,f0(1)) &
                 + B_plane(1) * ricker_Bielak_accel(t - sin(angleforce(1))*x/c_inc - cos(angleforce(1))*z/c_inc,f0(1)) &
                 + C_plane(1) * ricker_Bielak_accel(t - sin(angleforce_refl)*x/c_refl - cos(angleforce_refl)*z/c_refl,f0(1))
       accel_elastic(2,i) = A_plane(2) * ricker_Bielak_accel(t - sin(angleforce(1))*x/c_inc + cos(angleforce(1))*z/c_inc,f0(1)) &
                 + B_plane(2) * ricker_Bielak_accel(t - sin(angleforce(1))*x/c_inc - cos(angleforce(1))*z/c_inc,f0(1)) &
                 + C_plane(2) * ricker_Bielak_accel(t - sin(angleforce_refl)*x/c_refl - cos(angleforce_refl)*z/c_refl,f0(1))

         enddo

      else ! beyond critical angle

         if (myrank == 0) then
            if (source_type(1)/=3) write(IOUT,*) 'You are beyond the critical angle ( > ',asin(c_inc/c_refl)*180d0/pi,')'

            write(IOUT,*)  '*************'
            write(IOUT,*)  'We have to compute the initial field in the frequency domain'
            write(IOUT,*)  'and then convert it to the time domain (can be long... be patient...)'
            write(IOUT,*)  '*************'
         endif

         allocate(left_bound(nelemabs*NGLLX))
         allocate(right_bound(nelemabs*NGLLX))
         allocate(bot_bound(nelemabs*NGLLZ))

         count_bot=0
         count_left=0
         count_right=0
         do ispecabs=1,nelemabs
            ispec=numabs(ispecabs)
            if(codeabs(ILEFT,ispecabs)) then
               i = 1
               do j = 1,NGLLZ
                  count_left=count_left+1
                  iglob = ibool(i,j,ispec)
                  left_bound(count_left)=iglob
               enddo
            endif
            if(codeabs(IRIGHT,ispecabs)) then
               i = NGLLX
               do j = 1,NGLLZ
                  count_right=count_right+1
                  iglob = ibool(i,j,ispec)
                  right_bound(count_right)=iglob
               enddo
            endif
            if(codeabs(IBOTTOM,ispecabs)) then
               j = 1
               ! exclude corners to make sure there is no contradiction regarding the normal
               ibegin = 1
               iend = NGLLX
               if(codeabs(ILEFT,ispecabs)) ibegin = 2
               if(codeabs(IRIGHT,ispecabs)) iend = NGLLX-1
               do i = ibegin,iend
                  count_bot=count_bot+1
                  iglob = ibool(i,j,ispec)
                  bot_bound(count_bot)=iglob
               enddo
            endif
         enddo

         allocate(v0x_left(count_left,NSTEP))
         allocate(v0z_left(count_left,NSTEP))
         allocate(t0x_left(count_left,NSTEP))
         allocate(t0z_left(count_left,NSTEP))

         allocate(v0x_right(count_right,NSTEP))
         allocate(v0z_right(count_right,NSTEP))
         allocate(t0x_right(count_right,NSTEP))
         allocate(t0z_right(count_right,NSTEP))

         allocate(v0x_bot(count_bot,NSTEP))
         allocate(v0z_bot(count_bot,NSTEP))
         allocate(t0x_bot(count_bot,NSTEP))
         allocate(t0z_bot(count_bot,NSTEP))

! call Paco's routine to compute in frequency and convert to time by Fourier transform
         call paco_beyond_critical(coord,npoin,deltat,NSTEP,angleforce(1),&
              f0(1),cploc,csloc,TURN_ATTENUATION_ON,Qp_attenuation,source_type(1),v0x_left,v0z_left,&
              v0x_right,v0z_right,v0x_bot,v0z_bot,t0x_left,t0z_left,t0x_right,t0z_right,&
              t0x_bot,t0z_bot,left_bound(1:count_left),right_bound(1:count_right),bot_bound(1:count_bot)&
              ,count_left,count_right,count_bot,displ_elastic,veloc_elastic,accel_elastic)

         deallocate(left_bound)
         deallocate(right_bound)
         deallocate(bot_bound)

         if (myrank == 0) then
            write(IOUT,*)  '***********'
            write(IOUT,*)  'done calculating the initial wave field'
            write(IOUT,*)  '***********'
         endif

      endif ! beyond critical angle

    write(IOUT,*) 'Max norm of initial elastic displacement = ',maxval(sqrt(displ_elastic(1,:)**2 + displ_elastic(2,:)**2))
  endif ! initialfield

  deltatsquare = deltat * deltat
  deltatcube = deltatsquare * deltat
  deltatfourth = deltatsquare * deltatsquare

  twelvedeltat = 12.d0 * deltat
  fourdeltatsquare = 4.d0 * deltatsquare

! compute the source time function and store it in a text file
  if(.not. initialfield) then

    allocate(source_time_function(NSOURCE,NSTEP))

    if (myrank == 0) then
      write(IOUT,*)
      write(IOUT,*) 'Saving the source time function in a text file...'
      write(IOUT,*)
      open(unit=55,file='OUTPUT_FILES/source.txt',status='unknown')
    endif

! loop on all the sources
    do i_source=1,NSOURCE

! loop on all the time steps
    do it = 1,NSTEP

! compute current time
      time = (it-1)*deltat

! Ricker (second derivative of a Gaussian) source time function
      if(time_function_type(i_source) == 1) then
!        source_time_function(i_source,it) = - factor(i_source) * (ONE-TWO*aval(i_source)*(time-t0(i_source))**2) * &
!                                           exp(-aval(i_source)*(time-t0(i_source))**2)
        source_time_function(i_source,it) = - factor(i_source) * TWO*aval(i_source)*sqrt(aval(i_source))*&
                                            (time-t0(i_source))/pi * exp(-aval(i_source)*(time-t0(i_source))**2)

! first derivative of a Gaussian source time function
      else if(time_function_type(i_source) == 2) then
        source_time_function(i_source,it) = - factor(i_source) * TWO*aval(i_source)*(time-t0(i_source)) * &
                                           exp(-aval(i_source)*(time-t0(i_source))**2)

! Gaussian or Dirac (we use a very thin Gaussian instead) source time function
      else if(time_function_type(i_source) == 3 .or. time_function_type(i_source) == 4) then
        source_time_function(i_source,it) = factor(i_source) * exp(-aval(i_source)*(time-t0(i_source))**2)

! Heaviside source time function (we use a very thin error function instead)
      else if(time_function_type(i_source) == 5) then
        hdur(i_source) = 1.d0 / f0(i_source)
        hdur_gauss(i_source) = hdur(i_source) * 5.d0 / 3.d0
        source_time_function(i_source,it) = factor(i_source) * 0.5d0*(1.0d0 + &
                                           netlib_specfun_erf(SOURCE_DECAY_MIMIC_TRIANGLE*(time-t0(i_source))/hdur_gauss(i_source)))

      else
        call exit_MPI('unknown source time function')
      endif

! output absolute time in third column, in case user wants to check it as well
      if (myrank == 0 .and. i_source==1 ) write(55,*) sngl(time),real(source_time_function(1,it),4),sngl(time-t0(1))
   enddo
   enddo ! i_source=1,NSOURCE 

   if (myrank == 0) close(55)

! nb_proc_source is the number of processes that own the source (the nearest point). It can be greater
! than one if the nearest point is on the interface between several partitions with an explosive source.
! since source contribution is linear, the source_time_function is cut down by that number (it would have been similar
! if we just had elected one of those processes).
   do i_source=1,NSOURCE
    source_time_function(i_source,:) = source_time_function(i_source,:) / nb_proc_source(i_source)
   enddo

  else

    allocate(source_time_function(1,1))

  endif

! determine if coupled fluid-solid simulation
  coupled_acoustic_elastic = any_acoustic .and. any_elastic
  coupled_acoustic_poroelastic = any_acoustic .and. any_poroelastic

! fluid/solid (elastic) edge detection
! the two elements (fluid and solid) forming an edge are already known (computed in meshfem2D),
! the common nodes forming the edge are computed here
  if(coupled_acoustic_elastic) then

    if (myrank == 0) then
      print *
      print *,'Mixed acoustic/elastic simulation'
      print *
      print *,'Beginning of fluid/solid edge detection'
    endif

! define the edges of a given element
    i_begin(IBOTTOM) = 1
    j_begin(IBOTTOM) = 1
    i_end(IBOTTOM) = NGLLX
    j_end(IBOTTOM) = 1

    i_begin(IRIGHT) = NGLLX
    j_begin(IRIGHT) = 1
    i_end(IRIGHT) = NGLLX
    j_end(IRIGHT) = NGLLZ

    i_begin(ITOP) = NGLLX
    j_begin(ITOP) = NGLLZ
    i_end(ITOP) = 1
    j_end(ITOP) = NGLLZ

    i_begin(ILEFT) = 1
    j_begin(ILEFT) = NGLLZ
    i_end(ILEFT) = 1
    j_end(ILEFT) = 1

! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = ipoin1D
      ivalue_inverse(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,IBOTTOM) = 1
      jvalue_inverse(ipoin1D,IBOTTOM) = 1

      ivalue(ipoin1D,IRIGHT) = NGLLX
      ivalue_inverse(ipoin1D,IRIGHT) = NGLLX
      jvalue(ipoin1D,IRIGHT) = ipoin1D
      jvalue_inverse(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1

      ivalue(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,ITOP) = ipoin1D
      jvalue(ipoin1D,ITOP) = NGLLZ
      jvalue_inverse(ipoin1D,ITOP) = NGLLZ

      ivalue(ipoin1D,ILEFT) = 1
      ivalue_inverse(ipoin1D,ILEFT) = 1
      jvalue(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,ILEFT) = ipoin1D

    enddo

    do inum = 1, num_fluid_solid_edges
       ispec_acoustic =  fluid_solid_acoustic_ispec(inum)
       ispec_elastic =  fluid_solid_elastic_ispec(inum)

! one element must be acoustic and the other must be elastic
        if(ispec_acoustic /= ispec_elastic .and. .not. elastic(ispec_acoustic) .and. &
             .not. poroelastic(ispec_acoustic) .and. elastic(ispec_elastic)) then

! loop on the four edges of the two elements
          do iedge_acoustic = 1,NEDGES
            do iedge_elastic = 1,NEDGES

! store the matching topology if the two edges match in inverse order
              if(ibool(i_begin(iedge_acoustic),j_begin(iedge_acoustic),ispec_acoustic) == &
                   ibool(i_end(iedge_elastic),j_end(iedge_elastic),ispec_elastic) .and. &
                   ibool(i_end(iedge_acoustic),j_end(iedge_acoustic),ispec_acoustic) == &
                   ibool(i_begin(iedge_elastic),j_begin(iedge_elastic),ispec_elastic)) then
                 fluid_solid_acoustic_iedge(inum) = iedge_acoustic
                 fluid_solid_elastic_iedge(inum) = iedge_elastic
!                  print *,'edge ',iedge_acoustic,' of acoustic element ',ispec_acoustic, &
!                          ' is in contact with edge ',iedge_elastic,' of elastic element ',ispec_elastic
                endif

             enddo
          enddo

       endif

    enddo

! make sure fluid/solid (elastic) matching has been perfectly detected: check that the grid points
! have the same physical coordinates
! loop on all the coupling edges

    if(myrank == 0) print *,'Checking fluid/solid edge topology...'

    do inum = 1,num_fluid_solid_edges

! get the edge of the acoustic element
      ispec_acoustic = fluid_solid_acoustic_ispec(inum)
      iedge_acoustic = fluid_solid_acoustic_iedge(inum)

! get the corresponding edge of the elastic element
      ispec_elastic = fluid_solid_elastic_ispec(inum)
      iedge_elastic = fluid_solid_elastic_iedge(inum)

! implement 1D coupling along the edge
      do ipoin1D = 1,NGLLX

! get point values for the elastic side, which matches our side in the inverse direction
        i = ivalue_inverse(ipoin1D,iedge_elastic)
        j = jvalue_inverse(ipoin1D,iedge_elastic)
        iglob = ibool(i,j,ispec_elastic)

! get point values for the acoustic side
        i = ivalue(ipoin1D,iedge_acoustic)
        j = jvalue(ipoin1D,iedge_acoustic)
        iglob2 = ibool(i,j,ispec_acoustic)

! if distance between the two points is not negligible, there is an error, since it should be zero
        if(sqrt((coord(1,iglob) - coord(1,iglob2))**2 + (coord(2,iglob) - coord(2,iglob2))**2) > TINYVAL) &
            call exit_MPI( 'error in fluid/solid coupling buffer')

      enddo

    enddo

    if (myrank == 0) then
      print *,'End of fluid/solid edge detection'
      print *
    endif

  endif

! fluid/solid (poroelastic) edge detection
! the two elements (fluid and solid) forming an edge are already known (computed in meshfem2D),
! the common nodes forming the edge are computed here
  if(coupled_acoustic_poroelastic) then
    if ( myrank == 0 ) then
    print *
    print *,'Mixed acoustic/poroelastic simulation'
    print *
    print *,'Beginning of fluid/solid (poroelastic) edge detection'
    endif

! define the edges of a given element
    i_begin(IBOTTOM) = 1
    j_begin(IBOTTOM) = 1
    i_end(IBOTTOM) = NGLLX
    j_end(IBOTTOM) = 1

    i_begin(IRIGHT) = NGLLX
    j_begin(IRIGHT) = 1
    i_end(IRIGHT) = NGLLX
    j_end(IRIGHT) = NGLLZ

    i_begin(ITOP) = NGLLX
    j_begin(ITOP) = NGLLZ
    i_end(ITOP) = 1
    j_end(ITOP) = NGLLZ

    i_begin(ILEFT) = 1
    j_begin(ILEFT) = NGLLZ
    i_end(ILEFT) = 1
    j_end(ILEFT) = 1

! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = ipoin1D
      ivalue_inverse(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,IBOTTOM) = 1
      jvalue_inverse(ipoin1D,IBOTTOM) = 1

      ivalue(ipoin1D,IRIGHT) = NGLLX
      ivalue_inverse(ipoin1D,IRIGHT) = NGLLX
      jvalue(ipoin1D,IRIGHT) = ipoin1D
      jvalue_inverse(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1

      ivalue(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,ITOP) = ipoin1D
      jvalue(ipoin1D,ITOP) = NGLLZ
      jvalue_inverse(ipoin1D,ITOP) = NGLLZ

      ivalue(ipoin1D,ILEFT) = 1
      ivalue_inverse(ipoin1D,ILEFT) = 1
      jvalue(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,ILEFT) = ipoin1D

    enddo

    do inum = 1, num_fluid_poro_edges
       ispec_acoustic =  fluid_poro_acoustic_ispec(inum)
       ispec_poroelastic =  fluid_poro_poroelastic_ispec(inum)

! one element must be acoustic and the other must be poroelastic
        if(ispec_acoustic /= ispec_poroelastic .and. .not. poroelastic(ispec_acoustic) .and. &
                 .not. elastic(ispec_acoustic) .and. poroelastic(ispec_poroelastic)) then

! loop on the four edges of the two elements
          do iedge_acoustic = 1,NEDGES
            do iedge_poroelastic = 1,NEDGES

! store the matching topology if the two edges match in inverse order
              if(ibool(i_begin(iedge_acoustic),j_begin(iedge_acoustic),ispec_acoustic) == &
                   ibool(i_end(iedge_poroelastic),j_end(iedge_poroelastic),ispec_poroelastic) .and. &
                   ibool(i_end(iedge_acoustic),j_end(iedge_acoustic),ispec_acoustic) == &
                   ibool(i_begin(iedge_poroelastic),j_begin(iedge_poroelastic),ispec_poroelastic)) then
                 fluid_poro_acoustic_iedge(inum) = iedge_acoustic
                 fluid_poro_poroelastic_iedge(inum) = iedge_poroelastic
                endif

             enddo
          enddo

       endif

    enddo


! make sure fluid/solid (poroelastic) matching has been perfectly detected: check that the grid points
! have the same physical coordinates
! loop on all the coupling edges

    if ( myrank == 0 ) then
    print *,'Checking fluid/solid (poroelastic) edge topology...'
    endif

    do inum = 1,num_fluid_poro_edges

! get the edge of the acoustic element
      ispec_acoustic = fluid_poro_acoustic_ispec(inum)
      iedge_acoustic = fluid_poro_acoustic_iedge(inum)

! get the corresponding edge of the poroelastic element
      ispec_poroelastic = fluid_poro_poroelastic_ispec(inum)
      iedge_poroelastic = fluid_poro_poroelastic_iedge(inum)

! implement 1D coupling along the edge
      do ipoin1D = 1,NGLLX

! get point values for the poroelastic side, which matches our side in the inverse direction
        i = ivalue_inverse(ipoin1D,iedge_poroelastic)
        j = jvalue_inverse(ipoin1D,iedge_poroelastic)
        iglob = ibool(i,j,ispec_poroelastic)

! get point values for the acoustic side
        i = ivalue(ipoin1D,iedge_acoustic)
        j = jvalue(ipoin1D,iedge_acoustic)
        iglob2 = ibool(i,j,ispec_acoustic)

! if distance between the two points is not negligible, there is an error, since it should be zero
        if(sqrt((coord(1,iglob) - coord(1,iglob2))**2 + (coord(2,iglob) - coord(2,iglob2))**2) > TINYVAL) &
            call exit_MPI( 'error in fluid/solid (poroelastic) coupling buffer')

      enddo

    enddo

    if ( myrank == 0 ) then
    print *,'End of fluid/solid (poroelastic) edge detection'
    print *
    endif

  endif

! exclude common points between acoustic absorbing edges and acoustic/elastic matching interfaces
  if(coupled_acoustic_elastic .and. anyabs) then

    if (myrank == 0) &
      print *,'excluding common points between acoustic absorbing edges and acoustic/elastic matching interfaces, if any'

! loop on all the absorbing elements
    do ispecabs = 1,nelemabs

      ispec = numabs(ispecabs)

! loop on all the coupling edges
      do inum = 1,num_fluid_solid_edges

! get the edge of the acoustic element
        ispec_acoustic = fluid_solid_acoustic_ispec(inum)
        iedge_acoustic = fluid_solid_acoustic_iedge(inum)

! if acoustic absorbing element and acoustic/elastic coupled element is the same
        if(ispec_acoustic == ispec) then

          if(iedge_acoustic == IBOTTOM) then
            jbegin_left(ispecabs) = 2
            jbegin_right(ispecabs) = 2
          endif

          if(iedge_acoustic == ITOP) then
            jend_left(ispecabs) = NGLLZ - 1
            jend_right(ispecabs) = NGLLZ - 1
          endif

          if(iedge_acoustic == ILEFT) then
            ibegin_bottom(ispecabs) = 2
            ibegin_top(ispecabs) = 2
          endif

          if(iedge_acoustic == IRIGHT) then
            iend_bottom(ispecabs) = NGLLX - 1
            iend_top(ispecabs) = NGLLX - 1
          endif

        endif

      enddo

    enddo

  endif

! exclude common points between acoustic absorbing edges and acoustic/poroelastic matching interfaces
  if(coupled_acoustic_poroelastic .and. anyabs) then

    if (myrank == 0) &
      print *,'excluding common points between acoustic absorbing edges and acoustic/poroelastic matching interfaces, if any'

! loop on all the absorbing elements
    do ispecabs = 1,nelemabs

      ispec = numabs(ispecabs)

! loop on all the coupling edges
      do inum = 1,num_fluid_poro_edges

! get the edge of the acoustic element
        ispec_acoustic = fluid_poro_acoustic_ispec(inum)
        iedge_acoustic = fluid_poro_acoustic_iedge(inum)

! if acoustic absorbing element and acoustic/poroelastic coupled element is the same
        if(ispec_acoustic == ispec) then

          if(iedge_acoustic == IBOTTOM) then
            jbegin_left(ispecabs) = 2
            jbegin_right(ispecabs) = 2
          endif

          if(iedge_acoustic == ITOP) then
            jend_left(ispecabs) = NGLLZ - 1
            jend_right(ispecabs) = NGLLZ - 1
          endif

          if(iedge_acoustic == ILEFT) then
            ibegin_bottom(ispecabs) = 2
            ibegin_top(ispecabs) = 2
          endif

          if(iedge_acoustic == IRIGHT) then
            iend_bottom(ispecabs) = NGLLX - 1
            iend_top(ispecabs) = NGLLX - 1
          endif

        endif

      enddo

    enddo

  endif


! determine if coupled elastic-poroelastic simulation
  coupled_elastic_poroelastic = any_elastic .and. any_poroelastic

! solid/porous edge detection
! the two elements forming an edge are already known (computed in meshfem2D),
! the common nodes forming the edge are computed here
  if(coupled_elastic_poroelastic) then

    if(TURN_ATTENUATION_ON .or. TURN_VISCATTENUATION_ON) &
                   stop 'Attenuation not supported for mixed elastic/poroelastic simulations'

    if ( myrank == 0 ) then
    print *
    print *,'Mixed elastic/poroelastic simulation'
    print *
    print *,'Beginning of solid/porous edge detection'
    endif

! define the edges of a given element
    i_begin(IBOTTOM) = 1
    j_begin(IBOTTOM) = 1
    i_end(IBOTTOM) = NGLLX
    j_end(IBOTTOM) = 1

    i_begin(IRIGHT) = NGLLX
    j_begin(IRIGHT) = 1
    i_end(IRIGHT) = NGLLX
    j_end(IRIGHT) = NGLLZ

    i_begin(ITOP) = NGLLX
    j_begin(ITOP) = NGLLZ
    i_end(ITOP) = 1
    j_end(ITOP) = NGLLZ

    i_begin(ILEFT) = 1
    j_begin(ILEFT) = NGLLZ
    i_end(ILEFT) = 1
    j_end(ILEFT) = 1

! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = ipoin1D
      ivalue_inverse(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,IBOTTOM) = 1
      jvalue_inverse(ipoin1D,IBOTTOM) = 1

      ivalue(ipoin1D,IRIGHT) = NGLLX
      ivalue_inverse(ipoin1D,IRIGHT) = NGLLX
      jvalue(ipoin1D,IRIGHT) = ipoin1D
      jvalue_inverse(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1

      ivalue(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,ITOP) = ipoin1D
      jvalue(ipoin1D,ITOP) = NGLLZ
      jvalue_inverse(ipoin1D,ITOP) = NGLLZ

      ivalue(ipoin1D,ILEFT) = 1
      ivalue_inverse(ipoin1D,ILEFT) = 1
      jvalue(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,ILEFT) = ipoin1D

    enddo


    do inum = 1, num_solid_poro_edges
       ispec_elastic =  solid_poro_elastic_ispec(inum)
       ispec_poroelastic =  solid_poro_poroelastic_ispec(inum)

! one element must be elastic and the other must be poroelastic
        if(ispec_elastic /= ispec_poroelastic .and. elastic(ispec_elastic) .and. &
                 poroelastic(ispec_poroelastic)) then

! loop on the four edges of the two elements
          do iedge_poroelastic = 1,NEDGES
            do iedge_elastic = 1,NEDGES

! store the matching topology if the two edges match in inverse order
              if(ibool(i_begin(iedge_poroelastic),j_begin(iedge_poroelastic),ispec_poroelastic) == &
                   ibool(i_end(iedge_elastic),j_end(iedge_elastic),ispec_elastic) .and. &
                   ibool(i_end(iedge_poroelastic),j_end(iedge_poroelastic),ispec_poroelastic) == &
                   ibool(i_begin(iedge_elastic),j_begin(iedge_elastic),ispec_elastic)) then
                 solid_poro_elastic_iedge(inum) = iedge_elastic
                 solid_poro_poroelastic_iedge(inum) = iedge_poroelastic
                endif

             enddo
          enddo

       endif

    enddo

! make sure solid/porous matching has been perfectly detected: check that the grid points
! have the same physical coordinates
! loop on all the coupling edges

    if ( myrank == 0 ) then
    print *,'Checking solid/porous edge topology...'
    endif

    do inum = 1,num_solid_poro_edges

! get the edge of the elastic element
      ispec_elastic = solid_poro_elastic_ispec(inum)
      iedge_elastic = solid_poro_elastic_iedge(inum)

! get the corresponding edge of the poroelastic element
      ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
      iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

! implement 1D coupling along the edge
      do ipoin1D = 1,NGLLX

! get point values for the poroelastic side, which matches our side in the inverse direction
        i = ivalue_inverse(ipoin1D,iedge_elastic)
        j = jvalue_inverse(ipoin1D,iedge_elastic)
        iglob = ibool(i,j,ispec_elastic)

! get point values for the elastic side
        i = ivalue(ipoin1D,iedge_poroelastic)
        j = jvalue(ipoin1D,iedge_poroelastic)
        iglob2 = ibool(i,j,ispec_poroelastic)

! if distance between the two points is not negligible, there is an error, since it should be zero
        if(sqrt((coord(1,iglob) - coord(1,iglob2))**2 + (coord(2,iglob) - coord(2,iglob2))**2) > TINYVAL) &
            call exit_MPI( 'error in solid/porous coupling buffer')

      enddo

    enddo

    if ( myrank == 0 ) then
    print *,'End of solid/porous edge detection'
    print *
    endif

  endif

! initiation
 if(any_poroelastic .and. anyabs) then
! loop on all the absorbing elements
    do ispecabs = 1,nelemabs
            jbegin_left_poro(ispecabs) = 1
            jbegin_right_poro(ispecabs) = 1

            jend_left_poro(ispecabs) = NGLLZ 
            jend_right_poro(ispecabs) = NGLLZ 

            ibegin_bottom_poro(ispecabs) = 1
            ibegin_top_poro(ispecabs) = 1

            iend_bottom_poro(ispecabs) = NGLLX
            iend_top_poro(ispecabs) = NGLLX
    enddo
 endif

! exclude common points between poroelastic absorbing edges and elastic/poroelastic matching interfaces
  if(coupled_elastic_poroelastic .and. anyabs) then

    if (myrank == 0) &
      print *,'excluding common points between poroelastic absorbing edges and elastic/poroelastic matching interfaces, if any'

! loop on all the absorbing elements
    do ispecabs = 1,nelemabs

      ispec = numabs(ispecabs)

! loop on all the coupling edges
      do inum = 1,num_solid_poro_edges

! get the edge of the acoustic element
        ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

! if poroelastic absorbing element and elastic/poroelastic coupled element is the same
        if(ispec_poroelastic == ispec) then

          if(iedge_poroelastic == IBOTTOM) then
            jbegin_left_poro(ispecabs) = 2
            jbegin_right_poro(ispecabs) = 2
          endif

          if(iedge_poroelastic == ITOP) then
            jend_left_poro(ispecabs) = NGLLZ - 1
            jend_right_poro(ispecabs) = NGLLZ - 1
          endif

          if(iedge_poroelastic == ILEFT) then
            ibegin_bottom_poro(ispecabs) = 2
            ibegin_top_poro(ispecabs) = 2
          endif

          if(iedge_poroelastic == IRIGHT) then
            iend_bottom_poro(ispecabs) = NGLLX - 1
            iend_top_poro(ispecabs) = NGLLX - 1
          endif

        endif

      enddo

    enddo

  endif


#ifdef USE_MPI
  if(OUTPUT_ENERGY) stop 'energy calculation only currently serial only, should add an MPI_REDUCE in parallel'
#endif
! open the file in which we will store the energy curve
  if(OUTPUT_ENERGY) open(unit=IENERGY,file='energy.gnu',status='unknown')

!
!----          s t a r t   t i m e   i t e r a t i o n s
!
  if (myrank == 0) write(IOUT,400)

! count elapsed wall-clock time
  call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! this fails if we cross the end of the month
  time_start = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
               60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0

  if(output_color_image) then
! to display the P-velocity model in background on color images
  allocate(vp_display(npoin))
  do ispec = 1,nspec

   if(poroelastic(ispec)) then
!get parameters of current spectral element
    phil = porosity(kmato(ispec))
    tortl = tortuosity(kmato(ispec))
!solid properties
    mul_s = poroelastcoef(2,1,kmato(ispec))
    kappal_s = poroelastcoef(3,1,kmato(ispec)) - 4._CUSTOM_REAL*mul_s/3._CUSTOM_REAL
    rhol_s = density(1,kmato(ispec))
!fluid properties
    kappal_f = poroelastcoef(1,2,kmato(ispec))
    rhol_f = density(2,kmato(ispec))
!frame properties
    mul_fr = poroelastcoef(2,3,kmato(ispec))
    kappal_fr = poroelastcoef(3,3,kmato(ispec)) - 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
    rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f
!Biot coefficients for the input phi
      D_biot = kappal_s*(1._CUSTOM_REAL + phil*(kappal_s/kappal_f - 1._CUSTOM_REAL))
      H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) + kappal_fr + 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
      C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
      M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)
      B_biot = H_biot - 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
! Approximated velocities (no viscous dissipation)
      afactor = rhol_bar - phil/tortl*rhol_f
      bfactor = H_biot + phil*rhol_bar/(tortl*rhol_f)*M_biot - TWO*phil/tortl*C_biot
      cfactor = phil/(tortl*rhol_f)*(H_biot*M_biot - C_biot*C_biot)
      cpIsquare = (bfactor + sqrt(bfactor*bfactor - 4._CUSTOM_REAL*afactor*cfactor))/(2._CUSTOM_REAL*afactor)
      cpIIsquare = (bfactor - sqrt(bfactor*bfactor - 4._CUSTOM_REAL*afactor*cfactor))/(2._CUSTOM_REAL*afactor)
      cssquare = mul_fr/afactor

! Approximated ratio r = amplitude "w" field/amplitude "s" field (no viscous dissipation)
! used later for kernels calculation
      gamma1 = H_biot - phil/tortl*C_biot
      gamma2 = C_biot - phil/tortl*M_biot
      gamma3 = phil/tortl*( M_biot*(afactor/rhol_f + phil/tortl) - C_biot)
      gamma4 = phil/tortl*( C_biot*(afactor/rhol_f + phil/tortl) - H_biot)
      ratio = HALF*(gamma1 - gamma3)/gamma4 + HALF*sqrt((gamma1-gamma3)**2/gamma4**2 + 4._CUSTOM_REAL * gamma2/gamma4)

    do j = 1,NGLLZ
      do i = 1,NGLLX
            vp_display(ibool(i,j,ispec)) = sqrt(cpIsquare)
      enddo
    enddo

   else
! get relaxed elastic parameters of current spectral element
    rhol = density(1,kmato(ispec))
    lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
    mul_relaxed = poroelastcoef(2,1,kmato(ispec))
    do j = 1,NGLLZ
      do i = 1,NGLLX
!--- if external medium, get elastic parameters of current grid point
          if(assign_external_model) then
            vp_display(ibool(i,j,ispec)) = vpext(i,j,ispec)
          else
            vp_display(ibool(i,j,ispec)) = sqrt((lambdal_relaxed + 2.d0*mul_relaxed) / rhol)
          endif
      enddo
    enddo
   endif !if(poroelastic(ispec)) then
  enddo

! getting velocity for each local pixels
  image_color_vp_display(:,:) = 0.d0

  do k = 1, nb_pixel_loc
    j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
    i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color
    image_color_vp_display(i,j) = vp_display(iglob_image_color(i,j))
  enddo

! assembling array image_color_vp_display on process zero for color output
#ifdef USE_MPI
  if (nproc > 1) then
    if (myrank == 0) then
      do iproc = 1, nproc-1
        call MPI_RECV(data_pixel_recv(1),nb_pixel_per_proc(iproc+1), MPI_DOUBLE_PRECISION, &
                iproc, 43, MPI_COMM_WORLD, request_mpi_status, ier)

        do k = 1, nb_pixel_per_proc(iproc+1)
          j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
          i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color
          image_color_vp_display(i,j) = data_pixel_recv(k)
        enddo
      enddo

    else
      do k = 1, nb_pixel_loc
        j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
        i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color
        data_pixel_send(k) = vp_display(iglob_image_color(i,j))
      enddo

      call MPI_SEND(data_pixel_send(1),nb_pixel_loc,MPI_DOUBLE_PRECISION, 0, 43, MPI_COMM_WORLD, ier)

    endif
  endif

#endif
  endif

! dummy allocation of plane wave arrays if they are unused (but still need to exist because
! they are used as arguments to subroutines)
  if(.not. over_critical_angle) then
    allocate(v0x_left(1,NSTEP))
    allocate(v0z_left(1,NSTEP))
    allocate(t0x_left(1,NSTEP))
    allocate(t0z_left(1,NSTEP))

    allocate(v0x_right(1,NSTEP))
    allocate(v0z_right(1,NSTEP))
    allocate(t0x_right(1,NSTEP))
    allocate(t0z_right(1,NSTEP))

    allocate(v0x_bot(1,NSTEP))
    allocate(v0z_bot(1,NSTEP))
    allocate(t0x_bot(1,NSTEP))
    allocate(t0z_bot(1,NSTEP))
  endif

! initialize variables for writing seismograms
  seismo_offset = 0
  seismo_current = 0

! Precompute Runge Kutta coefficients if viscous attenuation
  if(TURN_VISCATTENUATION_ON) then
    theta_e = (sqrt(Q0**2+1.d0) +1.d0)/(2.d0*pi*freq0*Q0)
    theta_s = (sqrt(Q0**2+1.d0) -1.d0)/(2.d0*pi*freq0*Q0)

    thetainv = - 1.d0 / theta_s
    alphaval = 1.d0 + deltat*thetainv + deltat**2*thetainv**2 / 2.d0 + &
      deltat**3*thetainv**3 / 6.d0 + deltat**4*thetainv**4 / 24.d0
    betaval = deltat / 2.d0 + deltat**2*thetainv / 3.d0 + deltat**3*thetainv**2 / 8.d0 + deltat**4*thetainv**3 / 24.d0
    gammaval = deltat / 2.d0 + deltat**2*thetainv / 6.d0 + deltat**3*thetainv**2 / 24.d0
   print*,'************************************************************'
   print*,'****** Visco attenuation coefficients (poroelastic) ********'
   print*,'theta_e = ', theta_e
   print*,'theta_s = ', theta_s
   print*,'alpha = ', alphaval
   print*,'beta = ', betaval
   print*,'gamma = ', gammaval
   print*,'************************************************************'

! initialize memory variables for attenuation
    viscox(:,:,:) = 0.d0
    viscoz(:,:,:) = 0.d0
    rx_viscous(:,:,:) = 0.d0
    rz_viscous(:,:,:) = 0.d0

  endif

! allocate arrays for postscript output
#ifdef USE_MPI
  if(modelvect) then
  d1_coorg_recv_ps_velocity_model=2
  call mpi_allreduce(nspec,d2_coorg_recv_ps_velocity_model,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
  d2_coorg_recv_ps_velocity_model=d2_coorg_recv_ps_velocity_model*((NGLLX-subsamp)/subsamp)*((NGLLX-subsamp)/subsamp)*4
  d1_RGB_recv_ps_velocity_model=1
  call mpi_allreduce(nspec,d2_RGB_recv_ps_velocity_model,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
  d2_RGB_recv_ps_velocity_model=d2_RGB_recv_ps_velocity_model*((NGLLX-subsamp)/subsamp)*((NGLLX-subsamp)/subsamp)*4
  else
  d1_coorg_recv_ps_velocity_model=1
  d2_coorg_recv_ps_velocity_model=1
  d1_RGB_recv_ps_velocity_model=1
  d2_RGB_recv_ps_velocity_model=1
  endif

  d1_coorg_send_ps_element_mesh=2
  if ( ngnod == 4 ) then
    if ( numbers == 1 ) then
      d2_coorg_send_ps_element_mesh=nspec*5
      if ( colors == 1 ) then
        d1_color_send_ps_element_mesh=2*nspec
      else
        d1_color_send_ps_element_mesh=1*nspec
      endif
    else
      d2_coorg_send_ps_element_mesh=nspec*6
      if ( colors == 1 ) then
        d1_color_send_ps_element_mesh=1*nspec
      endif
    endif
  else
    if ( numbers == 1 ) then
      d2_coorg_send_ps_element_mesh=nspec*((pointsdisp-1)*3+max(0,pointsdisp-2)+1+1)
      if ( colors == 1 ) then
        d1_color_send_ps_element_mesh=2*nspec
      else
        d1_color_send_ps_element_mesh=1*nspec
      endif
    else
      d2_coorg_send_ps_element_mesh=nspec*((pointsdisp-1)*3+max(0,pointsdisp-2)+1)
      if ( colors == 1 ) then
        d1_color_send_ps_element_mesh=1*nspec
      endif
    endif
  endif

call mpi_allreduce(d1_coorg_send_ps_element_mesh,d1_coorg_recv_ps_element_mesh,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
call mpi_allreduce(d2_coorg_send_ps_element_mesh,d2_coorg_recv_ps_element_mesh,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
call mpi_allreduce(d1_color_send_ps_element_mesh,d1_color_recv_ps_element_mesh,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

  d1_coorg_send_ps_abs=4
  d2_coorg_send_ps_abs=4*nelemabs
call mpi_allreduce(d1_coorg_send_ps_abs,d1_coorg_recv_ps_abs,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
call mpi_allreduce(d2_coorg_send_ps_abs,d2_coorg_recv_ps_abs,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

  d1_coorg_send_ps_free_surface=4
  d2_coorg_send_ps_free_surface=4*nelem_acoustic_surface
call mpi_allreduce(d1_coorg_send_ps_free_surface,d1_coorg_recv_ps_free_surface,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
call mpi_allreduce(d2_coorg_send_ps_free_surface,d2_coorg_recv_ps_free_surface,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

  d1_coorg_send_ps_vector_field=8
  if(interpol) then
    if(plot_lowerleft_corner_only) then
      d2_coorg_send_ps_vector_field=nspec*1*1
    else
      d2_coorg_send_ps_vector_field=nspec*pointsdisp*pointsdisp
    endif
  else
    d2_coorg_send_ps_vector_field=npoin
  endif
call mpi_allreduce(d1_coorg_send_ps_vector_field,d1_coorg_recv_ps_vector_field,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
call mpi_allreduce(d2_coorg_send_ps_vector_field,d2_coorg_recv_ps_vector_field,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)


#else
  d1_coorg_recv_ps_velocity_model=1
  d2_coorg_recv_ps_velocity_model=1
  d1_RGB_recv_ps_velocity_model=1
  d2_RGB_recv_ps_velocity_model=1

  d1_coorg_send_ps_element_mesh=1
  d2_coorg_send_ps_element_mesh=1
  d1_coorg_recv_ps_element_mesh=1
  d2_coorg_recv_ps_element_mesh=1
  d1_color_send_ps_element_mesh=1
  d1_color_recv_ps_element_mesh=1

  d1_coorg_send_ps_abs=1
  d2_coorg_send_ps_abs=1
  d1_coorg_recv_ps_abs=1
  d2_coorg_recv_ps_abs=1
  d1_coorg_send_ps_free_surface=1
  d2_coorg_send_ps_free_surface=1
  d1_coorg_recv_ps_free_surface=1
  d2_coorg_recv_ps_free_surface=1

  d1_coorg_send_ps_vector_field=1
  d2_coorg_send_ps_vector_field=1
  d1_coorg_recv_ps_vector_field=1
  d2_coorg_recv_ps_vector_field=1

#endif
  d1_coorg_send_ps_velocity_model=2
  d2_coorg_send_ps_velocity_model=nspec*((NGLLX-subsamp)/subsamp)*((NGLLX-subsamp)/subsamp)*4
  d1_RGB_send_ps_velocity_model=1
  d2_RGB_send_ps_velocity_model=nspec*((NGLLX-subsamp)/subsamp)*((NGLLX-subsamp)/subsamp)

  allocate(coorg_send_ps_velocity_model(d1_coorg_send_ps_velocity_model,d2_coorg_send_ps_velocity_model))
  allocate(RGB_send_ps_velocity_model(d1_RGB_send_ps_velocity_model,d2_RGB_send_ps_velocity_model))

  allocate(coorg_recv_ps_velocity_model(d1_coorg_recv_ps_velocity_model,d2_coorg_recv_ps_velocity_model))
  allocate(RGB_recv_ps_velocity_model(d1_RGB_recv_ps_velocity_model,d2_RGB_recv_ps_velocity_model))

  allocate(coorg_send_ps_element_mesh(d1_coorg_send_ps_element_mesh,d2_coorg_send_ps_element_mesh))
  allocate(coorg_recv_ps_element_mesh(d1_coorg_recv_ps_element_mesh,d2_coorg_recv_ps_element_mesh))
  allocate(color_send_ps_element_mesh(d1_color_send_ps_element_mesh))
  allocate(color_recv_ps_element_mesh(d1_color_recv_ps_element_mesh))

  allocate(coorg_send_ps_abs(d1_coorg_send_ps_abs,d2_coorg_send_ps_abs))
  allocate(coorg_recv_ps_abs(d1_coorg_recv_ps_abs,d2_coorg_recv_ps_abs))

  allocate(coorg_send_ps_free_surface(d1_coorg_send_ps_free_surface,d2_coorg_send_ps_free_surface))
  allocate(coorg_recv_ps_free_surface(d1_coorg_recv_ps_free_surface,d2_coorg_recv_ps_free_surface))

  allocate(coorg_send_ps_vector_field(d1_coorg_send_ps_vector_field,d2_coorg_send_ps_vector_field))
  allocate(coorg_recv_ps_vector_field(d1_coorg_recv_ps_vector_field,d2_coorg_recv_ps_vector_field))

! *********************************************************
! ************* MAIN LOOP OVER THE TIME STEPS *************
! *********************************************************

#ifdef USE_MPI
! add a barrier if we generate traces of the run for analysis with "ParaVer"
  if(GENERATE_PARAVER_TRACES) call MPI_BARRIER(MPI_COMM_WORLD,ier)
#endif

  do it = 1,NSTEP

! update position in seismograms
    seismo_current = seismo_current + 1

! compute current time
    time = (it-1)*deltat

! update displacement using finite-difference time scheme (Newmark)
    if(any_elastic) then
      displ_elastic = displ_elastic + deltat*veloc_elastic + deltatsquareover2*accel_elastic
      veloc_elastic = veloc_elastic + deltatover2*accel_elastic
      accel_elastic = ZERO

     if(isolver == 2) then ! Adjoint calculation
      b_displ_elastic = b_displ_elastic + b_deltat*b_veloc_elastic + b_deltatsquareover2*b_accel_elastic
      b_veloc_elastic = b_veloc_elastic + b_deltatover2*b_accel_elastic
      b_accel_elastic = ZERO
     endif
    endif

    if(any_poroelastic) then
!for the solid
      displs_poroelastic = displs_poroelastic + deltat*velocs_poroelastic + deltatsquareover2*accels_poroelastic
      velocs_poroelastic = velocs_poroelastic + deltatover2*accels_poroelastic
      accels_poroelastic = ZERO
!for the fluid 
      displw_poroelastic = displw_poroelastic + deltat*velocw_poroelastic + deltatsquareover2*accelw_poroelastic
      velocw_poroelastic = velocw_poroelastic + deltatover2*accelw_poroelastic
      accelw_poroelastic = ZERO

     if(isolver == 2) then ! Adjoint calculation
!for the solid
      b_displs_poroelastic = b_displs_poroelastic + b_deltat*b_velocs_poroelastic + b_deltatsquareover2*b_accels_poroelastic
      b_velocs_poroelastic = b_velocs_poroelastic + b_deltatover2*b_accels_poroelastic
      b_accels_poroelastic = ZERO
!for the fluid 
      b_displw_poroelastic = b_displw_poroelastic + b_deltat*b_velocw_poroelastic + b_deltatsquareover2*b_accelw_poroelastic
      b_velocw_poroelastic = b_velocw_poroelastic + b_deltatover2*b_accelw_poroelastic
      b_accelw_poroelastic = ZERO
     endif
    endif

!--------------------------------------------------------------------------------------------
! implement viscous attenuation for poroelastic media
!
  if(TURN_VISCATTENUATION_ON .and. any_poroelastic) then
! update memory variables with fourth-order Runge-Kutta time scheme for attenuation
! loop over spectral elements

  do ispec = 1,nspec

    etal_f = poroelastcoef(2,2,kmato(ispec))
    permlxx = permeability(1,kmato(ispec))
    permlxz = permeability(2,kmato(ispec))
    permlzz = permeability(3,kmato(ispec))

 ! calcul of the inverse of k

    detk = permlxx*permlzz - permlxz*permlxz

    if(detk /= ZERO) then
     invpermlxx = permlzz/detk
     invpermlxz = -permlxz/detk
     invpermlzz = permlxx/detk
    else
      stop 'Permeability matrix is not inversible'
    endif

! relaxed viscous coef
          bl_relaxed(1) = etal_f*invpermlxx
          bl_relaxed(2) = etal_f*invpermlxz
          bl_relaxed(3) = etal_f*invpermlzz

  do j=1,NGLLZ
  do i=1,NGLLX

    iglob = ibool(i,j,ispec)

   viscox_loc(i,j) = velocw_poroelastic(1,iglob)*bl_relaxed(1) + &
                               velocw_poroelastic(2,iglob)*bl_relaxed(2)
   viscoz_loc(i,j) = velocw_poroelastic(1,iglob)*bl_relaxed(2) + &
                               velocw_poroelastic(2,iglob)*bl_relaxed(3)

! evolution rx_viscous
  Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscox(i,j,ispec)
  Snp1 = - (1.d0 - theta_e/theta_s)/theta_s*viscox_loc(i,j)
  rx_viscous(i,j,ispec) = alphaval * rx_viscous(i,j,ispec) + betaval * Sn + gammaval * Snp1

! evolution rz_viscous
  Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscoz(i,j,ispec)
  Snp1 = - (1.d0 - theta_e/theta_s)/theta_s*viscoz_loc(i,j)
  rz_viscous(i,j,ispec) = alphaval * rz_viscous(i,j,ispec) + betaval * Sn + gammaval * Snp1


  enddo
  enddo

! save visco for Runge-Kutta scheme
     viscox(:,:,ispec) = viscox_loc(:,:)
     viscoz(:,:,ispec) = viscoz_loc(:,:)

  enddo   ! end of spectral element loop
  endif ! end of viscous attenuation for porous media

!-----------------------------------------

    if(any_acoustic) then

      potential_acoustic = potential_acoustic + deltat*potential_dot_acoustic + deltatsquareover2*potential_dot_dot_acoustic
      potential_dot_acoustic = potential_dot_acoustic + deltatover2*potential_dot_dot_acoustic
      potential_dot_dot_acoustic = ZERO

    if(isolver == 2) then ! Adjoint calculation
      b_potential_acoustic = b_potential_acoustic + b_deltat*b_potential_dot_acoustic + &
                             b_deltatsquareover2*b_potential_dot_dot_acoustic
      b_potential_dot_acoustic = b_potential_dot_acoustic + b_deltatover2*b_potential_dot_dot_acoustic
      b_potential_dot_dot_acoustic = ZERO
    endif

! free surface for an acoustic medium
      if ( nelem_acoustic_surface > 0 ) then
        call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
           potential_acoustic,acoustic_surface,ibool,nelem_acoustic_surface,npoin,nspec)

   if(isolver == 2) then ! Adjoint calculation
    call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                b_potential_acoustic,acoustic_surface,ibool,nelem_acoustic_surface,npoin,nspec)
   endif
      endif

! *********************************************************
! ************* compute forces for the acoustic elements
! *********************************************************

    call compute_forces_acoustic(npoin,nspec,nelemabs,numat,it,NSTEP, &
               anyabs,assign_external_model,ibool,kmato,numabs, &
               elastic,poroelastic,codeabs,potential_dot_dot_acoustic,potential_dot_acoustic, &
               potential_acoustic,b_potential_dot_dot_acoustic,b_potential_acoustic, &
               density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,rhoext,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
               jbegin_left,jend_left,jbegin_right,jend_right,isolver,save_forward,b_absorb_acoustic_left,&
               b_absorb_acoustic_right,b_absorb_acoustic_bottom,&
               b_absorb_acoustic_top,nspec_xmin,nspec_xmax,&
               nspec_zmin,nspec_zmax,ib_xmin,ib_xmax,ib_zmin,ib_zmax,kappa_ac_k)

    if(anyabs .and. save_forward .and. isolver == 1) then

!--- left absorbing boundary
      if(nspec_xmin >0) then
      do ispec = 1,nspec_xmin
         do i=1,NGLLZ
     write(65) b_absorb_acoustic_left(i,ispec,it)
         enddo
      enddo
      endif

!--- right absorbing boundary
      if(nspec_xmax >0) then
      do ispec = 1,nspec_xmax
         do i=1,NGLLZ
     write(66) b_absorb_acoustic_right(i,ispec,it)
         enddo
      enddo
      endif

!--- bottom absorbing boundary
      if(nspec_zmin >0) then
      do ispec = 1,nspec_zmin
         do i=1,NGLLX
     write(67) b_absorb_acoustic_bottom(i,ispec,it)
         enddo
      enddo
      endif

!--- top absorbing boundary
      if(nspec_zmax >0) then
      do ispec = 1,nspec_zmax
         do i=1,NGLLX
     write(68) b_absorb_acoustic_top(i,ispec,it)
         enddo
      enddo
      endif

    endif ! if(anyabs .and. save_forward .and. isolver == 1)

 endif ! end of test if any acoustic element

! *********************************************************
! ************* add coupling with the elastic side
! *********************************************************

  if(coupled_acoustic_elastic) then

! loop on all the coupling edges
      do inum = 1,num_fluid_solid_edges

! get the edge of the acoustic element
        ispec_acoustic = fluid_solid_acoustic_ispec(inum)
        iedge_acoustic = fluid_solid_acoustic_iedge(inum)

! get the corresponding edge of the elastic element
        ispec_elastic = fluid_solid_elastic_ispec(inum)
        iedge_elastic = fluid_solid_elastic_iedge(inum)

! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

! get point values for the elastic side, which matches our side in the inverse direction
          i = ivalue_inverse(ipoin1D,iedge_elastic)
          j = jvalue_inverse(ipoin1D,iedge_elastic)
          iglob = ibool(i,j,ispec_elastic)

          displ_x = displ_elastic(1,iglob)
          displ_z = displ_elastic(2,iglob)

          if(isolver == 2) then
          b_displ_x = b_displ_elastic(1,iglob)
          b_displ_z = b_displ_elastic(2,iglob)
          endif

! get point values for the acoustic side
          i = ivalue(ipoin1D,iedge_acoustic)
          j = jvalue(ipoin1D,iedge_acoustic)
          iglob = ibool(i,j,ispec_acoustic)

! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
! Blackwell Science, page 110, equation (4.60).
          if(iedge_acoustic == ITOP)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D
          weight = jacobian1D * wxgll(i)
          elseif(iedge_acoustic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
          weight = jacobian1D * wxgll(i)
          elseif(iedge_acoustic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
          weight = jacobian1D * wzgll(j)
          elseif(iedge_acoustic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
          weight = jacobian1D * wzgll(j)
          endif

! compute dot product
          displ_n = displ_x*nx + displ_z*nz

! formulation with generalized potential
          weight = jacobian1D * wxgll(i)

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

          if(isolver == 2) then
          b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) +&
                      weight*(b_displ_x*nx + b_displ_z*nz)
          endif !if(isolver == 2) then

        enddo

      enddo

   endif

! *********************************************************
! ************* add coupling with the poroelastic side
! *********************************************************

    if(coupled_acoustic_poroelastic) then

! loop on all the coupling edges
      do inum = 1,num_fluid_poro_edges

! get the edge of the acoustic element
        ispec_acoustic = fluid_poro_acoustic_ispec(inum)
        iedge_acoustic = fluid_poro_acoustic_iedge(inum)

! get the corresponding edge of the poroelastic element
        ispec_poroelastic = fluid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = fluid_poro_poroelastic_iedge(inum)

! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

! get point values for the poroelastic side, which matches our side in the inverse direction
          i = ivalue_inverse(ipoin1D,iedge_poroelastic)
          j = jvalue_inverse(ipoin1D,iedge_poroelastic)
          iglob = ibool(i,j,ispec_poroelastic)

          displ_x = displs_poroelastic(1,iglob)
          displ_z = displs_poroelastic(2,iglob)

          phil = porosity(kmato(ispec_poroelastic))
          displw_x = displw_poroelastic(1,iglob)
          displw_z = displw_poroelastic(2,iglob)

          if(isolver == 2) then
          b_displ_x = b_displs_poroelastic(1,iglob)
          b_displ_z = b_displs_poroelastic(2,iglob)

          b_displw_x = b_displw_poroelastic(1,iglob)
          b_displw_z = b_displw_poroelastic(2,iglob)
          endif

! get point values for the acoustic side
! get point values for the acoustic side
          i = ivalue(ipoin1D,iedge_acoustic)
          j = jvalue(ipoin1D,iedge_acoustic)
          iglob = ibool(i,j,ispec_acoustic)

! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
! Blackwell Science, page 110, equation (4.60).
          if(iedge_acoustic == ITOP)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D
          weight = jacobian1D * wxgll(i)
          elseif(iedge_acoustic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
          weight = jacobian1D * wxgll(i)
          elseif(iedge_acoustic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
          weight = jacobian1D * wzgll(j)
          elseif(iedge_acoustic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
          weight = jacobian1D * wzgll(j)
          endif

! compute dot product [u_s + w]*n
          displ_n = (displ_x + displw_x)*nx + (displ_z + displw_z)*nz

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

          if(isolver == 2) then
          b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) +&
                    weight*((b_displ_x + b_displw_x)*nx + (b_displ_z + b_displw_z)*nz)
          endif !if(isolver == 2) then

        enddo

      enddo

    endif


! assembling potential_dot_dot for acoustic elements
#ifdef USE_MPI
  if ( nproc > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_ac(potential_dot_dot_acoustic,npoin, &
           ninterface, ninterface_acoustic,inum_interfaces_acoustic, &
           max_interface_size, max_ibool_interfaces_size_ac,&
           ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
           tab_requests_send_recv_acoustic,buffer_send_faces_vector_ac, &
           buffer_recv_faces_vector_ac, my_neighbours)
  endif

    if ( nproc > 1 .and. any_acoustic .and. ninterface_acoustic > 0 .and. isolver == 2) then
    call assemble_MPI_vector_ac(b_potential_dot_dot_acoustic,npoin, &
           ninterface, ninterface_acoustic,inum_interfaces_acoustic, &
           max_interface_size, max_ibool_interfaces_size_ac,&
           ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
           tab_requests_send_recv_acoustic,buffer_send_faces_vector_ac, &
           buffer_recv_faces_vector_ac, my_neighbours)
  endif
#endif


! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

  if(any_acoustic) then

! --- add the source
    if(.not. initialfield) then

    do i_source=1,NSOURCE
! if this processor carries the source and the source element is acoustic
      if (is_proc_source(i_source) == 1 .and. .not. elastic(ispec_selected_source(i_source)) .and. &
        .not. poroelastic(ispec_selected_source(i_source))) then
! collocated force
! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
! the sign is negative because pressure p = - Chi_dot_dot therefore we need
! to add minus the source to Chi_dot_dot to get plus the source in pressure
        if(source_type(i_source) == 1) then
      if(isolver == 1) then  ! forward wavefield
          potential_dot_dot_acoustic(iglob_source(i_source)) = potential_dot_dot_acoustic(iglob_source(i_source)) &
                 - source_time_function(i_source,it)
      else                   ! backward wavefield
      b_potential_dot_dot_acoustic(iglob_source(i_source)) = b_potential_dot_dot_acoustic(iglob_source(i_source)) &
                 - source_time_function(i_source,NSTEP-it+1)
      endif
! moment tensor
        else if(source_type(i_source) == 2) then
          call exit_MPI('cannot have moment tensor source in acoustic element')
        endif
      endif ! if this processor carries the source and the source element is acoustic
    enddo ! do i_source=1,NSOURCE

    if(isolver == 2) then   ! adjoint wavefield
      irec_local = 0
      do irec = 1,nrec
!   add the source (only if this proc carries the source)
      if (myrank == which_proc_receiver(irec) .and. .not. elastic(ispec_selected_rec(irec)) .and. &
         .not. poroelastic(ispec_selected_rec(irec))) then
      irec_local = irec_local + 1
! add source array
      do j=1,NGLLZ
        do i=1,NGLLX
      iglob = ibool(i,j,ispec_selected_rec(irec))
      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
          adj_sourcearrays(irec_local,NSTEP-it+1,1,i,j)
        enddo
      enddo
      endif ! if this processor carries the adjoint source
      enddo ! irec = 1,nrec
    endif ! isolver == 2 adjoint wavefield

    endif ! if not using an initial field

    potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic
    potential_dot_acoustic = potential_dot_acoustic + deltatover2*potential_dot_dot_acoustic

! free surface for an acoustic medium
    if ( nelem_acoustic_surface > 0 ) then
    call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                potential_acoustic,acoustic_surface,ibool,nelem_acoustic_surface,npoin,nspec)

   if(isolver == 2) then
    call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                b_potential_acoustic,acoustic_surface,ibool,nelem_acoustic_surface,npoin,nspec)
   endif

    endif
  endif

  if(any_acoustic .and. isolver == 2) then ! kernels calculation
      do iglob = 1,npoin
            rho_ac_k(iglob) = potential_dot_dot_acoustic(iglob)*b_potential_acoustic(iglob)
      enddo
  endif


! ****************************************************************************************
!   If coupling elastic/poroelastic domain, average some arrays at the interface first
! ****************************************************************************************
    if(coupled_elastic_poroelastic) then

! loop on all the coupling edges
      do inum = 1,num_solid_poro_edges

! get the edge of the elastic element
        ispec_elastic = solid_poro_elastic_ispec(inum)
        iedge_elastic = solid_poro_elastic_iedge(inum)

! get the corresponding edge of the poroelastic element
        ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = solid_poro_poroelastic_iedge(inum)
    phil = porosity(kmato(ispec_poroelastic))
! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

! get point values for the poroelastic side, which matches our side in the inverse direction
          i = ivalue(ipoin1D,iedge_poroelastic)
          j = jvalue(ipoin1D,iedge_poroelastic)
          iglob = ibool(i,j,ispec_poroelastic)

! get point values for the elastic side
          ii2 = ivalue_inverse(ipoin1D,iedge_elastic)
          jj2 = jvalue_inverse(ipoin1D,iedge_elastic)
          iglob2 = ibool(ii2,jj2,ispec_elastic)

     if(iglob /= iglob2) &
            call exit_MPI( 'error in solid/porous iglob detection')

           displ(1,iglob)=(displs_poroelastic(1,iglob) &
                              +displ_elastic(1,iglob2))/2.d0
           displ(2,iglob)=(displs_poroelastic(2,iglob) &
                              +displ_elastic(2,iglob2))/2.d0

           veloc(1,iglob)=(velocs_poroelastic(1,iglob) &
                              +veloc_elastic(1,iglob2))/2.d0
           veloc(2,iglob)=(velocs_poroelastic(2,iglob) &
                              +veloc_elastic(2,iglob2))/2.d0

        enddo
      enddo

! loop on all the coupling edges
      do inum = 1,num_solid_poro_edges

! get the corresponding edge of the poroelastic element
        ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = solid_poro_poroelastic_iedge(inum)
! imnplement 1D coupling along the edge
        do ipoin1D = 1,NGLLX
! get point values for the poroelastic side, which matches our side in the inverse direction
          i = ivalue(ipoin1D,iedge_poroelastic)
          j = jvalue(ipoin1D,iedge_poroelastic)
          iglob = ibool(i,j,ispec_poroelastic)

           displs_poroelastic(1,iglob)=displ(1,iglob)
           displs_poroelastic(2,iglob)=displ(2,iglob)

           displ_elastic(1,iglob)=displ(1,iglob)
           displ_elastic(2,iglob)=displ(2,iglob)

           velocs_poroelastic(1,iglob)=veloc(1,iglob)
           velocs_poroelastic(2,iglob)=veloc(2,iglob)

           veloc_elastic(1,iglob)=veloc(1,iglob)
           veloc_elastic(2,iglob)=veloc(2,iglob)

           displw_poroelastic(1,iglob)=ZERO
           displw_poroelastic(2,iglob)=ZERO

           velocw_poroelastic(1,iglob)=ZERO
           velocw_poroelastic(2,iglob)=ZERO
        enddo
      enddo

    endif

! *********************************************************
! ************* main solver for the elastic elements
! *********************************************************

 if(any_elastic) then
    call compute_forces_elastic(npoin,nspec,myrank,nelemabs,numat, &
               ispec_selected_source,ispec_selected_rec,is_proc_source,which_proc_receiver, &
               source_type,it,NSTEP,anyabs,assign_external_model, &
               initialfield,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,angleforce,deltatcube, &
               deltatfourth,twelvedeltat,fourdeltatsquare,ibool,kmato,numabs,elastic,codeabs, &
               accel_elastic,veloc_elastic,displ_elastic,b_accel_elastic,b_displ_elastic, &
               density,poroelastcoef,xix,xiz,gammax,gammaz, &
               jacobian,vpext,vsext,rhoext,source_time_function,sourcearray,adj_sourcearrays, &
               e1,e11,e13,dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n, &
               dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll,inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,Mu_nu1,Mu_nu2,N_SLS, &
               deltat,coord,add_Bielak_conditions, x0_source, z0_source, &
               A_plane, B_plane, C_plane, angleforce_refl, c_inc, c_refl, time_offset, f0(1),&
               v0x_left(1,it),v0z_left(1,it),v0x_right(1,it),v0z_right(1,it),v0x_bot(1,it),v0z_bot(1,it), &
               t0x_left(1,it),t0z_left(1,it),t0x_right(1,it),t0z_right(1,it),t0x_bot(1,it),t0z_bot(1,it), &
               count_left,count_right,count_bot,over_critical_angle,NSOURCE,nrec,isolver,save_forward,b_absorb_elastic_left,&
               b_absorb_elastic_right,b_absorb_elastic_bottom,b_absorb_elastic_top,nspec_xmin,nspec_xmax,&
               nspec_zmin,nspec_zmax,ib_xmin,ib_xmax,ib_zmin,ib_zmax,mu_k,kappa_k)

    if(anyabs .and. save_forward .and. isolver == 1) then
!--- left absorbing boundary
      if(nspec_xmin >0) then
      do ispec = 1,nspec_xmin
       do id =1,2
         do i=1,NGLLZ
     write(35) b_absorb_elastic_left(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

!--- right absorbing boundary
      if(nspec_xmax >0) then
      do ispec = 1,nspec_xmax
       do id =1,2
         do i=1,NGLLZ
     write(36) b_absorb_elastic_right(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

!--- bottom absorbing boundary
      if(nspec_zmin >0) then
      do ispec = 1,nspec_zmin
       do id =1,2
         do i=1,NGLLX
     write(37) b_absorb_elastic_bottom(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

!--- top absorbing boundary
      if(nspec_zmax >0) then
      do ispec = 1,nspec_zmax
       do id =1,2
         do i=1,NGLLX
     write(38) b_absorb_elastic_top(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

    endif ! if(anyabs .and. save_forward .and. isolver == 1)

  endif !if(any_elastic)

! *********************************************************
! ************* add coupling with the acoustic side
! *********************************************************

    if(coupled_acoustic_elastic) then

! loop on all the coupling edges
      do inum = 1,num_fluid_solid_edges

! get the edge of the acoustic element
        ispec_acoustic = fluid_solid_acoustic_ispec(inum)
        iedge_acoustic = fluid_solid_acoustic_iedge(inum)

! get the corresponding edge of the elastic element
        ispec_elastic = fluid_solid_elastic_ispec(inum)
        iedge_elastic = fluid_solid_elastic_iedge(inum)

! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

! get point values for the acoustic side, which matches our side in the inverse direction
          i = ivalue_inverse(ipoin1D,iedge_acoustic)
          j = jvalue_inverse(ipoin1D,iedge_acoustic)
          iglob = ibool(i,j,ispec_acoustic)

! compute pressure on the fluid/solid edge
          pressure = - potential_dot_dot_acoustic(iglob)
          if(isolver == 2) then
          b_pressure = - b_potential_dot_dot_acoustic(iglob)
          endif
! get point values for the elastic side
          ii2 = ivalue(ipoin1D,iedge_elastic)
          jj2 = jvalue(ipoin1D,iedge_elastic)
          iglob = ibool(ii2,jj2,ispec_elastic)

! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
! Blackwell Science, page 110, equation (4.60).
          if(iedge_acoustic == ITOP)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D
          weight = jacobian1D * wxgll(i)
          elseif(iedge_acoustic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
          weight = jacobian1D * wxgll(i)
          elseif(iedge_acoustic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
          weight = jacobian1D * wzgll(j)
          elseif(iedge_acoustic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
          weight = jacobian1D * wzgll(j)
          endif

          accel_elastic(1,iglob) = accel_elastic(1,iglob) + weight*nx*pressure
          accel_elastic(2,iglob) = accel_elastic(2,iglob) + weight*nz*pressure

          if(isolver == 2) then
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) + weight*nx*b_pressure
          b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) + weight*nz*b_pressure
          endif !if(isolver == 2) then

        enddo

      enddo

    endif

! ****************************************************************************
! ************* add coupling with the poroelastic side
! ****************************************************************************
    if(coupled_elastic_poroelastic) then

! loop on all the coupling edges
      do inum = 1,num_solid_poro_edges

! get the edge of the elastic element
        ispec_elastic = solid_poro_elastic_ispec(inum)
        iedge_elastic = solid_poro_elastic_iedge(inum)

! get the corresponding edge of the poroelastic element
        ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

! get point values for the poroelastic side, which matches our side in the inverse direction
          i = ivalue_inverse(ipoin1D,iedge_poroelastic)
          j = jvalue_inverse(ipoin1D,iedge_poroelastic)
          iglob = ibool(i,j,ispec_poroelastic)

! get poroelastic domain paramters
    phil = porosity(kmato(ispec_poroelastic))
    tortl = tortuosity(kmato(ispec_poroelastic))
!solid properties
    mul_s = poroelastcoef(2,1,kmato(ispec_poroelastic))
    kappal_s = poroelastcoef(3,1,kmato(ispec_poroelastic)) - 4._CUSTOM_REAL*mul_s/3._CUSTOM_REAL
    rhol_s = density(1,kmato(ispec_poroelastic))
!fluid properties
    kappal_f = poroelastcoef(1,2,kmato(ispec_poroelastic))
    rhol_f = density(2,kmato(ispec_poroelastic))
!frame properties
    mul_fr = poroelastcoef(2,3,kmato(ispec_poroelastic))
    kappal_fr = poroelastcoef(3,3,kmato(ispec_poroelastic)) - 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
    rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f
!Biot coefficients for the input phi
      D_biot = kappal_s*(1._CUSTOM_REAL + phil*(kappal_s/kappal_f - 1._CUSTOM_REAL))
      H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) + &
                kappal_fr + 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
      C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
      M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)
      mul_G = mul_fr
      lambdal_G = H_biot - 2._CUSTOM_REAL*mul_fr
      lambdalplus2mul_G = lambdal_G + TWO*mul_G

! derivative along x and along z for u_s and w
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          dwx_dxi = ZERO
          dwz_dxi = ZERO

          dwx_dgamma = ZERO
          dwz_dgamma = ZERO

          if(isolver == 2) then
          b_dux_dxi = ZERO
          b_duz_dxi = ZERO

          b_dux_dgamma = ZERO
          b_duz_dgamma = ZERO

          b_dwx_dxi = ZERO
          b_dwz_dxi = ZERO

          b_dwx_dgamma = ZERO
          b_dwz_dgamma = ZERO
          endif

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)

            dwx_dxi = dwx_dxi + displw_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            dwz_dxi = dwz_dxi + displw_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            dwx_dgamma = dwx_dgamma + displw_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            dwz_dgamma = dwz_dgamma + displw_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            if(isolver == 2) then
            b_dux_dxi = b_dux_dxi + b_displs_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            b_duz_dxi = b_duz_dxi + b_displs_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            b_dux_dgamma = b_dux_dgamma + b_displs_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            b_duz_dgamma = b_duz_dgamma + b_displs_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)

            b_dwx_dxi = b_dwx_dxi + b_displw_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            b_dwz_dxi = b_dwz_dxi + b_displw_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
            b_dwx_dgamma = b_dwx_dgamma + b_displw_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            b_dwz_dgamma = b_dwz_dgamma + b_displw_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
            endif
          enddo

          xixl = xix(i,j,ispec_poroelastic)
          xizl = xiz(i,j,ispec_poroelastic)
          gammaxl = gammax(i,j,ispec_poroelastic)
          gammazl = gammaz(i,j,ispec_poroelastic)

! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          dwx_dxl = dwx_dxi*xixl + dwx_dgamma*gammaxl
          dwx_dzl = dwx_dxi*xizl + dwx_dgamma*gammazl

          dwz_dxl = dwz_dxi*xixl + dwz_dgamma*gammaxl
          dwz_dzl = dwz_dxi*xizl + dwz_dgamma*gammazl

          if(isolver == 2) then
          b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
          b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

          b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
          b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl

          b_dwx_dxl = b_dwx_dxi*xixl + b_dwx_dgamma*gammaxl
          b_dwx_dzl = b_dwx_dxi*xizl + b_dwx_dgamma*gammazl

          b_dwz_dxl = b_dwz_dxi*xixl + b_dwz_dgamma*gammaxl
          b_dwz_dzl = b_dwz_dxi*xizl + b_dwz_dgamma*gammazl
          endif
! compute stress tensor (include attenuation or anisotropy if needed)

! no attenuation
    sigma_xx = lambdalplus2mul_G*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
    sigma_xz = mul_G*(duz_dxl + dux_dzl)
    sigma_zz = lambdalplus2mul_G*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

    if(isolver == 2) then
    b_sigma_xx = lambdalplus2mul_G*b_dux_dxl + lambdal_G*b_duz_dzl + C_biot*(b_dwx_dxl + b_dwz_dzl)
    b_sigma_xz = mul_G*(b_duz_dxl + b_dux_dzl)
    b_sigma_zz = lambdalplus2mul_G*b_duz_dzl + lambdal_G*b_dux_dxl + C_biot*(b_dwx_dxl + b_dwz_dzl)
    endif
! get point values for the elastic domain, which matches our side in the inverse direction
          ii2 = ivalue(ipoin1D,iedge_elastic)
          jj2 = jvalue(ipoin1D,iedge_elastic)
          iglob = ibool(ii2,jj2,ispec_elastic)

! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
! Blackwell Science, page 110, equation (4.60).
          if(iedge_poroelastic == ITOP)then
            xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D
          weight = jacobian1D * wxgll(i)
          elseif(iedge_poroelastic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
          weight = jacobian1D * wxgll(i)
          elseif(iedge_poroelastic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
          weight = jacobian1D * wzgll(j)
          elseif(iedge_poroelastic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
          weight = jacobian1D * wzgll(j)
          endif

          accel_elastic(1,iglob) = accel_elastic(1,iglob) - weight* &
                (sigma_xx*nx + sigma_xz*nz)

          accel_elastic(2,iglob) = accel_elastic(2,iglob) - weight* &
                (sigma_xz*nx + sigma_zz*nz)

          if(isolver == 2) then
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - weight*( &
                b_sigma_xx*nx + b_sigma_xz*nz)

          b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - weight*( &
                b_sigma_xz*nx + b_sigma_zz*nz)
          endif !if(isolver == 2) then

        enddo

      enddo

    endif


! assembling accel_elastic for elastic elements
#ifdef USE_MPI
  if (nproc > 1 .and. any_elastic .and. ninterface_elastic > 0) then
    call assemble_MPI_vector_el(accel_elastic,npoin, &
      ninterface, ninterface_elastic,inum_interfaces_elastic, &
      max_interface_size, max_ibool_interfaces_size_el,&
      ibool_interfaces_elastic, nibool_interfaces_elastic, &
      tab_requests_send_recv_elastic,buffer_send_faces_vector_el, &
      buffer_recv_faces_vector_el, my_neighbours)
  endif

  if (nproc > 1 .and. any_elastic .and. ninterface_elastic > 0 .and. isolver == 2) then
    call assemble_MPI_vector_el(b_accel_elastic,npoin, &
      ninterface, ninterface_elastic,inum_interfaces_elastic, &
      max_interface_size, max_ibool_interfaces_size_el,&
      ibool_interfaces_elastic, nibool_interfaces_elastic, &
      tab_requests_send_recv_elastic,buffer_send_faces_vector_el, &
      buffer_recv_faces_vector_el, my_neighbours)
  endif
#endif


! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

  if(any_elastic) then

! --- add the source if it is a collocated force
    if(.not. initialfield) then

    do i_source=1,NSOURCE
! if this processor carries the source and the source element is elastic
      if (is_proc_source(i_source) == 1 .and. elastic(ispec_selected_source(i_source))) then

! collocated force
        if(source_type(i_source) == 1) then
       if(isolver == 1) then  ! forward wavefield
          accel_elastic(1,iglob_source(i_source)) = accel_elastic(1,iglob_source(i_source)) &
                            - sin(angleforce(i_source))*source_time_function(i_source,it)
          accel_elastic(2,iglob_source(i_source)) = accel_elastic(2,iglob_source(i_source)) &
                            + cos(angleforce(i_source))*source_time_function(i_source,it)
       else                   ! backward wavefield
      b_accel_elastic(1,iglob_source(i_source)) = b_accel_elastic(1,iglob_source(i_source)) &
                            - sin(angleforce(i_source))*source_time_function(i_source,NSTEP-it+1)
      b_accel_elastic(2,iglob_source(i_source)) = b_accel_elastic(2,iglob_source(i_source)) &
                            + cos(angleforce(i_source))*source_time_function(i_source,NSTEP-it+1)
       endif  !endif isolver == 1
        endif

      endif ! if this processor carries the source and the source element is elastic
    enddo ! do i_source=1,NSOURCE

    endif ! if not using an initial field

    accel_elastic(1,:) = accel_elastic(1,:) * rmass_inverse_elastic
    accel_elastic(2,:) = accel_elastic(2,:) * rmass_inverse_elastic

    veloc_elastic = veloc_elastic + deltatover2*accel_elastic

   if(isolver == 2) then
    b_accel_elastic(1,:) = b_accel_elastic(1,:) * rmass_inverse_elastic(:)
    b_accel_elastic(2,:) = b_accel_elastic(2,:) * rmass_inverse_elastic(:)

    b_veloc_elastic = b_veloc_elastic + b_deltatover2*b_accel_elastic
   endif

  endif

  if(any_elastic .and. isolver == 2) then ! kernels calculation
      do iglob = 1,npoin
            rho_k(iglob) =  accel_elastic(1,iglob)*b_displ_elastic(1,iglob) +&
                            accel_elastic(2,iglob)*b_displ_elastic(2,iglob)
      enddo
  endif

! ******************************************************************************************************************
! ************* main solver for the poroelastic elements: first the solid (u_s) than the fluid (w)
! ******************************************************************************************************************

  if(any_poroelastic) then

    if(isolver == 2) then
! if inviscid fluid, comment the reading and uncomment the zeroing
     read(23,rec=NSTEP-it+1) b_viscodampx
     read(24,rec=NSTEP-it+1) b_viscodampz
!     b_viscodampx(:) = ZERO
!     b_viscodampz(:) = ZERO
    endif

    call compute_forces_solid(npoin,nspec,myrank,nelemabs,numat,iglob_source, &
               ispec_selected_source,ispec_selected_rec,is_proc_source,which_proc_receiver,&
               source_type,it,NSTEP,anyabs, &
               initialfield,TURN_ATTENUATION_ON,TURN_VISCATTENUATION_ON,angleforce,deltatcube, &
               deltatfourth,twelvedeltat,fourdeltatsquare,ibool,kmato,numabs,poroelastic,codeabs, &
               accels_poroelastic,velocs_poroelastic,velocw_poroelastic,displs_poroelastic,displw_poroelastic,&
               b_accels_poroelastic,b_displs_poroelastic,b_velocw_poroelastic,b_displw_poroelastic,&
               density,porosity,tortuosity,permeability,poroelastcoef,xix,xiz,gammax,gammaz, &
               jacobian,source_time_function,sourcearray,adj_sourcearrays,e11, &
               e13,dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n, &
               dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll,inv_tau_sigma_nu2,&
               phi_nu2,Mu_nu2,N_SLS, &
               rx_viscous,rz_viscous,theta_e,theta_s,&
               b_viscodampx,b_viscodampz,&
               ibegin_bottom_poro,iend_bottom_poro,ibegin_top_poro,iend_top_poro, &
               jbegin_left_poro,jend_left_poro,jbegin_right_poro,jend_right_poro,&
               mufr_k,B_k,NSOURCE,nrec,isolver,save_forward,&
               b_absorb_poro_s_left,b_absorb_poro_s_right,b_absorb_poro_s_bottom,b_absorb_poro_s_top,&
               nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax,ib_xmin,ib_xmax,ib_zmin,ib_zmax)



    call compute_forces_fluid(npoin,nspec,myrank,nelemabs,numat,iglob_source, &
               ispec_selected_source,ispec_selected_rec,is_proc_source,which_proc_receiver,&
               source_type,it,NSTEP,anyabs, &
               initialfield,TURN_ATTENUATION_ON,TURN_VISCATTENUATION_ON,angleforce,deltatcube, &
               deltatfourth,twelvedeltat,fourdeltatsquare,ibool,kmato,numabs,poroelastic,codeabs, &
               accelw_poroelastic,velocw_poroelastic,displw_poroelastic,velocs_poroelastic,displs_poroelastic,&
               b_accelw_poroelastic,b_velocw_poroelastic,b_displw_poroelastic,b_displs_poroelastic,&
               density,porosity,tortuosity,permeability,poroelastcoef,xix,xiz,gammax,gammaz, &
               jacobian,source_time_function,sourcearray,adj_sourcearrays,e11, &
               e13,dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n, &
               dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll,inv_tau_sigma_nu2,&
               phi_nu2,Mu_nu2,N_SLS, &
               rx_viscous,rz_viscous,theta_e,theta_s,&
               b_viscodampx,b_viscodampz,&
               ibegin_bottom_poro,iend_bottom_poro,ibegin_top_poro,iend_top_poro, &
               jbegin_left_poro,jend_left_poro,jbegin_right_poro,jend_right_poro,&
               C_k,M_k,NSOURCE,nrec,isolver,save_forward,&
               b_absorb_poro_w_left,b_absorb_poro_w_right,b_absorb_poro_w_bottom,b_absorb_poro_w_top,&
               nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax,ib_xmin,ib_xmax,ib_zmin,ib_zmax)


    if(save_forward .and. isolver == 1) then
! if inviscid fluid, comment
     write(23,rec=it) b_viscodampx
     write(24,rec=it) b_viscodampz
    endif

    if(anyabs .and. save_forward .and. isolver == 1) then

!--- left absorbing boundary
      if(nspec_xmin >0) then
      do ispec = 1,nspec_xmin
       do id =1,2
         do i=1,NGLLZ
     write(45) b_absorb_poro_s_left(id,i,ispec,it)
     write(25) b_absorb_poro_w_left(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

!--- right absorbing boundary
      if(nspec_xmax >0) then
      do ispec = 1,nspec_xmax
       do id =1,2
         do i=1,NGLLZ
     write(46) b_absorb_poro_s_right(id,i,ispec,it)
     write(26) b_absorb_poro_w_right(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

!--- bottom absorbing boundary
      if(nspec_zmin >0) then
      do ispec = 1,nspec_zmin
       do id =1,2
         do i=1,NGLLX
     write(47) b_absorb_poro_s_bottom(id,i,ispec,it)
     write(29) b_absorb_poro_w_bottom(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

!--- top absorbing boundary
      if(nspec_zmax >0) then
      do ispec = 1,nspec_zmax
       do id =1,2
         do i=1,NGLLX
     write(48) b_absorb_poro_s_top(id,i,ispec,it)
     write(28) b_absorb_poro_w_top(id,i,ispec,it)
         enddo
       enddo
      enddo
      endif

    endif ! if(anyabs .and. save_forward .and. isolver == 1)

  endif !if(any_poroelastic) then

! *********************************************************
! ************* add coupling with the acoustic side
! *********************************************************

    if(coupled_acoustic_poroelastic) then

! loop on all the coupling edges
      do inum = 1,num_fluid_poro_edges

! get the edge of the acoustic element
        ispec_acoustic = fluid_poro_acoustic_ispec(inum)
        iedge_acoustic = fluid_poro_acoustic_iedge(inum)

! get the corresponding edge of the poroelastic element
        ispec_poroelastic = fluid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = fluid_poro_poroelastic_iedge(inum)

! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

! get point values for the acoustic side, which matches our side in the inverse direction
          i = ivalue_inverse(ipoin1D,iedge_acoustic)
          j = jvalue_inverse(ipoin1D,iedge_acoustic)
          iglob = ibool(i,j,ispec_acoustic)

! get poroelastic parameters
            phil = porosity(kmato(ispec_poroelastic))
            tortl = tortuosity(kmato(ispec_poroelastic))
            rhol_f = density(2,kmato(ispec_poroelastic))
            rhol_s = density(1,kmato(ispec_poroelastic))
            rhol_bar = (1._CUSTOM_REAL-phil)*rhol_s + phil*rhol_f

! compute pressure on the fluid/porous medium edge
          pressure = - potential_dot_dot_acoustic(iglob)
          if(isolver == 2) then
          b_pressure = - b_potential_dot_dot_acoustic(iglob)
          endif

! get point values for the poroelastic side
          ii2 = ivalue(ipoin1D,iedge_poroelastic)
          jj2 = jvalue(ipoin1D,iedge_poroelastic)
          iglob = ibool(ii2,jj2,ispec_poroelastic)

! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
! Blackwell Science, page 110, equation (4.60).
          if(iedge_acoustic == ITOP)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D
          weight = jacobian1D * wxgll(i)
          elseif(iedge_acoustic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
          weight = jacobian1D * wxgll(i)
          elseif(iedge_acoustic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
          weight = jacobian1D * wzgll(j)
          elseif(iedge_acoustic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
          weight = jacobian1D * wzgll(j)
          endif

! contribution to the solid phase
          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) + weight*nx*pressure*(1._CUSTOM_REAL-phil/tortl)
          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + weight*nz*pressure*(1._CUSTOM_REAL-phil/tortl)

! contribution to the fluid phase
          accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) + weight*nx*pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
          accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) + weight*nz*pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)

          if(isolver == 2) then
! contribution to the solid phase
          b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) + weight*nx*b_pressure*(1._CUSTOM_REAL-phil/tortl)
          b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + weight*nz*b_pressure*(1._CUSTOM_REAL-phil/tortl)

! contribution to the fluid phase
          b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) + weight*nx*b_pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
          b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) + weight*nz*b_pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
          endif !if(isolver == 2) then

        enddo ! do ipoin1D = 1,NGLLX

      enddo ! do inum = 1,num_fluid_poro_edges

    endif ! if(coupled_acoustic_poroelastic)

! ****************************************************************************
! ************* add coupling with the elastic side
! ****************************************************************************

    if(coupled_elastic_poroelastic) then

! loop on all the coupling edges
      do inum = 1,num_solid_poro_edges

! get the edge of the elastic element
        ispec_elastic = solid_poro_elastic_ispec(inum)
        iedge_elastic = solid_poro_elastic_iedge(inum)

! get the corresponding edge of the poroelastic element
        ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

! get point values for the elastic side, which matches our side in the inverse direction
          i = ivalue_inverse(ipoin1D,iedge_elastic)
          j = jvalue_inverse(ipoin1D,iedge_elastic)
          iglob = ibool(i,j,ispec_elastic)

! get poroelastic medium properties
    phil = porosity(kmato(ispec_poroelastic))
    tortl = tortuosity(kmato(ispec_poroelastic))
!
    rhol_s = density(1,kmato(ispec_poroelastic))
    rhol_f = density(2,kmato(ispec_poroelastic))
    rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f

! get elastic properties
    lambdal_relaxed = poroelastcoef(1,1,kmato(ispec_elastic))
    mul_relaxed = poroelastcoef(2,1,kmato(ispec_elastic))
    lambdalplus2mul_relaxed = poroelastcoef(3,1,kmato(ispec_elastic))

! derivative along x and along z for u_s and w
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          if(isolver == 2) then
          b_dux_dxi = ZERO
          b_duz_dxi = ZERO

          b_dux_dgamma = ZERO
          b_duz_dgamma = ZERO
          endif

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec_elastic))*hprime_zz(j,k)

            if(isolver == 2) then
            b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
            b_duz_dxi = b_duz_dxi + b_displ_elastic(2,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
            b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
            b_duz_dgamma = b_duz_dgamma + b_displ_elastic(2,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
            endif
          enddo

          xixl = xix(i,j,ispec_elastic)
          xizl = xiz(i,j,ispec_elastic)
          gammaxl = gammax(i,j,ispec_elastic)
          gammazl = gammaz(i,j,ispec_elastic)

! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          if(isolver == 2) then
          b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
          b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

          b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
          b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl
          endif
! compute stress tensor

! no attenuation
    sigma_xx = lambdalplus2mul_relaxed*dux_dxl + lambdal_relaxed*duz_dzl
    sigma_xz = mul_relaxed*(duz_dxl + dux_dzl)
    sigma_zz = lambdalplus2mul_relaxed*duz_dzl + lambdal_relaxed*dux_dxl

! full anisotropy
  if(TURN_ANISOTROPY_ON) then
! implement anisotropy in 2D
     sigma_xx = c11val*dux_dxl + c15val*(duz_dxl + dux_dzl) + c13val*duz_dzl
     sigma_zz = c13val*dux_dxl + c35val*(duz_dxl + dux_dzl) + c33val*duz_dzl
     sigma_xz = c15val*dux_dxl + c55val*(duz_dxl + dux_dzl) + c35val*duz_dzl
  endif

    if(isolver == 2) then
    b_sigma_xx = lambdalplus2mul_relaxed*b_dux_dxl + lambdal_relaxed*b_duz_dzl
    b_sigma_xz = mul_relaxed*(b_duz_dxl + b_dux_dzl)
    b_sigma_zz = lambdalplus2mul_relaxed*b_duz_dzl + lambdal_relaxed*b_dux_dxl
    endif ! if(isolver == 2)

! get point values for the poroelastic side
          i = ivalue(ipoin1D,iedge_poroelastic)
          j = jvalue(ipoin1D,iedge_poroelastic)
          iglob = ibool(i,j,ispec_poroelastic)

! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
! Blackwell Science, page 110, equation (4.60).
          if(iedge_poroelastic == ITOP)then
            xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D
          weight = jacobian1D * wxgll(i)
          elseif(iedge_poroelastic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
          weight = jacobian1D * wxgll(i)
          elseif(iedge_poroelastic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
          weight = jacobian1D * wzgll(j)
          elseif(iedge_poroelastic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
          weight = jacobian1D * wzgll(j)
          endif

! contribution to the solid phase
          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) + &
                weight*(sigma_xx*nx + sigma_xz*nz)*(1._CUSTOM_REAL - phil/tortl )

          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + &
                weight*(sigma_xz*nx + sigma_zz*nz)*(1._CUSTOM_REAL - phil/tortl )

! contribution to the fluid phase
          accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob)  - &
                weight*(rhol_f/rhol_bar - 1._CUSTOM_REAL)*(sigma_xx*nx+sigma_xz*nz)

          accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - &
                weight*(rhol_f/rhol_bar - 1._CUSTOM_REAL)*(sigma_xz*nx+sigma_zz*nz)

          if(isolver == 2) then
! contribution to the solid phase
          b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) + &
                weight*(b_sigma_xx*nx + b_sigma_xz*nz)*(1._CUSTOM_REAL - phil/tortl)

          b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + &
                weight*(b_sigma_xz*nx + b_sigma_zz*nz)*(1._CUSTOM_REAL - phil/tortl)

! contribution to the fluid phase
          b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob)  - &
                weight*(rhol_f/rhol_bar - 1._CUSTOM_REAL)*(b_sigma_xx*nx + b_sigma_xz*nz)

          b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) - &
                weight*(rhol_f/rhol_bar - 1._CUSTOM_REAL)*(b_sigma_xz*nx + b_sigma_zz*nz)
          endif !if(isolver == 2) then

        enddo

      enddo

    endif ! if(coupled_elastic_poroelastic)


! assembling accels_proelastic & accelw_poroelastic for poroelastic elements
#ifdef USE_MPI
  if (nproc > 1 .and. any_poroelastic .and. ninterface_poroelastic > 0) then
    call assemble_MPI_vector_po(accels_poroelastic,accelw_poroelastic,npoin, &
      ninterface, ninterface_poroelastic,inum_interfaces_poroelastic, &
      max_interface_size, max_ibool_interfaces_size_po,&
      ibool_interfaces_poroelastic, nibool_interfaces_poroelastic, &
      tab_requests_send_recv_poroelastic,buffer_send_faces_vector_pos,buffer_send_faces_vector_pow, &
      buffer_recv_faces_vector_pos,buffer_recv_faces_vector_pow, &
      my_neighbours)
  endif

  if (nproc > 1 .and. any_poroelastic .and. ninterface_poroelastic > 0 .and. isolver == 2) then
    call assemble_MPI_vector_po(b_accels_poroelastic,b_accelw_poroelastic,npoin, &
      ninterface, ninterface_poroelastic,inum_interfaces_poroelastic, &
      max_interface_size, max_ibool_interfaces_size_po,&
      ibool_interfaces_poroelastic, nibool_interfaces_poroelastic, &
      tab_requests_send_recv_poroelastic,buffer_send_faces_vector_pos,buffer_send_faces_vector_pow, &
      buffer_recv_faces_vector_pos,buffer_recv_faces_vector_pow, &
      my_neighbours)
   endif
#endif


! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

 if(any_poroelastic) then


! --- add the source if it is a collocated force
    if(.not. initialfield) then

    do i_source=1,NSOURCE
! if this processor carries the source and the source element is elastic
      if (is_proc_source(i_source) == 1 .and. poroelastic(ispec_selected_source(i_source))) then

    phil = porosity(kmato(ispec_selected_source(i_source)))
    tortl = tortuosity(kmato(ispec_selected_source(i_source)))
    rhol_s = density(1,kmato(ispec_selected_source(i_source)))
    rhol_f = density(2,kmato(ispec_selected_source(i_source)))
    rhol_bar = (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f

! collocated force
        if(source_type(i_source) == 1) then
       if(isolver == 1) then  ! forward wavefield
! s
      accels_poroelastic(1,iglob_source(i_source)) = accels_poroelastic(1,iglob_source(i_source)) - &
                               (1._CUSTOM_REAL - phil/tortl)*sin(angleforce(i_source))*source_time_function(i_source,it)
      accels_poroelastic(2,iglob_source(i_source)) = accels_poroelastic(2,iglob_source(i_source)) + &
                               (1._CUSTOM_REAL - phil/tortl)*cos(angleforce(i_source))*source_time_function(i_source,it)
! w
      accelw_poroelastic(1,iglob_source(i_source)) = accelw_poroelastic(1,iglob_source(i_source)) - &
         (1._CUSTOM_REAL - rhol_f/rhol_bar)*sin(angleforce(i_source))*source_time_function(i_source,it)
      accelw_poroelastic(2,iglob_source(i_source)) = accelw_poroelastic(2,iglob_source(i_source)) + &
         (1._CUSTOM_REAL - rhol_f/rhol_bar)*cos(angleforce(i_source))*source_time_function(i_source,it)

       else                   ! backward wavefield
! b_s
      b_accels_poroelastic(1,iglob_source(i_source)) = b_accels_poroelastic(1,iglob_source(i_source)) - &
                               (1._CUSTOM_REAL - phil/tortl)*sin(angleforce(i_source))*source_time_function(i_source,NSTEP-it+1)
      b_accels_poroelastic(2,iglob_source(i_source)) = b_accels_poroelastic(2,iglob_source(i_source)) + &
                               (1._CUSTOM_REAL - phil/tortl)*cos(angleforce(i_source))*source_time_function(i_source,NSTEP-it+1)
!b_w
      b_accelw_poroelastic(1,iglob_source(i_source)) = b_accelw_poroelastic(1,iglob_source(i_source)) - &
         (1._CUSTOM_REAL - rhol_f/rhol_bar)*sin(angleforce(i_source))*source_time_function(i_source,NSTEP-it+1)
      b_accelw_poroelastic(2,iglob_source(i_source)) = b_accelw_poroelastic(2,iglob_source(i_source)) + &
         (1._CUSTOM_REAL - rhol_f/rhol_bar)*cos(angleforce(i_source))*source_time_function(i_source,NSTEP-it+1)
       endif !endif isolver == 1
        endif

      endif ! if this processor carries the source and the source element is elastic
    enddo ! do i_source=1,NSOURCE

    endif ! if not using an initial field

    accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
    accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
    velocs_poroelastic = velocs_poroelastic + deltatover2*accels_poroelastic

    accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
    accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)
    velocw_poroelastic = velocw_poroelastic + deltatover2*accelw_poroelastic

   if(isolver == 2) then
    b_accels_poroelastic(1,:) = b_accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
    b_accels_poroelastic(2,:) = b_accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
    b_velocs_poroelastic = b_velocs_poroelastic + b_deltatover2*b_accels_poroelastic

    b_accelw_poroelastic(1,:) = b_accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
    b_accelw_poroelastic(2,:) = b_accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)
    b_velocw_poroelastic = b_velocw_poroelastic + b_deltatover2*b_accelw_poroelastic
   endif

  endif

  if(any_poroelastic .and. isolver ==2) then
   do iglob =1,npoin
            rhot_k(iglob) = accels_poroelastic(1,iglob) * b_displs_poroelastic(1,iglob) + &
                  accels_poroelastic(2,iglob) * b_displs_poroelastic(2,iglob)
            rhof_k(iglob) = accelw_poroelastic(1,iglob) * b_displs_poroelastic(1,iglob) + &
                  accelw_poroelastic(2,iglob) * b_displs_poroelastic(2,iglob) + &
                  accels_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                  accels_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)
            sm_k(iglob) =  accelw_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                  accelw_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)
            eta_k(iglob) = velocw_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                  velocw_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)
   enddo
  endif

!----  compute kinetic and potential energy
  if(OUTPUT_ENERGY) &
     call compute_energy(displ_elastic,veloc_elastic,displs_poroelastic,velocs_poroelastic, &
         displw_poroelastic,velocw_poroelastic, &
         xix,xiz,gammax,gammaz,jacobian,ibool,elastic,poroelastic,hprime_xx,hprime_zz, &
         nspec,npoin,assign_external_model,it,deltat,t0(1),kmato,poroelastcoef,density, &
         porosity,tortuosity, &
         vpext,vsext,rhoext,wxgll,wzgll,numat, &
         pressure_element,vector_field_element,e1,e11, &
         potential_dot_acoustic,potential_dot_dot_acoustic,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,Mu_nu1,Mu_nu2,N_SLS)

!----  display time step and max of norm of displacement
  if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then

    if (myrank == 0) then
      write(IOUT,*)
      if(time >= 1.d-3 .and. time < 1000.d0) then
        write(IOUT,"('Time step number ',i7,'   t = ',f9.4,' s out of ',i7)") it,time,NSTEP
      else
        write(IOUT,"('Time step number ',i7,'   t = ',1pe12.6,' s out of ',i7)") it,time,NSTEP
      endif
      write(IOUT,*) 'We have done ',sngl(100.d0*dble(it-1)/dble(NSTEP-1)),'% of the total'
    endif

    if(any_elastic_glob) then
      if(any_elastic) then
        displnorm_all = maxval(sqrt(displ_elastic(1,:)**2 + displ_elastic(2,:)**2))
      else
        displnorm_all = 0.d0
      endif
      displnorm_all_glob = displnorm_all
#ifdef USE_MPI
      call MPI_ALLREDUCE (displnorm_all, displnorm_all_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
#endif
      if (myrank == 0) write(IOUT,*) 'Max norm of vector field in solid (elastic) = ',displnorm_all_glob
! check stability of the code in solid, exit if unstable
! negative values can occur with some compilers when the unstable value is greater
! than the greatest possible floating-point number of the machine
      if(displnorm_all_glob > STABILITY_THRESHOLD .or. displnorm_all_glob < 0) &
        call exit_MPI('code became unstable and blew up in solid')
    endif

    if(any_poroelastic_glob) then
      if(any_poroelastic) then
        displnorm_all = maxval(sqrt(displs_poroelastic(1,:)**2 + displs_poroelastic(2,:)**2))
      else
        displnorm_all = 0.d0
      endif
      displnorm_all_glob = displnorm_all
#ifdef USE_MPI
      call MPI_ALLREDUCE (displnorm_all, displnorm_all_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
#endif
      if (myrank == 0) write(IOUT,*) 'Max norm of vector field in solid (poroelastic) = ',displnorm_all_glob
! check stability of the code in solid, exit if unstable
! negative values can occur with some compilers when the unstable value is greater
! than the greatest possible floating-point number of the machine
      if(displnorm_all_glob > STABILITY_THRESHOLD .or. displnorm_all_glob < 0) &
        call exit_MPI('code became unstable and blew up in solid (poroelastic)')

      if(any_poroelastic) then
        displnorm_all = maxval(sqrt(displw_poroelastic(1,:)**2 + displw_poroelastic(2,:)**2))
      else
        displnorm_all = 0.d0
      endif
      displnorm_all_glob = displnorm_all
#ifdef USE_MPI
      call MPI_ALLREDUCE (displnorm_all, displnorm_all_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
#endif
      if (myrank == 0) write(IOUT,*) 'Max norm of vector field in fluid (poroelastic) = ',displnorm_all_glob
! check stability of the code in solid, exit if unstable
! negative values can occur with some compilers when the unstable value is greater
! than the greatest possible floating-point number of the machine
      if(displnorm_all_glob > STABILITY_THRESHOLD .or. displnorm_all_glob < 0) &
        call exit_MPI('code became unstable and blew up in fluid (poroelastic)')
    endif

    if(any_acoustic_glob) then
      if(any_acoustic) then
        displnorm_all = maxval(abs(potential_acoustic(:)))
      else
        displnorm_all = 0.d0
      endif
      displnorm_all_glob = displnorm_all
#ifdef USE_MPI
      call MPI_ALLREDUCE (displnorm_all, displnorm_all_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
#endif
      if (myrank == 0) write(IOUT,*) 'Max absolute value of scalar field in fluid (acoustic) = ',displnorm_all_glob
! check stability of the code in fluid, exit if unstable
! negative values can occur with some compilers when the unstable value is greater
! than the greatest possible floating-point number of the machine
      if(displnorm_all_glob > STABILITY_THRESHOLD .or. displnorm_all_glob < 0) &
        call exit_MPI('code became unstable and blew up in fluid')
    endif
    if (myrank == 0) write(IOUT,*)
  endif

! loop on all the receivers to compute and store the seismograms
  do irecloc = 1,nrecloc

    irec = recloc(irecloc)

    ispec = ispec_selected_rec(irec)

! compute pressure in this element if needed
    if(seismotype == 4) then

       call compute_pressure_one_element(pressure_element,potential_dot_dot_acoustic,displ_elastic,&
            displs_poroelastic,displw_poroelastic,elastic,poroelastic,&
            xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,assign_external_model, &
            numat,kmato,density,porosity,tortuosity,poroelastcoef,vpext,vsext,rhoext,ispec,e1,e11, &
            TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,Mu_nu1,Mu_nu2,N_SLS)

    else if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then

! for acoustic medium, compute vector field from gradient of potential for seismograms
       if(seismotype == 1) then
          call compute_vector_one_element(vector_field_element,potential_acoustic,displ_elastic,displs_poroelastic,&
               elastic,poroelastic,xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,ispec,numat,kmato,&
               density,rhoext,assign_external_model)
       else if(seismotype == 2) then
          call compute_vector_one_element(vector_field_element,potential_dot_acoustic,veloc_elastic,velocs_poroelastic, &
               elastic,poroelastic,xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,ispec,numat,kmato,&
               density,rhoext,assign_external_model)
       else if(seismotype == 3) then
          call compute_vector_one_element(vector_field_element,potential_dot_dot_acoustic,accel_elastic,accels_poroelastic, &
               elastic,poroelastic,xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,ispec,numat,kmato,&
               density,rhoext,assign_external_model)
       endif

    else if(seismotype == 5) then
       call compute_curl_one_element(curl_element,displ_elastic,displs_poroelastic,elastic,poroelastic, &
            xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin, ispec)
    endif

! perform the general interpolation using Lagrange polynomials
    valux = ZERO
    valuz = ZERO
    valcurl = ZERO

    do j = 1,NGLLZ
      do i = 1,NGLLX

        iglob = ibool(i,j,ispec)

        hlagrange = hxir_store(irec,i)*hgammar_store(irec,j)

        dcurld=ZERO

        if(seismotype == 4) then

          dxd = pressure_element(i,j)
          dzd = ZERO

        else if(.not. elastic(ispec) .and. .not. poroelastic(ispec) .and.  seismotype /= 6) then

          dxd = vector_field_element(1,i,j)
          dzd = vector_field_element(2,i,j)

        else if(seismotype == 6) then

          dxd = potential_acoustic(iglob)
          dzd = ZERO

        else if(seismotype == 1) then

             if(poroelastic(ispec)) then
          dxd = displs_poroelastic(1,iglob)
          dzd = displs_poroelastic(2,iglob)
             elseif(elastic(ispec)) then
          dxd = displ_elastic(1,iglob)
          dzd = displ_elastic(2,iglob)
             endif

        else if(seismotype == 2) then

             if(poroelastic(ispec)) then
          dxd = velocs_poroelastic(1,iglob)
          dzd = velocs_poroelastic(2,iglob)
             elseif(elastic(ispec)) then
          dxd = veloc_elastic(1,iglob)
          dzd = veloc_elastic(2,iglob)
             endif

        else if(seismotype == 3) then

             if(poroelastic(ispec)) then
          dxd = accels_poroelastic(1,iglob)
          dzd = accels_poroelastic(2,iglob)
             elseif(elastic(ispec)) then
          dxd = accel_elastic(1,iglob)
          dzd = accel_elastic(2,iglob)
             endif

        else if(seismotype == 5) then

             if(poroelastic(ispec)) then
          dxd = displs_poroelastic(1,iglob)
          dzd = displs_poroelastic(2,iglob)
             elseif(elastic(ispec)) then
          dxd = displ_elastic(1,iglob)
          dzd = displ_elastic(2,iglob)
             endif
          dcurld = curl_element(i,j)

        endif

! compute interpolated field
        valux = valux + dxd*hlagrange
        valuz = valuz + dzd*hlagrange
        valcurl = valcurl + dcurld*hlagrange

      enddo
    enddo

! rotate seismogram components if needed, except if recording pressure, which is a scalar
    if(seismotype /= 4 .and. seismotype /= 6) then
      sisux(seismo_current,irecloc) =   cosrot_irec(irecloc)*valux + sinrot_irec(irecloc)*valuz
      sisuz(seismo_current,irecloc) = - sinrot_irec(irecloc)*valux + cosrot_irec(irecloc)*valuz
    else
      sisux(seismo_current,irecloc) = valux
      sisuz(seismo_current,irecloc) = ZERO
    endif
    siscurl(seismo_current,irecloc) = valcurl

 enddo


!
!----- ecriture des kernels
!
! kernels output
  if(isolver == 2) then

   if(any_acoustic) then

    do ispec = 1, nspec
     if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then
      do k = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,k,ispec)
    kappal_ac_global(iglob) = poroelastcoef(1,2,kmato(ispec))
    rhol_ac_global(iglob) = density(2,kmato(ispec))
          enddo
      enddo
     endif
    enddo

          do iglob =1,npoin
            rho_ac_kl(iglob) = rho_ac_kl(iglob) - rhol_ac_global(iglob)  * rho_ac_k(iglob) * deltat
            kappa_ac_kl(iglob) = kappa_ac_kl(iglob) - kappal_ac_global(iglob) * kappa_ac_k(iglob) * deltat
!
            rhop_ac_kl(iglob) = rho_ac_kl(iglob) + kappa_ac_kl(iglob)
            alpha_ac_kl(iglob) = TWO *  kappa_ac_kl(iglob)
          enddo

    endif !if(any_acoustic)

   if(any_elastic) then

    do ispec = 1, nspec
     if(elastic(ispec)) then
      do k = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,k,ispec)
    mul_global(iglob) = poroelastcoef(2,1,kmato(ispec))
    kappal_global(iglob) = poroelastcoef(3,1,kmato(ispec)) - 4._CUSTOM_REAL*mul_global(iglob)/3._CUSTOM_REAL
    rhol_global(iglob) = density(1,kmato(ispec))
          enddo
      enddo
     endif
    enddo

          do iglob =1,npoin
            rho_kl(iglob) = rho_kl(iglob) - rhol_global(iglob)  * rho_k(iglob) * deltat
            mu_kl(iglob) = mu_kl(iglob) - TWO * mul_global(iglob) * mu_k(iglob) * deltat
            kappa_kl(iglob) = kappa_kl(iglob) - kappal_global(iglob) * kappa_k(iglob) * deltat
!
            rhop_kl(iglob) = rho_kl(iglob) + kappa_kl(iglob) + mu_kl(iglob)
            beta_kl(iglob) = TWO * (mu_kl(iglob) - 4._CUSTOM_REAL * mul_global(iglob) &
                  / (3._CUSTOM_REAL * kappal_global(iglob)) * kappa_kl(iglob))
            alpha_kl(iglob) = TWO * (1._CUSTOM_REAL + 4._CUSTOM_REAL * mul_global(iglob)/&
                   (3._CUSTOM_REAL * kappal_global(iglob))) * kappa_kl(iglob)
          enddo

   endif !if(any_elastic)

  if(any_poroelastic) then

    do ispec = 1, nspec
     if(poroelastic(ispec)) then
      do k = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,k,ispec)
    phil_global(iglob) = porosity(kmato(ispec))
    tortl_global(iglob) = tortuosity(kmato(ispec))
    rhol_s_global(iglob) = density(1,kmato(ispec))
    rhol_f_global(iglob) = density(2,kmato(ispec))
    rhol_bar_global(iglob) =  (1._CUSTOM_REAL - phil_global(iglob))*rhol_s_global(iglob) &
              + phil_global(iglob)*rhol_f_global(iglob)
    etal_f_global(iglob) = poroelastcoef(2,2,kmato(ispec))
    permlxx_global(iglob) = permeability(1,kmato(ispec))
    permlxz_global(iglob) = permeability(2,kmato(ispec))
    permlzz_global(iglob) = permeability(3,kmato(ispec))
    mulfr_global(iglob) = poroelastcoef(2,3,kmato(ispec))
          enddo
       enddo
     endif
    enddo

      do iglob =1,npoin
            rhot_kl(iglob) = rhot_kl(iglob) - deltat * rhol_bar_global(iglob) * rhot_k(iglob)
            rhof_kl(iglob) = rhof_kl(iglob) - deltat * rhol_f_global(iglob) * rhof_k(iglob)
            sm_kl(iglob) = sm_kl(iglob) - deltat * rhol_f_global(iglob)*tortl_global(iglob)/phil_global(iglob) * sm_k(iglob)
!at the moment works with constant permeability
            eta_kl(iglob) = eta_kl(iglob) - deltat * etal_f_global(iglob)/permlxx_global(iglob) * eta_k(iglob)
            B_kl(iglob) = B_kl(iglob) - deltat * B_k(iglob)
            C_kl(iglob) = C_kl(iglob) - deltat * C_k(iglob)
            M_kl(iglob) = M_kl(iglob) - deltat * M_k(iglob)
            mufr_kl(iglob) = mufr_kl(iglob) - TWO * deltat * mufr_k(iglob)
! density kernels
            rholb = rhol_bar_global(iglob) - phil_global(iglob)*rhol_f_global(iglob)/tortl_global(iglob)
            rhob_kl(iglob) = rhot_kl(iglob) + B_kl(iglob) + mufr_kl(iglob)
            rhofb_kl(iglob) = rhof_kl(iglob) + C_kl(iglob) + M_kl(iglob) + sm_kl(iglob)
            Bb_kl(iglob) = B_kl(iglob)
            Cb_kl(iglob) = C_kl(iglob)
            Mb_kl(iglob) = M_kl(iglob)
            mufrb_kl(iglob) = mufr_kl(iglob)
            phi_kl(iglob) = - sm_kl(iglob) - M_kl(iglob)
! wave speed kernels
            dd1 = (1._CUSTOM_REAL+rholb/rhol_f_global(iglob))*ratio**2 + 2._CUSTOM_REAL*ratio +&
                tortl_global(iglob)/phil_global(iglob)
            rhobb_kl(iglob) = rhob_kl(iglob) - &
                phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*B_biot) * &
                   (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1+&
                   (rhol_bar_global(iglob)**2*ratio**2/rhol_f_global(iglob)**2*(phil_global(iglob)/tortl_global(iglob)*&
                   ratio+1)*(phil_global(iglob)/tortl_global(iglob)*ratio+phil_global(iglob)/tortl_global(iglob)*&
                   (1+rhol_f_global(iglob)/rhol_bar_global(iglob))-1) )/dd1**2 )- FOUR_THIRDS*cssquare )*&
                   Bb_kl(iglob) - &
                rhol_bar_global(iglob)*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                   (phil_global(iglob)/tortl_global(iglob)*ratio + 1._CUSTOM_REAL)**2/dd1**2*Mb_kl(iglob) + &
                rhol_bar_global(iglob)*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                   phil_global(iglob)*ratio/tortl_global(iglob)*(phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                   (1+rhol_bar_global(iglob)*ratio/rhol_f_global(iglob))/dd1**2)*Cb_kl(iglob)+ &
                phil_global(iglob)*rhol_f_global(iglob)*cssquare/(tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(iglob)
           rhofbb_kl(iglob) = rhofb_kl(iglob) + &
                phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*B_biot) * &
                   (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1+&
                   (rhol_bar_global(iglob)**2*ratio**2/rhol_f_global(iglob)**2*(phil_global(iglob)/tortl_global(iglob)*&
                   ratio+1)*(phil_global(iglob)/tortl_global(iglob)*ratio+phil_global(iglob)/tortl_global(iglob)*&
                   (1+rhol_f_global(iglob)/rhol_bar_global(iglob))-1) )/dd1**2 )- FOUR_THIRDS*cssquare )*&
                   Bb_kl(iglob) + &
                rhol_bar_global(iglob)*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                   (phil_global(iglob)/tortl_global(iglob)*ratio + 1._CUSTOM_REAL)**2/dd1**2*Mb_kl(iglob) - &
                rhol_bar_global(iglob)*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                   phil_global(iglob)*ratio/tortl_global(iglob)*(phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                   (1+rhol_bar_global(iglob)*ratio/rhol_f_global(iglob))/dd1**2)*Cb_kl(iglob)- &
                phil_global(iglob)*rhol_f_global(iglob)*cssquare/(tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(iglob)
           phib_kl(iglob) = phi_kl(iglob) - &
                phil_global(iglob)*rhol_bar_global(iglob)/(tortl_global(iglob)*B_biot) * ( cpIsquare - rhol_f_global(iglob)/&
                   rhol_bar_global(iglob)*cpIIsquare- &
                   (cpIsquare-cpIIsquare)*( (TWO*ratio**2*phil_global(iglob)/tortl_global(iglob) + (1._CUSTOM_REAL+&
                   rhol_f_global(iglob)/rhol_bar_global(iglob))*(TWO*ratio*phil_global(iglob)/tortl_global(iglob)+&
                   1._CUSTOM_REAL))/dd1 + (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(phil_global(iglob)*&
                   ratio/tortl_global(iglob)+phil_global(iglob)/tortl_global(iglob)*(1._CUSTOM_REAL+rhol_f_global(iglob)/&
                   rhol_bar_global(iglob))-1._CUSTOM_REAL)*((1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-&
                   TWO*phil_global(iglob)/tortl_global(iglob))*ratio**2+TWO*ratio)/dd1**2 ) - &
                   FOUR_THIRDS*rhol_f_global(iglob)*cssquare/rhol_bar_global(iglob) )*Bb_kl(iglob) + &
                rhol_f_global(iglob)/M_biot * (cpIsquare-cpIIsquare)*(&
                   TWO*ratio*(phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2*((1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                   rhol_f_global(iglob)-TWO*phil_global(iglob)/tortl_global(iglob))*ratio**2+TWO*ratio)/dd1**2&
                   )*Mb_kl(iglob) + &
                phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*C_biot)*(cpIsquare-cpIIsquare)*ratio* (&
                   (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob)*ratio)/dd1 - &
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                   rhol_f_global(iglob)*ratio)*((1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-TWO*&
                   phil_global(iglob)/tortl_global(iglob))*ratio+TWO)/dd1**2&
                    )*Cb_kl(iglob) -&
                phil_global(iglob)*rhol_f_global(iglob)*cssquare/(tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(iglob)
           cpI_kl(iglob) = 2._CUSTOM_REAL*cpIsquare/B_biot*rhol_bar_global(iglob)*( &
                   1._CUSTOM_REAL-phil_global(iglob)/tortl_global(iglob) + &
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*&
                   ratio+phil_global(iglob)/tortl_global(iglob)*(1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-&
                   1._CUSTOM_REAL)/dd1 &
                    )* Bb_kl(iglob) +&
                2._CUSTOM_REAL*cpIsquare*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot) *&
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2/dd1*Mb_kl(iglob)+&
                2._CUSTOM_REAL*cpIsquare*rhol_f_global(iglob)/C_biot * &
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                   rhol_f_global(iglob)*ratio)/dd1*Cb_kl(iglob)
           cpII_kl(iglob) = 2._CUSTOM_REAL*cpIIsquare*rhol_bar_global(iglob)/B_biot * (&
                   phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*rhol_bar_global(iglob)) - &
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*&
                   ratio+phil_global(iglob)/tortl_global(iglob)*(1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-&
                   1._CUSTOM_REAL)/dd1  ) * Bb_kl(iglob) +&
                2._CUSTOM_REAL*cpIIsquare*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot) * (&
                   1._CUSTOM_REAL - (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2/dd1  )*Mb_kl(iglob) + &
                2._CUSTOM_REAL*cpIIsquare*rhol_f_global(iglob)/C_biot * (&
                   1._CUSTOM_REAL - (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(1._CUSTOM_REAL+&
                   rhol_bar_global(iglob)/rhol_f_global(iglob)*ratio)/dd1  )*Cb_kl(iglob)
           cs_kl(iglob) = - 8._CUSTOM_REAL/3._CUSTOM_REAL*cssquare*rhol_bar_global(iglob)/B_biot*(1._CUSTOM_REAL-&
                   phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*rhol_bar_global(iglob)))*Bb_kl(iglob) + &
                2._CUSTOM_REAL*(rhol_bar_global(iglob)-rhol_f_global(iglob)*phil_global(iglob)/tortl_global(iglob))/&
                   mulfr_global(iglob)*cssquare*mufrb_kl(iglob)
           ratio_kl(iglob) = ratio*rhol_bar_global(iglob)*phil_global(iglob)/(tortl_global(iglob)*B_biot) * &
                   (cpIsquare-cpIIsquare) * ( &
                   phil_global(iglob)/tortl_global(iglob)*(2._CUSTOM_REAL*ratio+1._CUSTOM_REAL+rhol_f_global(iglob)/ &
                   rhol_bar_global(iglob))/dd1 - (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                   (phil_global(iglob)/tortl_global(iglob)*ratio+phil_global(iglob)/tortl_global(iglob)*(&
                   1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-1._CUSTOM_REAL)*(2._CUSTOM_REAL*ratio*(&
                   1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob)) +&
                   2._CUSTOM_REAL)/dd1**2  )*Bb_kl(iglob) + &
                ratio*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot)*(cpIsquare-cpIIsquare) * &
                   2._CUSTOM_REAL*phil_global(iglob)/tortl_global(iglob) * (&
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2*((1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                   rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob))*ratio+1._CUSTOM_REAL)/dd1**2 )*Mb_kl(iglob) +&
                ratio*rhol_f_global(iglob)/C_biot*(cpIsquare-cpIIsquare) * (&
                   (2._CUSTOM_REAL*phil_global(iglob)*rhol_bar_global(iglob)*ratio/(tortl_global(iglob)*rhol_f_global(iglob))+&
                   phil_global(iglob)/tortl_global(iglob)+rhol_bar_global(iglob)/rhol_f_global(iglob))/dd1 - &
                   2._CUSTOM_REAL*phil_global(iglob)/tortl_global(iglob)*(phil_global(iglob)/tortl_global(iglob)*ratio+&
                   1._CUSTOM_REAL)*(1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)*ratio)*((1._CUSTOM_REAL+&
                   rhol_bar_global(iglob)/rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob))*ratio+1._CUSTOM_REAL)/&
                   dd1**2 )*Cb_kl(iglob)
      enddo

   endif ! if(any_poroelastic)

   endif ! if(isolver == 2)

!
!----  display results at given time steps
!
  if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then

!
! kernels output files
!

   if(isolver == 2 .and. it == NSTEP) then

  if ( myrank == 0 ) then
  write(IOUT,*) 'Writing Kernels file'
  endif

    if(any_acoustic) then
        write(filename,'(a,i7.7)') 'OUTPUT_FILES/snapshot_rho_kappa_',it
        write(filename2,'(a,i7.7)') 'OUTPUT_FILES/snapshot_rhop_c_',it

        open(unit = 97, file = trim(filename),status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing snapshot to disk'
        open(unit = 98, file = trim(filename2),status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing snapshot to disk'

     do iglob =1,npoin
        xx = coord(1,iglob)/maxval(coord(1,:))
        zz = coord(2,iglob)/maxval(coord(1,:))
         write(97,'(5e12.3)')xx,zz,rho_ac_kl(iglob),kappa_ac_kl(iglob)
         write(98,'(5e12.3)')xx,zz,rhop_ac_kl(iglob),alpha_ac_kl(iglob)
     enddo
    close(97)
    close(98)
    endif

    if(any_elastic) then
        write(filename,'(a,i7.7)') 'OUTPUT_FILES/snapshot_rho_kappa_mu_',it
        write(filename2,'(a,i7.7)') 'OUTPUT_FILES/snapshot_rhop_alpha_beta_',it

        open(unit = 97, file = trim(filename),status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing snapshot to disk'
        open(unit = 98, file = trim(filename2),status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing snapshot to disk'

     do iglob =1,npoin
        xx = coord(1,iglob)/maxval(coord(1,:))
        zz = coord(2,iglob)/maxval(coord(1,:))
         write(97,'(5e12.3)')xx,zz,rho_kl(iglob),kappa_kl(iglob),mu_kl(iglob)
         write(98,'(5e12.3)')xx,zz,rhop_kl(iglob),alpha_kl(iglob),beta_kl(iglob)
     enddo
    close(97)
    close(98)
    endif

    if(any_poroelastic) then
! Primary kernels
        write(filename,'(a,i7.7)') 'OUTPUT_FILES/snapshot_mu_B_C_',it
        write(filename2,'(a,i7.7)') 'OUTPUT_FILES/snapshot_M_rho_rhof_',it
        write(filename3,'(a,i7.7)') 'OUTPUT_FILES/snapshot_m_eta_',it

        open(unit = 97, file = trim(filename),status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing snapshot to disk'

        open(unit = 98, file = trim(filename2),status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing snapshot to disk'

        open(unit = 99, file = trim(filename3),status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing snapshot to disk'
! Wavespeed kernels
!        write(filename,'(a,i7.7)') 'OUTPUT_FILES/snapshot_cpI_cpII_cs_',it
!        write(filename2,'(a,i7.7)') 'OUTPUT_FILES/snapshot_rhobb_rhofbb_ratio_',it
!        write(filename3,'(a,i7.7)') 'OUTPUT_FILES/snapshot_phib_eta_',it
! Density normalized kernels
        write(filename,'(a,i7.7)') 'OUTPUT_FILES/snapshot_mub_Bb_Cb_',it
        write(filename2,'(a,i7.7)') 'OUTPUT_FILES/snapshot_Mb_rhob_rhofb_',it
        write(filename3,'(a,i7.7)') 'OUTPUT_FILES/snapshot_mb_etab_',it

        open(unit = 17, file = trim(filename),status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing snapshot to disk'

        open(unit = 18, file = trim(filename2),status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing snapshot to disk'

        open(unit = 19, file = trim(filename3),status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing snapshot to disk'

     do iglob =1,npoin
        xx = coord(1,iglob)/maxval(coord(1,:))
        zz = coord(2,iglob)/maxval(coord(1,:))
         write(97,'(5e12.3)')xx,zz,mufr_kl(iglob),B_kl(iglob),C_kl(iglob)
         write(98,'(5e12.3)')xx,zz,M_kl(iglob),rhot_kl(iglob),rhof_kl(iglob)
         write(99,'(5e12.3)')xx,zz,sm_kl(iglob),eta_kl(iglob)
         write(17,'(5e12.3)')xx,zz,mufrb_kl(iglob),Bb_kl(iglob),Cb_kl(iglob)
         write(18,'(5e12.3)')xx,zz,Mb_kl(iglob),rhob_kl(iglob),rhofb_kl(iglob)
         write(19,'(5e12.3)')xx,zz,phi_kl(iglob),eta_kl(iglob)
!         write(17,'(5e12.3)')xx,zz,cpI_kl(iglob),cpII_kl(iglob),cs_kl(iglob)
!         write(18,'(5e12.3)')xx,zz,rhobb_kl(iglob),rhofbb_kl(iglob),ratio_kl(iglob)
!         write(19,'(5e12.3)')xx,zz,phib_kl(iglob),eta_kl(iglob)
     enddo
    close(97)
    close(98)
    close(99)
    close(17)
    close(18)
    close(19)
    endif

    endif

!
!----  PostScript display
!
  if(output_postscript_snapshot) then

  if (myrank == 0) write(IOUT,*) 'Writing PostScript file'

  if(imagetype == 1) then

    if (myrank == 0) write(IOUT,*) 'drawing displacement vector as small arrows...'

    call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic,&
          elastic,poroelastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,numat,kmato,density,rhoext,assign_external_model)

    call plotpost(vector_field_display,coord,vpext,x_source,z_source,st_xval,st_zval, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,porosity,tortuosity,&
          poroelastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,nelem_acoustic_surface,acoustic_edges, &
          simulation_title,npoin,npgeo,vpImin,vpImax,nrec,NSOURCE, &
          colors,numbers,subsamp,imagetype,interpol,meshvect,modelvect, &
          boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,numat,pointsdisp, &
          nspec,ngnod,coupled_acoustic_elastic,coupled_acoustic_poroelastic,coupled_elastic_poroelastic, &
          any_acoustic,any_poroelastic,plot_lowerleft_corner_only, &
          fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,num_fluid_solid_edges,&
          fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge,num_fluid_poro_edges, &
          solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge,num_solid_poro_edges, &
          myrank,nproc,ier,&
          d1_coorg_send_ps_velocity_model,d2_coorg_send_ps_velocity_model, &
          d1_coorg_recv_ps_velocity_model,d2_coorg_recv_ps_velocity_model, &
          d1_RGB_send_ps_velocity_model,d2_RGB_send_ps_velocity_model,d1_RGB_recv_ps_velocity_model,d2_RGB_recv_ps_velocity_model, &
          coorg_send_ps_velocity_model,RGB_send_ps_velocity_model,coorg_recv_ps_velocity_model,RGB_recv_ps_velocity_model, &
          d1_coorg_send_ps_element_mesh,d2_coorg_send_ps_element_mesh,d1_coorg_recv_ps_element_mesh,d2_coorg_recv_ps_element_mesh, &
          d1_color_send_ps_element_mesh,d1_color_recv_ps_element_mesh, &
          coorg_send_ps_element_mesh,color_send_ps_element_mesh,coorg_recv_ps_element_mesh,color_recv_ps_element_mesh, &
          d1_coorg_send_ps_abs,d1_coorg_recv_ps_abs,d2_coorg_send_ps_abs,d2_coorg_recv_ps_abs, &
          coorg_send_ps_abs,coorg_recv_ps_abs, &
          d1_coorg_send_ps_free_surface,d1_coorg_recv_ps_free_surface,d2_coorg_send_ps_free_surface,d2_coorg_recv_ps_free_surface, &
          coorg_send_ps_free_surface,coorg_recv_ps_free_surface, &
          d1_coorg_send_ps_vector_field,d1_coorg_recv_ps_vector_field,d2_coorg_send_ps_vector_field,d2_coorg_recv_ps_vector_field, &
          coorg_send_ps_vector_field,coorg_recv_ps_vector_field)

  else if(imagetype == 2) then

    if (myrank == 0) write(IOUT,*) 'drawing velocity vector as small arrows...'

    call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic,&
          elastic,poroelastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,numat,kmato,density,rhoext,assign_external_model)

    call plotpost(vector_field_display,coord,vpext,x_source,z_source,st_xval,st_zval, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,porosity,tortuosity,&
          poroelastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,nelem_acoustic_surface,acoustic_edges, &
          simulation_title,npoin,npgeo,vpImin,vpImax,nrec,NSOURCE, &
          colors,numbers,subsamp,imagetype,interpol,meshvect,modelvect, &
          boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,numat,pointsdisp, &
          nspec,ngnod,coupled_acoustic_elastic,coupled_acoustic_poroelastic,coupled_elastic_poroelastic, &
          any_acoustic,any_poroelastic,plot_lowerleft_corner_only, &
          fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,num_fluid_solid_edges,&
          fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge,num_fluid_poro_edges, &
          solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge,num_solid_poro_edges, &
          myrank,nproc,ier,&
          d1_coorg_send_ps_velocity_model,d2_coorg_send_ps_velocity_model, &
          d1_coorg_recv_ps_velocity_model,d2_coorg_recv_ps_velocity_model, &
          d1_RGB_send_ps_velocity_model,d2_RGB_send_ps_velocity_model,d1_RGB_recv_ps_velocity_model,d2_RGB_recv_ps_velocity_model, &
          coorg_send_ps_velocity_model,RGB_send_ps_velocity_model,coorg_recv_ps_velocity_model,RGB_recv_ps_velocity_model, &
          d1_coorg_send_ps_element_mesh,d2_coorg_send_ps_element_mesh,d1_coorg_recv_ps_element_mesh,d2_coorg_recv_ps_element_mesh, &
          d1_color_send_ps_element_mesh,d1_color_recv_ps_element_mesh, &
          coorg_send_ps_element_mesh,color_send_ps_element_mesh,coorg_recv_ps_element_mesh,color_recv_ps_element_mesh, &
          d1_coorg_send_ps_abs,d1_coorg_recv_ps_abs,d2_coorg_send_ps_abs,d2_coorg_recv_ps_abs, &
          coorg_send_ps_abs,coorg_recv_ps_abs, &
          d1_coorg_send_ps_free_surface,d1_coorg_recv_ps_free_surface,d2_coorg_send_ps_free_surface,d2_coorg_recv_ps_free_surface, &
          coorg_send_ps_free_surface,coorg_recv_ps_free_surface, &
          d1_coorg_send_ps_vector_field,d1_coorg_recv_ps_vector_field,d2_coorg_send_ps_vector_field,d2_coorg_recv_ps_vector_field, &
          coorg_send_ps_vector_field,coorg_recv_ps_vector_field)

  else if(imagetype == 3) then

    if (myrank == 0) write(IOUT,*) 'drawing acceleration vector as small arrows...'

    call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,accels_poroelastic,&
          elastic,poroelastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,numat,kmato,density,rhoext,assign_external_model)

    call plotpost(vector_field_display,coord,vpext,x_source,z_source,st_xval,st_zval, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,porosity,tortuosity,&
          poroelastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,nelem_acoustic_surface,acoustic_edges, &
          simulation_title,npoin,npgeo,vpImin,vpImax,nrec,NSOURCE, &
          colors,numbers,subsamp,imagetype,interpol,meshvect,modelvect, &
          boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,numat,pointsdisp, &
          nspec,ngnod,coupled_acoustic_elastic,coupled_acoustic_poroelastic,coupled_elastic_poroelastic, &
          any_acoustic,any_poroelastic,plot_lowerleft_corner_only, &
          fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,num_fluid_solid_edges, &
          fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge,num_fluid_poro_edges, &
          solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge,num_solid_poro_edges, &
          myrank,nproc,ier,&
          d1_coorg_send_ps_velocity_model,d2_coorg_send_ps_velocity_model, &
          d1_coorg_recv_ps_velocity_model,d2_coorg_recv_ps_velocity_model, &
          d1_RGB_send_ps_velocity_model,d2_RGB_send_ps_velocity_model,d1_RGB_recv_ps_velocity_model,d2_RGB_recv_ps_velocity_model, &
          coorg_send_ps_velocity_model,RGB_send_ps_velocity_model,coorg_recv_ps_velocity_model,RGB_recv_ps_velocity_model, &
          d1_coorg_send_ps_element_mesh,d2_coorg_send_ps_element_mesh,d1_coorg_recv_ps_element_mesh,d2_coorg_recv_ps_element_mesh, &
          d1_color_send_ps_element_mesh,d1_color_recv_ps_element_mesh, &
          coorg_send_ps_element_mesh,color_send_ps_element_mesh,coorg_recv_ps_element_mesh,color_recv_ps_element_mesh, &
          d1_coorg_send_ps_abs,d1_coorg_recv_ps_abs,d2_coorg_send_ps_abs,d2_coorg_recv_ps_abs, &
          coorg_send_ps_abs,coorg_recv_ps_abs, &
          d1_coorg_send_ps_free_surface,d1_coorg_recv_ps_free_surface,d2_coorg_send_ps_free_surface,d2_coorg_recv_ps_free_surface, &
          coorg_send_ps_free_surface,coorg_recv_ps_free_surface, &
          d1_coorg_send_ps_vector_field,d1_coorg_recv_ps_vector_field,d2_coorg_send_ps_vector_field,d2_coorg_recv_ps_vector_field, &
          coorg_send_ps_vector_field,coorg_recv_ps_vector_field)

  else if(imagetype == 4) then

    if (myrank == 0) write(IOUT,*) 'cannot draw scalar pressure field as a vector plot, skipping...'

  else
    call exit_MPI('wrong type for snapshots')
  endif

  if (myrank == 0 .and. imagetype /= 4) write(IOUT,*) 'PostScript file written'

  endif

!
!----  display color image
!
  if(output_color_image) then

  if (myrank == 0) write(IOUT,*) 'Creating color image of size ',NX_IMAGE_color,' x ',NZ_IMAGE_color

  if(imagetype == 1) then

    if (myrank == 0) write(IOUT,*) 'drawing image of vertical component of displacement vector...'

    call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic,&
          elastic,poroelastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,numat,kmato,density,rhoext,assign_external_model)

  else if(imagetype == 2) then

    if (myrank == 0) write(IOUT,*) 'drawing image of vertical component of velocity vector...'

    call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic,&
          elastic,poroelastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,numat,kmato,density,rhoext,assign_external_model)

  else if(imagetype == 3) then

    if (myrank == 0) write(IOUT,*) 'drawing image of vertical component of acceleration vector...'

    call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,accels_poroelastic,&
          elastic,poroelastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,numat,kmato,density,rhoext,assign_external_model)

  else if(imagetype == 4) then

    if (myrank == 0) write(IOUT,*) 'drawing image of pressure field...'

    call compute_pressure_whole_medium(potential_dot_dot_acoustic,displ_elastic,&
         displs_poroelastic,displw_poroelastic,elastic,poroelastic,vector_field_display, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,assign_external_model, &
         numat,kmato,density,porosity,tortuosity,poroelastcoef,vpext,vsext,rhoext,e1,e11, &
         TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,Mu_nu1,Mu_nu2,N_SLS)

  else
    call exit_MPI('wrong type for snapshots')
  endif

  image_color_data(:,:) = 0.d0

  do k = 1, nb_pixel_loc
     j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
     i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color
     image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))
  enddo

! assembling array image_color_data on process zero for color output
#ifdef USE_MPI
  if (nproc > 1) then
     if (myrank == 0) then

        do iproc = 1, nproc-1
           call MPI_RECV(data_pixel_recv(1),nb_pixel_per_proc(iproc+1), MPI_DOUBLE_PRECISION, &
                iproc, 43, MPI_COMM_WORLD, request_mpi_status, ier)

           do k = 1, nb_pixel_per_proc(iproc+1)
              j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
              i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color
              image_color_data(i,j) = data_pixel_recv(k)
           enddo
        enddo

     else
        do k = 1, nb_pixel_loc
           j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
           i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color
           data_pixel_send(k) = vector_field_display(2,iglob_image_color(i,j))
        enddo

        call MPI_SEND(data_pixel_send(1),nb_pixel_loc,MPI_DOUBLE_PRECISION, 0, 43, MPI_COMM_WORLD, ier)

     endif
  endif

#endif

  if (myrank == 0) then
     call create_color_image(image_color_data,iglob_image_color,NX_IMAGE_color,NZ_IMAGE_color,it,cutsnaps,image_color_vp_display)
     write(IOUT,*) 'Color image created'
  endif

  endif

!----  save temporary or final seismograms
! suppress seismograms if we generate traces of the run for analysis with "ParaVer", because time consuming
  if(.not. GENERATE_PARAVER_TRACES) call write_seismograms(sisux,sisuz,siscurl,station_name,network_name,NSTEP, &
        nrecloc,which_proc_receiver,nrec,myrank,deltat,seismotype,st_xval,t0, &
        NTSTEP_BETWEEN_OUTPUT_SEISMO,seismo_offset,seismo_current)

  seismo_offset = seismo_offset + seismo_current
  seismo_current = 0

! count elapsed wall-clock time
  call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! this fails if we cross the end of the month
  time_end = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
             60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0

! elapsed time since beginning of the simulation
  tCPU = time_end - time_start
  int_tCPU = int(tCPU)
  ihours = int_tCPU / 3600
  iminutes = (int_tCPU - 3600*ihours) / 60
  iseconds = int_tCPU - 3600*ihours - 60*iminutes
  if (myrank == 0) then
    write(IOUT,*) 'Elapsed time in seconds = ',tCPU
    write(IOUT,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(IOUT,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
    write(IOUT,*)
  endif

  endif

#ifdef USE_MPI
! add a barrier if we generate traces of the run for analysis with "ParaVer"
  if(GENERATE_PARAVER_TRACES) call MPI_BARRIER(MPI_COMM_WORLD,ier)
#endif

  enddo ! end of the main time loop

  if((save_forward .and. isolver==1) .or. isolver ==2) then
   if(any_acoustic) then
  close(65)
  close(66)
  close(67)
  close(68)
   endif
   if(any_elastic) then
  close(35)
  close(36)
  close(37)
  close(38)
   endif
   if(any_poroelastic) then
  close(25)
  close(45)
  close(26)
  close(46)
  close(29)
  close(47)
  close(28)
  close(48)
   endif
  endif

!
!--- save last frame 
!
  if(save_forward .and. isolver ==1 .and. any_elastic) then
  if ( myrank == 0 ) then
    write(IOUT,*)
    write(IOUT,*) 'Saving elastic last frame...'
    write(IOUT,*)
  endif
    open(unit=55,file='OUTPUT_FILES/lastframe_elastic.bin',status='unknown',form='unformatted')
       do j=1,npoin
      write(55) (displ_elastic(i,j), i=1,NDIM), &
                  (veloc_elastic(i,j), i=1,NDIM), &
                  (accel_elastic(i,j), i=1,NDIM)
       enddo
    close(55)
  endif

  if(save_forward .and. isolver ==1 .and. any_poroelastic) then
  if ( myrank == 0 ) then
    write(IOUT,*)
    write(IOUT,*) 'Saving poroelastic last frame...'
    write(IOUT,*)
  endif
    open(unit=55,file='OUTPUT_FILES/lastframe_poroelastic_s.bin',status='unknown',form='unformatted')
    open(unit=56,file='OUTPUT_FILES/lastframe_poroelastic_w.bin',status='unknown',form='unformatted')
       do j=1,npoin
      write(55) (displs_poroelastic(i,j), i=1,NDIM), &
                  (velocs_poroelastic(i,j), i=1,NDIM), &
                  (accels_poroelastic(i,j), i=1,NDIM)
      write(56) (displw_poroelastic(i,j), i=1,NDIM), &
                  (velocw_poroelastic(i,j), i=1,NDIM), &
                  (accelw_poroelastic(i,j), i=1,NDIM)
       enddo
    close(55)
    close(56)
  endif

  if(save_forward .and. isolver ==1 .and. any_acoustic) then
  if ( myrank == 0 ) then
    write(IOUT,*)
    write(IOUT,*) 'Saving acoustic last frame...'
    write(IOUT,*)
  endif
    open(unit=55,file='OUTPUT_FILES/lastframe_acoustic.bin',status='unknown',form='unformatted')
       do j=1,npoin
      write(55) potential_acoustic(j),&
               potential_dot_acoustic(j),&
               potential_dot_dot_acoustic(j)
       enddo
    close(55)
  endif


  deallocate(v0x_left)
  deallocate(v0z_left)
  deallocate(t0x_left)
  deallocate(t0z_left)

  deallocate(v0x_right)
  deallocate(v0z_right)
  deallocate(t0x_right)
  deallocate(t0z_right)

  deallocate(v0x_bot)
  deallocate(v0z_bot)
  deallocate(t0x_bot)
  deallocate(t0z_bot)

!----  close energy file and create a gnuplot script to display it
  if(OUTPUT_ENERGY .and. myrank == 0) then
    close(IENERGY)
    open(unit=IENERGY,file='plotenergy',status='unknown')
    write(IENERGY,*) 'set term postscript landscape color solid "Helvetica" 22'
    write(IENERGY,*) 'set output "energy.ps"'
    write(IENERGY,*) 'set xlabel "Time (s)"'
    write(IENERGY,*) 'set ylabel "Energy (J)"'
    write(IENERGY,*) 'plot "energy.gnu" us 1:4 t ''Total Energy'' w l 1, "energy.gnu" us 1:3 t ''Potential Energy'' w l 2'
    close(IENERGY)
  endif

! print exit banner
  if (myrank == 0) call datim(simulation_title)

!
!----  close output file
!
  if(IOUT /= ISTANDARD_OUTPUT) close(IOUT)

!
!----  end MPI
!
#ifdef USE_MPI
  call MPI_FINALIZE(ier)
#endif

!
!----  formats
!

 400 format(/1x,41('=')/,' =  T i m e  e v o l u t i o n  l o o p  ='/1x,41('=')/)

 200 format(//1x,'C o n t r o l',/1x,13('='),//5x,&
  'Number of spectral element control nodes. . .(npgeo) =',i8/5x, &
  'Number of space dimensions. . . . . . . . . . (NDIM) =',i8)

 600 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Display frequency . . . (NTSTEP_BETWEEN_OUTPUT_INFO) = ',i6/ 5x, &
  'Color display . . . . . . . . . . . . . . . (colors) = ',i6/ 5x, &
  '        ==  0     black and white display              ',  / 5x, &
  '        ==  1     color display                        ',  /5x, &
  'Numbered mesh . . . . . . . . . . . . . . .(numbers) = ',i6/ 5x, &
  '        ==  0     do not number the mesh               ',  /5x, &
  '        ==  1     number the mesh                      ')

 700 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Seismograms recording type . . . . . . .(seismotype) = ',i6/5x, &
  'Angle for first line of receivers. . . . .(anglerec) = ',f6.2)

 750 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Read external initial field. . . . . .(initialfield) = ',l6/5x, &
  'Add Bielak conditions . . . .(add_Bielak_conditions) = ',l6/5x, &
  'Assign external model . . . .(assign_external_model) = ',l6/5x, &
  'Turn anisotropy on or off. . . .(TURN_ANISOTROPY_ON) = ',l6/5x, &
  'Turn attenuation on or off. . .(TURN_ATTENUATION_ON) = ',l6/5x, &
  'Save grid in external file or not. . . .(outputgrid) = ',l6)

 800 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Vector display type . . . . . . . . . . .(imagetype) = ',i6/5x, &
  'Percentage of cut for vector plots . . . .(cutsnaps) = ',f6.2/5x, &
  'Subsampling for velocity model display. . .(subsamp) = ',i6)

 703 format(//' I t e r a t i o n s '/1x,19('='),//5x, &
      'Number of time iterations . . . . .(NSTEP) =',i8,/5x, &
      'Time step increment. . . . . . . .(deltat) =',1pe15.6,/5x, &
      'Total simulation duration . . . . . (ttot) =',1pe15.6)

 107 format(/5x,'--> Isoparametric Spectral Elements <--',//)

 207 format(5x,'Number of spectral elements . . . . .  (nspec) =',i7,/5x, &
               'Number of control nodes per element .  (ngnod) =',i7,/5x, &
               'Number of points in X-direction . . .  (NGLLX) =',i7,/5x, &
               'Number of points in Y-direction . . .  (NGLLZ) =',i7,/5x, &
               'Number of points per element. . .(NGLLX*NGLLZ) =',i7,/5x, &
               'Number of points for display . . .(pointsdisp) =',i7,/5x, &
               'Number of element material sets . . .  (numat) =',i7,/5x, &
               'Number of absorbing elements . . . .(nelemabs) =',i7)

 212 format(//,5x,'Source Type. . . . . . . . . . . . . . = Collocated Force',/5x, &
                  'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
                  'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
                  'Angle from vertical direction (deg). . =',1pe20.10,/5x)

 222 format(//,5x,'Source Type. . . . . . . . . . . . . . = Moment-tensor',/5x, &
                  'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
                  'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mxx. . . . . . . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mzz. . . . . . . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mxz. . . . . . . . . . . . . . . . . . =',1pe20.10)

  end program specfem2D


subroutine tri_quad(n, n1, nnodes)

      implicit none

      integer  :: n1, nnodes
      integer, dimension(4)  :: n


      n(2) = n1

      if ( n1 == 1 ) then
         n(1) = nnodes
      else
         n(1) = n1-1
      endif

      if ( n1 == nnodes ) then
         n(3) = 1
      else
         n(3) = n1+1
      endif

      if ( n(3) == nnodes ) then
         n(4) = 1
      else
         n(4) = n(3)+1
      endif


end subroutine tri_quad


subroutine calcul_normale( angle, n1_x, n2_x, n3_x, n4_x, n1_z, n2_z, n3_z, n4_z )

      implicit none

      include 'constants.h'

      double precision :: angle
      double precision :: n1_x, n2_x, n3_x, n4_x, n1_z, n2_z, n3_z, n4_z

      double precision  :: theta1, theta2, theta3
      double precision  :: costheta1, costheta2, costheta3

      if ( abs(n2_z - n1_z) < TINYVAL ) then
         costheta1 = 0
      else
         costheta1 = (n2_z - n1_z) / sqrt((n2_x - n1_x)**2 + (n2_z - n1_z)**2)
      endif
      if ( abs(n3_z - n2_z) < TINYVAL ) then
         costheta2 = 0
      else
         costheta2 = (n3_z - n2_z) / sqrt((n3_x - n2_x)**2 + (n3_z - n2_z)**2)
      endif
      if ( abs(n4_z - n3_z) < TINYVAL ) then
         costheta3 = 0
      else
        costheta3 = (n4_z - n3_z) / sqrt((n4_x - n3_x)**2 + (n4_z - n3_z)**2)
      endif

      theta1 = - sign(1.d0,n2_x - n1_x) * acos(costheta1)
      theta2 = - sign(1.d0,n3_x - n2_x) * acos(costheta2)
      theta3 = - sign(1.d0,n4_x - n3_x) * acos(costheta3)


      angle = angle + ( theta1 + theta2 + theta3 ) / 3.d0 + PI/2.d0
      !angle = theta2

end subroutine calcul_normale
