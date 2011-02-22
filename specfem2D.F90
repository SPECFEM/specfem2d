
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
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
!   for the anelastic anisotropic or poroelastic wave equation.
!
!====================================================================================

! If you use this code for your own research, please cite at least one article
! written by the developers of the package, for instance:
!
! @ARTICLE{TrKoLi08,
! author = {Jeroen Tromp and Dimitri Komatitsch and Qinya Liu},
! title = {Spectral-Element and Adjoint Methods in Seismology},
! journal = {Communications in Computational Physics},
! year = {2008},
! volume = {3},
! pages = {1-32},
! number = {1}}
!
! or
!
! @ARTICLE{VaCaSaKoVi99,
! author = {R. Vai and J. M. Castillo-Covarrubias and F. J. S\'anchez-Sesma and
! D. Komatitsch and J. P. Vilotte},
! title = {Elastic wave propagation in an irregularly layered medium},
! journal = {Soil Dynamics and Earthquake Engineering},
! year = {1999},
! volume = {18},
! pages = {11-18},
! number = {1},
! doi = {10.1016/S0267-7261(98)00027-X}}
!
! @ARTICLE{LeChKoHuTr09,
! author = {Shiann Jong Lee and Yu Chang Chan and Dimitri Komatitsch and Bor
! Shouh Huang and Jeroen Tromp},
! title = {Effects of realistic surface topography on seismic ground motion
! in the {Y}angminshan region of {T}aiwan based upon the spectral-element
! method and {LiDAR DTM}},
! journal = {Bull. Seismol. Soc. Am.},
! year = {2009},
! volume = {99},
! pages = {681-693},
! number = {2A},
! doi = {10.1785/0120080264}}
!
! @ARTICLE{LeChLiKoHuTr08,
! author = {Shiann Jong Lee and How Wei Chen and Qinya Liu and Dimitri Komatitsch
! and Bor Shouh Huang and Jeroen Tromp},
! title = {Three-Dimensional Simulations of Seismic Wave Propagation in the
! {T}aipei Basin with Realistic Topography Based upon the Spectral-Element Method},
! journal = {Bull. Seismol. Soc. Am.},
! year = {2008},
! volume = {98},
! pages = {253-264},
! number = {1},
! doi = {10.1785/0120070033}}
!
! @ARTICLE{LeKoHuTr09,
! author = {S. J. Lee and Dimitri Komatitsch and B. S. Huang and J. Tromp},
! title = {Effects of topography on seismic wave propagation: An example from
! northern {T}aiwan},
! journal = {Bull. Seismol. Soc. Am.},
! year = {2009},
! volume = {99},
! pages = {314-325},
! number = {1},
! doi = {10.1785/0120080020}}
!
! @ARTICLE{KoErGoMi10,
! author = {Dimitri Komatitsch and Gordon Erlebacher and Dominik G\"oddeke and
! David Mich\'ea},
! title = {High-order finite-element seismic wave propagation modeling with
! {MPI} on a large {GPU} cluster},
! journal = {J. Comput. Phys.},
! year = {2010},
! volume = {229},
! pages = {7692-7714},
! number = {20},
! doi = {10.1016/j.jcp.2010.06.024}}
!
! @ARTICLE{KoGoErMi10,
! author = {Dimitri Komatitsch and Dominik G\"oddeke and Gordon Erlebacher and
! David Mich\'ea},
! title = {Modeling the propagation of elastic waves using spectral elements
! on a cluster of 192 {GPU}s},
! journal = {Computer Science Research and Development},
! year = {2010},
! volume = {25},
! pages = {75-82},
! number = {1-2},
! doi = {10.1007/s00450-010-0109-1}}
!
! @ARTICLE{KoMiEr09,
! author = {Dimitri Komatitsch and David Mich\'ea and Gordon Erlebacher},
! title = {Porting a high-order finite-element earthquake modeling application
! to {NVIDIA} graphics cards using {CUDA}},
! journal = {Journal of Parallel and Distributed Computing},
! year = {2009},
! volume = {69},
! pages = {451-460},
! number = {5},
! doi = {10.1016/j.jpdc.2009.01.006}}
!
! @ARTICLE{LiPoKoTr04,
! author = {Qinya Liu and Jascha Polet and Dimitri Komatitsch and Jeroen Tromp},
! title = {Spectral-element moment tensor inversions for earthquakes in {S}outhern {C}alifornia},
! journal={Bull. Seismol. Soc. Am.},
! year = {2004},
! volume = {94},
! pages = {1748-1761},
! number = {5},
! doi = {10.1785/012004038}}
!
! @INCOLLECTION{ChKoViCaVaFe07,
! author = {Emmanuel Chaljub and Dimitri Komatitsch and Jean-Pierre Vilotte and
! Yann Capdeville and Bernard Valette and Gaetano Festa},
! title = {Spectral Element Analysis in Seismology},
! booktitle = {Advances in Wave Propagation in Heterogeneous Media},
! publisher = {Elsevier - Academic Press},
! year = {2007},
! editor = {Ru-Shan Wu and Val\'erie Maupin},
! volume = {48},
! series = {Advances in Geophysics},
! pages = {365-419}}
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
! @ARTICLE{KoTr99,
! author={D. Komatitsch and J. Tromp},
! year=1999,
! title={Introduction to the spectral-element method for 3-{D} seismic wave propagation},
! journal={Geophys. J. Int.},
! volume=139,
! number=3,
! pages={806-822},
! doi={10.1046/j.1365-246x.1999.00967.x}}
!
! @ARTICLE{KoLiTrSuStSh04,
! author={Dimitri Komatitsch and Qinya Liu and Jeroen Tromp and Peter S\"{u}ss
!   and Christiane Stidham and John H. Shaw},
! year=2004,
! title={Simulations of Ground Motion in the {L}os {A}ngeles {B}asin
!   based upon the Spectral-Element Method},
! journal={Bull. Seism. Soc. Am.},
! volume=94,
! number=1,
! pages={187-206}}
!
! @ARTICLE{MoTr08,
! author={C. Morency and J. Tromp},
! title={Spectral-element simulations of wave propagation in poroelastic media},
! journal={Geophys. J. Int.},
! year=2008,
! volume=175,
! pages={301-345}}
!
! and/or other articles from http://web.univ-pau.fr/~dkomati1/publications.html
!
! If you use the kernel capabilities of the code, please cite at least one article
! written by the developers of the package, for instance:
!
! @ARTICLE{TrKoLi08,
! author = {Jeroen Tromp and Dimitri Komatitsch and Qinya Liu},
! title = {Spectral-Element and Adjoint Methods in Seismology},
! journal = {Communications in Computational Physics},
! year = {2008},
! volume = {3},
! pages = {1-32},
! number = {1}}
!
! or
!
! @ARTICLE{LiTr06,
! author={Qinya Liu and Jeroen Tromp},
! title={Finite-frequency kernels based on adjoint methods},
! journal={Bull. Seismol. Soc. Am.},
! year=2006,
! volume=96,
! number=6,
! pages={2383-2397},
! doi={10.1785/0120060041}}
!
! @ARTICLE{MoLuTr09,
! author={C. Morency and Y. Luo and J. Tromp},
! title={Finite-frequency kernels for wave propagation in porous media based upon adjoint methods},
! year=2009,
! journal={Geophys. J. Int.},
! doi={10.1111/j.1365-246X.2009.04332}}
!
! If you use the METIS / SCOTCH / CUBIT non-structured capabilities, please also cite:
!
! @ARTICLE{MaKoBlLe08,
! author = {R. Martin and D. Komatitsch and C. Blitz and N. {Le Goff}},
! title = {Simulation of seismic wave propagation in an asteroid based upon
! an unstructured {MPI} spectral-element method: blocking and non-blocking
! communication strategies},
! journal = {Lecture Notes in Computer Science},
! year = {2008},
! volume = {5336},
! pages = {350-363}}
!
! version 6.1, Christina Morency and Pieyre Le Loher, March 2010:
!               - added SH (membrane) waves calculation for elastic media
!               - added support for external fully anisotropic media
!               - fixed some bugs in acoustic kernels
!
! version 6.0, Christina Morency and Yang Luo, August 2009:
!               - support for poroelastic media
!               - adjoint method for acoustic/elastic/poroelastic
!
! version 5.2, Dimitri Komatitsch, Nicolas Le Goff and Roland Martin, February 2008:
!               - support for CUBIT and GiD meshes
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

!  character(len=80) datlin

  integer NSOURCES,i_source
  integer, dimension(:), allocatable :: source_type,time_function_type
  double precision, dimension(:), allocatable :: x_source,z_source,xi_source,gamma_source,&
                  Mxx,Mzz,Mxz,f0,tshift_src,factor,angleforce 
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: sourcearray
  double precision :: t0

  double precision, dimension(:,:), allocatable :: coorg

! for P-SV or SH (membrane) waves calculation
  logical :: p_sv

! receiver information
  integer :: nrec,ios
  integer, dimension(:), allocatable :: ispec_selected_rec
  double precision, dimension(:), allocatable :: xi_receiver,gamma_receiver,st_xval,st_zval
  character(len=150) dummystring

! for seismograms
  double precision, dimension(:,:), allocatable :: sisux,sisuz,siscurl
  integer :: seismo_offset, seismo_current

! vector field in an element
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLX) :: vector_field_element

! pressure in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: pressure_element

! curl in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: curl_element

  integer :: i,j,k,it,irec,id,n,ispec,npoin,npgeo,iglob 
  logical :: anyabs
  double precision :: dxd,dyd,dzd,dcurld,valux,valuy,valuz,valcurl,hlagrange,rhol,xi,gamma,x,z

! coefficients of the explicit Newmark time scheme
  integer NSTEP
  double precision :: deltatover2,deltatsquareover2,time
  double precision :: deltat

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
  double precision, dimension(:,:), allocatable :: &
    coord, flagrange,xinterp,zinterp,Uxinterp,Uzinterp,vector_field_display

! material properties of the poroelastic medium (solid phase:s and fluid phase [defined as w=phi(u_f-u_s)]: w)
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    accels_poroelastic,velocs_poroelastic,displs_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    accelw_poroelastic,velocw_poroelastic,displw_poroelastic
  double precision, dimension(:), allocatable :: porosity,tortuosity
  double precision, dimension(:,:), allocatable :: density,permeability

! poroelastic and elastic coefficients
  double precision, dimension(:,:,:), allocatable :: poroelastcoef

! anisotropy parameters
  logical :: all_anisotropic
  double precision ::  c11,c13,c15,c33,c35,c55
  logical, dimension(:), allocatable :: anisotropic
  double precision, dimension(:,:), allocatable :: anisotropy

! for acoustic medium
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic

! inverse mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_elastic,rmass_inverse_acoustic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic

! to evaluate cpI, cpII, and cs, and rI (poroelastic medium)
  real(kind=CUSTOM_REAL) :: rhol_s,rhol_f,rhol_bar,phil,tortl
  real(kind=CUSTOM_REAL) :: mul_s,kappal_s
  real(kind=CUSTOM_REAL) :: kappal_f
  real(kind=CUSTOM_REAL) :: mul_fr,kappal_fr
  real(kind=CUSTOM_REAL) :: D_biot,H_biot,C_biot,M_biot,B_biot,cpIsquare,cpIIsquare,cssquare
  real(kind=CUSTOM_REAL) :: ratio,dd1 

  double precision, dimension(:,:,:), allocatable :: vpext,vsext,rhoext
  double precision, dimension(:,:,:), allocatable :: Qp_attenuationext,Qs_attenuationext
  double precision, dimension(:,:,:), allocatable :: c11ext,c13ext,c15ext,c33ext,c35ext,c55ext

  double precision, dimension(:,:,:), allocatable :: shape2D,shape2D_display
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: xix,xiz,gammax,gammaz,jacobian

  double precision, dimension(:,:,:,:), allocatable :: dershape2D,dershape2D_display

  integer, dimension(:,:,:), allocatable :: ibool,ibool_outer,ibool_inner
  integer, dimension(:,:), allocatable  :: knods
  integer, dimension(:), allocatable :: kmato,numabs, &
     ibegin_bottom,iend_bottom,ibegin_top,iend_top,jbegin_left,jend_left,jbegin_right,jend_right

  integer, dimension(:), allocatable :: ispec_selected_source,iglob_source,&
                                        is_proc_source,nb_proc_source
  double precision, dimension(:), allocatable :: aval
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: source_time_function
  double precision, external :: netlib_specfun_erf

  double precision :: vpImin,vpImax,vpIImin,vpIImax

  integer :: colors,numbers,subsamp,imagetype, &
    NTSTEP_BETWEEN_OUTPUT_INFO,NTSTEP_BETWEEN_OUTPUT_SEISMO,seismotype
  integer :: numat,ngnod,nspec,pointsdisp, &
    nelemabs,nelem_acoustic_surface,ispecabs,UPPER_LIMIT_DISPLAY

  logical interpol,meshvect,modelvect,boundvect,assign_external_model,initialfield, &
    outputgrid,gnuplot,TURN_ATTENUATION_ON,output_postscript_snapshot,output_color_image, &
    plot_lowerleft_corner_only,add_Bielak_conditions,OUTPUT_ENERGY,READ_EXTERNAL_SEP_FILE

  double precision :: cutsnaps,sizemax_arrows,anglerec,xirec,gammarec

! for absorbing and acoustic free surface conditions
  integer :: ispec_acoustic_surface,inum 
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
  double precision :: theta_e,theta_s
  double precision :: Q0,freq0
  double precision :: alphaval,betaval,gammaval,thetainv
  logical :: TURN_VISCATTENUATION_ON
  double precision, dimension(NGLLX,NGLLZ) :: viscox_loc,viscoz_loc
  double precision :: Sn,Snp1,etal_f
  double precision, dimension(3):: bl_relaxed
  double precision :: permlxx,permlxz,permlzz,invpermlxx,invpermlxz,invpermlzz,detk
! adjoint
  double precision, dimension(:), allocatable :: b_viscodampx,b_viscodampz
  integer reclen

! for fluid/solid coupling and edge detection
  logical, dimension(:), allocatable :: elastic
  integer, dimension(NEDGES) :: i_begin,j_begin,i_end,j_end
  integer, dimension(NGLLX,NEDGES) :: ivalue,jvalue,ivalue_inverse,jvalue_inverse
  integer, dimension(:), allocatable :: fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
                                        fluid_solid_elastic_ispec,fluid_solid_elastic_iedge
  integer :: num_fluid_solid_edges,ispec_acoustic,ispec_elastic, &
             iedge_acoustic,iedge_elastic,ipoin1D,iglob2
  logical :: any_acoustic,any_acoustic_glob,any_elastic,any_elastic_glob,coupled_acoustic_elastic
  real(kind=CUSTOM_REAL) :: displ_x,displ_z,displ_n,displw_x,displw_z,zxi,xgamma,jacobian1D,pressure
  real(kind=CUSTOM_REAL) :: b_displ_x,b_displ_z,b_displw_x,b_displw_z,b_pressure
  logical :: any_fluid_solid_edges

! for fluid/porous medium coupling and edge detection
  logical, dimension(:), allocatable :: poroelastic
  logical :: any_poroelastic,any_poroelastic_glob
  integer, dimension(:), allocatable :: fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge, &
                                        fluid_poro_poroelastic_ispec,fluid_poro_poroelastic_iedge
  integer :: num_fluid_poro_edges,iedge_poroelastic
  logical :: coupled_acoustic_poro
  double precision :: mul_G,lambdal_G,lambdalplus2mul_G
  double precision :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  double precision :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  double precision :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  double precision :: dwx_dxl,dwz_dxl,dwx_dzl,dwz_dzl
  double precision :: b_dux_dxi,b_dux_dgamma,b_duz_dxi,b_duz_dgamma
  double precision :: b_dwx_dxi,b_dwx_dgamma,b_dwz_dxi,b_dwz_dgamma
  double precision :: b_dux_dxl,b_duz_dxl,b_dux_dzl,b_duz_dzl
  double precision :: b_dwx_dxl,b_dwz_dxl,b_dwx_dzl,b_dwz_dzl
  logical :: any_fluid_poro_edges

! for solid/porous medium coupling and edge detection
  integer, dimension(:), allocatable :: solid_poro_elastic_ispec,solid_poro_elastic_iedge, &
                                        solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge
  integer :: num_solid_poro_edges,ispec_poroelastic,ii2,jj2
  logical :: coupled_elastic_poro
  integer, dimension(:), allocatable :: icount
  double precision :: sigma_xx,sigma_xz,sigma_zz,sigmap
  double precision :: b_sigma_xx,b_sigma_xz,b_sigma_zz,b_sigmap
  integer, dimension(:), allocatable :: ibegin_bottom_poro,iend_bottom_poro,ibegin_top_poro,&
            iend_top_poro,jbegin_left_poro,jend_left_poro,jbegin_right_poro,jend_right_poro
  logical :: any_solid_poro_edges
  
! for adjoint method
  logical :: SAVE_FORWARD ! whether or not the last frame is saved to reconstruct the forward field
  integer :: SIMULATION_TYPE      ! 1 = forward wavefield, 2 = backward and adjoint wavefields and kernels
  double precision :: b_deltatover2,b_deltatsquareover2,b_deltat ! coefficients of the explicit Newmark time scheme
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_accels_poroelastic,b_velocs_poroelastic,b_displs_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_accelw_poroelastic,b_velocw_poroelastic,b_displw_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_accel_elastic,b_veloc_elastic,b_displ_elastic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_ac,b_displ_ac,b_accel_ac
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rho_kl, mu_kl, kappa_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhol_global, mul_global, kappal_global
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: mu_k, kappa_k,rho_k
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhop_kl, beta_kl, alpha_kl
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rho_ac_kl, kappa_ac_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhol_ac_global, kappal_ac_global
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhop_ac_kl, alpha_ac_kl
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhot_kl, rhof_kl, sm_kl, eta_kl, mufr_kl, B_kl, &
    C_kl, M_kl, rhob_kl, rhofb_kl, phi_kl, Bb_kl, Cb_kl, Mb_kl, mufrb_kl, &
    rhobb_kl, rhofbb_kl, phib_kl, cpI_kl, cpII_kl, cs_kl, ratio_kl
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhot_k, rhof_k, sm_k, eta_k, mufr_k, B_k, &
    C_k, M_k
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: phil_global,etal_f_global,rhol_s_global,rhol_f_global,rhol_bar_global, &
    tortl_global,mulfr_global
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: permlxx_global,permlxz_global,permlzz_global
  character(len=150) :: adj_source_file
  integer :: irec_local,nadj_rec_local
  double precision :: xx,zz,rholb,tempx1l,tempx2l,b_tempx1l,b_tempx2l,bb_tempx1l,bb_tempx2l
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: adj_sourcearray
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: adj_sourcearrays
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_absorb_elastic_left,b_absorb_poro_s_left,b_absorb_poro_w_left
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_absorb_elastic_right,b_absorb_poro_s_right,b_absorb_poro_w_right
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_absorb_elastic_bottom,b_absorb_poro_s_bottom,b_absorb_poro_w_bottom
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_absorb_elastic_top,b_absorb_poro_s_top,b_absorb_poro_w_top
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable ::  b_absorb_acoustic_left,b_absorb_acoustic_right,&
                      b_absorb_acoustic_bottom, b_absorb_acoustic_top
  integer :: nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax
  integer, dimension(:), allocatable :: ib_left,ib_right,ib_bottom,ib_top

! for color images
  integer :: NX_IMAGE_color,NZ_IMAGE_color
  double precision :: xmin_color_image,xmax_color_image, &
    zmin_color_image,zmax_color_image 
  integer, dimension(:,:), allocatable :: iglob_image_color,copy_iglob_image_color
  double precision, dimension(:,:), allocatable :: image_color_data
  double precision, dimension(:,:), allocatable :: image_color_vp_display
  integer  :: nb_pixel_loc
  integer, dimension(:), allocatable  :: num_pixel_loc

#ifdef USE_MPI
  integer, dimension(MPI_STATUS_SIZE)  :: request_mpi_status
  integer, dimension(:), allocatable  :: nb_pixel_per_proc
  integer, dimension(:,:), allocatable  :: num_pixel_recv
  double precision, dimension(:), allocatable  :: data_pixel_recv
  double precision, dimension(:), allocatable  :: data_pixel_send
#endif

! timing information for the stations
  character(len=MAX_LENGTH_STATION_NAME), allocatable, dimension(:) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), allocatable, dimension(:) :: network_name

! title of the plot
  character(len=60) simulation_title

! Lagrange interpolators at receivers
  double precision, dimension(:), allocatable :: hxir,hgammar,hpxir,hpgammar
  double precision, dimension(:,:), allocatable :: hxir_store,hgammar_store

! Lagrange interpolators at sources
  double precision, dimension(:), allocatable :: hxis,hgammas,hpxis,hpgammas
  double precision, dimension(:,:), allocatable :: hxis_store,hgammas_store

! for Lagrange interpolants
  double precision, external :: hgll

! timer to count elapsed time
  double precision :: time_start 
  integer :: year_start,month_start 

  ! to determine date and time at which the run will finish
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values

! for MPI and partitioning
  integer  :: ier
  integer  :: nproc
  integer  :: myrank
  character(len=150) :: outputname,outputname2

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
  integer, dimension(:), allocatable  :: tab_requests_send_recv_poro
  integer :: max_ibool_interfaces_size_ac, max_ibool_interfaces_size_el, max_ibool_interfaces_size_po
  integer :: iproc
#endif

! for overlapping MPI communications with computation
  integer  :: nspec_outer, nspec_inner, num_ispec_outer, num_ispec_inner
  integer, dimension(:), allocatable  :: ispec_outer_to_glob, ispec_inner_to_glob
  logical, dimension(:), allocatable  :: mask_ispec_inner_outer

  integer, dimension(:,:), allocatable  :: acoustic_surface
  integer, dimension(:,:), allocatable  :: acoustic_edges
  logical :: any_acoustic_edges
  
  integer  :: ixmin, ixmax, izmin, izmax

  integer  :: nrecloc, irecloc
  integer, dimension(:), allocatable :: recloc, which_proc_receiver

! mask to sort ibool
  integer, dimension(:), allocatable :: mask_ibool
  integer, dimension(:,:,:), allocatable :: copy_ibool_ori
  integer :: inumber

! to compute analytical initial plane wave field
  double precision :: angleforce_refl, c_inc, c_refl, cploc, csloc 
  double precision, dimension(2) :: A_plane, B_plane, C_plane
  double precision :: z0_source, x0_source, time_offset 

! beyond critical angle
  integer , dimension(:), allocatable :: left_bound,right_bound,bot_bound
  double precision , dimension(:,:), allocatable :: v0x_left,v0z_left,v0x_right,v0z_right,v0x_bot,v0z_bot
  double precision , dimension(:,:), allocatable :: t0x_left,t0z_left,t0x_right,t0z_right,t0x_bot,t0z_bot
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_paco,veloc_paco,displ_paco
  integer count_left,count_right,count_bottom 
  logical :: over_critical_angle

! further reduce cache misses inner/outer in two passes in the case of an MPI simulation
  integer :: ipass,ispec_inner,ispec_outer,NUMBER_OF_PASSES
  integer :: npoin_outer,npoin_inner
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

  integer, dimension(:), allocatable :: source_courbe_eros

  integer  :: nnodes_tangential_curve
  double precision, dimension(:,:), allocatable  :: nodes_tangential_curve
  logical  :: any_tangential_curve

  integer  :: n1_tangential_detection_curve
  integer, dimension(4)  :: n_tangential_detection_curve
  integer, dimension(:), allocatable  :: rec_tangential_detection_curve
  double precision :: distmin, dist_current, angleforce_recv
  double precision, dimension(:), allocatable :: dist_tangential_detection_curve
  double precision :: x_final_receiver_dummy, z_final_receiver_dummy
!!!!!!!!!!
  double precision, dimension(:,:,:),allocatable:: rho_local,vp_local,vs_local
!!!! hessian
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhorho_el_hessian_final1, rhorho_el_hessian_final2
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhorho_el_hessian_temp1, rhorho_el_hessian_temp2
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhorho_ac_hessian_final1, rhorho_ac_hessian_final2
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: weight_line_x, weight_line_z, weight_surface,weight_jacobian
  integer, dimension(:), allocatable :: weight_gll
  real(kind=CUSTOM_REAL) :: zmin_yang, zmax_yang, xmin_yang, xmax_yang

! to help locate elements with a negative Jacobian using OpenDX
  logical :: found_a_negative_jacobian

!! DK DK Feb 2010 for periodic conditions: detect common points between left and right edges
  logical, parameter :: ADD_PERIODIC_CONDITIONS = .false.

!! DK DK the periodic conditions below are currently specific to a Gmsh model designed by Paul Cristini

!! DK DK the horizontal periodicity distance is:
  double precision, parameter :: PERIODIC_horiz_dist =   0.3597d0

!! DK DK the length of an edge is about 1d-003, thus use e.g. 1/300 of that
  double precision, parameter :: PERIODIC_DETECT_TOL = 1d-003 / 300.d0

  integer, parameter :: NSPEC_PERIO = 670 / 2  ! 414 / 2

  integer, dimension(NSPEC_PERIO) :: numperio_left
  integer, dimension(NSPEC_PERIO) :: numperio_right

  logical, dimension(4,NSPEC_PERIO) :: codeabs_perio_left
  logical, dimension(4,NSPEC_PERIO) :: codeabs_perio_right

  integer :: idummy1, idummy2, idummy3, idummy4, idummy5, idummy6, idummy7, idummy8
  integer :: ispecperio, ispecperio2, ispec2, i2, j2
  integer :: iglob_target_to_replace, ispec3, i3, j3

!! DK DK Feb 2010 for periodic conditions: detect common points between left and right edges

!***********************************************************************
!
!             i n i t i a l i z a t i o n    p h a s e
!
!***********************************************************************
  call initialize_simulation(nproc,myrank,NUMBER_OF_PASSES, &
                  ninterface_acoustic,ninterface_elastic,ninterface_poroelastic)


  ! reduction of cache misses inner/outer in two passes
  do ipass = 1,NUMBER_OF_PASSES

  ! starts reading in Database file
  call read_databases_init(myrank,ipass, &
                  simulation_title,SIMULATION_TYPE,SAVE_FORWARD,npgeo, &
                  gnuplot,interpol,NTSTEP_BETWEEN_OUTPUT_INFO, &
                  output_postscript_snapshot,output_color_image,colors,numbers, &
                  meshvect,modelvect,boundvect,cutsnaps,subsamp,sizemax_arrows, &
                  anglerec,initialfield,add_Bielak_conditions, &
                  seismotype,imagetype,assign_external_model,READ_EXTERNAL_SEP_FILE, &
                  outputgrid,OUTPUT_ENERGY,TURN_ATTENUATION_ON, &
                  TURN_VISCATTENUATION_ON,Q0,freq0,p_sv, &
                  NSTEP,deltat,NTSTEP_BETWEEN_OUTPUT_SEISMO,NSOURCES)

  !
  !--- source information
  !
  if(ipass == 1) then
    allocate( source_type(NSOURCES) )
    allocate( time_function_type(NSOURCES) )
    allocate( x_source(NSOURCES) )
    allocate( z_source(NSOURCES) )
    allocate( f0(NSOURCES) )
    allocate( tshift_src(NSOURCES) )
    allocate( factor(NSOURCES) )
    allocate( angleforce(NSOURCES) )
    allocate( Mxx(NSOURCES) )
    allocate( Mxz(NSOURCES) )
    allocate( Mzz(NSOURCES) )
    allocate( aval(NSOURCES) )
    allocate( ispec_selected_source(NSOURCES) )
    allocate( iglob_source(NSOURCES) )
    allocate( source_courbe_eros(NSOURCES) )
    allocate( xi_source(NSOURCES) )
    allocate( gamma_source(NSOURCES) )
    allocate( is_proc_source(NSOURCES) )
    allocate( nb_proc_source(NSOURCES) )
    allocate( sourcearray(NSOURCES,NDIM,NGLLX,NGLLZ) )
  endif

  ! reads in source infos
  call read_databases_sources(NSOURCES,source_type,time_function_type, &
                      x_source,z_source,Mxx,Mzz,Mxz,f0,tshift_src,factor,angleforce)

  ! sets source parameters
  call set_sources(myrank,NSOURCES,source_type,time_function_type, &
                      x_source,z_source,Mxx,Mzz,Mxz,f0,tshift_src,factor,angleforce,aval, &
                      t0,initialfield,ipass,deltat)

  !
  !----  read attenuation information
  !
  call read_databases_atten(N_SLS,f0_attenuation)
  
  ! if source is not a Dirac or Heavyside then f0_attenuation is f0 of the first source
  if(.not. (time_function_type(1) == 4 .or. time_function_type(1) == 5)) then
    f0_attenuation = f0(1)
  endif


  !
  !---- read the spectral macrobloc nodal coordinates
  !
  if(ipass == 1) allocate(coorg(NDIM,npgeo))

  ! reads the spectral macrobloc nodal coordinates
  ! and basic properties of the spectral elements
  call read_databases_coorg_elem(myrank,ipass,npgeo,coorg,numat,ngnod,nspec, &
                              pointsdisp,plot_lowerleft_corner_only, &
                              nelemabs,nelem_acoustic_surface, &
                              num_fluid_solid_edges,num_fluid_poro_edges, &
                              num_solid_poro_edges,nnodes_tangential_curve)


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
    allocate(anisotropy(6,numat))
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
    allocate(anisotropic(nspec))
    allocate(inv_tau_sigma_nu1(NGLLX,NGLLZ,nspec,N_SLS))
    allocate(inv_tau_sigma_nu2(NGLLX,NGLLZ,nspec,N_SLS))
    allocate(phi_nu1(NGLLX,NGLLZ,nspec,N_SLS))
    allocate(phi_nu2(NGLLX,NGLLZ,nspec,N_SLS))
    allocate(inv_tau_sigma_nu1_sent(N_SLS))
    allocate(inv_tau_sigma_nu2_sent(N_SLS))
    allocate(phi_nu1_sent(N_SLS))
    allocate(phi_nu2_sent(N_SLS))
  endif

  !
  !---- read the material properties
  !
  call gmat01(density,porosity,tortuosity,anisotropy,permeability,poroelastcoef,numat,&
              myrank,ipass,Qp_attenuation,Qs_attenuation,freq0,Q0,f0(1),TURN_VISCATTENUATION_ON)
  !
  !----  read spectral macrobloc data
  !
  if(ipass == 1) then
    allocate(antecedent_list(nspec))
    allocate(perm(nspec))    
  endif  
  call read_databases_mato(ipass,nspec,ngnod,kmato,knods, &
                                perm,antecedent_list)
  

!-------------------------------------------------------------------------------
!----  determine if each spectral element is elastic, poroelastic, or acoustic
!-------------------------------------------------------------------------------
  ! initializes
  any_acoustic = .false.
  any_elastic = .false.
  any_poroelastic = .false.
  
  anisotropic(:) = .false.
  elastic(:) = .false.
  poroelastic(:) = .false.

  ! loops over all elements
  do ispec = 1,nspec

    if( nint(porosity(kmato(ispec))) == 1 ) then  
      ! acoustic domain
      elastic(ispec) = .false.
      poroelastic(ispec) = .false.
      any_acoustic = .true.
    elseif( porosity(kmato(ispec)) < TINYVAL) then  
      ! elastic domain
      elastic(ispec) = .true.
      poroelastic(ispec) = .false.
      any_elastic = .true.
      if(any(anisotropy(:,kmato(ispec)) /= 0)) then
         anisotropic(ispec) = .true.
      end if
    else                                       
      ! poroelastic domain
      elastic(ispec) = .false.
      poroelastic(ispec) = .true.
      any_poroelastic = .true.
    endif

  enddo !do ispec = 1,nspec


  if(.not. p_sv .and. .not. any_elastic) then
    print*, '*************** WARNING ***************'
    print*, 'Surface (membrane) waves calculation needs an elastic medium'
    print*, '*************** WARNING ***************'
    stop
  endif
  if(.not. p_sv .and. (TURN_ATTENUATION_ON)) then
    print*, '*************** WARNING ***************'
    print*, 'Attenuation and anisotropy are not implemented for surface (membrane) waves calculation'
    print*, '*************** WARNING ***************'
    stop
  endif


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
            f0_attenuation,inv_tau_sigma_nu1_sent,phi_nu1_sent, &
            inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu1_sent,Mu_nu2_sent)
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
  call read_databases_ninterface(ninterface,max_interface_size)    
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
   call read_databases_interfaces(ipass,ninterface,nspec,max_interface_size, &
                              my_neighbours,my_nelmnts_neighbours,my_interfaces, &                              
                              perm,antecedent_list)  

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

    allocate(ib_left(nelemabs))
    allocate(ib_right(nelemabs))
    allocate(ib_bottom(nelemabs))
    allocate(ib_top(nelemabs))
    
  endif

  !
  !----  read absorbing boundary data
  !
  call read_databases_absorbing(myrank,ipass,nelemabs,nspec,anyabs, &
                            ibegin_bottom,iend_bottom,jbegin_right,jend_right, &
                            ibegin_top,iend_top,jbegin_left,jend_left, &
                            numabs,codeabs,perm,antecedent_list, &
                            nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax, &
                            ib_right,ib_left,ib_bottom,ib_top)

  
  if( anyabs ) then

!    nspec_xmin = 0
!    nspec_xmax = 0
!    nspec_zmin = 0
!    nspec_zmax = 0
!    if(ipass == 1) then
!      allocate(ib_left(nelemabs))
!      allocate(ib_right(nelemabs))
!      allocate(ib_bottom(nelemabs))
!      allocate(ib_top(nelemabs))
!    endif
    
!    do inum = 1,nelemabs
!      if (codeabs(IBOTTOM,inum)) then
!        nspec_zmin = nspec_zmin + 1
!        ib_bottom(inum) =  nspec_zmin
!      endif
!      if (codeabs(IRIGHT,inum)) then
!        nspec_xmax = nspec_xmax + 1
!        ib_right(inum) =  nspec_xmax
!      endif
!      if (codeabs(ITOP,inum)) then
!        nspec_zmax = nspec_zmax + 1
!        ib_top(inum) = nspec_zmax
!      endif
!      if (codeabs(ILEFT,inum)) then
!        nspec_xmin = nspec_xmin + 1
!        ib_left(inum) =  nspec_xmin
!      endif
!    enddo

! Files to save absorbed waves needed to reconstruct backward wavefield for adjoint method
    if(ipass == 1) then
      if(any_elastic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 2)) then
        allocate(b_absorb_elastic_left(3,NGLLZ,nspec_xmin,NSTEP))
        allocate(b_absorb_elastic_right(3,NGLLZ,nspec_xmax,NSTEP))
        allocate(b_absorb_elastic_bottom(3,NGLLX,nspec_zmin,NSTEP))
        allocate(b_absorb_elastic_top(3,NGLLX,nspec_zmax,NSTEP))
      else
        allocate(b_absorb_elastic_left(1,1,1,1))
        allocate(b_absorb_elastic_right(1,1,1,1))
        allocate(b_absorb_elastic_bottom(1,1,1,1))
        allocate(b_absorb_elastic_top(1,1,1,1))
      endif
      if(any_poroelastic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 2)) then
        allocate(b_absorb_poro_s_left(NDIM,NGLLZ,nspec_xmin,NSTEP))
        allocate(b_absorb_poro_s_right(NDIM,NGLLZ,nspec_xmax,NSTEP))
        allocate(b_absorb_poro_s_bottom(NDIM,NGLLX,nspec_zmin,NSTEP))
        allocate(b_absorb_poro_s_top(NDIM,NGLLX,nspec_zmax,NSTEP))
        allocate(b_absorb_poro_w_left(NDIM,NGLLZ,nspec_xmin,NSTEP))
        allocate(b_absorb_poro_w_right(NDIM,NGLLZ,nspec_xmax,NSTEP))
        allocate(b_absorb_poro_w_bottom(NDIM,NGLLX,nspec_zmin,NSTEP))
        allocate(b_absorb_poro_w_top(NDIM,NGLLX,nspec_zmax,NSTEP))
      else
        allocate(b_absorb_poro_s_left(1,1,1,1))
        allocate(b_absorb_poro_s_right(1,1,1,1))
        allocate(b_absorb_poro_s_bottom(1,1,1,1))
        allocate(b_absorb_poro_s_top(1,1,1,1))
        allocate(b_absorb_poro_w_left(1,1,1,1))
        allocate(b_absorb_poro_w_right(1,1,1,1))
        allocate(b_absorb_poro_w_bottom(1,1,1,1))
        allocate(b_absorb_poro_w_top(1,1,1,1))
      endif
      if(any_acoustic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 2)) then
        allocate(b_absorb_acoustic_left(NGLLZ,nspec_xmin,NSTEP))
        allocate(b_absorb_acoustic_right(NGLLZ,nspec_xmax,NSTEP))
        allocate(b_absorb_acoustic_bottom(NGLLX,nspec_zmin,NSTEP))
        allocate(b_absorb_acoustic_top(NGLLX,nspec_zmax,NSTEP))
      else
        allocate(b_absorb_acoustic_left(1,1,1))
        allocate(b_absorb_acoustic_right(1,1,1))
        allocate(b_absorb_acoustic_bottom(1,1,1))
        allocate(b_absorb_acoustic_top(1,1,1))
      endif
    endif

!    if (myrank == 0 ) then
!      write(IOUT,*)
!      write(IOUT,*) 'nspec_xmin = ',nspec_xmin
!      write(IOUT,*) 'nspec_xmax = ',nspec_xmax
!      write(IOUT,*) 'nspec_zmin = ',nspec_zmin
!      write(IOUT,*) 'nspec_zmax = ',nspec_zmax
!    endif
    
  else

!    if(.not. allocated(ib_left)) then
!      allocate(ib_left(1))
!      allocate(ib_right(1))
!      allocate(ib_bottom(1))
!      allocate(ib_top(1))
!    endif

    if(.not. allocated(b_absorb_elastic_left)) then
      allocate(b_absorb_elastic_left(1,1,1,1))
      allocate(b_absorb_elastic_right(1,1,1,1))
      allocate(b_absorb_elastic_bottom(1,1,1,1))
      allocate(b_absorb_elastic_top(1,1,1,1))
    endif

    if(.not. allocated(b_absorb_poro_s_left)) then
      allocate(b_absorb_poro_s_left(1,1,1,1))
      allocate(b_absorb_poro_s_right(1,1,1,1))
      allocate(b_absorb_poro_s_bottom(1,1,1,1))
      allocate(b_absorb_poro_s_top(1,1,1,1))
      allocate(b_absorb_poro_w_left(1,1,1,1))
      allocate(b_absorb_poro_w_right(1,1,1,1))
      allocate(b_absorb_poro_w_bottom(1,1,1,1))
      allocate(b_absorb_poro_w_top(1,1,1,1))
    endif

    if(.not. allocated(b_absorb_acoustic_left)) then
      allocate(b_absorb_acoustic_left(1,1,1))
      allocate(b_absorb_acoustic_right(1,1,1))
      allocate(b_absorb_acoustic_bottom(1,1,1))
      allocate(b_absorb_acoustic_top(1,1,1))
    endif

  endif

!
!----  read acoustic free surface data
!
  if(nelem_acoustic_surface > 0) then
    any_acoustic_edges = .true.
  else
    any_acoustic_edges = .false.
    nelem_acoustic_surface = 1
  endif
  if( ipass == 1 ) then
    allocate(acoustic_edges(4,nelem_acoustic_surface))
    allocate(acoustic_surface(5,nelem_acoustic_surface))
  endif
  call read_databases_free_surf(ipass,nelem_acoustic_surface,nspec, &
                            acoustic_edges,perm,antecedent_list,any_acoustic_edges)
  ! resets nelem_acoustic_surface
  if( any_acoustic_edges .eqv. .false. ) nelem_acoustic_surface = 0

  ! constructs acoustic surface
  if(nelem_acoustic_surface > 0) then
    call construct_acoustic_surface ( nspec, ngnod, knods, nelem_acoustic_surface, &
                                     acoustic_edges, acoustic_surface)        
    if (myrank == 0 .and. ipass == 1) then
      write(IOUT,*)
      write(IOUT,*) 'Number of free surface elements: ',nelem_acoustic_surface
    endif
  endif


  !
  !---- read coupled edges
  !
  if( num_fluid_solid_edges > 0 ) then
    any_fluid_solid_edges = .true.
  else
    any_fluid_solid_edges = .false.
    num_fluid_solid_edges = 1
  endif
  if(ipass == 1) then
    allocate(fluid_solid_acoustic_ispec(num_fluid_solid_edges))
    allocate(fluid_solid_acoustic_iedge(num_fluid_solid_edges))
    allocate(fluid_solid_elastic_ispec(num_fluid_solid_edges))
    allocate(fluid_solid_elastic_iedge(num_fluid_solid_edges))
  endif
  if( num_fluid_poro_edges > 0 ) then  
    any_fluid_poro_edges = .true.
  else
    any_fluid_poro_edges = .false.
    num_fluid_poro_edges = 1
  endif
  if(ipass == 1) then
    allocate(fluid_poro_acoustic_ispec(num_fluid_poro_edges))
    allocate(fluid_poro_acoustic_iedge(num_fluid_poro_edges))
    allocate(fluid_poro_poroelastic_ispec(num_fluid_poro_edges))
    allocate(fluid_poro_poroelastic_iedge(num_fluid_poro_edges))
  endif
  if ( num_solid_poro_edges > 0 ) then  
    any_solid_poro_edges = .true.
  else
    any_solid_poro_edges = .false.
    num_solid_poro_edges = 1
  endif
  if(ipass == 1) then
    allocate(solid_poro_elastic_ispec(num_solid_poro_edges))
    allocate(solid_poro_elastic_iedge(num_solid_poro_edges))
    allocate(solid_poro_poroelastic_ispec(num_solid_poro_edges))
    allocate(solid_poro_poroelastic_iedge(num_solid_poro_edges))
  endif
  
  call read_databases_coupled(ipass,nspec,num_fluid_solid_edges,any_fluid_solid_edges, &
                            fluid_solid_acoustic_ispec,fluid_solid_elastic_ispec, &
                            num_fluid_poro_edges,any_fluid_poro_edges, &
                            fluid_poro_acoustic_ispec,fluid_poro_poroelastic_ispec, &
                            num_solid_poro_edges,any_solid_poro_edges, &
                            solid_poro_elastic_ispec,solid_poro_poroelastic_ispec, &
                            perm,antecedent_list)  
  
  ! resets counters
  if( any_fluid_solid_edges .eqv. .false. ) num_fluid_solid_edges = 0
  if( any_fluid_poro_edges .eqv. .false. ) num_fluid_poro_edges = 0
  if( any_solid_poro_edges .eqv. .false. ) num_solid_poro_edges = 0

  
  !
  !---- read tangential detection curve
  !      and close Database file
  !
  if (nnodes_tangential_curve > 0) then
    any_tangential_curve = .true.
  else
    any_tangential_curve = .false.
    nnodes_tangential_curve = 1
  endif  
  if (ipass == 1) then
    allocate(nodes_tangential_curve(2,nnodes_tangential_curve))
    allocate(dist_tangential_detection_curve(nnodes_tangential_curve))
  endif
  call read_databases_final(nnodes_tangential_curve,nodes_tangential_curve, &
                                force_normal_to_surface,rec_normal_to_surface, &
                                any_tangential_curve)                                
  ! resets nnode_tangential_curve
  if( any_tangential_curve .eqv. .false. ) nnodes_tangential_curve = 0
  
!
!---- compute shape functions and their derivatives for SEM grid
!

! set up Gauss-Lobatto-Legendre derivation matrices
  call define_derivation_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,hprimewgll_xx,hprimewgll_zz)

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

    print *
    print *,'Xmin,Xmax of the whole mesh = ',minval(coord(1,:)),maxval(coord(1,:))
    print *,'Zmin,Zmax of the whole mesh = ',minval(coord(2,:)),maxval(coord(2,:))
    print *

!! DK DK Feb 2010 for periodic conditions: detect common points between left and right edges

    if(ADD_PERIODIC_CONDITIONS) then

#ifdef USE_MPI
  stop 'periodic conditions currently implemented for a serial simulation only (due e.g. to mass matrix rebuilding)'
#endif

  if(any_poroelastic .or. any_acoustic) stop 'periodic conditions currently implemented for purely elastic models only'

  if(ACTUALLY_IMPLEMENT_PERM_OUT .or. ACTUALLY_IMPLEMENT_PERM_INN .or. ACTUALLY_IMPLEMENT_PERM_WHOLE) &
    stop 'currently, all permutations should be off for periodic conditions'

print *
open(unit=123,file='Database00000_left_edge_only',status='old')
do ispecperio = 1,NSPEC_PERIO
  read(123,*) numperio_left(ispecperio), &
     codeabs_perio_left(IBOTTOM,ispecperio), &
     codeabs_perio_left(IRIGHT,ispecperio), &
     codeabs_perio_left(ITOP,ispecperio), &
     codeabs_perio_left(ILEFT,ispecperio), &
     idummy1, idummy2, idummy3, idummy4, idummy5, idummy6, idummy7, idummy8
enddo
close(123)
print *,'read ',NSPEC_PERIO,' elements for left periodic edge'

open(unit=123,file='Database00000_right_edge_only',status='old')
do ispecperio = 1,NSPEC_PERIO
  read(123,*) numperio_right(ispecperio), &
     codeabs_perio_right(IBOTTOM,ispecperio), &
     codeabs_perio_right(IRIGHT,ispecperio), &
     codeabs_perio_right(ITOP,ispecperio), &
     codeabs_perio_right(ILEFT,ispecperio), &
     idummy1, idummy2, idummy3, idummy4, idummy5, idummy6, idummy7, idummy8
enddo
close(123)
print *,'read ',NSPEC_PERIO,' elements for right periodic edge'
print *

print *,'because of periodic conditions, values computed by checkgrid() are not reliable'
print *

!---------------------------------------------------------------------------

         do ispecperio = 1,NSPEC_PERIO

            ispec = numperio_left(ispecperio)

! print *,'dist of edge is ',sqrt((coord(2,ibool(1,1,ispec)) - coord(2,ibool(1,NGLLZ,ispec))) ** 2 + &
!                                 (coord(1,ibool(1,1,ispec)) - coord(1,ibool(1,NGLLZ,ispec))) ** 2)

            if(codeabs_perio_left(ILEFT,ispecperio)) then
               i = 1
               do j = 1,NGLLZ
                  iglob = ibool(i,j,ispec)
!----------------------------------------------------------------------
                  include "include_for_periodic_conditions.f90"
!----------------------------------------------------------------------
               enddo
            endif

            if(codeabs_perio_left(IRIGHT,ispecperio)) then
               i = NGLLX
               do j = 1,NGLLZ
                  iglob = ibool(i,j,ispec)
!----------------------------------------------------------------------
                  include "include_for_periodic_conditions.f90"
!----------------------------------------------------------------------
               enddo
            endif

            if(codeabs_perio_left(IBOTTOM,ispecperio)) then
               j = 1
               do i = 1,NGLLX
                  iglob = ibool(i,j,ispec)
!----------------------------------------------------------------------
                  include "include_for_periodic_conditions.f90"
!----------------------------------------------------------------------
               enddo
            endif

            if(codeabs_perio_left(ITOP,ispecperio)) then
               j = NGLLZ
               do i = 1,NGLLX
                  iglob = ibool(i,j,ispec)
!----------------------------------------------------------------------
                  include "include_for_periodic_conditions.f90"
!----------------------------------------------------------------------
               enddo
            endif

         enddo

! rebuild the mass matrix based on this new numbering
!
!---- build the global mass matrix and invert it once and for all
!
      rmass_inverse_elastic(:) = ZERO
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
              kappal = lambdal_relaxed + 2.d0/3.d0*mul_relaxed
            endif

             rmass_inverse_elastic(iglob) = rmass_inverse_elastic(iglob) &
                                + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)

          enddo
        enddo
      enddo ! do ispec = 1,nspec

! invert the mass matrix once and for all
! set entries that are equal to zero to something else, e.g. 1, to avoid division by zero
! these degrees of freedom correspond to points that have been replaced with their periodic counterpart
! and thus are not used any more
      where(rmass_inverse_elastic == ZERO) rmass_inverse_elastic = 1._CUSTOM_REAL
      rmass_inverse_elastic(:) = 1._CUSTOM_REAL / rmass_inverse_elastic(:)

    endif ! of if(ADD_PERIODIC_CONDITIONS)

!! DK DK Feb 2010 for periodic conditions: detect common points between left and right edges

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

  endif ! ipass

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
  open(unit=IIN,file='DATA/STATIONS_target',iostat=ios,status='old',action='read')
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
    allocate(hxis(NGLLX))
    allocate(hpxir(NGLLX))
    allocate(hpxis(NGLLX))
    allocate(hgammar(NGLLZ))
    allocate(hgammas(NGLLZ))
    allocate(hpgammar(NGLLZ))
    allocate(hpgammas(NGLLZ))

! allocate Lagrange interpolators for receivers
    allocate(hxir_store(nrec,NGLLX))
    allocate(hgammar_store(nrec,NGLLZ))

! allocate Lagrange interpolators for sources
    allocate(hxis_store(NSOURCES,NGLLX))
    allocate(hgammas_store(NSOURCES,NGLLZ))

! allocate other global arrays
    allocate(coord(NDIM,npoin))

! to display acoustic elements
    allocate(vector_field_display(3,npoin))

!    if(assign_external_model) then

! note: so far, full external array needed/defined in subroutine calls
      allocate(vpext(NGLLX,NGLLZ,nspec))
      allocate(vsext(NGLLX,NGLLZ,nspec))
      allocate(rhoext(NGLLX,NGLLZ,nspec))
      allocate(Qp_attenuationext(NGLLX,NGLLZ,nspec))
      allocate(Qs_attenuationext(NGLLX,NGLLZ,nspec))
      allocate(c11ext(NGLLX,NGLLZ,nspec))
      allocate(c13ext(NGLLX,NGLLZ,nspec))
      allocate(c15ext(NGLLX,NGLLZ,nspec))
      allocate(c33ext(NGLLX,NGLLZ,nspec))
      allocate(c35ext(NGLLX,NGLLZ,nspec))
      allocate(c55ext(NGLLX,NGLLZ,nspec))
!    else
!      allocate(vpext(1,1,1))
!      allocate(vsext(1,1,1))
!      allocate(rhoext(1,1,1))
!      allocate(c11ext(1,1,1))
!      allocate(c13ext(1,1,1))
!      allocate(c15ext(1,1,1))
!      allocate(c33ext(1,1,1))
!      allocate(c35ext(1,1,1))
!      allocate(c55ext(1,1,1))
!    endif

  endif

!
!----  set the coordinates of the points of the global grid
!
  found_a_negative_jacobian = .false.
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX

        xi = xigll(i)
        gamma = zigll(j)

        call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl, &
                        jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo, &
                        .false.)

        if(jacobianl <= ZERO) found_a_negative_jacobian = .true.

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

! create an OpenDX file containing all the negative elements displayed in red, if any
! this allows users to locate problems in a mesh based on the OpenDX file created at the second iteration
! do not create OpenDX files if no negative Jacobian has been found, or if we are running in parallel
! (because writing OpenDX routines is much easier in serial)
  if(found_a_negative_jacobian .and. nproc == 1) then
    call save_openDX_jacobian(nspec,npgeo,ngnod,knods,coorg,xigll,zigll)
  endif

! stop the code at the first negative element found, because such a mesh cannot be computed
  if(found_a_negative_jacobian) then

    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX

          xi = xigll(i)
          gamma = zigll(j)

          call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl, &
                          jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo, &
                          .true.)

        enddo
      enddo
    enddo

  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! yang  output weights for line, surface integrals !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!define_derivation_matrices(xigll(NGLLX),zigll(NGLLZ),wxgll(NGLLX),wzgll(NGLLZ),hprime_xx(NGLLX,NGLLX),hprime_zz(NGLLZ,NGLLZ),&
!                           hprimewgll_xx(NGLLX,NGLLX),hprimewgll_zz(NGLLZ,NGLLZ))
!xix(NGLLX,NGLLZ,nspec),xiz,gammax,gammaz,jacobian
!recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl,jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo, &
!          .true.)
  allocate(weight_line_x(npoin))
  allocate(weight_line_z(npoin))
  allocate(weight_surface(npoin))
  allocate(weight_jacobian(npoin))
  allocate(weight_gll(npoin))
  weight_line_x=0.0
  weight_line_z=0.0
  weight_surface=0.0
  zmin_yang=minval(coord(2,:))
  xmin_yang=minval(coord(1,:))
  zmax_yang=maxval(coord(2,:))
  xmax_yang=maxval(coord(1,:))
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
            iglob=ibool(i,j,ispec)
            z=coord(2,ibool(i,j,ispec))
            xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
            zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
            if ((j==1 .OR. j==NGLLZ) .AND. ( (abs(z-zmin_yang).GE.1) .AND. (abs(z-zmax_yang)).GE.1) )    xxi=xxi/2.0
            if ((i==1 .OR. i==NGLLZ) .AND. ( (abs(x-xmin_yang).GE.1) .AND. (abs(x-xmax_yang)).GE.1) )    zgamma=zgamma/2.0
            weight_line_x(iglob) =  weight_line_x(iglob) + xxi * wxgll(i)
            weight_line_z(iglob) =  weight_line_z(iglob) + zgamma * wzgll(j)
            weight_surface(iglob) = weight_surface(iglob) + wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
            weight_jacobian(iglob) = jacobian(i,j,ispec)
            weight_gll(iglob) = 10*j+i
      enddo
    enddo
  enddo
  open(unit=55,file='OUTPUT_FILES/x_z_weightLineX_weightLineZ_weightSurface',status='unknown')
  do n = 1,npoin
    write(55,*) coord(1,n), coord(2,n), weight_line_x(n), weight_line_z(n), weight_surface(n),weight_jacobian(n),weight_gll(n)
  enddo
  close(55)
  deallocate(weight_line_x)
  deallocate(weight_line_z)
  deallocate(weight_surface)
  deallocate(weight_jacobian)
  deallocate(weight_gll)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!--- save the grid of points in a file
!
  if(outputgrid .and. myrank == 0 .and. ipass == 1) then
     write(IOUT,*)
     write(IOUT,*) 'Saving the grid in a text file...'
     write(IOUT,*)
     open(unit=55,file='OUTPUT_FILES/grid_points_and_model.txt',status='unknown')
     write(55,*) npoin
     do n = 1,npoin
        write(55,*) (coord(i,n), i=1,NDIM)
     enddo
     close(55)
  endif

!
!-----   plot the GLL mesh in a Gnuplot file
!
  if(gnuplot .and. myrank == 0 .and. ipass == 1)  &
    call plotgll(knods,ibool,coorg,coord,npoin,npgeo,ngnod,nspec)

  if(myrank == 0 .and. ipass == 1)  &
    write(IOUT,*) 'assign_external_model = ', assign_external_model

!if ( assign_external_model .and. ipass == 1) then
  if ( assign_external_model) then
    call read_external_model(any_acoustic,any_elastic,any_poroelastic, &
                elastic,poroelastic,anisotropic,nspec,npoin,N_SLS,ibool, &
                f0_attenuation,inv_tau_sigma_nu1_sent,phi_nu1_sent, &
                inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu1_sent,Mu_nu2_sent, &
                inv_tau_sigma_nu1,inv_tau_sigma_nu2,phi_nu1,phi_nu2,Mu_nu1,Mu_nu2,&
                coord,kmato,myrank,rhoext,vpext,vsext, &
                Qp_attenuationext,Qs_attenuationext, &
                c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,READ_EXTERNAL_SEP_FILE)
  end if

!
!----  perform basic checks on parameters read
!
  all_anisotropic = .false.
  if(count(anisotropic(:) .eqv. .true.) == nspec) all_anisotropic = .true.
  
  if(all_anisotropic .and. anyabs) &
    call exit_MPI('Cannot put absorbing boundaries if anisotropic materials along edges')
    
  if(TURN_ATTENUATION_ON .and. all_anisotropic) then
    call exit_MPI('Cannot turn attenuation on in anisotropic materials')
  end if

  ! global domain flags
  any_elastic_glob = any_elastic
#ifdef USE_MPI
  call MPI_ALLREDUCE(any_elastic, any_elastic_glob, 1, MPI_LOGICAL, &
                    MPI_LOR, MPI_COMM_WORLD, ier)
#endif

  any_poroelastic_glob = any_poroelastic
#ifdef USE_MPI
  call MPI_ALLREDUCE(any_poroelastic, any_poroelastic_glob, 1, MPI_LOGICAL, &
                    MPI_LOR, MPI_COMM_WORLD, ier)
#endif

  any_acoustic_glob = any_acoustic
#ifdef USE_MPI
  call MPI_ALLREDUCE(any_acoustic, any_acoustic_glob, 1, MPI_LOGICAL, &
                    MPI_LOR, MPI_COMM_WORLD, ier)
#endif

  ! for acoustic
  if(TURN_ATTENUATION_ON .and. .not. any_elastic_glob) &
    call exit_MPI('currently cannot have attenuation if acoustic/poroelastic simulation only')

!
!----   define coefficients of the Newmark time scheme
!
  deltatover2 = HALF*deltat
  deltatsquareover2 = HALF*deltat*deltat

  if(SIMULATION_TYPE == 2) then
!  define coefficients of the Newmark time scheme for the backward wavefield
    b_deltat = - deltat
    b_deltatover2 = HALF*b_deltat
    b_deltatsquareover2 = HALF*b_deltat*b_deltat
  endif

!---- define actual location of source and receivers

  call setup_sources_receivers(NSOURCES,initialfield,source_type,&
     coord,ibool,npoin,nspec,nelem_acoustic_surface,acoustic_surface,elastic,poroelastic, &
     x_source,z_source,ispec_selected_source,ispec_selected_rec, &
     is_proc_source,nb_proc_source,ipass,&
     sourcearray,Mxx,Mzz,Mxz,xix,xiz,gammax,gammaz,xigll,zigll,npgeo,&
     nproc,myrank,xi_source,gamma_source,coorg,knods,ngnod, &
     nrec,nrecloc,recloc,which_proc_receiver,st_xval,st_zval, &
     xi_receiver,gamma_receiver,station_name,network_name,x_final_receiver,z_final_receiver,iglob_source)

! compute source array for adjoint source
  if(SIMULATION_TYPE == 2) then  ! adjoint calculation
    nadj_rec_local = 0
    do irec = 1,nrec
      if(myrank == which_proc_receiver(irec))then
!   check that the source proc number is okay
        if(which_proc_receiver(irec) < 0 .or. which_proc_receiver(irec) > NPROC-1) &
              call exit_MPI('something is wrong with the source proc number in adjoint simulation')
        nadj_rec_local = nadj_rec_local + 1
      endif
    enddo
    if(ipass == 1) allocate(adj_sourcearray(NSTEP,3,NGLLX,NGLLZ))
    if (nadj_rec_local > 0 .and. ipass == 1)  then
      allocate(adj_sourcearrays(nadj_rec_local,NSTEP,3,NGLLX,NGLLZ))
    else if (ipass == 1) then
      allocate(adj_sourcearrays(1,1,1,1,1))
    endif

    irec_local = 0
    do irec = 1, nrec
!   compute only adjoint source arrays in the local proc
      if(myrank == which_proc_receiver(irec))then
        irec_local = irec_local + 1
        adj_source_file = trim(station_name(irec))//'.'//trim(network_name(irec))
        call compute_arrays_adj_source(adj_source_file, &
                            xi_receiver(irec), gamma_receiver(irec), &
                            adj_sourcearray, xigll,zigll,NSTEP)
        adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
      endif
    enddo
  else if (ipass == 1) then
     allocate(adj_sourcearrays(1,1,1,1,1))
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

    if (rec_normal_to_surface .and. abs(anglerec) > 1.d-6) &
      stop 'anglerec should be zero when receivers are normal to the topography'

    anglerec_irec(:) = anglerec * pi / 180.d0
    cosrot_irec(:) = cos(anglerec_irec(:))
    sinrot_irec(:) = sin(anglerec_irec(:))
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
         call tri_quad(n_tangential_detection_curve, n1_tangential_detection_curve, &
                      nnodes_tangential_curve)

         call compute_normal_vector( anglerec_irec(irecloc), &
           nodes_tangential_curve(1,n_tangential_detection_curve(1)), &
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

      do i_source=1,NSOURCES
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

          call tri_quad(n_tangential_detection_curve, n1_tangential_detection_curve, &
                       nnodes_tangential_curve)

          ! in the case of a source force vector
          ! users can give an angle with respect to the normal to the topography surface,
          ! in which case we must compute the normal to the topography
          ! and add it the existing rotation angle
          call compute_normal_vector( angleforce(i_source), &
                            nodes_tangential_curve(1,n_tangential_detection_curve(1)), &
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
              call MPI_recv(source_courbe_eros(i_source),1,MPI_INTEGER, &
                          MPI_ANY_SOURCE,42,MPI_COMM_WORLD,request_mpi_status,ier)
              call MPI_recv(angleforce_recv,1,MPI_DOUBLE_PRECISION, &
                          MPI_ANY_SOURCE,43,MPI_COMM_WORLD,request_mpi_status,ier)
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
      enddo ! do i_source=1,NSOURCES
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
        open(unit=11,file='OUTPUT_FILES/dist_rec_tangential_detection_curve', &
              form='formatted', status='unknown')
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

    endif ! force_normal_to_surface
    
  endif ! ipass

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
            call exit_MPI('an acoustic pressure receiver cannot be located exactly '// &
                            'on the free surface because pressure is zero there')
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

! define and store Lagrange interpolators at all the sources
  do i = 1,NSOURCES
    call lagrange_any(xi_source(i),NGLLX,xigll,hxis,hpxis)
    call lagrange_any(gamma_source(i),NGLLZ,zigll,hgammas,hpgammas)
    hxis_store(i,:) = hxis(:)
    hgammas_store(i,:) = hgammas(:)
  enddo

! displacement, velocity, acceleration and inverse of the mass matrix for elastic elements
  if(ipass == 1) then

    if(any_elastic) then
      allocate(displ_elastic(3,npoin))
      allocate(veloc_elastic(3,npoin))
      allocate(accel_elastic(3,npoin))
      allocate(rmass_inverse_elastic(npoin))
    else
    ! allocate unused arrays with fictitious size
      allocate(displ_elastic(1,1))
      allocate(veloc_elastic(1,1))
      allocate(accel_elastic(1,1))
      allocate(rmass_inverse_elastic(1))
    endif
    ! extra array if adjoint and kernels calculation
    if(SIMULATION_TYPE == 2 .and. any_elastic) then
      allocate(b_displ_elastic(3,npoin))
      allocate(b_veloc_elastic(3,npoin))
      allocate(b_accel_elastic(3,npoin))
      allocate(rho_kl(NGLLX,NGLLZ,nspec))
      allocate(rho_k(npoin))
      allocate(rhol_global(npoin))
      allocate(mu_kl(NGLLX,NGLLZ,nspec))
      allocate(mu_k(npoin))
      allocate(mul_global(npoin))
      allocate(kappa_kl(NGLLX,NGLLZ,nspec))
      allocate(kappa_k(npoin))
      allocate(kappal_global(npoin))
      allocate(rhop_kl(NGLLX,NGLLZ,nspec))
      allocate(alpha_kl(NGLLX,NGLLZ,nspec))
      allocate(beta_kl(NGLLX,NGLLZ,nspec))
      allocate(rhorho_el_hessian_final2(NGLLX,NGLLZ,nspec))
      allocate(rhorho_el_hessian_temp2(npoin))
      allocate(rhorho_el_hessian_final1(NGLLX,NGLLZ,nspec))
      allocate(rhorho_el_hessian_temp1(npoin))
    else
      allocate(b_displ_elastic(1,1))
      allocate(b_veloc_elastic(1,1))
      allocate(b_accel_elastic(1,1))
      allocate(rho_kl(1,1,1))
      allocate(rho_k(1))
      allocate(rhol_global(1))
      allocate(mu_kl(1,1,1))
      allocate(mu_k(1))
      allocate(mul_global(1))
      allocate(kappa_kl(1,1,1))
      allocate(kappa_k(1))
      allocate(kappal_global(1))
      allocate(rhop_kl(1,1,1))
      allocate(alpha_kl(1,1,1))
      allocate(beta_kl(1,1,1))
      allocate(rhorho_el_hessian_final2(1,1,1))
      allocate(rhorho_el_hessian_temp2(1))
      allocate(rhorho_el_hessian_final1(1,1,1))
      allocate(rhorho_el_hessian_temp1(1))
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
    if(SIMULATION_TYPE == 2 .and. any_poroelastic) then
      allocate(b_displs_poroelastic(NDIM,npoin))
      allocate(b_velocs_poroelastic(NDIM,npoin))
      allocate(b_accels_poroelastic(NDIM,npoin))
      allocate(b_displw_poroelastic(NDIM,npoin))
      allocate(b_velocw_poroelastic(NDIM,npoin))
      allocate(b_accelw_poroelastic(NDIM,npoin))
      allocate(rhot_kl(NGLLX,NGLLZ,nspec))
      allocate(rhot_k(npoin))
      allocate(rhof_kl(NGLLX,NGLLZ,nspec))
      allocate(rhof_k(npoin))
      allocate(sm_kl(NGLLX,NGLLZ,nspec))
      allocate(sm_k(npoin))
      allocate(eta_kl(NGLLX,NGLLZ,nspec))
      allocate(eta_k(npoin))
      allocate(mufr_kl(NGLLX,NGLLZ,nspec))
      allocate(mufr_k(npoin))
      allocate(B_kl(NGLLX,NGLLZ,nspec))
      allocate(B_k(npoin))
      allocate(C_kl(NGLLX,NGLLZ,nspec))
      allocate(C_k(npoin))
      allocate(M_kl(NGLLX,NGLLZ,nspec))
      allocate(M_k(npoin))
      allocate(rhob_kl(NGLLX,NGLLZ,nspec))
      allocate(rhofb_kl(NGLLX,NGLLZ,nspec))
      allocate(phi_kl(NGLLX,NGLLZ,nspec))
      allocate(Bb_kl(NGLLX,NGLLZ,nspec))
      allocate(Cb_kl(NGLLX,NGLLZ,nspec))
      allocate(Mb_kl(NGLLX,NGLLZ,nspec))
      allocate(mufrb_kl(NGLLX,NGLLZ,nspec))
      allocate(rhobb_kl(NGLLX,NGLLZ,nspec))
      allocate(rhofbb_kl(NGLLX,NGLLZ,nspec))
      allocate(phib_kl(NGLLX,NGLLZ,nspec))
      allocate(cpI_kl(NGLLX,NGLLZ,nspec))
      allocate(cpII_kl(NGLLX,NGLLZ,nspec))
      allocate(cs_kl(NGLLX,NGLLZ,nspec))
      allocate(ratio_kl(NGLLX,NGLLZ,nspec))
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
      allocate(rhot_kl(1,1,1))
      allocate(rhot_k(1))
      allocate(rhof_kl(1,1,1))
      allocate(rhof_k(1))
      allocate(sm_kl(1,1,1))
      allocate(sm_k(1))
      allocate(eta_kl(1,1,1))
      allocate(eta_k(1))
      allocate(mufr_kl(1,1,1))
      allocate(mufr_k(1))
      allocate(B_kl(1,1,1))
      allocate(B_k(1))
      allocate(C_kl(1,1,1))
      allocate(C_k(1))
      allocate(M_kl(1,1,1))
      allocate(M_k(1))
      allocate(rhob_kl(1,1,1))
      allocate(rhofb_kl(1,1,1))
      allocate(phi_kl(1,1,1))
      allocate(Bb_kl(1,1,1))
      allocate(Cb_kl(1,1,1))
      allocate(Mb_kl(1,1,1))
      allocate(mufrb_kl(1,1,1))
      allocate(rhobb_kl(1,1,1))
      allocate(rhofbb_kl(1,1,1))
      allocate(phib_kl(1,1,1))
      allocate(cpI_kl(1,1,1))
      allocate(cpII_kl(1,1,1))
      allocate(cs_kl(1,1,1))
      allocate(ratio_kl(1,1,1))
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
      allocate(icount(npoin))
    else
      allocate(icount(1))
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
    if(SIMULATION_TYPE == 2 .and. any_acoustic) then
      allocate(b_potential_acoustic(npoin))
      allocate(b_potential_dot_acoustic(npoin))
      allocate(b_potential_dot_dot_acoustic(npoin))
      allocate(b_displ_ac(2,npoin))
      allocate(b_accel_ac(2,npoin))
      allocate(accel_ac(2,npoin))
      allocate(rho_ac_kl(NGLLX,NGLLZ,nspec))
      allocate(rhol_ac_global(npoin))
      allocate(kappa_ac_kl(NGLLX,NGLLZ,nspec))
      allocate(kappal_ac_global(npoin))
      allocate(rhop_ac_kl(NGLLX,NGLLZ,nspec))
      allocate(alpha_ac_kl(NGLLX,NGLLZ,nspec))
      allocate(rhorho_ac_hessian_final2(NGLLX,NGLLZ,nspec))
      allocate(rhorho_ac_hessian_final1(NGLLX,NGLLZ,nspec))
    else
    ! allocate unused arrays with fictitious size
      allocate(b_potential_acoustic(1))
      allocate(b_potential_dot_acoustic(1))
      allocate(b_potential_dot_dot_acoustic(1))
      allocate(b_displ_ac(1,1))
      allocate(b_accel_ac(1,1))
      allocate(accel_ac(1,1))
      allocate(rho_ac_kl(1,1,1))
      allocate(rhol_ac_global(1))
      allocate(kappa_ac_kl(1,1,1))
      allocate(kappal_ac_global(1))
      allocate(rhop_ac_kl(1,1,1))
      allocate(alpha_ac_kl(1,1,1))
      allocate(rhorho_ac_hessian_final2(1,1,1))
      allocate(rhorho_ac_hessian_final1(1,1,1))
    endif

  endif ! ipass == 1

  !
  !---- build the global mass matrix 
  !
  call invert_mass_matrix_init(any_elastic,any_acoustic,any_poroelastic,npoin, &
                                rmass_inverse_elastic,&
                                rmass_inverse_acoustic, &
                                rmass_s_inverse_poroelastic, &
                                rmass_w_inverse_poroelastic, &
                                nspec,ibool,kmato,wxgll,wzgll,jacobian, &
                                elastic,poroelastic, &
                                assign_external_model,numat, &
                                density,poroelastcoef,porosity,tortuosity, &
                                vpext,rhoext)
  


#ifdef USE_MPI
  if ( nproc > 1 ) then

    ! preparing for MPI communications
    if(ipass == 1) allocate(mask_ispec_inner_outer(nspec))
    mask_ispec_inner_outer(:) = .false.

    call get_MPI(nspec,ibool,knods,ngnod,npoin,elastic,poroelastic, &
                    ninterface, max_interface_size, &
                    my_nelmnts_neighbours,my_interfaces,my_neighbours, &
                    ibool_interfaces_acoustic, ibool_interfaces_elastic, &
                    ibool_interfaces_poroelastic, &
                    nibool_interfaces_acoustic, nibool_interfaces_elastic, &
                    nibool_interfaces_poroelastic, &
                    inum_interfaces_acoustic, inum_interfaces_elastic, &
                    inum_interfaces_poroelastic, &
                    ninterface_acoustic, ninterface_elastic, ninterface_poroelastic, &
                    mask_ispec_inner_outer, &
                    myrank,ipass,coord)


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

    ! buffers for MPI communications
    max_ibool_interfaces_size_ac = maxval(nibool_interfaces_acoustic(:))
    max_ibool_interfaces_size_el = 3*maxval(nibool_interfaces_elastic(:))
    max_ibool_interfaces_size_po = NDIM*maxval(nibool_interfaces_poroelastic(:))
    if(ipass == 1) then
      allocate(tab_requests_send_recv_acoustic(ninterface_acoustic*2))
      allocate(buffer_send_faces_vector_ac(max_ibool_interfaces_size_ac,ninterface_acoustic))
      allocate(buffer_recv_faces_vector_ac(max_ibool_interfaces_size_ac,ninterface_acoustic))
      allocate(tab_requests_send_recv_elastic(ninterface_elastic*2))
      allocate(buffer_send_faces_vector_el(max_ibool_interfaces_size_el,ninterface_elastic))
      allocate(buffer_recv_faces_vector_el(max_ibool_interfaces_size_el,ninterface_elastic))
      allocate(tab_requests_send_recv_poro(ninterface_poroelastic*4))
      allocate(buffer_send_faces_vector_pos(max_ibool_interfaces_size_po,ninterface_poroelastic))
      allocate(buffer_recv_faces_vector_pos(max_ibool_interfaces_size_po,ninterface_poroelastic))
      allocate(buffer_send_faces_vector_pow(max_ibool_interfaces_size_po,ninterface_poroelastic))
      allocate(buffer_recv_faces_vector_pow(max_ibool_interfaces_size_po,ninterface_poroelastic))
    endif

! assembling the mass matrix
    call assemble_MPI_scalar(rmass_inverse_acoustic, &
                            rmass_inverse_elastic, &
                            rmass_s_inverse_poroelastic, &
                            rmass_w_inverse_poroelastic,npoin, &
                            ninterface, max_interface_size, max_ibool_interfaces_size_ac, &
                            max_ibool_interfaces_size_el, &
                            max_ibool_interfaces_size_po, &
                            ibool_interfaces_acoustic,ibool_interfaces_elastic, &
                            ibool_interfaces_poroelastic, &
                            nibool_interfaces_acoustic,nibool_interfaces_elastic, &
                            nibool_interfaces_poroelastic,my_neighbours)

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

    !  allocate(antecedent_list(nspec))

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

    !allocate(perm(nspec))

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

!============================================
!
!            end inner/outer passes
!
!============================================

!---
!---  end of section performed in two passes
!---

  call invert_mass_matrix(any_elastic,any_acoustic,any_poroelastic,npoin,rmass_inverse_elastic,&
              rmass_inverse_acoustic,rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic)

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
  call checkgrid(vpext,vsext,rhoext,density,poroelastcoef, &
                 porosity,tortuosity,permeability,ibool,kmato, &
                 coord,npoin,vpImin,vpImax,vpIImin,vpIImax, &
                 assign_external_model,nspec,UPPER_LIMIT_DISPLAY,numat,deltat, &
                 f0,initialfield,time_function_type, &
                 coorg,xinterp,zinterp,shape2D_display,knods,simulation_title, &
                 npgeo,pointsdisp,ngnod,any_elastic,any_poroelastic,all_anisotropic, &
                 myrank,nproc,NSOURCES,poroelastic, &
                 freq0,Q0,TURN_VISCATTENUATION_ON)

! convert receiver angle to radians
  anglerec = anglerec * pi / 180.d0

!
!---- for color images
!

  if(output_color_image) then
    ! prepares dimension of image
    call prepare_color_image_init(NX_IMAGE_color,NZ_IMAGE_color, &
                            xmin_color_image,xmax_color_image, &
                            zmin_color_image,zmax_color_image, &
                            coord,npoin,npgeo)

    ! allocate an array for image data
    allocate(image_color_data(NX_IMAGE_color,NZ_IMAGE_color))
    allocate(image_color_vp_display(NX_IMAGE_color,NZ_IMAGE_color))

    ! allocate an array for the grid point that corresponds to a given image data point
    allocate(iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color))
    allocate(copy_iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color))

    ! creates pixels indexing
    call prepare_color_image_pixels(myrank,NX_IMAGE_color,NZ_IMAGE_color, &
                            xmin_color_image,xmax_color_image, &
                            zmin_color_image,zmax_color_image, &
                            coord,npoin,coorg,npgeo,nspec,ngnod,knods,ibool, &
                            nb_pixel_loc,iglob_image_color)
  

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

    call MPI_GATHER( nb_pixel_loc, 1, MPI_INTEGER, nb_pixel_per_proc(1), &
                    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

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
#endif

    if (myrank == 0) write(IOUT,*) 'done locating all the pixels of color images'

  endif ! color_image

!
!---- initialize seismograms
!
  sisux = ZERO
  sisuz = ZERO

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
!----- Files where viscous damping are saved during forward wavefield calculation
!
  if(any_poroelastic .and. (SAVE_FORWARD .or. SIMULATION_TYPE .eq. 2)) then
    allocate(b_viscodampx(npoin))
    allocate(b_viscodampz(npoin))
    write(outputname,'(a,i6.6,a)') 'viscodampingx',myrank,'.bin'
    write(outputname2,'(a,i6.6,a)') 'viscodampingz',myrank,'.bin'
    if(SIMULATION_TYPE == 2) then
      reclen = CUSTOM_REAL * npoin
      open(unit=23,file='OUTPUT_FILES/'//outputname,status='old',&
            action='read',form='unformatted',access='direct',&
            recl=reclen)
      open(unit=24,file='OUTPUT_FILES/'//outputname2,status='old',&
            action='read',form='unformatted',access='direct',&
            recl=reclen)
    else
      reclen = CUSTOM_REAL * npoin
      open(unit=23,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',access='direct',&
            recl=reclen)
      open(unit=24,file='OUTPUT_FILES/'//outputname2,status='unknown',&
            form='unformatted',access='direct',&
            recl=reclen)
    endif
  else
    allocate(b_viscodampx(1))
    allocate(b_viscodampz(1))
  endif

!
!----- Files where absorbing signal are saved during forward wavefield calculation
!

  if( ((SAVE_FORWARD .and. SIMULATION_TYPE ==1) .or. SIMULATION_TYPE == 2) .and. anyabs ) then
    ! opens files for absorbing boundary data
    call prepare_absorb_files(myrank,any_elastic,any_poroelastic,any_acoustic, &
                      nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax,SIMULATION_TYPE)
  endif 

  if(anyabs .and. SIMULATION_TYPE == 2) then

    ! reads in absorbing bounday data
    if(any_elastic) then
      call prepare_absorb_elastic(NSTEP,p_sv, &
                      nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax, &
                      b_absorb_elastic_left,b_absorb_elastic_right, &
                      b_absorb_elastic_bottom,b_absorb_elastic_top)

    endif 
    if(any_poroelastic) then
      call prepare_absorb_poroelastic(NSTEP, &
                      nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax, &
                      b_absorb_poro_s_left,b_absorb_poro_w_left, &
                      b_absorb_poro_s_right,b_absorb_poro_w_right, &
                      b_absorb_poro_s_bottom,b_absorb_poro_w_bottom, &
                      b_absorb_poro_s_top,b_absorb_poro_w_top)

    endif 
    if(any_acoustic) then
      call prepare_absorb_acoustic(NSTEP, &
                      nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax, &
                      b_absorb_acoustic_left,b_absorb_acoustic_right, &
                      b_absorb_acoustic_bottom,b_absorb_acoustic_top)
    endif 

  endif ! if(anyabs .and. SIMULATION_TYPE == 2)



!
!----- Read last frame for backward wavefield calculation
!

  if(SIMULATION_TYPE == 2) then

    if(any_elastic) then
      write(outputname,'(a,i6.6,a)') 'snapshot_rho_kappa_mu_',myrank
      open(unit = 97, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'
      write(outputname,'(a,i6.6,a)') 'snapshot_rhop_alpha_beta_',myrank
      open(unit = 98, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'

      rho_kl(:,:,:) = ZERO
      mu_kl(:,:,:) = ZERO
      kappa_kl(:,:,:) = ZERO

      rhop_kl(:,:,:) = ZERO
      beta_kl(:,:,:) = ZERO
      alpha_kl(:,:,:) = ZERO
      rhorho_el_hessian_final2(:,:,:) = ZERO
      rhorho_el_hessian_temp2(:) = ZERO
      rhorho_el_hessian_final1(:,:,:) = ZERO
      rhorho_el_hessian_temp1(:) = ZERO
    endif

    if(any_poroelastic) then

      ! Primary kernels
      write(outputname,'(a,i6.6,a)') 'snapshot_mu_B_C_',myrank
      open(unit = 144, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'
      write(outputname,'(a,i6.6,a)') 'snapshot_M_rho_rhof_',myrank
      open(unit = 155, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'
      write(outputname,'(a,i6.6,a)') 'snapshot_m_eta_',myrank
      open(unit = 16, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'
      ! Wavespeed kernels
      write(outputname,'(a,i6.6,a)') 'snapshot_cpI_cpII_cs_',myrank
      open(unit = 20, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'
      write(outputname,'(a,i6.6,a)') 'snapshot_rhobb_rhofbb_ratio_',myrank
      open(unit = 21, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'
      write(outputname,'(a,i6.6,a)') 'snapshot_phib_eta_',myrank
      open(unit = 22, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'
      ! Density normalized kernels
      write(outputname,'(a,i6.6,a)') 'snapshot_mub_Bb_Cb_',myrank
      open(unit = 17, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'
      write(outputname,'(a,i6.6,a)') 'snapshot_Mb_rhob_rhofb_',myrank
      open(unit = 18, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'
      write(outputname,'(a,i6.6,a)') 'snapshot_mb_etab_',myrank
      open(unit = 19, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'

      rhot_kl(:,:,:) = ZERO
      rhof_kl(:,:,:) = ZERO
      eta_kl(:,:,:) = ZERO
      sm_kl(:,:,:) = ZERO
      mufr_kl(:,:,:) = ZERO
      B_kl(:,:,:) = ZERO
      C_kl(:,:,:) = ZERO
      M_kl(:,:,:) = ZERO

      rhob_kl(:,:,:) = ZERO
      rhofb_kl(:,:,:) = ZERO
      phi_kl(:,:,:) = ZERO
      mufrb_kl(:,:,:) = ZERO
      Bb_kl(:,:,:) = ZERO
      Cb_kl(:,:,:) = ZERO
      Mb_kl(:,:,:) = ZERO

      rhobb_kl(:,:,:) = ZERO
      rhofbb_kl(:,:,:) = ZERO
      phib_kl(:,:,:) = ZERO
      cs_kl(:,:,:) = ZERO
      cpI_kl(:,:,:) = ZERO
      cpII_kl(:,:,:) = ZERO
      ratio_kl(:,:,:) = ZERO
    endif

    if(any_acoustic) then
      write(outputname,'(a,i6.6,a)') 'snapshot_rho_kappa_',myrank
      open(unit = 95, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'
      write(outputname,'(a,i6.6,a)') 'snapshot_rhop_c_',myrank
      open(unit = 96, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing snapshot to disk'

      rho_ac_kl(:,:,:) = ZERO
      kappa_ac_kl(:,:,:) = ZERO

      rhop_ac_kl(:,:,:) = ZERO
      alpha_ac_kl(:,:,:) = ZERO
      rhorho_ac_hessian_final2(:,:,:) = ZERO
      rhorho_ac_hessian_final1(:,:,:) = ZERO
    endif

  endif ! if(SIMULATION_TYPE == 2)

!
!----  read initial fields from external file if needed
!

! if we are looking a plane wave beyond critical angle we use other method
  over_critical_angle = .false.

  if(initialfield) then
  
    ! Calculation of the initial field for a plane wave
    call prepare_initialfield(myrank,any_acoustic,any_poroelastic,over_critical_angle, &
                        NSOURCES,source_type,angleforce,x_source,z_source,f0, &
                        npoin,numat,poroelastcoef,density,coord, &
                        angleforce_refl,c_inc,c_refl,cploc,csloc,time_offset, &
                        A_plane, B_plane, C_plane, &
                        accel_elastic,veloc_elastic,displ_elastic)
    
    if( over_critical_angle ) then
    
      allocate(left_bound(nelemabs*NGLLX))
      allocate(right_bound(nelemabs*NGLLX))
      allocate(bot_bound(nelemabs*NGLLZ))

      call prepare_initialfield_paco(myrank,nelemabs,left_bound,right_bound,bot_bound, &
                                    numabs,codeabs,ibool,nspec, &
                                    source_type,NSOURCES,c_inc,c_refl, &
                                    count_bottom,count_left,count_right)

      allocate(v0x_left(count_left,NSTEP))
      allocate(v0z_left(count_left,NSTEP))
      allocate(t0x_left(count_left,NSTEP))
      allocate(t0z_left(count_left,NSTEP))

      allocate(v0x_right(count_right,NSTEP))
      allocate(v0z_right(count_right,NSTEP))
      allocate(t0x_right(count_right,NSTEP))
      allocate(t0z_right(count_right,NSTEP))

      allocate(v0x_bot(count_bottom,NSTEP))
      allocate(v0z_bot(count_bottom,NSTEP))
      allocate(t0x_bot(count_bottom,NSTEP))
      allocate(t0z_bot(count_bottom,NSTEP))

      allocate(displ_paco(NDIM,npoin))
      allocate(veloc_paco(NDIM,npoin))
      allocate(accel_paco(NDIM,npoin))

      ! call Paco's routine to compute in frequency and convert to time by Fourier transform
      call paco_beyond_critical(coord,npoin,deltat,NSTEP,angleforce(1),&
              f0(1),cploc,csloc,TURN_ATTENUATION_ON,Qp_attenuation,source_type(1),v0x_left,v0z_left,&
              v0x_right,v0z_right,v0x_bot,v0z_bot,t0x_left,t0z_left,t0x_right,t0z_right,&
              t0x_bot,t0z_bot,left_bound(1:count_left),right_bound(1:count_right),bot_bound(1:count_bottom)&
              ,count_left,count_right,count_bottom,displ_paco,veloc_paco,accel_paco)

      displ_elastic(1,:) = displ_paco(1,:)
      displ_elastic(3,:) = displ_paco(2,:)
      veloc_elastic(1,:) = veloc_paco(1,:)
      veloc_elastic(3,:) = veloc_paco(2,:)
      accel_elastic(1,:) = accel_paco(1,:)
      accel_elastic(3,:) = accel_paco(2,:)

      deallocate(left_bound)
      deallocate(right_bound)
      deallocate(bot_bound)

      deallocate(displ_paco)
      deallocate(veloc_paco)
      deallocate(accel_paco)

      if (myrank == 0) then
        write(IOUT,*)  '***********'
        write(IOUT,*)  'done calculating the initial wave field'
        write(IOUT,*)  '***********'
      endif

    endif ! beyond critical angle

    write(IOUT,*) 'Max norm of initial elastic displacement = ', &
      maxval(sqrt(displ_elastic(1,:)**2 + displ_elastic(3,:)**2))

  endif ! initialfield

  deltatsquare = deltat * deltat
  deltatcube = deltatsquare * deltat
  deltatfourth = deltatsquare * deltatsquare

  twelvedeltat = 12.d0 * deltat
  fourdeltatsquare = 4.d0 * deltatsquare

! compute the source time function and store it in a text file
  if(.not. initialfield) then

    allocate(source_time_function(NSOURCES,NSTEP))
    source_time_function(:,:) = 0.0_CUSTOM_REAL

    ! computes source time function array
    call prepare_source_time_function(myrank,NSTEP,NSOURCES,source_time_function, &
                          time_function_type,f0,tshift_src,factor,aval, &
                          t0,nb_proc_source,deltat)
  else
    ! uses an initialfield
    ! dummy allocation
    allocate(source_time_function(1,1))
  endif

! determine if coupled fluid-solid simulation
  coupled_acoustic_elastic = any_acoustic .and. any_elastic
  coupled_acoustic_poro = any_acoustic .and. any_poroelastic

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
  if(coupled_acoustic_poro) then
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
  if(coupled_acoustic_poro .and. anyabs) then

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
  coupled_elastic_poro = any_elastic .and. any_poroelastic

! solid/porous edge detection
! the two elements forming an edge are already known (computed in meshfem2D),
! the common nodes forming the edge are computed here
  if(coupled_elastic_poro) then

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
  if(coupled_elastic_poro .and. anyabs) then

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
  if(OUTPUT_ENERGY) stop 'energy calculation currently serial only, should add an MPI_REDUCE in parallel'
#endif
! open the file in which we will store the energy curve
  if(OUTPUT_ENERGY) open(unit=IOUT_ENERGY,file='energy.gnu',status='unknown')

!
!----          s t a r t   t i m e   i t e r a t i o n s
!
  if (myrank == 0) write(IOUT,400)

  ! count elapsed wall-clock time
  call date_and_time(datein,timein,zone,time_values)
  ! time_values(1): year
  ! time_values(2): month of the year
  ! time_values(3): day of the month
  ! time_values(5): hour of the day
  ! time_values(6): minutes of the hour
  ! time_values(7): seconds of the minute
  ! time_values(8): milliseconds of the second
  ! this fails if we cross the end of the month
  time_start = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
               60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0
  month_start = time_values(2)
  year_start = time_values(1)

  ! prepares image background
  if(output_color_image) then
    call prepare_color_image_vp(npoin,image_color_vp_display,iglob_image_color, &
                            NX_IMAGE_color,NZ_IMAGE_color,nb_pixel_loc, &
                            num_pixel_loc,nspec,poroelastic,ibool,kmato, &
                            numat,density,poroelastcoef,porosity,tortuosity, &
                            nproc,myrank,assign_external_model,vpext)

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
      displ_elastic = displ_elastic &
                    + deltat*veloc_elastic &
                    + deltatsquareover2*accel_elastic
      veloc_elastic = veloc_elastic + deltatover2*accel_elastic
      accel_elastic = ZERO

      if(SIMULATION_TYPE == 2) then ! Adjoint calculation
        b_displ_elastic = b_displ_elastic &
                        + b_deltat*b_veloc_elastic &
                        + b_deltatsquareover2*b_accel_elastic
        b_veloc_elastic = b_veloc_elastic + b_deltatover2*b_accel_elastic
        b_accel_elastic = ZERO
      endif
    endif

    if(any_poroelastic) then
      !for the solid
      displs_poroelastic = displs_poroelastic &
                         + deltat*velocs_poroelastic &
                         + deltatsquareover2*accels_poroelastic
      velocs_poroelastic = velocs_poroelastic + deltatover2*accels_poroelastic
      accels_poroelastic = ZERO
      !for the fluid
      displw_poroelastic = displw_poroelastic &
                         + deltat*velocw_poroelastic &
                         + deltatsquareover2*accelw_poroelastic
      velocw_poroelastic = velocw_poroelastic + deltatover2*accelw_poroelastic
      accelw_poroelastic = ZERO

      if(SIMULATION_TYPE == 2) then ! Adjoint calculation
        !for the solid
        b_displs_poroelastic = b_displs_poroelastic &
                             + b_deltat*b_velocs_poroelastic &
                             + b_deltatsquareover2*b_accels_poroelastic
        b_velocs_poroelastic = b_velocs_poroelastic + b_deltatover2*b_accels_poroelastic
        b_accels_poroelastic = ZERO
        !for the fluid
        b_displw_poroelastic = b_displw_poroelastic &
                             + b_deltat*b_velocw_poroelastic &
                             + b_deltatsquareover2*b_accelw_poroelastic
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
          stop 'Permeability matrix is not invertible'
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
            rx_viscous(i,j,ispec) = alphaval * rx_viscous(i,j,ispec) &
                                  + betaval * Sn + gammaval * Snp1

            ! evolution rz_viscous
            Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscoz(i,j,ispec)
            Snp1 = - (1.d0 - theta_e/theta_s)/theta_s*viscoz_loc(i,j)
            rz_viscous(i,j,ispec) = alphaval * rz_viscous(i,j,ispec) &
                                  + betaval * Sn + gammaval * Snp1


          enddo
        enddo

        ! save visco for Runge-Kutta scheme
        viscox(:,:,ispec) = viscox_loc(:,:)
        viscoz(:,:,ispec) = viscoz_loc(:,:)

      enddo   ! end of spectral element loop
    endif ! end of viscous attenuation for porous media

!-----------------------------------------
    if(any_acoustic) then

      ! Newmark time scheme
      potential_acoustic = potential_acoustic &
                          + deltat*potential_dot_acoustic &
                          + deltatsquareover2*potential_dot_dot_acoustic
      potential_dot_acoustic = potential_dot_acoustic &
                              + deltatover2*potential_dot_dot_acoustic
      potential_dot_dot_acoustic = ZERO

      if(SIMULATION_TYPE == 2) then ! Adjoint calculation
        b_potential_acoustic = b_potential_acoustic &
                            + b_deltat*b_potential_dot_acoustic &
                            + b_deltatsquareover2*b_potential_dot_dot_acoustic
        b_potential_dot_acoustic = b_potential_dot_acoustic &
                                  + b_deltatover2*b_potential_dot_dot_acoustic
        b_potential_dot_dot_acoustic = ZERO
      endif

      ! free surface for an acoustic medium
      if ( nelem_acoustic_surface > 0 ) then
        call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                          potential_acoustic,acoustic_surface, &
                                          ibool,nelem_acoustic_surface,npoin,nspec)

        if(SIMULATION_TYPE == 2) then ! Adjoint calculation
          call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                            b_potential_acoustic,acoustic_surface, &
                                            ibool,nelem_acoustic_surface,npoin,nspec)
        endif
      endif

! *********************************************************
! ************* compute forces for the acoustic elements
! *********************************************************

!      call compute_forces_acoustic(npoin,nspec,nelemabs,numat,it,NSTEP, &
!               anyabs,assign_external_model,ibool,kmato,numabs, &
!               elastic,poroelastic,codeabs,potential_dot_dot_acoustic,potential_dot_acoustic, &
!               potential_acoustic,b_potential_dot_dot_acoustic,b_potential_acoustic, &
!               density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
!               vpext,rhoext,hprime_xx,hprimewgll_xx, &
!               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
!               ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
!               jbegin_left,jend_left,jbegin_right,jend_right, &
!               SIMULATION_TYPE,SAVE_FORWARD,b_absorb_acoustic_left,&
!               b_absorb_acoustic_right,b_absorb_acoustic_bottom,&
!               b_absorb_acoustic_top,nspec_xmin,nspec_xmax,&
!               nspec_zmin,nspec_zmax,ib_left,ib_right,ib_bottom,ib_top)


      call compute_forces_acoustic_2(npoin,nspec,nelemabs,numat,it,NSTEP, &
               anyabs,assign_external_model,ibool,kmato,numabs, &
               elastic,poroelastic,codeabs,potential_dot_dot_acoustic,potential_dot_acoustic, &
               potential_acoustic, &
               density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,rhoext,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
               jbegin_left,jend_left,jbegin_right,jend_right, &
               SIMULATION_TYPE,SAVE_FORWARD,nspec_xmin,nspec_xmax,&
               nspec_zmin,nspec_zmax,ib_left,ib_right,ib_bottom,ib_top, &
               b_absorb_acoustic_left,b_absorb_acoustic_right, &
               b_absorb_acoustic_bottom,b_absorb_acoustic_top)
      if( SIMULATION_TYPE == 2 ) then
        call compute_forces_acoustic_2(npoin,nspec,nelemabs,numat,it,NSTEP, &
               anyabs,assign_external_model,ibool,kmato,numabs, &
               elastic,poroelastic,codeabs,b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
               b_potential_acoustic, &
               density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,rhoext,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
               jbegin_left,jend_left,jbegin_right,jend_right, &
               SIMULATION_TYPE,SAVE_FORWARD,nspec_xmin,nspec_xmax,&
               nspec_zmin,nspec_zmax,ib_left,ib_right,ib_bottom,ib_top, &
               b_absorb_acoustic_left,b_absorb_acoustic_right, &
               b_absorb_acoustic_bottom,b_absorb_acoustic_top)          
      endif


      ! stores absorbing boundary contributions into files
      if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
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
      endif ! if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1)

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
          displ_z = displ_elastic(3,iglob)

          if(SIMULATION_TYPE == 2) then
            b_displ_x = b_displ_elastic(1,iglob)
            b_displ_z = b_displ_elastic(3,iglob)
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

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

          if(SIMULATION_TYPE == 2) then
          b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) +&
                      weight*(b_displ_x*nx + b_displ_z*nz)
          endif !if(SIMULATION_TYPE == 2) then

        enddo

      enddo

    endif

! *********************************************************
! ************* add coupling with the poroelastic side
! *********************************************************

    if(coupled_acoustic_poro) then

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

          if(SIMULATION_TYPE == 2) then
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

          if(SIMULATION_TYPE == 2) then
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) &
                   + weight*((b_displ_x + b_displw_x)*nx + (b_displ_z + b_displw_z)*nz)
          endif 

        enddo

      enddo

    endif


! ************************************************************************************
! ************************************ add force source
! ************************************************************************************

    if(any_acoustic) then

! --- add the source
      if(.not. initialfield) then

        do i_source=1,NSOURCES
          ! if this processor carries the source and the source element is acoustic
          if (is_proc_source(i_source) == 1 .and. &
            .not. elastic(ispec_selected_source(i_source)) .and. &
            .not. poroelastic(ispec_selected_source(i_source))) then
            
! collocated force
! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
! the sign is negative because pressure p = - Chi_dot_dot therefore we need
! to add minus the source to Chi_dot_dot to get plus the source in pressure
            if(source_type(i_source) == 1) then

              if(SIMULATION_TYPE == 1) then  
                ! forward wavefield
                do j = 1,NGLLZ
                  do i = 1,NGLLX
                    iglob = ibool(i,j,ispec_selected_source(i_source))
                    hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                    potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                            - source_time_function(i_source,it)*hlagrange
                  enddo
                enddo
              else                   
                ! backward wavefield
                do j = 1,NGLLZ
                  do i = 1,NGLLX
                    iglob = ibool(i,j,ispec_selected_source(i_source))
                    hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                    b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) &
                                          - source_time_function(i_source,NSTEP-it+1)*hlagrange
                  enddo
                enddo
              endif

            ! moment tensor
            else if(source_type(i_source) == 2) then
              call exit_MPI('cannot have moment tensor source in acoustic element')

            endif
          endif ! if this processor carries the source and the source element is acoustic
        enddo ! do i_source=1,NSOURCES

        if(SIMULATION_TYPE == 2) then   ! adjoint wavefield
          irec_local = 0
          do irec = 1,nrec
            !   add the source (only if this proc carries the source)
            if (myrank == which_proc_receiver(irec)) then

              irec_local = irec_local + 1
              if (.not. elastic(ispec_selected_rec(irec)) .and. &
                 .not. poroelastic(ispec_selected_rec(irec))) then
                ! add source array
                do j=1,NGLLZ
                  do i=1,NGLLX
                    iglob = ibool(i,j,ispec_selected_rec(irec))
                    potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                  - adj_sourcearrays(irec_local,NSTEP-it+1,1,i,j)
                  enddo
                enddo
              endif ! if element acoustic

            endif ! if this processor carries the adjoint source
          enddo ! irec = 1,nrec
        endif ! SIMULATION_TYPE == 2 adjoint wavefield

      endif ! if not using an initial field
      
    endif !if(any_acoustic)


! assembling potential_dot_dot for acoustic elements
#ifdef USE_MPI
  if ( nproc > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_ac(potential_dot_dot_acoustic,npoin, &
                    ninterface, ninterface_acoustic,inum_interfaces_acoustic, &
                    max_interface_size, max_ibool_interfaces_size_ac,&
                    ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
                    tab_requests_send_recv_acoustic,buffer_send_faces_vector_ac, &
                    buffer_recv_faces_vector_ac, my_neighbours)
           
    if ( SIMULATION_TYPE == 2) then
      call assemble_MPI_vector_ac(b_potential_dot_dot_acoustic,npoin, &
                     ninterface, ninterface_acoustic,inum_interfaces_acoustic, &
                     max_interface_size, max_ibool_interfaces_size_ac,&
                     ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
                     tab_requests_send_recv_acoustic,buffer_send_faces_vector_ac, &
                     buffer_recv_faces_vector_ac, my_neighbours)
          
    endif
           
  endif

!  if ( nproc > 1 .and. any_acoustic .and. ninterface_acoustic > 0 .and. SIMULATION_TYPE == 2) then
!    call assemble_MPI_vector_ac(b_potential_dot_dot_acoustic,npoin, &
!           ninterface, ninterface_acoustic,inum_interfaces_acoustic, &
!           max_interface_size, max_ibool_interfaces_size_ac,&
!           ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
!           tab_requests_send_recv_acoustic,buffer_send_faces_vector_ac, &
!           buffer_recv_faces_vector_ac, my_neighbours)
!  endif
#endif

! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

  if(any_acoustic) then

    potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic
    potential_dot_acoustic = potential_dot_acoustic + deltatover2*potential_dot_dot_acoustic

    if(SIMULATION_TYPE ==2)then
    b_potential_dot_dot_acoustic = b_potential_dot_dot_acoustic * rmass_inverse_acoustic
    b_potential_dot_acoustic = b_potential_dot_acoustic + b_deltatover2*b_potential_dot_dot_acoustic
    endif


! free surface for an acoustic medium
    if ( nelem_acoustic_surface > 0 ) then
      call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                        potential_acoustic,acoustic_surface, &
                                        ibool,nelem_acoustic_surface,npoin,nspec)

      if(SIMULATION_TYPE == 2) then
        call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                          b_potential_acoustic,acoustic_surface, &
                                          ibool,nelem_acoustic_surface,npoin,nspec)
      endif

    endif

  endif !if(any_acoustic)


! *********************************************************
! ************* main solver for the elastic elements
! *********************************************************

 if(any_elastic) then
    call compute_forces_viscoelastic(p_sv,npoin,nspec,myrank,nelemabs,numat, &
               ispec_selected_source,ispec_selected_rec,is_proc_source,which_proc_receiver, &
               source_type,it,NSTEP,anyabs,assign_external_model, &
               initialfield,TURN_ATTENUATION_ON,angleforce,deltatcube, &
               deltatfourth,twelvedeltat,fourdeltatsquare,ibool,kmato,numabs,elastic,codeabs, &
               accel_elastic,veloc_elastic,displ_elastic,b_accel_elastic,b_displ_elastic, &
               density,poroelastcoef,xix,xiz,gammax,gammaz, &
               jacobian,vpext,vsext,rhoext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,anisotropic,anisotropy, &
               source_time_function,sourcearray,adj_sourcearrays, &
               e1,e11,e13,dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n, &
               dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll,inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,Mu_nu1,Mu_nu2,N_SLS, &
               deltat,coord,add_Bielak_conditions, x0_source, z0_source, &
               A_plane, B_plane, C_plane, angleforce_refl, c_inc, c_refl, time_offset, f0(1),&
               v0x_left(1,it),v0z_left(1,it),v0x_right(1,it),v0z_right(1,it),v0x_bot(1,it),v0z_bot(1,it), &
               t0x_left(1,it),t0z_left(1,it),t0x_right(1,it),t0z_right(1,it),t0x_bot(1,it),t0z_bot(1,it), &
               count_left,count_right,count_bottom,over_critical_angle,NSOURCES,nrec,SIMULATION_TYPE,SAVE_FORWARD, &
               b_absorb_elastic_left,b_absorb_elastic_right,b_absorb_elastic_bottom,b_absorb_elastic_top, &
               nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax,ib_left,ib_right,ib_bottom,ib_top,mu_k,kappa_k)

    if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
!--- left absorbing boundary
      if(nspec_xmin >0) then
      do ispec = 1,nspec_xmin

      if(p_sv)then!P-SV waves
        do i=1,NGLLZ
          write(35) b_absorb_elastic_left(1,i,ispec,it)
        enddo
        do i=1,NGLLZ
          write(35) b_absorb_elastic_left(3,i,ispec,it)
        enddo
      else!SH (membrane) waves
        do i=1,NGLLZ
          write(35) b_absorb_elastic_left(2,i,ispec,it)
        enddo
      endif

      enddo
      endif

!--- right absorbing boundary
      if(nspec_xmax >0) then
      do ispec = 1,nspec_xmax


      if(p_sv)then!P-SV waves
         do i=1,NGLLZ
     write(36) b_absorb_elastic_right(1,i,ispec,it)
         enddo
         do i=1,NGLLZ
     write(36) b_absorb_elastic_right(3,i,ispec,it)
         enddo
      else!SH (membrane) waves
         do i=1,NGLLZ
     write(36) b_absorb_elastic_right(2,i,ispec,it)
         enddo
      endif

      enddo
      endif

!--- bottom absorbing boundary
      if(nspec_zmin >0) then
      do ispec = 1,nspec_zmin

      if(p_sv)then!P-SV waves
         do i=1,NGLLX
     write(37) b_absorb_elastic_bottom(1,i,ispec,it)
         enddo
         do i=1,NGLLX
     write(37) b_absorb_elastic_bottom(3,i,ispec,it)
         enddo
      else!SH (membrane) waves
         do i=1,NGLLX
     write(37) b_absorb_elastic_bottom(2,i,ispec,it)
         enddo
      endif

      enddo
      endif

!--- top absorbing boundary
      if(nspec_zmax >0) then
      do ispec = 1,nspec_zmax

      if(p_sv)then!P-SV waves
         do i=1,NGLLX
     write(38) b_absorb_elastic_top(1,i,ispec,it)
         enddo
         do i=1,NGLLX
     write(38) b_absorb_elastic_top(3,i,ispec,it)
         enddo
      else!SH (membrane) waves
         do i=1,NGLLX
     write(38) b_absorb_elastic_top(2,i,ispec,it)
         enddo
      endif

      enddo
      endif

    endif ! if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1)

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
          if(SIMULATION_TYPE == 2) then
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
          accel_elastic(3,iglob) = accel_elastic(3,iglob) + weight*nz*pressure

          if(SIMULATION_TYPE == 2) then
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) + weight*nx*b_pressure
          b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) + weight*nz*b_pressure
          endif !if(SIMULATION_TYPE == 2) then

        enddo

      enddo

    endif

! ****************************************************************************
! ************* add coupling with the poroelastic side
! ****************************************************************************
    if(coupled_elastic_poro) then

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

          if(SIMULATION_TYPE == 2) then
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
            if(SIMULATION_TYPE == 2) then
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

          if(SIMULATION_TYPE == 2) then
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

    if(SIMULATION_TYPE == 2) then
    b_sigma_xx = lambdalplus2mul_G*b_dux_dxl + lambdal_G*b_duz_dzl + C_biot*(b_dwx_dxl + b_dwz_dzl)
    b_sigma_xz = mul_G*(b_duz_dxl + b_dux_dzl)
    b_sigma_zz = lambdalplus2mul_G*b_duz_dzl + lambdal_G*b_dux_dxl + C_biot*(b_dwx_dxl + b_dwz_dzl)
    endif
! get point values for the elastic domain, which matches our side in the inverse direction
          ii2 = ivalue(ipoin1D,iedge_elastic)
          jj2 = jvalue(ipoin1D,iedge_elastic)
          iglob = ibool(ii2,jj2,ispec_elastic)

! get elastic properties
    lambdal_relaxed = poroelastcoef(1,1,kmato(ispec_elastic))
    mul_relaxed = poroelastcoef(2,1,kmato(ispec_elastic))
    lambdalplus2mul_relaxed = poroelastcoef(3,1,kmato(ispec_elastic))

! derivative along x and along z for u_s and w
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          if(SIMULATION_TYPE == 2) then
          b_dux_dxi = ZERO
          b_duz_dxi = ZERO

          b_dux_dgamma = ZERO
          b_duz_dgamma = ZERO
          endif

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
            duz_dxi = duz_dxi + displ_elastic(3,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)
            duz_dgamma = duz_dgamma + displ_elastic(3,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)

            if(SIMULATION_TYPE == 2) then
            b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
            b_duz_dxi = b_duz_dxi + b_displ_elastic(3,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
            b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)
            b_duz_dgamma = b_duz_dgamma + b_displ_elastic(3,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)
            endif
          enddo

          xixl = xix(ii2,jj2,ispec_elastic)
          xizl = xiz(ii2,jj2,ispec_elastic)
          gammaxl = gammax(ii2,jj2,ispec_elastic)
          gammazl = gammaz(ii2,jj2,ispec_elastic)

! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          if(SIMULATION_TYPE == 2) then
          b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
          b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

          b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
          b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl
          endif
! compute stress tensor
! full anisotropy
  if(kmato(ispec_elastic) == 2) then
! implement anisotropy in 2D
      if(assign_external_model) then
         c11 = c11ext(ii2,jj2,ispec_elastic)
         c13 = c13ext(ii2,jj2,ispec_elastic)
         c15 = c15ext(ii2,jj2,ispec_elastic)
         c33 = c33ext(ii2,jj2,ispec_elastic)
         c35 = c35ext(ii2,jj2,ispec_elastic)
         c55 = c55ext(ii2,jj2,ispec_elastic)
      else
         c11 = anisotropy(1,kmato(ispec_elastic))
         c13 = anisotropy(2,kmato(ispec_elastic))
         c15 = anisotropy(3,kmato(ispec_elastic))
         c33 = anisotropy(4,kmato(ispec_elastic))
         c35 = anisotropy(5,kmato(ispec_elastic))
         c55 = anisotropy(6,kmato(ispec_elastic))
      end if

     sigma_xx = sigma_xx + c11*dux_dxl + c15*(duz_dxl + dux_dzl) + c13*duz_dzl
     sigma_zz = sigma_zz + c13*dux_dxl + c35*(duz_dxl + dux_dzl) + c33*duz_dzl
     sigma_xz = sigma_xz + c15*dux_dxl + c55*(duz_dxl + dux_dzl) + c35*duz_dzl
  else
! no attenuation
    sigma_xx = sigma_xx + lambdalplus2mul_relaxed*dux_dxl + lambdal_relaxed*duz_dzl
    sigma_xz = sigma_xz + mul_relaxed*(duz_dxl + dux_dzl)
    sigma_zz = sigma_zz + lambdalplus2mul_relaxed*duz_dzl + lambdal_relaxed*dux_dxl
  endif

    if(SIMULATION_TYPE == 2) then
    b_sigma_xx = b_sigma_xx + lambdalplus2mul_relaxed*b_dux_dxl + lambdal_relaxed*b_duz_dzl
    b_sigma_xz = b_sigma_xz + mul_relaxed*(b_duz_dxl + b_dux_dzl)
    b_sigma_zz = b_sigma_zz + lambdalplus2mul_relaxed*b_duz_dzl + lambdal_relaxed*b_dux_dxl
    endif

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
                (sigma_xx*nx + sigma_xz*nz)/2.d0

          accel_elastic(3,iglob) = accel_elastic(3,iglob) - weight* &
                (sigma_xz*nx + sigma_zz*nz)/2.d0

          if(SIMULATION_TYPE == 2) then
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - weight* &
                (b_sigma_xx*nx + b_sigma_xz*nz)/2.d0

          b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) - weight* &
                (b_sigma_xz*nx + b_sigma_zz*nz)/2.d0
          endif !if(SIMULATION_TYPE == 2) then

        enddo

      enddo

    endif


! ************************************************************************************
! ************************************ add force source
! ************************************************************************************

  if(any_elastic) then

! --- add the source if it is a collocated force
    if(.not. initialfield) then

    do i_source=1,NSOURCES
! if this processor carries the source and the source element is elastic
      if (is_proc_source(i_source) == 1 .and. elastic(ispec_selected_source(i_source))) then

! collocated force
        if(source_type(i_source) == 1) then
       if(SIMULATION_TYPE == 1) then  ! forward wavefield

          if(p_sv) then ! P-SV calculation
          do j = 1,NGLLZ
           do i = 1,NGLLX
             iglob = ibool(i,j,ispec_selected_source(i_source))
             hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
          accel_elastic(1,iglob) = accel_elastic(1,iglob) &
              - sin(angleforce(i_source))*source_time_function(i_source,it)*hlagrange
          accel_elastic(3,iglob) = accel_elastic(3,iglob) &
              + cos(angleforce(i_source))*source_time_function(i_source,it)*hlagrange
           enddo
          enddo
          else    ! SH (membrane) calculation
          do j = 1,NGLLZ
           do i = 1,NGLLX
             iglob = ibool(i,j,ispec_selected_source(i_source))
             hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
          accel_elastic(2,iglob) = accel_elastic(2,iglob) &
                            + source_time_function(i_source,it)*hlagrange
           enddo
          enddo
          endif

       else                   ! backward wavefield

          if(p_sv) then ! P-SV calculation
          do j = 1,NGLLZ
           do i = 1,NGLLX
             iglob = ibool(i,j,ispec_selected_source(i_source))
             hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
      b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) &
            - sin(angleforce(i_source))*source_time_function(i_source,NSTEP-it+1) &
            *hlagrange
      b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) &
            + cos(angleforce(i_source))*source_time_function(i_source,NSTEP-it+1) &
            *hlagrange
           enddo
          enddo
          else    ! SH (membrane) calculation
          do j = 1,NGLLZ
           do i = 1,NGLLX
             iglob = ibool(i,j,ispec_selected_source(i_source))
             hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
      b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) &
                            + source_time_function(i_source,NSTEP-it+1)*hlagrange
           enddo
          enddo

          endif

       endif  !endif SIMULATION_TYPE == 1
        endif

      endif ! if this processor carries the source and the source element is elastic
    enddo ! do i_source=1,NSOURCES

    endif ! if not using an initial field
  endif !if(any_elastic)

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

  if (nproc > 1 .and. any_elastic .and. ninterface_elastic > 0 .and. SIMULATION_TYPE == 2) then
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
    accel_elastic(1,:) = accel_elastic(1,:) * rmass_inverse_elastic
    accel_elastic(2,:) = accel_elastic(2,:) * rmass_inverse_elastic
    accel_elastic(3,:) = accel_elastic(3,:) * rmass_inverse_elastic

    veloc_elastic = veloc_elastic + deltatover2*accel_elastic

   if(SIMULATION_TYPE == 2) then
    b_accel_elastic(1,:) = b_accel_elastic(1,:) * rmass_inverse_elastic(:)
    b_accel_elastic(2,:) = b_accel_elastic(2,:) * rmass_inverse_elastic(:)
    b_accel_elastic(3,:) = b_accel_elastic(3,:) * rmass_inverse_elastic(:)

    b_veloc_elastic = b_veloc_elastic + b_deltatover2*b_accel_elastic
   endif

  endif !if(any_elastic)


! ******************************************************************************************************************
! ************* main solver for the poroelastic elements: first the solid (u_s) than the fluid (w)
! ******************************************************************************************************************

  if(any_poroelastic) then

    if(SIMULATION_TYPE == 2) then
! if inviscid fluid, comment the reading and uncomment the zeroing
!     read(23,rec=NSTEP-it+1) b_viscodampx
!     read(24,rec=NSTEP-it+1) b_viscodampz
     b_viscodampx(:) = ZERO
     b_viscodampz(:) = ZERO
    endif

    call compute_forces_poro_solid(npoin,nspec,myrank,nelemabs,numat, &
               ispec_selected_source,ispec_selected_rec,is_proc_source,which_proc_receiver,&
               source_type,it,NSTEP,anyabs, &
               initialfield,TURN_ATTENUATION_ON,TURN_VISCATTENUATION_ON,deltatcube, &
               deltatfourth,twelvedeltat,fourdeltatsquare,ibool,kmato,numabs,poroelastic,codeabs, &
               accels_poroelastic,velocs_poroelastic,velocw_poroelastic,displs_poroelastic,displw_poroelastic,&
               b_accels_poroelastic,b_displs_poroelastic,b_displw_poroelastic,&
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
               mufr_k,B_k,NSOURCES,nrec,SIMULATION_TYPE,SAVE_FORWARD,&
               b_absorb_poro_s_left,b_absorb_poro_s_right,b_absorb_poro_s_bottom,b_absorb_poro_s_top,&
               nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax,ib_left,ib_right,ib_bottom,ib_top,f0(1),freq0,Q0)



    call compute_forces_poro_fluid(npoin,nspec,myrank,nelemabs,numat, &
               ispec_selected_source,ispec_selected_rec,is_proc_source,which_proc_receiver,&
               source_type,it,NSTEP,anyabs, &
               initialfield,TURN_ATTENUATION_ON,TURN_VISCATTENUATION_ON,deltatcube, &
               deltatfourth,twelvedeltat,fourdeltatsquare,ibool,kmato,numabs,poroelastic,codeabs, &
               accelw_poroelastic,velocw_poroelastic,displw_poroelastic,velocs_poroelastic,displs_poroelastic,&
               b_accelw_poroelastic,b_displw_poroelastic,b_displs_poroelastic,&
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
               C_k,M_k,NSOURCES,nrec,SIMULATION_TYPE,SAVE_FORWARD,&
               b_absorb_poro_w_left,b_absorb_poro_w_right,b_absorb_poro_w_bottom,b_absorb_poro_w_top,&
               nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax,ib_left,ib_right,ib_bottom,ib_top,f0(1),freq0,Q0)


    if(SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
! if inviscid fluid, comment
!     write(23,rec=it) b_viscodampx
!     write(24,rec=it) b_viscodampz
    endif

    if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1) then

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

    endif ! if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1)

  endif !if(any_poroelastic) then

! *********************************************************
! ************* add coupling with the acoustic side
! *********************************************************

    if(coupled_acoustic_poro) then

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
          if(SIMULATION_TYPE == 2) then
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

          if(SIMULATION_TYPE == 2) then
! contribution to the solid phase
          b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) + weight*nx*b_pressure*(1._CUSTOM_REAL-phil/tortl)
          b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + weight*nz*b_pressure*(1._CUSTOM_REAL-phil/tortl)

! contribution to the fluid phase
          b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) + weight*nx*b_pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
          b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) + weight*nz*b_pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
          endif !if(SIMULATION_TYPE == 2) then

        enddo ! do ipoin1D = 1,NGLLX

      enddo ! do inum = 1,num_fluid_poro_edges

    endif ! if(coupled_acoustic_poro)

! ****************************************************************************
! ************* add coupling with the elastic side
! ****************************************************************************

    if(coupled_elastic_poro) then

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

! get elastic properties
    lambdal_relaxed = poroelastcoef(1,1,kmato(ispec_elastic))
    mul_relaxed = poroelastcoef(2,1,kmato(ispec_elastic))
    lambdalplus2mul_relaxed = poroelastcoef(3,1,kmato(ispec_elastic))

! derivative along x and along z for u_s and w
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          if(SIMULATION_TYPE == 2) then
          b_dux_dxi = ZERO
          b_duz_dxi = ZERO

          b_dux_dgamma = ZERO
          b_duz_dgamma = ZERO
          endif

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec_elastic))*hprime_zz(j,k)

            if(SIMULATION_TYPE == 2) then
            b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
            b_duz_dxi = b_duz_dxi + b_displ_elastic(3,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
            b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
            b_duz_dgamma = b_duz_dgamma + b_displ_elastic(3,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
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

          if(SIMULATION_TYPE == 2) then
          b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
          b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

          b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
          b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl
          endif
! compute stress tensor
! full anisotropy
  if(kmato(ispec_elastic) == 2) then
! implement anisotropy in 2D
      if(assign_external_model) then
         c11 = c11ext(i,j,ispec_elastic)
         c13 = c13ext(i,j,ispec_elastic)
         c15 = c15ext(i,j,ispec_elastic)
         c33 = c33ext(i,j,ispec_elastic)
         c35 = c35ext(i,j,ispec_elastic)
         c55 = c55ext(i,j,ispec_elastic)
      else
         c11 = anisotropy(1,kmato(ispec_elastic))
         c13 = anisotropy(2,kmato(ispec_elastic))
         c15 = anisotropy(3,kmato(ispec_elastic))
         c33 = anisotropy(4,kmato(ispec_elastic))
         c35 = anisotropy(5,kmato(ispec_elastic))
         c55 = anisotropy(6,kmato(ispec_elastic))
      end if
     sigma_xx = c11*dux_dxl + c15*(duz_dxl + dux_dzl) + c13*duz_dzl
     sigma_zz = c13*dux_dxl + c35*(duz_dxl + dux_dzl) + c33*duz_dzl
     sigma_xz = c15*dux_dxl + c55*(duz_dxl + dux_dzl) + c35*duz_dzl
  else
! no attenuation
    sigma_xx = lambdalplus2mul_relaxed*dux_dxl + lambdal_relaxed*duz_dzl
    sigma_xz = mul_relaxed*(duz_dxl + dux_dzl)
    sigma_zz = lambdalplus2mul_relaxed*duz_dzl + lambdal_relaxed*dux_dxl
  endif

    if(SIMULATION_TYPE == 2) then
    b_sigma_xx = lambdalplus2mul_relaxed*b_dux_dxl + lambdal_relaxed*b_duz_dzl
    b_sigma_xz = mul_relaxed*(b_duz_dxl + b_dux_dzl)
    b_sigma_zz = lambdalplus2mul_relaxed*b_duz_dzl + lambdal_relaxed*b_dux_dxl
    endif ! if(SIMULATION_TYPE == 2)

! get point values for the poroelastic side
          i = ivalue(ipoin1D,iedge_poroelastic)
          j = jvalue(ipoin1D,iedge_poroelastic)
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

          if(SIMULATION_TYPE == 2) then
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
            if(SIMULATION_TYPE == 2) then
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

          if(SIMULATION_TYPE == 2) then
          b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
          b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

          b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
          b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl

          b_dwx_dxl = b_dwx_dxi*xixl + b_dwx_dgamma*gammaxl
          b_dwx_dzl = b_dwx_dxi*xizl + b_dwx_dgamma*gammazl

          b_dwz_dxl = b_dwz_dxi*xixl + b_dwz_dgamma*gammaxl
          b_dwz_dzl = b_dwz_dxi*xizl + b_dwz_dgamma*gammazl
          endif
! compute stress tensor

! no attenuation
    sigma_xx = sigma_xx + lambdalplus2mul_G*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
    sigma_xz = sigma_xz + mul_G*(duz_dxl + dux_dzl)
    sigma_zz = sigma_zz + lambdalplus2mul_G*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

    sigmap = C_biot*(dux_dxl + duz_dzl) + M_biot*(dwx_dxl + dwz_dzl)

    if(SIMULATION_TYPE == 2) then
    b_sigma_xx = b_sigma_xx + lambdalplus2mul_G*b_dux_dxl + lambdal_G*b_duz_dzl + C_biot*(b_dwx_dxl + b_dwz_dzl)
    b_sigma_xz = b_sigma_xz + mul_G*(b_duz_dxl + b_dux_dzl)
    b_sigma_zz = b_sigma_zz + lambdalplus2mul_G*b_duz_dzl + lambdal_G*b_dux_dxl + C_biot*(b_dwx_dxl + b_dwz_dzl)
    b_sigmap = C_biot*(b_dux_dxl + b_duz_dzl) + M_biot*(b_dwx_dxl + b_dwz_dzl)
    endif

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
                weight*((sigma_xx*nx + sigma_xz*nz)/2.d0 -phil/tortl*sigmap*nx)

          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + &
                weight*((sigma_xz*nx + sigma_zz*nz)/2.d0 -phil/tortl*sigmap*nz)

! contribution to the fluid phase
! w = 0

          if(SIMULATION_TYPE == 2) then
! contribution to the solid phase
          b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) + &
                weight*((b_sigma_xx*nx + b_sigma_xz*nz)/2.d0 -phil/tortl*b_sigmap*nx)

          b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + &
                weight*((b_sigma_xz*nx + b_sigma_zz*nz)/2.d0 -phil/tortl*b_sigmap*nz)

! contribution to the fluid phase
! w = 0
          endif !if(SIMULATION_TYPE == 2) then

        enddo

      enddo

    endif ! if(coupled_elastic_poro)


! ************************************************************************************
! ******************************** add force source
! ************************************************************************************

 if(any_poroelastic) then


! --- add the source if it is a collocated force
    if(.not. initialfield) then

    do i_source=1,NSOURCES
! if this processor carries the source and the source element is elastic
      if (is_proc_source(i_source) == 1 .and. poroelastic(ispec_selected_source(i_source))) then

    phil = porosity(kmato(ispec_selected_source(i_source)))
    tortl = tortuosity(kmato(ispec_selected_source(i_source)))
    rhol_s = density(1,kmato(ispec_selected_source(i_source)))
    rhol_f = density(2,kmato(ispec_selected_source(i_source)))
    rhol_bar = (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f

! collocated force
        if(source_type(i_source) == 1) then
       if(SIMULATION_TYPE == 1) then  ! forward wavefield
          do j = 1,NGLLZ
           do i = 1,NGLLX
             iglob = ibool(i,j,ispec_selected_source(i_source))
             hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
! s
      accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - hlagrange * &
                               (1._CUSTOM_REAL - phil/tortl)*sin(angleforce(i_source))*source_time_function(i_source,it)
      accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + hlagrange * &
                               (1._CUSTOM_REAL - phil/tortl)*cos(angleforce(i_source))*source_time_function(i_source,it)
! w
      accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - hlagrange * &
         (1._CUSTOM_REAL - rhol_f/rhol_bar)*sin(angleforce(i_source))*source_time_function(i_source,it)
      accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) + hlagrange * &
         (1._CUSTOM_REAL - rhol_f/rhol_bar)*cos(angleforce(i_source))*source_time_function(i_source,it)
           enddo
          enddo
       else                   ! backward wavefield
          do j = 1,NGLLZ
           do i = 1,NGLLX
             iglob = ibool(i,j,ispec_selected_source(i_source))
             hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
! b_s
      b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - hlagrange * &
                               (1._CUSTOM_REAL - phil/tortl)*sin(angleforce(i_source))*source_time_function(i_source,NSTEP-it+1)
      b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + hlagrange * &
                               (1._CUSTOM_REAL - phil/tortl)*cos(angleforce(i_source))*source_time_function(i_source,NSTEP-it+1)
!b_w
      b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) - hlagrange * &
         (1._CUSTOM_REAL - rhol_f/rhol_bar)*sin(angleforce(i_source))*source_time_function(i_source,NSTEP-it+1)
      b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) + hlagrange * &
         (1._CUSTOM_REAL - rhol_f/rhol_bar)*cos(angleforce(i_source))*source_time_function(i_source,NSTEP-it+1)
           enddo
          enddo
       endif !endif SIMULATION_TYPE == 1
        endif

      endif ! if this processor carries the source and the source element is elastic
    enddo ! do i_source=1,NSOURCES

    endif ! if not using an initial field
  endif !if(any_poroelastic)

! assembling accels_proelastic & accelw_poroelastic for poroelastic elements
#ifdef USE_MPI
  if (nproc > 1 .and. any_poroelastic .and. ninterface_poroelastic > 0) then
    call assemble_MPI_vector_po(accels_poroelastic,accelw_poroelastic,npoin, &
      ninterface, ninterface_poroelastic,inum_interfaces_poroelastic, &
      max_interface_size, max_ibool_interfaces_size_po,&
      ibool_interfaces_poroelastic, nibool_interfaces_poroelastic, &
      tab_requests_send_recv_poro,buffer_send_faces_vector_pos,buffer_send_faces_vector_pow, &
      buffer_recv_faces_vector_pos,buffer_recv_faces_vector_pow, &
      my_neighbours)
  endif

  if (nproc > 1 .and. any_poroelastic .and. ninterface_poroelastic > 0 .and. SIMULATION_TYPE == 2) then
    call assemble_MPI_vector_po(b_accels_poroelastic,b_accelw_poroelastic,npoin, &
      ninterface, ninterface_poroelastic,inum_interfaces_poroelastic, &
      max_interface_size, max_ibool_interfaces_size_po,&
      ibool_interfaces_poroelastic, nibool_interfaces_poroelastic, &
      tab_requests_send_recv_poro,buffer_send_faces_vector_pos,buffer_send_faces_vector_pow, &
      buffer_recv_faces_vector_pos,buffer_recv_faces_vector_pow, &
      my_neighbours)
   endif
#endif


! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

 if(any_poroelastic) then
    accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
    accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
    velocs_poroelastic = velocs_poroelastic + deltatover2*accels_poroelastic

    accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
    accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)
    velocw_poroelastic = velocw_poroelastic + deltatover2*accelw_poroelastic

   if(SIMULATION_TYPE == 2) then
    b_accels_poroelastic(1,:) = b_accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
    b_accels_poroelastic(2,:) = b_accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
    b_velocs_poroelastic = b_velocs_poroelastic + b_deltatover2*b_accels_poroelastic

    b_accelw_poroelastic(1,:) = b_accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
    b_accelw_poroelastic(2,:) = b_accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)
    b_velocw_poroelastic = b_velocw_poroelastic + b_deltatover2*b_accelw_poroelastic
   endif

  endif !if(any_poroelastic)

!*******************************************************************************
!         assembling the displacements on the elastic-poro boundaries
!*******************************************************************************
    if(coupled_elastic_poro) then
     icount(:)=ZERO

! loop on all the coupling edges
      do inum = 1,num_solid_poro_edges
! get the edge of the elastic element
        ispec_elastic = solid_poro_elastic_ispec(inum)
        iedge_elastic = solid_poro_elastic_iedge(inum)
! get the corresponding edge of the poroelastic element
        ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

        do ipoin1D = 1,NGLLX
! recovering original velocities and accelerations on boundaries (elastic side)
          i = ivalue(ipoin1D,iedge_poroelastic)
          j = jvalue(ipoin1D,iedge_poroelastic)
          iglob = ibool(i,j,ispec_poroelastic)
          icount(iglob) = icount(iglob) + 1

        if(icount(iglob) ==1)then
          veloc_elastic(1,iglob) = veloc_elastic(1,iglob) - deltatover2*accel_elastic(1,iglob)
          veloc_elastic(3,iglob) = veloc_elastic(3,iglob) - deltatover2*accel_elastic(3,iglob)
          accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic(iglob)
          accel_elastic(3,iglob) = accel_elastic(3,iglob) / rmass_inverse_elastic(iglob)
! recovering original velocities and accelerations on boundaries (poro side)
          velocs_poroelastic(1,iglob) = velocs_poroelastic(1,iglob) - deltatover2*accels_poroelastic(1,iglob)
          velocs_poroelastic(2,iglob) = velocs_poroelastic(2,iglob) - deltatover2*accels_poroelastic(2,iglob)
          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)
! assembling accelerations
          accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
          accel_elastic(3,iglob) = ( accel_elastic(3,iglob) + accels_poroelastic(2,iglob) ) / &
                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
          accels_poroelastic(1,iglob) = accel_elastic(1,iglob)
          accels_poroelastic(2,iglob) = accel_elastic(3,iglob)
! updating velocities
          velocs_poroelastic(1,iglob) = velocs_poroelastic(1,iglob) + deltatover2*accels_poroelastic(1,iglob)
          velocs_poroelastic(2,iglob) = velocs_poroelastic(2,iglob) + deltatover2*accels_poroelastic(2,iglob)
          veloc_elastic(1,iglob) = veloc_elastic(1,iglob) + deltatover2*accel_elastic(1,iglob)
          veloc_elastic(3,iglob) = veloc_elastic(3,iglob) + deltatover2*accel_elastic(3,iglob)
! zeros w
          accelw_poroelastic(1,iglob) = ZERO
          accelw_poroelastic(2,iglob) = ZERO
          velocw_poroelastic(1,iglob) = ZERO
          velocw_poroelastic(2,iglob) = ZERO

         if(SIMULATION_TYPE == 2) then
          b_veloc_elastic(1,iglob) = b_veloc_elastic(1,iglob) - b_deltatover2*b_accel_elastic(1,iglob)
          b_veloc_elastic(3,iglob) = b_veloc_elastic(3,iglob) - b_deltatover2*b_accel_elastic(3,iglob)
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) / rmass_inverse_elastic(iglob)
          b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) / rmass_inverse_elastic(iglob)
! recovering original velocities and accelerations on boundaries (poro side)
          b_velocs_poroelastic(1,iglob) = b_velocs_poroelastic(1,iglob) - b_deltatover2*b_accels_poroelastic(1,iglob)
          b_velocs_poroelastic(2,iglob) = b_velocs_poroelastic(2,iglob) - b_deltatover2*b_accels_poroelastic(2,iglob)
          b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
          b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)
! assembling accelerations
          b_accel_elastic(1,iglob) = ( b_accel_elastic(1,iglob) + b_accels_poroelastic(1,iglob) ) / &
                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
          b_accel_elastic(3,iglob) = ( b_accel_elastic(3,iglob) + b_accels_poroelastic(2,iglob) ) / &
                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
          b_accels_poroelastic(1,iglob) = b_accel_elastic(1,iglob)
          b_accels_poroelastic(2,iglob) = b_accel_elastic(3,iglob)
! updating velocities
          b_velocs_poroelastic(1,iglob) = b_velocs_poroelastic(1,iglob) + b_deltatover2*b_accels_poroelastic(1,iglob)
          b_velocs_poroelastic(2,iglob) = b_velocs_poroelastic(2,iglob) + b_deltatover2*b_accels_poroelastic(2,iglob)
          b_veloc_elastic(1,iglob) = b_veloc_elastic(1,iglob) + b_deltatover2*b_accel_elastic(1,iglob)
          b_veloc_elastic(3,iglob) = b_veloc_elastic(3,iglob) + b_deltatover2*b_accel_elastic(3,iglob)
! zeros w
          b_accelw_poroelastic(1,iglob) = ZERO
          b_accelw_poroelastic(2,iglob) = ZERO
          b_velocw_poroelastic(1,iglob) = ZERO
          b_velocw_poroelastic(2,iglob) = ZERO
         endif !if(SIMULATION_TYPE == 2)

        endif !if(icount(iglob) ==1)

        enddo

      enddo
    endif

! ********************************************************************************************
!                       reading lastframe for adjoint/kernels calculation
! ********************************************************************************************
   if(it == 1 .and. SIMULATION_TYPE == 2) then

! acoustic medium
    if(any_acoustic) then
      write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
      open(unit=55,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
      do j=1,npoin
        read(55) b_potential_acoustic(j),&
                b_potential_dot_acoustic(j),&
                b_potential_dot_dot_acoustic(j)
        enddo
      close(55)

! free surface for an acoustic medium
      if ( nelem_acoustic_surface > 0 ) then
        call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                          b_potential_acoustic,acoustic_surface, &
                                          ibool,nelem_acoustic_surface,npoin,nspec)
      endif
    endif

! elastic medium
    if(any_elastic) then
    write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
      if(p_sv)then !P-SV waves
       do j=1,npoin
      read(55) (b_displ_elastic(i,j), i=1,NDIM), &
                  (b_veloc_elastic(i,j), i=1,NDIM), &
                  (b_accel_elastic(i,j), i=1,NDIM)
       enddo
       b_displ_elastic(3,:) = b_displ_elastic(2,:)
       b_displ_elastic(2,:) = ZERO
       b_veloc_elastic(3,:) = b_veloc_elastic(2,:)
       b_veloc_elastic(2,:) = ZERO
       b_accel_elastic(3,:) = b_accel_elastic(2,:)
       b_accel_elastic(2,:) = ZERO
      else !SH (membrane) waves
       do j=1,npoin
      read(55) b_displ_elastic(2,j), &
                  b_veloc_elastic(2,j), &
                  b_accel_elastic(2,j)
       enddo
       b_displ_elastic(1,:) = ZERO
       b_displ_elastic(3,:) = ZERO
       b_veloc_elastic(1,:) = ZERO
       b_veloc_elastic(3,:) = ZERO
       b_accel_elastic(1,:) = ZERO
       b_accel_elastic(3,:) = ZERO
      endif
    close(55)
    endif

! poroelastic medium
    if(any_poroelastic) then
    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_s',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_w',myrank,'.bin'
    open(unit=56,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
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
    endif

  endif ! if(it == 1 .and. SIMULATION_TYPE == 2)

! ********************************************************************************************
!                                      kernels calculation
! ********************************************************************************************
  if(any_elastic .and. SIMULATION_TYPE == 2) then ! kernels calculation
      do iglob = 1,npoin
            rho_k(iglob) =  accel_elastic(1,iglob)*b_displ_elastic(1,iglob) +&
                            accel_elastic(2,iglob)*b_displ_elastic(2,iglob) +&
                            accel_elastic(3,iglob)*b_displ_elastic(3,iglob)
            rhorho_el_hessian_temp1(iglob) = accel_elastic(1,iglob)*accel_elastic(1,iglob) +&
                                            accel_elastic(2,iglob)*accel_elastic(2,iglob)  +&
                                            accel_elastic(3,iglob)*accel_elastic(3,iglob)
            rhorho_el_hessian_temp2(iglob) = accel_elastic(1,iglob)*b_accel_elastic(1,iglob) +&
                                            accel_elastic(2,iglob)*b_accel_elastic(2,iglob)  +&
                                            accel_elastic(3,iglob)*b_accel_elastic(3,iglob)
      enddo
  endif

  if(any_poroelastic .and. SIMULATION_TYPE ==2) then
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
                        nspec,npoin,assign_external_model,it,deltat,t0,kmato,poroelastcoef,density, &
                        porosity,tortuosity, &
                        vpext,vsext,rhoext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext, &
                        anisotropic,anisotropy,wxgll,wzgll,numat, &
                        pressure_element,vector_field_element,e1,e11, &
                        potential_dot_acoustic,potential_dot_dot_acoustic, &
                        TURN_ATTENUATION_ON,Mu_nu1,Mu_nu2,N_SLS)

!----  display time step and max of norm of displacement
  if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
    call check_stability(myrank,time,it,NSTEP,npoin, &
                        any_elastic_glob,any_elastic,displ_elastic, &
                        any_poroelastic_glob,any_poroelastic, &
                        displs_poroelastic,displw_poroelastic, &
                        any_acoustic_glob,any_acoustic,potential_acoustic, &
                        year_start,month_start,time_start)
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
            numat,kmato,density,porosity,tortuosity,poroelastcoef,vpext,vsext,rhoext, &
            c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,anisotropic,anisotropy,ispec,e1,e11, &
            TURN_ATTENUATION_ON,Mu_nu1,Mu_nu2,N_SLS)

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
    valuy = ZERO
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
          dzd = vector_field_element(3,i,j)

        else if(seismotype == 6) then

          dxd = potential_acoustic(iglob)
          dzd = ZERO

        else if(seismotype == 1) then

             if(poroelastic(ispec)) then
          dxd = displs_poroelastic(1,iglob)
          dzd = displs_poroelastic(2,iglob)
             elseif(elastic(ispec)) then
          dxd = displ_elastic(1,iglob)
          dyd = displ_elastic(2,iglob)
          dzd = displ_elastic(3,iglob)
             endif

        else if(seismotype == 2) then

             if(poroelastic(ispec)) then
          dxd = velocs_poroelastic(1,iglob)
          dzd = velocs_poroelastic(2,iglob)
             elseif(elastic(ispec)) then
          dxd = veloc_elastic(1,iglob)
          dyd = veloc_elastic(2,iglob)
          dzd = veloc_elastic(3,iglob)
             endif

        else if(seismotype == 3) then

             if(poroelastic(ispec)) then
          dxd = accels_poroelastic(1,iglob)
          dzd = accels_poroelastic(2,iglob)
             elseif(elastic(ispec)) then
          dxd = accel_elastic(1,iglob)
          dyd = accel_elastic(2,iglob)
          dzd = accel_elastic(3,iglob)
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
        if(elastic(ispec))  valuy = valuy + dyd*hlagrange
        valuz = valuz + dzd*hlagrange
        valcurl = valcurl + dcurld*hlagrange

      enddo
   enddo

! rotate seismogram components if needed, except if recording pressure, which is a scalar
    if(seismotype /= 4 .and. seismotype /= 6) then
      if(p_sv) then
        sisux(seismo_current,irecloc) =   cosrot_irec(irecloc)*valux + sinrot_irec(irecloc)*valuz
        sisuz(seismo_current,irecloc) = - sinrot_irec(irecloc)*valux + cosrot_irec(irecloc)*valuz
      else
        sisux(seismo_current,irecloc) = valuy
        sisuz(seismo_current,irecloc) = ZERO
      endif
    else
      sisux(seismo_current,irecloc) = valux
      sisuz(seismo_current,irecloc) = ZERO
    endif
    siscurl(seismo_current,irecloc) = valcurl

 enddo


!----- writing the kernels
!
! kernels output
  if(SIMULATION_TYPE == 2) then

   if(any_acoustic) then

    do ispec = 1, nspec
     if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then
      do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
    kappal_ac_global(iglob) = poroelastcoef(3,1,kmato(ispec))
    rhol_ac_global(iglob) = density(1,kmato(ispec))

! calcul the displacement by computing the gradient of potential / rho
! and calcul the acceleration by computing the gradient of potential_dot_dot / rho
        tempx1l = ZERO
        tempx2l = ZERO
        b_tempx1l = ZERO
        b_tempx2l = ZERO
        do k = 1,NGLLX
! derivative along x
          tempx1l = tempx1l + potential_dot_dot_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
          b_tempx1l = b_tempx1l + b_potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
          bb_tempx1l = bb_tempx1l + b_potential_dot_dot_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
! derivative along z
          tempx2l = tempx2l + potential_dot_dot_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
          b_tempx2l = b_tempx2l + b_potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
          bb_tempx2l = bb_tempx2l + b_potential_dot_dot_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        if(assign_external_model) rhol_ac_global(iglob) = rhoext(i,j,ispec)

! derivatives of potential
        accel_ac(1,iglob) = (tempx1l*xixl + tempx2l*gammaxl) / rhol_ac_global(iglob)
        accel_ac(2,iglob) = (tempx1l*xizl + tempx2l*gammazl) / rhol_ac_global(iglob)
        b_displ_ac(1,iglob) = (b_tempx1l*xixl + b_tempx2l*gammaxl) / rhol_ac_global(iglob)
        b_displ_ac(2,iglob) = (b_tempx1l*xizl + b_tempx2l*gammazl) / rhol_ac_global(iglob)
        b_accel_ac(1,iglob) = (bb_tempx1l*xixl + bb_tempx2l*gammaxl) / rhol_ac_global(iglob)
        b_accel_ac(2,iglob) = (bb_tempx1l*xizl + bb_tempx2l*gammazl) / rhol_ac_global(iglob)

          enddo !i = 1, NGLLX
      enddo !j = 1, NGLLZ
     endif
    enddo

          do ispec = 1,nspec
     if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then
      do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) - rhol_ac_global(iglob)  * &
                           dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * deltat
            kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) - kappal_ac_global(iglob) * &
                           potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * &
                           b_potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob)&
                           * deltat
!
            rhop_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + kappa_ac_kl(i,j,ispec)
            alpha_ac_kl(i,j,ispec) = TWO *  kappa_ac_kl(i,j,ispec)
            rhorho_ac_hessian_final1(i,j,ispec) =  rhorho_ac_hessian_final1(i,j,ispec) + &
                             dot_product(accel_ac(:,iglob),accel_ac(:,iglob)) * deltat
            rhorho_ac_hessian_final2(i,j,ispec) =  rhorho_ac_hessian_final2(i,j,ispec) + &
                             dot_product(accel_ac(:,iglob),b_accel_ac(:,iglob)) * deltat
         enddo
       enddo
      endif
          enddo

    endif !if(any_acoustic)

   if(any_elastic) then

    do ispec = 1, nspec
     if(elastic(ispec)) then
      do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
    mul_global(iglob) = poroelastcoef(2,1,kmato(ispec))
    kappal_global(iglob) = poroelastcoef(3,1,kmato(ispec)) - 4._CUSTOM_REAL*mul_global(iglob)/3._CUSTOM_REAL
    rhol_global(iglob) = density(1,kmato(ispec))

            rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - rhol_global(iglob)  * rho_k(iglob) * deltat
            mu_kl(i,j,ispec) =  mu_kl(i,j,ispec) - TWO * mul_global(iglob) * mu_k(iglob) * deltat
            kappa_kl(i,j,ispec) = kappa_kl(i,j,ispec) - kappal_global(iglob) * kappa_k(iglob) * deltat
!
            rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + kappa_kl(i,j,ispec) + mu_kl(i,j,ispec)
            beta_kl(i,j,ispec) = TWO * (mu_kl(i,j,ispec) - 4._CUSTOM_REAL * mul_global(iglob) &
                  / (3._CUSTOM_REAL * kappal_global(iglob)) * kappa_kl(i,j,ispec))
            alpha_kl(i,j,ispec) = TWO * (1._CUSTOM_REAL + 4._CUSTOM_REAL * mul_global(iglob)/&
                   (3._CUSTOM_REAL * kappal_global(iglob))) * kappa_kl(i,j,ispec)
            rhorho_el_hessian_final1(i,j,ispec) = rhorho_el_hessian_final1(i,j,ispec) + rhorho_el_hessian_temp1(iglob) * deltat
            rhorho_el_hessian_final2(i,j,ispec) = rhorho_el_hessian_final2(i,j,ispec) + rhorho_el_hessian_temp2(iglob) * deltat

          enddo
      enddo
     endif
    enddo

   endif !if(any_elastic)

  if(any_poroelastic) then

    do ispec = 1, nspec
     if(poroelastic(ispec)) then
      do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
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

            rhot_kl(i,j,ispec) = rhot_kl(i,j,ispec) - deltat * rhol_bar_global(iglob) * rhot_k(iglob)
            rhof_kl(i,j,ispec) = rhof_kl(i,j,ispec) - deltat * rhol_f_global(iglob) * rhof_k(iglob)
            sm_kl(i,j,ispec) = sm_kl(i,j,ispec) - deltat * rhol_f_global(iglob)*tortl_global(iglob)/phil_global(iglob) * sm_k(iglob)
!at the moment works with constant permeability
            eta_kl(i,j,ispec) = eta_kl(i,j,ispec) - deltat * etal_f_global(iglob)/permlxx_global(iglob) * eta_k(iglob)
            B_kl(i,j,ispec) = B_kl(i,j,ispec) - deltat * B_k(iglob)
            C_kl(i,j,ispec) = C_kl(i,j,ispec) - deltat * C_k(iglob)
            M_kl(i,j,ispec) = M_kl(i,j,ispec) - deltat * M_k(iglob)
            mufr_kl(i,j,ispec) = mufr_kl(i,j,ispec) - TWO * deltat * mufr_k(iglob)
! density kernels
            rholb = rhol_bar_global(iglob) - phil_global(iglob)*rhol_f_global(iglob)/tortl_global(iglob)
            rhob_kl(i,j,ispec) = rhot_kl(i,j,ispec) + B_kl(i,j,ispec) + mufr_kl(i,j,ispec)
            rhofb_kl(i,j,ispec) = rhof_kl(i,j,ispec) + C_kl(i,j,ispec) + M_kl(i,j,ispec) + sm_kl(i,j,ispec)
            Bb_kl(i,j,ispec) = B_kl(i,j,ispec)
            Cb_kl(i,j,ispec) = C_kl(i,j,ispec)
            Mb_kl(i,j,ispec) = M_kl(i,j,ispec)
            mufrb_kl(i,j,ispec) = mufr_kl(i,j,ispec)
            phi_kl(i,j,ispec) = - sm_kl(i,j,ispec) - M_kl(i,j,ispec)
! wave speed kernels
            dd1 = (1._CUSTOM_REAL+rholb/rhol_f_global(iglob))*ratio**2 + 2._CUSTOM_REAL*ratio +&
                tortl_global(iglob)/phil_global(iglob)
            rhobb_kl(i,j,ispec) = rhob_kl(i,j,ispec) - &
                phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*B_biot) * &
                   (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1+&
                   (rhol_bar_global(iglob)**2*ratio**2/rhol_f_global(iglob)**2*(phil_global(iglob)/tortl_global(iglob)*&
                   ratio+1)*(phil_global(iglob)/tortl_global(iglob)*ratio+phil_global(iglob)/tortl_global(iglob)*&
                   (1+rhol_f_global(iglob)/rhol_bar_global(iglob))-1) )/dd1**2 )- FOUR_THIRDS*cssquare )*&
                   Bb_kl(i,j,ispec) - &
                rhol_bar_global(iglob)*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                   (phil_global(iglob)/tortl_global(iglob)*ratio + 1._CUSTOM_REAL)**2/dd1**2*Mb_kl(i,j,ispec) + &
                rhol_bar_global(iglob)*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                   phil_global(iglob)*ratio/tortl_global(iglob)*(phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                   (1+rhol_bar_global(iglob)*ratio/rhol_f_global(iglob))/dd1**2)*Cb_kl(i,j,ispec)+ &
                phil_global(iglob)*rhol_f_global(iglob)*cssquare/(tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
           rhofbb_kl(i,j,ispec) = rhofb_kl(i,j,ispec) + &
                phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*B_biot) * &
                   (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1+&
                   (rhol_bar_global(iglob)**2*ratio**2/rhol_f_global(iglob)**2*(phil_global(iglob)/tortl_global(iglob)*&
                   ratio+1)*(phil_global(iglob)/tortl_global(iglob)*ratio+phil_global(iglob)/tortl_global(iglob)*&
                   (1+rhol_f_global(iglob)/rhol_bar_global(iglob))-1) )/dd1**2 )- FOUR_THIRDS*cssquare )*&
                   Bb_kl(i,j,ispec) + &
                rhol_bar_global(iglob)*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                   (phil_global(iglob)/tortl_global(iglob)*ratio + 1._CUSTOM_REAL)**2/dd1**2*Mb_kl(i,j,ispec) - &
                rhol_bar_global(iglob)*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                   phil_global(iglob)*ratio/tortl_global(iglob)*(phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                   (1+rhol_bar_global(iglob)*ratio/rhol_f_global(iglob))/dd1**2)*Cb_kl(i,j,ispec)- &
                phil_global(iglob)*rhol_f_global(iglob)*cssquare/(tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
           phib_kl(i,j,ispec) = phi_kl(i,j,ispec) - &
                phil_global(iglob)*rhol_bar_global(iglob)/(tortl_global(iglob)*B_biot) * ( cpIsquare - rhol_f_global(iglob)/&
                   rhol_bar_global(iglob)*cpIIsquare- &
                   (cpIsquare-cpIIsquare)*( (TWO*ratio**2*phil_global(iglob)/tortl_global(iglob) + (1._CUSTOM_REAL+&
                   rhol_f_global(iglob)/rhol_bar_global(iglob))*(TWO*ratio*phil_global(iglob)/tortl_global(iglob)+&
                   1._CUSTOM_REAL))/dd1 + (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(phil_global(iglob)*&
                   ratio/tortl_global(iglob)+phil_global(iglob)/tortl_global(iglob)*(1._CUSTOM_REAL+rhol_f_global(iglob)/&
                   rhol_bar_global(iglob))-1._CUSTOM_REAL)*((1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-&
                   TWO*phil_global(iglob)/tortl_global(iglob))*ratio**2+TWO*ratio)/dd1**2 ) - &
                   FOUR_THIRDS*rhol_f_global(iglob)*cssquare/rhol_bar_global(iglob) )*Bb_kl(i,j,ispec) + &
                rhol_f_global(iglob)/M_biot * (cpIsquare-cpIIsquare)*(&
                   TWO*ratio*(phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2*((1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                   rhol_f_global(iglob)-TWO*phil_global(iglob)/tortl_global(iglob))*ratio**2+TWO*ratio)/dd1**2&
                   )*Mb_kl(i,j,ispec) + &
                phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*C_biot)*(cpIsquare-cpIIsquare)*ratio* (&
                   (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob)*ratio)/dd1 - &
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                   rhol_f_global(iglob)*ratio)*((1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-TWO*&
                   phil_global(iglob)/tortl_global(iglob))*ratio+TWO)/dd1**2&
                    )*Cb_kl(i,j,ispec) -&
                phil_global(iglob)*rhol_f_global(iglob)*cssquare/(tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
           cpI_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIsquare/B_biot*rhol_bar_global(iglob)*( &
                   1._CUSTOM_REAL-phil_global(iglob)/tortl_global(iglob) + &
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*&
                   ratio+phil_global(iglob)/tortl_global(iglob)*(1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-&
                   1._CUSTOM_REAL)/dd1 &
                    )* Bb_kl(i,j,ispec) +&
                2._CUSTOM_REAL*cpIsquare*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot) *&
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2/dd1*Mb_kl(i,j,ispec)+&
                2._CUSTOM_REAL*cpIsquare*rhol_f_global(iglob)/C_biot * &
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                   rhol_f_global(iglob)*ratio)/dd1*Cb_kl(i,j,ispec)
           cpII_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIIsquare*rhol_bar_global(iglob)/B_biot * (&
                   phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*rhol_bar_global(iglob)) - &
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*&
                   ratio+phil_global(iglob)/tortl_global(iglob)*(1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-&
                   1._CUSTOM_REAL)/dd1  ) * Bb_kl(i,j,ispec) +&
                2._CUSTOM_REAL*cpIIsquare*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot) * (&
                   1._CUSTOM_REAL - (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2/dd1  )*Mb_kl(i,j,ispec) + &
                2._CUSTOM_REAL*cpIIsquare*rhol_f_global(iglob)/C_biot * (&
                   1._CUSTOM_REAL - (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(1._CUSTOM_REAL+&
                   rhol_bar_global(iglob)/rhol_f_global(iglob)*ratio)/dd1  )*Cb_kl(i,j,ispec)
           cs_kl(i,j,ispec) = - 8._CUSTOM_REAL/3._CUSTOM_REAL*cssquare*rhol_bar_global(iglob)/B_biot*(1._CUSTOM_REAL-&
                   phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*rhol_bar_global(iglob)))*Bb_kl(i,j,ispec) + &
                2._CUSTOM_REAL*(rhol_bar_global(iglob)-rhol_f_global(iglob)*phil_global(iglob)/tortl_global(iglob))/&
                   mulfr_global(iglob)*cssquare*mufrb_kl(i,j,ispec)
           ratio_kl(i,j,ispec) = ratio*rhol_bar_global(iglob)*phil_global(iglob)/(tortl_global(iglob)*B_biot) * &
                   (cpIsquare-cpIIsquare) * ( &
                   phil_global(iglob)/tortl_global(iglob)*(2._CUSTOM_REAL*ratio+1._CUSTOM_REAL+rhol_f_global(iglob)/ &
                   rhol_bar_global(iglob))/dd1 - (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                   (phil_global(iglob)/tortl_global(iglob)*ratio+phil_global(iglob)/tortl_global(iglob)*(&
                   1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-1._CUSTOM_REAL)*(2._CUSTOM_REAL*ratio*(&
                   1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob)) +&
                   2._CUSTOM_REAL)/dd1**2  )*Bb_kl(i,j,ispec) + &
                ratio*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot)*(cpIsquare-cpIIsquare) * &
                   2._CUSTOM_REAL*phil_global(iglob)/tortl_global(iglob) * (&
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                   (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2*((1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                   rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob))*ratio+1._CUSTOM_REAL)/dd1**2 )*Mb_kl(i,j,ispec) +&
                ratio*rhol_f_global(iglob)/C_biot*(cpIsquare-cpIIsquare) * (&
                   (2._CUSTOM_REAL*phil_global(iglob)*rhol_bar_global(iglob)*ratio/(tortl_global(iglob)*rhol_f_global(iglob))+&
                   phil_global(iglob)/tortl_global(iglob)+rhol_bar_global(iglob)/rhol_f_global(iglob))/dd1 - &
                   2._CUSTOM_REAL*phil_global(iglob)/tortl_global(iglob)*(phil_global(iglob)/tortl_global(iglob)*ratio+&
                   1._CUSTOM_REAL)*(1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)*ratio)*((1._CUSTOM_REAL+&
                   rhol_bar_global(iglob)/rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob))*ratio+1._CUSTOM_REAL)/&
                   dd1**2 )*Cb_kl(i,j,ispec)

          enddo
       enddo
     endif
    enddo

   endif ! if(any_poroelastic)

   endif ! if(SIMULATION_TYPE == 2)

!
!----  display results at given time steps
!
  if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then

!
! kernels output files
!

   if(SIMULATION_TYPE == 2 .and. it == NSTEP) then

  if ( myrank == 0 ) then
  write(IOUT,*) 'Writing Kernels file'
  endif

    if(any_acoustic) then
    do ispec = 1, nspec
      do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
        xx = coord(1,iglob)
        zz = coord(2,iglob)
         write(95,'(5e11.3)')xx,zz,rho_ac_kl(i,j,ispec),kappa_ac_kl(i,j,ispec)
         write(96,'(5e11.3)')rhorho_ac_hessian_final1(i,j,ispec), rhorho_ac_hessian_final2(i,j,ispec),&
                             rhop_ac_kl(i,j,ispec),alpha_ac_kl(i,j,ispec)
          enddo
      enddo
    enddo
    close(95)
    close(96)
    endif

    if(any_elastic) then
    do ispec = 1, nspec
      do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
        xx = coord(1,iglob)
        zz = coord(2,iglob)
         write(97,'(5e11.3)')xx,zz,rho_kl(i,j,ispec),kappa_kl(i,j,ispec),mu_kl(i,j,ispec)
         write(98,'(5e11.3)')xx,zz,rhop_kl(i,j,ispec),alpha_kl(i,j,ispec),beta_kl(i,j,ispec)
         !write(98,'(5e11.3)')rhorho_el_hessian_final1(i,j,ispec), rhorho_el_hessian_final2(i,j,ispec),&
         !                    rhop_kl(i,j,ispec),alpha_kl(i,j,ispec),beta_kl(i,j,ispec)
          enddo
      enddo
    enddo
    close(97)
    close(98)
    endif

    if(any_poroelastic) then
    do ispec = 1, nspec
      do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
        xx = coord(1,iglob)
        zz = coord(2,iglob)
         write(144,'(5e11.3)')xx,zz,mufr_kl(i,j,ispec),B_kl(i,j,ispec),C_kl(i,j,ispec)
         write(155,'(5e11.3)')xx,zz,M_kl(i,j,ispec),rhot_kl(i,j,ispec),rhof_kl(i,j,ispec)
         write(16,'(5e11.3)')xx,zz,sm_kl(i,j,ispec),eta_kl(i,j,ispec)
         write(17,'(5e11.3)')xx,zz,mufrb_kl(i,j,ispec),Bb_kl(i,j,ispec),Cb_kl(i,j,ispec)
         write(18,'(5e11.3)')xx,zz,Mb_kl(i,j,ispec),rhob_kl(i,j,ispec),rhofb_kl(i,j,ispec)
         write(19,'(5e11.3)')xx,zz,phi_kl(i,j,ispec),eta_kl(i,j,ispec)
         write(20,'(5e11.3)')xx,zz,cpI_kl(i,j,ispec),cpII_kl(i,j,ispec),cs_kl(i,j,ispec)
         write(21,'(5e11.3)')xx,zz,rhobb_kl(i,j,ispec),rhofbb_kl(i,j,ispec),ratio_kl(i,j,ispec)
         write(22,'(5e11.3)')xx,zz,phib_kl(i,j,ispec),eta_kl(i,j,ispec)
          enddo
      enddo
    enddo
    close(144)
    close(155)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)
    endif

    endif

!
!----  PostScript display
!
  if(output_postscript_snapshot) then

  if (myrank == 0) write(IOUT,*) 'Writing PostScript file'

  if(imagetype == 1 .and. p_sv) then

    if (myrank == 0) write(IOUT,*) 'drawing displacement vector as small arrows...'

    call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic,&
          elastic,poroelastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,numat,kmato,density,rhoext,assign_external_model)

    call plotpost(vector_field_display,coord,vpext,x_source,z_source,x_final_receiver,z_final_receiver, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,porosity,tortuosity,&
          poroelastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,nelem_acoustic_surface,acoustic_edges, &
          simulation_title,npoin,npgeo,vpImin,vpImax,nrec,NSOURCES, &
          colors,numbers,subsamp,imagetype,interpol,meshvect,modelvect, &
          boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,numat,pointsdisp, &
          nspec,ngnod,coupled_acoustic_elastic,coupled_acoustic_poro,coupled_elastic_poro, &
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

  else if(imagetype == 2 .and. p_sv) then

    if (myrank == 0) write(IOUT,*) 'drawing velocity vector as small arrows...'

    call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic,&
          elastic,poroelastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,numat,kmato,density,rhoext,assign_external_model)

    call plotpost(vector_field_display,coord,vpext,x_source,z_source,x_final_receiver,z_final_receiver, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,porosity,tortuosity,&
          poroelastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,nelem_acoustic_surface,acoustic_edges, &
          simulation_title,npoin,npgeo,vpImin,vpImax,nrec,NSOURCES, &
          colors,numbers,subsamp,imagetype,interpol,meshvect,modelvect, &
          boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,numat,pointsdisp, &
          nspec,ngnod,coupled_acoustic_elastic,coupled_acoustic_poro,coupled_elastic_poro, &
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

  else if(imagetype == 3 .and. p_sv) then

    if (myrank == 0) write(IOUT,*) 'drawing acceleration vector as small arrows...'

    call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,accels_poroelastic,&
          elastic,poroelastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,numat,kmato,density,rhoext,assign_external_model)

    call plotpost(vector_field_display,coord,vpext,x_source,z_source,x_final_receiver,z_final_receiver, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,porosity,tortuosity,&
          poroelastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,nelem_acoustic_surface,acoustic_edges, &
          simulation_title,npoin,npgeo,vpImin,vpImax,nrec,NSOURCES, &
          colors,numbers,subsamp,imagetype,interpol,meshvect,modelvect, &
          boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,numat,pointsdisp, &
          nspec,ngnod,coupled_acoustic_elastic,coupled_acoustic_poro,coupled_elastic_poro, &
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

  else if(imagetype == 4 .or. .not. p_sv) then

    if (myrank == 0) write(IOUT,*) 'cannot draw scalar pressure field or y-component field as a vector plot, skipping...'

  else
    call exit_MPI('wrong type for snapshots')
  endif

  if (myrank == 0 .and. imagetype /= 4 .and. p_sv) write(IOUT,*) 'PostScript file written'

  endif

!
!----  display color image
!
  if(output_color_image) then

  if (myrank == 0) write(IOUT,*) 'Creating color image of size ',NX_IMAGE_color,' x ',NZ_IMAGE_color,' for time step ',it

  if(imagetype == 1) then

    if (myrank == 0) write(IOUT,*) 'drawing image of z (if P-SV) or y (if SH) component of displacement vector...'

    call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic,&
          elastic,poroelastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,numat,kmato,density,rhoext,assign_external_model)

  else if(imagetype == 2) then

    if (myrank == 0) write(IOUT,*) 'drawing image of z (if P-SV) or y (if SH) component of velocity vector...'

    call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic,&
          elastic,poroelastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,numat,kmato,density,rhoext,assign_external_model)

  else if(imagetype == 3) then

    if (myrank == 0) write(IOUT,*) 'drawing image of z (if P-SV) or y (if SH) component of acceleration vector...'

    call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,accels_poroelastic,&
          elastic,poroelastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,numat,kmato,density,rhoext,assign_external_model)

  else if(imagetype == 4 .and. p_sv) then

    if (myrank == 0) write(IOUT,*) 'drawing image of pressure field...'

    call compute_pressure_whole_medium(potential_dot_dot_acoustic,displ_elastic,&
         displs_poroelastic,displw_poroelastic,elastic,poroelastic,vector_field_display, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,assign_external_model, &
         numat,kmato,density,porosity,tortuosity,poroelastcoef,vpext,vsext,rhoext, &
         c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,anisotropic,anisotropy,e1,e11, &
         TURN_ATTENUATION_ON,Mu_nu1,Mu_nu2,N_SLS)

  else if(imagetype == 4 .and. .not. p_sv) then
    call exit_MPI('cannot draw pressure field for SH (membrane) waves')
  else
    call exit_MPI('wrong type for snapshots')
  endif

  image_color_data(:,:) = 0.d0

  do k = 1, nb_pixel_loc
     j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
     i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color
    if(p_sv) then !P-SH waves, plot vertical component or pressure
     image_color_data(i,j) = vector_field_display(3,iglob_image_color(i,j))
    else !SH (membrane) waves, plot y-component
     image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))
    endif
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
           if(p_sv) then !P-SH waves, plot vertical component or pressure
             data_pixel_send(k) = vector_field_display(3,iglob_image_color(i,j))
           else !SH (membrane) waves, plot y-component
             data_pixel_send(k) = vector_field_display(2,iglob_image_color(i,j))
           endif
        enddo

        call MPI_SEND(data_pixel_send(1),nb_pixel_loc,MPI_DOUBLE_PRECISION, 0, 43, MPI_COMM_WORLD, ier)

     endif
  endif

#endif

  if (myrank == 0) then
     call create_color_image(image_color_data,iglob_image_color, &
                NX_IMAGE_color,NZ_IMAGE_color,it,cutsnaps,image_color_vp_display)
     write(IOUT,*) 'Color image created'
  endif

  endif

!----  save temporary or final seismograms
! suppress seismograms if we generate traces of the run for analysis with "ParaVer", because time consuming
  if(.not. GENERATE_PARAVER_TRACES) &
    call write_seismograms(sisux,sisuz,siscurl,station_name,network_name,NSTEP, &
                          nrecloc,which_proc_receiver,nrec,myrank,deltat,seismotype,st_xval,t0, &
                          NTSTEP_BETWEEN_OUTPUT_SEISMO,seismo_offset,seismo_current,p_sv)

  seismo_offset = seismo_offset + seismo_current
  seismo_current = 0

  endif

#ifdef USE_MPI
! add a barrier if we generate traces of the run for analysis with "ParaVer"
  if(GENERATE_PARAVER_TRACES) call MPI_BARRIER(MPI_COMM_WORLD,ier)
#endif

  enddo ! end of the main time loop

  if((SAVE_FORWARD .and. SIMULATION_TYPE==1) .or. SIMULATION_TYPE ==2) then
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
  if(SAVE_FORWARD .and. SIMULATION_TYPE ==1 .and. any_elastic) then
  if ( myrank == 0 ) then
    write(IOUT,*)
    write(IOUT,*) 'Saving elastic last frame...'
    write(IOUT,*)
  endif
    write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
      if(p_sv)then !P-SV waves
       do j=1,npoin
      write(55) displ_elastic(1,j), displ_elastic(3,j), &
                  veloc_elastic(1,j), veloc_elastic(3,j), &
                  accel_elastic(1,j), accel_elastic(3,j)
       enddo
      else !SH (membrane) waves
       do j=1,npoin
      write(55) displ_elastic(2,j), &
                  veloc_elastic(2,j), &
                  accel_elastic(2,j)
       enddo
      endif
    close(55)
  endif

  if(SAVE_FORWARD .and. SIMULATION_TYPE ==1 .and. any_poroelastic) then
  if ( myrank == 0 ) then
    write(IOUT,*)
    write(IOUT,*) 'Saving poroelastic last frame...'
    write(IOUT,*)
  endif
    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_s',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_w',myrank,'.bin'
    open(unit=56,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
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

  if(SAVE_FORWARD .and. SIMULATION_TYPE ==1 .and. any_acoustic) then
  if ( myrank == 0 ) then
    write(IOUT,*)
    write(IOUT,*) 'Saving acoustic last frame...'
    write(IOUT,*)
  endif
    write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
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
    close(IOUT_ENERGY)
    open(unit=IOUT_ENERGY,file='plotenergy',status='unknown')
    write(IOUT_ENERGY,*) 'set term postscript landscape color solid "Helvetica" 22'
    write(IOUT_ENERGY,*) 'set output "energy.ps"'
    write(IOUT_ENERGY,*) 'set xlabel "Time (s)"'
    write(IOUT_ENERGY,*) 'set ylabel "Energy (J)"'
    write(IOUT_ENERGY,*) 'plot "energy.gnu" us 1:4 t ''Total Energy'' w l 1, "energy.gnu" us 1:3 t ''Potential Energy'' w l 2'
    close(IOUT_ENERGY)
  endif

   if (.not. any_poroelastic) then
open(unit=1001,file='DATA/model_velocity.dat_output',status='unknown')
   if ( .NOT. assign_external_model) then
allocate(rho_local(ngllx,ngllz,nspec)); rho_local=0.
allocate(vp_local(ngllx,ngllz,nspec)); vp_local=0.
allocate(vs_local(ngllx,ngllz,nspec)); vs_local=0.
!!      write(1001,*) npoin
!!      do iglob = 1,npoin
!!         write(1001,*) coord(1,iglob),coord(2,iglob),rho_global(iglob),vp_global(iglob),vs_global(iglob)
!!      end do
    do ispec = 1,nspec
       do j = 1,NGLLZ
       do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          rho_local(i,j,ispec) = density(1,kmato(ispec))
          vp_local(i,j,ispec) = sqrt(poroelastcoef(3,1,kmato(ispec))/density(1,kmato(ispec)))
          vs_local(i,j,ispec) = sqrt(poroelastcoef(2,1,kmato(ispec))/density(1,kmato(ispec)))
          write(1001,'(I10, 5F13.4)') iglob, coord(1,iglob),coord(2,iglob),&
                                      rho_local(i,j,ispec),vp_local(i,j,ispec),vs_local(i,j,ispec)
       end do
       end do
    end do
   else
!!     write(1001,*) npoin
!!  do iglob = 1,npoin
!!     write(1001,*) coord(1,iglob),coord(2,iglob),rhoext_global(iglob),vpext_global(iglob),vsext_global(iglob)
!!  end do
     do ispec = 1,nspec
        do j = 1,NGLLZ
        do i = 1,NGLLX
           iglob = ibool(i,j,ispec)
           write(1001,'(I10,5F13.4)') iglob, coord(1,iglob),coord(2,iglob),&
                                       rhoext(i,j,ispec),vpext(i,j,ispec),vsext(i,j,ispec)
        end do
        end do
     end do
   endif
close(1001)
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

  end program specfem2D

