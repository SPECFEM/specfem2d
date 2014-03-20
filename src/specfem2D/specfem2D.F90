
  program specfem2D

!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, Inria and University of Pau, France,
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
! @ARTICLE{PeKoLuMaLeCaLeMaLiBlNiBaTr11,
! author = {Daniel Peter and Dimitri Komatitsch and Yang Luo and Roland Martin
!     and Nicolas {Le Goff} and Emanuele Casarotti and Pieyre {Le Loher}
!     and Federica Magnoni and Qinya Liu and C\'eline Blitz and Tarje Nissen-Meyer
!     and Piero Basini and Jeroen Tromp},
! title = {Forward and adjoint simulations of seismic wave propagation on fully
!     unstructured hexahedral meshes},
! journal={Geophys. J. Int.},
! year = {2011},
! volume = {186},
! pages = {721-739},
! number = {2},
! doi = {10.1111/j.1365-246X.2011.05044.x}}
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
! @ARTICLE{PeKoLuMaLeCaLeMaLiBlNiBaTr11,
! author = {Daniel Peter and Dimitri Komatitsch and Yang Luo and Roland Martin
!     and Nicolas {Le Goff} and Emanuele Casarotti and Pieyre {Le Loher}
!     and Federica Magnoni and Qinya Liu and C\'eline Blitz and Tarje Nissen-Meyer
!     and Piero Basini and Jeroen Tromp},
! title = {Forward and adjoint simulations of seismic wave propagation on fully
!     unstructured hexahedral meshes},
! journal={Geophys. J. Int.},
! year = {2011},
! volume = {186},
! pages = {721-739},
! number = {2},
! doi = {10.1111/j.1365-246X.2011.05044.x}}
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
! version 7.0, Dimitri Komatitsch, Zhinan Xie, Paul Cristini, Roland Martin and Rene Matzen, July 2012:
!               - added support for Convolution PML absorbing layers
!               - added higher-order time schemes (4th order Runge-Kutta and LDDRK4-6)
!               - many small or moderate bug fixes
!
! version 6.2, many developers, April 2011:
!               - restructured package source code into separate src/ directories
!               - added configure & Makefile scripts and a PDF manual in doc/
!               - added user examples in EXAMPLES/
!               - added a USER_T0 parameter to fix the onset time in simulation
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

#ifdef USE_MPI
  use :: mpi
#endif

  implicit none

  include "constants.h"
#ifdef USE_MPI
  include "precision.h"
#endif

!! DK DK uncomment this in order to force vectorization of the loops
!! DK DK using a trick that goes out of the array bounds
!! DK DK (then array bound checking cannot be used, thus for instance do NOT use -check all in Intel ifort)
! #define FORCE_VECTORIZATION

  integer NSOURCES,i_source
  integer, dimension(:), allocatable :: source_type,time_function_type
  double precision, dimension(:), allocatable :: x_source,z_source,xi_source,gamma_source,&
                  Mxx,Mzz,Mxz,f0,tshift_src,factor,anglesource
  integer, dimension(:), allocatable :: ix_image_color_source,iy_image_color_source
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: sourcearray
  double precision :: t0

  double precision, dimension(:,:), allocatable :: coorg

! for P-SV or SH (membrane) waves calculation
  logical :: p_sv

! factor to subsample color images output by the code (useful for very large models)
  double precision :: factor_subsample_image

! use snapshot number in the file name of JPG color snapshots instead of the time step
  logical :: USE_SNAPSHOT_NUMBER_IN_FILENAME

! display acoustic layers as constant blue, because they likely correspond to water in the case of ocean acoustics
! or in the case of offshore oil industry experiments.
! (if off, display them as greyscale, as for elastic or poroelastic elements,
!  for instance for acoustic-only oil industry models of solid media)
  logical :: DRAW_WATER_IN_BLUE

! US letter paper or European A4
  logical :: US_LETTER

! non linear display to enhance small amplitudes in color images
  double precision :: POWER_DISPLAY_COLOR

! output seismograms in Seismic Unix format (adjoint traces will be read in the same format)
  logical :: SU_FORMAT

! use this t0 as earliest starting time rather than the automatically calculated one
! (must be positive and bigger than the automatically one to be effective;
!  simulation will start at t = - t0)
  double precision :: USER_T0

! value of time_stepping_scheme to decide which time scheme will be used
! 1 = Newmark (2nd order), 2 = LDDRK4-6 (4th-order 6-stage low storage Runge-Kutta)
! 3 = classical 4th-order 4-stage Runge-Kutta
  integer :: time_stepping_scheme

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

  integer :: i,j,k,it,irec,id,n,ispec,ispec2,nglob,npgeo,iglob
  integer :: nglob_acoustic
  integer :: nglob_elastic
  integer :: nglob_poroelastic
  logical :: anyabs
  double precision :: dxd,dyd,dzd,dcurld,valux,valuy,valuz,valcurl,hlagrange,xi,gamma,x,z

! add a small crack (discontinuity) in the medium manually
  logical, parameter :: ADD_A_SMALL_CRACK_IN_THE_MEDIUM = .false.
!! must be set equal to the number of spectral elements on one vertical side of the crack
  integer :: NB_POINTS_TO_ADD_TO_NPGEO = 3
  integer :: check_nb_points_to_add_to_npgeo,current_last_point,npgeo_ori,original_value,ignod
  logical :: already_found_a_crack_element

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
  double precision :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: veloc_elastic_LDDRK,displ_elastic_LDDRK,&
                                                         veloc_elastic_LDDRK_temp
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: accel_elastic_rk,veloc_elastic_rk
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: veloc_elastic_initial_rk,displ_elastic_initial_rk
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_elastic_adj_coupling
  double precision, dimension(:,:), allocatable :: &
    coord, flagrange,xinterp,zinterp,Uxinterp,Uzinterp,vector_field_display

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

! poroelastic and elastic coefficients
  double precision, dimension(:,:,:), allocatable :: poroelastcoef

! anisotropy parameters
  logical :: all_anisotropic
  double precision ::  c11,c13,c15,c33,c35,c55,c12,c23,c25
  logical, dimension(:), allocatable :: anisotropic
  double precision, dimension(:,:), allocatable :: anisotropy

! for acoustic medium
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic,potential_acoustic_old
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_dot_acoustic_LDDRK, potential_acoustic_LDDRK
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_dot_acoustic_temp
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic_init_rk, potential_dot_acoustic_init_rk
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: potential_dot_dot_acoustic_rk, potential_dot_acoustic_rk
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_acoustic_adj_coupling

! inverse mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_elastic_one
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_elastic_three

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_acoustic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic

! to evaluate cpI, cpII, and cs, and rI (poroelastic medium)
  double precision :: rhol_s,rhol_f,rhol_bar,phil,tortl
  double precision :: mul_s,kappal_s
  double precision :: kappal_f
  double precision :: mul_fr,kappal_fr
  double precision :: D_biot,H_biot,C_biot,M_biot,B_biot,cpIsquare,cpIIsquare,cssquare
  real(kind=CUSTOM_REAL) :: ratio,dd1

  double precision, dimension(:,:,:), allocatable :: vpext,vsext,rhoext
  double precision, dimension(:,:,:), allocatable :: QKappa_attenuationext,Qmu_attenuationext
  double precision, dimension(:,:,:), allocatable :: c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext

  double precision, dimension(:,:,:), allocatable :: shape2D,shape2D_display
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: xix,xiz,gammax,gammaz,jacobian

  double precision, dimension(:,:,:,:), allocatable :: dershape2D,dershape2D_display

  integer, dimension(:,:,:), allocatable :: ibool,ibool_outer,ibool_inner
  integer, dimension(:,:), allocatable  :: knods
  integer, dimension(:), allocatable :: kmato,numabs, &
     ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3,ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2

  integer, dimension(:), allocatable :: ispec_selected_source,iglob_source,&
                                        is_proc_source,nb_proc_source
  double precision, dimension(:), allocatable :: aval
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: source_time_function
  double precision, external :: netlib_specfun_erf

  double precision :: vpImin,vpImax,vpIImin,vpIImax

  real(kind=CUSTOM_REAL) :: kinetic_energy,potential_energy,kinetic_energy_total,potential_energy_total

  integer :: colors,numbers,subsamp_postscript,imagetype_postscript, &
    NSTEP_BETWEEN_OUTPUT_INFO,seismotype,NSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP_BETWEEN_OUTPUT_IMAGES, &
    NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS,subsamp_seismos,imagetype_JPEG,imagetype_wavefield_dumps
  integer :: numat,ngnod,nspec,pointsdisp, &
    nelemabs,nelem_acoustic_surface,ispecabs,UPPER_LIMIT_DISPLAY,NELEM_PML_THICKNESS

  logical interpol,meshvect,modelvect,boundvect,assign_external_model,initialfield, &
    output_grid_ASCII,output_grid_Gnuplot,ATTENUATION_VISCOELASTIC_SOLID,output_postscript_snapshot,output_color_image, &
    plot_lowerleft_corner_only,add_Bielak_conditions,output_energy,READ_EXTERNAL_SEP_FILE, &
    output_wavefield_dumps,use_binary_for_wavefield_dumps,PML_BOUNDARY_CONDITIONS,ROTATE_PML_ACTIVATE,STACEY_BOUNDARY_CONDITIONS

  double precision :: ROTATE_PML_ANGLE

! AXISYM parameters                                                                                                  !axisym
                                                                                                                     !axisym
  logical :: AXISYM ! .true. if we are performing a 2.5D simulation                                                  !axisym
  ! Nember of elements on the symmetry axis                                                                          !axisym
  integer :: nelem_on_the_axis                                                                                       !axisym
  ! Flag to know if an element is on the axis                                                                        !axisym
  logical, dimension(:), allocatable :: is_on_the_axis                                                               !axisym
  integer, dimension(:), allocatable  ::ispec_of_axial_elements                                                      !axisym
  ! Gauss-Lobatto-Jacobi points and weights                                                                          !axisym
  double precision, dimension(NGLJ) :: xiglj                                                                         !axisym
  real(kind=CUSTOM_REAL), dimension(NGLJ) :: wxglj                                                                   !axisym
  ! derivatives of GLJ polynomials                                                                                   !axisym
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLJ) :: hprimeBar_xx,hprimeBarwglj_xx                                      !axisym
  ! Shape functions (and their derivatives) evaluated at the GLJ points                                              !axisym
  double precision, dimension(:,:), allocatable :: flagrange_GLJ                                                     !axisym
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1                                                         !axisym

! for CPML_element_file
  logical :: read_external_mesh

  double precision :: cutsnaps,sizemax_arrows,anglerec,xirec,gammarec

! for absorbing and acoustic free surface conditions
  integer :: ispec_acoustic_surface,inum
  real(kind=CUSTOM_REAL) :: nx,nz,weight,xxi,zgamma

  logical, dimension(:,:), allocatable  :: codeabs
  integer, dimension(:), allocatable  :: typeabs

! for attenuation
  integer  :: N_SLS
  double precision, dimension(:), allocatable  :: QKappa_attenuation
  double precision, dimension(:), allocatable  :: Qmu_attenuation
  double precision  :: f0_attenuation
  integer nspec_allocate

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: e1,e11,e13
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: e1_LDDRK,e11_LDDRK,e13_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: e1_initial_rk,e11_initial_rk,e13_initial_rk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: e1_force_rk,e11_force_rk,e13_force_rk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: inv_tau_sigma_nu1_sent,phi_nu1_sent,inv_tau_sigma_nu2_sent,phi_nu2_sent
  real(kind=CUSTOM_REAL), dimension(:,:,:) , allocatable :: Mu_nu1,Mu_nu2
  real(kind=CUSTOM_REAL) :: Mu_nu1_sent,Mu_nu2_sent

! for viscous attenuation
  double precision, dimension(:,:,:), allocatable :: &
    rx_viscous,rz_viscous,viscox,viscoz

  double precision, dimension(:,:,:), allocatable :: &
    rx_viscous_LDDRK,rz_viscous_LDDRK

  double precision, dimension(:,:,:), allocatable :: &
    rx_viscous_initial_rk,rz_viscous_initial_rk

  double precision, dimension(:,:,:,:), allocatable :: &
    rx_viscous_force_RK,rz_viscous_force_RK

  double precision :: theta_e,theta_s
  double precision :: Q0,freq0
  double precision :: alphaval,betaval,gammaval,thetainv
  logical :: ATTENUATION_PORO_FLUID_PART,save_ASCII_seismograms,save_binary_seismograms_single,save_binary_seismograms_double, &
             DRAW_SOURCES_AND_RECEIVERS, save_ASCII_kernels
  double precision, dimension(NGLLX,NGLLZ) :: viscox_loc,viscoz_loc
  double precision :: Sn,Snp1,etal_f
  double precision, dimension(3):: bl_unrelaxed_elastic
  double precision :: permlxx,permlxz,permlzz,invpermlxx,invpermlxz,invpermlzz,detk
! adjoint
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: b_viscodampx,b_viscodampz
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
  integer, dimension(:), allocatable :: ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro,&
            iend_edge3_poro,ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro
  logical :: any_solid_poro_edges

! for adjoint method
  logical :: SAVE_FORWARD ! whether or not the last frame is saved to reconstruct the forward field
  integer :: SIMULATION_TYPE      ! 1 = forward wavefield, 3 = backward and adjoint wavefields and kernels
  double precision :: b_deltatover2,b_deltatsquareover2,b_deltat ! coefficients of the explicit Newmark time scheme
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    b_accels_poroelastic,b_velocs_poroelastic,b_displs_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    b_accelw_poroelastic,b_velocw_poroelastic,b_displw_poroelastic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    b_accel_elastic,b_veloc_elastic,b_displ_elastic,b_displ_elastic_old
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic,b_potential_acoustic_old
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

! for color images
  integer :: NX_IMAGE_color,NZ_IMAGE_color
  double precision :: xmin_color_image,xmax_color_image, &
    zmin_color_image,zmax_color_image
  integer, dimension(:,:), allocatable :: iglob_image_color,copy_iglob_image_color
  double precision, dimension(:,:), allocatable :: image_color_data
  double precision, dimension(:,:), allocatable :: image_color_vp_display
  integer  :: nb_pixel_loc
  integer, dimension(:), allocatable  :: num_pixel_loc

! name of wavefield snapshot file
  character(len=150) :: wavefield_file

#ifdef USE_MPI
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

  ! to determine date and time at which the run will finish
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  integer :: year,mon,day,hr,minutes,timestamp
  double precision :: timestamp_seconds_start

! for MPI and partitioning
  integer  :: ier
  integer  :: nproc,nproc_read_from_database
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

! to compute analytical initial plane wave field
  double precision :: anglesource_refl, c_inc, c_refl, cploc, csloc
  double precision, dimension(2) :: A_plane, B_plane, C_plane
  double precision :: time_offset

! beyond critical angle
  integer , dimension(:), allocatable :: left_bound,right_bound,bot_bound
  double precision , dimension(:,:), allocatable :: v0x_left,v0z_left,v0x_right,v0z_right,v0x_bot,v0z_bot
  double precision , dimension(:,:), allocatable :: t0x_left,t0z_left,t0x_right,t0z_right,t0x_bot,t0z_bot
  integer count_left,count_right,count_bottom
  logical :: over_critical_angle

! inner/outer elements in the case of an MPI simulation
  integer :: ispec_inner,ispec_outer
  integer :: nglob_outer,nglob_inner

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
  integer, dimension(:), allocatable :: ix_image_color_receiver,iy_image_color_receiver
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

! to dump the wave field
  integer :: icounter,nb_of_values_to_save
  logical :: this_is_the_first_time_we_dump
  logical, dimension(:), allocatable  :: mask_ibool,mask_ibool_pml

! to create a sorted version of the indirect addressing to reduce cache misses
  integer, dimension(:,:,:), allocatable :: copy_ibool_ori
  integer, dimension(:), allocatable :: integer_mask_ibool

  double precision, dimension(:,:,:),allocatable:: rho_local,vp_local,vs_local
!!!! hessian
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhorho_el_hessian_final1, rhorho_el_hessian_final2
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rhorho_el_hessian_temp1, rhorho_el_hessian_temp2
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: rhorho_ac_hessian_final1, rhorho_ac_hessian_final2

! to help locate elements with a negative Jacobian using OpenDX
  logical :: found_a_negative_jacobian

! add spring to Stacey absorbing boundary condition
  logical :: ADD_SPRING_TO_STACEY
  double precision :: x_center_spring,z_center_spring
  double precision :: xmin,xmax,zmin,zmax
  double precision :: xmin_local,xmax_local,zmin_local,zmax_local

! for horizontal periodic conditions
  logical :: ADD_PERIODIC_CONDITIONS
  logical, dimension(:), allocatable :: this_ibool_is_a_periodic_edge

! horizontal periodicity distance for periodic conditions
  double precision :: PERIODIC_HORIZ_DIST

  double precision :: xmaxval,xminval,ymaxval,yminval,xtol,xtypdist
  integer :: counter

  integer :: isnapshot_number = 0

!<SU_FORMAT
  integer(kind=4) :: r4head(60)
  character(len=512) :: filename
  real(kind=4),dimension(:,:),allocatable :: adj_src_s
  integer(kind=2) :: header2(2)
!>SU_FORMAT

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

! for LDDRK46
  integer :: i_stage,stage_time_scheme
  real(kind=CUSTOM_REAL), dimension(Nstages):: alpha_LDDRK,beta_LDDRK,c_LDDRK

  ! parameters used in LDDRK scheme, from equation (2) of
  ! Berland, J., Bogey, C., & Bailly, C.
  ! Low-dissipation and low-dispersion fourth-order Runge–Kutta algorithm, Computers & Fluids, 35(10), 1459-1463.
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

! to count the number of degrees of freedom
  integer :: count_nspec_acoustic,count_nspec_acoustic_total,nspec_total,nglob_total,nb_acoustic_DOFs,nb_elastic_DOFs
  double precision :: ratio_1DOF,ratio_2DOFs

! PML parameters
  logical, dimension(:), allocatable :: is_PML
  integer, dimension(:), allocatable :: region_CPML
  double precision, dimension(:,:,:), allocatable :: &
                    K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
   rmemory_dux_dx,rmemory_duz_dx,rmemory_dux_dz,rmemory_duz_dz

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_fsb_displ_elastic

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
   rmemory_dux_dx_prime,rmemory_duz_dx_prime,rmemory_dux_dz_prime,rmemory_duz_dz_prime

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
   rmemory_dux_dx_LDDRK,rmemory_duz_dx_LDDRK,rmemory_dux_dz_LDDRK,rmemory_duz_dz_LDDRK

  integer, dimension(:), allocatable :: spec_to_PML
  logical, dimension(:,:), allocatable :: which_PML_elem

! defined for using PML in elastic simulation
  logical, dimension(:,:), allocatable :: PML_interior_interface
  integer, dimension(:), allocatable :: point_interface
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: pml_interface_history_displ
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: pml_interface_history_veloc
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: pml_interface_history_accel
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: pml_interface_history_potential
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: pml_interface_history_potential_dot
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: pml_interface_history_potential_dot_dot
  integer :: nglob_interface

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_displ_elastic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rmemory_displ_elastic_LDDRK

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_potential_acoustic
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                          rmemory_acoustic_dux_dx,rmemory_acoustic_dux_dz

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rmemory_potential_acoustic_LDDRK
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                          rmemory_acoustic_dux_dx_LDDRK,rmemory_acoustic_dux_dz_LDDRK

  logical :: anyabs_glob
  integer :: nspec_PML

!!  logical :: backward_simulation for seprate adjoint simulation from backward simulation

!***********************************************************************
!
!             i n i t i a l i z a t i o n    p h a s e
!
!***********************************************************************

  ! force Flush-To-Zero if available to avoid very slow Gradual Underflow trapping
  call force_ftz()

  call initialize_simulation(nproc,myrank,ninterface_acoustic,ninterface_elastic,ninterface_poroelastic)
  if(nproc < 1) stop 'should have nproc >= 1'

  ! starts reading in Database file
  call read_databases_init(myrank, &
                  simulation_title,AXISYM,SIMULATION_TYPE,NOISE_TOMOGRAPHY,SAVE_FORWARD,npgeo,nproc_read_from_database, &
                  output_grid_Gnuplot,interpol,NSTEP_BETWEEN_OUTPUT_INFO,NSTEP_BETWEEN_OUTPUT_SEISMOS, &
                  NSTEP_BETWEEN_OUTPUT_IMAGES,PML_BOUNDARY_CONDITIONS,ROTATE_PML_ACTIVATE,ROTATE_PML_ANGLE,NELEM_PML_THICKNESS, &
                  NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS,subsamp_seismos, &
                  imagetype_JPEG,imagetype_wavefield_dumps, &
                  output_postscript_snapshot,output_color_image,colors,numbers, &
                  meshvect,modelvect,boundvect,cutsnaps,subsamp_postscript,sizemax_arrows, &
                  anglerec,initialfield,add_Bielak_conditions, &
                  seismotype,imagetype_postscript,assign_external_model,READ_EXTERNAL_SEP_FILE, &
                  output_grid_ASCII,output_energy,output_wavefield_dumps,use_binary_for_wavefield_dumps, &
                  ATTENUATION_VISCOELASTIC_SOLID,ATTENUATION_PORO_FLUID_PART,save_ASCII_seismograms, &
                  save_binary_seismograms_single,save_binary_seismograms_double,DRAW_SOURCES_AND_RECEIVERS, &
                  Q0,freq0,p_sv,NSTEP,deltat,NSOURCES, &
                  factor_subsample_image,USE_SNAPSHOT_NUMBER_IN_FILENAME,DRAW_WATER_IN_BLUE,US_LETTER, &
                  POWER_DISPLAY_COLOR,SU_FORMAT,USER_T0, time_stepping_scheme, &
                  ADD_SPRING_TO_STACEY,ADD_PERIODIC_CONDITIONS,PERIODIC_HORIZ_DIST, &
                  read_external_mesh,save_ASCII_kernels)

  if(nproc_read_from_database < 1) stop 'should have nproc_read_from_database >= 1'
  if(SIMULATION_TYPE == 3 .and.(time_stepping_scheme == 2 .or. time_stepping_scheme == 3)) &
                                  stop 'RK and LDDRK time scheme not supported for adjoint inversion'
  if(nproc /= nproc_read_from_database) stop 'must always have nproc == nproc_read_from_database'

! add a small crack (discontinuity) in the medium manually
  npgeo_ori = npgeo
  if(ADD_A_SMALL_CRACK_IN_THE_MEDIUM) npgeo = npgeo + NB_POINTS_TO_ADD_TO_NPGEO

  !
  !--- source information
  !
    allocate( source_type(NSOURCES) )
    allocate( time_function_type(NSOURCES) )
    allocate( x_source(NSOURCES) )
    allocate( z_source(NSOURCES) )
    allocate( ix_image_color_source(NSOURCES) )
    allocate( iy_image_color_source(NSOURCES) )
    allocate( f0(NSOURCES) )
    allocate( tshift_src(NSOURCES) )
    allocate( factor(NSOURCES) )
    allocate( anglesource(NSOURCES) )
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

  ! reads in source infos
  call read_databases_sources(NSOURCES,source_type,time_function_type, &
                      x_source,z_source,Mxx,Mzz,Mxz,f0,tshift_src,factor,anglesource)

  !if(AXISYM) factor = factor/(TWO*PI)                                                                         !axisym TODO verify

  ! sets source parameters
  call set_sources(myrank,NSOURCES,source_type,time_function_type, &
                      x_source,z_source,Mxx,Mzz,Mxz,f0,tshift_src,factor,anglesource,aval, &
                      t0,initialfield,deltat,USER_T0)

  !----  define time stepping scheme
  if(time_stepping_scheme == 1)then
    stage_time_scheme=1
  else if(time_stepping_scheme == 2)then
    stage_time_scheme=Nstages
  else if(time_stepping_scheme == 3)then
    stage_time_scheme=4
  endif

  !----  read attenuation information
  call read_databases_atten(N_SLS,f0_attenuation)

  ! if source is not a Dirac or Heavyside then f0_attenuation is f0 of the first source
  if(.not. (time_function_type(1) == 4 .or. time_function_type(1) == 5)) then
    f0_attenuation = f0(1)
  endif

  !---- read the spectral macrobloc nodal coordinates
  allocate(coorg(NDIM,npgeo))

  ! reads the spectral macrobloc nodal coordinates
  ! and basic properties of the spectral elements
  !! DK DK  call read_databases_coorg_elem(myrank,npgeo,coorg,numat,ngnod,nspec, &
  !! DK DK  added a crack manually
  call read_databases_coorg_elem(myrank,npgeo_ori,coorg,numat,ngnod,nspec, &
                              pointsdisp,plot_lowerleft_corner_only, &
                              nelemabs,nelem_acoustic_surface, &
                              num_fluid_solid_edges,num_fluid_poro_edges, &
                              num_solid_poro_edges,nnodes_tangential_curve, &
                              nelem_on_the_axis)                                                                     !axisym

  !---- allocate arrays
    allocate(shape2D(ngnod,NGLLX,NGLLZ))
    allocate(dershape2D(NDIM,ngnod,NGLLX,NGLLZ))
    allocate(shape2D_display(ngnod,pointsdisp,pointsdisp))
    allocate(dershape2D_display(NDIM,ngnod,pointsdisp,pointsdisp))
    if(AXISYM) then                                                                  !axisym
      allocate(flagrange_GLJ(NGLJ,pointsdisp))                                       !axisym
    else                                                                             !axisym
      allocate(flagrange_GLJ(1,1))                                                   !axisym
    endif                                                                            !axisym
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
    allocate(anisotropy(9,numat))
    allocate(porosity(numat))
    allocate(tortuosity(numat))
    allocate(permeability(3,numat))
    allocate(poroelastcoef(4,3,numat))
    allocate(QKappa_attenuation(numat))
    allocate(Qmu_attenuation(numat))
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

  !
  !---- read the material properties
  !
  call gmat01(density,porosity,tortuosity,anisotropy,permeability,poroelastcoef,numat,&
              myrank,QKappa_attenuation,Qmu_attenuation,freq0,Q0,f0(1),ATTENUATION_PORO_FLUID_PART)
  !
  !----  read spectral macrobloc data
  !

! add support for using PML in MPI mode with external mesh
  allocate(region_CPML(nspec))
  call read_databases_mato(nspec,ngnod,kmato,knods,region_CPML)

! add a small crack (discontinuity) in the medium manually
  if(ADD_A_SMALL_CRACK_IN_THE_MEDIUM) then

#ifdef USE_MPI
  stop 'currently only serial runs are handled when adding a crack manually'
#endif
!! DK DK material number 2 indicates the spectral elements that form the left vertical side of the crack
  check_nb_points_to_add_to_npgeo = count(kmato == 2)
  print *
  print *,'adding a crack manually'
  print *,'need to add ',nb_points_to_add_to_npgeo,' npgeo mesh points to do that'

  if(check_nb_points_to_add_to_npgeo /= NB_POINTS_TO_ADD_TO_NPGEO) &
    stop 'must have check_nb_points_to_add_to_npgeo == NB_POINTS_TO_ADD_TO_NPGEO when adding a crack manually'

  if(ngnod /= 4) stop 'must currently have ngnod == 4 when adding a crack manually'

  if(FAST_NUMBERING) stop 'must not have FAST_NUMBERING when adding a crack manually'

!! DK DK modify arrays "knods" and "coorg" to introduce the crack manually by duplicating and splitting the nodes
  already_found_a_crack_element = .false.
  current_last_point = npgeo_ori

  do ispec = 1,nspec-1
!! DK DK my convention is to introduce a vertical crack between two elements with material numbers 2 and 3
    if(kmato(ispec) == 2 .and. kmato(ispec+1) == 3) then

      print *,'adding a crack between elements ',ispec,' and ',ispec+1

!! DK DK duplicate and split the lower-right corner of this element,
!! DK DK except if it is the first crack element found, because then it is the crack
!! DK DK tip and thus it should be assembled rather than split.
!! DK DK Lower-right corner of an element is local npgeo point #2
      if(already_found_a_crack_element .and. knods(2,ispec) <= npgeo_ori) then
        current_last_point = current_last_point + 1
        original_value = knods(2,ispec)
!! DK DK split this point number in all the elements in which it appears
        do ispec2 = 1,nspec
! do this only for elements that define the left vertical edge of the crack
          if(kmato(ispec2) /= 2) cycle
          do ignod = 1,ngnod
            if(knods(ignod,ispec2) == original_value) then
              knods(ignod,ispec2) = current_last_point
              coorg(:,current_last_point) = coorg(:,original_value)
            endif
          enddo
        enddo
      endif

!! DK DK duplicate and split the upper-right corner of this element
      already_found_a_crack_element = .true.

!! DK DK Upper-right corner of an element is local npgeo point #3
      if(knods(3,ispec) <= npgeo_ori) then

        current_last_point = current_last_point + 1
        original_value = knods(3,ispec)
!! DK DK split this point number in all the elements in which it appears
        do ispec2 = 1,nspec
! do this only for elements that define the left vertical edge of the crack
          if(kmato(ispec2) /= 2) cycle
          do ignod = 1,ngnod
            if(knods(ignod,ispec2) == original_value) then
              knods(ignod,ispec2) = current_last_point
              coorg(:,current_last_point) = coorg(:,original_value)
            endif
          enddo
        enddo
      endif

    endif ! of if(kmato(ispec) == 2 .and. kmato(ispec+1) == 3)

  enddo

  if(current_last_point /= npgeo) then
    print *,'current_last_point = ',current_last_point
    print *,'npgeo_new = ',npgeo
    stop 'did not find the right total number of points, should have current_last_point == npgeo_new'
  endif

  endif ! of if(ADD_A_SMALL_CRACK_IN_THE_MEDIUM) then

!-------------------------------------------------------------------------------
!----  determine if each spectral element is elastic, poroelastic, or acoustic
!-------------------------------------------------------------------------------
  call initialize_simulation_domains(any_acoustic,any_elastic,any_poroelastic, &
                                anisotropic,elastic,poroelastic,porosity,anisotropy,kmato,numat, &
                                nspec,nspec_allocate,p_sv,ATTENUATION_VISCOELASTIC_SOLID,count_nspec_acoustic)

  if(PML_BOUNDARY_CONDITIONS .and. any_poroelastic) then
    stop 'PML boundary conditions not implemented for poroelastic simulations yet'
  endif

  if(PML_BOUNDARY_CONDITIONS .and. any_elastic .and. (.not. p_sv)) then
    stop 'PML boundary conditions not implemented for SH simulations yet'
  endif

  if(PML_BOUNDARY_CONDITIONS .and. time_stepping_scheme == 3) then
    stop 'PML boundary conditions not implemented with standard Runge Kutta scheme'
  endif

#ifdef USE_MPI
  if(myrank == 0)then
   if(time_stepping_scheme == 3) then
    stop 'mpi support for standard Runge Kutta scheme is not implemented'
   endif
  endif
#endif

  ! allocate memory variables for attenuation
    allocate(e1(NGLLX,NGLLZ,nspec_allocate,N_SLS))
    allocate(e11(NGLLX,NGLLZ,nspec_allocate,N_SLS))
    allocate(e13(NGLLX,NGLLZ,nspec_allocate,N_SLS))

    e1(:,:,:,:) = 0._CUSTOM_REAL
    e11(:,:,:,:) = 0._CUSTOM_REAL
    e13(:,:,:,:) = 0._CUSTOM_REAL

    if(time_stepping_scheme == 2)then
      allocate(e1_LDDRK(NGLLX,NGLLZ,nspec_allocate,N_SLS))
      allocate(e11_LDDRK(NGLLX,NGLLZ,nspec_allocate,N_SLS))
      allocate(e13_LDDRK(NGLLX,NGLLZ,nspec_allocate,N_SLS))
    else
      allocate(e1_LDDRK(1,1,1,1))
      allocate(e11_LDDRK(1,1,1,1))
      allocate(e13_LDDRK(1,1,1,1))
    endif
    e1_LDDRK(:,:,:,:) = 0._CUSTOM_REAL
    e11_LDDRK(:,:,:,:) = 0._CUSTOM_REAL
    e13_LDDRK(:,:,:,:) = 0._CUSTOM_REAL

    if(time_stepping_scheme == 3)then
      allocate(e1_initial_rk(NGLLX,NGLLZ,nspec_allocate,N_SLS))
      allocate(e11_initial_rk(NGLLX,NGLLZ,nspec_allocate,N_SLS))
      allocate(e13_initial_rk(NGLLX,NGLLZ,nspec_allocate,N_SLS))
      allocate(e1_force_rk(NGLLX,NGLLZ,nspec_allocate,N_SLS,stage_time_scheme))
      allocate(e11_force_rk(NGLLX,NGLLZ,nspec_allocate,N_SLS,stage_time_scheme))
      allocate(e13_force_rk(NGLLX,NGLLZ,nspec_allocate,N_SLS,stage_time_scheme))
    else
      allocate(e1_initial_rk(1,1,1,1))
      allocate(e11_initial_rk(1,1,1,1))
      allocate(e13_initial_rk(1,1,1,1))
      allocate(e1_force_rk(1,1,1,1,1))
      allocate(e11_force_rk(1,1,1,1,1))
      allocate(e13_force_rk(1,1,1,1,1))
    endif
    e1_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
    e11_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
    e13_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
    e1_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL
    e11_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL
    e13_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL
    allocate(Mu_nu1(NGLLX,NGLLZ,nspec))
    allocate(Mu_nu2(NGLLX,NGLLZ,nspec))

! define the attenuation quality factors.
! they can be different for each element.
!! DK DK if needed in the future, here the quality factor could be different for each point
  do ispec = 1,nspec
    call attenuation_model(N_SLS,QKappa_attenuation(kmato(ispec)),Qmu_attenuation(kmato(ispec)), &
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
    if(ATTENUATION_PORO_FLUID_PART) then
      allocate(rx_viscous(NGLLX,NGLLZ,nspec))
      allocate(rz_viscous(NGLLX,NGLLZ,nspec))
      allocate(viscox(NGLLX,NGLLZ,nspec))
      allocate(viscoz(NGLLX,NGLLZ,nspec))

      if(time_stepping_scheme == 2) then
      allocate(rx_viscous_LDDRK(NGLLX,NGLLZ,nspec))
      allocate(rz_viscous_LDDRK(NGLLX,NGLLZ,nspec))
      endif

      if(time_stepping_scheme == 3) then
      allocate(rx_viscous_initial_rk(NGLLX,NGLLZ,nspec))
      allocate(rz_viscous_initial_rk(NGLLX,NGLLZ,nspec))
      allocate(rx_viscous_force_RK(NGLLX,NGLLZ,nspec,stage_time_scheme))
      allocate(rz_viscous_force_RK(NGLLX,NGLLZ,nspec,stage_time_scheme))
      endif

    else
      allocate(rx_viscous(NGLLX,NGLLZ,1))
      allocate(rz_viscous(NGLLX,NGLLZ,1))
      allocate(viscox(NGLLX,NGLLZ,1))
      allocate(viscoz(NGLLX,NGLLZ,1))
    endif

  !
  !----  read interfaces data
  !
  call read_databases_ninterface(ninterface,max_interface_size)
  if ( ninterface > 0 ) then
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
   call read_databases_interfaces(ninterface,max_interface_size,my_neighbours,my_nelmnts_neighbours,my_interfaces)

  else
       allocate(my_neighbours(1))
       allocate(my_nelmnts_neighbours(1))
       allocate(my_interfaces(1,1,1))
       allocate(ibool_interfaces_acoustic(1,1))
       allocate(ibool_interfaces_elastic(1,1))
       allocate(ibool_interfaces_poroelastic(1,1))
       allocate(nibool_interfaces_acoustic(1))
       allocate(nibool_interfaces_elastic(1))
       allocate(nibool_interfaces_poroelastic(1))
       allocate(inum_interfaces_acoustic(1))
       allocate(inum_interfaces_elastic(1))
       allocate(inum_interfaces_poroelastic(1))
  endif


! --- allocate arrays for absorbing boundary conditions

  if(nelemabs <= 0) then
    nelemabs = 1
    anyabs = .false.
  else
    anyabs = .true.
  endif

    allocate(numabs(nelemabs))
    allocate(codeabs(4,nelemabs))
    allocate(typeabs(nelemabs))

    allocate(ibegin_edge1(nelemabs))
    allocate(iend_edge1(nelemabs))
    allocate(ibegin_edge3(nelemabs))
    allocate(iend_edge3(nelemabs))

    allocate(ibegin_edge4(nelemabs))
    allocate(iend_edge4(nelemabs))
    allocate(ibegin_edge2(nelemabs))
    allocate(iend_edge2(nelemabs))

    allocate(ibegin_edge1_poro(nelemabs))
    allocate(iend_edge1_poro(nelemabs))
    allocate(ibegin_edge3_poro(nelemabs))
    allocate(iend_edge3_poro(nelemabs))

    allocate(ibegin_edge4_poro(nelemabs))
    allocate(iend_edge4_poro(nelemabs))
    allocate(ibegin_edge2_poro(nelemabs))
    allocate(iend_edge2_poro(nelemabs))

    allocate(ib_left(nelemabs))
    allocate(ib_right(nelemabs))
    allocate(ib_bottom(nelemabs))
    allocate(ib_top(nelemabs))

  !
  !----  read absorbing boundary data
  !
  call read_databases_absorbing(myrank,nelemabs,nspec,anyabs, &
                            ibegin_edge1,iend_edge1,ibegin_edge2,iend_edge2, &
                            ibegin_edge3,iend_edge3,ibegin_edge4,iend_edge4, &
                            numabs,codeabs,typeabs, &
                            nspec_left,nspec_right,nspec_bottom,nspec_top, &
                            ib_right,ib_left,ib_bottom,ib_top,PML_BOUNDARY_CONDITIONS)

  if(anyabs .and. (.not. PML_BOUNDARY_CONDITIONS))then
    STACEY_BOUNDARY_CONDITIONS = .true.
  else
    STACEY_BOUNDARY_CONDITIONS = .false.
  endif


  if( anyabs ) then
    ! files to save absorbed waves needed to reconstruct backward wavefield for adjoint method
      if(any_elastic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3).and. (.not. PML_BOUNDARY_CONDITIONS)) then
        allocate(b_absorb_elastic_left(3,NGLLZ,nspec_left,NSTEP))
        allocate(b_absorb_elastic_right(3,NGLLZ,nspec_right,NSTEP))
        allocate(b_absorb_elastic_bottom(3,NGLLX,nspec_bottom,NSTEP))
        allocate(b_absorb_elastic_top(3,NGLLX,nspec_top,NSTEP))
      else
        allocate(b_absorb_elastic_left(1,1,1,1))
        allocate(b_absorb_elastic_right(1,1,1,1))
        allocate(b_absorb_elastic_bottom(1,1,1,1))
        allocate(b_absorb_elastic_top(1,1,1,1))
      endif
      if(any_poroelastic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3).and. (.not. PML_BOUNDARY_CONDITIONS)) then
        allocate(b_absorb_poro_s_left(NDIM,NGLLZ,nspec_left,NSTEP))
        allocate(b_absorb_poro_s_right(NDIM,NGLLZ,nspec_right,NSTEP))
        allocate(b_absorb_poro_s_bottom(NDIM,NGLLX,nspec_bottom,NSTEP))
        allocate(b_absorb_poro_s_top(NDIM,NGLLX,nspec_top,NSTEP))
        allocate(b_absorb_poro_w_left(NDIM,NGLLZ,nspec_left,NSTEP))
        allocate(b_absorb_poro_w_right(NDIM,NGLLZ,nspec_right,NSTEP))
        allocate(b_absorb_poro_w_bottom(NDIM,NGLLX,nspec_bottom,NSTEP))
        allocate(b_absorb_poro_w_top(NDIM,NGLLX,nspec_top,NSTEP))
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
      if(any_acoustic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3) .and. (.not. PML_BOUNDARY_CONDITIONS)) then
        allocate(b_absorb_acoustic_left(NGLLZ,nspec_left,NSTEP))
        allocate(b_absorb_acoustic_right(NGLLZ,nspec_right,NSTEP))
        allocate(b_absorb_acoustic_bottom(NGLLX,nspec_bottom,NSTEP))
        allocate(b_absorb_acoustic_top(NGLLX,nspec_top,NSTEP))
      else
        allocate(b_absorb_acoustic_left(1,1,1))
        allocate(b_absorb_acoustic_right(1,1,1))
        allocate(b_absorb_acoustic_bottom(1,1,1))
        allocate(b_absorb_acoustic_top(1,1,1))
      endif

  else

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
  allocate(acoustic_edges(4,nelem_acoustic_surface))
  allocate(acoustic_surface(5,nelem_acoustic_surface))
  call read_databases_free_surf(nelem_acoustic_surface,acoustic_edges,any_acoustic_edges)
  ! resets nelem_acoustic_surface
  if( any_acoustic_edges .eqv. .false. ) nelem_acoustic_surface = 0

  ! constructs acoustic surface
  if(nelem_acoustic_surface > 0) then
    call construct_acoustic_surface ( nspec, ngnod, knods, nelem_acoustic_surface, &
                                     acoustic_edges, acoustic_surface)
    if (myrank == 0) then
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
  allocate(fluid_solid_acoustic_ispec(num_fluid_solid_edges))
  allocate(fluid_solid_acoustic_iedge(num_fluid_solid_edges))
  allocate(fluid_solid_elastic_ispec(num_fluid_solid_edges))
  allocate(fluid_solid_elastic_iedge(num_fluid_solid_edges))
  if( num_fluid_poro_edges > 0 ) then
    any_fluid_poro_edges = .true.
  else
    any_fluid_poro_edges = .false.
    num_fluid_poro_edges = 1
  endif
  allocate(fluid_poro_acoustic_ispec(num_fluid_poro_edges))
  allocate(fluid_poro_acoustic_iedge(num_fluid_poro_edges))
  allocate(fluid_poro_poroelastic_ispec(num_fluid_poro_edges))
  allocate(fluid_poro_poroelastic_iedge(num_fluid_poro_edges))
  if ( num_solid_poro_edges > 0 ) then
    any_solid_poro_edges = .true.
  else
    any_solid_poro_edges = .false.
    num_solid_poro_edges = 1
  endif
  allocate(solid_poro_elastic_ispec(num_solid_poro_edges))
  allocate(solid_poro_elastic_iedge(num_solid_poro_edges))
  allocate(solid_poro_poroelastic_ispec(num_solid_poro_edges))
  allocate(solid_poro_poroelastic_iedge(num_solid_poro_edges))

  call read_databases_coupled(num_fluid_solid_edges,any_fluid_solid_edges, &
                            fluid_solid_acoustic_ispec,fluid_solid_elastic_ispec, &
                            num_fluid_poro_edges,any_fluid_poro_edges, &
                            fluid_poro_acoustic_ispec,fluid_poro_poroelastic_ispec, &
                            num_solid_poro_edges,any_solid_poro_edges, &
                            solid_poro_elastic_ispec,solid_poro_poroelastic_ispec)

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
  allocate(nodes_tangential_curve(2,nnodes_tangential_curve))
  allocate(dist_tangential_detection_curve(nnodes_tangential_curve))
  call read_tangential_detection_curve(nnodes_tangential_curve,nodes_tangential_curve, &
                                force_normal_to_surface,rec_normal_to_surface, &
                                any_tangential_curve)
  ! resets nnode_tangential_curve
  if( any_tangential_curve .eqv. .false. ) nnodes_tangential_curve = 0

!                                                                                                                    !axisym
!----  read axial elements data                                                                                      !axisym
!                                                                                                                    !axisym
  allocate(is_on_the_axis(nspec),stat=ier)                                                                           !axisym
  if(ier /= 0) stop 'error: not enough memory to allocate array is_on_the_axis'                                      !axisym
  is_on_the_axis(:) = .false.                                                                                        !axisym
  if(nelem_on_the_axis == 0) then                                                                                    !axisym
    allocate(ispec_of_axial_elements(1))                                                                             !axisym
  else                                                                                                               !axisym
    allocate(ispec_of_axial_elements(nelem_on_the_axis))                                                             !axisym
    call read_databases_axial_elements(nelem_on_the_axis,ispec_of_axial_elements)                                    !axisym
    call build_is_on_the_axis(nspec, nelem_on_the_axis, ispec_of_axial_elements, is_on_the_axis)                     !axisym
  endif                                                                                                              !axisym
  if (myrank == 0 .and. AXISYM) then                                                                                 !axisym
    write(IOUT,*)                                                                                                    !axisym
    write(IOUT,*) 'Number of elements on the axis: ',nelem_on_the_axis                                               !axisym
  endif                                                                                                              !axisym

!
!----  end of reading
!

! closes input Database file
 close(IIN)

!
!---- compute shape functions and their derivatives for SEM grid
!

! set up Gauss-Lobatto-Legendre derivation matrices
  call define_derivation_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,hprimewgll_xx,hprimewgll_zz)

  if (AXISYM) then                                                                                                   !axisym
    ! set up Gauss-Lobatto-Jacobi derivation matrices                                                                !axisym
    call define_GLJ_derivation_matrix(xiglj,wxglj,hprimeBar_xx,hprimeBarwglj_xx)                                     !axisym
  endif                                                                                                              !axisym

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
    call createnum_fast(knods,ibool,shape2D,coorg,nglob,npgeo,nspec,ngnod,myrank)
  else
    call createnum_slow(knods,ibool,nglob,nspec,ngnod,myrank)
  endif

#ifdef USE_MPI
  call MPI_REDUCE(count_nspec_acoustic, count_nspec_acoustic_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
  call MPI_REDUCE(nspec, nspec_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
  call MPI_REDUCE(nglob, nglob_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
#else
  count_nspec_acoustic_total = count_nspec_acoustic
  nspec_total = nspec
  nglob_total = nglob
#endif
  if (myrank == 0) then
    write(IOUT,*)
    write(IOUT,*) 'Total number of elements: ',nspec_total
    write(IOUT,*) 'decomposed as follows:'
    write(IOUT,*)
    write(IOUT,*) 'Total number of elastic/visco/poro elements: ',nspec_total - count_nspec_acoustic_total
    write(IOUT,*) 'Total number of acoustic elements: ',count_nspec_acoustic_total
    write(IOUT,*)
#ifdef USE_MPI
    write(IOUT,*) 'Approximate total number of grid points in the mesh'
    write(IOUT,*) '(with a few duplicates coming from MPI buffers): ',nglob_total
#else
    write(IOUT,*) 'Exact total number of grid points in the mesh: ',nglob_total
#endif

! percentage of elements with 2 degrees of freedom per point
    ratio_2DOFs = (nspec_total - count_nspec_acoustic_total) / dble(nspec_total)
    ratio_1DOF  = count_nspec_acoustic_total / dble(nspec_total)
    nb_acoustic_DOFs = nint(nglob_total*ratio_1DOF)
! elastic elements have two degrees of freedom per point
    nb_elastic_DOFs  = nint(nglob_total*ratio_2DOFs*2)

    if(p_sv) then
      write(IOUT,*)
      write(IOUT,*) 'Approximate number of acoustic degrees of freedom in the mesh: ',nb_acoustic_DOFs
      write(IOUT,*) 'Approximate number of elastic degrees of freedom in the mesh: ',nb_elastic_DOFs
      write(IOUT,*) '  (there are 2 degrees of freedom per point for elastic elements)'
      write(IOUT,*)
      write(IOUT,*) 'Approximate total number of degrees of freedom in the mesh'
      write(IOUT,*) '(sum of the two values above): ',nb_acoustic_DOFs + nb_elastic_DOFs
      write(IOUT,*)
      write(IOUT,*) ' (for simplicity viscoelastic or poroelastic elements, if any,'
      write(IOUT,*) '  are counted as elastic in the above three estimates;'
      write(IOUT,*) '  in reality they have more degrees of freedom)'
      write(IOUT,*)
    endif
  endif

    ! allocate temporary arrays
    allocate(integer_mask_ibool(nglob),stat=ier)
    if( ier /= 0 ) stop 'error allocating mask_ibool'
    allocate(copy_ibool_ori(NGLLX,NGLLZ,nspec),stat=ier)
    if( ier /= 0 ) stop 'error allocating copy_ibool_ori'

    ! reduce cache misses by sorting the global numbering in the order in which it is accessed in the time loop.
    ! this speeds up the calculations significantly on modern processors
    call get_global(nspec,nglob,ibool,copy_ibool_ori,integer_mask_ibool)

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
      if(AXISYM) flagrange_GLJ(j,i) = hgll(j-1,xirec,xiglj,NGLJ)                                               !axisym
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

  if (myrank == 0) then
    write(IOUT,*)
    write(IOUT,*) 'Total number of receivers = ',nrec
    write(IOUT,*)
  endif

  if(nrec < 1) call exit_MPI('need at least one receiver')

! receiver information
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
    allocate(ix_image_color_receiver(nrec))
    allocate(iy_image_color_receiver(nrec))

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
    allocate(coord(NDIM,nglob))

! to display the whole vector field (it needs to be computed from the potential in acoustic elements,
! thus it does not exist as a whole it case of simulations that contain some acoustic elements
! and it thus needs to be computed specifically for display purposes)
    allocate(vector_field_display(3,nglob))

! when periodic boundary conditions are on, some global degrees of freedom are going to be removed,
! thus we need to set this array to zero otherwise some of its locations may contain random values
! if the memory is not cleaned
    vector_field_display(:,:) = 0.d0

    if(assign_external_model) then
      allocate(vpext(NGLLX,NGLLZ,nspec))
      allocate(vsext(NGLLX,NGLLZ,nspec))
      allocate(rhoext(NGLLX,NGLLZ,nspec))
      allocate(QKappa_attenuationext(NGLLX,NGLLZ,nspec))
      allocate(Qmu_attenuationext(NGLLX,NGLLZ,nspec))
      allocate(c11ext(NGLLX,NGLLZ,nspec))
      allocate(c13ext(NGLLX,NGLLZ,nspec))
      allocate(c15ext(NGLLX,NGLLZ,nspec))
      allocate(c33ext(NGLLX,NGLLZ,nspec))
      allocate(c35ext(NGLLX,NGLLZ,nspec))
      allocate(c55ext(NGLLX,NGLLZ,nspec))
      allocate(c12ext(NGLLX,NGLLZ,nspec))
      allocate(c23ext(NGLLX,NGLLZ,nspec))
      allocate(c25ext(NGLLX,NGLLZ,nspec))
    else
      allocate(vpext(1,1,1))
      allocate(vsext(1,1,1))
      allocate(rhoext(1,1,1))
      allocate(QKappa_attenuationext(1,1,1))
      allocate(Qmu_attenuationext(1,1,1))
      allocate(c11ext(1,1,1))
      allocate(c13ext(1,1,1))
      allocate(c15ext(1,1,1))
      allocate(c33ext(1,1,1))
      allocate(c35ext(1,1,1))
      allocate(c55ext(1,1,1))
      allocate(c12ext(1,1,1))
      allocate(c23ext(1,1,1))
      allocate(c25ext(1,1,1))
    endif

!
!----  set the coordinates of the points of the global grid
!
    found_a_negative_jacobian = .false.
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if(AXISYM) then                                                            !axisym
            if (is_on_the_axis(ispec)) then                                          !axisym
              xi = xiglj(i)                                                          !axisym
            else                                                                     !axisym
              xi = xigll(i)                                                          !axisym
            endif                                                                    !axisym
          else                                                                       !axisym
            xi = xigll(i)
          endif                                                                      !axisym
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
    call save_openDX_jacobian(nspec,npgeo,ngnod,knods,coorg,xigll,zigll, &
                                  AXISYM,is_on_the_axis,xiglj)                                                      !axisym)
  endif

! stop the code at the first negative element found, because such a mesh cannot be computed
  if(found_a_negative_jacobian) then

    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if(AXISYM) then                                                            !axisym
            if (is_on_the_axis(ispec)) then                                          !axisym
              xi = xiglj(i)                                                          !axisym
            else                                                                     !axisym
              xi = xigll(i)                                                          !axisym
            endif                                                                    !axisym
          else                                                                       !axisym
            xi = xigll(i)
          endif                                                                      !axisym
          gamma = zigll(j)

          call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl, &
                          jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo, &
                          .true.)

        enddo
      enddo
    enddo

  endif

  xmin_local = minval(coord(1,:))
  xmax_local = maxval(coord(1,:))
  zmin_local = minval(coord(2,:))
  zmax_local = maxval(coord(2,:))

#ifdef USE_MPI
  call MPI_REDUCE(xmin_local, xmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ier)
  call MPI_REDUCE(xmax_local, xmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ier)
  call MPI_REDUCE(zmin_local, zmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ier)
  call MPI_REDUCE(zmax_local, zmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ier)
#else
  xmin = xmin_local
  xmax = xmax_local
  zmin = zmin_local
  zmax = zmax_local
#endif

  if(myrank == 0) then
    write(IOUT,*)
    write(IOUT,*) 'Xmin,Xmax of the whole mesh = ',xmin,xmax
    write(IOUT,*) 'Zmin,Zmax of the whole mesh = ',zmin,zmax
    write(IOUT,*)

! check that no source is located outside the mesh
    do i = 1,NSOURCES
      if(x_source(i) < xmin) stop 'error: at least one source has x < xmin of the mesh'
      if(x_source(i) > xmax) stop 'error: at least one source has x > xmax of the mesh'

      if(z_source(i) < zmin) stop 'error: at least one source has z < zmin of the mesh'
      if(z_source(i) > zmax) stop 'error: at least one source has z > zmax of the mesh'
    enddo

  endif

! use a spring to improve the stability of the Stacey condition
  x_center_spring = (xmax + xmin)/2.d0
  z_center_spring = (zmax + zmin)/2.d0

! allocate an array to make sure that an acoustic free surface is not enforced on periodic edges
  allocate(this_ibool_is_a_periodic_edge(NGLOB))
  this_ibool_is_a_periodic_edge(:) = .false.

! periodic conditions: detect common points between left and right edges and replace one of them with the other
    if(ADD_PERIODIC_CONDITIONS) then

      if (myrank == 0) then
        write(IOUT,*)
        write(IOUT,*) 'implementing periodic boundary conditions'
        write(IOUT,*) 'in the horizontal direction with a periodicity distance of ',PERIODIC_HORIZ_DIST,' m'
        if(PERIODIC_HORIZ_DIST <= 0.d0) stop 'PERIODIC_HORIZ_DIST should be greater than zero when using ADD_PERIODIC_CONDITIONS'
        write(IOUT,*)
        write(IOUT,*) '*****************************************************************'
        write(IOUT,*) '*****************************************************************'
        write(IOUT,*) '**** BEWARE: because of periodic conditions, values computed ****'
        write(IOUT,*) '****         by checkgrid() below will not be reliable       ****'
        write(IOUT,*) '*****************************************************************'
        write(IOUT,*) '*****************************************************************'
        write(IOUT,*)
      endif

! set up a local geometric tolerance

  xtypdist = +HUGEVAL

  do ispec = 1,nspec

  xminval = +HUGEVAL
  yminval = +HUGEVAL
  xmaxval = -HUGEVAL
  ymaxval = -HUGEVAL

! only loop on the four corners of each element to get a typical size
  do j = 1,NGLLZ,NGLLZ-1
    do i = 1,NGLLX,NGLLX-1
      iglob = ibool(i,j,ispec)
      xmaxval = max(coord(1,iglob),xmaxval)
      xminval = min(coord(1,iglob),xminval)
      ymaxval = max(coord(2,iglob),ymaxval)
      yminval = min(coord(2,iglob),yminval)
    enddo
  enddo

! compute the minimum typical "size" of an element in the mesh
  xtypdist = min(xtypdist,xmaxval-xminval)
  xtypdist = min(xtypdist,ymaxval-yminval)

  enddo

! define a tolerance, small with respect to the minimum size
  xtol = 1.d-4 * xtypdist

! detect the points that are on the same horizontal line (i.e. at the same height Z)
! and that have a value of the horizontal coordinate X that differs by exactly the periodicity length;
! if so, make them all have the same global number, which will then implement periodic boundary conditions automatically.
! We select the smallest value of iglob and assign it to all the points that are the same due to periodicity,
! this way the maximum value of the ibool() array will remain as small as possible.
!
! *** IMPORTANT: this simple algorithm will be slow for large meshes because it has a cost of NGLOB^2 / 2
! (where NGLOB is the number of points per MPI slice, not of the whole mesh though). This could be
! reduced to O(NGLOB log(NGLOB)) by using a quicksort algorithm on the coordinates of the points to detect the multiples
! (as implemented in routine createnum_fast() elsewhere in the code). This could be done one day if needed instead
! of the very simple double loop below.
      if (myrank == 0) write(IOUT,*) &
        'start detecting points for periodic boundary conditions (the current algorithm can be slow and could be improved)...'
      counter = 0
      do iglob = 1,NGLOB-1
        do iglob2 = iglob + 1,NGLOB
          ! check if the two points have the exact same Z coordinate
          if(abs(coord(2,iglob2) - coord(2,iglob)) < xtol) then
            ! if so, check if their X coordinate differs by exactly the periodicity distance
            if(abs(abs(coord(1,iglob2) - coord(1,iglob)) - PERIODIC_HORIZ_DIST) < xtol) then
              ! if so, they are the same point, thus replace the highest value of ibool with the lowest
              ! to make them the same global point and thus implement periodicity automatically
              counter = counter + 1
              this_ibool_is_a_periodic_edge(iglob) = .true.
              this_ibool_is_a_periodic_edge(iglob2) = .true.
              do ispec = 1,nspec
                do j = 1,NGLLZ
                  do i = 1,NGLLX
                    if(ibool(i,j,ispec) == iglob2) ibool(i,j,ispec) = iglob
                  enddo
                enddo
              enddo
            endif
          endif
        enddo
      enddo

#ifdef USE_MPI
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
#endif

      if (myrank == 0) write(IOUT,*) 'done detecting points for periodic boundary conditions.'

      if(counter > 0) write(IOUT,*) 'implemented periodic conditions on ',counter,' grid points on proc ',myrank

    endif ! of if(ADD_PERIODIC_CONDITIONS)

!
!--- save the grid of points in a file
!
  if(output_grid_ASCII .and. myrank == 0) then
     write(IOUT,*)
     write(IOUT,*) 'Saving the grid in an ASCII text file...'
     write(IOUT,*)
     open(unit=55,file='OUTPUT_FILES/ASCII_dump_of_grid_points.txt',status='unknown')
     write(55,*) nglob
     do n = 1,nglob
        write(55,*) (coord(i,n), i=1,NDIM)
     enddo
     close(55)
  endif

!
!-----   plot the GLL mesh in a Gnuplot file
!
  if(output_grid_Gnuplot .and. myrank == 0)  &
    call plotgll(knods,ibool,coorg,coord,nglob,npgeo,ngnod,nspec)

  if (assign_external_model) then
    if(myrank == 0) write(IOUT,*) 'Assigning an external velocity and density model...'
    call read_external_model(any_acoustic,any_elastic,any_poroelastic, &
                elastic,poroelastic,anisotropic,nspec,nglob,N_SLS,ibool, &
                f0_attenuation,inv_tau_sigma_nu1_sent,phi_nu1_sent, &
                inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu1_sent,Mu_nu2_sent, &
                inv_tau_sigma_nu1,inv_tau_sigma_nu2,phi_nu1,phi_nu2,Mu_nu1,Mu_nu2,&
                coord,kmato,rhoext,vpext,vsext, &
                QKappa_attenuationext,Qmu_attenuationext, &
                c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,READ_EXTERNAL_SEP_FILE)
  endif

!
!----  perform basic checks on parameters read
!
  all_anisotropic = .false.
  if(count(anisotropic(:) .eqv. .true.) == nspec) all_anisotropic = .true.

  if(all_anisotropic .and. anyabs) &
    call exit_MPI('Cannot put absorbing boundaries if anisotropic materials along edges')

  if(ATTENUATION_VISCOELASTIC_SOLID .and. all_anisotropic) then
    call exit_MPI('Cannot turn attenuation on in anisotropic materials')
  endif

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
  if(ATTENUATION_VISCOELASTIC_SOLID .and. .not. any_elastic_glob) &
    call exit_MPI('currently cannot have attenuation if acoustic/poroelastic simulation only')

!
!----   define coefficients of the Newmark time scheme
!
  deltatover2 = HALF*deltat
  deltatsquareover2 = HALF*deltat*deltat

  if(SIMULATION_TYPE == 3) then
!  define coefficients of the Newmark time scheme for the backward wavefield
    b_deltat = - deltat
    b_deltatover2 = HALF*b_deltat
    b_deltatsquareover2 = HALF*b_deltat*b_deltat
  endif

!---- define actual location of source and receivers

  call setup_sources_receivers(NSOURCES,initialfield,source_type, &
     coord,ibool,nglob,nspec,nelem_acoustic_surface,acoustic_surface,elastic,poroelastic, &
     x_source,z_source,ispec_selected_source,ispec_selected_rec, &
     is_proc_source,nb_proc_source, &
     sourcearray,Mxx,Mzz,Mxz,xix,xiz,gammax,gammaz,xigll,zigll,npgeo, &
     nproc,myrank,xi_source,gamma_source,coorg,knods,ngnod, &
     nrec,nrecloc,recloc,which_proc_receiver,st_xval,st_zval, &
     xi_receiver,gamma_receiver,station_name,network_name,x_final_receiver,z_final_receiver,iglob_source)

! compute source array for adjoint source
  nadj_rec_local = 0
  if(SIMULATION_TYPE == 3) then  ! adjoint calculation

    do irec = 1,nrec
      if(myrank == which_proc_receiver(irec)) then
        ! check that the source proc number is okay
        if(which_proc_receiver(irec) < 0 .or. which_proc_receiver(irec) > NPROC-1) &
              call exit_MPI('something is wrong with the source proc number in adjoint simulation')
        nadj_rec_local = nadj_rec_local + 1
      endif
    enddo
    allocate(adj_sourcearray(NSTEP,3,NGLLX,NGLLZ))
    if (nadj_rec_local > 0)  then
      allocate(adj_sourcearrays(nadj_rec_local,NSTEP,3,NGLLX,NGLLZ))
    else
      allocate(adj_sourcearrays(1,1,1,1,1))
    endif

    if (.not. SU_FORMAT) then
       irec_local = 0
       do irec = 1, nrec
         ! compute only adjoint source arrays in the local proc
         if(myrank == which_proc_receiver(irec))then
           irec_local = irec_local + 1
           adj_source_file = trim(station_name(irec))//'.'//trim(network_name(irec))
           call compute_arrays_adj_source(adj_source_file, &
                               xi_receiver(irec), gamma_receiver(irec), &
                               adj_sourcearray, xigll,zigll,NSTEP)
           adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
         endif
       enddo
    else
       irec_local = 0
       write(filename, "('./SEM/Ux_file_single.bin.adj')")
       open(111,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ios)
               if (ios /= 0) call exit_MPI(' file '//trim(filename)//'does not exist')
       write(filename, "('./SEM/Uz_file_single.bin.adj')")
       open(113,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ios)
               if (ios /= 0) call exit_MPI(' file '//trim(filename)//'does not exist')

       allocate(adj_src_s(NSTEP,3))

       do irec = 1, nrec
         if(myrank == which_proc_receiver(irec))then
          irec_local = irec_local + 1
          adj_sourcearray(:,:,:,:) = 0.0
          read(111,rec=irec,iostat=ios) r4head, adj_src_s(:,1)
               if (ios /= 0) call exit_MPI(' file '//trim(filename)//' read error')
          read(113,rec=irec,iostat=ios) r4head, adj_src_s(:,3)
               if (ios /= 0) call exit_MPI(' file '//trim(filename)//' read error')
          header2=int(r4head(29), kind=2)
          if (irec==1) print*, r4head(1),r4head(19),r4head(20),r4head(21),r4head(22),header2(2)
          call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
          call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
          do k = 1, NGLLZ
              do i = 1, NGLLX
                adj_sourcearray(:,:,i,k) = hxir(i) * hgammar(k) * adj_src_s(:,:)
              enddo
          enddo
          adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
         endif !  if(myrank == which_proc_receiver(irec))
       enddo ! irec
       close(111)
       close(113)
       deallocate(adj_src_s)
    endif
  else
     allocate(adj_sourcearrays(1,1,1,1,1))
  endif

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

!
!--- tangential computation
!

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
          call compute_normal_vector( anglesource(i_source), &
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
            anglesource_recv = anglesource(i_source)
#ifdef USE_MPI
          else if ( myrank == 0 ) then
            do i = 1, nb_proc_source(i_source) - is_proc_source(i_source)
              call MPI_recv(source_courbe_eros(i_source),1,MPI_INTEGER, &
                          MPI_ANY_SOURCE,42,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
              call MPI_recv(anglesource_recv,1,MPI_DOUBLE_PRECISION, &
                          MPI_ANY_SOURCE,43,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
            enddo
          else if ( is_proc_source(i_source) == 1 ) then
            call MPI_send(n1_tangential_detection_curve,1,MPI_INTEGER,0,42,MPI_COMM_WORLD,ier)
            call MPI_send(anglesource(i_source),1,MPI_DOUBLE_PRECISION,0,43,MPI_COMM_WORLD,ier)
#endif
          endif

#ifdef USE_MPI
          call MPI_bcast(anglesource_recv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
          anglesource(i_source) = anglesource_recv
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
               which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
            call MPI_RECV(x_final_receiver_dummy,1,MPI_DOUBLE_PRECISION,&
               which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
            call MPI_RECV(z_final_receiver_dummy,1,MPI_DOUBLE_PRECISION,&
               which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

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

!
!---
!

! allocate seismogram arrays
  allocate(sisux(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrecloc))
  allocate(sisuz(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrecloc))
  allocate(siscurl(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrecloc))

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
    if (AXISYM) then
      if((is_proc_source(i) == 1) .and. is_on_the_axis(ispec_selected_source(i))) then                    !axisym
        call lagrange_any(xi_source(i),NGLJ,xiglj,hxis,hpxis)                                             !axisym
      else                                                                                                !axisym
        call lagrange_any(xi_source(i),NGLLX,xigll,hxis,hpxis)                                            !axisym
      endif                                                                                               !axisym
    else
      call lagrange_any(xi_source(i),NGLLX,xigll,hxis,hpxis)
    endif                                                                                                 !axisym

    call lagrange_any(gamma_source(i),NGLLZ,zigll,hgammas,hpgammas)
    hxis_store(i,:) = hxis(:)
    hgammas_store(i,:) = hgammas(:)
  enddo

! displacement, velocity, acceleration and inverse of the mass matrix for elastic elements
    if(any_elastic) then
      nglob_elastic = nglob
    else
      ! allocate unused arrays with fictitious size
      nglob_elastic = 1
    endif
    allocate(displ_elastic(3,nglob_elastic))
    allocate(displ_elastic_old(3,nglob_elastic))
    allocate(veloc_elastic(3,nglob_elastic))
    allocate(accel_elastic(3,nglob_elastic))
    allocate(accel_elastic_adj_coupling(3,nglob_elastic))

    allocate(rmass_inverse_elastic_one(nglob_elastic))
    allocate(rmass_inverse_elastic_three(nglob_elastic))

    if(time_stepping_scheme==2)then
      allocate(displ_elastic_LDDRK(3,nglob_elastic))
      allocate(veloc_elastic_LDDRK(3,nglob_elastic))
      allocate(veloc_elastic_LDDRK_temp(3,nglob_elastic))
    endif

    if(time_stepping_scheme == 3)then
      allocate(accel_elastic_rk(3,nglob_elastic,stage_time_scheme))
      allocate(veloc_elastic_rk(3,nglob_elastic,stage_time_scheme))
      allocate(veloc_elastic_initial_rk(3,nglob_elastic))
      allocate(displ_elastic_initial_rk(3,nglob_elastic))
    endif

    ! extra array if adjoint and kernels calculation
    if(SIMULATION_TYPE == 3 .and. any_elastic) then
      allocate(b_displ_elastic(3,nglob))
      allocate(b_displ_elastic_old(3,nglob))
      allocate(b_veloc_elastic(3,nglob))
      allocate(b_accel_elastic(3,nglob))
      allocate(rho_kl(NGLLX,NGLLZ,nspec))
      allocate(rho_k(nglob))
      allocate(rhol_global(nglob))
      allocate(mu_kl(NGLLX,NGLLZ,nspec))
      allocate(mu_k(nglob))
      allocate(mul_global(nglob))
      allocate(kappa_kl(NGLLX,NGLLZ,nspec))
      allocate(kappa_k(nglob))
      allocate(kappal_global(nglob))
      allocate(rhop_kl(NGLLX,NGLLZ,nspec))
      allocate(alpha_kl(NGLLX,NGLLZ,nspec))
      allocate(beta_kl(NGLLX,NGLLZ,nspec))
      allocate(rhorho_el_hessian_final2(NGLLX,NGLLZ,nspec))
      allocate(rhorho_el_hessian_temp2(nglob))
      allocate(rhorho_el_hessian_final1(NGLLX,NGLLZ,nspec))
      allocate(rhorho_el_hessian_temp1(nglob))
    else
      allocate(b_displ_elastic(1,1))
      allocate(b_displ_elastic_old(1,1))
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
      nglob_poroelastic = nglob
    else
      ! allocate unused arrays with fictitious size
      nglob_poroelastic = 1
    endif
    allocate(displs_poroelastic(NDIM,nglob_poroelastic))
    allocate(displs_poroelastic_old(NDIM,nglob_poroelastic))
    allocate(velocs_poroelastic(NDIM,nglob_poroelastic))
    allocate(accels_poroelastic(NDIM,nglob_poroelastic))
    allocate(accels_poroelastic_adj_coupling(NDIM,nglob_poroelastic))
    allocate(rmass_s_inverse_poroelastic(nglob_poroelastic))
    allocate(displw_poroelastic(NDIM,nglob_poroelastic))
    allocate(velocw_poroelastic(NDIM,nglob_poroelastic))
    allocate(accelw_poroelastic(NDIM,nglob_poroelastic))
    allocate(accelw_poroelastic_adj_coupling(NDIM,nglob_poroelastic))
    allocate(rmass_w_inverse_poroelastic(nglob_poroelastic))

    if(time_stepping_scheme == 2)then
    allocate(displs_poroelastic_LDDRK(NDIM,nglob_poroelastic))
    allocate(velocs_poroelastic_LDDRK(NDIM,nglob_poroelastic))
    allocate(displw_poroelastic_LDDRK(NDIM,nglob_poroelastic))
    allocate(velocw_poroelastic_LDDRK(NDIM,nglob_poroelastic))
    endif

    if(time_stepping_scheme == 3)then
    allocate(accels_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
    allocate(velocs_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
    allocate(accelw_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
    allocate(velocw_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
    allocate(displs_poroelastic_initial_rk(NDIM,nglob_poroelastic))
    allocate(velocs_poroelastic_initial_rk(NDIM,nglob_poroelastic))
    allocate(displw_poroelastic_initial_rk(NDIM,nglob_poroelastic))
    allocate(velocw_poroelastic_initial_rk(NDIM,nglob_poroelastic))
    endif

    ! extra array if adjoint and kernels calculation
    if(SIMULATION_TYPE == 3 .and. any_poroelastic) then
      allocate(b_displs_poroelastic(NDIM,nglob))
      allocate(b_velocs_poroelastic(NDIM,nglob))
      allocate(b_accels_poroelastic(NDIM,nglob))
      allocate(b_displw_poroelastic(NDIM,nglob))
      allocate(b_velocw_poroelastic(NDIM,nglob))
      allocate(b_accelw_poroelastic(NDIM,nglob))
      allocate(rhot_kl(NGLLX,NGLLZ,nspec))
      allocate(rhot_k(nglob))
      allocate(rhof_kl(NGLLX,NGLLZ,nspec))
      allocate(rhof_k(nglob))
      allocate(sm_kl(NGLLX,NGLLZ,nspec))
      allocate(sm_k(nglob))
      allocate(eta_kl(NGLLX,NGLLZ,nspec))
      allocate(eta_k(nglob))
      allocate(mufr_kl(NGLLX,NGLLZ,nspec))
      allocate(mufr_k(nglob))
      allocate(B_kl(NGLLX,NGLLZ,nspec))
      allocate(B_k(nglob))
      allocate(C_kl(NGLLX,NGLLZ,nspec))
      allocate(C_k(nglob))
      allocate(M_kl(NGLLX,NGLLZ,nspec))
      allocate(M_k(nglob))
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
      allocate(phil_global(nglob))
      allocate(mulfr_global(nglob))
      allocate(etal_f_global(nglob))
      allocate(rhol_s_global(nglob))
      allocate(rhol_f_global(nglob))
      allocate(rhol_bar_global(nglob))
      allocate(tortl_global(nglob))
      allocate(permlxx_global(nglob))
      allocate(permlxz_global(nglob))
      allocate(permlzz_global(nglob))
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
      allocate(icount(nglob))
    else
      allocate(icount(1))
    endif

    ! potential, its first and second derivative, and inverse of the mass matrix for acoustic elements
    if(any_acoustic) then
      nglob_acoustic = nglob
    else
      ! allocate unused arrays with fictitious size
      nglob_acoustic = 1
    endif
    allocate(potential_acoustic(nglob_acoustic))
    allocate(potential_acoustic_old(nglob_acoustic))
    allocate(potential_acoustic_adj_coupling(nglob_acoustic))
    allocate(potential_dot_acoustic(nglob_acoustic))
    allocate(potential_dot_dot_acoustic(nglob_acoustic))
    allocate(rmass_inverse_acoustic(nglob_acoustic))
    if(time_stepping_scheme == 2) then
    allocate(potential_acoustic_LDDRK(nglob_acoustic))
    allocate(potential_dot_acoustic_LDDRK(nglob_acoustic))
    allocate(potential_dot_acoustic_temp(nglob_acoustic))
    endif

    if(time_stepping_scheme == 3) then
    allocate(potential_acoustic_init_rk(nglob_acoustic))
    allocate(potential_dot_acoustic_init_rk(nglob_acoustic))
    allocate(potential_dot_dot_acoustic_rk(nglob_acoustic,stage_time_scheme))
    allocate(potential_dot_acoustic_rk(nglob_acoustic,stage_time_scheme))
    endif


    if(SIMULATION_TYPE == 3 .and. any_acoustic) then
      allocate(b_potential_acoustic(nglob))
      allocate(b_potential_acoustic_old(nglob))
      allocate(b_potential_dot_acoustic(nglob))
      allocate(b_potential_dot_dot_acoustic(nglob))
      allocate(b_displ_ac(2,nglob))
      allocate(b_accel_ac(2,nglob))
      allocate(accel_ac(2,nglob))
      allocate(rho_ac_kl(NGLLX,NGLLZ,nspec))
      allocate(rhol_ac_global(nglob))
      allocate(kappa_ac_kl(NGLLX,NGLLZ,nspec))
      allocate(kappal_ac_global(nglob))
      allocate(rhop_ac_kl(NGLLX,NGLLZ,nspec))
      allocate(alpha_ac_kl(NGLLX,NGLLZ,nspec))
      allocate(rhorho_ac_hessian_final2(NGLLX,NGLLZ,nspec))
      allocate(rhorho_ac_hessian_final1(NGLLX,NGLLZ,nspec))
    else
    ! allocate unused arrays with fictitious size
      allocate(b_potential_acoustic(1))
      allocate(b_potential_acoustic_old(1))
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

  ! PML absorbing conditionds
    anyabs_glob=anyabs
#ifdef USE_MPI
    call MPI_ALLREDUCE(anyabs, anyabs_glob, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ier)
#endif

    ! allocate this in all cases, even if PML is not set, because we use it for PostScript display as well
    allocate(is_PML(nspec),stat=ier)
    if(ier /= 0) stop 'error: not enough memory to allocate array is_PML'
    is_PML(:) = .false.

    if(PML_BOUNDARY_CONDITIONS .and. anyabs_glob ) then
      allocate(spec_to_PML(nspec),stat=ier)
      if(ier /= 0) stop 'error: not enough memory to allocate array spec_to_PML'

      allocate(which_PML_elem(4,nspec),stat=ier)
      if(ier /= 0) stop 'error: not enough memory to allocate array which_PML_elem'
      which_PML_elem(:,:) = .false.

      if(SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))then
        allocate(PML_interior_interface(4,nspec),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array PML_interior_interface'
        PML_interior_interface = .false.
      else
        allocate(PML_interior_interface(4,1))
      endif

!   DK DK add support for using pml in mpi mode with external mesh
      if(read_external_mesh)then
        allocate(mask_ibool_pml(nglob))
      else
        allocate(mask_ibool_pml(1))
      endif

      call pml_init(myrank,SIMULATION_TYPE,SAVE_FORWARD,nspec,nglob,ibool,anyabs,nelemabs,codeabs,numabs,&
                    NELEM_PML_THICKNESS,nspec_PML,is_PML,which_PML_elem,spec_to_PML,region_CPML,&
                    PML_interior_interface,nglob_interface,mask_ibool_pml,read_external_mesh, &
                    AXISYM)                                                                                          !axisym

      if((SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD)) .and. PML_BOUNDARY_CONDITIONS)then

        if(nglob_interface > 0) then
          allocate(point_interface(nglob_interface),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array point_interface'
        endif

        if(any_elastic .and. nglob_interface > 0)then
          allocate(pml_interface_history_displ(3,nglob_interface,NSTEP),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_displ'
          allocate(pml_interface_history_veloc(3,nglob_interface,NSTEP),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_veloc'
          allocate(pml_interface_history_accel(3,nglob_interface,NSTEP),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_accel'
        endif

        if(any_acoustic .and. nglob_interface > 0)then
          allocate(pml_interface_history_potential(nglob_interface,NSTEP),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_potential'
          allocate(pml_interface_history_potential_dot(nglob_interface,NSTEP),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_potential_dot'
          allocate(pml_interface_history_potential_dot_dot(nglob_interface,NSTEP),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_potential_dot_dot'
        endif

        if(nglob_interface > 0) then
          call determin_interface_pml_interior(nglob_interface,nspec,ibool,PML_interior_interface,&
                                               which_PML_elem,point_interface,read_external_mesh,&
                                               mask_ibool_pml,region_CPML,nglob)
          deallocate(PML_interior_interface)
          deallocate(mask_ibool_pml)
        endif

        if(any_elastic .and. nglob_interface > 0)then
          write(outputname,'(a,i6.6,a)') 'pml_interface_elastic',myrank,'.bin'
          open(unit=71,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
        endif

        if(any_acoustic .and. nglob_interface > 0)then
          write(outputname,'(a,i6.6,a)') 'pml_interface_acoustic',myrank,'.bin'
          open(unit=72,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
        endif
      else
        allocate(point_interface(1))
        allocate(pml_interface_history_displ(3,1,1))
        allocate(pml_interface_history_veloc(3,1,1))
        allocate(pml_interface_history_accel(3,1,1))
        allocate(pml_interface_history_potential(1,1))
        allocate(pml_interface_history_potential_dot(1,1))
        allocate(pml_interface_history_potential_dot_dot(1,1))
      endif

      if(SIMULATION_TYPE == 3 .and. PML_BOUNDARY_CONDITIONS)then

        if(any_elastic .and. nglob_interface > 0)then
          do it = 1,NSTEP
            do i = 1, nglob_interface
              read(71)pml_interface_history_accel(1,i,it),pml_interface_history_accel(2,i,it),&
                      pml_interface_history_accel(3,i,it),&
                      pml_interface_history_veloc(1,i,it),pml_interface_history_veloc(2,i,it),&
                      pml_interface_history_veloc(3,i,it),&
                      pml_interface_history_displ(1,i,it),pml_interface_history_displ(2,i,it),&
                      pml_interface_history_displ(3,i,it)
          enddo
         enddo
       endif

       if(any_acoustic .and. nglob_interface > 0)then
         do it = 1,NSTEP
           do i = 1, nglob_interface
             read(72)pml_interface_history_potential_dot_dot(i,it),pml_interface_history_potential_dot(i,it),&
                     pml_interface_history_potential(i,it)
           enddo
         enddo
       endif
     endif

      deallocate(which_PML_elem)

      if (nspec_PML==0) nspec_PML=1 ! DK DK added this

      if (nspec_PML > 0) then
        allocate(K_x_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array K_x_store'
        allocate(K_z_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array K_z_store'
        allocate(d_x_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array d_x_store'
        allocate(d_z_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array d_z_store'
        allocate(alpha_x_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array alpha_x_store'
        allocate(alpha_z_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array alpha_z_store'
        K_x_store(:,:,:) = ZERO
        K_z_store(:,:,:) = ZERO
        d_x_store(:,:,:) = ZERO
        d_z_store(:,:,:) = ZERO
        alpha_x_store(:,:,:) = ZERO
        alpha_z_store(:,:,:) = ZERO
        call define_PML_coefficients(nglob,nspec,kmato,density,poroelastcoef,numat,f0(1),&
                  ibool,coord,is_PML,region_CPML,spec_to_PML,nspec_PML,&
                  K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store)
      else
        allocate(K_x_store(NGLLX,NGLLZ,1),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array K_x_store'
        allocate(K_z_store(NGLLX,NGLLZ,1),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array K_z_store'
        allocate(d_x_store(NGLLX,NGLLZ,1),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array d_x_store'
        allocate(d_z_store(NGLLX,NGLLZ,1),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array d_z_store'
        allocate(alpha_x_store(NGLLX,NGLLZ,1),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array alpha_x_store'
        allocate(alpha_z_store(NGLLX,NGLLZ,1),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array alpha_z_store'
        K_x_store(:,:,:) = ZERO
        K_z_store(:,:,:) = ZERO
        d_x_store(:,:,:) = ZERO
        d_z_store(:,:,:) = ZERO
        alpha_x_store(:,:,:) = ZERO
        alpha_z_store(:,:,:) = ZERO
      endif

      !elastic PML memory variables
      if (any_elastic .and. nspec_PML>0) then
        allocate(rmemory_displ_elastic(2,3,NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_displ_elastic'
        allocate(rmemory_dux_dx(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dx'
        allocate(rmemory_dux_dz(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dz'
        allocate(rmemory_duz_dx(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dx'
        allocate(rmemory_duz_dz(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dz'
        if(any_acoustic .and. num_fluid_solid_edges > 0)then
          allocate(rmemory_fsb_displ_elastic(1,3,NGLLX,NGLLZ,num_fluid_solid_edges),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_fsb_displ_elastic'
        endif

        if(ROTATE_PML_ACTIVATE)then
          allocate(rmemory_dux_dx_prime(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dx_prime'
          allocate(rmemory_dux_dz_prime(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dz_prime'
          allocate(rmemory_duz_dx_prime(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dx_prime'
          allocate(rmemory_duz_dz_prime(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dz_prime'
        else
          allocate(rmemory_dux_dx_prime(1,1,1,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dx_prime'
          allocate(rmemory_dux_dz_prime(1,1,1,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dz_prime'
          allocate(rmemory_duz_dx_prime(1,1,1,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dx_prime'
          allocate(rmemory_duz_dz_prime(1,1,1,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dz_prime'
        endif

        if(time_stepping_scheme == 2)then
          allocate(rmemory_displ_elastic_LDDRK(2,3,NGLLX,NGLLZ,nspec_PML),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_displ_elastic'
          allocate(rmemory_dux_dx_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dx'
          allocate(rmemory_dux_dz_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dz'
          allocate(rmemory_duz_dx_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dx'
          allocate(rmemory_duz_dz_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dz'
        else
          allocate(rmemory_displ_elastic_LDDRK(1,1,1,1,1),stat=ier)
          allocate(rmemory_dux_dx_LDDRK(1,1,1,2),stat=ier)
          allocate(rmemory_dux_dz_LDDRK(1,1,1,2),stat=ier)
          allocate(rmemory_duz_dx_LDDRK(1,1,1,2),stat=ier)
          allocate(rmemory_duz_dz_LDDRK(1,1,1,2),stat=ier)
        endif

        rmemory_displ_elastic(:,:,:,:,:) = ZERO
        rmemory_dux_dx(:,:,:,:) = ZERO
        rmemory_dux_dz(:,:,:,:) = ZERO
        rmemory_duz_dx(:,:,:,:) = ZERO
        rmemory_duz_dz(:,:,:,:) = ZERO

        if(any_acoustic .and. num_fluid_solid_edges > 0)then
          rmemory_fsb_displ_elastic(:,:,:,:,:) = ZERO
        endif

        if(ROTATE_PML_ACTIVATE)then
          rmemory_dux_dx_prime(:,:,:,:) = ZERO
          rmemory_dux_dz_prime(:,:,:,:) = ZERO
          rmemory_duz_dx_prime(:,:,:,:) = ZERO
          rmemory_duz_dz_prime(:,:,:,:) = ZERO
        endif

        if(time_stepping_scheme == 2)then
          rmemory_displ_elastic_LDDRK(:,:,:,:,:) = ZERO
          rmemory_dux_dx_LDDRK(:,:,:,:) = ZERO
          rmemory_dux_dz_LDDRK(:,:,:,:) = ZERO
          rmemory_duz_dx_LDDRK(:,:,:,:) = ZERO
          rmemory_duz_dz_LDDRK(:,:,:,:) = ZERO
        endif

      else

        allocate(rmemory_displ_elastic(1,1,1,1,1))
        allocate(rmemory_dux_dx(1,1,1,1))
        allocate(rmemory_dux_dz(1,1,1,1))
        allocate(rmemory_duz_dx(1,1,1,1))
        allocate(rmemory_duz_dz(1,1,1,1))
        if(any_acoustic .and. num_fluid_solid_edges > 0)then
          allocate(rmemory_fsb_displ_elastic(1,3,NGLLX,NGLLZ,1))
        endif

        allocate(rmemory_dux_dx_prime(1,1,1,1))
        allocate(rmemory_dux_dz_prime(1,1,1,1))
        allocate(rmemory_duz_dx_prime(1,1,1,1))
        allocate(rmemory_duz_dz_prime(1,1,1,1))

        allocate(rmemory_displ_elastic_LDDRK(1,1,1,1,1))
        allocate(rmemory_dux_dx_LDDRK(1,1,1,1))
        allocate(rmemory_dux_dz_LDDRK(1,1,1,1))
        allocate(rmemory_duz_dx_LDDRK(1,1,1,1))
        allocate(rmemory_duz_dz_LDDRK(1,1,1,1))
      endif

      if (any_acoustic .and. nspec_PML>0) then
        allocate(rmemory_potential_acoustic(2,NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_potential_acoustic'
        allocate(rmemory_acoustic_dux_dx(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_acoustic_dux_dx'
        allocate(rmemory_acoustic_dux_dz(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_acoustic_dux_dz'

        rmemory_potential_acoustic = ZERO
        rmemory_acoustic_dux_dx = ZERO
        rmemory_acoustic_dux_dz = ZERO

        if(time_stepping_scheme == 2)then
          allocate(rmemory_potential_acoustic_LDDRK(2,NGLLX,NGLLZ,nspec_PML),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_potential_acoustic'
          allocate(rmemory_acoustic_dux_dx_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_acoustic_dux_dx'
          allocate(rmemory_acoustic_dux_dz_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_acoustic_dux_dz'
        else
          allocate(rmemory_potential_acoustic_LDDRK(1,1,1,1),stat=ier)
          allocate(rmemory_acoustic_dux_dx_LDDRK(1,1,1,1),stat=ier)
          allocate(rmemory_acoustic_dux_dz_LDDRK(1,1,1,1),stat=ier)
        endif

        rmemory_potential_acoustic_LDDRK = ZERO
        rmemory_acoustic_dux_dx_LDDRK = ZERO
        rmemory_acoustic_dux_dz_LDDRK = ZERO

      else
        allocate(rmemory_potential_acoustic(1,1,1,1))
        allocate(rmemory_acoustic_dux_dx(1,1,1,1))
        allocate(rmemory_acoustic_dux_dz(1,1,1,1))
      endif

    else
      allocate(rmemory_dux_dx(1,1,1,1))
      allocate(rmemory_dux_dz(1,1,1,1))
      allocate(rmemory_duz_dx(1,1,1,1))
      allocate(rmemory_duz_dz(1,1,1,1))
      allocate(rmemory_fsb_displ_elastic(1,3,NGLLX,NGLLZ,1))

      allocate(rmemory_dux_dx_prime(1,1,1,1))
      allocate(rmemory_dux_dz_prime(1,1,1,1))
      allocate(rmemory_duz_dx_prime(1,1,1,1))
      allocate(rmemory_duz_dz_prime(1,1,1,1))

      allocate(rmemory_displ_elastic(1,1,1,1,1))

      allocate(rmemory_displ_elastic_LDDRK(1,1,1,1,1))
      allocate(rmemory_dux_dx_LDDRK(1,1,1,1))
      allocate(rmemory_dux_dz_LDDRK(1,1,1,1))
      allocate(rmemory_duz_dx_LDDRK(1,1,1,1))
      allocate(rmemory_duz_dz_LDDRK(1,1,1,1))

      allocate(rmemory_potential_acoustic(1,1,1,1))
      allocate(rmemory_acoustic_dux_dx(1,1,1,1))
      allocate(rmemory_acoustic_dux_dz(1,1,1,1))

      allocate(rmemory_potential_acoustic_LDDRK(1,1,1,1))
      allocate(rmemory_acoustic_dux_dx_LDDRK(1,1,1,1))
      allocate(rmemory_acoustic_dux_dz_LDDRK(1,1,1,1))

      allocate(spec_to_PML(1))

      allocate(K_x_store(1,1,1))
      allocate(K_z_store(1,1,1))
      allocate(d_x_store(1,1,1))
      allocate(d_z_store(1,1,1))
      allocate(alpha_x_store(1,1,1))
      allocate(alpha_z_store(1,1,1))
    endif ! PML_BOUNDARY_CONDITIONS

  ! avoid a potential side effect owing to the "if" statements above: this array may be unallocated,
  ! if so we need to allocate a dummy version in order to be able to use that array as an argument
  ! in some subroutine calls below
  if(.not. allocated(rmemory_fsb_displ_elastic)) allocate(rmemory_fsb_displ_elastic(1,3,NGLLX,NGLLZ,1))

! Test compatibility with axisym.                                                                                    !axisym
  if(AXISYM) then                                                                                                    !axisym
    call check_compatibility_axisym(any_poroelastic, ATTENUATION_VISCOELASTIC_SOLID, anisotropic, &                  !axisym
    ROTATE_PML_ACTIVATE, STACEY_BOUNDARY_CONDITIONS, SIMULATION_TYPE, SAVE_FORWARD, &                                !axisym
    time_stepping_scheme, ADD_PERIODIC_CONDITIONS, NOISE_TOMOGRAPHY, NSOURCES, source_type, &                        !axisym
    ispec_selected_source, xi_source, anglesource, nrec, ispec_selected_rec, xi_receiver, nspec, is_on_the_axis, &   !axisym
    elastic, myrank, is_proc_source, which_proc_receiver)                                                            !axisym
  endif                                                                                                              !axisym

  !
  !---- build the global mass matrix
  !
  call invert_mass_matrix_init(any_elastic,any_acoustic,any_poroelastic, &
                                rmass_inverse_elastic_one,nglob_elastic, &
                                rmass_inverse_acoustic,nglob_acoustic, &
                                rmass_s_inverse_poroelastic, &
                                rmass_w_inverse_poroelastic,nglob_poroelastic, &
                                nspec,ibool,kmato,wxgll,wzgll,jacobian, &
                                elastic,poroelastic, &
                                assign_external_model,numat, &
                                density,poroelastcoef,porosity,tortuosity, &
                                vpext,rhoext,&
                                anyabs,numabs,deltat,codeabs,&
                                ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
                                ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
                                rmass_inverse_elastic_three,&
                                nelemabs,vsext,xix,xiz,gammaz,gammax, &
                                K_x_store,K_z_store,is_PML,&
                                AXISYM,nglob,is_on_the_axis,coord,wxglj,xiglj, &                                     !axisym
                                d_x_store,d_z_store,PML_BOUNDARY_CONDITIONS,region_CPML, &
                                nspec_PML,spec_to_PML,time_stepping_scheme)

#ifdef USE_MPI
  if ( nproc > 1 ) then

    ! preparing for MPI communications
    allocate(mask_ispec_inner_outer(nspec))
    mask_ispec_inner_outer(:) = .false.

    call get_MPI(nspec,ibool,knods,ngnod,nglob,elastic,poroelastic, &
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
                    myrank,coord)

    nspec_outer = count(mask_ispec_inner_outer)
    nspec_inner = nspec - nspec_outer

    allocate(ispec_outer_to_glob(nspec_outer))
    allocate(ispec_inner_to_glob(nspec_inner))

    ! building of corresponding arrays between inner/outer elements and their global number
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

    ! buffers for MPI communications
    max_ibool_interfaces_size_ac = maxval(nibool_interfaces_acoustic(:))
    max_ibool_interfaces_size_el = 3*maxval(nibool_interfaces_elastic(:))
    max_ibool_interfaces_size_po = NDIM*maxval(nibool_interfaces_poroelastic(:))
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

! assembling the mass matrix
    call assemble_MPI_scalar(rmass_inverse_acoustic,nglob_acoustic, &
                            rmass_inverse_elastic_one,rmass_inverse_elastic_three,nglob_elastic, &
                            rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic,nglob_poroelastic, &
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
    allocate(mask_ispec_inner_outer(1))

    nspec_outer = 0
    nspec_inner = nspec

    allocate(ispec_inner_to_glob(nspec_inner))
    do ispec = 1, nspec
      ispec_inner_to_glob(ispec) = ispec
    enddo

  endif ! end of test on whether there is more than one process (nproc > 1)

#else
  num_ispec_outer = 0
  num_ispec_inner = 0
  allocate(mask_ispec_inner_outer(1))

  nspec_outer = 0
  nspec_inner = nspec

  allocate(ispec_outer_to_glob(1))
  allocate(ispec_inner_to_glob(nspec_inner))
  do ispec = 1, nspec
     ispec_inner_to_glob(ispec) = ispec
  enddo

#endif

    ! loop over spectral elements
    do ispec_outer = 1,nspec_outer
    ! get global numbering for inner or outer elements
      ispec = ispec_outer_to_glob(ispec_outer)
    enddo

    ! loop over spectral elements
    do ispec_inner = 1,nspec_inner
    ! get global numbering for inner or outer elements
      ispec = ispec_inner_to_glob(ispec_inner)
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

    ! reduces cache misses for outer elements
    call get_global_indirect_addressing(nspec_outer,nglob,ibool_outer,copy_ibool_ori,integer_mask_ibool)

    ! the total number of points without multiples in this region is now known
    nglob_outer = maxval(ibool_outer)

    ! reduces cache misses for inner elements
    call get_global_indirect_addressing(nspec_inner,nglob,ibool_inner,copy_ibool_ori,integer_mask_ibool)

    ! the total number of points without multiples in this region is now known
    nglob_inner = maxval(ibool_inner)

  call invert_mass_matrix(any_elastic,any_acoustic,any_poroelastic,&
              rmass_inverse_elastic_one,rmass_inverse_elastic_three,nglob_elastic, &
              rmass_inverse_acoustic,nglob_acoustic, &
              rmass_s_inverse_poroelastic, &
              rmass_w_inverse_poroelastic,nglob_poroelastic)

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
                 coord,nglob,vpImin,vpImax,vpIImin,vpIImax, &
                 assign_external_model,nspec,UPPER_LIMIT_DISPLAY,numat,deltat, &
                 f0,initialfield,time_function_type, &
                 coorg,xinterp,zinterp,shape2D_display,knods,simulation_title, &
                 npgeo,pointsdisp,ngnod,any_elastic,any_poroelastic,all_anisotropic, &
                 myrank,nproc,NSOURCES,poroelastic, &
                 freq0,Q0,ATTENUATION_PORO_FLUID_PART,US_LETTER,output_postscript_snapshot)

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
                            coord,nglob,npgeo,factor_subsample_image)

    ! allocate an array for image data
    allocate(image_color_data(NX_IMAGE_color,NZ_IMAGE_color),stat=ier); if(ier /= 0) stop 'error in an allocate statement 1'
    allocate(image_color_vp_display(NX_IMAGE_color,NZ_IMAGE_color),stat=ier); if(ier /= 0) stop 'error in an allocate statement 2'

    ! allocate an array for the grid point that corresponds to a given image data point
    allocate(iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier); if(ier /= 0) stop 'error in an allocate statement 3'
    allocate(copy_iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier); if(ier /= 0) stop 'error in an allocate statement 4'

    ! creates pixels indexing
    call prepare_color_image_pixels(myrank,NX_IMAGE_color,NZ_IMAGE_color, &
                            xmin_color_image,xmax_color_image, &
                            zmin_color_image,zmax_color_image, &
                            coord,nglob,coorg,npgeo,nspec,ngnod,knods,ibool, &
                            nb_pixel_loc,iglob_image_color, &
                            DRAW_SOURCES_AND_RECEIVERS,NSOURCES,nrec,x_source,z_source,st_xval,st_zval, &
                            ix_image_color_source,iy_image_color_source,ix_image_color_receiver,iy_image_color_receiver)

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
                  iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
             do k = 1, nb_pixel_per_proc(iproc+1)
                j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
                i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color

                ! avoid edge effects
                if(i < 1) i = 1
                if(j < 1) j = 1

                if(i > NX_IMAGE_color) i = NX_IMAGE_color
                if(j > NZ_IMAGE_color) j = NZ_IMAGE_color

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
  sisux = ZERO ! double precision zero
  sisuz = ZERO

! initialize arrays to zero
  displ_elastic = 0._CUSTOM_REAL
  displ_elastic_old = 0._CUSTOM_REAL
  veloc_elastic = 0._CUSTOM_REAL
  accel_elastic = 0._CUSTOM_REAL

    if(SIMULATION_TYPE == 3 .and. any_elastic) then
      b_displ_elastic_old = 0._CUSTOM_REAL
      b_displ_elastic = 0._CUSTOM_REAL
      b_veloc_elastic = 0._CUSTOM_REAL
      b_accel_elastic = 0._CUSTOM_REAL
    endif

  if(time_stepping_scheme == 2) then
  displ_elastic_LDDRK = 0._CUSTOM_REAL
  veloc_elastic_LDDRK = 0._CUSTOM_REAL
  veloc_elastic_LDDRK_temp = 0._CUSTOM_REAL
  endif

  if(time_stepping_scheme == 3)then
   accel_elastic_rk = 0._CUSTOM_REAL
   veloc_elastic_rk = 0._CUSTOM_REAL
   veloc_elastic_initial_rk = 0._CUSTOM_REAL
   displ_elastic_initial_rk = 0._CUSTOM_REAL
  endif

  displs_poroelastic = 0._CUSTOM_REAL
  displs_poroelastic_old = 0._CUSTOM_REAL
  velocs_poroelastic = 0._CUSTOM_REAL
  accels_poroelastic = 0._CUSTOM_REAL
  displw_poroelastic = 0._CUSTOM_REAL
  velocw_poroelastic = 0._CUSTOM_REAL
  accelw_poroelastic = 0._CUSTOM_REAL

  if(time_stepping_scheme == 2) then
    displs_poroelastic_LDDRK = 0._CUSTOM_REAL
    velocs_poroelastic_LDDRK = 0._CUSTOM_REAL
    displw_poroelastic_LDDRK = 0._CUSTOM_REAL
    velocw_poroelastic_LDDRK = 0._CUSTOM_REAL
  endif

  if(time_stepping_scheme == 3) then
    accels_poroelastic_rk = 0._CUSTOM_REAL
    velocs_poroelastic_rk = 0._CUSTOM_REAL

    accelw_poroelastic_rk = 0._CUSTOM_REAL
    velocw_poroelastic_rk = 0._CUSTOM_REAL

    velocs_poroelastic_initial_rk = 0._CUSTOM_REAL
    displs_poroelastic_initial_rk = 0._CUSTOM_REAL

    velocw_poroelastic_initial_rk = 0._CUSTOM_REAL
    displw_poroelastic_initial_rk = 0._CUSTOM_REAL

  endif

  potential_acoustic = 0._CUSTOM_REAL
  potential_acoustic_old = 0._CUSTOM_REAL
  potential_dot_acoustic = 0._CUSTOM_REAL
  potential_dot_dot_acoustic = 0._CUSTOM_REAL

  if(time_stepping_scheme == 2 )then
    potential_acoustic_LDDRK = 0._CUSTOM_REAL
    potential_dot_acoustic_LDDRK = 0._CUSTOM_REAL
    potential_dot_acoustic_temp = 0._CUSTOM_REAL
  endif

  if(time_stepping_scheme == 3 )then
    potential_acoustic_init_rk = 0._CUSTOM_REAL
    potential_dot_acoustic_init_rk = 0._CUSTOM_REAL
    potential_dot_dot_acoustic_rk = 0._CUSTOM_REAL
    potential_dot_acoustic_rk = 0._CUSTOM_REAL
  endif

!
!----- Files where viscous damping are saved during forward wavefield calculation
!
  if(any_poroelastic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3)) then
    allocate(b_viscodampx(nglob))
    allocate(b_viscodampz(nglob))
    write(outputname,'(a,i6.6,a)') 'viscodampingx',myrank,'.bin'
    write(outputname2,'(a,i6.6,a)') 'viscodampingz',myrank,'.bin'
    if(SIMULATION_TYPE == 3) then
      reclen = CUSTOM_REAL * nglob
      open(unit=23,file='OUTPUT_FILES/'//outputname,status='old',&
            action='read',form='unformatted',access='direct',&
            recl=reclen)
      open(unit=24,file='OUTPUT_FILES/'//outputname2,status='old',&
            action='read',form='unformatted',access='direct',&
            recl=reclen)
    else
      reclen = CUSTOM_REAL * nglob
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

  if( ((SAVE_FORWARD .and. SIMULATION_TYPE ==1) .or. SIMULATION_TYPE == 3) .and. anyabs &
      .and. (.not. PML_BOUNDARY_CONDITIONS) ) then
    ! opens files for absorbing boundary data
    call prepare_absorb_files(myrank,any_elastic,any_poroelastic,any_acoustic, &
                      nspec_left,nspec_right,nspec_bottom,nspec_top,SIMULATION_TYPE)
  endif

  if(anyabs .and. SIMULATION_TYPE == 3 .and. (.not. PML_BOUNDARY_CONDITIONS)) then

    ! reads in absorbing boundary data
    if(any_elastic) then
      call prepare_absorb_elastic(NSTEP,p_sv, &
                      nspec_left,nspec_right,nspec_bottom,nspec_top, &
                      b_absorb_elastic_left,b_absorb_elastic_right, &
                      b_absorb_elastic_bottom,b_absorb_elastic_top)

    endif
    if(any_poroelastic) then
      call prepare_absorb_poroelastic(NSTEP, &
                      nspec_left,nspec_right,nspec_bottom,nspec_top, &
                      b_absorb_poro_s_left,b_absorb_poro_w_left, &
                      b_absorb_poro_s_right,b_absorb_poro_w_right, &
                      b_absorb_poro_s_bottom,b_absorb_poro_w_bottom, &
                      b_absorb_poro_s_top,b_absorb_poro_w_top)

    endif
    if(any_acoustic) then
      call prepare_absorb_acoustic(NSTEP, &
                      nspec_left,nspec_right,nspec_bottom,nspec_top, &
                      b_absorb_acoustic_left,b_absorb_acoustic_right, &
                      b_absorb_acoustic_bottom,b_absorb_acoustic_top)
    endif

  endif ! if(anyabs .and. SIMULATION_TYPE == 3)



!
!----- Read last frame for backward wavefield calculation
!

  if(SIMULATION_TYPE == 3) then

    if(any_elastic) then

      if(.not. save_ASCII_kernels)then
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_mu_kernel.bin'
        open(unit = 97, file = 'OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      else
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_mu_kernel.dat'
        open(unit = 97, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      endif

      if(.not. save_ASCII_kernels)then
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_alpha_beta_kernel.bin'
        open(unit = 98, file = 'OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      else
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_alpha_beta_kernel.dat'
        open(unit = 98, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      endif

      rho_kl(:,:,:) = 0._CUSTOM_REAL
      mu_kl(:,:,:) = 0._CUSTOM_REAL
      kappa_kl(:,:,:) = 0._CUSTOM_REAL

      rhop_kl(:,:,:) = 0._CUSTOM_REAL
      beta_kl(:,:,:) = 0._CUSTOM_REAL
      alpha_kl(:,:,:) = 0._CUSTOM_REAL
      rhorho_el_hessian_final2(:,:,:) = 0._CUSTOM_REAL
      rhorho_el_hessian_temp2(:) = 0._CUSTOM_REAL
      rhorho_el_hessian_final1(:,:,:) = 0._CUSTOM_REAL
      rhorho_el_hessian_temp1(:) = 0._CUSTOM_REAL
    endif

    if(any_poroelastic) then

      ! Primary kernels
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mu_B_C_kernel.dat'
      open(unit = 144, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_M_rho_rhof_kernel.dat'
      open(unit = 155, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_m_eta_kernel.dat'
      open(unit = 16, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      ! Wavespeed kernels
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_cpI_cpII_cs_kernel.dat'
      open(unit = 20, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhobb_rhofbb_ratio_kernel.dat'
      open(unit = 21, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_phib_eta_kernel.dat'
      open(unit = 22, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      ! Density normalized kernels
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mub_Bb_Cb_kernel.dat'
      open(unit = 17, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Mb_rhob_rhofb_kernel.dat'
      open(unit = 18, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mb_etab_kernel.dat'
      open(unit = 19, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'

      rhot_kl(:,:,:) = 0._CUSTOM_REAL
      rhof_kl(:,:,:) = 0._CUSTOM_REAL
      eta_kl(:,:,:) = 0._CUSTOM_REAL
      sm_kl(:,:,:) = 0._CUSTOM_REAL
      mufr_kl(:,:,:) = 0._CUSTOM_REAL
      B_kl(:,:,:) = 0._CUSTOM_REAL
      C_kl(:,:,:) = 0._CUSTOM_REAL
      M_kl(:,:,:) = 0._CUSTOM_REAL

      rhob_kl(:,:,:) = 0._CUSTOM_REAL
      rhofb_kl(:,:,:) = 0._CUSTOM_REAL
      phi_kl(:,:,:) = 0._CUSTOM_REAL
      mufrb_kl(:,:,:) = 0._CUSTOM_REAL
      Bb_kl(:,:,:) = 0._CUSTOM_REAL
      Cb_kl(:,:,:) = 0._CUSTOM_REAL
      Mb_kl(:,:,:) = 0._CUSTOM_REAL

      rhobb_kl(:,:,:) = 0._CUSTOM_REAL
      rhofbb_kl(:,:,:) = 0._CUSTOM_REAL
      phib_kl(:,:,:) = 0._CUSTOM_REAL
      cs_kl(:,:,:) = 0._CUSTOM_REAL
      cpI_kl(:,:,:) = 0._CUSTOM_REAL
      cpII_kl(:,:,:) = 0._CUSTOM_REAL
      ratio_kl(:,:,:) = 0._CUSTOM_REAL
    endif

    if(any_acoustic) then

      if(.not. save_ASCII_kernels)then
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_kernel.bin'
        open(unit = 95, file = 'OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',&
            iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      else
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_kernel.dat'
        open(unit = 95, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      endif

      if(.not. save_ASCII_kernels)then
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_c_kernel.bin'
        open(unit = 96, file = 'OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',&
            iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      else
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_c_kernel.dat'
        open(unit = 96, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      endif

      rho_ac_kl(:,:,:) = 0._CUSTOM_REAL
      kappa_ac_kl(:,:,:) = 0._CUSTOM_REAL

      rhop_ac_kl(:,:,:) = 0._CUSTOM_REAL
      alpha_ac_kl(:,:,:) = 0._CUSTOM_REAL
      rhorho_ac_hessian_final2(:,:,:) = 0._CUSTOM_REAL
      rhorho_ac_hessian_final1(:,:,:) = 0._CUSTOM_REAL
    endif

  endif ! if(SIMULATION_TYPE == 3)

!
!----  read initial fields from external file if needed
!

! if we are looking a plane wave beyond critical angle we use other method
  over_critical_angle = .false.

  if(initialfield) then

    ! Calculation of the initial field for a plane wave
    if( any_elastic ) then
      call prepare_initialfield(myrank,any_acoustic,any_poroelastic,over_critical_angle, &
                        NSOURCES,source_type,anglesource,x_source,z_source,f0,t0, &
                        nglob,numat,poroelastcoef,density,coord, &
                        anglesource_refl,c_inc,c_refl,cploc,csloc,time_offset, &
                        A_plane, B_plane, C_plane, &
                        accel_elastic,veloc_elastic,displ_elastic)
    endif

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

      ! call Paco's routine to compute in frequency and convert to time by Fourier transform
      call paco_beyond_critical(coord,nglob,deltat,NSTEP,anglesource(1),&
              f0(1),cploc,csloc,ATTENUATION_VISCOELASTIC_SOLID,QKappa_attenuation(1),source_type(1),v0x_left,v0z_left, &
              v0x_right,v0z_right,v0x_bot,v0z_bot,t0x_left,t0z_left,t0x_right,t0z_right, &
              t0x_bot,t0z_bot,left_bound(1:count_left),right_bound(1:count_right),bot_bound(1:count_bottom), &
              count_left,count_right,count_bottom,displ_elastic,veloc_elastic,accel_elastic,x_source(1))

      deallocate(left_bound)
      deallocate(right_bound)
      deallocate(bot_bound)

      if (myrank == 0) then
        write(IOUT,*)  '***********'
        write(IOUT,*)  'done calculating the initial wave field'
        write(IOUT,*)  '***********'
      endif

    endif ! beyond critical angle

    write(IOUT,*) 'Max norm of initial elastic displacement = ', &
      maxval(sqrt(displ_elastic(1,:)**2 + displ_elastic(3,:)**2))

  endif ! initialfield

! compute the source time function and store it in a text file
  if(.not. initialfield) then

    allocate(source_time_function(NSOURCES,NSTEP,stage_time_scheme))
    source_time_function(:,:,:) = 0._CUSTOM_REAL

    ! computes source time function array
    call prepare_source_time_function(myrank,NSTEP,NSOURCES,source_time_function, &
                          time_function_type,f0,tshift_src,factor,aval, &
                          t0,nb_proc_source,deltat,stage_time_scheme,c_LDDRK)
  else
    ! uses an initialfield
    ! dummy allocation
    allocate(source_time_function(1,1,1))
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
            ibegin_edge4(ispecabs) = 2
            ibegin_edge2(ispecabs) = 2
          endif

          if(iedge_acoustic == ITOP) then
            iend_edge4(ispecabs) = NGLLZ - 1
            iend_edge2(ispecabs) = NGLLZ - 1
          endif

          if(iedge_acoustic == ILEFT) then
            ibegin_edge1(ispecabs) = 2
            ibegin_edge3(ispecabs) = 2
          endif

          if(iedge_acoustic == IRIGHT) then
            iend_edge1(ispecabs) = NGLLX - 1
            iend_edge3(ispecabs) = NGLLX - 1
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
            ibegin_edge4(ispecabs) = 2
            ibegin_edge2(ispecabs) = 2
          endif

          if(iedge_acoustic == ITOP) then
            iend_edge4(ispecabs) = NGLLZ - 1
            iend_edge2(ispecabs) = NGLLZ - 1
          endif

          if(iedge_acoustic == ILEFT) then
            ibegin_edge1(ispecabs) = 2
            ibegin_edge3(ispecabs) = 2
          endif

          if(iedge_acoustic == IRIGHT) then
            iend_edge1(ispecabs) = NGLLX - 1
            iend_edge3(ispecabs) = NGLLX - 1
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

if(ATTENUATION_PORO_FLUID_PART .and. any_poroelastic .and. &
(time_stepping_scheme == 2.or. time_stepping_scheme == 3)) &
    stop 'RK and LDDRK time scheme not supported poroelastic simulations with attenuation'

if(coupled_elastic_poro) then

    if(ATTENUATION_VISCOELASTIC_SOLID .or. ATTENUATION_PORO_FLUID_PART) &
                   stop 'Attenuation not supported for mixed elastic/poroelastic simulations'

    if(time_stepping_scheme == 2.or. time_stepping_scheme == 3) &
                   stop 'RK and LDDRK time scheme not supported for mixed elastic/poroelastic simulations'

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
            ibegin_edge4_poro(ispecabs) = 1
            ibegin_edge2_poro(ispecabs) = 1

            iend_edge4_poro(ispecabs) = NGLLZ
            iend_edge2_poro(ispecabs) = NGLLZ

            ibegin_edge1_poro(ispecabs) = 1
            ibegin_edge3_poro(ispecabs) = 1

            iend_edge1_poro(ispecabs) = NGLLX
            iend_edge3_poro(ispecabs) = NGLLX
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
            ibegin_edge4_poro(ispecabs) = 2
            ibegin_edge2_poro(ispecabs) = 2
          endif

          if(iedge_poroelastic == ITOP) then
            iend_edge4_poro(ispecabs) = NGLLZ - 1
            iend_edge2_poro(ispecabs) = NGLLZ - 1
          endif

          if(iedge_poroelastic == ILEFT) then
            ibegin_edge1_poro(ispecabs) = 2
            ibegin_edge3_poro(ispecabs) = 2
          endif

          if(iedge_poroelastic == IRIGHT) then
            iend_edge1_poro(ispecabs) = NGLLX - 1
            iend_edge3_poro(ispecabs) = NGLLX - 1
          endif

        endif

      enddo

    enddo

  endif

!----  create a Gnuplot script to display the energy curve in log scale
  if(output_energy .and. myrank == 0) then
    close(IOUT_ENERGY)
    open(unit=IOUT_ENERGY,file='plot_energy.gnu',status='unknown',action='write')
    write(IOUT_ENERGY,*) 'set term wxt'
    write(IOUT_ENERGY,*) '#set term postscript landscape color solid "Helvetica" 22'
    write(IOUT_ENERGY,*) '#set output "energy.ps"'
    write(IOUT_ENERGY,*) '# set xrange [0:60]'
    write(IOUT_ENERGY,*) 'set logscale y'
    write(IOUT_ENERGY,*) 'set xlabel "Time (s)"'
    write(IOUT_ENERGY,*) 'set ylabel "Energy (J)"'
    write(IOUT_ENERGY,'(A)') &
      'plot "energy.dat" us 1:4 t ''Total Energy'' w l lc 1, "energy.dat" us 1:3 t ''Potential Energy'' w l lc 2'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    close(IOUT_ENERGY)
  endif

! open the file in which we will store the energy curve
  if(output_energy .and. myrank == 0) open(unit=IOUT_ENERGY,file='energy.dat',status='unknown',action='write')

!<NOISE_TOMOGRAPHY

  if (NOISE_TOMOGRAPHY /= 0) then

    !allocate arrays for noise tomography
    allocate(time_function_noise(NSTEP))
    allocate(source_array_noise(3,NGLLX,NGLLZ,NSTEP))
    allocate(mask_noise(nglob))
    allocate(surface_movie_x_noise(nglob))
    allocate(surface_movie_y_noise(nglob))
    allocate(surface_movie_z_noise(nglob))

    !read in parameters for noise tomography
    call read_parameters_noise(NOISE_TOMOGRAPHY,SIMULATION_TYPE,SAVE_FORWARD, &
                                 any_acoustic,any_poroelastic,p_sv,irec_master, &
                                 Mxx,Mxz,Mzz,factor,NSOURCES, &
                                 xi_receiver,gamma_receiver,ispec_selected_rec,nrec, &
                                 xi_noise,gamma_noise,ispec_noise,angle_noise)

  endif ! NOISE_TOMOGRAPHY /= 0


  if (NOISE_TOMOGRAPHY == 1) then
    call compute_source_array_noise(p_sv,NSTEP,deltat,nglob,ibool,ispec_noise, &
                                 xi_noise,gamma_noise,xigll,zigll, &
                                 time_function_noise,source_array_noise)

    !write out coordinates of mesh
    open(unit=504,file='OUTPUT_FILES/mesh_spec',status='unknown',action='write')
      do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            write(504,'(1pe11.3,1pe11.3,2i3,i7)') &
              coord(1,iglob), coord(2,iglob), i, j, ispec
         enddo
        enddo
      enddo
    close(504)

    open(unit=504,file='OUTPUT_FILES/mesh_glob',status='unknown',action='write')
      do iglob = 1, nglob
        write(504,'(1pe11.3,1pe11.3,i7)') &
          coord(1,iglob), coord(2,iglob), iglob
      enddo
    close(504)

    !write out spatial distribution of noise sources
    call create_mask_noise(nglob,coord,mask_noise)
    open(unit=504,file='OUTPUT_FILES/mask_noise',status='unknown',action='write')
      do iglob = 1, nglob
            write(504,'(1pe11.3,1pe11.3,1pe11.3)') &
              coord(1,iglob), coord(2,iglob), mask_noise(iglob)
      enddo
    close(504)

    !write out velocity model
    if(assign_external_model) then
      open(unit=504,file='OUTPUT_FILES/model_rho_vp_vs',status='unknown',action='write')
        do ispec = 1, nspec
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob = ibool(i,j,ispec)
              write(504,'(1pe11.3,1pe11.3,1pe11.3,1pe11.3,1pe11.3)') &
                coord(1,iglob), coord(2,iglob), &
                rhoext(i,j,ispec), vpext(i,j,ispec), vsext(i,j,ispec)
            enddo
          enddo
        enddo
      close(504)
    else
      open(unit=504,file='OUTPUT_FILES/model_rho_kappa_mu',status='unknown',action='write')
        do ispec = 1, nspec
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob = ibool(i,j,ispec)
              write(504,'(1pe11.3,1pe11.3,1pe11.3,1pe11.3,1pe11.3)') &
                coord(1,iglob), coord(2,iglob), density(1,kmato(ispec)), &
                poroelastcoef(1,1,kmato(ispec)) + 2.d0/3.d0*poroelastcoef(2,1,kmato(ispec)), &
                poroelastcoef(2,1,kmato(ispec))

            enddo
          enddo
        enddo
      close(504)
    endif

  else if (NOISE_TOMOGRAPHY == 2) then
    call create_mask_noise(nglob,coord,mask_noise)

  else if (NOISE_TOMOGRAPHY == 3) then

    if (output_wavefields_noise) then
      call create_mask_noise(nglob,coord,mask_noise)

      !prepare array that will hold wavefield snapshots
      noise_output_ncol = 5
      allocate(noise_output_array(noise_output_ncol,nglob))
      allocate(noise_output_rhokl(nglob))
    endif

  endif

!>NOISE_TOMOGRAPHY

  ! prepares image background
  if(output_color_image) then
    call prepare_color_image_vp(nglob,image_color_vp_display,iglob_image_color, &
                            NX_IMAGE_color,NZ_IMAGE_color,nb_pixel_loc, &
                            num_pixel_loc,nspec,elastic,poroelastic,ibool,kmato, &
                            numat,density,poroelastcoef,porosity,tortuosity, &
                            nproc,myrank,assign_external_model,vpext,DRAW_WATER_IN_BLUE)
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
  if(ATTENUATION_PORO_FLUID_PART) then
! viscous attenuation is implemented following the memory variable formulation of
! J. M. Carcione Wave fields in real media: wave propagation in anisotropic,
! anelastic and porous media, Elsevier, p. 304-305, 2007
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
    if(time_stepping_scheme == 2) then
     rx_viscous_LDDRK = 0.d0
     rz_viscous_LDDRK = 0.d0
    endif

    if(time_stepping_scheme == 3) then
     rx_viscous_initial_rk = 0.d0
     rz_viscous_initial_rk = 0.d0
     rx_viscous_force_RK = 0.d0
     rz_viscous_force_RK = 0.d0
    endif

  endif

! allocate arrays for postscript output
#ifdef USE_MPI
  if(modelvect) then
  d1_coorg_recv_ps_velocity_model=2
  call mpi_allreduce(nspec,d2_coorg_recv_ps_velocity_model,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
  d2_coorg_recv_ps_velocity_model=d2_coorg_recv_ps_velocity_model*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
       ((NGLLX-subsamp_postscript)/subsamp_postscript)*4
  d1_RGB_recv_ps_velocity_model=1
  call mpi_allreduce(nspec,d2_RGB_recv_ps_velocity_model,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
  d2_RGB_recv_ps_velocity_model=d2_RGB_recv_ps_velocity_model*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
       ((NGLLX-subsamp_postscript)/subsamp_postscript)*4
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
    d2_coorg_send_ps_vector_field=nglob
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
  d2_coorg_send_ps_velocity_model=nspec*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
                                        ((NGLLX-subsamp_postscript)/subsamp_postscript)*4
  d1_RGB_send_ps_velocity_model=1
  d2_RGB_send_ps_velocity_model=nspec*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
                                      ((NGLLX-subsamp_postscript)/subsamp_postscript)

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

! to dump the wave field
  this_is_the_first_time_we_dump = .true.

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

  ! get timestamp in minutes of current date and time
  year = time_values(1)
  mon = time_values(2)
  day = time_values(3)
  hr = time_values(5)
  minutes = time_values(6)
  call convtime(timestamp,year,mon,day,hr,minutes)

  ! convert to seconds instead of minutes, to be more precise for 2D runs, which can be fast
  timestamp_seconds_start = timestamp*60.d0 + time_values(7) + time_values(8)/1000.d0

! *********************************************************
! ************* MAIN LOOP OVER THE TIME STEPS *************
! *********************************************************

  do it = 1,NSTEP

! compute current time
    time = (it-1)*deltat

    do i_stage=1, stage_time_scheme

      call update_displacement_precondition_newmark(time_stepping_scheme,SIMULATION_TYPE,&
                                                    nglob_acoustic,nglob_elastic,nglob_poroelastic,&
                                                    any_acoustic,any_elastic,any_poroelastic,deltat,deltatover2,&
                                                    deltatsquareover2,potential_acoustic,potential_dot_acoustic,&
                                                    potential_dot_dot_acoustic,potential_acoustic_old,&
                                                    displ_elastic,veloc_elastic,accel_elastic,displ_elastic_old,&
                                                    displs_poroelastic,velocs_poroelastic,accels_poroelastic,&
                                                    displs_poroelastic_old,displw_poroelastic,velocw_poroelastic,&
                                                    accelw_poroelastic,b_deltat,b_deltatover2,b_deltatsquareover2,&
                                                    b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic,&
                                                    b_potential_acoustic_old,&
                                                    b_displ_elastic,b_veloc_elastic,b_accel_elastic,b_displ_elastic_old,&
                                                    accel_elastic_adj_coupling,&
                                                    b_displs_poroelastic,b_velocs_poroelastic,b_accels_poroelastic,&
                                                    accels_poroelastic_adj_coupling,&
                                                    b_displw_poroelastic,b_velocw_poroelastic,b_accelw_poroelastic,&
                                                    accelw_poroelastic_adj_coupling)

      if (AXISYM) then                                                    !axisym
        do ispec=1,nspec                                                  !axisym
          if (elastic(ispec) .and. is_on_the_axis(ispec)) then            !axisym
            do j = 1,NGLLZ                                                !axisym
              do i = 1,NGLJ                                               !axisym
                if (abs(coord(1,ibool(i,j,ispec))) < TINYVAL) then        !axisym
                  displ_elastic(1,ibool(i,j,ispec))=ZERO                  !axisym
                endif                                                     !axisym
              enddo                                                       !axisym
            enddo                                                         !axisym
          endif                                                           !axisym
        enddo                                                             !axisym
      endif                                                               !axisym

    if(any_acoustic) then
      ! free surface for an acoustic medium
      if ( nelem_acoustic_surface > 0 ) then
        call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                          potential_acoustic,acoustic_surface, &
                                          ibool,nelem_acoustic_surface,nglob,nspec,this_ibool_is_a_periodic_edge)

        if(SIMULATION_TYPE == 3) then ! Adjoint calculation
          call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                            b_potential_acoustic,acoustic_surface, &
                                            ibool,nelem_acoustic_surface,nglob,nspec,this_ibool_is_a_periodic_edge)
        endif
      endif

! *********************************************************
! ************* compute forces for the acoustic elements
! *********************************************************

      call compute_forces_acoustic(nglob,nspec,nelemabs,numat,it,NSTEP, &
               anyabs,assign_external_model,ibool,kmato,numabs, &
               elastic,poroelastic,codeabs,potential_dot_dot_acoustic,potential_dot_acoustic, &
               potential_acoustic,potential_acoustic_old,stage_time_scheme,i_stage, &
               density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,rhoext,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               AXISYM,coord, is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj, &                             !axisym
               ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
               ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
               SIMULATION_TYPE,SAVE_FORWARD,nspec_left,nspec_right,&
               nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top, &
               b_absorb_acoustic_left,b_absorb_acoustic_right, &
               b_absorb_acoustic_bottom,b_absorb_acoustic_top,.false.,&
               is_PML,nspec_PML,spec_to_PML,region_CPML, &
               K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store,&
               rmemory_potential_acoustic,&
               rmemory_acoustic_dux_dx,rmemory_acoustic_dux_dz,&
               rmemory_potential_acoustic_LDDRK,alpha_LDDRK,beta_LDDRK,c_LDDRK, &
               rmemory_acoustic_dux_dx_LDDRK,rmemory_acoustic_dux_dz_LDDRK,&
               deltat,PML_BOUNDARY_CONDITIONS,STACEY_BOUNDARY_CONDITIONS)

      if( SIMULATION_TYPE == 3 ) then

       if(PML_BOUNDARY_CONDITIONS)then
          do ispec = 1,nspec
            do i = 1, NGLLX
              do j = 1, NGLLZ
                if(.not. elastic(ispec) .and. .not. poroelastic(ispec) .and. is_pml(ispec))then
                  b_potential_dot_dot_acoustic(ibool(i,j,ispec)) = 0.
                  b_potential_dot_acoustic(ibool(i,j,ispec)) = 0.
                  b_potential_acoustic(ibool(i,j,ispec)) = 0.
                endif
               enddo
            enddo
          enddo
       endif

       if(PML_BOUNDARY_CONDITIONS)then
         if(any_acoustic .and. nglob_interface > 0)then
           do i = 1, nglob_interface
             b_potential_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot(i,NSTEP-it+1)
             b_potential_acoustic(point_interface(i)) = pml_interface_history_potential(i,NSTEP-it+1)
           enddo
         endif
       endif

        call compute_forces_acoustic(nglob,nspec,nelemabs,numat,it,NSTEP, &
               anyabs,assign_external_model,ibool,kmato,numabs, &
               elastic,poroelastic,codeabs,b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
               b_potential_acoustic,b_potential_acoustic_old,stage_time_scheme, i_stage, &
               density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,rhoext,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               AXISYM,coord, is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj, &                             !axisym
               ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
               ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
               SIMULATION_TYPE,SAVE_FORWARD,nspec_left,nspec_right,&
               nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top, &
               b_absorb_acoustic_left,b_absorb_acoustic_right, &
               b_absorb_acoustic_bottom,b_absorb_acoustic_top,.true.,&
               is_PML,nspec_PML,spec_to_PML,region_CPML, &
               K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store,&
               rmemory_potential_acoustic,&
               rmemory_acoustic_dux_dx,rmemory_acoustic_dux_dz,&
               rmemory_potential_acoustic_LDDRK,alpha_LDDRK,beta_LDDRK,c_LDDRK, &
               rmemory_acoustic_dux_dx_LDDRK,rmemory_acoustic_dux_dz_LDDRK,&
!               deltat,PML_BOUNDARY_CONDITIONS)
               deltat,.false.,STACEY_BOUNDARY_CONDITIONS)

       if(PML_BOUNDARY_CONDITIONS)then
          do ispec = 1,nspec
            do i = 1, NGLLX
              do j = 1, NGLLZ
                if(.not. elastic(ispec) .and. .not. poroelastic(ispec) .and. is_pml(ispec))then
                  b_potential_dot_acoustic(ibool(i,j,ispec)) = 0.
                  b_potential_acoustic(ibool(i,j,ispec)) = 0.
                endif
               enddo
            enddo
          enddo
       endif

       if(PML_BOUNDARY_CONDITIONS)then
         if(any_acoustic .and. nglob_interface > 0)then
           do i = 1, nglob_interface
             b_potential_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot(i,NSTEP-it+1)
             b_potential_acoustic(point_interface(i)) = pml_interface_history_potential(i,NSTEP-it+1)
           enddo
         endif
       endif

      endif


      ! stores absorbing boundary contributions into files
      if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. PML_BOUNDARY_CONDITIONS)) then
        !--- left absorbing boundary
        if(nspec_left >0) then
          do ispec = 1,nspec_left
            do i=1,NGLLZ
              write(65) b_absorb_acoustic_left(i,ispec,it)
            enddo
          enddo
        endif
        !--- right absorbing boundary
        if(nspec_right >0) then
          do ispec = 1,nspec_right
            do i=1,NGLLZ
              write(66) b_absorb_acoustic_right(i,ispec,it)
            enddo
          enddo
        endif
        !--- bottom absorbing boundary
        if(nspec_bottom >0) then
          do ispec = 1,nspec_bottom
            do i=1,NGLLX
              write(67) b_absorb_acoustic_bottom(i,ispec,it)
            enddo
          enddo
        endif
        !--- top absorbing boundary
        if(nspec_top >0) then
          do ispec = 1,nspec_top
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

      if(SIMULATION_TYPE == 1)then
        call compute_coupling_acoustic_el(nspec,nglob_elastic,nglob_acoustic,num_fluid_solid_edges,ibool,wxgll,wzgll,xix,xiz,&
               gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse,displ_elastic,displ_elastic_old,&
               potential_dot_dot_acoustic,fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
               fluid_solid_elastic_ispec,fluid_solid_elastic_iedge,&
               AXISYM,nglob,coord,is_on_the_axis,xiglj,wxglj, &                                               !axisym
               PML_BOUNDARY_CONDITIONS,nspec_PML,K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,&
               alpha_z_store,is_PML,spec_to_PML,region_CPML,rmemory_fsb_displ_elastic,time,deltat)
      endif

      if(SIMULATION_TYPE == 3)then
        call compute_coupling_acoustic_el(nspec,nglob_elastic,nglob_acoustic,num_fluid_solid_edges,ibool,wxgll,wzgll,xix,xiz,&
               gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse,-accel_elastic_adj_coupling,displ_elastic_old,&
               potential_dot_dot_acoustic,fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
               fluid_solid_elastic_ispec,fluid_solid_elastic_iedge,&
               AXISYM,nglob,coord,is_on_the_axis,xiglj,wxglj, &                                               !axisym
               PML_BOUNDARY_CONDITIONS,nspec_PML,K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,&
               alpha_z_store,is_PML,spec_to_PML,region_CPML,rmemory_fsb_displ_elastic,time,deltat)

        call compute_coupling_acoustic_el(nspec,nglob_elastic,nglob_acoustic,num_fluid_solid_edges,ibool,wxgll,wzgll,xix,xiz,&
               gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse,b_displ_elastic,b_displ_elastic_old,&
               b_potential_dot_dot_acoustic,fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
               fluid_solid_elastic_ispec,fluid_solid_elastic_iedge,&
               AXISYM,nglob,coord,is_on_the_axis,xiglj,wxglj, &                                               !axisym
               PML_BOUNDARY_CONDITIONS,nspec_PML,K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,&
               alpha_z_store,is_PML,spec_to_PML,region_CPML,rmemory_fsb_displ_elastic,time,deltat)
      endif

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

          if(SIMULATION_TYPE == 3) then
            b_displ_x = b_displs_poroelastic(1,iglob)
            b_displ_z = b_displs_poroelastic(2,iglob)

            b_displw_x = b_displw_poroelastic(1,iglob)
            b_displw_z = b_displw_poroelastic(2,iglob)

            ! new definition of adjoint displacement and adjoint potential
            displ_x = -accels_poroelastic_adj_coupling(1,iglob)
            displ_z = -accels_poroelastic_adj_coupling(2,iglob)

            displw_x = -accelw_poroelastic_adj_coupling(1,iglob)
            displw_z = -accelw_poroelastic_adj_coupling(2,iglob)
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
          else if(iedge_acoustic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_acoustic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          else if(iedge_acoustic ==IRIGHT)then
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

          if(SIMULATION_TYPE == 3) then
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
          ! if this processor core carries the source and the source element is acoustic
          if (is_proc_source(i_source) == 1 .and. (.not. elastic(ispec_selected_source(i_source))) .and. &
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
                                            - source_time_function(i_source,it,i_stage)*hlagrange
                  enddo
                enddo
              else
                ! backward wavefield
                do j = 1,NGLLZ
                  do i = 1,NGLLX
                    iglob = ibool(i,j,ispec_selected_source(i_source))
                    hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                    b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) &
                                          - source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)*hlagrange
                  enddo
                enddo
              endif

            ! moment tensor
            else if(source_type(i_source) == 2) then
              call exit_MPI('cannot have moment tensor source in acoustic element')

            endif
          endif ! if this processor core carries the source and the source element is acoustic
        enddo ! do i_source=1,NSOURCES

        if(SIMULATION_TYPE == 3) then   ! adjoint wavefield
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
                                  + adj_sourcearrays(irec_local,NSTEP-it+1,1,i,j)
                  enddo
                enddo
              endif ! if element acoustic

            endif ! if this processor core carries the adjoint source
          enddo ! irec = 1,nrec
        endif ! SIMULATION_TYPE == 3 adjoint wavefield

      endif ! if not using an initial field

    endif !if(any_acoustic)


! assembling potential_dot_dot for acoustic elements
#ifdef USE_MPI
    if ( nproc > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
      call assemble_MPI_vector_ac(potential_dot_dot_acoustic,nglob, &
                    ninterface, ninterface_acoustic,inum_interfaces_acoustic, &
                    max_interface_size, max_ibool_interfaces_size_ac,&
                    ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
                    tab_requests_send_recv_acoustic,buffer_send_faces_vector_ac, &
                    buffer_recv_faces_vector_ac, my_neighbours)

     if(time_stepping_scheme == 2)then
      if(i_stage==1 .and. it == 1 .and. (.not. initialfield))then
       potential_dot_acoustic_temp = potential_dot_acoustic
       call assemble_MPI_vector_ac(potential_dot_acoustic,nglob, &
                    ninterface, ninterface_acoustic,inum_interfaces_acoustic, &
                    max_interface_size, max_ibool_interfaces_size_ac,&
                    ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
                    tab_requests_send_recv_acoustic,buffer_send_faces_vector_ac, &
                    buffer_recv_faces_vector_ac, my_neighbours)
      endif
     endif

      if ( SIMULATION_TYPE == 3) then
        call assemble_MPI_vector_ac(b_potential_dot_dot_acoustic,nglob, &
                     ninterface, ninterface_acoustic,inum_interfaces_acoustic, &
                     max_interface_size, max_ibool_interfaces_size_ac,&
                     ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
                     tab_requests_send_recv_acoustic,buffer_send_faces_vector_ac, &
                     buffer_recv_faces_vector_ac, my_neighbours)

      endif

    endif

#endif

      if(PML_BOUNDARY_CONDITIONS .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1)then
       if(any_acoustic .and. nglob_interface > 0)then
        do i = 1, nglob_interface
          write(72)potential_dot_dot_acoustic(point_interface(i)),&
                   potential_dot_acoustic(point_interface(i)),&
                   potential_acoustic(point_interface(i))
        enddo
       endif
      endif

     if(SIMULATION_TYPE == 3)then
       if(PML_BOUNDARY_CONDITIONS)then
         if(any_acoustic .and. nglob_interface > 0)then
           do i = 1, nglob_interface
             b_potential_dot_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot_dot(i,NSTEP-it+1)
           enddo
         endif
       endif
     endif

! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

    if(any_acoustic) then
      if(time_stepping_scheme == 1)then
!! DK DK this should be vectorized
      potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic
      potential_dot_acoustic = potential_dot_acoustic + deltatover2*potential_dot_dot_acoustic
      endif

      if(time_stepping_scheme == 2)then

!! DK DK this should be vectorized
        potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic

        potential_dot_acoustic_LDDRK = alpha_LDDRK(i_stage) * potential_dot_acoustic_LDDRK &
                                       + deltat * potential_dot_dot_acoustic

        potential_acoustic_LDDRK = alpha_LDDRK(i_stage) * potential_acoustic_LDDRK &
                                       +deltat*potential_dot_acoustic

        if(i_stage==1 .and. it == 1 .and. (.not. initialfield))then
!! DK DK this should be vectorized
        potential_dot_acoustic_temp = potential_dot_acoustic_temp &
                                      + beta_LDDRK(i_stage) * potential_dot_acoustic_LDDRK
        potential_dot_acoustic = potential_dot_acoustic_temp
        else
        potential_dot_acoustic = potential_dot_acoustic + beta_LDDRK(i_stage) * potential_dot_acoustic_LDDRK
        endif

!! DK DK this should be vectorized
        potential_acoustic = potential_acoustic + beta_LDDRK(i_stage) * potential_acoustic_LDDRK

      endif

      if(time_stepping_scheme == 3)then

!! DK DK this should be vectorized
        potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic

        potential_dot_dot_acoustic_rk(:,i_stage) = deltat * potential_dot_dot_acoustic(:)
        potential_dot_acoustic_rk(:,i_stage) = deltat * potential_dot_acoustic(:)

        if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

        if(i_stage == 1)weight_rk = 0.5d0
        if(i_stage == 2)weight_rk = 0.5d0
        if(i_stage == 3)weight_rk = 1.0d0

        if(i_stage==1)then

!! DK DK this should be vectorized
        potential_dot_acoustic_init_rk = potential_dot_acoustic
        potential_acoustic_init_rk = potential_acoustic

        endif

!! DK DK this should be vectorized
        potential_dot_acoustic(:) = potential_dot_acoustic_init_rk(:) + weight_rk * potential_dot_dot_acoustic_rk(:,i_stage)
        potential_acoustic(:) = potential_acoustic_init_rk(:) + weight_rk * potential_dot_acoustic_rk(:,i_stage)

        else if(i_stage==4)then

!! DK DK this should be vectorized
        potential_dot_acoustic(:) = potential_dot_acoustic_init_rk(:) + 1.0d0 / 6.0d0 * &
        (potential_dot_dot_acoustic_rk(:,1) + 2.0d0 * potential_dot_dot_acoustic_rk(:,2) + &
         2.0d0 * potential_dot_dot_acoustic_rk(:,3) + potential_dot_dot_acoustic_rk(:,4))

!! DK DK this should be vectorized
        potential_acoustic(:) = potential_acoustic_init_rk(:) + 1.0d0 / 6.0d0 * &
        (potential_dot_acoustic_rk(:,1) + 2.0d0 * potential_dot_acoustic_rk(:,2) + &
         2.0d0 * potential_dot_acoustic_rk(:,3) + potential_dot_acoustic_rk(:,4))

        endif

      endif

      if(SIMULATION_TYPE == 3)then
!! DK DK this should be vectorized
        b_potential_dot_dot_acoustic = b_potential_dot_dot_acoustic * rmass_inverse_acoustic
        b_potential_dot_acoustic = b_potential_dot_acoustic + b_deltatover2*b_potential_dot_dot_acoustic
      endif


      ! free surface for an acoustic medium
      if ( nelem_acoustic_surface > 0 ) then
        call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                        potential_acoustic,acoustic_surface, &
                                        ibool,nelem_acoustic_surface,nglob,nspec,this_ibool_is_a_periodic_edge)

        if(SIMULATION_TYPE == 3) then
          call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                          b_potential_acoustic,acoustic_surface, &
                                          ibool,nelem_acoustic_surface,nglob,nspec,this_ibool_is_a_periodic_edge)
        endif

      endif

      ! update the potential field (use a new array here) for coupling terms
      potential_acoustic_adj_coupling = potential_acoustic &
                          + deltat*potential_dot_acoustic &
                          + deltatsquareover2*potential_dot_dot_acoustic

    endif !if(any_acoustic)


! *********************************************************
! ************* main solver for the elastic elements
! *********************************************************

    if(any_elastic) then

      call compute_forces_viscoelastic(p_sv,nglob,nspec,myrank,nelemabs,numat, &
               ispec_selected_source,ispec_selected_rec,is_proc_source,which_proc_receiver, &
               source_type,it,NSTEP,anyabs,assign_external_model, &
               initialfield,ATTENUATION_VISCOELASTIC_SOLID,anglesource, &
               ibool,kmato,numabs,elastic,codeabs, &
               accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old, &
               density,poroelastcoef,xix,xiz,gammax,gammaz, &
               jacobian,vpext,vsext,rhoext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,anisotropic,anisotropy, &
               source_time_function,sourcearray,adj_sourcearrays, &
               e1,e11,e13,e1_LDDRK,e11_LDDRK,e13_LDDRK,alpha_LDDRK,beta_LDDRK,c_LDDRK, &
               e1_initial_rk,e11_initial_rk,e13_initial_rk,e1_force_rk, e11_force_rk, e13_force_rk, &
               hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               AXISYM,is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj, &                                   !axisym
               inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,Mu_nu1,Mu_nu2,N_SLS, &
               deltat,coord,add_Bielak_conditions, x_source(1), z_source(1), &
               A_plane, B_plane, C_plane, anglesource_refl, c_inc, c_refl, time_offset, f0(1),&
               v0x_left(1,it),v0z_left(1,it),v0x_right(1,it),v0z_right(1,it),v0x_bot(1,it),v0z_bot(1,it), &
               t0x_left(1,it),t0z_left(1,it),t0x_right(1,it),t0z_right(1,it),t0x_bot(1,it),t0z_bot(1,it), &
               count_left,count_right,count_bottom,over_critical_angle,NSOURCES,nrec,SIMULATION_TYPE,SAVE_FORWARD, &
               b_absorb_elastic_left,b_absorb_elastic_right,b_absorb_elastic_bottom,b_absorb_elastic_top, &
               nspec_left,nspec_right,nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top, &
               stage_time_scheme,i_stage,ADD_SPRING_TO_STACEY,x_center_spring,z_center_spring,max(1,nadj_rec_local), &
               is_PML,nspec_PML,spec_to_PML,region_CPML, &
               K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store, &
               rmemory_displ_elastic,rmemory_dux_dx,rmemory_dux_dz,rmemory_duz_dx,rmemory_duz_dz, &
               rmemory_dux_dx_prime,rmemory_dux_dz_prime,rmemory_duz_dx_prime,rmemory_duz_dz_prime, &
               rmemory_displ_elastic_LDDRK,rmemory_dux_dx_LDDRK,rmemory_dux_dz_LDDRK,&
               rmemory_duz_dx_LDDRK,rmemory_duz_dz_LDDRK, &
               PML_BOUNDARY_CONDITIONS,ROTATE_PML_ACTIVATE,ROTATE_PML_ANGLE,.false.,STACEY_BOUNDARY_CONDITIONS)

      if(SIMULATION_TYPE == 3)then
       if(PML_BOUNDARY_CONDITIONS)then
          do ispec = 1,nspec
            do i = 1, NGLLX
              do j = 1, NGLLZ
                if(elastic(ispec) .and. is_pml(ispec))then
                  b_veloc_elastic(:,ibool(i,j,ispec)) = 0.
                  b_displ_elastic(:,ibool(i,j,ispec)) = 0.
                endif
               enddo
            enddo
          enddo
       endif

       if(PML_BOUNDARY_CONDITIONS)then
         if(any_elastic .and. nglob_interface > 0)then
           do i = 1, nglob_interface
             b_veloc_elastic(1,point_interface(i)) = pml_interface_history_veloc(1,i,NSTEP-it+1)
             b_veloc_elastic(2,point_interface(i)) = pml_interface_history_veloc(2,i,NSTEP-it+1)
             b_veloc_elastic(3,point_interface(i)) = pml_interface_history_veloc(3,i,NSTEP-it+1)
             b_displ_elastic(1,point_interface(i)) = pml_interface_history_displ(1,i,NSTEP-it+1)
             b_displ_elastic(2,point_interface(i)) = pml_interface_history_displ(2,i,NSTEP-it+1)
             b_displ_elastic(3,point_interface(i)) = pml_interface_history_displ(3,i,NSTEP-it+1)
           enddo
         endif
       endif

      call compute_forces_viscoelastic(p_sv,nglob,nspec,myrank,nelemabs,numat, &
               ispec_selected_source,ispec_selected_rec,is_proc_source,which_proc_receiver, &
               source_type,it,NSTEP,anyabs,assign_external_model, &
               initialfield,ATTENUATION_VISCOELASTIC_SOLID,anglesource, &
               ibool,kmato,numabs,elastic,codeabs, &
               b_accel_elastic,b_veloc_elastic,b_displ_elastic,b_displ_elastic_old, &
               density,poroelastcoef,xix,xiz,gammax,gammaz, &
               jacobian,vpext,vsext,rhoext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,anisotropic,anisotropy, &
               source_time_function,sourcearray,adj_sourcearrays, &
               e1,e11,e13,e1_LDDRK,e11_LDDRK,e13_LDDRK,alpha_LDDRK,beta_LDDRK,c_LDDRK, &
               e1_initial_rk,e11_initial_rk,e13_initial_rk,e1_force_rk, e11_force_rk, e13_force_rk, &
               hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               AXISYM,is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj, &                                   !axisym
               inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,Mu_nu1,Mu_nu2,N_SLS, &
               deltat,coord,add_Bielak_conditions, x_source(1), z_source(1), &
               A_plane, B_plane, C_plane, anglesource_refl, c_inc, c_refl, time_offset, f0(1),&
               v0x_left(1,it),v0z_left(1,it),v0x_right(1,it),v0z_right(1,it),v0x_bot(1,it),v0z_bot(1,it), &
               t0x_left(1,it),t0z_left(1,it),t0x_right(1,it),t0z_right(1,it),t0x_bot(1,it),t0z_bot(1,it), &
               count_left,count_right,count_bottom,over_critical_angle, &
               NSOURCES,nrec,SIMULATION_TYPE,SAVE_FORWARD, &
               b_absorb_elastic_left,b_absorb_elastic_right,b_absorb_elastic_bottom,b_absorb_elastic_top, &
               nspec_left,nspec_right,nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top, &
               stage_time_scheme,i_stage,ADD_SPRING_TO_STACEY,x_center_spring,z_center_spring,max(1,nadj_rec_local), &
               is_PML,nspec_PML,spec_to_PML,region_CPML, &
               K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store, &
               rmemory_displ_elastic,rmemory_dux_dx,rmemory_dux_dz,rmemory_duz_dx,rmemory_duz_dz, &
               rmemory_dux_dx_prime,rmemory_dux_dz_prime,rmemory_duz_dx_prime,rmemory_duz_dz_prime, &
               rmemory_displ_elastic_LDDRK,rmemory_dux_dx_LDDRK,rmemory_dux_dz_LDDRK,&
               rmemory_duz_dx_LDDRK,rmemory_duz_dz_LDDRK, &
               .false.,ROTATE_PML_ACTIVATE,ROTATE_PML_ANGLE,.true.,STACEY_BOUNDARY_CONDITIONS)

       if(PML_BOUNDARY_CONDITIONS)then
          do ispec = 1,nspec
            do i = 1, NGLLX
              do j = 1, NGLLZ
                if(elastic(ispec) .and. is_pml(ispec))then
                  b_veloc_elastic(:,ibool(i,j,ispec)) = 0.
                  b_displ_elastic(:,ibool(i,j,ispec)) = 0.
                endif
               enddo
            enddo
          enddo
       endif

       if(PML_BOUNDARY_CONDITIONS)then
         if(any_elastic .and. nglob_interface > 0)then
           do i = 1, nglob_interface
             b_veloc_elastic(1,point_interface(i)) = pml_interface_history_veloc(1,i,NSTEP-it+1)
             b_veloc_elastic(2,point_interface(i)) = pml_interface_history_veloc(2,i,NSTEP-it+1)
             b_veloc_elastic(3,point_interface(i)) = pml_interface_history_veloc(3,i,NSTEP-it+1)
             b_displ_elastic(1,point_interface(i)) = pml_interface_history_displ(1,i,NSTEP-it+1)
             b_displ_elastic(2,point_interface(i)) = pml_interface_history_displ(2,i,NSTEP-it+1)
             b_displ_elastic(3,point_interface(i)) = pml_interface_history_displ(3,i,NSTEP-it+1)
           enddo
         endif
       endif


      call compute_forces_viscoelastic_pre_kernel(p_sv,nglob,nspec,displ_elastic,b_displ_elastic,&
              mu_k,kappa_k,elastic,ibool,hprime_xx,hprime_zz,xix,xiz,gammax,gammaz,SIMULATION_TYPE)
      endif


      if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. PML_BOUNDARY_CONDITIONS)) then
        !--- left absorbing boundary
        if(nspec_left >0) then
          do ispec = 1,nspec_left
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
        if(nspec_right >0) then
          do ispec = 1,nspec_right
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
        if(nspec_bottom >0) then
          do ispec = 1,nspec_bottom
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
        if(nspec_top >0) then
          do ispec = 1,nspec_top
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
          if(SIMULATION_TYPE == 3) then
            b_pressure = - b_potential_dot_dot_acoustic(iglob)
            !<YANGL
            ! new definition of adjoint displacement and adjoint potential
            pressure = potential_acoustic_adj_coupling(iglob)
            !>YANGL
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

          if (AXISYM) then                                                                                           !axisym
            if (abs(coord(1,ibool(i,j,ispec_acoustic))) < TINYVAL) then                                              !axisym
              xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)                                      !axisym
              r_xiplus1(i,j) = xxi                                                                                   !axisym
            else if (is_on_the_axis(ispec_acoustic)) then                                                            !axisym
               r_xiplus1(i,j) = coord(1,ibool(i,j,ispec_acoustic))/(xiglj(i)+ONE)                                    !axisym
            endif                                                                                                    !axisym
          endif                                                                                                      !axisym

          if(iedge_acoustic == ITOP)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D
            if (AXISYM) then                                                                                         !axisym
              if (is_on_the_axis(ispec_acoustic)) then                                                               !axisym
                weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)                                                      !axisym
              else                                                                                                   !axisym
                weight = jacobian1D * wxgll(i) * coord(1,ibool(i,j,ispec_acoustic))                                  !axisym
              endif                                                                                                  !axisym
            else                                                                                                     !axisym
              weight = jacobian1D * wxgll(i)
            endif                                                                                                    !axisym
          else if(iedge_acoustic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
            if (AXISYM) then                                                                                         !axisym
              if (is_on_the_axis(ispec_acoustic)) then                                                               !axisym
                weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)                                                      !axisym
              else                                                                                                   !axisym
                weight = jacobian1D * wxgll(i) * coord(1,ibool(i,j,ispec_acoustic))                                  !axisym
              endif                                                                                                  !axisym
            else                                                                                                     !axisym
              weight = jacobian1D * wxgll(i)
            endif                                                                                                    !axisym
          else if(iedge_acoustic == ILEFT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          else if(iedge_acoustic == IRIGHT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          endif

          accel_elastic(1,iglob) = accel_elastic(1,iglob) + weight*nx*pressure
          accel_elastic(3,iglob) = accel_elastic(3,iglob) + weight*nz*pressure

          if(SIMULATION_TYPE == 3) then
            b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) + weight*nx*b_pressure
            b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) + weight*nz*b_pressure
          endif !if(SIMULATION_TYPE == 3) then

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

          if(SIMULATION_TYPE == 3) then
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
            if(SIMULATION_TYPE == 3) then
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

          if(SIMULATION_TYPE == 3) then
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

          if(SIMULATION_TYPE == 3) then
            b_sigma_xx = lambdalplus2mul_G*b_dux_dxl + lambdal_G*b_duz_dzl + C_biot*(b_dwx_dxl + b_dwz_dzl)
            b_sigma_xz = mul_G*(b_duz_dxl + b_dux_dzl)
            b_sigma_zz = lambdalplus2mul_G*b_duz_dzl + lambdal_G*b_dux_dxl + C_biot*(b_dwx_dxl + b_dwz_dzl)
          endif
          ! get point values for the elastic domain, which matches our side in the inverse direction
          ii2 = ivalue(ipoin1D,iedge_elastic)
          jj2 = jvalue(ipoin1D,iedge_elastic)
          iglob = ibool(ii2,jj2,ispec_elastic)

          ! get elastic properties
          lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec_elastic))
          mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec_elastic))
          lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec_elastic))

          ! derivative along x and along z for u_s and w
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          if(SIMULATION_TYPE == 3) then
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

            if(SIMULATION_TYPE == 3) then
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

          if(SIMULATION_TYPE == 3) then
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
              c12 = c12ext(ii2,jj2,ispec_elastic)
              c23 = c23ext(ii2,jj2,ispec_elastic)
              c25 = c25ext(ii2,jj2,ispec_elastic)
            else
              c11 = anisotropy(1,kmato(ispec_elastic))
              c13 = anisotropy(2,kmato(ispec_elastic))
              c15 = anisotropy(3,kmato(ispec_elastic))
              c33 = anisotropy(4,kmato(ispec_elastic))
              c35 = anisotropy(5,kmato(ispec_elastic))
              c55 = anisotropy(6,kmato(ispec_elastic))
              c12 = anisotropy(7,kmato(ispec_elastic))
              c23 = anisotropy(8,kmato(ispec_elastic))
              c25 = anisotropy(9,kmato(ispec_elastic))
            endif

            sigma_xx = sigma_xx + c11*dux_dxl + c15*(duz_dxl + dux_dzl) + c13*duz_dzl
            sigma_zz = sigma_zz + c13*dux_dxl + c35*(duz_dxl + dux_dzl) + c33*duz_dzl
            sigma_xz = sigma_xz + c15*dux_dxl + c55*(duz_dxl + dux_dzl) + c35*duz_dzl
          else
            ! no attenuation
            sigma_xx = sigma_xx + lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
            sigma_xz = sigma_xz + mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
            sigma_zz = sigma_zz + lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
          endif

          if(SIMULATION_TYPE == 3) then
            b_sigma_xx = b_sigma_xx + lambdaplus2mu_unrelaxed_elastic*b_dux_dxl + lambdal_unrelaxed_elastic*b_duz_dzl
            b_sigma_xz = b_sigma_xz + mul_unrelaxed_elastic*(b_duz_dxl + b_dux_dzl)
            b_sigma_zz = b_sigma_zz + lambdaplus2mu_unrelaxed_elastic*b_duz_dzl + lambdal_unrelaxed_elastic*b_dux_dxl
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
          else if(iedge_poroelastic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_poroelastic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          else if(iedge_poroelastic ==IRIGHT)then
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

          if(SIMULATION_TYPE == 3) then
            b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - weight* &
                (b_sigma_xx*nx + b_sigma_xz*nz)/2.d0

            b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) - weight* &
                (b_sigma_xz*nx + b_sigma_zz*nz)/2.d0
          endif !if(SIMULATION_TYPE == 3) then

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
          ! if this processor core carries the source and the source element is elastic
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
                        - sin(anglesource(i_source))*source_time_function(i_source,it,i_stage)*hlagrange
                      accel_elastic(3,iglob) = accel_elastic(3,iglob) &
                        + cos(anglesource(i_source))*source_time_function(i_source,it,i_stage)*hlagrange
                    enddo
                  enddo
                else    ! SH (membrane) calculation
                  do j = 1,NGLLZ
                    do i = 1,NGLLX
                      iglob = ibool(i,j,ispec_selected_source(i_source))
                      hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                      accel_elastic(2,iglob) = accel_elastic(2,iglob) &
                            + source_time_function(i_source,it,i_stage)*hlagrange
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
                        - sin(anglesource(i_source))*source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1) &
                          *hlagrange
                      b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) &
                        + cos(anglesource(i_source))*source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1) &
                          *hlagrange
                    enddo
                  enddo
                else    ! SH (membrane) calculation
                  do j = 1,NGLLZ
                    do i = 1,NGLLX
                      iglob = ibool(i,j,ispec_selected_source(i_source))
                      hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                      b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) &
                            + source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)*hlagrange
                    enddo
                  enddo
                endif

              endif  !endif SIMULATION_TYPE == 1
            endif

          endif ! if this processor core carries the source and the source element is elastic
        enddo ! do i_source=1,NSOURCES

!<NOISE_TOMOGRAPHY

        ! inject wavefield sources for noise simulations

        if (NOISE_TOMOGRAPHY == 1) then
          call  add_point_source_noise(p_sv,it,NSTEP,nglob,ibool,ispec_noise, &
                            accel_elastic,angle_noise,source_array_noise)

        else if (NOISE_TOMOGRAPHY == 2) then
          call add_surface_movie_noise(p_sv,NOISE_TOMOGRAPHY,it,NSTEP,nspec,nglob,ibool,accel_elastic, &
                            surface_movie_x_noise,surface_movie_y_noise, &
                            surface_movie_z_noise,mask_noise,jacobian,wxgll,wzgll)

        else if (NOISE_TOMOGRAPHY == 3) then
          if (.not. save_everywhere) then
            call add_surface_movie_noise(p_sv,NOISE_TOMOGRAPHY,it,NSTEP,nspec,nglob,ibool,b_accel_elastic, &
                              surface_movie_x_noise,surface_movie_y_noise, &
                              surface_movie_z_noise,mask_noise,jacobian,wxgll,wzgll)
          endif
        endif

!>NOISE_TOMOGRAPHY


      endif ! if not using an initial field
    endif !if(any_elastic)

! assembling accel_elastic for elastic elements
#ifdef USE_MPI

    if(time_stepping_scheme == 2)then
    if(i_stage==1 .and. it == 1 .and. (.not. initialfield))then
    veloc_elastic_LDDRK_temp = veloc_elastic
      if (nproc > 1 .and. any_elastic .and. ninterface_elastic > 0) then
       call assemble_MPI_vector_el(veloc_elastic,nglob, &
            ninterface, ninterface_elastic,inum_interfaces_elastic, &
            max_interface_size, max_ibool_interfaces_size_el,&
            ibool_interfaces_elastic, nibool_interfaces_elastic, &
            tab_requests_send_recv_elastic,buffer_send_faces_vector_el, &
            buffer_recv_faces_vector_el, my_neighbours)
       endif
    endif
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,ier)

    if (nproc > 1 .and. any_elastic .and. ninterface_elastic > 0) then
      call assemble_MPI_vector_el(accel_elastic,nglob, &
            ninterface, ninterface_elastic,inum_interfaces_elastic, &
            max_interface_size, max_ibool_interfaces_size_el,&
            ibool_interfaces_elastic, nibool_interfaces_elastic, &
            tab_requests_send_recv_elastic,buffer_send_faces_vector_el, &
            buffer_recv_faces_vector_el, my_neighbours)
    endif


    if (nproc > 1 .and. any_elastic .and. ninterface_elastic > 0 .and. SIMULATION_TYPE == 3) then
      call assemble_MPI_vector_el(b_accel_elastic,nglob, &
            ninterface, ninterface_elastic,inum_interfaces_elastic, &
            max_interface_size, max_ibool_interfaces_size_el,&
            ibool_interfaces_elastic, nibool_interfaces_elastic, &
            tab_requests_send_recv_elastic,buffer_send_faces_vector_el, &
            buffer_recv_faces_vector_el, my_neighbours)
    endif
#endif

      if(PML_BOUNDARY_CONDITIONS .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1)then
       if(any_elastic .and. nglob_interface > 0)then
        do i = 1, nglob_interface
          write(71)accel_elastic(1,point_interface(i)),accel_elastic(2,point_interface(i)),&
                   accel_elastic(3,point_interface(i)),&
                   veloc_elastic(1,point_interface(i)),veloc_elastic(2,point_interface(i)),&
                   veloc_elastic(3,point_interface(i)),&
                   displ_elastic(1,point_interface(i)),displ_elastic(2,point_interface(i)),&
                   displ_elastic(3,point_interface(i))
        enddo
       endif
      endif

      if(SIMULATION_TYPE == 3)then
        if(PML_BOUNDARY_CONDITIONS)then
          if(any_elastic .and. nglob_interface > 0)then
            do i = 1, nglob_interface
              b_accel_elastic(1,point_interface(i)) = pml_interface_history_accel(1,i,NSTEP-it+1)
              b_accel_elastic(2,point_interface(i)) = pml_interface_history_accel(2,i,NSTEP-it+1)
              b_accel_elastic(3,point_interface(i)) = pml_interface_history_accel(3,i,NSTEP-it+1)
            enddo
          endif
        endif
      endif


! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

    if(any_elastic) then

!! DK DK this should be vectorized
      accel_elastic(1,:) = accel_elastic(1,:) * rmass_inverse_elastic_one
      accel_elastic(2,:) = accel_elastic(2,:) * rmass_inverse_elastic_one
      accel_elastic(3,:) = accel_elastic(3,:) * rmass_inverse_elastic_three

      if(time_stepping_scheme == 1)then
!! DK DK this should be vectorized
        veloc_elastic = veloc_elastic + deltatover2 * accel_elastic
      endif

      if(time_stepping_scheme == 2)then

!! DK DK this should be vectorized
        veloc_elastic_LDDRK = alpha_LDDRK(i_stage) * veloc_elastic_LDDRK + deltat * accel_elastic
        displ_elastic_LDDRK = alpha_LDDRK(i_stage) * displ_elastic_LDDRK + deltat * veloc_elastic
        if(i_stage==1 .and. it == 1 .and. (.not. initialfield))then
        veloc_elastic_LDDRK_temp = veloc_elastic_LDDRK_temp + beta_LDDRK(i_stage) * veloc_elastic_LDDRK
        veloc_elastic = veloc_elastic_LDDRK_temp
        else
        veloc_elastic = veloc_elastic + beta_LDDRK(i_stage) * veloc_elastic_LDDRK
        endif
        displ_elastic = displ_elastic + beta_LDDRK(i_stage) * displ_elastic_LDDRK

      endif

      if(time_stepping_scheme == 3)then

!! DK DK this should be vectorized
        accel_elastic_rk(1,:,i_stage) = deltat * accel_elastic(1,:)
        accel_elastic_rk(2,:,i_stage) = deltat * accel_elastic(2,:)
        accel_elastic_rk(3,:,i_stage) = deltat * accel_elastic(3,:)

        veloc_elastic_rk(1,:,i_stage) = deltat * veloc_elastic(1,:)
        veloc_elastic_rk(2,:,i_stage) = deltat * veloc_elastic(2,:)
        veloc_elastic_rk(3,:,i_stage) = deltat * veloc_elastic(3,:)

        if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

        if(i_stage == 1)weight_rk = 0.5d0
        if(i_stage == 2)weight_rk = 0.5d0
        if(i_stage == 3)weight_rk = 1.0d0

        if(i_stage==1)then

!! DK DK this should be vectorized
        veloc_elastic_initial_rk(1,:) = veloc_elastic(1,:)
        veloc_elastic_initial_rk(2,:) = veloc_elastic(2,:)
        veloc_elastic_initial_rk(3,:) = veloc_elastic(3,:)

        displ_elastic_initial_rk(1,:) = displ_elastic(1,:)
        displ_elastic_initial_rk(2,:) = displ_elastic(2,:)
        displ_elastic_initial_rk(3,:) = displ_elastic(3,:)

        endif

!! DK DK this should be vectorized
        veloc_elastic(1,:) = veloc_elastic_initial_rk(1,:) + weight_rk * accel_elastic_rk(1,:,i_stage)
        veloc_elastic(2,:) = veloc_elastic_initial_rk(2,:) + weight_rk * accel_elastic_rk(2,:,i_stage)
        veloc_elastic(3,:) = veloc_elastic_initial_rk(3,:) + weight_rk * accel_elastic_rk(3,:,i_stage)

        displ_elastic(1,:) = displ_elastic_initial_rk(1,:) + weight_rk * veloc_elastic_rk(1,:,i_stage)
        displ_elastic(2,:) = displ_elastic_initial_rk(2,:) + weight_rk * veloc_elastic_rk(2,:,i_stage)
        displ_elastic(3,:) = displ_elastic_initial_rk(3,:) + weight_rk * veloc_elastic_rk(3,:,i_stage)

        else if(i_stage==4)then

!! DK DK this should be vectorized
        veloc_elastic(1,:) = veloc_elastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
        (accel_elastic_rk(1,:,1) + 2.0d0 * accel_elastic_rk(1,:,2) + &
         2.0d0 * accel_elastic_rk(1,:,3) + accel_elastic_rk(1,:,4))

        veloc_elastic(2,:) = veloc_elastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
        (accel_elastic_rk(2,:,1) + 2.0d0 * accel_elastic_rk(2,:,2) + &
         2.0d0 * accel_elastic_rk(2,:,3) + accel_elastic_rk(2,:,4))

         veloc_elastic(3,:) = veloc_elastic_initial_rk(3,:) + 1.0d0 / 6.0d0 * &
        (accel_elastic_rk(3,:,1) + 2.0d0 * accel_elastic_rk(3,:,2) + &
         2.0d0 * accel_elastic_rk(3,:,3) + accel_elastic_rk(3,:,4))

        displ_elastic(1,:) = displ_elastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
        (veloc_elastic_rk(1,:,1) + 2.0d0 * veloc_elastic_rk(1,:,2) + &
         2.0d0 * veloc_elastic_rk(1,:,3) + veloc_elastic_rk(1,:,4))

        displ_elastic(2,:) = displ_elastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
        (veloc_elastic_rk(2,:,1) + 2.0d0 * veloc_elastic_rk(2,:,2) + &
         2.0d0 * veloc_elastic_rk(2,:,3) + veloc_elastic_rk(2,:,4))

        displ_elastic(3,:) = displ_elastic_initial_rk(3,:) + 1.0d0 / 6.0d0 * &
        (veloc_elastic_rk(3,:,1) + 2.0d0 * veloc_elastic_rk(3,:,2) + &
         2.0d0 * veloc_elastic_rk(3,:,3) + veloc_elastic_rk(3,:,4))

        endif

      endif

      if(SIMULATION_TYPE == 3) then
!! DK DK this should be vectorized
        b_accel_elastic(1,:) = b_accel_elastic(1,:) * rmass_inverse_elastic_one(:)
        b_accel_elastic(2,:) = b_accel_elastic(2,:) * rmass_inverse_elastic_one(:)
        b_accel_elastic(3,:) = b_accel_elastic(3,:) * rmass_inverse_elastic_three(:)

        b_veloc_elastic = b_veloc_elastic + b_deltatover2*b_accel_elastic
      endif

    endif !if(any_elastic)


! ******************************************************************************************************************
! ************* main solver for the poroelastic elements: first the solid (u_s) then the fluid (w)
! ******************************************************************************************************************

    if(any_poroelastic) then

!--------------------------------------------------------------------------------------------
! implement viscous attenuation for poroelastic media
!
    if(ATTENUATION_PORO_FLUID_PART) then
      ! update memory variables with fourth-order Runge-Kutta time scheme for attenuation
      ! loop over spectral elements
      do ispec = 1,nspec

       if(poroelastic(ispec) .and. poroelastcoef(2,2,kmato(ispec)) >0.d0) then

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
        bl_unrelaxed_elastic(1) = etal_f*invpermlxx
        bl_unrelaxed_elastic(2) = etal_f*invpermlxz
        bl_unrelaxed_elastic(3) = etal_f*invpermlzz

        do j=1,NGLLZ
          do i=1,NGLLX

            iglob = ibool(i,j,ispec)
            viscox_loc(i,j) = velocw_poroelastic(1,iglob)*bl_unrelaxed_elastic(1) + &
                               velocw_poroelastic(2,iglob) * bl_unrelaxed_elastic(2)
            viscoz_loc(i,j) = velocw_poroelastic(1,iglob)*bl_unrelaxed_elastic(2) + &
                               velocw_poroelastic(2,iglob)*bl_unrelaxed_elastic(3)

            if(time_stepping_scheme == 1) then
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
            endif

            if(time_stepping_scheme == 2) then
              Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscox(i,j,ispec)
              rx_viscous_LDDRK(i,j,ispec) = alpha_LDDRK(i_stage) * rx_viscous_LDDRK(i,j,ispec) + &
                                            deltat * (Sn + thetainv * rx_viscous(i,j,ispec))
              rx_viscous(i,j,ispec)= rx_viscous(i,j,ispec)+beta_LDDRK(i_stage) * rx_viscous_LDDRK(i,j,ispec)

              Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscoz(i,j,ispec)
              rz_viscous_LDDRK(i,j,ispec)= alpha_LDDRK(i_stage) * rz_viscous_LDDRK(i,j,ispec)+&
                                           deltat * (Sn + thetainv * rz_viscous(i,j,ispec))
              rz_viscous(i,j,ispec)= rz_viscous(i,j,ispec)+beta_LDDRK(i_stage) * rz_viscous_LDDRK(i,j,ispec)
            endif

            if(time_stepping_scheme == 3) then

              Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscox(i,j,ispec)
              rx_viscous_force_RK(i,j,ispec,i_stage) = deltat * (Sn + thetainv * rx_viscous(i,j,ispec))

              if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then
                if(i_stage == 1)weight_rk = 0.5d0
                if(i_stage == 2)weight_rk = 0.5d0
                if(i_stage == 3)weight_rk = 1.0d0

                if(i_stage==1)then
                  rx_viscous_initial_rk(i,j,ispec) = rx_viscous(i,j,ispec)
                endif
                  rx_viscous(i,j,ispec) = rx_viscous_initial_rk(i,j,ispec) + &
                                          weight_rk * rx_viscous_force_RK(i,j,ispec,i_stage)
              else if(i_stage==4)then

                rx_viscous(i,j,ispec) = rx_viscous_initial_rk(i,j,ispec) + &
                                        1.0d0 / 6.0d0 * (rx_viscous_force_RK(i,j,ispec,i_stage) + &
                                        2.0d0 * rx_viscous_force_RK(i,j,ispec,i_stage) + &
                                        2.0d0 * rx_viscous_force_RK(i,j,ispec,i_stage) + &
                                        rx_viscous_force_RK(i,j,ispec,i_stage))
              endif

              Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscoz(i,j,ispec)
              rz_viscous_force_RK(i,j,ispec,i_stage) = deltat * (Sn + thetainv * rz_viscous(i,j,ispec))

              if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then
                if(i_stage == 1)weight_rk = 0.5d0
                if(i_stage == 2)weight_rk = 0.5d0
                if(i_stage == 3)weight_rk = 1.0d0
                if(i_stage==1)then
                  rz_viscous_initial_rk(i,j,ispec) = rz_viscous(i,j,ispec)
                endif
                rz_viscous(i,j,ispec) = rz_viscous_initial_rk(i,j,ispec) + &
                                        weight_rk * rz_viscous_force_RK(i,j,ispec,i_stage)
              else if(i_stage==4)then
                rz_viscous(i,j,ispec) = rz_viscous_initial_rk(i,j,ispec) + &
                                        1.0d0 / 6.0d0 * (rz_viscous_force_RK(i,j,ispec,i_stage) + &
                                        2.0d0 * rz_viscous_force_RK(i,j,ispec,i_stage) + &
                                        2.0d0 * rz_viscous_force_RK(i,j,ispec,i_stage) + &
                                        rz_viscous_force_RK(i,j,ispec,i_stage))
              endif
            endif
          enddo
        enddo

        if(stage_time_scheme == 1) then
        ! save visco for Runge-Kutta scheme when used together with Newmark
        viscox(:,:,ispec) = viscox_loc(:,:)
        viscoz(:,:,ispec) = viscoz_loc(:,:)
        endif

       endif  ! end of poroelastic element loop

      enddo   ! end of spectral element loop
    endif ! end of viscous attenuation for porous media

!-----------------------------------------

      if(SIMULATION_TYPE == 3) then
        ! if inviscid fluid, comment the reading and uncomment the zeroing
        !     read(23,rec=NSTEP-it+1) b_viscodampx
        !     read(24,rec=NSTEP-it+1) b_viscodampz
        b_viscodampx(:) = ZERO
        b_viscodampz(:) = ZERO
      endif

      call compute_forces_poro_solid(nglob,nspec,myrank,nelemabs,numat, &
               ispec_selected_source,ispec_selected_rec,is_proc_source,which_proc_receiver,&
               source_type,it,NSTEP,anyabs, &
               initialfield,ATTENUATION_VISCOELASTIC_SOLID,ATTENUATION_PORO_FLUID_PART,deltat, &
               ibool,kmato,numabs,poroelastic,codeabs, &
               accels_poroelastic,velocs_poroelastic,velocw_poroelastic,displs_poroelastic,displs_poroelastic_old,&
               displw_poroelastic,b_accels_poroelastic,b_displs_poroelastic,b_displw_poroelastic,&
               density,porosity,tortuosity,permeability,poroelastcoef,xix,xiz,gammax,gammaz, &
               jacobian,source_time_function,sourcearray,adj_sourcearrays, &
               e11,e13,hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll,&
               inv_tau_sigma_nu2,phi_nu2,Mu_nu2,N_SLS, &
               rx_viscous,rz_viscous,theta_e,theta_s,b_viscodampx,b_viscodampz,&
               ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro,iend_edge3_poro, &
               ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro,&
               mufr_k,B_k,NSOURCES,nrec,SIMULATION_TYPE,SAVE_FORWARD,&
               b_absorb_poro_s_left,b_absorb_poro_s_right,b_absorb_poro_s_bottom,b_absorb_poro_s_top,&
               nspec_left,nspec_right,nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top,f0(1),freq0,Q0,&
               e11_LDDRK,e13_LDDRK,alpha_LDDRK,beta_LDDRK, &
               e11_initial_rk,e13_initial_rk,e11_force_RK, e13_force_RK, &
               stage_time_scheme,i_stage)

      call compute_forces_poro_fluid(nglob,nspec,myrank,nelemabs,numat, &
               ispec_selected_source,ispec_selected_rec,is_proc_source,which_proc_receiver,&
               source_type,it,NSTEP,anyabs, &
               initialfield,ATTENUATION_VISCOELASTIC_SOLID,ATTENUATION_PORO_FLUID_PART,deltat, &
               ibool,kmato,numabs,poroelastic,codeabs, &
               accelw_poroelastic,velocw_poroelastic,displw_poroelastic,velocs_poroelastic,displs_poroelastic,&
               displs_poroelastic_old,b_accelw_poroelastic,b_displw_poroelastic,b_displs_poroelastic,&
               density,porosity,tortuosity,permeability,poroelastcoef,xix,xiz,gammax,gammaz, &
               jacobian,source_time_function,sourcearray,adj_sourcearrays, &
               e11,e13,hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll,&
               inv_tau_sigma_nu2,phi_nu2,Mu_nu2,N_SLS, &
               rx_viscous,rz_viscous,theta_e,theta_s,b_viscodampx,b_viscodampz,&
               ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro,iend_edge3_poro, &
               ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro,&
               C_k,M_k,NSOURCES,nrec,SIMULATION_TYPE,SAVE_FORWARD,&
               b_absorb_poro_w_left,b_absorb_poro_w_right,b_absorb_poro_w_bottom,b_absorb_poro_w_top,&
               nspec_left,nspec_right,nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top,f0(1),freq0,Q0,&
               e11_LDDRK,e13_LDDRK,alpha_LDDRK,beta_LDDRK, &
               e11_initial_rk,e13_initial_rk,e11_force_RK, e13_force_RK, &
               stage_time_scheme,i_stage)

      if(SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
        ! if inviscid fluid, comment
        !     write(23,rec=it) b_viscodampx
        !     write(24,rec=it) b_viscodampz
      endif

      if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1) then

        !--- left absorbing boundary
        if(nspec_left >0) then
          do ispec = 1,nspec_left
            do id =1,2
              do i=1,NGLLZ
                write(45) b_absorb_poro_s_left(id,i,ispec,it)
                write(25) b_absorb_poro_w_left(id,i,ispec,it)
              enddo
            enddo
          enddo
        endif

        !--- right absorbing boundary
        if(nspec_right >0) then
          do ispec = 1,nspec_right
            do id =1,2
              do i=1,NGLLZ
                write(46) b_absorb_poro_s_right(id,i,ispec,it)
                write(26) b_absorb_poro_w_right(id,i,ispec,it)
              enddo
            enddo
          enddo
        endif

        !--- bottom absorbing boundary
        if(nspec_bottom >0) then
          do ispec = 1,nspec_bottom
            do id =1,2
              do i=1,NGLLX
                write(47) b_absorb_poro_s_bottom(id,i,ispec,it)
                write(29) b_absorb_poro_w_bottom(id,i,ispec,it)
              enddo
            enddo
          enddo
        endif

        !--- top absorbing boundary
        if(nspec_top >0) then
          do ispec = 1,nspec_top
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
          if(SIMULATION_TYPE == 3) then
            b_pressure = - b_potential_dot_dot_acoustic(iglob)
            ! new definition of adjoint displacement and adjoint potential
            pressure = potential_acoustic_adj_coupling(iglob)
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
          else if(iedge_acoustic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_acoustic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          else if(iedge_acoustic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          endif

          ! contribution to the solid phase
          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) &
            + weight*nx*pressure*(1._CUSTOM_REAL-phil/tortl)
          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) &
            + weight*nz*pressure*(1._CUSTOM_REAL-phil/tortl)

          ! contribution to the fluid phase
          accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) &
            + weight*nx*pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
          accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) &
            + weight*nz*pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)

          if(SIMULATION_TYPE == 3) then
            ! contribution to the solid phase
            b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) &
              + weight*nx*b_pressure*(1._CUSTOM_REAL-phil/tortl)
            b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) &
              + weight*nz*b_pressure*(1._CUSTOM_REAL-phil/tortl)

            ! contribution to the fluid phase
            b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) &
              + weight*nx*b_pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
            b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) &
              + weight*nz*b_pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
          endif !if(SIMULATION_TYPE == 3) then

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
          lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec_elastic))
          mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec_elastic))
          lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec_elastic))

          ! derivative along x and along z for u_s and w
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          if(SIMULATION_TYPE == 3) then
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

            if(SIMULATION_TYPE == 3) then
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

          if(SIMULATION_TYPE == 3) then
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
              c12 = c12ext(i,j,ispec_elastic)
              c23 = c23ext(i,j,ispec_elastic)
              c25 = c25ext(i,j,ispec_elastic)
            else
              c11 = anisotropy(1,kmato(ispec_elastic))
              c13 = anisotropy(2,kmato(ispec_elastic))
              c15 = anisotropy(3,kmato(ispec_elastic))
              c33 = anisotropy(4,kmato(ispec_elastic))
              c35 = anisotropy(5,kmato(ispec_elastic))
              c55 = anisotropy(6,kmato(ispec_elastic))
              c12 = anisotropy(7,kmato(ispec_elastic))
              c23 = anisotropy(8,kmato(ispec_elastic))
              c25 = anisotropy(9,kmato(ispec_elastic))
            endif
            sigma_xx = c11*dux_dxl + c15*(duz_dxl + dux_dzl) + c13*duz_dzl
            sigma_zz = c13*dux_dxl + c35*(duz_dxl + dux_dzl) + c33*duz_dzl
            sigma_xz = c15*dux_dxl + c55*(duz_dxl + dux_dzl) + c35*duz_dzl
          else
            ! no attenuation
            sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
            sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
            sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
          endif

          if(SIMULATION_TYPE == 3) then
            b_sigma_xx = lambdaplus2mu_unrelaxed_elastic*b_dux_dxl + lambdal_unrelaxed_elastic*b_duz_dzl
            b_sigma_xz = mul_unrelaxed_elastic*(b_duz_dxl + b_dux_dzl)
            b_sigma_zz = lambdaplus2mu_unrelaxed_elastic*b_duz_dzl + lambdal_unrelaxed_elastic*b_dux_dxl
          endif ! if(SIMULATION_TYPE == 3)

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

          if(SIMULATION_TYPE == 3) then
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
            if(SIMULATION_TYPE == 3) then
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

          if(SIMULATION_TYPE == 3) then
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

          if(SIMULATION_TYPE == 3) then
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
          else if(iedge_poroelastic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_poroelastic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          else if(iedge_poroelastic ==IRIGHT)then
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

          if(SIMULATION_TYPE == 3) then
            ! contribution to the solid phase
            b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) + &
                weight*((b_sigma_xx*nx + b_sigma_xz*nz)/2.d0 -phil/tortl*b_sigmap*nx)

            b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + &
                weight*((b_sigma_xz*nx + b_sigma_zz*nz)/2.d0 -phil/tortl*b_sigmap*nz)

            ! contribution to the fluid phase
            ! w = 0
          endif !if(SIMULATION_TYPE == 3) then

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
          ! if this processor core carries the source and the source element is elastic
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
                      (1._CUSTOM_REAL - phil/tortl)*sin(anglesource(i_source))*source_time_function(i_source,it,i_stage)
                    accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + hlagrange * &
                      (1._CUSTOM_REAL - phil/tortl)*cos(anglesource(i_source))*source_time_function(i_source,it,i_stage)
                    ! w
                    accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - hlagrange * &
                      (1._CUSTOM_REAL - rhol_f/rhol_bar)*sin(anglesource(i_source))*source_time_function(i_source,it,i_stage)
                    accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) + hlagrange * &
                      (1._CUSTOM_REAL - rhol_f/rhol_bar)*cos(anglesource(i_source))*source_time_function(i_source,it,i_stage)
                  enddo
                enddo
              else                   ! backward wavefield
                do j = 1,NGLLZ
                  do i = 1,NGLLX
                    iglob = ibool(i,j,ispec_selected_source(i_source))
                    hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
                    ! b_s
                    b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - hlagrange * &
                      (1._CUSTOM_REAL - phil/tortl)*sin(anglesource(i_source))* &
                      source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)
                    b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + hlagrange * &
                      (1._CUSTOM_REAL - phil/tortl)*cos(anglesource(i_source))* &
                      source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)
                    !b_w
                    b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) - hlagrange * &
                      (1._CUSTOM_REAL - rhol_f/rhol_bar)*sin(anglesource(i_source))* &
                      source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)
                    b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) + hlagrange * &
                      (1._CUSTOM_REAL - rhol_f/rhol_bar)*cos(anglesource(i_source))* &
                      source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)
                  enddo
                enddo
              endif !endif SIMULATION_TYPE == 1
            endif

          endif ! if this processor core carries the source and the source element is elastic
        enddo ! do i_source=1,NSOURCES

      endif ! if not using an initial field
    endif !if(any_poroelastic)

! assembling accels_proelastic & accelw_poroelastic for poroelastic elements
#ifdef USE_MPI
    if (nproc > 1 .and. any_poroelastic .and. ninterface_poroelastic > 0) then
      call assemble_MPI_vector_po(accels_poroelastic,accelw_poroelastic,nglob, &
            ninterface, ninterface_poroelastic,inum_interfaces_poroelastic, &
            max_interface_size, max_ibool_interfaces_size_po,&
            ibool_interfaces_poroelastic, nibool_interfaces_poroelastic, &
            tab_requests_send_recv_poro,buffer_send_faces_vector_pos,buffer_send_faces_vector_pow, &
            buffer_recv_faces_vector_pos,buffer_recv_faces_vector_pow, &
            my_neighbours)
    endif

    if (nproc > 1 .and. any_poroelastic .and. ninterface_poroelastic > 0 .and. SIMULATION_TYPE == 3) then
      call assemble_MPI_vector_po(b_accels_poroelastic,b_accelw_poroelastic,nglob, &
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

      if(time_stepping_scheme == 1)then

      accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
      accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)
      velocs_poroelastic = velocs_poroelastic + deltatover2*accels_poroelastic

      accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
      accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)
      velocw_poroelastic = velocw_poroelastic + deltatover2*accelw_poroelastic

      endif

      if(time_stepping_scheme == 2)then

        accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
        accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)

  velocs_poroelastic_LDDRK = alpha_LDDRK(i_stage) * velocs_poroelastic_LDDRK + deltat * accels_poroelastic
  displs_poroelastic_LDDRK = alpha_LDDRK(i_stage) * displs_poroelastic_LDDRK + deltat * velocs_poroelastic

  velocs_poroelastic = velocs_poroelastic + beta_LDDRK(i_stage) * velocs_poroelastic_LDDRK
  displs_poroelastic = displs_poroelastic + beta_LDDRK(i_stage) * displs_poroelastic_LDDRK

        accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
        accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)

  velocw_poroelastic_LDDRK = alpha_LDDRK(i_stage) * velocw_poroelastic_LDDRK + deltat * accelw_poroelastic
  displw_poroelastic_LDDRK = alpha_LDDRK(i_stage) * displw_poroelastic_LDDRK + deltat * velocw_poroelastic

  velocw_poroelastic = velocw_poroelastic + beta_LDDRK(i_stage) * velocw_poroelastic_LDDRK
  displw_poroelastic = displw_poroelastic + beta_LDDRK(i_stage) * displw_poroelastic_LDDRK

      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(time_stepping_scheme == 3)then

        accels_poroelastic(1,:) = accels_poroelastic(1,:) * rmass_s_inverse_poroelastic(:)
        accels_poroelastic(2,:) = accels_poroelastic(2,:) * rmass_s_inverse_poroelastic(:)

        accels_poroelastic_rk(1,:,i_stage) = deltat * accels_poroelastic(1,:)
        accels_poroelastic_rk(2,:,i_stage) = deltat * accels_poroelastic(2,:)
        velocs_poroelastic_rk(1,:,i_stage) = deltat * velocs_poroelastic(1,:)
        velocs_poroelastic_rk(2,:,i_stage) = deltat * velocs_poroelastic(2,:)

        accelw_poroelastic(1,:) = accelw_poroelastic(1,:) * rmass_w_inverse_poroelastic(:)
        accelw_poroelastic(2,:) = accelw_poroelastic(2,:) * rmass_w_inverse_poroelastic(:)

        accelw_poroelastic_rk(1,:,i_stage) = deltat * accelw_poroelastic(1,:)
        accelw_poroelastic_rk(2,:,i_stage) = deltat * accelw_poroelastic(2,:)
        velocw_poroelastic_rk(1,:,i_stage) = deltat * velocw_poroelastic(1,:)
        velocw_poroelastic_rk(2,:,i_stage) = deltat * velocw_poroelastic(2,:)

        if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

        if(i_stage == 1)weight_rk = 0.5d0
        if(i_stage == 2)weight_rk = 0.5d0
        if(i_stage == 3)weight_rk = 1.0d0

        if(i_stage==1)then

        velocs_poroelastic_initial_rk(1,:) = velocs_poroelastic(1,:)
        velocs_poroelastic_initial_rk(2,:) = velocs_poroelastic(2,:)
        displs_poroelastic_initial_rk(1,:) = displs_poroelastic(1,:)
        displs_poroelastic_initial_rk(2,:) = displs_poroelastic(2,:)

        velocw_poroelastic_initial_rk(1,:) = velocw_poroelastic(1,:)
        velocw_poroelastic_initial_rk(2,:) = velocw_poroelastic(2,:)
        displw_poroelastic_initial_rk(1,:) = displw_poroelastic(1,:)
        displw_poroelastic_initial_rk(2,:) = displw_poroelastic(2,:)

        endif

        velocs_poroelastic(1,:) = velocs_poroelastic_initial_rk(1,:) + weight_rk * accels_poroelastic_rk(1,:,i_stage)
  velocs_poroelastic(2,:) = velocs_poroelastic_initial_rk(2,:) + weight_rk * accels_poroelastic_rk(2,:,i_stage)
        displs_poroelastic(1,:) = displs_poroelastic_initial_rk(1,:) + weight_rk * velocs_poroelastic_rk(1,:,i_stage)
  displs_poroelastic(2,:) = displs_poroelastic_initial_rk(2,:) + weight_rk * velocs_poroelastic_rk(2,:,i_stage)

        velocw_poroelastic(1,:) = velocw_poroelastic_initial_rk(1,:) + weight_rk * accelw_poroelastic_rk(1,:,i_stage)
  velocw_poroelastic(2,:) = velocw_poroelastic_initial_rk(2,:) + weight_rk * accelw_poroelastic_rk(2,:,i_stage)
        displw_poroelastic(1,:) = displw_poroelastic_initial_rk(1,:) + weight_rk * velocw_poroelastic_rk(1,:,i_stage)
  displw_poroelastic(2,:) = displw_poroelastic_initial_rk(2,:) + weight_rk * velocw_poroelastic_rk(2,:,i_stage)


        else if(i_stage==4)then

        velocs_poroelastic(1,:) = velocs_poroelastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
        (accels_poroelastic_rk(1,:,1) + 2.0d0 * accels_poroelastic_rk(1,:,2) + &
         2.0d0 * accels_poroelastic_rk(1,:,3) + accels_poroelastic_rk(1,:,4))

        velocs_poroelastic(2,:) = velocs_poroelastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
        (accels_poroelastic_rk(2,:,1) + 2.0d0 * accels_poroelastic_rk(2,:,2) + &
         2.0d0 * accels_poroelastic_rk(2,:,3) + accels_poroelastic_rk(2,:,4))

        displs_poroelastic(1,:) = displs_poroelastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
        (velocs_poroelastic_rk(1,:,1) + 2.0d0 * velocs_poroelastic_rk(1,:,2) + &
         2.0d0 * velocs_poroelastic_rk(1,:,3) + velocs_poroelastic_rk(1,:,4))

        displs_poroelastic(2,:) = displs_poroelastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
        (velocs_poroelastic_rk(2,:,1) + 2.0d0 * velocs_poroelastic_rk(2,:,2) + &
         2.0d0 * velocs_poroelastic_rk(2,:,3) + velocs_poroelastic_rk(2,:,4))

        velocw_poroelastic(1,:) = velocw_poroelastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
        (accelw_poroelastic_rk(1,:,1) + 2.0d0 * accelw_poroelastic_rk(1,:,2) + &
         2.0d0 * accelw_poroelastic_rk(1,:,3) + accelw_poroelastic_rk(1,:,4))

        velocw_poroelastic(2,:) = velocw_poroelastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
        (accelw_poroelastic_rk(2,:,1) + 2.0d0 * accelw_poroelastic_rk(2,:,2) + &
         2.0d0 * accelw_poroelastic_rk(2,:,3) + accelw_poroelastic_rk(2,:,4))

        displw_poroelastic(1,:) = displw_poroelastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
        (velocw_poroelastic_rk(1,:,1) + 2.0d0 * velocw_poroelastic_rk(1,:,2) + &
         2.0d0 * velocw_poroelastic_rk(1,:,3) + velocw_poroelastic_rk(1,:,4))

        displw_poroelastic(2,:) = displw_poroelastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
        (velocw_poroelastic_rk(2,:,1) + 2.0d0 * velocw_poroelastic_rk(2,:,2) + &
         2.0d0 * velocw_poroelastic_rk(2,:,3) + velocw_poroelastic_rk(2,:,4))

        endif

      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(SIMULATION_TYPE == 3) then
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

!
! Explanation of the code below, from Christina Morency and Yang Luo, January 2012:
!
! Coupled elastic-poroelastic simulations imply continuity of traction and
! displacement at the interface.
! For the traction we pass on both sides n*(T + Te)/2 , that is, the average
! between the total stress (from the poroelastic part) and the elastic stress.
! For the displacement, we enforce its continuity in the assembling stage,
! realizing that continuity of displacement correspond to the continuity of
! the acceleration we have:
!
! accel_elastic = rmass_inverse_elastic * force_elastic
! accels_poroelastic = rmass_s_inverse_poroelastic * force_poroelastic
!
! Therefore, continuity of acceleration gives
!
! accel = (force_elastic + force_poroelastic)/
!     (1/rmass_inverse_elastic + 1/rmass_inverse_poroelastic)
!
! Then
!
! accel_elastic = accel
! accels_poroelastic = accel
! accelw_poroelastic = 0
!
! From there, the velocity and displacement are updated.
! Note that force_elastic and force_poroelastic are the right hand sides of
! the equations we solve, that is, the acceleration terms before the
! division by the inverse of the mass matrices. This is why in the code below
! we first need to recover the accelerations (which are then
! the right hand sides forces) and the velocities before the update.
!
! This implementation highly helped stability especially with unstructured meshes.
!

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

            if(time_stepping_scheme == 1)then

            veloc_elastic(1,iglob) = veloc_elastic(1,iglob) - deltatover2*accel_elastic(1,iglob)
            veloc_elastic(3,iglob) = veloc_elastic(3,iglob) - deltatover2*accel_elastic(3,iglob)
            accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic_one(iglob)
            accel_elastic(3,iglob) = accel_elastic(3,iglob) / rmass_inverse_elastic_three(iglob)
            ! recovering original velocities and accelerations on boundaries (poro side)
            velocs_poroelastic(1,iglob) = velocs_poroelastic(1,iglob) - deltatover2*accels_poroelastic(1,iglob)
            velocs_poroelastic(2,iglob) = velocs_poroelastic(2,iglob) - deltatover2*accels_poroelastic(2,iglob)
            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)
            ! assembling accelerations
            accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
                                   ( 1.0/rmass_inverse_elastic_one(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
            accel_elastic(3,iglob) = ( accel_elastic(3,iglob) + accels_poroelastic(2,iglob) ) / &
                                   ( 1.0/rmass_inverse_elastic_three(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
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

            endif

!            if(time_stepping_scheme == 2)then
            ! recovering original velocities and accelerations on boundaries (elastic side)
!      veloc_elastic = veloc_elastic - beta_LDDRK(i_stage) * veloc_elastic_LDDRK
!      displ_elastic = displ_elastic - beta_LDDRK(i_stage) * displ_elastic_LDDRK
!      veloc_elastic_LDDRK = (veloc_elastic_LDDRK - deltat * accel_elastic) / alpha_LDDRK(i_stage)
!      displ_elastic_LDDRK = (displ_elastic_LDDRK - deltat * veloc_elastic) / alpha_LDDRK(i_stage)
!            accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic(iglob)
!            accel_elastic(3,iglob) = accel_elastic(3,iglob) / rmass_inverse_elastic(iglob)

            ! recovering original velocities and accelerations on boundaries (poro side)
!      velocs_poroelastic = velocs_poroelastic - beta_LDDRK(i_stage) * velocs_poroelastic_LDDRK
!      displs_poroelastic = displs_poroelastic - beta_LDDRK(i_stage) * displs_poroelastic_LDDRK
!      velocs_poroelastic_LDDRK = (velocs_poroelastic_LDDRK - deltat * accels_poroelastic) / alpha_LDDRK(i_stage)
!      displs_poroelastic_LDDRK = (velocs_poroelastic_LDDRK - deltat * velocs_poroelastic) / alpha_LDDRK(i_stage)
!            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
!            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)

            ! assembling accelerations
!            accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accel_elastic(3,iglob) = ( accel_elastic(3,iglob) + accels_poroelastic(2,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accels_poroelastic(1,iglob) = accel_elastic(1,iglob)
!            accels_poroelastic(2,iglob) = accel_elastic(3,iglob)

      ! updating velocities
            ! updating velocities(elastic side)
!      veloc_elastic_LDDRK = alpha_LDDRK(i_stage) * veloc_elastic_LDDRK + deltat * accel_elastic
!      displ_elastic_LDDRK = alpha_LDDRK(i_stage) * displ_elastic_LDDRK + deltat * veloc_elastic
!      veloc_elastic = veloc_elastic + beta_LDDRK(i_stage) * veloc_elastic_LDDRK
!      displ_elastic = displ_elastic + beta_LDDRK(i_stage) * displ_elastic_LDDRK
            ! updating velocities(poro side)
!      velocs_poroelastic_LDDRK = alpha_LDDRK(i_stage) * velocs_poroelastic_LDDRK + deltat * accels_poroelastic
!      displs_poroelastic_LDDRK = alpha_LDDRK(i_stage) * displs_poroelastic_LDDRK + deltat * velocs_poroelastic
!      velocs_poroelastic = velocs_poroelastic + beta_LDDRK(i_stage) * velocs_poroelastic_LDDRK
!      displs_poroelastic = displs_poroelastic + beta_LDDRK(i_stage) * displs_poroelastic_LDDRK

            ! zeros w
!            accelw_poroelastic(1,iglob) = ZERO
!            accelw_poroelastic(2,iglob) = ZERO
!            velocw_poroelastic(1,iglob) = ZERO
!            velocw_poroelastic(2,iglob) = ZERO
!            endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if(time_stepping_scheme == 3)then

        ! recovering original velocities and accelerations on boundaries (elastic side)
!        if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

!        if(i_stage == 1)weight_rk = 0.5d0
!        if(i_stage == 2)weight_rk = 0.5d0
!        if(i_stage == 3)weight_rk = 1.0d0

!        veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) - weight_rk * accel_elastic_rk(1,iglob,i_stage)
!  veloc_elastic(3,iglob) = veloc_elastic_initial_rk(3,iglob) - weight_rk * accel_elastic_rk(3,iglob,i_stage)
!        displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) - weight_rk * veloc_elastic_rk(1,iglob,i_stage)
!  displ_elastic(3,iglob) = displ_elastic_initial_rk(3,iglob) - weight_rk * veloc_elastic_rk(3,iglob,i_stage)


!        else if(i_stage==4)then

!        veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (accel_elastic_rk(1,iglob,1) + 2.0d0 * accel_elastic_rk(1,iglob,2) + &
!         2.0d0 * accel_elastic_rk(1,iglob,3) + accel_elastic_rk(1,iglob,4))

!        veloc_elastic(3,iglob) = veloc_elastic_initial_rk(3,iglob) - 1.0d0 / 6.0d0 * &
!        (accel_elastic_rk(3,iglob,1) + 2.0d0 * accel_elastic_rk(3,iglob,2) + &
!         2.0d0 * accel_elastic_rk(3,iglob,3) + accel_elastic_rk(3,iglob,4))

!        displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (veloc_elastic_rk(1,iglob,1) + 2.0d0 * veloc_elastic_rk(1,iglob,2) + &
!         2.0d0 * veloc_elastic_rk(1,iglob,3) + veloc_elastic_rk(1,iglob,4))

!        displ_elastic(3,iglob) = displ_elastic_initial_rk(3,iglob) - 1.0d0 / 6.0d0 * &
!        (veloc_elastic_rk(3,iglob,1) + 2.0d0 * veloc_elastic_rk(3,iglob,2) + &
!         2.0d0 * veloc_elastic_rk(3,iglob,3) + veloc_elastic_rk(3,iglob,4))

!        endif

!        accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic(iglob)
!        accel_elastic(3,iglob) = accel_elastic(3,iglob) / rmass_inverse_elastic(iglob)

!        accel_elastic_rk(1,iglob,i_stage) = accel_elastic(1,iglob) / deltat
!        accel_elastic_rk(3,iglob,i_stage) = accel_elastic(3,iglob) / deltat
!        veloc_elastic_rk(1,iglob,i_stage) = veloc_elastic(1,iglob) / deltat
!        veloc_elastic_rk(3,iglob,i_stage) = veloc_elastic(3,iglob) / deltat


        ! recovering original velocities and accelerations on boundaries (poro side)
!        if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

!        if(i_stage == 1)weight_rk = 0.5d0
!        if(i_stage == 2)weight_rk = 0.5d0
!        if(i_stage == 3)weight_rk = 1.0d0

!        velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) - weight_rk * accels_poroelastic_rk(1,iglob,i_stage)
!  velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) - weight_rk * accels_poroelastic_rk(2,iglob,i_stage)
!        displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) - weight_rk * velocs_poroelastic_rk(1,iglob,i_stage)
!  displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) - weight_rk * velocs_poroelastic_rk(2,iglob,i_stage)


!        else if(i_stage==4)then

!        velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (accels_poroelastic_rk(1,iglob,1) + 2.0d0 * accels_poroelastic_rk(1,iglob,2) + &
!         2.0d0 * accels_poroelastic_rk(1,iglob,3) + accels_poroelastic_rk(1,iglob,4))

!        velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) - 1.0d0 / 6.0d0 * &
!        (accels_poroelastic_rk(2,iglob,1) + 2.0d0 * accels_poroelastic_rk(2,iglob,2) + &
!         2.0d0 * accels_poroelastic_rk(2,iglob,3) + accels_poroelastic_rk(2,iglob,4))

!        displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (velocs_poroelastic_rk(1,iglob,1) + 2.0d0 * velocs_poroelastic_rk(1,iglob,2) + &
!         2.0d0 * velocs_poroelastic_rk(1,iglob,3) + velocs_poroelastic_rk(1,iglob,4))

!        displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) - 1.0d0 / 6.0d0 * &
!        (velocs_poroelastic_rk(2,iglob,1) + 2.0d0 * velocs_poroelastic_rk(2,iglob,2) + &
!         2.0d0 * velocs_poroelastic_rk(2,iglob,3) + velocs_poroelastic_rk(2,iglob,4))

!        endif

!        accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
!        accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)

!        accels_poroelastic_rk(1,iglob,i_stage) = accels_poroelastic(1,iglob) / deltat
!        accels_poroelastic_rk(2,iglob,i_stage) = accels_poroelastic(2,iglob) / deltat
!        velocs_poroelastic_rk(1,iglob,i_stage) = velocs_poroelastic(1,iglob) / deltat
!        velocs_poroelastic_rk(2,iglob,i_stage) = velocs_poroelastic(2,iglob) / deltat


        ! assembling accelerations
!            accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accel_elastic(3,iglob) = ( accel_elastic(3,iglob) + accels_poroelastic(2,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accels_poroelastic(1,iglob) = accel_elastic(1,iglob)
!            accels_poroelastic(2,iglob) = accel_elastic(3,iglob)

   ! updating velocities
        ! updating velocities(elastic side)

 !       accel_elastic_rk(1,iglob,i_stage) = accel_elastic(1,iglob) * deltat
 !       accel_elastic_rk(3,iglob,i_stage) = accel_elastic(3,iglob) * deltat

 !       if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

 !       if(i_stage == 1)weight_rk = 0.5d0
 !       if(i_stage == 2)weight_rk = 0.5d0
 !       if(i_stage == 3)weight_rk = 1.0d0

 !       veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) + weight_rk * accel_elastic_rk(1,iglob,i_stage)
 ! veloc_elastic(3,iglob) = veloc_elastic_initial_rk(3,iglob) + weight_rk * accel_elastic_rk(3,iglob,i_stage)
 !       displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) + weight_rk * veloc_elastic_rk(1,iglob,i_stage)
 ! displ_elastic(3,iglob) = displ_elastic_initial_rk(3,iglob) + weight_rk * veloc_elastic_rk(3,iglob,i_stage)


 !       else if(i_stage==4)then

 !       veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (accel_elastic_rk(1,iglob,1) + 2.0d0 * accel_elastic_rk(1,iglob,2) + &
 !        2.0d0 * accel_elastic_rk(1,iglob,3) + accel_elastic_rk(1,iglob,4))
!
 !       veloc_elastic(3,iglob) = veloc_elastic_initial_rk(3,iglob) + 1.0d0 / 6.0d0 * &
 !       (accel_elastic_rk(3,iglob,1) + 2.0d0 * accel_elastic_rk(3,iglob,2) + &
 !        2.0d0 * accel_elastic_rk(3,iglob,3) + accel_elastic_rk(3,iglob,4))

 !       displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (veloc_elastic_rk(1,iglob,1) + 2.0d0 * veloc_elastic_rk(1,iglob,2) + &
 !        2.0d0 * veloc_elastic_rk(1,iglob,3) + veloc_elastic_rk(1,iglob,4))

 !       displ_elastic(3,iglob) = displ_elastic_initial_rk(3,iglob) + 1.0d0 / 6.0d0 * &
 !       (veloc_elastic_rk(3,iglob,1) + 2.0d0 * veloc_elastic_rk(3,iglob,2) + &
 !        2.0d0 * veloc_elastic_rk(3,iglob,3) + veloc_elastic_rk(3,iglob,4))

 !       endif
        ! updating velocities(poro side)

 !       accels_poroelastic_rk(1,iglob,i_stage) = deltat * accels_poroelastic(1,iglob)
 !       accels_poroelastic_rk(2,iglob,i_stage) = deltat * accels_poroelastic(2,iglob)
 !       velocs_poroelastic_rk(1,iglob,i_stage) = deltat * velocs_poroelastic(1,iglob)
 !       velocs_poroelastic_rk(2,iglob,i_stage) = deltat * velocs_poroelastic(2,iglob)


 !       if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then

 !       if(i_stage == 1)weight_rk = 0.5d0
 !       if(i_stage == 2)weight_rk = 0.5d0
 !       if(i_stage == 3)weight_rk = 1.0d0

 !       velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) + weight_rk * accels_poroelastic_rk(1,iglob,i_stage)
 ! velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) + weight_rk * accels_poroelastic_rk(2,iglob,i_stage)
 !       displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) + weight_rk * velocs_poroelastic_rk(1,iglob,i_stage)
 ! displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) + weight_rk * velocs_poroelastic_rk(2,iglob,i_stage)


 !       else if(i_stage==4)then

 !       velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (accels_poroelastic_rk(1,iglob,1) + 2.0d0 * accels_poroelastic_rk(1,iglob,2) + &
 !        2.0d0 * accels_poroelastic_rk(1,iglob,3) + accels_poroelastic_rk(1,iglob,4))

 !       velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) + 1.0d0 / 6.0d0 * &
 !       (accels_poroelastic_rk(2,iglob,1) + 2.0d0 * accels_poroelastic_rk(2,iglob,2) + &
 !        2.0d0 * accels_poroelastic_rk(2,iglob,3) + accels_poroelastic_rk(2,iglob,4))
!
 !       displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (velocs_poroelastic_rk(1,iglob,1) + 2.0d0 * velocs_poroelastic_rk(1,iglob,2) + &
 !        2.0d0 * velocs_poroelastic_rk(1,iglob,3) + velocs_poroelastic_rk(1,iglob,4))
!
 !       displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) + 1.0d0 / 6.0d0 * &
 !       (velocs_poroelastic_rk(2,iglob,1) + 2.0d0 * velocs_poroelastic_rk(2,iglob,2) + &
 !        2.0d0 * velocs_poroelastic_rk(2,iglob,3) + velocs_poroelastic_rk(2,iglob,4))

 !       endif

 !     endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if(SIMULATION_TYPE == 3) then
              b_veloc_elastic(1,iglob) = b_veloc_elastic(1,iglob) - b_deltatover2*b_accel_elastic(1,iglob)
              b_veloc_elastic(3,iglob) = b_veloc_elastic(3,iglob) - b_deltatover2*b_accel_elastic(3,iglob)
              b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) / rmass_inverse_elastic_one(iglob)
              b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) / rmass_inverse_elastic_three(iglob)
              ! recovering original velocities and accelerations on boundaries (poro side)
              b_velocs_poroelastic(1,iglob) = b_velocs_poroelastic(1,iglob) - b_deltatover2*b_accels_poroelastic(1,iglob)
              b_velocs_poroelastic(2,iglob) = b_velocs_poroelastic(2,iglob) - b_deltatover2*b_accels_poroelastic(2,iglob)
              b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
              b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)
              ! assembling accelerations
              b_accel_elastic(1,iglob) = ( b_accel_elastic(1,iglob) + b_accels_poroelastic(1,iglob) ) / &
                                   ( 1.0/rmass_inverse_elastic_one(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
              b_accel_elastic(3,iglob) = ( b_accel_elastic(3,iglob) + b_accels_poroelastic(2,iglob) ) / &
                                   ( 1.0/rmass_inverse_elastic_three(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
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
            endif !if(SIMULATION_TYPE == 3)

          endif !if(icount(iglob) ==1)

        enddo

      enddo
    endif

   enddo !LDDRK or RK

! ********************************************************************************************
!                       reading lastframe for adjoint/kernels calculation
! ********************************************************************************************
    if(it == 1 .and. SIMULATION_TYPE == 3) then

      ! acoustic medium
      if(any_acoustic) then
        write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
        open(unit=55,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
        do j=1,nglob
          read(55) b_potential_acoustic(j),&
                  b_potential_dot_acoustic(j),&
                  b_potential_dot_dot_acoustic(j)
          enddo
        close(55)

        ! free surface for an acoustic medium
        if ( nelem_acoustic_surface > 0 ) then
          call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                            b_potential_acoustic,acoustic_surface, &
                                            ibool,nelem_acoustic_surface,nglob,nspec,this_ibool_is_a_periodic_edge)
        endif
      endif

      ! elastic medium
      if(any_elastic) then
        write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
        open(unit=55,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
        if(p_sv)then !P-SV waves
          do j=1,nglob
            read(55) (b_displ_elastic(i,j), i=1,NDIM), &
                      (b_veloc_elastic(i,j), i=1,NDIM), &
                      (b_accel_elastic(i,j), i=1,NDIM)
          enddo
          b_displ_elastic(3,:) = b_displ_elastic(2,:)
          b_displ_elastic(2,:) = 0._CUSTOM_REAL
          b_veloc_elastic(3,:) = b_veloc_elastic(2,:)
          b_veloc_elastic(2,:) = 0._CUSTOM_REAL
          b_accel_elastic(3,:) = b_accel_elastic(2,:)
          b_accel_elastic(2,:) = 0._CUSTOM_REAL
        else !SH (membrane) waves
          do j=1,nglob
            read(55) b_displ_elastic(2,j), &
                      b_veloc_elastic(2,j), &
                      b_accel_elastic(2,j)
          enddo
          b_displ_elastic(1,:) = 0._CUSTOM_REAL
          b_displ_elastic(3,:) = 0._CUSTOM_REAL
          b_veloc_elastic(1,:) = 0._CUSTOM_REAL
          b_veloc_elastic(3,:) = 0._CUSTOM_REAL
          b_accel_elastic(1,:) = 0._CUSTOM_REAL
          b_accel_elastic(3,:) = 0._CUSTOM_REAL
        endif
        close(55)
      endif

      ! poroelastic medium
      if(any_poroelastic) then
        write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_s',myrank,'.bin'
        open(unit=55,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
        write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_w',myrank,'.bin'
        open(unit=56,file='OUTPUT_FILES/'//outputname,status='old',action='read',form='unformatted')
        do j=1,nglob
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

    endif ! if(it == 1 .and. SIMULATION_TYPE == 3)

!<NOISE_TOMOGRAPHY

  if ( NOISE_TOMOGRAPHY == 1 ) then
    call save_surface_movie_noise(NOISE_TOMOGRAPHY,p_sv,it,NSTEP,nglob,displ_elastic)

  else if ( NOISE_TOMOGRAPHY == 2 .and. save_everywhere ) then
    call save_surface_movie_noise(NOISE_TOMOGRAPHY,p_sv,it,NSTEP,nglob,displ_elastic)

  else if ( NOISE_TOMOGRAPHY == 3 .and. save_everywhere ) then
    if (it==1) &
      open(unit=500,file='OUTPUT_FILES/NOISE_TOMOGRAPHY/phi',access='direct', &
      recl=nglob*CUSTOM_REAL,action='write',iostat=ios)
    if( ios /= 0) write(*,*) 'Error retrieving ensemble forward wavefield.'
    if(p_sv) then
      call exit_mpi('P-SV case not yet implemented.')
    else
      read(unit=500,rec=NSTEP-it+1) b_displ_elastic(2,:)
    endif

  endif



!>NOISE_TOMOGRAPHY

! ********************************************************************************************
!                                      kernels calculation
! ********************************************************************************************
    if(any_elastic .and. SIMULATION_TYPE == 3) then ! kernels calculation
      do iglob = 1,nglob
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

    if(any_poroelastic .and. SIMULATION_TYPE == 3) then
      do iglob =1,nglob
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
    if(output_energy) then

      call compute_energy(displ_elastic,veloc_elastic, &
                        displs_poroelastic,velocs_poroelastic, &
                        displw_poroelastic,velocw_poroelastic, &
                        xix,xiz,gammax,gammaz,jacobian,ibool,elastic,poroelastic,hprime_xx,hprime_zz, &
                        AXISYM,nglob,coord,is_on_the_axis,hprimeBar_xx, &                                            !axisym
                        nspec,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                        assign_external_model,kmato,poroelastcoef,density,porosity,tortuosity, &
                        vpext,vsext,rhoext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext, &
                        anisotropic,anisotropy,wxgll,wzgll,numat, &
                        pressure_element,vector_field_element,e1,e11, &
                        potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ATTENUATION_VISCOELASTIC_SOLID,Mu_nu1,Mu_nu2,N_SLS,p_sv,kinetic_energy,potential_energy)

#ifdef USE_MPI
      call MPI_REDUCE(kinetic_energy, kinetic_energy_total, 1, CUSTOM_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ier)
      call MPI_REDUCE(potential_energy, potential_energy_total, 1, CUSTOM_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD, ier)
#else
      kinetic_energy_total = kinetic_energy
      potential_energy_total = potential_energy
#endif

! save kinetic, potential and total energy for this time step in external file
      if(myrank == 0) write(IOUT_ENERGY,*) real(dble(it-1)*deltat - t0,4),real(kinetic_energy_total,4), &
                     real(potential_energy_total,4),real(kinetic_energy_total + potential_energy_total,4)

    endif

!----  display time step and max of norm of displacement
    if(mod(it,NSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
      call check_stability(myrank,time,it,NSTEP,NOISE_TOMOGRAPHY, &
                        nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                        any_elastic_glob,any_elastic,displ_elastic, &
                        any_poroelastic_glob,any_poroelastic, &
                        displs_poroelastic,displw_poroelastic, &
                        any_acoustic_glob,any_acoustic,potential_acoustic, &
                        timestamp_seconds_start)

    endif

!---- loop on all the receivers to compute and store the seismograms
    if(mod(it-1,subsamp_seismos) == 0) then

! update position in seismograms
    seismo_current = seismo_current + 1

    do irecloc = 1,nrecloc

      irec = recloc(irecloc)

      ispec = ispec_selected_rec(irec)

      ! compute pressure in this element if needed
      if(seismotype == 4) then

        call compute_pressure_one_element(pressure_element,potential_dot_dot_acoustic,displ_elastic,&
              displs_poroelastic,displw_poroelastic,elastic,poroelastic,&
              xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec, &
              AXISYM,nglob,coord,jacobian,is_on_the_axis,hprimeBar_xx, &                                             !axisym
              nglob_acoustic,nglob_elastic,nglob_poroelastic,assign_external_model, &
              numat,kmato,density,porosity,tortuosity,poroelastcoef,vpext,vsext,rhoext, &
              c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,anisotropic,anisotropy,ispec,e1,e11, &
              ATTENUATION_VISCOELASTIC_SOLID,Mu_nu1,Mu_nu2,N_SLS)

      else if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then

        ! for acoustic medium, compute vector field from gradient of potential for seismograms
        if(seismotype == 1) then
          call compute_vector_one_element(vector_field_element,potential_acoustic, &
                              displ_elastic,displs_poroelastic,&
                              elastic,poroelastic,xix,xiz,gammax,gammaz, &
                              ibool,hprime_xx,hprime_zz, &
                              AXISYM,is_on_the_axis,hprimeBar_xx, &                                                  !axisym
                              nspec,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                              ispec,numat,kmato,density,rhoext,assign_external_model)
        else if(seismotype == 2) then
          call compute_vector_one_element(vector_field_element,potential_dot_acoustic, &
                              veloc_elastic,velocs_poroelastic, &
                              elastic,poroelastic,xix,xiz,gammax,gammaz, &
                              ibool,hprime_xx,hprime_zz, &
                              AXISYM,is_on_the_axis,hprimeBar_xx, &                                                  !axisym
                              nspec,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                              ispec,numat,kmato,density,rhoext,assign_external_model)
        else if(seismotype == 3) then
          call compute_vector_one_element(vector_field_element,potential_dot_dot_acoustic, &
                              accel_elastic,accels_poroelastic, &
                              elastic,poroelastic,xix,xiz,gammax,gammaz, &
                              ibool,hprime_xx,hprime_zz, &
                              AXISYM,is_on_the_axis,hprimeBar_xx, &                                                  !axisym
                              nspec,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                              ispec,numat,kmato,density,rhoext,assign_external_model)
        endif

      else if(seismotype == 5) then
        call compute_curl_one_element(curl_element,displ_elastic, &
                            displs_poroelastic,elastic,poroelastic, &
                            xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz, &
                            nspec,nglob_elastic,nglob_poroelastic,ispec)
      endif

      ! perform the general interpolation using Lagrange polynomials
      valux = ZERO
      valuy = ZERO
      valuz = ZERO

      valcurl = ZERO

      dxd = 0
      dyd = 0
      dzd = 0

      dcurld = 0

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
            else if(elastic(ispec)) then
              dxd = displ_elastic(1,iglob)
              dyd = displ_elastic(2,iglob)
              dzd = displ_elastic(3,iglob)
            endif

          else if(seismotype == 2) then

            if(poroelastic(ispec)) then
              dxd = velocs_poroelastic(1,iglob)
              dzd = velocs_poroelastic(2,iglob)
            else if(elastic(ispec)) then
              dxd = veloc_elastic(1,iglob)
              dyd = veloc_elastic(2,iglob)
              dzd = veloc_elastic(3,iglob)
            endif

          else if(seismotype == 3) then

            if(poroelastic(ispec)) then
              dxd = accels_poroelastic(1,iglob)
              dzd = accels_poroelastic(2,iglob)
            else if(elastic(ispec)) then
              dxd = accel_elastic(1,iglob)
              dyd = accel_elastic(2,iglob)
              dzd = accel_elastic(3,iglob)
            endif

          else if(seismotype == 5) then

            if(poroelastic(ispec)) then
              dxd = displs_poroelastic(1,iglob)
              dzd = displs_poroelastic(2,iglob)
            else if(elastic(ispec)) then
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

      ! check for edge effects
      if(seismo_current < 1 .or. seismo_current > NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos) &
        stop 'error: seismo_current out of bounds in recording of seismograms'

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

    endif

!----- writing the kernels
    ! kernels output
    if(SIMULATION_TYPE == 3) then

      if(any_acoustic) then

        do ispec = 1, nspec
          if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then
            do j = 1, NGLLZ
              do i = 1, NGLLX
                iglob = ibool(i,j,ispec)
                if (.not. assign_external_model) then
                   kappal_ac_global(iglob) = poroelastcoef(3,1,kmato(ispec))
                   rhol_ac_global(iglob) = density(1,kmato(ispec))
                else
                   rhol_ac_global(iglob)   = rhoext(i,j,ispec)
                   kappal_ac_global(iglob) = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec)
                endif

! calcul the displacement by computing the gradient of potential / rho
! and calcul the acceleration by computing the gradient of potential_dot_dot / rho
                tempx1l = ZERO
                tempx2l = ZERO
                b_tempx1l = ZERO
                b_tempx2l = ZERO
                bb_tempx1l = ZERO
                bb_tempx2l = ZERO
                do k = 1,NGLLX
                  ! derivative along x
                  !tempx1l = tempx1l + potential_dot_dot_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
                  tempx1l = tempx1l + potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k) !!! YANGL
                  b_tempx1l = b_tempx1l + b_potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
                  bb_tempx1l = bb_tempx1l + b_potential_dot_dot_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
                  ! derivative along z
                  !tempx2l = tempx2l + potential_dot_dot_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
                  tempx2l = tempx2l + potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k) !!! YANGL
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
                !<YANGL
                !!!! old expression (from elastic kernels)
                !!!rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) - rhol_ac_global(iglob)  * &
                !!!           dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * deltat
                !!!kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) - kappal_ac_global(iglob) * &
                !!!           potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * &
                !!!           b_potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob)&
                !!!           * deltat
                !!!! new expression (from PDE-constrained optimization, coupling terms changed as well)
                rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + rhol_ac_global(iglob)  * &
                           dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * deltat
                kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) + kappal_ac_global(iglob) * &
                           potential_acoustic(iglob)/kappal_ac_global(iglob) * &
                           b_potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * deltat
                !>YANGL
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
                if (.not. assign_external_model) then
                   mul_global(iglob) = poroelastcoef(2,1,kmato(ispec))
                   kappal_global(iglob) = poroelastcoef(3,1,kmato(ispec)) &
                                       - 4._CUSTOM_REAL*mul_global(iglob)/3._CUSTOM_REAL
                   rhol_global(iglob) = density(1,kmato(ispec))
                else
                   rhol_global(iglob)   = rhoext(i,j,ispec)
                   mul_global(iglob)    = rhoext(i,j,ispec)*vsext(i,j,ispec)*vsext(i,j,ispec)
                   kappal_global(iglob) = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) &
                                       -4._CUSTOM_REAL*mul_global(iglob)/3._CUSTOM_REAL
                endif

                rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - rhol_global(iglob)  * rho_k(iglob) * deltat
                mu_kl(i,j,ispec) =  mu_kl(i,j,ispec) - TWO * mul_global(iglob) * mu_k(iglob) * deltat
                kappa_kl(i,j,ispec) = kappa_kl(i,j,ispec) - kappal_global(iglob) * kappa_k(iglob) * deltat
                !
                rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + kappa_kl(i,j,ispec) + mu_kl(i,j,ispec)
                beta_kl(i,j,ispec) = TWO * (mu_kl(i,j,ispec) - 4._CUSTOM_REAL * mul_global(iglob) &
                    / (3._CUSTOM_REAL * kappal_global(iglob)) * kappa_kl(i,j,ispec))
                alpha_kl(i,j,ispec) = TWO * (1._CUSTOM_REAL + 4._CUSTOM_REAL * mul_global(iglob)/&
                     (3._CUSTOM_REAL * kappal_global(iglob))) * kappa_kl(i,j,ispec)
                rhorho_el_hessian_final1(i,j,ispec) = rhorho_el_hessian_final1(i,j,ispec) &
                                                  + rhorho_el_hessian_temp1(iglob) * deltat
                rhorho_el_hessian_final2(i,j,ispec) = rhorho_el_hessian_final2(i,j,ispec) &
                                                  + rhorho_el_hessian_temp2(iglob) * deltat

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
                sm_kl(i,j,ispec) = sm_kl(i,j,ispec) - &
                        deltat * rhol_f_global(iglob)*tortl_global(iglob)/phil_global(iglob) * sm_k(iglob)
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
                dd1 = (1._CUSTOM_REAL+rholb/rhol_f_global(iglob))*ratio**2 &
                      + 2._CUSTOM_REAL*ratio &
                      + tortl_global(iglob)/phil_global(iglob)

                rhobb_kl(i,j,ispec) = rhob_kl(i,j,ispec) - &
                      phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*B_biot) * &
                      (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil_global(iglob) / &
                      tortl_global(iglob)*ratio +1._CUSTOM_REAL)/dd1 + &
                      (rhol_bar_global(iglob)**2*ratio**2/rhol_f_global(iglob)**2*(phil_global(iglob) / &
                      tortl_global(iglob)*ratio+1)*(phil_global(iglob)/tortl_global(iglob)*ratio + &
                      phil_global(iglob)/tortl_global(iglob) * &
                      (1+rhol_f_global(iglob)/rhol_bar_global(iglob))-1) )/dd1**2 ) - &
                      FOUR_THIRDS*cssquare )*Bb_kl(i,j,ispec) - &
                      rhol_bar_global(iglob)*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                      (phil_global(iglob)/tortl_global(iglob)*ratio + &
                      1._CUSTOM_REAL)**2/dd1**2*Mb_kl(i,j,ispec) + &
                      rhol_bar_global(iglob)*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                      (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                      phil_global(iglob)*ratio/tortl_global(iglob)*(phil_global(iglob) / &
                      tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                      (1+rhol_bar_global(iglob)*ratio/rhol_f_global(iglob))/dd1**2)*Cb_kl(i,j,ispec)+ &
                      phil_global(iglob)*rhol_f_global(iglob)*cssquare / &
                      (tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
                rhofbb_kl(i,j,ispec) = rhofb_kl(i,j,ispec) + &
                        phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*B_biot) * &
                       (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil_global(iglob)/ &
                       tortl_global(iglob)*ratio +1._CUSTOM_REAL)/dd1+&
                       (rhol_bar_global(iglob)**2*ratio**2/rhol_f_global(iglob)**2*(phil_global(iglob)/ &
                       tortl_global(iglob)*ratio+1)*(phil_global(iglob)/tortl_global(iglob)*ratio+ &
                       phil_global(iglob)/tortl_global(iglob)*&
                       (1+rhol_f_global(iglob)/rhol_bar_global(iglob))-1) )/dd1**2 )- &
                       FOUR_THIRDS*cssquare )*Bb_kl(i,j,ispec) + &
                        rhol_bar_global(iglob)*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                       (phil_global(iglob)/tortl_global(iglob)*ratio + &
                       1._CUSTOM_REAL)**2/dd1**2*Mb_kl(i,j,ispec) - &
                       rhol_bar_global(iglob)*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                       phil_global(iglob)*ratio/tortl_global(iglob)*(phil_global(iglob)/ &
                       tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                       (1+rhol_bar_global(iglob)*ratio/rhol_f_global(iglob))/dd1**2)*Cb_kl(i,j,ispec)- &
                       phil_global(iglob)*rhol_f_global(iglob)*cssquare/ &
                       (tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
                phib_kl(i,j,ispec) = phi_kl(i,j,ispec) - &
                       phil_global(iglob)*rhol_bar_global(iglob)/(tortl_global(iglob)*B_biot) &
                       * ( cpIsquare - rhol_f_global(iglob)/rhol_bar_global(iglob)*cpIIsquare- &
                       (cpIsquare-cpIIsquare)*( (TWO*ratio**2*phil_global(iglob)/ &
                       tortl_global(iglob) + (1._CUSTOM_REAL+&
                       rhol_f_global(iglob)/rhol_bar_global(iglob))* &
                       (TWO*ratio*phil_global(iglob)/tortl_global(iglob)+&
                       1._CUSTOM_REAL))/dd1 + (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                       1._CUSTOM_REAL)*(phil_global(iglob)*&
                       ratio/tortl_global(iglob)+phil_global(iglob)/tortl_global(iglob)* &
                       (1._CUSTOM_REAL+rhol_f_global(iglob)/&
                       rhol_bar_global(iglob))-1._CUSTOM_REAL)*((1._CUSTOM_REAL+ &
                       rhol_bar_global(iglob)/rhol_f_global(iglob)-&
                       TWO*phil_global(iglob)/tortl_global(iglob))*ratio**2+TWO*ratio)/dd1**2 ) - &
                       FOUR_THIRDS*rhol_f_global(iglob)*cssquare/rhol_bar_global(iglob) )*Bb_kl(i,j,ispec) + &
                       rhol_f_global(iglob)/M_biot * (cpIsquare-cpIIsquare)*(&
                       TWO*ratio*(phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2*( &
                       (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                       rhol_f_global(iglob)-TWO*phil_global(iglob)/tortl_global(iglob))*ratio**2+TWO*ratio)/dd1**2 &
                       )*Mb_kl(i,j,ispec) + &
                       phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*C_biot)* &
                       (cpIsquare-cpIIsquare)*ratio* (&
                       (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob)*ratio)/dd1 - &
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)* &
                       (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                       rhol_f_global(iglob)*ratio)*((1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-TWO*&
                       phil_global(iglob)/tortl_global(iglob))*ratio+TWO)/dd1**2&
                        )*Cb_kl(i,j,ispec) -&
                       phil_global(iglob)*rhol_f_global(iglob)*cssquare &
                       /(tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
                cpI_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIsquare/B_biot*rhol_bar_global(iglob)*( &
                       1._CUSTOM_REAL-phil_global(iglob)/tortl_global(iglob) + &
                       (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                       1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*&
                       ratio+phil_global(iglob)/tortl_global(iglob)* &
                       (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-&
                       1._CUSTOM_REAL)/dd1 &
                        )* Bb_kl(i,j,ispec) +&
                       2._CUSTOM_REAL*cpIsquare*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot) *&
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2/dd1*Mb_kl(i,j,ispec)+&
                       2._CUSTOM_REAL*cpIsquare*rhol_f_global(iglob)/C_biot * &
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)* &
                       (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                       rhol_f_global(iglob)*ratio)/dd1*Cb_kl(i,j,ispec)
                cpII_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIIsquare*rhol_bar_global(iglob)/B_biot * (&
                       phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*rhol_bar_global(iglob)) - &
                       (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                       1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*&
                       ratio+phil_global(iglob)/tortl_global(iglob)* &
                       (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-&
                       1._CUSTOM_REAL)/dd1  ) * Bb_kl(i,j,ispec) +&
                       2._CUSTOM_REAL*cpIIsquare*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot) * (&
                       1._CUSTOM_REAL - (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                       1._CUSTOM_REAL)**2/dd1  )*Mb_kl(i,j,ispec) + &
                       2._CUSTOM_REAL*cpIIsquare*rhol_f_global(iglob)/C_biot * (&
                       1._CUSTOM_REAL - (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                       1._CUSTOM_REAL)*(1._CUSTOM_REAL+&
                       rhol_bar_global(iglob)/rhol_f_global(iglob)*ratio)/dd1  )*Cb_kl(i,j,ispec)
                cs_kl(i,j,ispec) = - 8._CUSTOM_REAL/3._CUSTOM_REAL*cssquare* &
                       rhol_bar_global(iglob)/B_biot*(1._CUSTOM_REAL-&
                       phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)* &
                       rhol_bar_global(iglob)))*Bb_kl(i,j,ispec) + &
                       2._CUSTOM_REAL*(rhol_bar_global(iglob)-rhol_f_global(iglob)*&
                       phil_global(iglob)/tortl_global(iglob))/&
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
                       (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2*( &
                       (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                       rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob))*ratio+ &
                       1._CUSTOM_REAL)/dd1**2 )*Mb_kl(i,j,ispec) +&
                       ratio*rhol_f_global(iglob)/C_biot*(cpIsquare-cpIIsquare) * (&
                       (2._CUSTOM_REAL*phil_global(iglob)*rhol_bar_global(iglob)* &
                       ratio/(tortl_global(iglob)*rhol_f_global(iglob))+&
                       phil_global(iglob)/tortl_global(iglob)+rhol_bar_global(iglob)/rhol_f_global(iglob))/dd1 - &
                       2._CUSTOM_REAL*phil_global(iglob)/tortl_global(iglob)*(phil_global(iglob)/tortl_global(iglob)*ratio+&
                       1._CUSTOM_REAL)*(1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)*ratio)*((1._CUSTOM_REAL+&
                       rhol_bar_global(iglob)/rhol_f_global(iglob)- &
                       phil_global(iglob)/tortl_global(iglob))*ratio+1._CUSTOM_REAL)/&
                       dd1**2 )*Cb_kl(i,j,ispec)

              enddo
            enddo
          endif
        enddo

      endif ! if(any_poroelastic)

    endif ! if(SIMULATION_TYPE == 3)


!
!----  display results at given time steps
!
    if(mod(it,NSTEP_BETWEEN_OUTPUT_IMAGES) == 0 .or. it == 5 .or. it == NSTEP) then

!
! kernels output files
!

      if(SIMULATION_TYPE == 3 .and. it == NSTEP) then

        if ( myrank == 0 ) then
          write(IOUT,*) 'Writing Kernels file'
        endif

        if(any_acoustic) then
          if(.not. save_ASCII_kernels)then
             write(95)coord
             write(95)rho_ac_kl
             write(95)kappa_ac_kl
             write(96)coord
             write(96)rho_ac_kl
             write(96)alpha_ac_kl
          else
            do ispec = 1, nspec
              do j = 1, NGLLZ
                do i = 1, NGLLX
                  iglob = ibool(i,j,ispec)
                  xx = coord(1,iglob)
                  zz = coord(2,iglob)
                  write(95,'(4e15.5e4)')xx,zz,rho_ac_kl(i,j,ispec),kappa_ac_kl(i,j,ispec)
                  write(96,'(4e15.5e4)')xx,zz,rhop_ac_kl(i,j,ispec),alpha_ac_kl(i,j,ispec)
                  !write(96,'(4e15.5e4)')rhorho_ac_hessian_final1(i,j,ispec), rhorho_ac_hessian_final2(i,j,ispec),&
                  !                rhop_ac_kl(i,j,ispec),alpha_ac_kl(i,j,ispec)
                enddo
              enddo
            enddo
          endif
          close(95)
          close(96)
        endif

        if(any_elastic) then
          if(.not. save_ASCII_kernels)then
             write(97)coord
             write(97)rho_kl
             write(97)kappa_kl
             write(97)mu_kl
             write(98)coord
             write(98)rhop_kl
             write(98)alpha_kl
             write(98)beta_kl
          else
            do ispec = 1, nspec
              do j = 1, NGLLZ
                do i = 1, NGLLX
                  iglob = ibool(i,j,ispec)
                  xx = coord(1,iglob)
                  zz = coord(2,iglob)
                  write(97,'(5e15.5e4)')xx,zz,rho_kl(i,j,ispec),kappa_kl(i,j,ispec),mu_kl(i,j,ispec)
                  write(98,'(5e15.5e4)')xx,zz,rhop_kl(i,j,ispec),alpha_kl(i,j,ispec),beta_kl(i,j,ispec)
                  !write(98,'(5e15.5e4)')rhorho_el_hessian_final1(i,j,ispec), rhorho_el_hessian_final2(i,j,ispec),&
                  !                    rhop_kl(i,j,ispec),alpha_kl(i,j,ispec),beta_kl(i,j,ispec)
                enddo
              enddo
            enddo
          endif
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

!<NOISE_TOMOGRAPHY

      if (NOISE_TOMOGRAPHY == 3 .and. output_wavefields_noise) then

        !load ensemble forward source
        inquire(unit=500,exist=ex,opened=od)
        if (.not. od) &
          open(unit=500,file='OUTPUT_FILES/NOISE_TOMOGRAPHY/eta',access='direct', &
          recl=nglob*CUSTOM_REAL,action='write',iostat=ios)
        read(unit=500,rec=it) surface_movie_y_noise

        !load product of fwd, adj wavefields
        call spec2glob(nspec,nglob,ibool,rho_kl,noise_output_rhokl)

        !write text file
        noise_output_array(1,:) = surface_movie_y_noise(:) * mask_noise(:)
        noise_output_array(2,:) = b_displ_elastic(2,:)
        noise_output_array(3,:) = accel_elastic(2,:)
        noise_output_array(4,:) = rho_k(:)
        noise_output_array(5,:) = noise_output_rhokl(:)
        write(noise_output_file,"('OUTPUT_FILES/snapshot_all_',i6.6)") it
        call snapshots_noise(noise_output_ncol,nglob,noise_output_file,noise_output_array)

      endif

!>NOISE_TOMOGRAPHY


!
!----  PostScript display
!
      if(output_postscript_snapshot) then

        if (myrank == 0) then
          write(IOUT,*)
          write(IOUT,*) 'Writing PostScript vector plot for time step ',it
        endif

        if(imagetype_postscript == 1 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'drawing displacement vector as small arrows...'

          call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic,&
                          elastic,poroelastic,vector_field_display, &
                          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz, &
                          AXISYM,is_on_the_axis,hprimeBar_xx, &                                                      !axisym
                          nspec,nglob,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                          numat,kmato,density,rhoext,assign_external_model)

          call plotpost(vector_field_display,coord,vpext,x_source,z_source,st_xval,st_zval, &
                      it,deltat,coorg,xinterp,zinterp,shape2D_display, &
                      Uxinterp,Uzinterp,flagrange,density,porosity,tortuosity,&
                      AXISYM,is_on_the_axis,flagrange_GLJ, &                                                        !axisym
                      poroelastcoef,knods,kmato,ibool, &
                      numabs,codeabs,typeabs,anyabs,nelem_acoustic_surface,acoustic_edges, &
                      simulation_title,nglob,npgeo,vpImin,vpImax,nrec,NSOURCES, &
                      colors,numbers,subsamp_postscript,imagetype_postscript,interpol,meshvect,modelvect, &
                      boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,numat,pointsdisp, &
                      nspec,ngnod,coupled_acoustic_elastic,coupled_acoustic_poro,coupled_elastic_poro, &
                      any_acoustic,any_poroelastic,plot_lowerleft_corner_only, &
                      fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,num_fluid_solid_edges,&
                      fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge,num_fluid_poro_edges, &
                      solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge,num_solid_poro_edges, &
                      poroelastic,myrank,nproc,ier,&
                      d1_coorg_send_ps_velocity_model,d2_coorg_send_ps_velocity_model, &
                      d1_coorg_recv_ps_velocity_model,d2_coorg_recv_ps_velocity_model, &
                      d1_RGB_send_ps_velocity_model,d2_RGB_send_ps_velocity_model, &
                      d1_RGB_recv_ps_velocity_model,d2_RGB_recv_ps_velocity_model, &
                      coorg_send_ps_velocity_model,RGB_send_ps_velocity_model, &
                      coorg_recv_ps_velocity_model,RGB_recv_ps_velocity_model, &
                      d1_coorg_send_ps_element_mesh,d2_coorg_send_ps_element_mesh, &
                      d1_coorg_recv_ps_element_mesh,d2_coorg_recv_ps_element_mesh, &
                      d1_color_send_ps_element_mesh,d1_color_recv_ps_element_mesh, &
                      coorg_send_ps_element_mesh,color_send_ps_element_mesh, &
                      coorg_recv_ps_element_mesh,color_recv_ps_element_mesh, &
                      d1_coorg_send_ps_abs,d1_coorg_recv_ps_abs,d2_coorg_send_ps_abs,d2_coorg_recv_ps_abs, &
                      coorg_send_ps_abs,coorg_recv_ps_abs, &
                      d1_coorg_send_ps_free_surface,d1_coorg_recv_ps_free_surface, &
                      d2_coorg_send_ps_free_surface,d2_coorg_recv_ps_free_surface, &
                      coorg_send_ps_free_surface,coorg_recv_ps_free_surface, &
                      d1_coorg_send_ps_vector_field,d1_coorg_recv_ps_vector_field, &
                      d2_coorg_send_ps_vector_field,d2_coorg_recv_ps_vector_field, &
                      coorg_send_ps_vector_field,coorg_recv_ps_vector_field,US_LETTER,is_PML)

        else if(imagetype_postscript == 2 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'drawing velocity vector as small arrows...'

          call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic,&
                          elastic,poroelastic,vector_field_display, &
                          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz, &
                          AXISYM,is_on_the_axis,hprimeBar_xx, &                                                      !axisym
                          nspec,nglob,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                          numat,kmato,density,rhoext,assign_external_model)

          call plotpost(vector_field_display,coord,vpext,x_source,z_source,st_xval,st_zval, &
                      it,deltat,coorg,xinterp,zinterp,shape2D_display, &
                      Uxinterp,Uzinterp,flagrange,density,porosity,tortuosity,&
                      AXISYM,is_on_the_axis,flagrange_GLJ, &                                                        !axisym
                      poroelastcoef,knods,kmato,ibool, &
                      numabs,codeabs,typeabs,anyabs,nelem_acoustic_surface,acoustic_edges, &
                      simulation_title,nglob,npgeo,vpImin,vpImax,nrec,NSOURCES, &
                      colors,numbers,subsamp_postscript,imagetype_postscript,interpol,meshvect,modelvect, &
                      boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,numat,pointsdisp, &
                      nspec,ngnod,coupled_acoustic_elastic,coupled_acoustic_poro,coupled_elastic_poro, &
                      any_acoustic,any_poroelastic,plot_lowerleft_corner_only, &
                      fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,num_fluid_solid_edges,&
                      fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge,num_fluid_poro_edges, &
                      solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge,num_solid_poro_edges, &
                      poroelastic,myrank,nproc,ier,&
                      d1_coorg_send_ps_velocity_model,d2_coorg_send_ps_velocity_model, &
                      d1_coorg_recv_ps_velocity_model,d2_coorg_recv_ps_velocity_model, &
                      d1_RGB_send_ps_velocity_model,d2_RGB_send_ps_velocity_model, &
                      d1_RGB_recv_ps_velocity_model,d2_RGB_recv_ps_velocity_model, &
                      coorg_send_ps_velocity_model,RGB_send_ps_velocity_model, &
                      coorg_recv_ps_velocity_model,RGB_recv_ps_velocity_model, &
                      d1_coorg_send_ps_element_mesh,d2_coorg_send_ps_element_mesh, &
                      d1_coorg_recv_ps_element_mesh,d2_coorg_recv_ps_element_mesh, &
                      d1_color_send_ps_element_mesh,d1_color_recv_ps_element_mesh, &
                      coorg_send_ps_element_mesh,color_send_ps_element_mesh, &
                      coorg_recv_ps_element_mesh,color_recv_ps_element_mesh, &
                      d1_coorg_send_ps_abs,d1_coorg_recv_ps_abs,d2_coorg_send_ps_abs,d2_coorg_recv_ps_abs, &
                      coorg_send_ps_abs,coorg_recv_ps_abs, &
                      d1_coorg_send_ps_free_surface,d1_coorg_recv_ps_free_surface, &
                      d2_coorg_send_ps_free_surface,d2_coorg_recv_ps_free_surface, &
                      coorg_send_ps_free_surface,coorg_recv_ps_free_surface, &
                      d1_coorg_send_ps_vector_field,d1_coorg_recv_ps_vector_field, &
                      d2_coorg_send_ps_vector_field,d2_coorg_recv_ps_vector_field, &
                      coorg_send_ps_vector_field,coorg_recv_ps_vector_field,US_LETTER,is_PML)

        else if(imagetype_postscript == 3 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'drawing acceleration vector as small arrows...'

          call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,accels_poroelastic,&
                          elastic,poroelastic,vector_field_display, &
                          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz, &
                          AXISYM,is_on_the_axis,hprimeBar_xx, &                                                      !axisym
                          nspec,nglob,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                          numat,kmato,density,rhoext,assign_external_model)

          call plotpost(vector_field_display,coord,vpext,x_source,z_source,st_xval,st_zval, &
                      it,deltat,coorg,xinterp,zinterp,shape2D_display, &
                      Uxinterp,Uzinterp,flagrange,density,porosity,tortuosity,&
                      AXISYM,is_on_the_axis,flagrange_GLJ, &                                                        !axisym
                      poroelastcoef,knods,kmato,ibool, &
                      numabs,codeabs,typeabs,anyabs,nelem_acoustic_surface,acoustic_edges, &
                      simulation_title,nglob,npgeo,vpImin,vpImax,nrec,NSOURCES, &
                      colors,numbers,subsamp_postscript,imagetype_postscript,interpol,meshvect,modelvect, &
                      boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,numat,pointsdisp, &
                      nspec,ngnod,coupled_acoustic_elastic,coupled_acoustic_poro,coupled_elastic_poro, &
                      any_acoustic,any_poroelastic,plot_lowerleft_corner_only, &
                      fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,num_fluid_solid_edges, &
                      fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge,num_fluid_poro_edges, &
                      solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge,num_solid_poro_edges, &
                      poroelastic,myrank,nproc,ier,&
                      d1_coorg_send_ps_velocity_model,d2_coorg_send_ps_velocity_model, &
                      d1_coorg_recv_ps_velocity_model,d2_coorg_recv_ps_velocity_model, &
                      d1_RGB_send_ps_velocity_model,d2_RGB_send_ps_velocity_model, &
                      d1_RGB_recv_ps_velocity_model,d2_RGB_recv_ps_velocity_model, &
                      coorg_send_ps_velocity_model,RGB_send_ps_velocity_model, &
                      coorg_recv_ps_velocity_model,RGB_recv_ps_velocity_model, &
                      d1_coorg_send_ps_element_mesh,d2_coorg_send_ps_element_mesh, &
                      d1_coorg_recv_ps_element_mesh,d2_coorg_recv_ps_element_mesh, &
                      d1_color_send_ps_element_mesh,d1_color_recv_ps_element_mesh, &
                      coorg_send_ps_element_mesh,color_send_ps_element_mesh, &
                      coorg_recv_ps_element_mesh,color_recv_ps_element_mesh, &
                      d1_coorg_send_ps_abs,d1_coorg_recv_ps_abs,d2_coorg_send_ps_abs,d2_coorg_recv_ps_abs, &
                      coorg_send_ps_abs,coorg_recv_ps_abs, &
                      d1_coorg_send_ps_free_surface,d1_coorg_recv_ps_free_surface, &
                      d2_coorg_send_ps_free_surface,d2_coorg_recv_ps_free_surface, &
                      coorg_send_ps_free_surface,coorg_recv_ps_free_surface, &
                      d1_coorg_send_ps_vector_field,d1_coorg_recv_ps_vector_field, &
                      d2_coorg_send_ps_vector_field,d2_coorg_recv_ps_vector_field, &
                      coorg_send_ps_vector_field,coorg_recv_ps_vector_field,US_LETTER,is_PML)

        else if(.not. p_sv) then
          call exit_MPI('cannot draw a SH scalar field as a vector plot, turn PostScript plots off')

        else
          call exit_MPI('wrong type for PostScript snapshots')
        endif

        if (myrank == 0 .and. imagetype_postscript /= 4 .and. p_sv) write(IOUT,*) 'PostScript file written'

      endif

!
!----  display color image
!
      if(output_color_image) then

        if (myrank == 0) then
          write(IOUT,*)
          write(IOUT,*) 'Creating color image of size ',NX_IMAGE_color,' x ',NZ_IMAGE_color,' for time step ',it
        endif

        if(imagetype_JPEG >= 1 .and. imagetype_JPEG <= 3) then

          if (myrank == 0) write(IOUT,*) 'drawing scalar image of part of the displacement vector...'
          call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic,&
                          elastic,poroelastic,vector_field_display, &
                          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz, &
                          AXISYM,is_on_the_axis,hprimeBar_xx, &                                                      !axisym
                          nspec,nglob,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                          numat,kmato,density,rhoext,assign_external_model)

        else if(imagetype_JPEG >= 4 .and. imagetype_JPEG <= 6) then

          if (myrank == 0) write(IOUT,*) 'drawing scalar image of part of the velocity vector...'
          call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic,&
                          elastic,poroelastic,vector_field_display, &
                          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz, &
                          AXISYM,is_on_the_axis,hprimeBar_xx, &                                                      !axisym
                          nspec,nglob,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                          numat,kmato,density,rhoext,assign_external_model)

        else if(imagetype_JPEG >= 7 .and. imagetype_JPEG <= 9) then

          if (myrank == 0) write(IOUT,*) 'drawing scalar image of part of the acceleration vector...'
          call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,accels_poroelastic,&
                          elastic,poroelastic,vector_field_display, &
                          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz, &
                          AXISYM,is_on_the_axis,hprimeBar_xx, &                                                      !axisym
                          nspec,nglob,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                          numat,kmato,density,rhoext,assign_external_model)

        else if(imagetype_JPEG == 10 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'drawing image of pressure field...'
          call compute_pressure_whole_medium(potential_dot_dot_acoustic,displ_elastic,&
                     displs_poroelastic,displw_poroelastic,elastic,poroelastic,vector_field_display, &
                     xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec, &
                     AXISYM,coord,jacobian,is_on_the_axis,hprimeBar_xx, &                                            !axisym
                     nglob,nglob_acoustic,nglob_elastic,nglob_poroelastic,assign_external_model, &
                     numat,kmato,density,porosity,tortuosity,poroelastcoef,vpext,vsext,rhoext, &
                     c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,anisotropic,anisotropy,e1,e11, &
                     ATTENUATION_VISCOELASTIC_SOLID,Mu_nu1,Mu_nu2,N_SLS)

        else if(imagetype_JPEG == 10 .and. .not. p_sv) then
          call exit_MPI('cannot draw pressure field for SH (membrane) waves')

        else
          call exit_MPI('wrong type for JPEG snapshots')
        endif

        image_color_data(:,:) = 0.d0

        do k = 1, nb_pixel_loc
          j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
          i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

          ! avoid edge effects
          if(i < 1) i = 1
          if(j < 1) j = 1

          if(i > NX_IMAGE_color) i = NX_IMAGE_color
          if(j > NZ_IMAGE_color) j = NZ_IMAGE_color

          if(p_sv) then ! P-SH waves, plot a component of vector, its norm, or else pressure
            if(iglob_image_color(i,j) /= -1) then
              if(imagetype_JPEG == 1 .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7) then
                image_color_data(i,j) = vector_field_display(1,iglob_image_color(i,j))  ! draw the X component of the vector

              else if(imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8) then
                image_color_data(i,j) = vector_field_display(3,iglob_image_color(i,j))  ! draw the Z component of the vector

              else if(imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9) then
                image_color_data(i,j) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2 + &
                                             vector_field_display(3,iglob_image_color(i,j))**2)  ! draw the norm of the vector

              else if(imagetype_JPEG == 10) then
! by convention we have stored pressure in the third component of the array
                image_color_data(i,j) = vector_field_display(3,iglob_image_color(i,j))

              else
                call exit_MPI('wrong type for JPEG snapshots')
              endif
            endif

          else ! SH (membrane) waves, plot y-component
            if(iglob_image_color(i,j) /= -1) image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))
          endif
        enddo

! assembling array image_color_data on process zero for color output
#ifdef USE_MPI
        if (nproc > 1) then
          if (myrank == 0) then
            do iproc = 1, nproc-1
              call MPI_RECV(data_pixel_recv(1),nb_pixel_per_proc(iproc+1), MPI_DOUBLE_PRECISION, &
                  iproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

              do k = 1, nb_pixel_per_proc(iproc+1)
                j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
                i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color

                ! avoid edge effects
                if(i < 1) i = 1
                if(j < 1) j = 1

                if(i > NX_IMAGE_color) i = NX_IMAGE_color
                if(j > NZ_IMAGE_color) j = NZ_IMAGE_color

                image_color_data(i,j) = data_pixel_recv(k)
              enddo
            enddo
          else
            do k = 1, nb_pixel_loc
              j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
              i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

              ! avoid edge effects
              if(i < 1) i = 1
              if(j < 1) j = 1

              if(i > NX_IMAGE_color) i = NX_IMAGE_color
              if(j > NZ_IMAGE_color) j = NZ_IMAGE_color

              if(p_sv) then ! P-SH waves, plot a component of vector, its norm, or else pressure

                if(imagetype_JPEG == 1 .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7) then
                  data_pixel_send(k) = vector_field_display(1,iglob_image_color(i,j))  ! draw the X component of the vector

                else if(imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8) then
                  data_pixel_send(k) = vector_field_display(3,iglob_image_color(i,j))  ! draw the Z component of the vector

                else if(imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9) then
                  data_pixel_send(k) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2 + &
                                            vector_field_display(3,iglob_image_color(i,j))**2)  ! draw the norm of the vector

                else if(imagetype_JPEG == 10) then
! by convention we have stored pressure in the third component of the array
                  data_pixel_send(k) = vector_field_display(3,iglob_image_color(i,j))

                else
                  call exit_MPI('wrong type for JPEG snapshots')
                endif

              else ! SH (membrane) waves, plot y-component
                if(iglob_image_color(i,j) /= -1) data_pixel_send(k) = vector_field_display(2,iglob_image_color(i,j))
              endif
            enddo

            call MPI_SEND(data_pixel_send(1),nb_pixel_loc,MPI_DOUBLE_PRECISION, 0, 43, MPI_COMM_WORLD, ier)
          endif
        endif
#endif

        if (myrank == 0) then
          call create_color_image(image_color_data,iglob_image_color, &
                  NX_IMAGE_color,NZ_IMAGE_color,it,isnapshot_number,cutsnaps,image_color_vp_display, &
                  USE_SNAPSHOT_NUMBER_IN_FILENAME,POWER_DISPLAY_COLOR, &
                  DRAW_SOURCES_AND_RECEIVERS,NSOURCES,nrec, &
                  ix_image_color_source,iy_image_color_source,ix_image_color_receiver,iy_image_color_receiver)
          write(IOUT,*) 'Color image created'
        endif

      endif
    endif  ! of display images at a given time step

!----------------------------------------------

! dump the full (local) wavefield to a file
! note: in the case of MPI, in the future it would be more convenient to output a single file rather than one for each myrank

    if(output_wavefield_dumps .and. (mod(it,NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS) == 0 .or. it == 5 .or. it == NSTEP)) then

      if (myrank == 0) then
        write(IOUT,*)
        write(IOUT,*) 'Dumping the wave field to a file for time step ',it
      endif

      if(this_is_the_first_time_we_dump) then

        allocate(mask_ibool(nglob))

! save the grid separately once and for all
        if(use_binary_for_wavefield_dumps) then
          write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.bin')") myrank
          open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='unknown', &
               action='write',recl=2*SIZE_REAL)
        else
          write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.txt')") myrank
          open(unit=27,file=wavefield_file,status='unknown',action='write')
        endif

        icounter = 0
        mask_ibool(:) = .false.
        do ispec = 1,nspec
          do j = 1,NGLLZ
            do i = 1,NGLLX
               iglob = ibool(i,j,ispec)
               if(.not. mask_ibool(iglob)) then
                 icounter = icounter + 1
                 mask_ibool(iglob) = .true.
                 if(use_binary_for_wavefield_dumps) then
                   write(27,rec=icounter) sngl(coord(1,iglob)),sngl(coord(2,iglob))
                 else
                   write(27,'(2e16.6)') coord(1,iglob),coord(2,iglob)
                 endif
               endif
            enddo
          enddo
        enddo

        close(27)

! save nglob to a file once and for all
        write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_value_of_nglob_',i3.3,'.txt')") myrank
        open(unit=27,file=wavefield_file,status='unknown',action='write')
        write(27,*) icounter
        close(27)
        if(icounter /= nglob) stop 'error: should have icounter == nglob in wavefield dumps'

        this_is_the_first_time_we_dump = .false.

      endif

        if(imagetype_wavefield_dumps == 1) then

          if (myrank == 0) write(IOUT,*) 'dumping the displacement vector...'
          call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic,&
                          elastic,poroelastic,vector_field_display, &
                          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz, &
                          AXISYM,is_on_the_axis,hprimeBar_xx, &                                                      !axisym
                          nspec,nglob,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                          numat,kmato,density,rhoext,assign_external_model)

        else if(imagetype_wavefield_dumps == 2) then

          if (myrank == 0) write(IOUT,*) 'dumping the velocity vector...'
          call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic,&
                          elastic,poroelastic,vector_field_display, &
                          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz, &
                          AXISYM,is_on_the_axis,hprimeBar_xx, &                                                      !axisym
                          nspec,nglob,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                          numat,kmato,density,rhoext,assign_external_model)

        else if(imagetype_wavefield_dumps == 3) then

          if (myrank == 0) write(IOUT,*) 'dumping the acceleration vector...'
          call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,accels_poroelastic,&
                          elastic,poroelastic,vector_field_display, &
                          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz, &
                          AXISYM,is_on_the_axis,hprimeBar_xx, &                                                      !axisym
                          nspec,nglob,nglob_acoustic,nglob_elastic,nglob_poroelastic, &
                          numat,kmato,density,rhoext,assign_external_model)

        else if(imagetype_wavefield_dumps == 4 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'dumping the pressure field...'
          call compute_pressure_whole_medium(potential_dot_dot_acoustic,displ_elastic,&
                     displs_poroelastic,displw_poroelastic,elastic,poroelastic,vector_field_display, &
                     xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec, &
                     AXISYM,coord,jacobian,is_on_the_axis,hprimeBar_xx, &                                            !axisym
                     nglob,nglob_acoustic,nglob_elastic,nglob_poroelastic,assign_external_model, &
                     numat,kmato,density,porosity,tortuosity,poroelastcoef,vpext,vsext,rhoext, &
                     c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,anisotropic,anisotropy,e1,e11, &
                     ATTENUATION_VISCOELASTIC_SOLID,Mu_nu1,Mu_nu2,N_SLS)

        else if(imagetype_wavefield_dumps == 4 .and. .not. p_sv) then
          call exit_MPI('cannot dump the pressure field for SH (membrane) waves')

        else
          call exit_MPI('wrong type of flag for wavefield dumping')
        endif

        if(use_binary_for_wavefield_dumps) then
          if(p_sv .and. .not. imagetype_wavefield_dumps == 4) then
            nb_of_values_to_save = 2
          else
            nb_of_values_to_save = 1
          endif
          write(wavefield_file,"('OUTPUT_FILES/wavefield',i7.7,'_',i2.2,'_',i3.3,'.bin')") it,SIMULATION_TYPE,myrank
          open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='unknown', &
                       action='write',recl=nb_of_values_to_save*SIZE_REAL)
        else
          write(wavefield_file,"('OUTPUT_FILES/wavefield',i7.7,'_',i2.2,'_',i3.3,'.txt')") it,SIMULATION_TYPE,myrank
          open(unit=27,file=wavefield_file,status='unknown',action='write')
        endif

        icounter = 0
        mask_ibool(:) = .false.
        do ispec = 1,nspec
          do j = 1,NGLLZ
            do i = 1,NGLLX
               iglob = ibool(i,j,ispec)
               if(.not. mask_ibool(iglob)) then
                 icounter = icounter + 1
                 mask_ibool(iglob) = .true.
                 if(use_binary_for_wavefield_dumps) then

                   if(p_sv .and. .not. imagetype_wavefield_dumps == 4) then
                     write(27,rec=icounter) sngl(vector_field_display(1,iglob)),sngl(vector_field_display(3,iglob))
                   else if(p_sv .and. imagetype_wavefield_dumps == 4) then
! by convention we use the third component of the array to store the pressure above
                     write(27,rec=icounter) sngl(vector_field_display(3,iglob))
                   else ! SH case
                     write(27,rec=icounter) sngl(vector_field_display(2,iglob))
                   endif

                 else

                   if(p_sv .and. .not. imagetype_wavefield_dumps == 4) then
                     write(27,*) sngl(vector_field_display(1,iglob)),sngl(vector_field_display(3,iglob))
                   else if(p_sv .and. imagetype_wavefield_dumps == 4) then
! by convention we use the third component of the array to store the pressure above
                     write(27,*) sngl(vector_field_display(3,iglob))
                   else ! SH case
                     write(27,*) sngl(vector_field_display(2,iglob))
                   endif

                 endif
               endif
            enddo
          enddo
        enddo

        close(27)

        if(myrank ==0) write(IOUT,*) 'Wave field dumped'

    endif  ! of display wavefield dumps at a given time step

!----  save temporary or final seismograms
! suppress seismograms if we generate traces of the run for analysis with "ParaVer", because time consuming
    if(mod(it,NSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
        call write_seismograms(sisux,sisuz,siscurl,station_name,network_name,NSTEP, &
                            nrecloc,which_proc_receiver,nrec,myrank,deltat,seismotype,st_xval,t0, &
                            NSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current,p_sv, &
                            st_zval,x_source(1),z_source(1),SU_FORMAT,save_ASCII_seismograms, &
                            save_binary_seismograms_single,save_binary_seismograms_double,subsamp_seismos)

      seismo_offset = seismo_offset + seismo_current
      seismo_current = 0

    endif  ! of display images at a given time step

  enddo ! end of the main time loop

! *********************************************************
! ************* END MAIN LOOP OVER THE TIME STEPS *********
! *********************************************************

  if(output_wavefield_dumps) deallocate(mask_ibool)

  if((SAVE_FORWARD .and. SIMULATION_TYPE==1) .or. SIMULATION_TYPE == 3) then
    if(any_acoustic) then
      close(65)
      close(66)
      close(67)
      close(68)
      close(72)
    endif
    if(any_elastic) then
      close(35)
      close(36)
      close(37)
      close(38)
      close(71)

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
      do j=1,nglob
        write(55) displ_elastic(1,j), displ_elastic(3,j), &
                  veloc_elastic(1,j), veloc_elastic(3,j), &
                  accel_elastic(1,j), accel_elastic(3,j)
      enddo
    else !SH (membrane) waves
      do j=1,nglob
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
       do j=1,nglob
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
       do j=1,nglob
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

!----  close energy file
  if(output_energy .and. myrank == 0) close(IOUT_ENERGY)

  if (OUTPUT_MODEL_VELOCITY_FILE .and. .not. any_poroelastic) then
    open(unit=1001,file='DATA/model_velocity.dat_output',status='unknown')
    if ( .NOT. assign_external_model) then
      allocate(rho_local(ngllx,ngllz,nspec)); rho_local=0.
      allocate(vp_local(ngllx,ngllz,nspec)); vp_local=0.
      allocate(vs_local(ngllx,ngllz,nspec)); vs_local=0.
      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            rho_local(i,j,ispec) = density(1,kmato(ispec))
            vp_local(i,j,ispec) = sqrt(poroelastcoef(3,1,kmato(ispec))/density(1,kmato(ispec)))
            vs_local(i,j,ispec) = sqrt(poroelastcoef(2,1,kmato(ispec))/density(1,kmato(ispec)))
            write(1001,'(I10, 5F13.4)') iglob, coord(1,iglob),coord(2,iglob),&
                                      rho_local(i,j,ispec),vp_local(i,j,ispec),vs_local(i,j,ispec)
          enddo
        enddo
      enddo
    else
      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            write(1001,'(I10,5F13.4)') iglob, coord(1,iglob),coord(2,iglob),&
                                       rhoext(i,j,ispec),vpext(i,j,ispec),vsext(i,j,ispec)
          enddo
        enddo
      enddo
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

