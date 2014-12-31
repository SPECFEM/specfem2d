
  program specfem2D

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
  use mpi
#endif

  use specfem_par

  implicit none

#ifdef USE_MPI
  include "precision.h"
#endif

!! DK DK uncomment this in order to force vectorization of the loops
!! DK DK using a trick that goes out of the array bounds
!! DK DK (then array bound checking cannot be used, thus for instance do NOT use -check all in Intel ifort)
! #define FORCE_VECTORIZATION


!***********************************************************************
!
!             i n i t i a l i z a t i o n    p h a s e
!
!***********************************************************************

  ! force Flush-To-Zero if available to avoid very slow Gradual Underflow trapping
  call force_ftz()

  call initialize_simulation()
  if(nproc < 1) stop 'should have nproc >= 1'

  ! starts reading in Database file
  call read_databases_init()

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
  call read_databases_sources()

  !if(AXISYM) factor = factor/(TWO*PI)   !!!!!axisym TODO verify

  ! sets source parameters
  call set_sources()

  !----  define time stepping scheme
  if(time_stepping_scheme == 1)then
    stage_time_scheme=1
  else if(time_stepping_scheme == 2)then
    stage_time_scheme=Nstages
  else if(time_stepping_scheme == 3)then
    stage_time_scheme=4
  endif

  !----  read attenuation information
  call read_databases_atten()

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
  call read_databases_coorg_elem()

  !---- allocate arrays
    allocate(shape2D(ngnod,NGLLX,NGLLZ))
    allocate(dershape2D(NDIM,ngnod,NGLLX,NGLLZ))
    allocate(shape2D_display(ngnod,pointsdisp,pointsdisp))
    allocate(dershape2D_display(NDIM,ngnod,pointsdisp,pointsdisp))
    if(AXISYM) then
      allocate(flagrange_GLJ(NGLJ,pointsdisp))
    else
      allocate(flagrange_GLJ(1,1))
    endif
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
    allocate(already_shifted_velocity(numat))
    allocate(QKappa_attenuation(numat))
    allocate(Qmu_attenuation(numat))
    allocate(kmato(nspec))
    allocate(knods(ngnod,nspec))
    allocate(ibool(NGLLX,NGLLZ,nspec))
    allocate(elastic(nspec))
    allocate(acoustic(nspec))
    allocate(gravitoacoustic(nspec))
    allocate(poroelastic(nspec))
    allocate(anisotropic(nspec))
    allocate(inv_tau_sigma_nu1(NGLLX,NGLLZ,nspec,N_SLS))
    allocate(inv_tau_sigma_nu2(NGLLX,NGLLZ,nspec,N_SLS))
    allocate(phi_nu1(NGLLX,NGLLZ,nspec,N_SLS))
    allocate(phi_nu2(NGLLX,NGLLZ,nspec,N_SLS))
    allocate(tau_epsilon_nu1(N_SLS))
    allocate(tau_epsilon_nu2(N_SLS))
    allocate(inv_tau_sigma_nu1_sent(N_SLS))
    allocate(inv_tau_sigma_nu2_sent(N_SLS))
    allocate(phi_nu1_sent(N_SLS))
    allocate(phi_nu2_sent(N_SLS))

    already_shifted_velocity(:) = .false.

  !
  !---- read the material properties
  !
  call gmat01(f0(1))
  !
  !----  read spectral macrobloc data
  !

! add support for using PML in MPI mode with external mesh
  allocate(region_CPML(nspec))
  call read_databases_mato()

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
  call initialize_simulation_domains()

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
    stop 'MPI support for standard Runge-Kutta scheme is not implemented'
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

! initialize to dummy values
! convention to indicate that Q = 9999 in that element i.e. that there is no viscoelasticity in that element
  inv_tau_sigma_nu1(:,:,:,:) = -1._CUSTOM_REAL
  phi_nu1(:,:,:,:) = -1._CUSTOM_REAL
  inv_tau_sigma_nu2(:,:,:,:) = -1._CUSTOM_REAL
  phi_nu2(:,:,:,:) = -1._CUSTOM_REAL
  Mu_nu1(:,:,:) = -1._CUSTOM_REAL
  Mu_nu2(:,:,:) = -1._CUSTOM_REAL

! define the attenuation quality factors.
! they can be different for each element.
!! DK DK if needed in the future, here the quality factor could be different for each point
  do ispec = 1,nspec

!   attenuation is not implemented in acoustic (i.e. fluid) media for now, only in viscoelastic (i.e. solid) media
    if(acoustic(ispec)) cycle

!   check that attenuation values entered by the user make sense
    if((QKappa_attenuation(kmato(ispec)) <= 9998.999d0 .and. Qmu_attenuation(kmato(ispec)) >  9998.999d0) .or. &
       (QKappa_attenuation(kmato(ispec)) >  9998.999d0 .and. Qmu_attenuation(kmato(ispec)) <= 9998.999d0)) stop &
     'need to have Qkappa and Qmu both above or both below 9999 for a given material; trick: use 9998 if you want to turn off one'

!   if no attenuation in that elastic element
    if(QKappa_attenuation(kmato(ispec)) > 9998.999d0) cycle

    call attenuation_model(QKappa_attenuation(kmato(ispec)),Qmu_attenuation(kmato(ispec)))

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

    if(ATTENUATION_VISCOELASTIC_SOLID .and. READ_VELOCITIES_AT_F0 .and. .not. assign_external_model) then
      if(anisotropic(ispec) .or. poroelastic(ispec) .or. gravitoacoustic(ispec)) &
         stop 'READ_VELOCITIES_AT_F0 only implemented for non anisotropic, non poroelastic, non gravitoacoustic materials for now'
      n = kmato(ispec)
      if(.not. already_shifted_velocity(n)) then
        rho = density(1,n)
        lambda = poroelastcoef(1,1,n)
        mu = poroelastcoef(2,1,n)
        vp = dsqrt((lambda + TWO * mu) / rho)
        vs = dsqrt(mu / rho)
        call shift_velocities_from_f0(vp,vs,rho,mu,lambda)
        poroelastcoef(1,1,n) = lambda
        poroelastcoef(2,1,n) = mu
        poroelastcoef(3,1,n) = lambda + TWO*mu
        already_shifted_velocity(n) = .true.
      endif
    endif

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
  call read_databases_ninterface()
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
   call read_databases_interfaces()

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

!---codeabs_corner(1,nelemabs) denotes whether element is on bottom-left corner of absorbing boundary or not
!---codeabs_corner(2,nelemabs) denotes whether element is on bottom-right corner of absorbing boundary or not
!---codeabs_corner(3,nelemabs) denotes whether element is on top-left corner of absorbing boundary or not
!---codeabs_corner(4,nelemabs) denotes whether element is on top-right corner of absorbing boundary or not
    allocate(codeabs_corner(4,nelemabs))
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
  call read_databases_absorbing()

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

! --- allocate arrays for acoustic forcing boundary conditions

  if(.not. ACOUSTIC_FORCING) then
    nelem_acforcing = 1
  endif

    allocate(numacforcing(nelem_acforcing))
    allocate(codeacforcing(4,nelem_acforcing))
    allocate(typeacforcing(nelem_acforcing))

    allocate(ibegin_edge1_acforcing(nelem_acforcing))
    allocate(iend_edge1_acforcing(nelem_acforcing))
    allocate(ibegin_edge3_acforcing(nelem_acforcing))
    allocate(iend_edge3_acforcing(nelem_acforcing))

    allocate(ibegin_edge4_acforcing(nelem_acforcing))
    allocate(iend_edge4_acforcing(nelem_acforcing))
    allocate(ibegin_edge2_acforcing(nelem_acforcing))
    allocate(iend_edge2_acforcing(nelem_acforcing))

    allocate(ib_left_acforcing(nelem_acforcing))
    allocate(ib_right_acforcing(nelem_acforcing))
    allocate(ib_bottom_acforcing(nelem_acforcing))
    allocate(ib_top_acforcing(nelem_acforcing))

  !
  !----  read acoustic forcing boundary data
  !
  call read_databases_acoustic_forcing()


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
  call read_databases_free_surf()
  ! resets nelem_acoustic_surface
  if( any_acoustic_edges .eqv. .false. ) nelem_acoustic_surface = 0

  ! constructs acoustic surface
  if(nelem_acoustic_surface > 0) then
    call construct_acoustic_surface ()
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

  call read_databases_coupled()

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
  call read_tangential_detection_curve()
  ! resets nnode_tangential_curve
  if( any_tangential_curve .eqv. .false. ) nnodes_tangential_curve = 0

!
!----  read axial elements data
!
  allocate(is_on_the_axis(nspec),stat=ier)
  if(ier /= 0) stop 'error: not enough memory to allocate array is_on_the_axis'
  is_on_the_axis(:) = .false.
  if(nelem_on_the_axis == 0) then
    allocate(ispec_of_axial_elements(1))
  else
    allocate(ispec_of_axial_elements(nelem_on_the_axis))
    call read_databases_axial_elements()
    call build_is_on_the_axis()
  endif
  if (myrank == 0 .and. AXISYM) then
    write(IOUT,*)
    write(IOUT,*) 'Number of elements on the axis: ',nelem_on_the_axis
  endif

!
!----  end of reading
!

! closes input Database file
 close(IIN)

!
!---- compute shape functions and their derivatives for SEM grid
!

! set up Gauss-Lobatto-Legendre derivation matrices
  call define_derivation_matrices()

  if (AXISYM) then
    ! set up Gauss-Lobatto-Jacobi derivation matrices
    call define_GLJ_derivation_matrix()
  endif

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
    call createnum_fast()
  else
    call createnum_slow()
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
    call get_global()

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
      if(AXISYM) flagrange_GLJ(j,i) = hgll(j-1,xirec,xiglj,NGLJ)
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
      allocate(gravityext(NGLLX,NGLLZ,nspec))
      allocate(Nsqext(NGLLX,NGLLZ,nspec))
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
      allocate(gravityext(1,1,1))
      allocate(Nsqext(1,1,1))
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
          if(AXISYM) then
            if (is_on_the_axis(ispec)) then
              xi = xiglj(i)
            else
              xi = xigll(i)
            endif
          else
            xi = xigll(i)
          endif
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
    call save_openDX_jacobian(nspec,npgeo,ngnod,knods,coorg,xigll,zigll,AXISYM,is_on_the_axis,xiglj)
  endif

! stop the code at the first negative element found, because such a mesh cannot be computed
  if(found_a_negative_jacobian) then

    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if(AXISYM) then
            if (is_on_the_axis(ispec)) then
              xi = xiglj(i)
            else
              xi = xigll(i)
            endif
          else
            xi = xigll(i)
          endif
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
    call plotgll()

  if (assign_external_model) then
    if(myrank == 0) write(IOUT,*) 'Assigning an external velocity and density model...'
    call read_external_model()
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

   any_gravitoacoustic_glob = any_gravitoacoustic
#ifdef USE_MPI
  call MPI_ALLREDUCE(any_gravitoacoustic, any_gravitoacoustic_glob, 1, MPI_LOGICAL, &
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

  call setup_sources_receivers()

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
           adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
           call compute_arrays_adj_source(xi_receiver(irec), gamma_receiver(irec))
           adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
         endif
       enddo
    else ! (SU_FORMAT)
        call add_adjoint_sources_SU()
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
      if(acoustic(ispec) .and. ispec == ispec_selected_rec(irec)) then
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
      if((is_proc_source(i) == 1) .and. is_on_the_axis(ispec_selected_source(i))) then
        call lagrange_any(xi_source(i),NGLJ,xiglj,hxis,hpxis)
      else
        call lagrange_any(xi_source(i),NGLLX,xigll,hxis,hpxis)
      endif
    else
      call lagrange_any(xi_source(i),NGLLX,xigll,hxis,hpxis)
    endif

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
    allocate(accel_elastic_adj_coupling2(3,nglob_elastic))

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

    ! potential, its first and second derivative, and inverse of the mass matrix for gravitoacoustic elements
    if(any_gravitoacoustic) then
      nglob_gravitoacoustic = nglob
    else
      ! allocate unused arrays with fictitious size
      nglob_gravitoacoustic = 1
    endif
    allocate(potential_gravitoacoustic(nglob_gravitoacoustic))
    allocate(potential_dot_gravitoacoustic(nglob_gravitoacoustic))
    allocate(potential_dot_dot_gravitoacoustic(nglob_gravitoacoustic))
    allocate(rmass_inverse_gravitoacoustic(nglob_gravitoacoustic))
    allocate(potential_gravito(nglob_gravitoacoustic))
    allocate(potential_dot_gravito(nglob_gravitoacoustic))
    allocate(potential_dot_dot_gravito(nglob_gravitoacoustic))
    allocate(rmass_inverse_gravito(nglob_gravitoacoustic))

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

! add support for using PML in MPI mode with external mesh
      if(read_external_mesh)then
        allocate(mask_ibool_pml(nglob))
      else
        allocate(mask_ibool_pml(1))
      endif

      call pml_init()

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
          call determin_interface_pml_interior()
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
        call define_PML_coefficients(f0(1))
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

      ! elastic PML memory variables
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
          allocate(rmemory_sfb_potential_ddot_acoustic(1,NGLLX,NGLLZ,num_fluid_solid_edges),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_sfb_potential_ddot_acoustic'
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
          if(any_acoustic .and. num_fluid_solid_edges > 0)then
            allocate(rmemory_fsb_displ_elastic_LDDRK(1,3,NGLLX,NGLLZ,num_fluid_solid_edges),stat=ier)
            if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_fsb_displ_elastic'
            allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,num_fluid_solid_edges),stat=ier)
            if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_sfb_potential_ddot_acoustic'
          endif
        else
          allocate(rmemory_displ_elastic_LDDRK(1,1,1,1,1),stat=ier)
          allocate(rmemory_dux_dx_LDDRK(1,1,1,2),stat=ier)
          allocate(rmemory_dux_dz_LDDRK(1,1,1,2),stat=ier)
          allocate(rmemory_duz_dx_LDDRK(1,1,1,2),stat=ier)
          allocate(rmemory_duz_dz_LDDRK(1,1,1,2),stat=ier)
          if(any_acoustic .and. num_fluid_solid_edges > 0)then
            allocate(rmemory_fsb_displ_elastic_LDDRK(1,3,NGLLX,NGLLZ,1),stat=ier)
            if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_fsb_displ_elastic'
            allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,1),stat=ier)
            if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_sfb_potential_ddot_acoustic'
          endif
        endif

        rmemory_displ_elastic(:,:,:,:,:) = ZERO
        rmemory_dux_dx(:,:,:,:) = ZERO
        rmemory_dux_dz(:,:,:,:) = ZERO
        rmemory_duz_dx(:,:,:,:) = ZERO
        rmemory_duz_dz(:,:,:,:) = ZERO

        if(any_acoustic .and. num_fluid_solid_edges > 0)then
          rmemory_fsb_displ_elastic(:,:,:,:,:) = ZERO
          rmemory_sfb_potential_ddot_acoustic(:,:,:,:) = ZERO
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
          if(any_acoustic .and. num_fluid_solid_edges > 0)then
            rmemory_fsb_displ_elastic_LDDRK(:,:,:,:,:) = ZERO
            rmemory_sfb_potential_ddot_acoustic_LDDRK(:,:,:,:) = ZERO
          endif
        endif

      else

        allocate(rmemory_displ_elastic(1,1,1,1,1))
        allocate(rmemory_dux_dx(1,1,1,1))
        allocate(rmemory_dux_dz(1,1,1,1))
        allocate(rmemory_duz_dx(1,1,1,1))
        allocate(rmemory_duz_dz(1,1,1,1))
        if(any_acoustic .and. num_fluid_solid_edges > 0)then
          allocate(rmemory_fsb_displ_elastic(1,3,NGLLX,NGLLZ,1))
          allocate(rmemory_sfb_potential_ddot_acoustic(1,NGLLX,NGLLZ,1))
          allocate(rmemory_fsb_displ_elastic_LDDRK(1,3,NGLLX,NGLLZ,1))
          allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,1))
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
      allocate(rmemory_sfb_potential_ddot_acoustic(1,NGLLX,NGLLZ,1))
      allocate(rmemory_fsb_displ_elastic_LDDRK(1,3,NGLLX,NGLLZ,1))
      allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,1))

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
  if(.not. allocated(rmemory_sfb_potential_ddot_acoustic)) allocate(rmemory_sfb_potential_ddot_acoustic(1,NGLLX,NGLLZ,1))
  if(.not. allocated(rmemory_fsb_displ_elastic_LDDRK)) then
    allocate(rmemory_fsb_displ_elastic_LDDRK(1,3,NGLLX,NGLLZ,1))
  endif
  if(.not. allocated(rmemory_sfb_potential_ddot_acoustic_LDDRK)) then
    allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,1))
  endif

! Test compatibility with axisymmetric formulation
  if(AXISYM) call check_compatibility_axisym()


  !
  !---- build the global mass matrix
  !
  call invert_mass_matrix_init()

#ifdef USE_MPI
  if ( nproc > 1 ) then

    ! preparing for MPI communications
    allocate(mask_ispec_inner_outer(nspec))
    mask_ispec_inner_outer(:) = .false.

    call get_MPI()

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
                            rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic,nglob_poroelastic)

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

  call invert_mass_matrix()

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
  call checkgrid()

! convert receiver angle to radians
  anglerec = anglerec * pi / 180.d0

!
!---- for color images
!

  if(output_color_image) then
    ! prepares dimension of image
    call prepare_color_image_init()

    ! allocate an array for image data
    allocate(image_color_data(NX_IMAGE_color,NZ_IMAGE_color),stat=ier); if(ier /= 0) stop 'error in an allocate statement 1'
    allocate(image_color_vp_display(NX_IMAGE_color,NZ_IMAGE_color),stat=ier); if(ier /= 0) stop 'error in an allocate statement 2'

    ! allocate an array for the grid point that corresponds to a given image data point
    allocate(iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier); if(ier /= 0) stop 'error in an allocate statement 3'
    allocate(copy_iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier); if(ier /= 0) stop 'error in an allocate statement 4'

    ! creates pixels indexing
    call prepare_color_image_pixels()

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

  potential_gravitoacoustic = 0._CUSTOM_REAL
  potential_dot_gravitoacoustic = 0._CUSTOM_REAL
  potential_dot_dot_gravitoacoustic = 0._CUSTOM_REAL
  potential_gravito = 0._CUSTOM_REAL
  potential_dot_gravito = 0._CUSTOM_REAL
  potential_dot_dot_gravito = 0._CUSTOM_REAL

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
!----- Read last frame for forward wavefield reconstruction
!

  if( ((SAVE_FORWARD .and. SIMULATION_TYPE ==1) .or. SIMULATION_TYPE == 3) .and. anyabs &
      .and. (.not. PML_BOUNDARY_CONDITIONS) ) then
    ! opens files for absorbing boundary data
    call prepare_absorb_files()
  endif

  if(anyabs .and. SIMULATION_TYPE == 3 .and. (.not. PML_BOUNDARY_CONDITIONS)) then

    ! reads in absorbing boundary data

    if(any_elastic) call prepare_absorb_elastic()

    if(any_poroelastic) call prepare_absorb_poroelastic()

    if(any_acoustic) call prepare_absorb_acoustic()

  endif ! if(anyabs .and. SIMULATION_TYPE == 3)



!
!----- Allocate sensitivity kernel arrays
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
      call prepare_initialfield()
    endif

    if( over_critical_angle ) then

      allocate(left_bound(nelemabs*NGLLX))
      allocate(right_bound(nelemabs*NGLLX))
      allocate(bot_bound(nelemabs*NGLLZ))

      call prepare_initialfield_paco()

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
      call paco_beyond_critical(anglesource(1),&
                                f0(1),QKappa_attenuation(1),source_type(1),left_bound(1:count_left),&
                                right_bound(1:count_right),bot_bound(1:count_bottom), &
                                count_left,count_right,count_bottom,x_source(1))

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
    call prepare_source_time_function()
  else
    ! uses an initialfield
    ! dummy allocation
    allocate(source_time_function(1,1,1))
  endif

! acoustic forcing edge detection
! the elements forming an edge are already known (computed in meshfem2D),
! the common nodes forming the edge are computed here
  if(ACOUSTIC_FORCING) then

    if (myrank == 0) then
      print *
      print *,'Acoustic forcing simulation'
      print *
      print *,'Beginning of acoustic forcing edge detection'
    endif

! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,IBOTTOM) = ipoin1D
      jvalue(ipoin1D,IBOTTOM) = NGLLZ
      jvalue_inverse(ipoin1D,IBOTTOM) = NGLLZ

      ivalue(ipoin1D,IRIGHT) = 1
      ivalue_inverse(ipoin1D,IRIGHT) = 1
      jvalue(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,IRIGHT) = ipoin1D

      ivalue(ipoin1D,ITOP) = ipoin1D
      ivalue_inverse(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,ITOP) = 1
      jvalue_inverse(ipoin1D,ITOP) = 1

      ivalue(ipoin1D,ILEFT) = NGLLX
      ivalue_inverse(ipoin1D,ILEFT) = NGLLX
      jvalue(ipoin1D,ILEFT) = ipoin1D
      jvalue_inverse(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1

    enddo

  endif ! if(ACOUSTIC_FORCING)


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
    call read_parameters_noise()

  endif ! NOISE_TOMOGRAPHY /= 0


  if (NOISE_TOMOGRAPHY == 1) then
    call compute_source_array_noise()

    !write out coordinates of mesh
    open(unit=504,file='OUTPUT_FILES/mesh_spec',status='unknown',action='write')
      do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            write(504,'(1pe11.3,1pe11.3,2i3,i7)') coord(1,iglob), coord(2,iglob), i, j, ispec
         enddo
        enddo
      enddo
    close(504)

    open(unit=504,file='OUTPUT_FILES/mesh_glob',status='unknown',action='write')
      do iglob = 1, nglob
        write(504,'(1pe11.3,1pe11.3,i7)') coord(1,iglob), coord(2,iglob), iglob
      enddo
    close(504)

    !write out spatial distribution of noise sources
    call create_mask_noise()
    open(unit=504,file='OUTPUT_FILES/mask_noise',status='unknown',action='write')
      do iglob = 1, nglob
            write(504,'(1pe11.3,1pe11.3,1pe11.3)') coord(1,iglob), coord(2,iglob), mask_noise(iglob)
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
                coord(1,iglob), coord(2,iglob), rhoext(i,j,ispec), vpext(i,j,ispec), vsext(i,j,ispec)
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
    call create_mask_noise()

  else if (NOISE_TOMOGRAPHY == 3) then

    if (output_wavefields_noise) then
      call create_mask_noise()

      !prepare array that will hold wavefield snapshots
      noise_output_ncol = 5
      allocate(noise_output_array(noise_output_ncol,nglob))
      allocate(noise_output_rhokl(nglob))
    endif

  endif

!>NOISE_TOMOGRAPHY

  ! prepares image background
  if(output_color_image) then
    call prepare_color_image_vp()
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
    timeval = (it-1)*deltat

    do i_stage=1, stage_time_scheme

      call update_displacement_precondition_newmark()
      if (AXISYM) then
        do ispec=1,nspec
          if (elastic(ispec) .and. is_on_the_axis(ispec)) then
            do j = 1,NGLLZ
              do i = 1,NGLJ
                if (abs(coord(1,ibool(i,j,ispec))) < TINYVAL) then
                  displ_elastic(1,ibool(i,j,ispec))=ZERO
                endif
              enddo
            enddo
          endif
        enddo
      endif

    if(any_acoustic) then
      ! free surface for an acoustic medium
      if ( nelem_acoustic_surface > 0 ) then
        call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                          potential_acoustic)

        if(SIMULATION_TYPE == 3) then ! Adjoint calculation
          call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                            b_potential_acoustic)
        endif
      endif

! *********************************************************
! ************* compute forces for the acoustic elements
! *********************************************************

      call compute_forces_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic, &
               potential_acoustic,potential_acoustic_old,.false.,PML_BOUNDARY_CONDITIONS)

      if( SIMULATION_TYPE == 3 ) then

       if(PML_BOUNDARY_CONDITIONS)then
          do ispec = 1,nspec
            do i = 1, NGLLX
              do j = 1, NGLLZ
                if(acoustic(ispec) .and. is_pml(ispec))then
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

        call compute_forces_acoustic(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
               b_potential_acoustic,b_potential_acoustic_old,.true.,.false.)

       if(PML_BOUNDARY_CONDITIONS)then
          do ispec = 1,nspec
            do i = 1, NGLLX
              do j = 1, NGLLZ
                if(acoustic(ispec) .and. is_pml(ispec))then
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


    ! *********************************************************
    ! ************* add acoustic forcing at a rigid boundary
    ! *********************************************************
    if(ACOUSTIC_FORCING) then
      ! loop on all the forced edges

     do inum = 1,nelem_acforcing

        ispec = numacforcing(inum)

        !--- left acoustic forcing boundary
        if(codeacforcing(IEDGE4,inum)) then

           i = 1

           do j = 1,NGLLZ

              ! acoustic spectral element
              if(acoustic(ispec)) then
                 iglob = ibool(i,j,ispec)

                 xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
                 zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xgamma**2 + zgamma**2)
                 nx = - zgamma / jacobian1D
                 nz = + xgamma / jacobian1D

                 weight = jacobian1D * wzgll(j)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then
              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()
            endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

              endif  !end of acoustic
           enddo

        endif  !  end of left acoustic forcing boundary

        !--- right acoustic forcing boundary
        if(codeacforcing(IEDGE2,inum)) then

           i = NGLLX

           do j = 1,NGLLZ

              ! acoustic spectral element
              if(acoustic(ispec)) then
                 iglob = ibool(i,j,ispec)

                 xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
                 zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xgamma**2 + zgamma**2)
                 nx = + zgamma / jacobian1D
                 nz = - xgamma / jacobian1D

                 weight = jacobian1D * wzgll(j)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then

              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()

            endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

              endif  !end of acoustic
           enddo

        endif  !  end of right acoustic forcing boundary

        !--- bottom acoustic forcing boundary
        if(codeacforcing(IEDGE1,inum)) then

           j = 1

           do i = 1,NGLLX

              ! acoustic spectral element
              if(acoustic(ispec)) then
                 iglob = ibool(i,j,ispec)

                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 nx = + zxi / jacobian1D
                 nz = - xxi / jacobian1D

                 weight = jacobian1D * wxgll(i)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then
              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else
            call acoustic_forcing_boundary()

            endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

              endif  !end of acoustic
           enddo

        endif  !  end of bottom acoustic forcing boundary

        !--- top acoustic forcing boundary
        if(codeacforcing(IEDGE3,inum)) then

           j = NGLLZ

           do i = 1,NGLLX

              ! acoustic spectral element
              if(acoustic(ispec)) then
                 iglob = ibool(i,j,ispec)

                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 nx = - zxi / jacobian1D
                 nz = + xxi / jacobian1D

                 weight = jacobian1D * wxgll(i)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then

              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()

            endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

              endif  !end of acoustic
           enddo

        endif  !  end of top acoustic forcing boundary

     enddo

     endif ! of if ACOUSTIC_FORCING

    endif ! end of test if any acoustic element

! *********************************************************
! ************* add coupling with the elastic side
! *********************************************************

    if(coupled_acoustic_elastic) then

      if(SIMULATION_TYPE == 1)then
        call compute_coupling_acoustic_el(displ_elastic,displ_elastic_old,potential_dot_dot_acoustic, PML_BOUNDARY_CONDITIONS)
      endif

      if(SIMULATION_TYPE == 3)then

        accel_elastic_adj_coupling2 = - accel_elastic_adj_coupling

        call compute_coupling_acoustic_el(accel_elastic_adj_coupling2,displ_elastic_old,potential_dot_dot_acoustic,&
                                          PML_BOUNDARY_CONDITIONS)

        call compute_coupling_acoustic_el(b_displ_elastic,b_displ_elastic_old,b_potential_dot_dot_acoustic,.false.)

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
          if (is_proc_source(i_source) == 1 .and. acoustic(ispec_selected_source(i_source))) then

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
              if (acoustic(ispec_selected_rec(irec))) then
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
      call assemble_MPI_vector_ac(potential_dot_dot_acoustic)

     if(time_stepping_scheme == 2)then
      if(i_stage==1 .and. it == 1 .and. (.not. initialfield))then
       potential_dot_acoustic_temp = potential_dot_acoustic
       call assemble_MPI_vector_ac(potential_dot_acoustic)
      endif
     endif

      if ( SIMULATION_TYPE == 3) then
        call assemble_MPI_vector_ac(b_potential_dot_dot_acoustic)

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
                                   + deltat*potential_dot_acoustic

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
                                        potential_acoustic)

        if(SIMULATION_TYPE == 3) then
          call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                          b_potential_acoustic)
        endif

      endif

      ! update the potential field (use a new array here) for coupling terms
      potential_acoustic_adj_coupling = potential_acoustic &
                          + deltat*potential_dot_acoustic &
                          + deltatsquareover2*potential_dot_dot_acoustic

    endif ! of if(any_acoustic)


! *********************************************************
! ************* main solver for the gravitoacoustic elements
! *********************************************************
! only SIMULATION_TYPE == 1, time_stepping_scheme == 1, and no PML or STACEY yet
! NO MIX OF ACOUSTIC AND GRAVITOACOUTIC ELEMENTS
! NO COUPLING TO ELASTIC AND POROELASTIC SIDES
! *********************************************************
!-----------------------------------------
    if ((any_gravitoacoustic)) then

      if(time_stepping_scheme==1)then
      ! Newmark time scheme
!! DK DK this should be vectorized
      potential_gravitoacoustic = potential_gravitoacoustic &
                          + deltat*potential_dot_gravitoacoustic &
                          + deltatsquareover2*potential_dot_dot_gravitoacoustic
      potential_dot_gravitoacoustic = potential_dot_gravitoacoustic &
                              + deltatover2*potential_dot_dot_gravitoacoustic
      potential_gravito = potential_gravito &
                          + deltat*potential_dot_gravito &
                          + deltatsquareover2*potential_dot_dot_gravito
      potential_dot_gravito = potential_dot_gravito &
                              + deltatover2*potential_dot_dot_gravito
      else
      stop 'Only time_stepping_scheme=1 for gravitoacoustic'
      endif
      potential_dot_dot_gravitoacoustic = ZERO
      potential_dot_dot_gravito = ZERO

! Impose displacements from boundary forcing here
! because at this step the displacement (potentials) values
! are already equal to value at n+1
! equivalent to free surface condition
! the contour integral u.n is computed after compute_forces_gravitoacoustic
! *********************************************************
! ** impose displacement from acoustic forcing at a rigid boundary
! ** force potential_dot_dot_gravito by displacement
! *********************************************************
    if(ACOUSTIC_FORCING) then

      ! loop on all the forced edges

     do inum = 1,nelem_acforcing

        ispec = numacforcing(inum)

        !--- left acoustic forcing boundary
        if(codeacforcing(IEDGE4,inum)) then

           i = 1

           do j = 1,NGLLZ

              ! acoustic spectral element
              if(gravitoacoustic(ispec)) then
                 iglob = ibool(i,j,ispec)
                 xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
                 zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xgamma**2 + zgamma**2)
                 nx = - zgamma / jacobian1D
                 nz = + xgamma / jacobian1D

                 weight = jacobian1D * wzgll(j)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then

              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()

            endif

! compute displacement at this point
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

          ! if external density model
          if(assign_external_model)then
            if(CUSTOM_REAL == SIZE_REAL) then
              rhol = sngl(rhoext(i,j,ispec))
              gravityl = sngl(gravityext(i,j,ispec))
            else
              rhol = rhoext(i,j,ispec)
              gravityl = gravityext(i,j,ispec)
            endif
          endif

! impose potential_gravito in order to have z displacement equal to forced
! value
          iglob = ibool(i,j,ispec)
          displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
          potential_gravito(iglob) = (rhol*displ_n - &
          (tempx1l*xizl + tempx2l*gammazl)*nz - (tempx1l*xixl + tempx2l*gammaxl)*nx)/ &
          (0.0 - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz
          potential_dot_dot_gravito(iglob) = potential_dot_dot_gravito(iglob) - rhol*weight*displ_n

              endif  !end of gravitoacoustic
           enddo

        endif  !  end of left acoustic forcing boundary

        !--- right acoustic forcing boundary
        if(codeacforcing(IEDGE2,inum)) then

           i = NGLLX

           do j = 1,NGLLZ

              ! acoustic spectral element
              if(gravitoacoustic(ispec)) then
                 iglob = ibool(i,j,ispec)
                 xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
                 zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xgamma**2 + zgamma**2)
                 nx = + zgamma / jacobian1D
                 nz = - xgamma / jacobian1D

                 weight = jacobian1D * wzgll(j)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then

              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()

            endif

! compute displacement at this point
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

          ! if external density model
          if(assign_external_model)then
            if(CUSTOM_REAL == SIZE_REAL) then
              rhol = sngl(rhoext(i,j,ispec))
              gravityl = sngl(gravityext(i,j,ispec))
            else
              rhol = rhoext(i,j,ispec)
              gravityl = gravityext(i,j,ispec)
            endif
          endif

! impose potential_gravito in order to have z displacement equal to forced
! value
          iglob = ibool(i,j,ispec)
          displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
          potential_gravito(iglob) = (rhol*displ_n - &
          (tempx1l*xizl + tempx2l*gammazl)*nz - (tempx1l*xixl + tempx2l*gammaxl)*nx)/ &
          (0.0-gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz
          potential_dot_dot_gravito(iglob) = potential_dot_dot_gravito(iglob) - rhol*weight*displ_n

              endif  !end of gravitoacoustic
           enddo

        endif  !  end of right acoustic forcing boundary

        !--- bottom acoustic forcing boundary
        if(codeacforcing(IEDGE1,inum)) then

           j = 1

           do i = 1,NGLLX

              ! acoustic spectral element
              if(gravitoacoustic(ispec)) then
                 iglob = ibool(i,j,ispec)
                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 nx = + zxi / jacobian1D
                 nz = - xxi / jacobian1D

                 weight = jacobian1D * wxgll(i)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then

              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()

            endif

! compute displacement at this point
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

          ! if external density model
          if(assign_external_model)then
            if(CUSTOM_REAL == SIZE_REAL) then
              rhol = sngl(rhoext(i,j,ispec))
              gravityl = sngl(gravityext(i,j,ispec))
            else
              rhol = rhoext(i,j,ispec)
              gravityl = gravityext(i,j,ispec)
            endif
          endif

! impose potential_gravito in order to have z displacement equal to forced
! value
          iglob = ibool(i,j,ispec)
          displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
          potential_gravito(iglob) = (rhol*displ_n - &
          (tempx1l*xizl + tempx2l*gammazl)*nz - (tempx1l*xixl + tempx2l*gammaxl)*nx)/ &
          (0.0 - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz
          potential_dot_dot_gravito(iglob) = potential_dot_dot_gravito(iglob) - rhol*weight*displ_n

              endif  !end of gravitoacoustic
           enddo

        endif  !  end of bottom acoustic forcing boundary

        !--- top acoustic forcing boundary
        if(codeacforcing(IEDGE3,inum)) then

           j = NGLLZ

           do i = 1,NGLLX

              ! acoustic spectral element
              if(gravitoacoustic(ispec)) then
                 iglob = ibool(i,j,ispec)
                 xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                 zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
                 jacobian1D = sqrt(xxi**2 + zxi**2)
                 nx = - zxi / jacobian1D
                 nz = + xxi / jacobian1D

                 weight = jacobian1D * wxgll(i)

          ! define displacement components which will force the boundary

            if(PML_BOUNDARY_CONDITIONS) then

              if(is_PML(ispec)) then
              displ_x = 0
              displ_z = 0
              else
              call acoustic_forcing_boundary()
              endif

            else

            call acoustic_forcing_boundary()

            endif

! compute z displacement at this point
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

          ! if external density model
          if(assign_external_model)then
            if(CUSTOM_REAL == SIZE_REAL) then
              rhol = sngl(rhoext(i,j,ispec))
              gravityl = sngl(gravityext(i,j,ispec))
            else
              rhol = rhoext(i,j,ispec)
              gravityl = gravityext(i,j,ispec)
              Nsql = Nsqext(i,j,ispec)
            endif
          endif

! impose potential_gravito in order to have z displacement equal to forced
! value on the boundary
!!!! Passe deux fois sur le même iglob
!!!! Mais vrai pour tous les points partagés entre deux elements
          iglob = ibool(i,j,ispec)
          displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
          potential_gravito(iglob) = (rhol*displ_n - &
          (tempx1l*xizl + tempx2l*gammazl)*nz - (tempx1l*xixl + tempx2l*gammaxl)*nx)/ &
          (0.0 - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz
          potential_dot_dot_gravito(iglob) = potential_dot_dot_gravito(iglob) - rhol*weight*displ_n

              endif  !end of gravitoacoustic
           enddo

!       write(*,*) 'ispec detection =',ispec
!       if ((ispec==2000).and.(mod(it,100)==0)) then
       if ((ispec==800).and.(mod(it,100)==0)) then
!       if ((ispec==800)) then
       iglobzero=iglob
       write(*,*) ispec,it,Nsql,rhol,displ_n, &
       maxval(potential_dot_dot_gravito),potential_dot_dot_gravito(iglob), &
       maxval(potential_gravitoacoustic),potential_gravitoacoustic(iglob), &
       maxval(potential_gravito),potential_gravito(iglob)
       endif

        endif  !  end of top acoustic forcing boundary

     enddo

     endif ! end ACOUSTIC_FORCING !

      ! free surface for a gravitoacoustic medium
      !!! to be coded !!!
!      if ( nelem_acoustic_surface > 0 ) then
!        call enforce_acoustic_free_surface(potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
!                                          potential_gravitoacoustic)

!        if(SIMULATION_TYPE == 3) then ! Adjoint calculation
!          call enforce_acoustic_free_surface(b_potential_dot_dot_gravitoacoustic,b_potential_dot_gravitoacoustic, &
!                                            b_potential_gravitoacoustic)
!        endif
!      endif

! *********************************************************
! ************* compute forces for the gravitoacoustic elements
! *********************************************************

      call compute_forces_gravitoacoustic(potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
               potential_gravitoacoustic, potential_dot_dot_gravito, &
               potential_gravito,.false.,PML_BOUNDARY_CONDITIONS)

       if ((mod(it,100)==0)) then
         iglob=iglobzero
         write(*,*) it,Nsql,gravityl, &
         maxval(potential_dot_dot_gravito),potential_dot_dot_gravito(iglob), &
         maxval(potential_dot_dot_gravitoacoustic),potential_dot_dot_gravitoacoustic(iglob)
       endif

    endif ! end of test if any gravitoacoustic element

! *********************************************************
! ************* add coupling with the elastic side
! *********************************************************

! *********************************************************
! ************* add coupling with the poroelastic side
! *********************************************************

! ************************************************************************************
! ************************************ add force source
! ************************************************************************************

! assembling potential_dot_dot for gravitoacoustic elements
!#ifdef USE_MPI
!    if ( nproc > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
!      call assemble_MPI_vector_ac(potential_dot_dot_gravitoacoustic)
!
!    endif
!
!#endif

! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

    if((any_gravitoacoustic)) then
      if(time_stepping_scheme == 1)then
!! DK DK this should be vectorized

      potential_dot_dot_gravitoacoustic = potential_dot_dot_gravitoacoustic * rmass_inverse_gravitoacoustic
      potential_dot_gravitoacoustic = potential_dot_gravitoacoustic + deltatover2*potential_dot_dot_gravitoacoustic

!! line below already done in compute_forces_gravitoacoustic, because necessary
!! for the computation of potential_dot_dot_gravitoacoustic
!      potential_dot_dot_gravito = potential_dot_dot_gravito * rmass_inverse_gravito
      potential_dot_gravito = potential_dot_gravito + deltatover2*potential_dot_dot_gravito
      else
        stop 'Only time_stepping_scheme = 1 implemented for gravitoacoustic case'
      endif

      ! free surface for an acoustic medium
!      if ( nelem_acoustic_surface > 0 ) then
!        call enforce_acoustic_free_surface(potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
!                                        potential_gravitoacoustic)
!
!        if(SIMULATION_TYPE == 3) then
!          call enforce_acoustic_free_surface(b_potential_dot_dot_gravitoacoustic,b_potential_dot_gravitoacoustic, &
!                                          b_potential_gravitoacoustic)
!        endif
!
!      endif
!
      ! update the potential field (use a new array here) for coupling terms
!      potential_gravitoacoustic_adj_coupling = potential_gravitoacoustic &
!                          + deltat*potential_dot_gravitoacoustic &
!                          + deltatsquareover2*potential_dot_dot_gravitoacoustic

    endif ! of if(any_gravitoacoustic)



! *********************************************************
! ************* main solver for the elastic elements
! *********************************************************

    if(any_elastic) then

      call compute_forces_viscoelastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old,x_source(1),z_source(1), &
               f0(1),v0x_left(1,it),v0z_left(1,it),v0x_right(1,it),v0z_right(1,it),v0x_bot(1,it),v0z_bot(1,it), &
               t0x_left(1,it),t0z_left(1,it),t0x_right(1,it),t0z_right(1,it),t0x_bot(1,it),t0z_bot(1,it), &
               count_left,count_right,count_bottom,PML_BOUNDARY_CONDITIONS,.false.)

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

      call compute_forces_viscoelastic(b_accel_elastic,b_veloc_elastic,b_displ_elastic,&
               b_displ_elastic_old,x_source(1),z_source(1),f0(1),&
               v0x_left(1,it),v0z_left(1,it),v0x_right(1,it),v0z_right(1,it),v0x_bot(1,it),v0z_bot(1,it), &
               t0x_left(1,it),t0z_left(1,it),t0x_right(1,it),t0z_right(1,it),t0x_bot(1,it),t0z_bot(1,it), &
               count_left,count_right,count_bottom,.false.,.true.)

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


      call compute_forces_viscoelastic_pre_kernel()

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

    if(coupled_acoustic_elastic) call compute_coupling_viscoelastic_ac()

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
          call  add_point_source_noise()

        else if (NOISE_TOMOGRAPHY == 2) then
          call add_surface_movie_noise(accel_elastic)

        else if (NOISE_TOMOGRAPHY == 3) then
          if (.not. save_everywhere) then
            call add_surface_movie_noise(b_accel_elastic)
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
       call assemble_MPI_vector_el(veloc_elastic)
       endif
    endif
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,ier)

    if (nproc > 1 .and. any_elastic .and. ninterface_elastic > 0) then
      call assemble_MPI_vector_el(accel_elastic)
    endif


    if (nproc > 1 .and. any_elastic .and. ninterface_elastic > 0 .and. SIMULATION_TYPE == 3) then
      call assemble_MPI_vector_el(b_accel_elastic)
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

      call compute_forces_poro_solid(f0(1))

      call compute_forces_poro_fluid(f0(1))


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
      call assemble_MPI_vector_po(accels_poroelastic,accelw_poroelastic)
    endif

    if (nproc > 1 .and. any_poroelastic .and. ninterface_poroelastic > 0 .and. SIMULATION_TYPE == 3) then
      call assemble_MPI_vector_po(b_accels_poroelastic,b_accelw_poroelastic)
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
                                            b_potential_acoustic)
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
    call save_surface_movie_noise()

  else if ( NOISE_TOMOGRAPHY == 2 .and. save_everywhere ) then
    call save_surface_movie_noise()

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

      call compute_energy()

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
      call check_stability()

    endif

!---- loop on all the receivers to compute and store the seismograms
    if(mod(it-1,subsamp_seismos) == 0) then
      call write_seismograms()
    endif

!----- writing the kernels
    ! kernels output
    if(SIMULATION_TYPE == 3) then

      if(any_acoustic) then

        do ispec = 1, nspec
          if(acoustic(ispec)) then
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
          if(acoustic(ispec)) then
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
! write kernel files
!

      if(SIMULATION_TYPE == 3 .and. it == NSTEP) then
          call save_adjoint_kernels()
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

          call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                          potential_gravito,displ_elastic,displs_poroelastic)

          call plotpost()

        else if(imagetype_postscript == 2 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'drawing velocity vector as small arrows...'

          call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                          potential_dot_gravito,veloc_elastic,velocs_poroelastic)

          call plotpost()

        else if(imagetype_postscript == 3 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'drawing acceleration vector as small arrows...'

          call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
                          potential_dot_dot_gravito,accel_elastic,accels_poroelastic)

          call plotpost()

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
          call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                          potential_gravito,displ_elastic,displs_poroelastic)

        else if(imagetype_JPEG >= 4 .and. imagetype_JPEG <= 6) then

          if (myrank == 0) write(IOUT,*) 'drawing scalar image of part of the velocity vector...'
          call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                          potential_dot_gravito,veloc_elastic,velocs_poroelastic)

        else if(imagetype_JPEG >= 7 .and. imagetype_JPEG <= 9) then

          if (myrank == 0) write(IOUT,*) 'drawing scalar image of part of the acceleration vector...'
          call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
                          potential_dot_dot_gravito,accel_elastic,accels_poroelastic)

        else if(imagetype_JPEG >= 11 .and. imagetype_JPEG <= 13) then
! allocation for normalized representation in JPEG image
! for an atmosphere model

          if (myrank == 0) write(IOUT,*) 'drawing scalar image of part of normalized displacement vector...'

          call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                          potential_gravito,displ_elastic,displs_poroelastic)


          do ispec = 1,nspec
            do j = 1,NGLLZ
              do i = 1,NGLLX
                iglob = ibool(i,j,ispec)
                vector_field_display(1,iglob) = sqrt(rhoext(i,j,ispec)) * vector_field_display(1,iglob)
                vector_field_display(2,iglob) = sqrt(rhoext(i,j,ispec)) * vector_field_display(2,iglob)
                vector_field_display(3,iglob) = sqrt(rhoext(i,j,ispec)) * vector_field_display(3,iglob)
              enddo
            enddo
          enddo

        else if(imagetype_JPEG >= 14 .and. imagetype_JPEG <= 16) then
! allocation for normalized representation in JPEG image
! for an atmosphere model
          call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                          potential_dot_gravito,veloc_elastic,velocs_poroelastic)

          do ispec = 1,nspec
            do j = 1,NGLLZ
              do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            vector_field_display(1,iglob) = sqrt(rhoext(i,j,ispec)) * vector_field_display(1,iglob)
            vector_field_display(2,iglob) = sqrt(rhoext(i,j,ispec)) * vector_field_display(2,iglob)
            vector_field_display(3,iglob) = sqrt(rhoext(i,j,ispec)) * vector_field_display(3,iglob)
              enddo
            enddo
          enddo

        else if(imagetype_JPEG == 10 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'drawing image of pressure field...'
          call compute_pressure_whole_medium()

        else if(imagetype_JPEG == 10 .and. .not. p_sv) then
          call exit_MPI('cannot draw pressure field for SH (membrane) waves')

        else
          call exit_MPI('wrong type for JPEG snapshots')
        endif

!! DK DK quick hack to remove the PMLs from JPEG images if needed: set the vector field to zero there
        if(PML_BOUNDARY_CONDITIONS .and. REMOVE_PMLS_FROM_JPEG_IMAGES) then
          do ispec = 1,nspec
            if(is_PML(ispec)) then
              do j = 1,NGLLZ
                do i = 1,NGLLX
                  iglob = ibool(i,j,ispec)
                  vector_field_display(1,iglob) = 0.d0
                  vector_field_display(2,iglob) = 0.d0
                  vector_field_display(3,iglob) = 0.d0
                enddo
              enddo
            endif
          enddo
        endif
!! DK DK quick hack to remove the PMLs from JPEG images if needed

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
              if(imagetype_JPEG == 1 .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7 .or. imagetype_JPEG == 11 &
                                     .or. imagetype_JPEG == 14) then
                image_color_data(i,j) = vector_field_display(1,iglob_image_color(i,j))  ! draw the X component of the vector

              else if(imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8 .or. imagetype_JPEG == 12 &
                                          .or. imagetype_JPEG == 15) then
                image_color_data(i,j) = vector_field_display(3,iglob_image_color(i,j))  ! draw the Z component of the vector
              else if(imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9 .or. imagetype_JPEG == 13 &
                                          .or. imagetype_JPEG == 16) then
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

              if(imagetype_JPEG == 1 .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7 .or. imagetype_JPEG == 11 &
                                     .or. imagetype_JPEG == 14) then
                  data_pixel_send(k) = vector_field_display(1,iglob_image_color(i,j))  ! draw the X component of the vector

              else if(imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8 .or. imagetype_JPEG == 12 &
                                          .or. imagetype_JPEG == 15) then
                  data_pixel_send(k) = vector_field_display(3,iglob_image_color(i,j))  ! draw the Z component of the vector

              else if(imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9 .or. imagetype_JPEG == 13 &
                                          .or. imagetype_JPEG == 16) then
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
          call create_color_image()
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
          call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                            potential_gravito,displ_elastic,displs_poroelastic)
        else if(imagetype_wavefield_dumps == 2) then

          if (myrank == 0) write(IOUT,*) 'dumping the velocity vector...'
          call compute_vector_whole_medium(potential_dot_acoustic,potential_gravitoacoustic, &
                            potential_gravito,veloc_elastic,velocs_poroelastic)

        else if(imagetype_wavefield_dumps == 3) then

          if (myrank == 0) write(IOUT,*) 'dumping the acceleration vector...'
          call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_gravitoacoustic, &
                            potential_gravito,accel_elastic,accels_poroelastic)

        else if(imagetype_wavefield_dumps == 4 .and. p_sv) then

          if (myrank == 0) write(IOUT,*) 'dumping the pressure field...'
          call compute_pressure_whole_medium()

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

      call write_seismograms_to_file(x_source(1),z_source(1))

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

