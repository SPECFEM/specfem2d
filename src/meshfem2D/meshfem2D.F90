
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

!========================================================================
!
!  Basic mesh generator for SPECFEM2D
!
!========================================================================

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

program meshfem2D

  use part_unstruct
  use parameter_file
  use source_file
  use interfaces_file
  implicit none

  include "constants.h"

  ! coordinates of the grid points of the mesh
  double precision, dimension(:,:), allocatable :: x,z

  ! to compute the coordinate transformation
  integer :: ioffset
  double precision :: gamma,absx,a00,a01,bot0,top0

  ! to store density and velocity model
  integer, dimension(:), allocatable :: num_material

  ! to store the position of pml element in array region_pml_external_mesh
  ! this is only useful when using pml together with external mesh
  integer, dimension(:), allocatable :: region_pml_external_mesh
  integer :: nspec_cpml

  ! interface data
  integer :: max_npoints_interface,number_of_interfaces,npoints_interface_bottom, &
    npoints_interface_top
  integer :: number_of_layers
  integer :: nz,nxread,nzread

  integer :: ilayer,ipoint_current
  integer, dimension(:), pointer :: nz_layer
  double precision, dimension(:), allocatable :: &
       xinterface_bottom,zinterface_bottom,coefs_interface_bottom, &
       xinterface_top,zinterface_top,coefs_interface_top

  integer :: nspec
  integer :: nbregion

  ! external functions
  integer, external :: num_4, num_9
  double precision, external :: value_spline

  ! variables used for storing info about the mesh and partitions
  integer, dimension(:), allocatable  :: my_interfaces
  integer, dimension(:), allocatable  :: my_nb_interfaces

  integer  :: remove_min_to_start_at_zero
  integer  :: num_node

  ! variables used for tangential detection
  integer ::  nnodes_tangential_curve
  double precision, dimension(:,:), allocatable  :: nodes_tangential_curve

#ifdef USE_SCOTCH
  integer  :: edgecut
#endif

  integer :: iproc
  integer :: ix,iz,i,j
  integer :: imaterial_number,inumelem
  integer :: i_source
  double precision :: tang1,tangN

  ! ***
  ! *** read the parameter file
  ! ***

  print *,'Reading the parameter file...'
  print *

  ! opens file Par_file
  call open_parameter_file()

  ! reads in parameters in DATA/Par_file
  call read_parameter_file()

  ! reads in mesh elements
  if ( read_external_mesh ) then
     call read_external_mesh_file(mesh_file, remove_min_to_start_at_zero, ngnod)

  else
     call read_interfaces_file(interfacesfile,max_npoints_interface, &
                                number_of_interfaces,npoints_interface_bottom, &
                                number_of_layers,nz_layer,nx,nz,nxread,nzread,ngnod, &
                                nelmnts,elmnts)
  endif

  allocate(num_material(nelmnts))
  num_material(:) = 0

  allocate(region_pml_external_mesh(nelmnts))
  region_pml_external_mesh(:) = 0

  ! assigns materials to mesh elements
  if ( read_external_mesh ) then
     call read_mat(materials_file, num_material)
     if(PML_BOUNDARY_CONDITIONS)then
      call read_pml_element(CPML_element_file, region_pml_external_mesh, nspec_cpml)
     endif
  else
     call read_regions(nbregion,nb_materials,icodemat,cp,cs, &
                      rho_s,QKappa,Qmu,aniso3,aniso4,aniso5,aniso6,aniso7,aniso8,aniso9,aniso10,aniso11, &
                      nelmnts,num_material,nxread,nzread)
  endif

  ! closes file Par_file
  call close_parameter_file()

  print *
  print *,'Parameter file successfully read... '

  ! reads in source descriptions
  call read_source_file(NSOURCES)

  ! reads in tangential detection
  if (force_normal_to_surface .or. rec_normal_to_surface) then
     open(unit=IIN,file=tangential_detection_curve_file,status='old',action='read')
     read(IIN,*) nnodes_tangential_curve
     allocate(nodes_tangential_curve(2,nnodes_tangential_curve))
     do i = 1, nnodes_tangential_curve
        read(IIN,*) nodes_tangential_curve(1,i), nodes_tangential_curve(2,i)
     enddo
     close(IIN)
  else
     nnodes_tangential_curve = 1 ! dummy values instead of 0
     allocate(nodes_tangential_curve(2,1))
  endif


  !---

  if(ngnod /= 4 .and. ngnod /= 9) stop 'ngnod different from 4 or 9!'

  print *
  print *,'The mesh contains ',nelmnts,' elements'
  print *
  print *,'Control elements have ',ngnod,' nodes'
  print *

  !---

  if ( .not. read_external_mesh ) then
     ! allocate arrays for the grid
     allocate(x(0:nx,0:nz))
     allocate(z(0:nx,0:nz))

     x(:,:) = 0.d0
     z(:,:) = 0.d0

     ! get interface data from external file
     print *,'Reading interface data from file DATA/',interfacesfile(1:len_trim(interfacesfile))
     open(unit=IIN_INTERFACES,file='DATA/'//interfacesfile,status='old')

     allocate(xinterface_bottom(max_npoints_interface))
     allocate(zinterface_bottom(max_npoints_interface))
     allocate(coefs_interface_bottom(max_npoints_interface))

     allocate(xinterface_top(max_npoints_interface))
     allocate(zinterface_top(max_npoints_interface))
     allocate(coefs_interface_top(max_npoints_interface))

     ! read number of interfaces
     call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,number_of_interfaces)

     ! read bottom interface
     call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface_bottom)

     ! loop on all the points describing this interface
     do ipoint_current = 1,npoints_interface_bottom
        call read_two_interface_points(IIN_INTERFACES,DONT_IGNORE_JUNK, &
             xinterface_bottom(ipoint_current),zinterface_bottom(ipoint_current))
     enddo

     ! loop on all the layers
     do ilayer = 1,number_of_layers

        ! read top interface
        call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface_top)

        ! loop on all the points describing this interface
        do ipoint_current = 1,npoints_interface_top
           call read_two_interface_points(IIN_INTERFACES,DONT_IGNORE_JUNK, &
                xinterface_top(ipoint_current),zinterface_top(ipoint_current))
        enddo

        ! compute the spline for the bottom interface, impose the tangent on both edges
        tang1 = (zinterface_bottom(2)-zinterface_bottom(1)) / (xinterface_bottom(2)-xinterface_bottom(1))
        tangN = (zinterface_bottom(npoints_interface_bottom)-zinterface_bottom(npoints_interface_bottom-1)) / &
             (xinterface_bottom(npoints_interface_bottom)-xinterface_bottom(npoints_interface_bottom-1))
        call spline_construction(xinterface_bottom,zinterface_bottom,npoints_interface_bottom, &
                              tang1,tangN,coefs_interface_bottom)

        ! compute the spline for the top interface, impose the tangent on both edges
        tang1 = (zinterface_top(2)-zinterface_top(1)) / (xinterface_top(2)-xinterface_top(1))
        tangN = (zinterface_top(npoints_interface_top)-zinterface_top(npoints_interface_top-1)) / &
             (xinterface_top(npoints_interface_top)-xinterface_top(npoints_interface_top-1))
        call spline_construction(xinterface_top,zinterface_top,npoints_interface_top,tang1,tangN,coefs_interface_top)

        ! check if we are in the last layer, which contains topography,
        ! and modify the position of the source accordingly if it is located exactly at the surface
        do i_source=1,NSOURCES
           if(source_surf(i_source) .and. ilayer == number_of_layers ) then
                print *, 'source ', i_source
                print *, '  target (input) z: ', zs(i_source)
                zs(i_source) = value_spline(xs(i_source),xinterface_top,zinterface_top, &
                                            coefs_interface_top,npoints_interface_top)
                print *, '  surface (actual) z: ', zs(i_source)
           endif
        enddo

        ! compute the offset of this layer in terms of number of spectral elements below along Z
        if(ilayer > 1) then
           ioffset = sum(nz_layer(1:ilayer-1))
        else
           ioffset = 0
        endif

        !--- definition of the mesh

        do ix = 0,nx

           ! evenly spaced points along X
           absx = xmin + (xmax - xmin) * dble(ix) / dble(nx)

           ! value of the bottom and top splines
           bot0 = value_spline(absx,xinterface_bottom,zinterface_bottom,coefs_interface_bottom,npoints_interface_bottom)
           top0 = value_spline(absx,xinterface_top,zinterface_top,coefs_interface_top,npoints_interface_top)

           do iz = 0,nz_layer(ilayer)

              ! linear interpolation between bottom and top
              gamma = dble(iz) / dble(nz_layer(ilayer))
              a00 = 1.d0 - gamma
              a01 = gamma

              ! coordinates of the grid points
              x(ix,iz + ioffset) = absx
              z(ix,iz + ioffset) = a00*bot0 + a01*top0

           enddo

        enddo

        ! the top interface becomes the bottom interface before switching to the next layer
        npoints_interface_bottom = npoints_interface_top
        xinterface_bottom(:) = xinterface_top(:)
        zinterface_bottom(:) = zinterface_top(:)

     enddo

     close(IIN_INTERFACES)

     nnodes = (nz+1)*(nx+1)
     allocate(nodes_coords(2,nnodes))
     if ( ngnod == 4 ) then
        do j = 0, nz
           do i = 0, nx
              num_node = num_4(i,j,nxread)
              nodes_coords(1, num_node) = x(i,j)
              nodes_coords(2, num_node) = z(i,j)

           enddo
        enddo

     else
        do j = 0, nz
           do i = 0, nx
              num_node = num_9(i,j,nxread,nzread)
              nodes_coords(1, num_node) = x(i,j)
              nodes_coords(2, num_node) = z(i,j)
           enddo
        enddo

     endif
  else
     call read_nodes_coords(nodes_coords_file)
  endif


  if ( read_external_mesh ) then
     call read_acoustic_surface(free_surface_file, num_material, &
                        ANISOTROPIC_MATERIAL, nb_materials, icodemat, phi, remove_min_to_start_at_zero)

     if ( any_abs ) then
        call read_abs_surface(absorbing_surface_file, remove_min_to_start_at_zero)
! rotate the elements that are located on the edges of the mesh if needed
! otherwise the plane wave and Bielak conditions may not be applied correctly
        if(initialfield) call rotate_mesh_for_plane_wave(ngnod)
     endif

     if ( ACOUSTIC_FORCING ) then
        call read_acoustic_forcing_surface(acoustic_forcing_surface_file, remove_min_to_start_at_zero)
! rotate the elements that are located on the edges of the mesh if needed
! otherwise the plane wave and Bielak conditions may not be applied correctly
        if(initialfield) call rotate_mesh_for_plane_wave(ngnod)
     endif

  else

     ! count the number of acoustic free-surface elements
     nelem_acoustic_surface = 0

     ! if the surface is absorbing, it cannot be free at the same time
     if(.not. abstop) then
        j = nzread
        do i = 1,nxread
           imaterial_number = num_material((j-1)*nxread+i)
           if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. phi(imaterial_number) >= 1.d0 ) then
              nelem_acoustic_surface = nelem_acoustic_surface + 1
           endif
        enddo
     endif
     if(.not. absbottom) then
        j = 1
        do i = 1,nxread
           imaterial_number = num_material((j-1)*nxread+i)
           if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. phi(imaterial_number) >= 1.d0 ) then
              nelem_acoustic_surface = nelem_acoustic_surface + 1
           endif
        enddo
     endif
     ! in the axisymmetric case if xmin == 0 the axis is a symmetry axis and thus cannot be a free surface as well
     if(.not. absleft .and. .not. (AXISYM .and. abs(xmin) < TINYVAL)) then
        i = 1
        do j = 1,nzread
           imaterial_number = num_material((j-1)*nxread+i)
           if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. phi(imaterial_number) >= 1.d0 ) then
              nelem_acoustic_surface = nelem_acoustic_surface + 1
           endif
        enddo
     endif
     if(.not. absright) then
        i = nxread
        do j = 1,nzread
           imaterial_number = num_material((j-1)*nxread+i)
           if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. phi(imaterial_number) >= 1.d0 ) then
              nelem_acoustic_surface = nelem_acoustic_surface + 1
           endif
        enddo
     endif


     allocate(acoustic_surface(4,nelem_acoustic_surface))

     nelem_acoustic_surface = 0

     if(.not. abstop) then
        j = nzread
        do i = 1,nxread
           imaterial_number = num_material((j-1)*nxread+i)
           if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. phi(imaterial_number) >=1.d0 ) then
              nelem_acoustic_surface = nelem_acoustic_surface + 1
              acoustic_surface(1,nelem_acoustic_surface) = (j-1)*nxread + (i-1)
              acoustic_surface(2,nelem_acoustic_surface) = 2
              acoustic_surface(3,nelem_acoustic_surface) = elmnts(3+ngnod*((j-1)*nxread+i-1))
              acoustic_surface(4,nelem_acoustic_surface) = elmnts(2+ngnod*((j-1)*nxread+i-1))
           endif
        enddo
     endif
     if(.not. absbottom) then
        j = 1
        do i = 1,nxread
           imaterial_number = num_material((j-1)*nxread+i)
           if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. phi(imaterial_number) >=1.d0 ) then
              nelem_acoustic_surface = nelem_acoustic_surface + 1
              acoustic_surface(1,nelem_acoustic_surface) = (j-1)*nxread + (i-1)
              acoustic_surface(2,nelem_acoustic_surface) = 2
              acoustic_surface(3,nelem_acoustic_surface) = elmnts(0+ngnod*((j-1)*nxread+i-1))
              acoustic_surface(4,nelem_acoustic_surface) = elmnts(1+ngnod*((j-1)*nxread+i-1))
           endif
        enddo
     endif
     ! in the axisymmetric case if xmin == 0 the axis is a symmetry axis and thus cannot be a free surface as well
     if(.not. absleft .and. .not. (AXISYM .and. abs(xmin) < TINYVAL)) then
        i = 1
        do j = 1,nzread
           imaterial_number = num_material((j-1)*nxread+i)
           if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. phi(imaterial_number) >=1.d0 ) then
              nelem_acoustic_surface = nelem_acoustic_surface + 1
              acoustic_surface(1,nelem_acoustic_surface) = (j-1)*nxread + (i-1)
              acoustic_surface(2,nelem_acoustic_surface) = 2
              acoustic_surface(3,nelem_acoustic_surface) = elmnts(0+ngnod*((j-1)*nxread+i-1))
              acoustic_surface(4,nelem_acoustic_surface) = elmnts(3+ngnod*((j-1)*nxread+i-1))
           endif
        enddo
     endif
     if(.not. absright) then
        i = nxread
        do j = 1,nzread
           imaterial_number = num_material((j-1)*nxread+i)
           if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. phi(imaterial_number) >=1.d0 ) then
              nelem_acoustic_surface = nelem_acoustic_surface + 1
              acoustic_surface(1,nelem_acoustic_surface) = (j-1)*nxread + (i-1)
              acoustic_surface(2,nelem_acoustic_surface) = 2
              acoustic_surface(3,nelem_acoustic_surface) = elmnts(1+ngnod*((j-1)*nxread+i-1))
              acoustic_surface(4,nelem_acoustic_surface) = elmnts(2+ngnod*((j-1)*nxread+i-1))
           endif
        enddo
     endif

     !
     !--- definition of absorbing boundaries
     !
     nelemabs = 0
     if(absbottom) nelemabs = nelemabs + nxread
     if(abstop) nelemabs = nelemabs + nxread
     if(absleft) nelemabs = nelemabs + nzread
     if(absright) nelemabs = nelemabs + nzread

     allocate(abs_surface(5,nelemabs))

     ! generate the list of absorbing elements
     if(nelemabs > 0) then
        nelemabs = 0
        do iz = 1,nzread
           do ix = 1,nxread
              inumelem = (iz-1)*nxread + ix
              if(absbottom    .and. iz == 1) then
                 nelemabs = nelemabs + 1
                 abs_surface(1,nelemabs) = inumelem-1
                 abs_surface(2,nelemabs) = 2
                 abs_surface(3,nelemabs) = elmnts(0+ngnod*(inumelem-1))
                 abs_surface(4,nelemabs) = elmnts(1+ngnod*(inumelem-1))
                 abs_surface(5,nelemabs) = IBOTTOM
              endif
              if(absright .and. ix == nxread) then
                 nelemabs = nelemabs + 1
                 abs_surface(1,nelemabs) = inumelem-1
                 abs_surface(2,nelemabs) = 2
                 abs_surface(3,nelemabs) = elmnts(1+ngnod*(inumelem-1))
                 abs_surface(4,nelemabs) = elmnts(2+ngnod*(inumelem-1))
                 abs_surface(5,nelemabs) = IRIGHT
              endif
              if(abstop   .and. iz == nzread) then
                 nelemabs = nelemabs + 1
                 abs_surface(1,nelemabs) = inumelem-1
                 abs_surface(2,nelemabs) = 2
                 abs_surface(3,nelemabs) = elmnts(3+ngnod*(inumelem-1))
                 abs_surface(4,nelemabs) = elmnts(2+ngnod*(inumelem-1))
                 abs_surface(5,nelemabs) = ITOP
              endif
              if(absleft .and. ix == 1) then
                 nelemabs = nelemabs + 1
                 abs_surface(1,nelemabs) = inumelem-1
                 abs_surface(2,nelemabs) = 2
                 abs_surface(3,nelemabs) = elmnts(0+ngnod*(inumelem-1))
                 abs_surface(4,nelemabs) = elmnts(3+ngnod*(inumelem-1))
                 abs_surface(5,nelemabs) = ILEFT
              endif
           enddo
        enddo
     endif
  endif

  if(AXISYM) then
    if(read_external_mesh) then
      call read_axial_elements_file(axial_elements_file,remove_min_to_start_at_zero)
      ! the mesh can have elements that are rotated, but for our GLJ axisymmetric implementation
      ! we assume that the r axis is along the i direction;
      ! thus this routine fixes that by rotating the elements backwards if needed to make sure
      ! this assumption is always true
      call rotate_mesh_for_axisym(ngnod)

    else ! if the mesh has been made by the internal mesher

      remove_min_to_start_at_zero = 0

      if(xmin * xmax < 0) stop 'in axisymmetric mode xmin and xmax must have the same sign, they cannot cross the symmetry axis'
      if(xmin < 0) stop 'in axisymmetric mode, case of symmetry axis on the right edge instead of left not supported yet'

      ! count the number of axial elements
      nelem_on_the_axis = 0

      ! test if the left edge is on the symmetry axis
      if(abs(xmin) < TINYVAL) then

        ! if the surface is absorbing, it cannot be axial at the same time
        if(absleft) stop 'in axisymmetric mode, the left edge cannot be both axial and absorbing'
        !all the elements on the left edge are axial because that edge is vertical and located in x = 0
        nelem_on_the_axis = nzread
        allocate(ispec_of_axial_elements(nelem_on_the_axis))
        i = 1
        do j = 1,nzread
          ispec_of_axial_elements(j) = (j-1)*nxread + (i-1) + 1
        enddo

      else ! no elements on the symmetry axis
        allocate(ispec_of_axial_elements(1))
      endif

    endif ! of if(read_external_mesh) then

  else ! of AXISYM

    nelem_on_the_axis = 0
    allocate(ispec_of_axial_elements(1))

  endif

  ! compute min and max of X and Z in the grid
  print *
  print *,'Min and max value of X in the grid = ',minval(nodes_coords(1,:)),maxval(nodes_coords(1,:))
  print *,'Min and max value of Z in the grid = ',minval(nodes_coords(2,:)),maxval(nodes_coords(2,:))
  print *


  ! ***
  ! *** create a Gnuplot file that displays the grid
  ! ***
  if (output_grid_Gnuplot .and. .not. read_external_mesh) call save_gnuplot_file(ngnod,nx,nz,x,z)


  !*****************************
  ! partitioning
  !*****************************

  ! allocates and initializes partitioning of elements
  allocate(part(0:nelmnts-1))
  part(:) = -1

  if( nproc > 1 ) then
    allocate(xadj_g(0:nelmnts))
    allocate(adjncy_g(0:MAX_NEIGHBORS*nelmnts-1))
    xadj_g(:) = 0
    adjncy_g(:) = -1
  endif

  ! construction of the graph

  ! if ngnod == 9, we work on a subarray of elements that represents the elements with four nodes (four corners) only
  ! because the adjacency of the mesh elements can be entirely determined from the knowledge of the four corners only
  if ( ngnod == 9 ) then
     allocate(elmnts_bis(0:NCORNERS*nelmnts-1))
     do i = 0, nelmnts-1
       elmnts_bis(i*NCORNERS:i*NCORNERS+NCORNERS-1) = elmnts(i*ngnod:i*ngnod+NCORNERS-1)
     enddo

     if ( nproc > 1 ) then

!! DK DK fixed problem in the previous implementation by Nicolas Le Goff:
!! DK DK (nxread+1)*(nzread+1) is OK for a regular internal mesh only, not for non structured external meshes
!! DK DK      call mesh2dual_ncommonnodes(nelmnts, (nxread+1)*(nzread+1), &
!! DK DK                                    elmnts_bis, xadj, adjncy, nnodes_elmnts, nodes_elmnts,1)
!! DK DK the subset of element corners is not renumbered therefore we must still use the nnodes computed for 9 nodes here
        ! determines maximum neighbors based on 1 common node
        call mesh2dual_ncommonnodes(elmnts_bis,1,xadj_g,adjncy_g)
     endif

  else
     if ( nproc > 1 ) then
        ! determines maximum neighbors based on 1 common node
        call mesh2dual_ncommonnodes(elmnts,1,xadj_g,adjncy_g)
     endif

  endif


  if ( nproc == 1 ) then
     part(:) = 0 ! single process has rank 0
  else

     ! number of common edges
     nb_edges = xadj_g(nelmnts)

     ! giving weight to edges and vertices. Currently not used.
!! DK DK
!! DK DK could be used to define different weights for acoustic, elastic and poroelastic elements
!! DK DK and also to define different weights for acoustic PML and elastic PML elements
!! DK DK
     call read_weights()

     ! partitioning
     select case (partitioning_method)

     case(1)

        do iproc = 0, nproc-2
           part(iproc*floor(real(nelmnts)/real(nproc)):(iproc+1)*floor(real(nelmnts)/real(nproc))-1) = iproc
        enddo
        part(floor(real(nelmnts)/real(nproc))*(nproc-1):nelmnts-1) = nproc - 1

     case(2)

!#ifdef USE_METIS
!       call Part_metis(nelmnts, xadj, adjncy, vwgt, adjwgt, nproc, nb_edges, edgecut, part, metis_options)
!#else
!       print *, 'This version of SPECFEM was not compiled with support of METIS.'
!       print *, 'Please recompile with -DUSE_METIS in order to enable use of METIS.'
!       stop
!#endif
       stop 'support for the METIS graph partitioner has been discontinued, please use SCOTCH (option 3) instead'

     case(3)

#ifdef USE_SCOTCH
        call Part_scotch(nproc, edgecut)
#else
        print *, 'This version of SPECFEM was not compiled with support of SCOTCH.'
        print *, 'Please recompile with -DUSE_SCOTCH in order to enable use of SCOTCH.'
        stop
#endif

     end select

  endif

  ! fluid-solid edges: coupled elements are transferred to the same partition
  if ( ngnod == 9 ) then
     call acoustic_elastic_repartitioning(elmnts_bis, nb_materials, phi, num_material, nproc)
  else
     call acoustic_elastic_repartitioning(elmnts, nb_materials, phi, num_material, nproc)
  endif

  ! fluid-porous edges: coupled elements are transferred to the same partition
  if ( ngnod == 9 ) then
     call acoustic_poro_repartitioning(elmnts_bis, nb_materials, phi, num_material, nproc)
  else
     call acoustic_poro_repartitioning(elmnts, nb_materials, phi, num_material, nproc)
  endif

  ! porous-solid edges: coupled elements are transferred to the same partition
  if ( ngnod == 9 ) then
     call poro_elastic_repartitioning(elmnts_bis, nb_materials, phi, num_material, nproc)
  else
     call poro_elastic_repartitioning(elmnts, nb_materials, phi, num_material, nproc)
  endif

  ! periodic edges: coupled elements are transferred to the same partition
  if(ADD_PERIODIC_CONDITIONS .and. nproc > 1) then
    if ( ngnod == 9 ) then
       call periodic_edges_repartitioning(elmnts_bis,nnodes,nodes_coords,PERIODIC_HORIZ_DIST)
    else
       call periodic_edges_repartitioning(elmnts,nnodes,nodes_coords,PERIODIC_HORIZ_DIST)
    endif
  endif

  ! local number of each element for each partition
  call Construct_glob2loc_elmnts(nproc)

  if ( ngnod == 9 ) then
    if( allocated(nnodes_elmnts) ) deallocate(nnodes_elmnts)
    if( allocated(nodes_elmnts) ) deallocate(nodes_elmnts)
    allocate(nnodes_elmnts(0:nnodes-1))
    allocate(nodes_elmnts(0:nsize*nnodes-1))
    nnodes_elmnts(:) = 0
    nodes_elmnts(:) = 0
    do i = 0, ngnod*nelmnts-1
      nodes_elmnts(elmnts(i)*nsize+nnodes_elmnts(elmnts(i))) = i/ngnod
      nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1
    enddo
  else
    if ( nproc < 2 ) then
      if( .not. allocated(nnodes_elmnts) ) allocate(nnodes_elmnts(0:nnodes-1))
      if( .not. allocated(nodes_elmnts) ) allocate(nodes_elmnts(0:nsize*nnodes-1))
      nnodes_elmnts(:) = 0
      nodes_elmnts(:) = 0
      do i = 0, ngnod*nelmnts-1
        nodes_elmnts(elmnts(i)*nsize+nnodes_elmnts(elmnts(i))) = i/ngnod
        nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1
      enddo
    endif
  endif

  ! local number of each node for each partition
  call Construct_glob2loc_nodes(nproc)

  ! construct the interfaces between partitions (used for MPI assembly)
  if ( nproc /= 1 ) then
     if ( ngnod == 9 ) then
        call Construct_interfaces(nproc, elmnts_bis, &
                                  nb_materials, phi, num_material)
     else
        call Construct_interfaces(nproc, elmnts, &
                                  nb_materials, phi, num_material)
     endif
     allocate(my_interfaces(0:ninterfaces-1))
     allocate(my_nb_interfaces(0:ninterfaces-1))
  else
     ! dummy allocation
     ninterfaces=0
     allocate(my_interfaces(0:ninterfaces-1))
     allocate(my_nb_interfaces(0:ninterfaces-1))
  endif

  ! setting absorbing boundaries by elements instead of edges
  if ( any_abs ) then
     call merge_abs_boundaries(nb_materials, phi, num_material, ngnod)
  endif

  ! setting acoustic forcing boundaries by elements instead of edges
  if ( ACOUSTIC_FORCING ) then
     call merge_acoustic_forcing_boundaries(ngnod)
  endif

  ! *** generate the databases for the solver
  call save_databases(nspec,num_material, region_pml_external_mesh, &
                      my_interfaces,my_nb_interfaces, &
                      nnodes_tangential_curve,nodes_tangential_curve,remove_min_to_start_at_zero)

  ! print position of the source
  do i_source=1,NSOURCES
     print *
     print *,'Position (x,z) of the source = ',xs(i_source),zs(i_source)
     print *
  enddo

  !--- compute position of the receivers and write the STATIONS file
  if (.not. use_existing_STATIONS) then

!! DK DK for now we cannot use both record_at_surface_same_vertical and read_external_mesh
!! DK DK because we need to know splines to define the shape of the surface of the model
    if(any(record_at_surface_same_vertical) .and. read_external_mesh) &
      stop 'for now we cannot use both record_at_surface_same_vertical and read_external_mesh'

!! DK DK if we read an external mesh, the splines are not defined for the shape of the surface and of the interfaces
!! DK DK therefore let us allocate dummy arrays just to be able to call the "save_stations_file" subroutine
    if (read_external_mesh) then
      npoints_interface_top = 1
      max_npoints_interface = 1
      allocate(xinterface_top(1))
      allocate(zinterface_top(1))
      allocate(coefs_interface_top(1))
    endif

    call save_stations_file(nreceiversets,nrec,xdeb,zdeb,xfin,zfin,record_at_surface_same_vertical, &
                            xinterface_top,zinterface_top,coefs_interface_top, &
                            npoints_interface_top,max_npoints_interface)
  endif

  print *
  if (nproc == 1) then
     print *,'This will be a serial simulation'
  else
     print *,'This will be a parallel simulation on ',nproc,' processor cores'
  endif
  print *

  if(associated(nz_layer)) deallocate(nz_layer)
  if(associated(elmnts)) deallocate(elmnts)

end program meshfem2D


