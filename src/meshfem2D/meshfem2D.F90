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
! journal = {communications in Computational Physics},
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
! number= 2,
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
! number= 1,
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
! journal = {communications in Computational Physics},
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

  use constants, only: IMAIN,ISTANDARD_OUTPUT,TINYVAL,OUTPUT_FILES

  use shared_parameters
  use part_unstruct_par
  use source_file_par
  use compute_elements_load_par

  implicit none

  include 'version.fh'

  integer :: nspec_cpml
  integer :: i,j,i_source,ier,num_elmnt

  ! MPI initialization
  call init_mpi()
  call world_size(NPROC)
  call world_rank(myrank)

  ! open main output file, only written to by process 0
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) then
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'output_meshfem2D.txt',status='unknown',iostat=ier)
    if (ier /= 0) then
      print *,'Error could not open output file :',trim(OUTPUT_FILES)//'output_meshfem2D.txt'
      call stop_the_code('Error opening output file')
    endif
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
#ifdef USE_MPI
    write(IMAIN,*) '**********************************************'
    write(IMAIN,*) '*** Specfem 2-D Mesher - MPI version       ***'
    write(IMAIN,*) '**********************************************'
#else
    write(IMAIN,*) '**********************************************'
    write(IMAIN,*) '*** Specfem 2-D Mesher - serial version    ***'
    write(IMAIN,*) '**********************************************'
#endif
    write(IMAIN,*)
    write(IMAIN,*) 'Running Git version of the code corresponding to ', git_commit_version
    write(IMAIN,*) 'dating ', git_date_version
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! initializes
  remove_min_to_start_at_zero = 0

  ! mesher works only for single process
  ! (slave processes can sit idle)
  if (myrank == 0) then
    ! ***
    ! *** read the parameter file
    ! ***
    write(IMAIN,*) 'Reading the parameter file...'
    write(IMAIN,*)

    ! reads in parameters in DATA/Par_file
    call read_parameter_file(1,.false.)

    ! reads in additional files for mesh elements
    if (read_external_mesh) then
      ! external meshing
      ! user output
      write(IMAIN,*)
      write(IMAIN,*) 'Mesh from external meshing:'

      ! reads in mesh
      call read_external_mesh_file(mesh_file, remove_min_to_start_at_zero, ngnod)

      ! reads in material defined in external file
      call read_external_materials_file(materials_file)

    else
      ! internal meshing
      allocate(elmnts(0:ngnod*nelmnts-1),stat=ier)
      if (ier /= 0) call stop_the_code('Error allocating array elmnts')

      ! stores mesh point indices in array 'elmnts'
      if (ngnod == 4) then
        num_elmnt = 0
        do j = 1, nzread
           do i = 1, nxread
              elmnts(num_elmnt*ngnod)   = (j-1)*(nxread+1) + (i-1)
              elmnts(num_elmnt*ngnod+1) = (j-1)*(nxread+1) + (i-1) + 1
              elmnts(num_elmnt*ngnod+2) = j*(nxread+1) + (i-1) + 1
              elmnts(num_elmnt*ngnod+3) = j*(nxread+1) + (i-1)
              num_elmnt = num_elmnt + 1
           enddo
        enddo
      else if (ngnod == 9) then
        num_elmnt = 0
        do j = 1, nzread
           do i = 1, nxread
              elmnts(num_elmnt*ngnod)   = (j-1)*(nxread+1) + (i-1)
              elmnts(num_elmnt*ngnod+1) = (j-1)*(nxread+1) + (i-1) + 1
              elmnts(num_elmnt*ngnod+2) = j*(nxread+1) + (i-1) + 1
              elmnts(num_elmnt*ngnod+3) = j*(nxread+1) + (i-1)
              elmnts(num_elmnt*ngnod+4) = (nxread+1)*(nzread+1) + (j-1)*nxread + (i-1)
              elmnts(num_elmnt*ngnod+5) = (nxread+1)*(nzread+1) + nxread*(nzread+1) + (j-1)*(nxread*2+1) + (i-1)*2 + 2
              elmnts(num_elmnt*ngnod+6) = (nxread+1)*(nzread+1) + j*nxread + (i-1)
              elmnts(num_elmnt*ngnod+7) = (nxread+1)*(nzread+1) + nxread*(nzread+1) + (j-1)*(nxread*2+1) + (i-1)*2
              elmnts(num_elmnt*ngnod+8) = (nxread+1)*(nzread+1) + nxread*(nzread+1) + (j-1)*(nxread*2+1) + (i-1)*2 + 1
              num_elmnt = num_elmnt + 1
           enddo
        enddo
      else
        call stop_the_code('ngnod must be either 4 or 9')
      endif

      ! user output
      write(IMAIN,*) 'Total number of spectral elements         = ',nelmnts
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! PML mesh elements
    allocate(region_pml_external_mesh(nelmnts),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating array region_pml_external_mesh')
    region_pml_external_mesh(:) = 0

    allocate(is_pml(0:nelmnts-1),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating array is_pml')
    is_pml(:) = .false.

    if (read_external_mesh) then
      if (PML_BOUNDARY_CONDITIONS) then
        call read_external_pml_element(absorbing_cpml_file, region_pml_external_mesh, nspec_cpml)
      else
        nspec_cpml = 0
      endif
    endif

    ! user output
    write(IMAIN,*)
    write(IMAIN,*) 'Parameter file successfully read '
    write(IMAIN,*)
    write(IMAIN,*) 'The mesh contains ',nelmnts,' elements'
    write(IMAIN,*)
    write(IMAIN,*) 'Control elements have ',ngnod,' nodes'
    write(IMAIN,*)

    ! reads in source descriptions
    call read_source_file(NSOURCES)

    ! reads in tangential detection
    call read_mesh_tangential_curve_file()

    ! reads in node coordinates
    if (read_external_mesh) then
      call read_external_mesh_nodes_coords(nodes_coords_file)
    else
      ! reads interfaces and sets node coordinates
      call read_mesh_nodes_coords_from_interfaces()
    endif

    if (read_external_mesh) then
      call read_external_acoustic_surface(free_surface_file, num_material, &
                                          nbmodels, icodemat, phi_read, remove_min_to_start_at_zero)

      if (any_abs) then
        call read_external_abs_surface(absorbing_surface_file, remove_min_to_start_at_zero)

        ! rotate the elements that are located on the edges of the mesh if needed
        ! otherwise the plane wave and Bielak conditions may not be applied correctly
        if (initialfield) call rotate_mesh_for_plane_wave(ngnod)
      endif

      if (ACOUSTIC_FORCING) then
        call read_external_acoustic_forcing_surface(acoustic_forcing_surface_file, remove_min_to_start_at_zero)

        ! rotate the elements that are located on the edges of the mesh if needed
        ! otherwise the plane wave and Bielak conditions may not be applied correctly
        call rotate_mesh_for_acoustic_forcing(ngnod)
      endif

    else
      ! determines acoustic free surface
      call determine_acoustic_surface()

      ! determines absorbing boundary elements
      call determine_abs_surface()
    endif

    if (AXISYM) then
      if (read_external_mesh) then
        ! external meshing
        call read_external_axial_elements_file(axial_elements_file,remove_min_to_start_at_zero)
        ! the mesh can have elements that are rotated, but for our GLJ axisymmetric implementation
        ! we assume that the r axis is along the i direction;
        ! thus this routine fixes that by rotating the elements backwards if needed to make sure
        ! this assumption is always true
        call rotate_mesh_for_axisym(ngnod)

      else
        ! internal meshing
        ! if the mesh has been made by the internal mesher
        if (xmin_param * xmax_param < 0) &
          call stop_the_code('in axisymmetric mode xmin and xmax must have the same sign, they cannot cross the symmetry axis')
        if (xmin_param < 0) &
          call stop_the_code('in axisymmetric mode, case of symmetry axis on the right edge instead of left not supported yet')

        ! count the number of axial elements
        nelem_on_the_axis = 0

        ! test if the left edge is on the symmetry axis
        if (abs(xmin_param) < TINYVAL) then

          ! if the surface is absorbing, it cannot be axial at the same time
          if (absorbleft) call stop_the_code('in axisymmetric mode, the left edge cannot be both axial and absorbing')

          !all the elements on the left edge are axial because that edge is vertical and located in x = 0
          nelem_on_the_axis = nzread
          allocate(ispec_of_axial_elements(nelem_on_the_axis),stat=ier)
          if (ier /= 0) call stop_the_code('Error allocating array ispec_of_axial_elements')

          i = 1
          do j = 1,nzread
            ispec_of_axial_elements(j) = (j-1)*nxread + (i-1) + 1
          enddo

        else
          ! no elements on the symmetry axis
          allocate(ispec_of_axial_elements(nelem_on_the_axis),stat=ier)
          if (ier /= 0) call stop_the_code('Error allocating array ispec_of_axial_elements')
        endif

      endif ! of if (read_external_mesh) then
      ! user output
      write(IMAIN,*)
      write(IMAIN,*) 'Axial elements: ',nelem_on_the_axis
      write(IMAIN,*)
    else
      ! .not. AXISYM
      ! dummy allocation
      nelem_on_the_axis = 0
      allocate(ispec_of_axial_elements(nelem_on_the_axis),stat=ier)
      if (ier /= 0) call stop_the_code('Error allocating array ispec_of_axial_elements')
    endif

    ! compute min and max of X and Z in the grid
    write(IMAIN,*)
    write(IMAIN,*) 'Mesh dimensions: '
    write(IMAIN,*) '  Min and max value of X in the grid = ',minval(nodes_coords(1,:)),maxval(nodes_coords(1,:))
    write(IMAIN,*) '  Min and max value of Z in the grid = ',minval(nodes_coords(2,:)),maxval(nodes_coords(2,:))
    write(IMAIN,*)

    ! create a Gnuplot file that displays the grid
    if (output_grid_Gnuplot .and. .not. read_external_mesh) call save_gnuplot_file(ngnod,nx,nz,grid_point_x,grid_point_z)

    ! partitioning
    call decompose_mesh()

    ! setting absorbing boundaries by elements instead of edges
    if (any_abs) then
      call merge_abs_boundaries(nbmodels, phi_read, num_material, ngnod)
    endif

    ! setting acoustic forcing boundaries by elements instead of edges
    if (ACOUSTIC_FORCING) then
      call merge_acoustic_forcing_boundaries(ngnod)
    endif

    ! generate the databases for the solver
    call save_databases()

    ! print position of the source
    do i_source= 1,NSOURCES
      write(IMAIN,*)
      write(IMAIN,*) 'Position (x,z) of the source = ',xs(i_source),zs(i_source)
      write(IMAIN,*)
    enddo

    !--- compute position of the receivers and write the STATIONS file
    if (.not. use_existing_STATIONS) then

!! DK DK for now we cannot use both record_at_surface_same_vertical and read_external_mesh
!! DK DK because we need to know splines to define the shape of the surface of the model
      if (any(record_at_surface_same_vertical) .and. read_external_mesh) &
        call stop_the_code('for now we cannot use both record_at_surface_same_vertical and read_external_mesh')

!! DK DK if we read an external mesh, the splines are not defined for the shape of the surface and of the interfaces
!! DK DK therefore let us allocate dummy arrays just to be able to call the "save_stations_file" subroutine
      if (read_external_mesh) then
        npoints_interface_top = 1
        max_npoints_interface = 1
        allocate(xinterface_top(1))
        allocate(zinterface_top(1))
        allocate(coefs_interface_top(1))
      endif

      call save_stations_file(nreceiversets,nrec_line,xdeb,zdeb,xfin,zfin,record_at_surface_same_vertical, &
                              xinterface_top,zinterface_top,coefs_interface_top, &
                              npoints_interface_top,max_npoints_interface)
    endif

    ! user output
    if (NPROC == 1) then
      write(IMAIN,*)
      write(IMAIN,*) 'This will be a serial simulation'
      write(IMAIN,*)
    else
      write(IMAIN,*)
      write(IMAIN,*) 'This will be a parallel simulation on ',NPROC,' processor cores'
      write(IMAIN,*)
    endif

    ! frees memory
    if (allocated(nz_layer)) deallocate(nz_layer)
    if (allocated(elmnts)) deallocate(elmnts)

  ! mesher works only for single process
  endif ! myrank == 0

  ! close main output file
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

  ! slave processes wait
  call synchronize_all()

  ! MPI finish
  call finalize_mpi()

  end program meshfem2D


