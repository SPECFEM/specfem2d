!========================================================================
!
!                            S P E C F E M 2 D
!                            -----------------
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
!
! Please find in the header of specfem2D.F90 further code informations.
!
! ************** PROGRAM STARTS HERE **************

  program meshfem2D

  use constants, only: IMAIN,ISTANDARD_OUTPUT,TINYVAL,OUTPUT_FILES

  use shared_parameters
  use part_unstruct_par
  use compute_elements_load_par

  implicit none

  include 'version.fh'

  integer :: nspec_cpml
  integer :: i,j,ier,num_elmnt
  logical :: BROADCAST_AFTER_READ

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
#ifdef WITH_MPI
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
  ! (secondary processes can sit idle)
  if (myrank == 0) then
    ! ***
    ! *** read the parameter file
    ! ***
    ! reads in parameters in DATA/Par_file
    BROADCAST_AFTER_READ = .false.
    call read_parameter_file(.true.,BROADCAST_AFTER_READ)

    ! reads in additional files for mesh elements
    if (read_external_mesh) then
      ! external meshing
      ! user output
      write(IMAIN,*)
      write(IMAIN,*) 'Mesh from external meshing:'
      call flush_IMAIN()

      ! reads in mesh
      call read_external_mesh_file(mesh_file, remove_min_to_start_at_zero, NGNOD)

      ! reads in material defined in external file
      call read_external_materials_file(materials_file)

    else
      ! internal meshing
      ! user output
      write(IMAIN,*)
      write(IMAIN,*) 'Mesh from internal meshing:'
      call flush_IMAIN()

      allocate(elmnts(0:NGNOD*nelmnts-1),stat=ier)
      if (ier /= 0) call stop_the_code('Error allocating array elmnts')
      elmnts(:) = 0

      ! stores mesh point indices in array 'elmnts'
      if (NGNOD == 4) then
        num_elmnt = 0
        do j = 1, nzread
           do i = 1, nxread
              elmnts(num_elmnt*NGNOD)   = (j-1)*(nxread+1) + (i-1)
              elmnts(num_elmnt*NGNOD+1) = (j-1)*(nxread+1) + (i-1) + 1
              elmnts(num_elmnt*NGNOD+2) = j*(nxread+1) + (i-1) + 1
              elmnts(num_elmnt*NGNOD+3) = j*(nxread+1) + (i-1)
              num_elmnt = num_elmnt + 1
           enddo
        enddo
      else if (NGNOD == 9) then
        num_elmnt = 0
        do j = 1, nzread
           do i = 1, nxread
              elmnts(num_elmnt*NGNOD)   = (j-1)*(nxread+1) + (i-1)
              elmnts(num_elmnt*NGNOD+1) = (j-1)*(nxread+1) + (i-1) + 1
              elmnts(num_elmnt*NGNOD+2) = j*(nxread+1) + (i-1) + 1
              elmnts(num_elmnt*NGNOD+3) = j*(nxread+1) + (i-1)
              elmnts(num_elmnt*NGNOD+4) = (nxread+1)*(nzread+1) + (j-1)*nxread + (i-1)
              elmnts(num_elmnt*NGNOD+5) = (nxread+1)*(nzread+1) + nxread*(nzread+1) + (j-1)*(nxread*2+1) + (i-1)*2 + 2
              elmnts(num_elmnt*NGNOD+6) = (nxread+1)*(nzread+1) + j*nxread + (i-1)
              elmnts(num_elmnt*NGNOD+7) = (nxread+1)*(nzread+1) + nxread*(nzread+1) + (j-1)*(nxread*2+1) + (i-1)*2
              elmnts(num_elmnt*NGNOD+8) = (nxread+1)*(nzread+1) + nxread*(nzread+1) + (j-1)*(nxread*2+1) + (i-1)*2 + 1
              num_elmnt = num_elmnt + 1
           enddo
        enddo
      else
        call stop_the_code('NGNOD must be either 4 or 9')
      endif

      ! user output
      write(IMAIN,*) '  Total number of spectral elements         = ',nelmnts
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! PML mesh elements
    ! user output
    write(IMAIN,*) 'PML mesh elements:'
    call flush_IMAIN()

    allocate(region_pml_external_mesh(nelmnts),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating array region_pml_external_mesh')
    region_pml_external_mesh(:) = 0

    allocate(is_pml(0:nelmnts-1),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating array is_pml')
    is_pml(:) = .false.
    nspec_cpml = 0

    if (PML_BOUNDARY_CONDITIONS) then
      if (read_external_mesh) then
        call read_external_pml_element(absorbing_cpml_file, region_pml_external_mesh, nspec_cpml)
      else
        ! no need to read in PML values.
        ! the internal mesher will assign PML elements in routine pml_init() in the solver.
        nspec_cpml = 0
        write(IMAIN,*) '  using internal mesh, PML elements will be determined in solver run...'
        write(IMAIN,*)
      endif
    else
      ! no PML condition
      ! user output
      write(IMAIN,*) '  Total number of PML elements = ',nspec_cpml
      write(IMAIN,*)
    endif

    ! user output
    write(IMAIN,*) 'The mesh contains ',nelmnts,' elements'
    write(IMAIN,*)
    write(IMAIN,*) 'Control elements have ',NGNOD,' nodes'
    write(IMAIN,*)
    call flush_IMAIN()

    ! reads in source descriptions
    ! daniel todo: for now needed for ADD_PERTURBATION_AROUND_SOURCE_ONLY and if source_surf(..) is set to .true.
    !              will try to avoid dependency on sources in mesher in future...
    call read_source_file(NSOURCES,BROADCAST_AFTER_READ)

    ! reads in tangential detection
    call read_mesh_tangential_curve_file()

    ! reads in node coordinates
    ! user output
    write(IMAIN,*) 'Node coordinates:'
    call flush_IMAIN()

    if (read_external_mesh) then
      call read_external_mesh_nodes_coords(nodes_coords_file)
    else
      ! reads interfaces and sets node coordinates
      call read_mesh_nodes_coords_from_interfaces()
    endif

    ! user output
    write(IMAIN,*) 'Mesh surfaces:'
    call flush_IMAIN()

    if (read_external_mesh) then
      call read_external_acoustic_surface(free_surface_file, num_material, &
                                          nbmodels, icodemat, phi_read, remove_min_to_start_at_zero)

      if (any_abs) then
        call read_external_abs_surface(absorbing_surface_file, remove_min_to_start_at_zero)

        ! rotate the elements that are located on the edges of the mesh if needed
        ! otherwise the plane wave and Bielak conditions may not be applied correctly
        if (initialfield) call rotate_mesh_for_plane_wave(NGNOD)
      endif

      if (ACOUSTIC_FORCING) then
        call read_external_acoustic_forcing_surface(acoustic_forcing_surface_file, remove_min_to_start_at_zero)

        ! rotate the elements that are located on the edges of the mesh if needed
        ! otherwise the plane wave and Bielak conditions may not be applied correctly
        call rotate_mesh_for_acoustic_forcing(NGNOD)
      endif

    else
      ! determines acoustic free surface
      call determine_acoustic_surface()

      ! determines absorbing boundary elements
      call determine_abs_surface()
    endif

    ! axi-symmetric mesh
    if (AXISYM) then
      ! user output
      write(IMAIN,*) 'Axisymmetric mesh:'
      call flush_IMAIN()

      if (read_external_mesh) then
        ! external meshing
        call read_external_axial_elements_file(axial_elements_file,remove_min_to_start_at_zero)
        ! the mesh can have elements that are rotated, but for our GLJ axisymmetric implementation
        ! we assume that the r axis is along the i direction;
        ! thus this routine fixes that by rotating the elements backwards if needed to make sure
        ! this assumption is always true
        call rotate_mesh_for_axisym(NGNOD)

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
    if (output_grid_Gnuplot .and. .not. read_external_mesh) &
      call save_gnuplot_file(NGNOD,nx_elem_internal,nz_elem_internal,grid_point_x,grid_point_z)

    ! partitioning
    ! user output
    write(IMAIN,*) 'Mesh partitioning:'
    call flush_IMAIN()

    call decompose_mesh()

    ! setting absorbing boundaries by elements instead of edges
    if (any_abs) then
      ! user output
      write(IMAIN,*) 'Absorbing boundaries:'
      call flush_IMAIN()

      call merge_abs_boundaries(nbmodels, phi_read, num_material, NGNOD)
    endif

    ! setting acoustic forcing boundaries by elements instead of edges
    if (ACOUSTIC_FORCING) then
      ! user output
      write(IMAIN,*) 'Acoustic forcing boundaries:'
      call flush_IMAIN()

      call merge_acoustic_forcing_boundaries(NGNOD)
    endif

    ! generate the databases for the solver
    ! user output
    write(IMAIN,*) 'Saving databases:'
    call flush_IMAIN()

    call save_databases()

    !--- compute position of the receivers and write the STATIONS file
    if (.not. use_existing_STATIONS) then
      ! user output
      write(IMAIN,*) 'creating STATIONS file...'
      call flush_IMAIN()

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

      ! daniel todo: move to solver?
      !         note that the STATIONS file will be written out here and then read in the solver to locate them.
      !         in case we have record_at_surface_same_vertical(.) set to .true., the stations will be
      !         placed at the surface, using splines for the top interface.
      !
      !         in future, we might want to move this to the solver, to make the solver more independant
      !         and work in a similar way like the 3D versions.
      !         still, this would require to have these top splines in the solver...
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

  ! secondary processes wait
  call synchronize_all()

  ! MPI finish
  call finalize_mpi()

  end program meshfem2D


