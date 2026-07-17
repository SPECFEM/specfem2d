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

! In this file we gather some subroutines related to the AXISYM option (some others can be found
! in gll_library.f90 and define_derivation_matrices.f90)

  subroutine  build_is_on_the_axis()

! This subroutine set the vector of logicals is_on_the_axis. is_on_the_axis(ispec)
! should be .true. if the element number ispec is located near the axis.
! The vector build_is_on_the_axis has been initialized to .false. before.

  use specfem_par, only: nspec, ispec_of_axial_elements, is_on_the_axis

  implicit none

  ! local parameters
  integer :: ispec

  is_on_the_axis(:) = .false.

  do ispec = 1,nspec
    ! checks if ispec is contained in list ispec_of_axial_element
    if (any(ispec_of_axial_elements(:) == ispec)) is_on_the_axis(ispec) = .true.
  enddo

  end subroutine build_is_on_the_axis


!
!------------------------------------------------------------------------------------------------
!

  subroutine check_compatibility_axisym()

! This subroutine check the parameters of the run and stop the program (or display a warning) if the configuration is not
! compatible with axisymetric simulations.
! _Not compatible fundamentaly :
!   If it is on the axis the source need to be in the x direction (test if sourceangle =0 /180+-epsilon)
!  Warning if the source is not on the axis (circular source)
! _Not implemented yet (or not tested):
!   poroelasticity, anisotropy, Stacey absorbing boundaries, time stepping scheme /= 1, PML rotated, adjoint
!   simulations, periodic conditions, noise tomographies

  use constants, only: PI,TWO,TINYVAL,myrank,NGLLZ,NGLLX,IMAIN,MAX_STRING_LEN,OUTPUT_FILES

  use specfem_par, only: any_poroelastic, ROTATE_PML_ACTIVATE, &
                         STACEY_ABSORBING_CONDITIONS, &
                         NSOURCES, source_type, ispec_selected_source, ADD_PERIODIC_CONDITIONS, &
                         anglesource, is_on_the_axis, ispec_is_elastic, islice_selected_source, &
                         NOISE_TOMOGRAPHY

  use specfem_par, only: ibool,coord,nspec,nglob,SAVE_MESH_FILES

  implicit none

  ! Local parameters
  integer :: isource
  integer :: ispec,i,j,ier
  ! vtk output
  double precision,dimension(:),allocatable :: xstore,zstore
  character(len=MAX_STRING_LEN) :: filename,prname

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Checking axisymmetric simulation'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (any_poroelastic) &
    call exit_MPI(myrank,'Poroelasticity is not implemented for axisymmetric simulations')
  if (ROTATE_PML_ACTIVATE) &
    call exit_MPI(myrank,'ROTATE_PML_ACTIVATE is not implemented for axisymmetric simulations')
  if (STACEY_ABSORBING_CONDITIONS) &
    call exit_MPI(myrank,'Stacey boundary conditions are not implemented for axisymmetric simulations, use PML instead')
  if (ADD_PERIODIC_CONDITIONS) &
    call exit_MPI(myrank,'Periodic conditions (ADD_PERIODIC_CONDITIONS) are not implemented for axisymmetric simulations')
  if (NOISE_TOMOGRAPHY /= 0) &
    call exit_MPI(myrank,'Axisymmetric noise tomographies are not possible yet')

  ! synchronizes MPI processes
  call synchronize_all()

  ! VTK file output for visualization
  if (SAVE_MESH_FILES) then
    if (myrank == 0) then
      write(IMAIN,*) 'axial element output files:'
      call flush_IMAIN()
    endif
    write(prname,"(a,i5.5,a)") trim(OUTPUT_FILES)//'mesh',myrank,'_'

    ! allocate xstore/zstore
    allocate(xstore(nglob), &
             zstore(nglob),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating temporary xstore/zstore arrays')
    xstore(:) = coord(1,:)
    zstore(:) = coord(2,:)

    ! element flags
    filename = trim(prname) // 'is_on_the_axis'
    call write_VTK_data_elem_l(nspec,nglob,ibool,xstore,zstore,is_on_the_axis,filename)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  written file: ',trim(filename) // '.vtk'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! free temporary arrays
    deallocate(xstore,zstore)
  endif

  ! synchronizes MPI processes
  call synchronize_all()

  ! check axial elements
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        ! check axial element
        if (is_on_the_axis(ispec)) then
          ! uses same check as in compute forces routine
          ! not first GLJ point
          if (abs(coord(1,ibool(i,j,ispec))) > TINYVAL) then
            ! check first GLJ point (i==1) is on axis and excluded here, otherwise axial element is invalid
            if (i == 1) then
              print *
              print *,'Error: AXISYM simulation has an invalid GLJ point.'
              print *,'       slice rank    : ',myrank
              print *,'       element number: ',ispec
              print *,'       GLJ point     : i/j/ispec =',i,j,ispec,' iglob =',ibool(i,j,ispec)
              print *,'       coordinate    : ',coord(:,ibool(i,j,ispec))
              print *
              print *,'       for axial elements, first GLJ point x-coordinate',coord(1,ibool(i,j,ispec)),' should be zero!'
              print *,'       That is, coordinate must be at least <= TINYVAL ',TINYVAL
              print *
              print *,'       element corner coordinates:'
              print *,'         (1,1)       : x/z = ',coord(:,ibool(1,1,ispec))
              print *,'         (NGLLX,1)   : x/z = ',coord(:,ibool(NGLLX,1,ispec))
              print *,'         (NGLLX,NGLJ): x/z = ',coord(:,ibool(NGLLX,NGLLZ,ispec))
              print *,'         (1,NGLJ)    : x/z = ',coord(:,ibool(1,NGLLZ,ispec))
              print *
              print *,'Check that your axial elements are on the symmetry axis.'
              print *,'Maybe also have a look at doc/problematic_case_that_we_exclude_for_axisymmetric.pdf'
              print *
              call exit_MPI(myrank,'Error: an axial element is invalid or rotated.')
            endif
          endif
        endif  ! is_on_the_axis
      enddo
    enddo
  enddo

  ! synchronizes MPI processes
  call synchronize_all()

  ! Check sources
  ! Loop on the sources :
  do isource = 1,NSOURCES
    if (myrank == islice_selected_source(isource)) then
      !  If the source is not an elastic force or an acoustic pressure
      if (source_type(isource) /= 1) then
        ! CMT
        ! Error (?)
        !print *,'Error: Source ',isource,' has invalid source type ',source_type(isource),' for AXISYM case!'
        !print *,'       Force/Pressure sources only (type == 1) for AXISYM simulations are implemented for now.'
        !print *,'       Please modify the source type in DATA/SOURCE accordingly...'
        !print *
        !call exit_MPI(myrank,'Axisymmetry : just elastic force or acoustic pressure sources has been tested so far')
        ! Warning (?)
        print *, '***** WARNING *****'
        print *, 'Axisymmetry: Source ',isource,' has CMT source type for AXISYM case.'
        print *, '             Only explosion monopole source (Mxx == Mzz and Mxz == 0) is possible.'
        print *
      endif

      !   If the source is on an axial element
      if (is_on_the_axis(ispec_selected_source(isource))) then
        !  ... or if the source is (at r=0) on an elastic axial element.
        if (ispec_is_elastic(ispec_selected_source(isource))) then
          if (source_type(isource) == 1) then
            ! point force source
            ! for a force angle /= 0 or 180 degree the force becomes a ring source, where the radial component cancels out.
            ! note: anglesource has been converted to radians after reading in from SOURCE file
            if (((anglesource(isource) > TINYVAL) .and. (anglesource(isource) < PI) ) &    ! ... and has a radial component.
              .or. ( (anglesource(isource) > PI) .and. (anglesource(isource) < TWO*PI))) then
              print *, '***** WARNING *****'
              print *, 'Axisymmetry: U_r(r=0)=0, Radial component of axial source will be ignored (anglesource /= 0 modulo 180)'
              print *
            endif
          endif
        endif
      else
        !   If the source is not on an axial element
        print *, '***** WARNING *****'
        print *, 'Axisymmetry: physically a source off the symmetry axis becomes a circular ring source!'
        print *
      endif
    endif
  enddo

  ! synchronizes MPI processes
  call synchronize_all()

  end subroutine check_compatibility_axisym

!
!------------------------------------------------------------------------------------------------
!

  subroutine enforce_zero_radial_displacements_on_the_axis(displ_elastic,veloc_elastic,accel_elastic)

! This subroutine enforces zero displacement, velocity and acceleration on the axis for elastic elements;
! for acoustic elements we do not need to do anything, some gradient components
! will be set to zero on the axis later in the code, when they are computed

  use constants, only: NDIM,NGLJ,NGLLZ,ZERO,CUSTOM_REAL

  use specfem_par, only: nglob, ispec_is_elastic, ibool, nelem_on_the_axis, ispec_of_axial_elements, is_on_the_axis

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: displ_elastic,veloc_elastic,accel_elastic

  ! local variables
  integer :: i_on_the_axis,ispec_axis,i,j

  do i_on_the_axis = 1,nelem_on_the_axis ! Loop on the elements on the axis
    ispec_axis = ispec_of_axial_elements(i_on_the_axis)

    ! if the element is acoustic we do not need to do anything, some gradient components
    ! will be set to zero on the axis later in the code, when they are computed

    ! if the element is elastic
    if (ispec_is_elastic(ispec_axis)) then
      do j = 1,NGLLZ ! Loop on the GLL/GLJ points
        do i = 1,NGLJ
          if (is_on_the_axis(ispec_axis) .and. i == 1) then ! If the point scanned is on the axis
            displ_elastic(1,ibool(i,j,ispec_axis)) = ZERO ! enforce the radial displacement to zero
            veloc_elastic(1,ibool(i,j,ispec_axis)) = ZERO ! enforce the radial velocity to zero
            accel_elastic(1,ibool(i,j,ispec_axis)) = ZERO ! enforce the radial acceleration to zero
          endif
        enddo
      enddo
    endif
  enddo

  end subroutine enforce_zero_radial_displacements_on_the_axis

