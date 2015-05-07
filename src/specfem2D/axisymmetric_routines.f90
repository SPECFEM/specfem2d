
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

! On this file we gather some subroutines related to the AXISYM option (some others can be found
! in gll_library.f90 and define_derivation_matrices.f90)
!axisym!axisym!axisym!axisym!axisym!axisym!axisym!axisym!axisym!axisym!axisym

  module qsort_c_module

! Recursive Fortran 95 quicksort routines (used by enforce_zero_radial_displacements_on_the_axis)
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing
! Made F conformant by Walt Brainerd

    implicit none
    public :: qsortC
    private :: partition

    contains

    recursive subroutine qsortC(A)
      double precision, intent(in out), dimension(:) :: A
      integer :: iq

      if(size(A) > 1) then
         call partition(A, iq)
         call qsortC(A(:iq-1))
         call qsortC(A(iq:))
      endif
    end subroutine QsortC

    subroutine partition(A, marker)
      double precision, intent(in out), dimension(:) :: A
      integer, intent(out) :: marker
      integer :: i, j
      double precision :: temp
      double precision :: x      ! pivot point
      x = A(1)
      i = 0
      j = size(A) + 1

      do
         j = j-1
         do
            if (A(j) <= x) exit
            j = j-1
         enddo
         i = i+1
         do
            if (A(i) >= x) exit
            i = i+1
         enddo
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
         else if (i == j) then
            marker = i+1
            return
         else
            marker = i
            return
         endif
      enddo

    end subroutine partition

  end module qsort_c_module

subroutine  build_is_on_the_axis()
  ! This subroutine set the vector of logicals is_on_the_axis. is_on_the_axis(ispec)
  ! should be .true. if the element number ispec is located near the axis.
  ! The vector build_is_on_the_axis has been initialized to .false. before.

    use specfem_par, only: nspec, ispec_of_axial_elements, is_on_the_axis

    implicit none

    ! include "constants.h"

    ! local parameters
    integer :: ispec

    is_on_the_axis(:) = .false. ! Normally this has been done before but we can never be too careful.

    do ispec = 1,nspec
      if (any(ispec_of_axial_elements == ispec)) is_on_the_axis(ispec) = .true.
    enddo

  ! do ispec = 1,nspec
  !   counter = 0
  !   do j = 1,NGLLZ,NGLLZ-1    ! just check the corners
  !     do i = 1,NGLLX,NGLLX-1  ! just check the corners
  !       iglob = ibool(i,j,ispec)
  !       if (abs(coord(1,iglob)) < TINYVAL) counter = counter + 1
  !     enddo
  !   enddo
  !
  !if ( counter == 1 ) call exit_MPI('error: element with a single point on the symmetry axis found')
  !if ( counter > 2 ) call exit_MPI('error: element with more than two points on the symmetry axis found; this should never happen')
  ! TODO see how to check that with our new implementation
  !   if ( counter == 2 ) then
  !     is_on_the_axis(ispec)=.true.
  ! DK DK and AB AB: for now make sure that the left edge is the one that is on the axis;
  ! DK DK and AB AB: we could remove this constraint but we purposely do not do it for now
  !     iglob1 = ibool(1,1,ispec)
  !     iglob2 = ibool(1,NGLLZ,ispec)
  !     if ( abs(coord(1,iglob1)) > TINYVAL .or. abs(coord(1,iglob2)) > TINYVAL ) then
  !        call exit_MPI('DK DK and AB AB: axial element found is rotated; this is currently not implemented')
  ! TODO see how to check that with our new implementation
  !     endif
  !   endif
  ! enddo

  end subroutine build_is_on_the_axis



  subroutine check_compatibility_axisym()
  ! This subroutine check the parameters of the run and stop the program (or display a warning) if the configuration is not
  ! compatible with axisymetric simulations.
  ! _Not compatible fundamentaly :
  !   If it is on the axis the source need to be in the x direction (test if sourceangle =0 /180+-epsilon)
  !  Warning if the source is not on the axis (circular source)
  ! _Not implemented yet (or not tested):
  !   poroelasticity, anisotropy, attenuation, Stacey absorbing boundaries, time stepping scheme /= 1, PML rotated, adjoint
  !   simulations, periodic conditions, noise tomographies, source and receivers on the elements along the axis but not
  !   exactly on the axis.

    use specfem_par, only: any_poroelastic, anisotropic,ROTATE_PML_ACTIVATE, &
                           STACEY_BOUNDARY_CONDITIONS, SIMULATION_TYPE, SAVE_FORWARD,time_stepping_scheme, &
                           NOISE_TOMOGRAPHY, NSOURCES, source_type, ispec_selected_source, xi_source, ADD_PERIODIC_CONDITIONS, &
                           anglesource, nrec, ispec_selected_rec, xi_receiver, is_on_the_axis, elastic, myrank, &
                           is_proc_source, which_proc_receiver

    implicit none

    include "constants.h"

    ! Local parameters
    integer :: irec, isource

    if ( any_poroelastic ) &
      call exit_MPI('Poroelasticity is presently not implemented for axisymmetric simulations')
    !if ( ATTENUATION_VISCOELASTIC_SOLID) &
    !  call exit_MPI('Attenuation (ATTENUATION_VISCOELASTIC_SOLID) is presently not implemented for axisymmetric simulations')
    if ( any(anisotropic(:) .eqv. .true.) ) &
      call exit_MPI('Anisotropy is presently not implemented for axisymmetric simulations')
    if ( ROTATE_PML_ACTIVATE ) &
      call exit_MPI('ROTATE_PML_ACTIVATE is presently not implemented for axisymmetric simulations')
    if ( STACEY_BOUNDARY_CONDITIONS ) &
      call exit_MPI('Stacey boundary conditions are presently not implemented for axisymmetric simulations -> use PML instead')
    if ( SIMULATION_TYPE /= 1 ) &
      call exit_MPI('Just axisymmetric FORWARD simulations are possible so far')
    if ( SAVE_FORWARD ) &
      call exit_MPI('SAVE_FORWARD has presently not been tested with axisymmetric simulations')
    if ( time_stepping_scheme /= 1 ) &
      call exit_MPI('Just Newmark scheme is presently possible for axisymmetric simulation')
    if ( ADD_PERIODIC_CONDITIONS ) &
      call exit_MPI('Periodic conditions (ADD_PERIODIC_CONDITIONS) are presently not implemented for axisymmetric simulations')
    if ( NOISE_TOMOGRAPHY /= 0 ) &
      call exit_MPI('Axisymmetric noise tomographies are not possible yet')

    ! Check receivers
    do irec = 1,nrec                                                               ! Loop on the receivers :
      if ( myrank == which_proc_receiver(irec) ) then
        if ( is_on_the_axis(ispec_selected_rec(irec)) ) then                         !   If the receiver is on an axial element ...
          if ( (xi_receiver(irec) > -1.d0) .and. (xi_receiver(irec) < 1.d0 ) ) then  !   ... but not exactly on the edges
            call exit_MPI('Axisymmetry : a receiver has been set on the 1st element but not on the edges, this is not possible yet')
          endif
        endif
      endif
    enddo

    ! Check sources
    do isource = 1,NSOURCES                                      ! Loop on the sources :
      if ( is_proc_source(isource) == 1 ) then
        if ( source_type(isource) /= 1 ) then                       !  If the source is not an elastic force or an acoustic pressure
          call exit_MPI('Axisymmetry : just elastic force or acoustic pressure sources are implemented so far')
        endif
        if ( is_on_the_axis(ispec_selected_source(isource)) ) then  !   If the source is on an axial element
          if ( (xi_source(isource) > -1.d0)  .and. (xi_source(isource) < 1.d0 )) then  !   ... but not exactly on the edges
            call exit_MPI('Axisymmetry : a source has been set on the first element but not on the axis, this is not possible yet')
          endif
          if ( elastic(ispec_selected_source(isource)) ) then       !  ... or if the source is (at r=0) on an elastic axial element.
            if ( ((anglesource(isource) > TINYVAL) .and. (anglesource(isource) < PI) ) &    ! ... and has a radial component.
               .or. ( (anglesource(isource) > PI) .and. (anglesource(isource) < TWO*PI)) ) then
               print *, '***** WARNING *****'
               print *, 'Axisymmetry : U_r(r=0)=0, Radial component of axial source will be ignored (anglesource /= 0 modulo 180)'
               print *, ''
            endif
          endif
        else                                                      !   If the source is not on an axial element
           print *, '***** WARNING *****'
           print *, 'Axisymmetry : physically a non axial source is a circular source!'
           print *, ''
        endif
      endif
    enddo

  end subroutine check_compatibility_axisym

  subroutine enforce_zero_radial_displacements_on_the_axis()
    ! This subroutine enforce zero displacement on the axis
    ! _For elastic elements it is straight forward : displ_elastic(1,ibool(i,j,ispec))=ZERO
    ! _For acoustic elements rho*u = grad(potential) hence the two elements near the axis need to share the same values
    !  of potential

    use specfem_par, only: acoustic, elastic, coord, ibool, nelem_on_the_axis, ispec_of_axial_elements, is_on_the_axis, &
                           potential_acoustic, displ_elastic!, potential_dot_acoustic, potential_dot_dot_acoustic, veloc_elastic, &
                           !accel_elastic
    use qsort_c_module

    implicit none

    include "constants.h"

    ! Local :
    integer i_on_the_axis,ispec_axis,i,j,i_1,i_2
    double precision, dimension(NGLJ) :: coord_line, coord_line_sorted
    integer, dimension(NGLJ) :: indices

    do i_on_the_axis = 1,nelem_on_the_axis ! Loop on the elements on the axis
      ispec_axis = ispec_of_axial_elements(i_on_the_axis)
      if ( acoustic(ispec_axis) ) then ! If the element is acoustic
        ! TODO : For the moment we don't do anything
        do j = 1,NGLLZ ! For each depth
          ! For each depth we have NGLJ points (say 4) : *  *  *  *
          ! We have to know which on is the first one, and which one is at its side
          do i = 1,NGLJ
            coord_line(i) = coord(1,ibool(i,j,ispec_axis)) ! Ex : coord_line = (8.12 3.68 15.42 0.0)
          enddo
          coord_line_sorted = coord_line
          !print *, "coord_line :",coord_line
          !print *, "coord_line_sorted :",coord_line_sorted
          call qsortC(coord_line_sorted) ! Ex : coord_line_sorted = (0.0,3.68,8.12,15.42)
          do i = 1,NGLJ
           indices(i) = minloc(abs(coord_line - coord_line_sorted(i)), 1) ! Ex : Indices = (3,2,4,1)
           if ( indices(i) == 1 ) i_1 = i ! Ex : i_1 = 4 -> index of the first point (the one on the axis)
           if ( indices(i) == 2 ) i_2 = i ! Ex : i_2 = 2 -> index of the one at its side
          enddo
          !print *, "coord_line :",coord_line
          !print *, "coord_line_sorted :",coord_line_sorted
          !print *, "indices :", indices
          potential_acoustic(ibool(i_1,j,ispec_axis)) = potential_acoustic(ibool(i_2,j,ispec_axis))
          !potential_dot_acoustic(ibool(i_2,j,ispec_axis)) = potential_dot_acoustic(ibool(i_1,j,ispec_axis))
          !potential_dot_dot_acoustic(ibool(i_2,j,ispec_axis)) = potential_dot_dot_acoustic(ibool(i_1,j,ispec_axis))
        enddo

      else if ( elastic(ispec_axis) ) then ! Else if the element is elastic

        do j = 1,NGLLZ ! Loop on the GLL/GLJ points
          do i = 1,NGLJ
            if( is_on_the_axis(ispec_axis) .and. i == 1 ) then ! If the point scanned is on the axis
              displ_elastic(1,ibool(i,j,ispec_axis)) = ZERO ! We enforce the radial displacement to zero
              !veloc_elastic(1,ibool(i,j,ispec_axis)) = ZERO ! We enforce the radial displacement to zero
              !accel_elastic(1,ibool(i,j,ispec_axis)) = ZERO ! We enforce the radial displacement to zero
            endif
          enddo
        enddo

      endif

    enddo

  end subroutine enforce_zero_radial_displacements_on_the_axis

