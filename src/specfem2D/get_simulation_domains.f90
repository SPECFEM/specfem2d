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

  subroutine get_simulation_domains()

  use constants, only: TINYVAL

  use specfem_par, only: any_acoustic,any_elastic,any_poroelastic, &
    ispec_is_anisotropic,ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
    nspec,porosity,anisotropy,kmato

  implicit none

  ! local parameters
  integer :: ispec

  ! initializes
  any_acoustic = .false.
  any_elastic = .false.
  any_poroelastic = .false.

  ispec_is_acoustic(:) = .false.
  ispec_is_anisotropic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  ! loops over all elements
  do ispec = 1,nspec

    ! checks domain properties
    if (nint(porosity(kmato(ispec))) == 1) then
      ! assume acoustic domain
      ispec_is_acoustic(ispec) = .true.
      any_acoustic = .true.

    else if (porosity(kmato(ispec)) < TINYVAL) then
      ! assume elastic domain
      ispec_is_elastic(ispec) = .true.
      any_elastic = .true.
      if (any(anisotropy(:,kmato(ispec)) /= 0)) then
        ispec_is_anisotropic(ispec) = .true.
      endif

    else
      ! assume poroelastic domain
      ispec_is_poroelastic(ispec) = .true.
      any_poroelastic = .true.
    endif

  enddo ! of do ispec = 1,nspec

  ! safety checks
  call get_simulation_domain_check()

  ! sets domain numbers
  call get_simulation_domain_counts()

  end subroutine get_simulation_domains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_simulation_domains_from_external_models()

  use constants, only: TINYVAL,NGLLX,NGLLZ,CUSTOM_REAL

  use specfem_par, only: any_acoustic,any_elastic,any_poroelastic, &
    ispec_is_anisotropic,ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
    nspec,myrank,P_SV,MODEL

  ! external model parameters
  use specfem_par, only: vsext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext

  implicit none

  ! local parameters
  integer :: ispec,i,j
  real(kind=CUSTOM_REAL) :: previous_vsext

  ! re-assigns flags
  ! initializes
  any_acoustic = .false.
  any_elastic = .false.
  any_poroelastic = .false.

  ispec_is_acoustic(:) = .false.
  ispec_is_anisotropic(:) = .false.
  ispec_is_elastic(:) = .false.
  ispec_is_poroelastic(:) = .false.

  do ispec = 1,nspec
    ! value at corner
    previous_vsext = vsext(1,1,ispec)

    do j = 1,NGLLZ
      do i = 1,NGLLX
        if (P_SV .and. (.not. (i == 1 .and. j == 1)) .and. &
          ((vsext(i,j,ispec) >= TINYVAL .and. previous_vsext < TINYVAL) .or. &
           (vsext(i,j,ispec) < TINYVAL .and. previous_vsext >= TINYVAL)))  &
          call exit_MPI(myrank,'external velocity model cannot be both fluid and solid inside the same spectral element')

        if (vsext(i,j,ispec) < TINYVAL) then
          ! acoustic
          ispec_is_acoustic(ispec) = .true.
        else
          ! elastic
          ispec_is_elastic(ispec) = .true.
        endif

        if ( trim(MODEL) == 'external' .or. trim(MODEL) == 'tomo' .or. trim(MODEL) == 'binary_voigt' ) then

          ! sets element type
          if (c11ext(i,j,ispec) > TINYVAL .or. c13ext(i,j,ispec) > TINYVAL .or. c15ext(i,j,ispec) > TINYVAL .or. &
              c33ext(i,j,ispec) > TINYVAL .or. c35ext(i,j,ispec) > TINYVAL .or. c55ext(i,j,ispec) > TINYVAL) then
            ! anisotropic elastic
            ispec_is_anisotropic(ispec) = .true.
            ispec_is_elastic(ispec) = .true.
            ispec_is_acoustic(ispec) = .false.
          else if (vsext(i,j,ispec) < TINYVAL) then
            ispec_is_acoustic(ispec) = .true.
          endif

        endif

        ! sets new GLL point value to compare against
        previous_vsext = vsext(i,j,ispec)
      enddo
    enddo
  enddo ! ispec

  ! sets domain numbers
  call get_simulation_domain_counts()

  ! safety check
  call get_simulation_domain_check()

  end subroutine get_simulation_domains_from_external_models

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_simulation_domain_check()

  use specfem_par

  implicit none
  ! local parameters
  integer :: ispec

  ! checks simulation domain flags
  if (ANY(ispec_is_acoustic(:)) .neqv. any_acoustic) call stop_the_code('Error any_acoustic invalid')
  if (ANY(ispec_is_elastic(:)) .neqv. any_elastic) call stop_the_code('Error any_elastic invalid')
  if (ANY(ispec_is_poroelastic(:)) .neqv. any_poroelastic) call stop_the_code('Error any_poroelastic invalid')

  ! safety checks
  if (.not. P_SV .and. .not. any_elastic) then
    print *, '*************** WARNING ***************'
    print *, 'Surface (membrane) waves calculation needs an elastic medium'
    print *, '*************** WARNING ***************'
    call stop_the_code('Please set P_SV flag to .true. for acoustic simulations')
  endif

  if (PML_BOUNDARY_CONDITIONS .and. any_poroelastic) then
    call stop_the_code('PML boundary conditions not implemented for poroelastic simulations yet')
  endif

  if (PML_BOUNDARY_CONDITIONS .and. any_elastic .and. (.not. P_SV)) then
    call stop_the_code('PML boundary conditions not implemented for SH simulations yet')
  endif

  ! checks material domains
  do ispec = 1,nspec
    ! checks if at least one domain is set
    if ((.not. ispec_is_acoustic(ispec)) .and. (.not. ispec_is_elastic(ispec)) &
        .and. (.not. ispec_is_poroelastic(ispec))) then
      print *,'Error material domain not assigned to element:',ispec
      print *,'acoustic       : ',ispec_is_acoustic(ispec)
      print *,'elastic        : ',ispec_is_elastic(ispec)
      print *,'poroelastic    : ',ispec_is_poroelastic(ispec)
      call stop_the_code('Error material domain index element')
    endif

    ! checks if domain is unique
    if ((ispec_is_acoustic(ispec) .and. ispec_is_elastic(ispec)) .or. &
        (ispec_is_acoustic(ispec) .and. ispec_is_poroelastic(ispec)) .or. &
        (ispec_is_elastic(ispec) .and. ispec_is_poroelastic(ispec))) then
      print *,'Error material domain assigned twice to element:',ispec
      print *,'acoustic       : ',ispec_is_acoustic(ispec)
      print *,'elastic        : ',ispec_is_elastic(ispec)
      print *,'poroelastic    : ',ispec_is_poroelastic(ispec)
      call stop_the_code('Error material domain index element')
    endif
  enddo

  end subroutine get_simulation_domain_check

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_simulation_domain_counts()

! counts number of elements in different domains

  use specfem_par

  implicit none

  ! acoustic
  ! number of acoustic elements in this partition
  nspec_acoustic = count(ispec_is_acoustic(:))
  if (nspec_acoustic > 0 ) any_acoustic = .true.

  ! elastic
  ! number of elastic elements in this partition
  nspec_elastic = count(ispec_is_elastic(:))
  if (nspec_elastic > 0 ) any_elastic = .true.

  ! poroelastic
  ! number of elastic elements in this partition
  nspec_poroelastic = count(ispec_is_poroelastic(:))
  if (nspec_poroelastic > 0 ) any_poroelastic = .true.

  end subroutine get_simulation_domain_counts

