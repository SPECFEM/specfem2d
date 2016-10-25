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
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

  subroutine compute_forces_gravitoacoustic_main()

  use constants, only: ZERO
  use specfem_par

  implicit none

  ! local parameters
  integer :: iglob

  ! main solver for the gravitoacoustic elements

! *********************************************************
! only SIMULATION_TYPE == 1, time_stepping_scheme == 1, and no PML or STACEY yet
! NO MIX OF ACOUSTIC AND GRAVITOACOUTIC ELEMENTS
! NO COUPLING TO ELASTIC AND POROELASTIC SIDES
! *********************************************************

  if ((any_gravitoacoustic)) then
    if (time_stepping_scheme == 1) then
      ! Newmark time scheme
      !! DK DK this should be vectorized
      potential_gravitoacoustic = potential_gravitoacoustic + deltat*potential_dot_gravitoacoustic + &
                                  deltatsquareover2*potential_dot_dot_gravitoacoustic
      potential_dot_gravitoacoustic = potential_dot_gravitoacoustic + &
                                      deltatover2*potential_dot_dot_gravitoacoustic
      potential_gravito = potential_gravito + deltat*potential_dot_gravito + &
                          deltatsquareover2*potential_dot_dot_gravito
      potential_dot_gravito = potential_dot_gravito + deltatover2*potential_dot_dot_gravito
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
    if (ACOUSTIC_FORCING) then
      call add_acoustic_forcing_at_rigid_boundary_gravitoacoustic()
    endif ! end ACOUSTIC_FORCING !

! free surface for a gravitoacoustic medium
!!! to be coded !!!
!        call enforce_acoustic_free_surface(potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
!                                          potential_gravitoacoustic)

!        if (SIMULATION_TYPE == 3) then ! Adjoint calculation
!          call enforce_acoustic_free_surface(b_potential_dot_dot_gravitoacoustic,b_potential_dot_gravitoacoustic, &
!                                            b_potential_gravitoacoustic)
!        endif

! *********************************************************
! ************* compute forces for the gravitoacoustic elements
! *********************************************************

    call compute_forces_gravitoacoustic(potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
                                        potential_gravitoacoustic, potential_dot_dot_gravito, &
                                        potential_gravito,.false.)

    ! debugging
    if ((mod(it,100) == 0)) then
      iglob=iglobzero
      write(*,*)it, & ! Nsql,gravityl,
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
!    if (NPROC > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
!      call assemble_MPI_scalar_ac_blocking(potential_dot_dot_gravitoacoustic)
!
!    endif
!
!#endif

! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

  if ((any_gravitoacoustic)) then
    if (time_stepping_scheme == 1) then
      !! DK DK this should be vectorized
      potential_dot_dot_gravitoacoustic(:) = potential_dot_dot_gravitoacoustic(:) * rmass_inverse_gravitoacoustic(:)
      potential_dot_gravitoacoustic(:) = potential_dot_gravitoacoustic(:) + &
                                         deltatover2 * potential_dot_dot_gravitoacoustic(:)

!! line below already done in compute_forces_gravitoacoustic, because necessary
!! for the computation of potential_dot_dot_gravitoacoustic
!      potential_dot_dot_gravito = potential_dot_dot_gravito * rmass_inverse_gravito
      potential_dot_gravito(:) = potential_dot_gravito(:) + deltatover2 * potential_dot_dot_gravito(:)
    else
      stop 'Only time_stepping_scheme = 1 implemented for gravitoacoustic case'
    endif

! free surface for an acoustic medium
!        call enforce_acoustic_free_surface(potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
!                                        potential_gravitoacoustic)
!
!        if (SIMULATION_TYPE == 3) then
!          call enforce_acoustic_free_surface(b_potential_dot_dot_gravitoacoustic,b_potential_dot_gravitoacoustic, &
!                                          b_potential_gravitoacoustic)
!        endif
!
      ! update the potential field (use a new array here) for coupling terms
!      potential_gravitoacoustic_adj_coupling = potential_gravitoacoustic &
!                          + deltat*potential_dot_gravitoacoustic &
!                          + deltatsquareover2*potential_dot_dot_gravitoacoustic

  endif ! of if (any_gravitoacoustic)

  end subroutine compute_forces_gravitoacoustic_main

