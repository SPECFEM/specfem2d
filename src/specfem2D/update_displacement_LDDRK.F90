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

!------------------------------------------------------------------------------------------------
!
! Low-dissipation and low-dispersion fourth-order Runge-Kutta algorithm
!
! reference:
! J. Berland, C. Bogey, and C. Bailly.
! Low-dissipation and low-dispersion fourth-order Runge-Kutta algorithm.
! Computers and Fluids, 35:1459-1463, 2006
!
! see: http://www.sciencedirect.com/science/article/pii/S0045793005000575?np=y
!
!------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------
!
! elastic domains
!
!-------------------------------------------------------------------------------------------------


  subroutine update_veloc_elastic_LDDRK()

! updates acceleration,velocity and displacement in elastic regions (crust/mantle,inner core)

  use constants, only: ALPHA_LDDRK,BETA_LDDRK
  use specfem_par

  implicit none

  ! local parameters

  !! DK DK this should be vectorized
  veloc_elastic_LDDRK(:,:) = ALPHA_LDDRK(i_stage) * veloc_elastic_LDDRK(:,:) + deltat * accel_elastic(:,:)
  displ_elastic_LDDRK(:,:) = ALPHA_LDDRK(i_stage) * displ_elastic_LDDRK(:,:) + deltat * veloc_elastic(:,:)

  if (i_stage == 1 .and. it == 1 .and. (.not. initialfield)) then
    veloc_elastic_LDDRK_temp(:,:) = veloc_elastic_LDDRK_temp(:,:) + BETA_LDDRK(i_stage) * veloc_elastic_LDDRK(:,:)
    veloc_elastic(:,:) = veloc_elastic_LDDRK_temp(:,:)
  else
    veloc_elastic(:,:) = veloc_elastic(:,:) + BETA_LDDRK(i_stage) * veloc_elastic_LDDRK(:,:)
  endif

  displ_elastic(:,:) = displ_elastic(:,:) + BETA_LDDRK(i_stage) * displ_elastic_LDDRK(:,:)

  end subroutine update_veloc_elastic_lddrk


!-------------------------------------------------------------------------------------------------
!
! acoustic domains
!
!------------------------------------------------------------------------------------------------

  subroutine update_veloc_acoustic_LDDRK()

! updates velocity potential (corrector)

  use constants, only: ALPHA_LDDRK,BETA_LDDRK
  use specfem_par

  implicit none

  !! DK DK this should be vectorized
  potential_dot_acoustic_LDDRK(:) = ALPHA_LDDRK(i_stage) * potential_dot_acoustic_LDDRK(:) + &
                                    deltat * potential_dot_dot_acoustic(:)
  potential_acoustic_LDDRK(:) = ALPHA_LDDRK(i_stage) * potential_acoustic_LDDRK(:) + &
                                deltat * potential_dot_acoustic(:)

  if (i_stage == 1 .and. it == 1 .and. (.not. initialfield)) then
    !! DK DK this should be vectorized
    potential_dot_acoustic_temp(:) = potential_dot_acoustic_temp(:) + &
                                     BETA_LDDRK(i_stage) * potential_dot_acoustic_LDDRK(:)
    potential_dot_acoustic(:) = potential_dot_acoustic_temp(:)
  else
    potential_dot_acoustic(:) = potential_dot_acoustic(:) + BETA_LDDRK(i_stage) * potential_dot_acoustic_LDDRK(:)
  endif

  !! DK DK this should be vectorized
  potential_acoustic(:) = potential_acoustic(:) + BETA_LDDRK(i_stage) * potential_acoustic_LDDRK(:)

  end subroutine update_veloc_acoustic_LDDRK

!------------------------------------------------------------------------------------------------
!
! poroelastic domains
!
!------------------------------------------------------------------------------------------------

  subroutine update_veloc_poroelastic_LDDRK()

! updates velocity potential (corrector)

  use constants, only: ALPHA_LDDRK,BETA_LDDRK
  use specfem_par

  implicit none

  ! solid
  velocs_poroelastic_LDDRK(:,:) = ALPHA_LDDRK(i_stage) * velocs_poroelastic_LDDRK(:,:) + deltat * accels_poroelastic(:,:)
  displs_poroelastic_LDDRK(:,:) = ALPHA_LDDRK(i_stage) * displs_poroelastic_LDDRK(:,:) + deltat * velocs_poroelastic(:,:)

  velocs_poroelastic(:,:) = velocs_poroelastic(:,:) + BETA_LDDRK(i_stage) * velocs_poroelastic_LDDRK(:,:)
  displs_poroelastic(:,:) = displs_poroelastic(:,:) + BETA_LDDRK(i_stage) * displs_poroelastic_LDDRK(:,:)

  ! fluid
  velocw_poroelastic_LDDRK(:,:) = ALPHA_LDDRK(i_stage) * velocw_poroelastic_LDDRK(:,:) + deltat * accelw_poroelastic(:,:)
  displw_poroelastic_LDDRK(:,:) = ALPHA_LDDRK(i_stage) * displw_poroelastic_LDDRK(:,:) + deltat * velocw_poroelastic(:,:)

  velocw_poroelastic(:,:) = velocw_poroelastic(:,:) + BETA_LDDRK(i_stage) * velocw_poroelastic_LDDRK(:,:)
  displw_poroelastic(:,:) = displw_poroelastic(:,:) + BETA_LDDRK(i_stage) * displw_poroelastic_LDDRK(:,:)

  end subroutine update_veloc_poroelastic_LDDRK

