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
! 4th-order Runge-Kutta time scheme
!
!------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
!
! elastic domains
!
!-------------------------------------------------------------------------------------------------


  subroutine update_veloc_elastic_RK()

! updates acceleration,velocity and displacement in elastic regions (crust/mantle,inner core)

  use specfem_par

  implicit none

  ! local parameters
  double precision :: weight_rk

  !! DK DK this should be vectorized
  accel_elastic_rk(1,:,i_stage) = deltat * accel_elastic(1,:)
  accel_elastic_rk(2,:,i_stage) = deltat * accel_elastic(2,:)

  veloc_elastic_rk(1,:,i_stage) = deltat * veloc_elastic(1,:)
  veloc_elastic_rk(2,:,i_stage) = deltat * veloc_elastic(2,:)

  select case (i_stage)
  case (1,2,3)
    ! weight
    if (i_stage == 1) then
      weight_rk = 0.5d0
    else if (i_stage == 2) then
      weight_rk = 0.5d0
    else if (i_stage == 3) then
      weight_rk = 1.0d0
    endif

    if (i_stage == 1) then
      !! DK DK this should be vectorized
      veloc_elastic_initial_rk(1,:) = veloc_elastic(1,:)
      veloc_elastic_initial_rk(2,:) = veloc_elastic(2,:)

      displ_elastic_initial_rk(1,:) = displ_elastic(1,:)
      displ_elastic_initial_rk(2,:) = displ_elastic(2,:)
    endif

    !! DK DK this should be vectorized
    veloc_elastic(1,:) = veloc_elastic_initial_rk(1,:) + weight_rk * accel_elastic_rk(1,:,i_stage)
    veloc_elastic(2,:) = veloc_elastic_initial_rk(2,:) + weight_rk * accel_elastic_rk(2,:,i_stage)

    displ_elastic(1,:) = displ_elastic_initial_rk(1,:) + weight_rk * veloc_elastic_rk(1,:,i_stage)
    displ_elastic(2,:) = displ_elastic_initial_rk(2,:) + weight_rk * veloc_elastic_rk(2,:,i_stage)

  case (4)
    !! DK DK this should be vectorized
    veloc_elastic(1,:) = veloc_elastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
                         ( accel_elastic_rk(1,:,1) + 2.0d0 * accel_elastic_rk(1,:,2) + &
                           2.0d0 * accel_elastic_rk(1,:,3) + accel_elastic_rk(1,:,4) )
    veloc_elastic(2,:) = veloc_elastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
                         ( accel_elastic_rk(2,:,1) + 2.0d0 * accel_elastic_rk(2,:,2) + &
                           2.0d0 * accel_elastic_rk(2,:,3) + accel_elastic_rk(2,:,4) )

    displ_elastic(1,:) = displ_elastic_initial_rk(1,:) + 1.0d0 / 6.0d0 * &
                         ( veloc_elastic_rk(1,:,1) + 2.0d0 * veloc_elastic_rk(1,:,2) + &
                           2.0d0 * veloc_elastic_rk(1,:,3) + veloc_elastic_rk(1,:,4) )
    displ_elastic(2,:) = displ_elastic_initial_rk(2,:) + 1.0d0 / 6.0d0 * &
                         ( veloc_elastic_rk(2,:,1) + 2.0d0 * veloc_elastic_rk(2,:,2) + &
                           2.0d0 * veloc_elastic_rk(2,:,3) + veloc_elastic_rk(2,:,4))
  end select

  end subroutine update_veloc_elastic_RK


!-------------------------------------------------------------------------------------------------
!
! acoustic domains
!
!------------------------------------------------------------------------------------------------

  subroutine update_veloc_acoustic_RK()

! updates velocity potential (corrector)

  use specfem_par

  implicit none

  ! local parameters
  double precision :: weight_rk

  !! DK DK this should be vectorized
  potential_dot_dot_acoustic_rk(:,i_stage) = deltat * potential_dot_dot_acoustic(:)
  potential_dot_acoustic_rk(:,i_stage) = deltat * potential_dot_acoustic(:)

  select case (i_stage)
  case (1,2,3)
    ! selects weight
    if (i_stage == 1) then
      weight_rk = 0.5d0
    else if (i_stage == 2) then
      weight_rk = 0.5d0
    else if (i_stage == 3) then
      weight_rk = 1.0d0
    endif

    if (i_stage == 1) then
!! DK DK this should be vectorized
      potential_dot_acoustic_init_rk(:) = potential_dot_acoustic(:)
      potential_acoustic_init_rk(:) = potential_acoustic(:)
    endif

!! DK DK this should be vectorized
    potential_dot_acoustic(:) = potential_dot_acoustic_init_rk(:) + &
                                weight_rk * potential_dot_dot_acoustic_rk(:,i_stage)
    potential_acoustic(:) = potential_acoustic_init_rk(:) + weight_rk * potential_dot_acoustic_rk(:,i_stage)

  case(4)
!! DK DK this should be vectorized
    potential_dot_acoustic(:) = potential_dot_acoustic_init_rk(:) + &
                                1.0d0 / 6.0d0 * ( potential_dot_dot_acoustic_rk(:,1) + &
                                                  2.0d0 * potential_dot_dot_acoustic_rk(:,2) + &
                                                  2.0d0 * potential_dot_dot_acoustic_rk(:,3) + &
                                                  potential_dot_dot_acoustic_rk(:,4) )

!! DK DK this should be vectorized
    potential_acoustic(:) = potential_acoustic_init_rk(:) + &
                            1.0d0 / 6.0d0 * ( potential_dot_acoustic_rk(:,1) + &
                                              2.0d0 * potential_dot_acoustic_rk(:,2) + &
                                              2.0d0 * potential_dot_acoustic_rk(:,3) + &
                                              potential_dot_acoustic_rk(:,4) )
  end select

  end subroutine update_veloc_acoustic_RK


!------------------------------------------------------------------------------------------------
!
! poroelastic domains
!
!------------------------------------------------------------------------------------------------

  subroutine update_veloc_poroelastic_RK()

! updates velocity potential (corrector)

  use specfem_par

  implicit none

  ! local parameters
  double precision :: weight_rk

  ! solid
  accels_poroelastic_rk(1,:,i_stage) = deltat * accels_poroelastic(1,:)
  accels_poroelastic_rk(2,:,i_stage) = deltat * accels_poroelastic(2,:)
  velocs_poroelastic_rk(1,:,i_stage) = deltat * velocs_poroelastic(1,:)
  velocs_poroelastic_rk(2,:,i_stage) = deltat * velocs_poroelastic(2,:)

  ! fluid
  accelw_poroelastic_rk(1,:,i_stage) = deltat * accelw_poroelastic(1,:)
  accelw_poroelastic_rk(2,:,i_stage) = deltat * accelw_poroelastic(2,:)
  velocw_poroelastic_rk(1,:,i_stage) = deltat * velocw_poroelastic(1,:)
  velocw_poroelastic_rk(2,:,i_stage) = deltat * velocw_poroelastic(2,:)

  select case (i_stage)
  case (1,2,3)
    if (i_stage == 1) then
      weight_rk = 0.5d0
    else if (i_stage == 2) then
      weight_rk = 0.5d0
    else if (i_stage == 3) then
      weight_rk = 1.0d0
    endif

    if (i_stage == 1) then
      ! solid
      velocs_poroelastic_initial_rk(1,:) = velocs_poroelastic(1,:)
      velocs_poroelastic_initial_rk(2,:) = velocs_poroelastic(2,:)
      displs_poroelastic_initial_rk(1,:) = displs_poroelastic(1,:)
      displs_poroelastic_initial_rk(2,:) = displs_poroelastic(2,:)
      ! fluid
      velocw_poroelastic_initial_rk(1,:) = velocw_poroelastic(1,:)
      velocw_poroelastic_initial_rk(2,:) = velocw_poroelastic(2,:)
      displw_poroelastic_initial_rk(1,:) = displw_poroelastic(1,:)
      displw_poroelastic_initial_rk(2,:) = displw_poroelastic(2,:)
    endif

    ! solid
    velocs_poroelastic(1,:) = velocs_poroelastic_initial_rk(1,:) + weight_rk * accels_poroelastic_rk(1,:,i_stage)
    velocs_poroelastic(2,:) = velocs_poroelastic_initial_rk(2,:) + weight_rk * accels_poroelastic_rk(2,:,i_stage)
    displs_poroelastic(1,:) = displs_poroelastic_initial_rk(1,:) + weight_rk * velocs_poroelastic_rk(1,:,i_stage)
    displs_poroelastic(2,:) = displs_poroelastic_initial_rk(2,:) + weight_rk * velocs_poroelastic_rk(2,:,i_stage)
    ! fluid
    velocw_poroelastic(1,:) = velocw_poroelastic_initial_rk(1,:) + weight_rk * accelw_poroelastic_rk(1,:,i_stage)
    velocw_poroelastic(2,:) = velocw_poroelastic_initial_rk(2,:) + weight_rk * accelw_poroelastic_rk(2,:,i_stage)
    displw_poroelastic(1,:) = displw_poroelastic_initial_rk(1,:) + weight_rk * velocw_poroelastic_rk(1,:,i_stage)
    displw_poroelastic(2,:) = displw_poroelastic_initial_rk(2,:) + weight_rk * velocw_poroelastic_rk(2,:,i_stage)

  case (4)
    ! solid
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
  end select

  end subroutine update_veloc_poroelastic_RK


