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
!------------------------------------------------------------------------------------------------
!
! symplectic time scheme
!
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

! wrapper functions

  subroutine update_displ_symplectic()

  use specfem_par

  implicit none

  ! time marching
  ! symplectic PEFRL
  if (.not. time_stepping_scheme == 4) return

  ! acoustic domain
  if (ACOUSTIC_SIMULATION) call update_displ_symplectic_acoustic_forward()

  ! elastic domain
  if (ELASTIC_SIMULATION) call update_displ_symplectic_elastic_forward()

  ! poroelastic domain
  if (POROELASTIC_SIMULATION) call update_displ_symplectic_poroelastic_forward()

  end subroutine update_displ_symplectic

!------------------------------------------------------------------------------------------------

  subroutine update_displ_symplectic_backward()

  use specfem_par

  implicit none

  ! time marching
  ! symplectic PEFRL
  if (.not. time_stepping_scheme == 4) return

  ! acoustic domain
  if (ACOUSTIC_SIMULATION) call update_displ_symplectic_acoustic_backward()

  ! elastic domain
  if (ELASTIC_SIMULATION) call update_displ_symplectic_elastic_backward()

  ! poroelastic domain
  if (POROELASTIC_SIMULATION) call update_displ_symplectic_poroelastic_backward()

  end subroutine update_displ_symplectic_backward


!------------------------------------------------------------------------------------------------


  subroutine update_displ_symplectic_acoustic_forward()

! acoustic domains

  use specfem_par

  implicit none

  ! checks if anything to do in this slice
  if (.not. any_acoustic) return

  if (.not. GPU_MODE) then
    ! on CPU
    call update_displacement_symplectic_acoustic(deltat, &
                                                 potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic)
  else
    ! GPU
    call stop_the_code('Symplectic time scheme on GPU for forward acoustic not implemented yet')
  endif

  end subroutine update_displ_symplectic_acoustic_forward

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displ_symplectic_acoustic_backward()

! acoustic domains

  use specfem_par

  implicit none

  ! checks
  if (SIMULATION_TYPE /= 3) return
  ! checks if anything to do in this slice
  if (.not. any_acoustic) return

  if (.not. GPU_MODE) then
    !Since we do not do anything in PML region in case of backward simulation, thus we set
    !PML_BOUNDARY_CONDITIONS = .false.
    call update_displacement_symplectic_acoustic(b_deltat, &
                                                 b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic)
  else
    ! on GPU
    call stop_the_code('Symplectic time scheme on GPU for backward acoustic not implemented yet')
  endif

  end subroutine update_displ_symplectic_acoustic_backward

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displ_symplectic_elastic_forward()

! visco-elastic domains

  use specfem_par

  implicit none

  ! checks if anything to do in this slice
  if (.not. any_elastic) return

  if (.not. GPU_MODE) then
    ! for coupling with adjoint wavefield, stores old (at time t_n) wavefield
    ! not needed anymore, taking care of by re-ordering domain updates
    !if (SIMULATION_TYPE == 3) then
    !  ! handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
    !  ! adjoint definition: \partial_t^2 \bfs^\dagger = - \frac{1}{\rho} \bfnabla \phi^\dagger
    !  if (coupled_acoustic_elastic) then
    !    accel_elastic_adj_coupling(:,:) = - accel_elastic(:,:)
    !  endif
    !endif

    ! updates elastic wavefields
    call update_displacement_symplectic_elastic(deltat,accel_elastic,veloc_elastic,displ_elastic)
  else
    ! on GPU
    call stop_the_code('Symplectic time scheme on GPU for forward elastic not implemented yet')
  endif

  end subroutine update_displ_symplectic_elastic_forward

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displ_symplectic_elastic_backward()

! visco-elastic domains

  use specfem_par

  implicit none

  ! checks
  if (SIMULATION_TYPE /= 3) return
  ! checks if anything to do in this slice
  if (.not. any_elastic) return

  if (.not. GPU_MODE) then
    !Since we do not do anything in PML region in case of backward simulation, thus we set
    !PML_BOUNDARY_CONDITIONS = .false.
    call update_displacement_symplectic_elastic(b_deltat,b_accel_elastic,b_veloc_elastic,b_displ_elastic)
  else
    ! on GPU
    call stop_the_code('Symplectic time scheme on GPU for backward elastic not implemented yet')
  endif

  end subroutine update_displ_symplectic_elastic_backward

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displ_symplectic_poroelastic_forward()

! poro-elastic domains

  use specfem_par

  implicit none

  ! checks if anything to do in this slice
  if (.not. any_poroelastic) return

  if (.not. GPU_MODE) then
    ! for coupling with adjoint wavefield, stores old (at time t_n) wavefield
    ! not needed anymore, taking care of by re-ordering domain updates
    !if (SIMULATION_TYPE == 3) then
    !  ! handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
    !  ! adjoint definition: \partial_t^2 \bfs^\dagger = - \frac{1}{\rho} \bfnabla \phi^\dagger
    !  ! not needed anymore...
    !  accels_poroelastic_adj_coupling(:,:) = - accels_poroelastic(:,:)
    !  accelw_poroelastic_adj_coupling(:,:) = - accelw_poroelastic(:,:)
    !endif

    ! updates poroelastic wavefields
    call update_displacement_symplectic_poroelastic(deltat, &
                                                    accels_poroelastic,velocs_poroelastic,displs_poroelastic, &
                                                    accelw_poroelastic,velocw_poroelastic,displw_poroelastic)
  else
    ! on GPU
    call stop_the_code('Symplectic time scheme on GPU for forward poroelastic not implemented yet')
  endif

  end subroutine update_displ_symplectic_poroelastic_forward


!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displ_symplectic_poroelastic_backward()

! poro-elastic domains

  use specfem_par

  implicit none

  ! checks
  if (SIMULATION_TYPE /= 3) return
  ! checks if anything to do in this slice
  if (.not. any_poroelastic) return

  if (.not. GPU_MODE) then
    !PML not implemented for poroelastic simulation
    call update_displacement_symplectic_poroelastic(b_deltat, &
                                                    b_accels_poroelastic,b_velocs_poroelastic,b_displs_poroelastic, &
                                                    b_accelw_poroelastic,b_velocw_poroelastic,b_displw_poroelastic)
  else
    ! on GPU
    call stop_the_code('Symplectic time scheme on GPU for backward poroelastic not implemented yet')
  endif

  end subroutine update_displ_symplectic_poroelastic_backward


!------------------------------------------------------------------------------------------------
!
! symplectic schemes
!
!------------------------------------------------------------------------------------------------
! Omelyan,  I.M. Mryglod and R. Folk, 2002.
! Optimized Forest-Ruth- and Suzuki-like algorithms for integration of motion in many-body systems,
! Computer Physics communications 146, 188
! http://arxiv.org/abs/cond-mat/0110585
!
! uses PEFRL scheme coefficients from Omelyan et al, 2002, eq. (20)
! together with algorithm in eq. (22) which update position r before velocity v:
!   r1 = r(t) + v(t) xi h
!   v1 = v(t) + 1/m f[r] (1 - 2 lambda) h/2
!
!   r2 = r1 + v1 chi h
!   v2 = v1 + 1/m f[r2] lambda h
!
!   r3 = r2 + v2 (1 - 2 (chi + xi)) h
!   v3 = v2 + 1/m f[r3] lambda h
!
!   r4 = r3 + v3 chi h
!   v(t+h) = v3 + 1/m f[r4] (1 - 2 lambda) h/2
!
!   r(t+h) = r4 + v(t+h) xi h <- final update

  subroutine update_displacement_symplectic_acoustic(deltat, &
                                                     potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic)

! acoustic wavefields

  use constants, only: CUSTOM_REAL,ALPHA_SYMPLECTIC
  use specfem_par, only: nglob_acoustic,i_stage

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: deltat
  real(kind=CUSTOM_REAL), dimension(nglob_acoustic),intent(inout) :: potential_acoustic,potential_dot_acoustic, &
                                                                     potential_dot_dot_acoustic

  ! symplectic
  potential_acoustic(:) = potential_acoustic(:) + ALPHA_SYMPLECTIC(i_stage) * deltat * potential_dot_acoustic(:)
  potential_dot_dot_acoustic(:) = 0.0_CUSTOM_REAL

  end subroutine update_displacement_symplectic_acoustic

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_symplectic_elastic(deltat,accel_elastic,veloc_elastic,displ_elastic)

! elastic wavefields

  use constants, only: CUSTOM_REAL,NDIM,ALPHA_SYMPLECTIC
  use specfem_par, only: nglob_elastic,i_stage

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: deltat
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic),intent(inout) :: accel_elastic,veloc_elastic,displ_elastic

  ! symplectic
  displ_elastic(:,:) = displ_elastic(:,:) + ALPHA_SYMPLECTIC(i_stage) * deltat * veloc_elastic(:,:)
  accel_elastic(:,:) = 0.0_CUSTOM_REAL

  ! PML not implemented yet

  end subroutine update_displacement_symplectic_elastic

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_symplectic_poroelastic(deltat, &
                                                        accels_poroelastic,velocs_poroelastic,displs_poroelastic, &
                                                        accelw_poroelastic,velocw_poroelastic,displw_poroelastic)

! poroelastic wavefields

  use constants, only: CUSTOM_REAL,NDIM,ALPHA_SYMPLECTIC
  use specfem_par, only: nglob_poroelastic,i_stage

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: deltat
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic),intent(inout) :: &
    accels_poroelastic,velocs_poroelastic,displs_poroelastic, &
    accelw_poroelastic,velocw_poroelastic,displw_poroelastic

  ! for the solid
  displs_poroelastic(:,:) = displs_poroelastic(:,:) + ALPHA_SYMPLECTIC(i_stage) * deltat * velocs_poroelastic(:,:)
  accels_poroelastic(:,:) = 0.0_CUSTOM_REAL

  ! for the fluid
  displw_poroelastic(:,:) = displw_poroelastic(:,:) + ALPHA_SYMPLECTIC(i_stage) * deltat * velocw_poroelastic(:,:)
  accelw_poroelastic(:,:) = 0.0_CUSTOM_REAL

  !PML not implemented yet

  end subroutine update_displacement_symplectic_poroelastic

!------------------------------------------------------------------------------------------------
!
! velocity updates & final stage update
!
!------------------------------------------------------------------------------------------------

  subroutine update_veloc_elastic_symplectic()

! elastic wavefields
! updates velocity vector

  use constants, only: ALPHA_SYMPLECTIC,BETA_SYMPLECTIC,NSTAGE_SYMPLECTIC
  use specfem_par

  implicit none

  ! time marching
  veloc_elastic(:,:) = veloc_elastic(:,:) + BETA_SYMPLECTIC(i_stage) * deltat * accel_elastic(:,:)

  ! final update
  if (i_stage == NSTAGE_SYMPLECTIC) then
    ! note: final update uses alpha coefficient xi, same as ALPHA_SYMPLECTIC(1)
    displ_elastic(:,:) = displ_elastic(:,:) + ALPHA_SYMPLECTIC(1) * deltat * veloc_elastic(:,:)
  endif

  end subroutine update_veloc_elastic_symplectic

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_veloc_elastic_symplectic_backward()

! elastic wavefields
! updates backward velocity vector

  use constants, only: ALPHA_SYMPLECTIC,BETA_SYMPLECTIC,NSTAGE_SYMPLECTIC
  use specfem_par

  implicit none

  b_veloc_elastic(:,:) = b_veloc_elastic(:,:) + BETA_SYMPLECTIC(i_stage) * b_deltat * b_accel_elastic(:,:)

  ! final update
  if (i_stage == NSTAGE_SYMPLECTIC) then
    ! note: final update uses alpha coefficient xi, same as ALPHA_SYMPLECTIC(1)
    b_displ_elastic(:,:) = b_displ_elastic(:,:) + ALPHA_SYMPLECTIC(1) * b_deltat * b_veloc_elastic(:,:)
  endif

  end subroutine update_veloc_elastic_symplectic_backward

!------------------------------------------------------------------------------------------------

  subroutine update_veloc_acoustic_symplectic()

! acoustic wavefields
! updates velocity potential

  use constants, only: ALPHA_SYMPLECTIC,BETA_SYMPLECTIC,NSTAGE_SYMPLECTIC
  use specfem_par

  implicit none

  potential_dot_acoustic(:) = potential_dot_acoustic(:) + BETA_SYMPLECTIC(i_stage) * deltat * potential_dot_dot_acoustic(:)

  ! update the potential field (use a new array here) for coupling terms
  if (SIMULATION_TYPE == 3 .and. coupled_acoustic_elastic) then
    stop 'adjoint coupling with symplectic not implemented yet'
  endif

  ! final update
  if (i_stage == NSTAGE_SYMPLECTIC) then
    ! note: final update uses alpha coefficient xi, same as ALPHA_SYMPLECTIC(1)
    potential_acoustic(:) = potential_acoustic(:) + ALPHA_SYMPLECTIC(1) * deltat * potential_dot_acoustic(:)
  endif


  end subroutine update_veloc_acoustic_symplectic


!------------------------------------------------------------------------------------------------

  subroutine update_veloc_acoustic_symplectic_backward()

! acoustic wavefields
! updates velocity potential

  use constants, only: ALPHA_SYMPLECTIC,BETA_SYMPLECTIC,NSTAGE_SYMPLECTIC
  use specfem_par

  implicit none

  b_potential_dot_acoustic(:) = b_potential_dot_acoustic(:) &
                                + BETA_SYMPLECTIC(i_stage) * b_deltat * b_potential_dot_dot_acoustic(:)

  ! final update
  if (i_stage == NSTAGE_SYMPLECTIC) then
    ! note: final update uses alpha coefficient xi, same as ALPHA_SYMPLECTIC(1)
    b_potential_acoustic(:) = b_potential_acoustic(:) + ALPHA_SYMPLECTIC(1) * b_deltat * b_potential_dot_acoustic(:)
  endif

  end subroutine update_veloc_acoustic_symplectic_backward

!------------------------------------------------------------------------------------------------

  subroutine update_veloc_poroelastic_symplectic()

! poroelastic wavefields
! updates velocity

  use constants, only: ALPHA_SYMPLECTIC,BETA_SYMPLECTIC,NSTAGE_SYMPLECTIC
  use specfem_par

  implicit none

  ! solid
  velocs_poroelastic(:,:) = velocs_poroelastic(:,:) + BETA_SYMPLECTIC(i_stage) * deltat * accels_poroelastic(:,:)
  ! fluid
  velocw_poroelastic(:,:) = velocw_poroelastic(:,:) + BETA_SYMPLECTIC(i_stage) * deltat * accelw_poroelastic(:,:)

  ! final update
  if (i_stage == NSTAGE_SYMPLECTIC) then
    ! note: final update uses alpha coefficient xi, same as ALPHA_SYMPLECTIC(1)
    displs_poroelastic(:,:) = displs_poroelastic(:,:) + ALPHA_SYMPLECTIC(1) * deltat * velocs_poroelastic(:,:)
    displw_poroelastic(:,:) = displw_poroelastic(:,:) + ALPHA_SYMPLECTIC(1) * deltat * velocw_poroelastic(:,:)
  endif

  end subroutine update_veloc_poroelastic_symplectic


!------------------------------------------------------------------------------------------------

  subroutine update_veloc_poroelastic_symplectic_backward()

! poroelastic wavefields
! updates velocity

  use constants, only: ALPHA_SYMPLECTIC,BETA_SYMPLECTIC,NSTAGE_SYMPLECTIC
  use specfem_par

  implicit none

  ! solid
  b_velocs_poroelastic(:,:) = b_velocs_poroelastic(:,:) + BETA_SYMPLECTIC(i_stage) * b_deltat * b_accels_poroelastic(:,:)
  ! fluid
  b_velocw_poroelastic(:,:) = b_velocw_poroelastic(:,:) + BETA_SYMPLECTIC(i_stage) * b_deltat * b_accelw_poroelastic(:,:)

  ! final update
  if (i_stage == NSTAGE_SYMPLECTIC) then
    ! note: final update uses alpha coefficient xi, same as ALPHA_SYMPLECTIC(1)
    b_displs_poroelastic(:,:) = b_displs_poroelastic(:,:) + ALPHA_SYMPLECTIC(1) * b_deltat * b_velocs_poroelastic(:,:)
    b_displw_poroelastic(:,:) = b_displw_poroelastic(:,:) + ALPHA_SYMPLECTIC(1) * b_deltat * b_velocw_poroelastic(:,:)
  endif

  end subroutine update_veloc_poroelastic_symplectic_backward

