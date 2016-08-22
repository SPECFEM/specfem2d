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

!------------------------------------------------------------------------------------------------
!
! Newmark time scheme
!
!------------------------------------------------------------------------------------------------
!
! explicit Newmark time scheme with acoustic & elastic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! chi(t+delta_t) = chi(t) + delta_t chi_dot(t) + 1/2 delta_t**2 chi_dot_dot(t)
! chi_dot(t+delta_t) = chi_dot(t) + 1/2 delta_t chi_dot_dot(t) + 1/2 delta_t chi_dot_dot(t+delta_t)
! chi_dot_dot(t+delta_t) = 1/M_acoustic( -K_acoustic chi(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! u(t+delta_t) = u(t) + delta_t  v(t) + 1/2  delta_t**2 a(t)
! v(t+delta_t) = v(t) + 1/2 delta_t a(t) + 1/2 delta_t a(t+delta_t)
! a(t+delta_t) = 1/M_elastic ( -K_elastic u(t+delta) + B_elastic chi_dot_dot(t+delta_t) + f( t+delta_t) )
!
! where
!   chi, chi_dot, chi_dot_dot are acoustic (fluid) potentials ( dotted with respect to time)
!   u, v, a are displacement,velocity & acceleration
!   M is mass matrix, K stiffness matrix and B boundary term for acoustic/elastic domains
!   f denotes a source term (acoustic/elastic)
!
! note that this stage calculates the predictor terms
!
!   for
!   potential chi_dot(t+delta) requires + 1/2 delta_t chi_dot_dot(t+delta_t)
!                                   at a later stage (corrector) once where chi_dot_dot(t+delta) is calculated
!   and similar,
!   velocity v(t+delta_t) requires  + 1/2 delta_t a(t+delta_t)
!                                   at a later stage once where a(t+delta) is calculated
! also:
!   boundary term B_elastic requires chi_dot_dot(t+delta)
!                                   thus chi_dot_dot has to be updated first before the elastic boundary term is considered


!------------------------------------------------------------------------------------------------

! predictor phase

!------------------------------------------------------------------------------------------------

  subroutine update_displ_Newmark()

  use specfem_par

  implicit none

  ! time marching

  ! acoustic domain
  if (ACOUSTIC_SIMULATION) call update_displ_acoustic_forward()

  ! elastic domain
  if (ELASTIC_SIMULATION) call update_displ_elastic_forward()

  ! poroelastic domain
  if (POROELASTIC_SIMULATION) call update_displ_poroelastic_forward()

  end subroutine update_displ_Newmark

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displ_Newmark_backward()

  use specfem_par

  implicit none

  ! time marching

  ! acoustic domain
  if (ACOUSTIC_SIMULATION) call update_displ_acoustic_backward()

  ! elastic domain
  if (ELASTIC_SIMULATION) call update_displ_elastic_backward()

  ! poroelastic domain
  if (POROELASTIC_SIMULATION) call update_displ_poroelastic_backward()

  end subroutine update_displ_Newmark_backward

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displ_acoustic_forward()

! acoustic domains

  use specfem_par

  implicit none

#ifdef FORCE_VECTORIZATION
  integer :: i
#endif

  ! checks if anything to do in this slice
  if (.not. any_acoustic) return

  if (.not. GPU_MODE) then

    if (time_stepping_scheme == 1) then
      call update_displacement_newmark_acoustic(deltat,deltatover2,deltatsquareover2, &
                                                potential_dot_dot_acoustic,potential_dot_acoustic, &
                                                potential_acoustic, &
                                                PML_BOUNDARY_CONDITIONS,potential_acoustic_old)
    else
#ifdef FORCE_VECTORIZATION
      do i = 1,nglob_acoustic
        potential_dot_dot_acoustic(i) = 0._CUSTOM_REAL
      enddo
#else
      potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
#endif
    endif

  else
    ! on GPU
    ! handles both forward and backward
    call update_displacement_newmark_GPU_acoustic()
  endif

  end subroutine update_displ_acoustic_forward

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displ_acoustic_backward()

! acoustic domains

  use specfem_par

  implicit none

#ifdef FORCE_VECTORIZATION
  integer :: i
#endif

  ! checks
  if (SIMULATION_TYPE /= 3) return

  ! checks if anything to do in this slice
  if (.not. any_acoustic) return

  if (.not. GPU_MODE) then
    !Since we do not do anything in PML region in case of backward simulation, thus we set
    !PML_BOUNDARY_CONDITIONS = .false.
    if (time_stepping_scheme == 1) then
      call update_displacement_newmark_acoustic(b_deltat,b_deltatover2,b_deltatsquareover2, &
                                                b_potential_dot_dot_acoustic,b_potential_dot_acoustic, &
                                                b_potential_acoustic, &
                                                .false.,b_potential_acoustic_old)
    else
#ifdef FORCE_VECTORIZATION
      do i = 1,nglob_acoustic
        b_potential_dot_dot_acoustic(i) = 0._CUSTOM_REAL
      enddo
#else
      b_potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
#endif
    endif
  else
    ! on GPU
    ! already done in forward call...
    continue
  endif

  end subroutine update_displ_acoustic_backward

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displ_elastic_forward()

! visco-elastic domains

  use specfem_par

  implicit none

#ifdef FORCE_VECTORIZATION
  integer :: i
#endif

  ! checks if anything to do in this slice
  if (.not. any_elastic) return

  if (.not. GPU_MODE) then

    ! for coupling with adjoint wavefield, stores old (at time t_n) wavefield
    if (SIMULATION_TYPE == 3) then
      if (time_stepping_scheme == 1) then
        ! handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
        ! adjoint definition: \partial_t^2 \bfs^\dagger = - \frac{1}{\rho} \bfnabla \phi^\dagger
        if (coupled_acoustic_elastic) then
#ifdef FORCE_VECTORIZATION
          do i = 1,NDIM*nglob_elastic
            accel_elastic_adj_coupling(i,1) = - accel_elastic(i,1)
          enddo
#else
          accel_elastic_adj_coupling(:,:) = - accel_elastic(:,:)
#endif
        endif
      endif
    endif

    ! updates elastic wavefields
    if (time_stepping_scheme == 1) then
      call update_displacement_newmark_elastic(deltat,deltatover2,deltatsquareover2, &
                                               accel_elastic,veloc_elastic, &
                                               displ_elastic,displ_elastic_old, &
                                               PML_BOUNDARY_CONDITIONS)
    else
#ifdef FORCE_VECTORIZATION
      do i = 1,NDIM*nglob_elastic
        accel_elastic(i,1) = 0._CUSTOM_REAL
      enddo
#else
      accel_elastic(:,:) = 0._CUSTOM_REAL
#endif
    endif
  else
    ! on GPU
    ! handles both forward and backward
    call update_displacement_newmark_GPU_elastic()
  endif

  end subroutine update_displ_elastic_forward

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displ_elastic_backward()

! visco-elastic domains

  use specfem_par

  implicit none

#ifdef FORCE_VECTORIZATION
  integer :: i
#endif

  ! checks
  if (SIMULATION_TYPE /= 3) return

  ! checks if anything to do in this slice
  if (.not. any_elastic) return

  if (.not. GPU_MODE) then
    !Since we do not do anything in PML region in case of backward simulation, thus we set
    !PML_BOUNDARY_CONDITIONS = .false.
    if (time_stepping_scheme == 1) then
      call update_displacement_newmark_elastic(b_deltat,b_deltatover2,b_deltatsquareover2, &
                                               b_accel_elastic,b_veloc_elastic, &
                                               b_displ_elastic,b_displ_elastic_old, &
                                               .false.)
    else
#ifdef FORCE_VECTORIZATION
      do i = 1,NDIM*nglob_elastic
        b_accel_elastic(i,1) = 0._CUSTOM_REAL
      enddo
#else
      b_accel_elastic(:,:) = 0._CUSTOM_REAL
#endif
    endif
  else
    ! on GPU
    ! already done in forward call..
    continue
  endif

  end subroutine update_displ_elastic_backward

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displ_poroelastic_forward()

! poro-elastic domains

  use specfem_par

  implicit none

  ! checks if anything to do in this slice
  if (.not. any_poroelastic) return

  if (.not. GPU_MODE) then

    ! for coupling with adjoint wavefield, stores old (at time t_n) wavefield
    if (SIMULATION_TYPE == 3) then
      if (time_stepping_scheme == 1) then
        ! handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
        ! adjoint definition: \partial_t^2 \bfs^\dagger = - \frac{1}{\rho} \bfnabla \phi^\dagger
        accels_poroelastic_adj_coupling(:,:) = - accels_poroelastic(:,:)
        accelw_poroelastic_adj_coupling(:,:) = - accelw_poroelastic(:,:)
      endif
    endif

    ! updates poroelastic wavefields
    if (time_stepping_scheme == 1) then
      call update_displacement_newmark_poroelastic(deltat,deltatover2,deltatsquareover2, &
                                                   accels_poroelastic,velocs_poroelastic, &
                                                   displs_poroelastic,accelw_poroelastic, &
                                                   velocw_poroelastic,displw_poroelastic)
    else
      accels_poroelastic(:,:) = 0._CUSTOM_REAL
      accelw_poroelastic(:,:) = 0._CUSTOM_REAL
    endif
  else
    ! on GPU
    ! safety stop
    call exit_MPI(myrank,'poroelastic time marching scheme on GPU not implemented yet...')
  endif

  end subroutine update_displ_poroelastic_forward


!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displ_poroelastic_backward()

! poro-elastic domains

  use specfem_par

  implicit none

  ! checks
  if (SIMULATION_TYPE /= 3) return

  ! checks if anything to do in this slice
  if (.not. any_poroelastic) return

  if (.not. GPU_MODE) then

    if (time_stepping_scheme == 1) then
      !PML not implemented for poroelastic simulation
      call update_displacement_newmark_poroelastic(b_deltat,b_deltatover2,b_deltatsquareover2, &
                                                   b_accels_poroelastic,b_velocs_poroelastic, &
                                                   b_displs_poroelastic,b_accelw_poroelastic, &
                                                   b_velocw_poroelastic,b_displw_poroelastic)
    else
      b_accels_poroelastic(:,:) = 0._CUSTOM_REAL
      b_accelw_poroelastic(:,:) = 0._CUSTOM_REAL
    endif
  else
    ! on GPU
    ! safety stop
    call exit_MPI(myrank,'poroelastic time marching scheme on GPU not implemented yet...')
  endif

  end subroutine update_displ_poroelastic_backward


!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_newmark_acoustic(deltat,deltatover2,deltatsquareover2, &
                                                  potential_dot_dot_acoustic,potential_dot_acoustic, &
                                                  potential_acoustic, &
                                                  PML_BOUNDARY_CONDITIONS,potential_acoustic_old)

  use specfem_par, only: nglob_acoustic,CUSTOM_REAL

  implicit none

  double precision,intent(in) :: deltat,deltatover2,deltatsquareover2
  real(kind=CUSTOM_REAL), dimension(nglob_acoustic),intent(inout) :: potential_acoustic,potential_dot_acoustic, &
                                                                     potential_dot_dot_acoustic

  logical,intent(in) :: PML_BOUNDARY_CONDITIONS
  real(kind=CUSTOM_REAL), dimension(nglob_acoustic),intent(inout) :: potential_acoustic_old

  ! local parameters
#ifdef FORCE_VECTORIZATION
  integer :: i
#endif

  ! PML simulations
  if (PML_BOUNDARY_CONDITIONS) then

    ! note: todo - for elastic, there is an additional factor 1/TWO to the default deltasquareover2 for the acceleration term
    !       find explanations where?
#ifdef FORCE_VECTORIZATION
    do i = 1,nglob_acoustic
      potential_acoustic_old(i) = potential_acoustic(i) + deltatsquareover2 * potential_dot_dot_acoustic(i)
    enddo
#else
    potential_acoustic_old(:) = potential_acoustic(:) + deltatsquareover2 * potential_dot_dot_acoustic(:)
#endif

  endif ! PML_BOUNDARY_CONDITIONS


#ifdef FORCE_VECTORIZATION
  do i = 1,nglob_acoustic
    potential_acoustic(i) = potential_acoustic(i) + deltat * potential_dot_acoustic(i) &
                                                  + deltatsquareover2 * potential_dot_dot_acoustic(i)
    potential_dot_acoustic(i) = potential_dot_acoustic(i) + deltatover2 * potential_dot_dot_acoustic(i)
    potential_dot_dot_acoustic(i) = 0._CUSTOM_REAL
  enddo
#else
  potential_acoustic(:) = potential_acoustic(:) + deltat * potential_dot_acoustic &
                                                + deltatsquareover2 * potential_dot_dot_acoustic(:)
  potential_dot_acoustic(:) = potential_dot_acoustic(:) + deltatover2 * potential_dot_dot_acoustic(:)
  potential_dot_dot_acoustic(:) = 0._CUSTOM_REAL
#endif

  end subroutine update_displacement_newmark_acoustic

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_newmark_elastic(deltat,deltatover2,deltatsquareover2, &
                                                 accel_elastic,veloc_elastic, &
                                                 displ_elastic,displ_elastic_old, &
                                                 PML_BOUNDARY_CONDITIONS)

  use constants, only: CUSTOM_REAL,NDIM,TWO

  use specfem_par, only: nglob_elastic,ATTENUATION_VISCOELASTIC

#ifndef FORCE_VECTORIZATION
  use constants, only: USE_ENFORCE_FIELDS
  use specfem_par, only: iglob_is_forced,it
#endif

  implicit none

  double precision,intent(in) :: deltat,deltatover2,deltatsquareover2
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic),intent(inout) :: accel_elastic,veloc_elastic, &
                                                                      displ_elastic,displ_elastic_old

  logical,intent(in) :: PML_BOUNDARY_CONDITIONS

  ! local parameters
#ifdef FORCE_VECTORIZATION
  integer :: i
#else
  integer :: iglob
#endif

  ! attenuation/PML simulations
  if (PML_BOUNDARY_CONDITIONS .or. ATTENUATION_VISCOELASTIC) then

    ! note: todo - there is an additional factor 1/TWO to the default deltasquareover2 for the acceleration term
    !       find explanations where?
#ifdef FORCE_VECTORIZATION
    do i = 1,NDIM*nglob_elastic
      displ_elastic_old(i,1) = displ_elastic(i,1) + deltatsquareover2/TWO * accel_elastic(i,1)
    enddo
#else
    displ_elastic_old(:,:) = displ_elastic(:,:) + deltatsquareover2/TWO * accel_elastic(:,:)
#endif

  endif ! PML

#ifdef FORCE_VECTORIZATION
  !! DK DK this allows for full vectorization by using a trick to see the 2D array as a 1D array
  !! DK DK whose dimension is the product of the two dimensions, the second dimension being equal to 1
  ! vectorized updates
  do i = 1,NDIM*nglob_elastic
    displ_elastic(i,1) = displ_elastic(i,1) + deltat * veloc_elastic(i,1) + deltatsquareover2 * accel_elastic(i,1)
    veloc_elastic(i,1) = veloc_elastic(i,1) + deltatover2 * accel_elastic(i,1)
    accel_elastic(i,1) = 0._CUSTOM_REAL
  enddo
#else
  ! non-vectorized updates
  if (USE_ENFORCE_FIELDS) then
    do iglob = 1,nglob_elastic
      if (iglob_is_forced(iglob)) then
        call enforce_fields(iglob,it)
      else
        ! big loop over all the global points (not elements) in the mesh to update
        ! the displacement and velocity vectors and clear the acceleration vector
        displ_elastic(:,iglob) = displ_elastic(:,iglob) &
                                 + deltat * veloc_elastic(:,iglob) + deltatsquareover2 * accel_elastic(:,iglob)
        veloc_elastic(:,iglob) = veloc_elastic(:,iglob) + deltatover2 * accel_elastic(:,iglob) ! Predictor
        accel_elastic(:,iglob) = 0._CUSTOM_REAL
      endif
    enddo
  else
    displ_elastic(:,:) = displ_elastic(:,:) + deltat * veloc_elastic(:,:) + deltatsquareover2 * accel_elastic(:,:)
    veloc_elastic(:,:) = veloc_elastic(:,:) + deltatover2 * accel_elastic(:,:)
    accel_elastic(:,:) = 0._CUSTOM_REAL
  endif

#endif

  end subroutine update_displacement_newmark_elastic

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_newmark_poroelastic(deltat,deltatover2,deltatsquareover2, &
                                                     accels_poroelastic,velocs_poroelastic, &
                                                     displs_poroelastic,accelw_poroelastic, &
                                                     velocw_poroelastic,displw_poroelastic)

  use constants, only: CUSTOM_REAL,NDIM

  use specfem_par, only: nglob_poroelastic,PML_BOUNDARY_CONDITIONS

  implicit none

  double precision,intent(in) :: deltat,deltatover2,deltatsquareover2
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic),intent(inout) :: &
    accels_poroelastic,velocs_poroelastic,displs_poroelastic, &
    accelw_poroelastic,velocw_poroelastic,displw_poroelastic

  ! local parameters
#ifdef FORCE_VECTORIZATION
  integer :: i
#endif

  ! safety check
  !PML did not implemented for poroelastic simulation
  if (PML_BOUNDARY_CONDITIONS) stop 'Updating displacement for PML on poroelastic domain not implemented yet'

  !for the solid
#ifdef FORCE_VECTORIZATION
  do i = 1, NDIM * nglob_poroelastic
    displs_poroelastic(i,1) = displs_poroelastic(i,1) + deltat*velocs_poroelastic(i,1) &
                              + deltatsquareover2*accels_poroelastic(i,1)
    velocs_poroelastic(i,1) = velocs_poroelastic(i,1) + deltatover2*accels_poroelastic(i,1)
    accels_poroelastic(i,1) = 0._CUSTOM_REAL
  enddo
#else
  displs_poroelastic(:,:) = displs_poroelastic(:,:) + deltat*velocs_poroelastic(:,:) &
                            + deltatsquareover2*accels_poroelastic(:,:)
  velocs_poroelastic(:,:) = velocs_poroelastic(:,:) + deltatover2*accels_poroelastic(:,:)
  accels_poroelastic(:,:) = 0._CUSTOM_REAL
#endif

  !for the fluid
#ifdef FORCE_VECTORIZATION
  do i = 1, NDIM * nglob_poroelastic
    displw_poroelastic(i,1) = displw_poroelastic(i,1) + deltat*velocw_poroelastic(i,1) &
                              + deltatsquareover2*accelw_poroelastic(i,1)
    velocw_poroelastic(i,1) = velocw_poroelastic(i,1) + deltatover2*accelw_poroelastic(i,1)
    accelw_poroelastic(i,1) = 0._CUSTOM_REAL
  enddo
#else
  displw_poroelastic(:,:) = displw_poroelastic(:,:) + deltat*velocw_poroelastic(:,:) &
                            + deltatsquareover2*accelw_poroelastic(:,:)
  velocw_poroelastic(:,:) = velocw_poroelastic(:,:) + deltatover2*accelw_poroelastic(:,:)
  accelw_poroelastic(:,:) = 0._CUSTOM_REAL
#endif

  end subroutine update_displacement_newmark_poroelastic

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_newmark_GPU()

  use specfem_par, only: any_acoustic,any_elastic,any_poroelastic,myrank

  implicit none

  ! update displacement using finite-difference time scheme (Newmark)

  if (any_acoustic) then
    ! wavefields on GPU
    call update_displacement_newmark_GPU_acoustic()
  endif

  if (any_elastic) then
    ! wavefields on GPU
    call update_displacement_newmark_GPU_elastic()
  endif

  if (any_poroelastic) then
    ! safety stop
    call exit_MPI(myrank,'poroelastic time marching scheme on GPU not implemented yet...')
  endif

  end subroutine update_displacement_newmark_GPU

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_newmark_GPU_acoustic()

  use specfem_par, only: SIMULATION_TYPE,PML_BOUNDARY_CONDITIONS,myrank

  use specfem_par_gpu, only: Mesh_pointer,deltatf,deltatover2f,deltatsquareover2f,b_deltatf,b_deltatover2f, &
    b_deltatsquareover2f

  implicit none

  ! update displacement using finite-difference time scheme (Newmark)

  ! wavefields on GPU
  ! check
  if (SIMULATION_TYPE == 3) then
    if (PML_BOUNDARY_CONDITIONS) then
      call exit_MPI(myrank,'acoustic time marching scheme with PML_CONDITIONS on GPU not implemented yet...')
    endif
  endif

  ! updates acoustic potentials
  call update_displacement_ac_cuda(Mesh_pointer,deltatf,deltatsquareover2f,deltatover2f, &
                                   b_deltatf,b_deltatsquareover2f,b_deltatover2f)

  end subroutine update_displacement_newmark_GPU_acoustic

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_newmark_GPU_elastic()

  use specfem_par, only: SIMULATION_TYPE,PML_BOUNDARY_CONDITIONS,myrank

  use specfem_par_gpu, only: Mesh_pointer,deltatf,deltatover2f,deltatsquareover2f,b_deltatf,b_deltatover2f, &
    b_deltatsquareover2f

  implicit none

  ! update displacement using finite-difference time scheme (Newmark)

  ! wavefields on GPU
  ! check
  if (SIMULATION_TYPE == 3) then
    if (PML_BOUNDARY_CONDITIONS) then
      call exit_MPI(myrank,'elastic time marching scheme with PML_CONDITIONS on GPU not implemented yet...')
    endif
  endif

  ! updates elastic displacement and velocity
  ! Includes SIM_TYPE 1 & 3 (for noise tomography)
  call update_displacement_cuda(Mesh_pointer,deltatf,deltatsquareover2f,deltatover2f, &
                                b_deltatf,b_deltatsquareover2f,b_deltatover2f)

  end subroutine update_displacement_newmark_GPU_elastic


!------------------------------------------------------------------------------------------------
!
! corrector phase
!
!------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
!
! elastic domains
!
!------------------------------------------------------------------------------------------------

  subroutine update_veloc_elastic_Newmark()

! updates velocity vector (corrector)

  use constants, only: USE_ENFORCE_FIELDS
  use specfem_par

  implicit none

  ! local parameters
  integer :: iglob

  ! time marching
  !! DK DK this should be vectorized
  if (USE_ENFORCE_FIELDS) then
    do iglob = 1,nglob_elastic
      ! big loop over all the global points (not elements) in the mesh to update velocity vector (corrector).
      if (.not. iglob_is_forced(iglob)) then
        veloc_elastic(:,iglob) = veloc_elastic(:,iglob) + deltatover2 * accel_elastic(:,iglob)
      endif
    enddo
  else
    veloc_elastic(:,:) = veloc_elastic(:,:) + deltatover2 * accel_elastic(:,:)
  endif

  end subroutine update_veloc_elastic_Newmark

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_veloc_elastic_Newmark_backward()

! updates backward velocity vector (corrector)

  use constants, only: USE_ENFORCE_FIELDS
  use specfem_par

  implicit none

  ! local parameters
  integer :: iglob

  ! time marching
  !! DK DK this should be vectorized
  if (USE_ENFORCE_FIELDS) then
    do iglob = 1,nglob_elastic
      ! big loop over all the global points (not elements) in the mesh to update velocity vector (corrector).
      if (.not. iglob_is_forced(iglob)) then
        b_veloc_elastic(:,iglob) = b_veloc_elastic(:,iglob) + b_deltatover2 * b_accel_elastic(:,iglob)
      endif
    enddo
  else
    b_veloc_elastic(:,:) = b_veloc_elastic(:,:) + b_deltatover2 * b_accel_elastic(:,:)
  endif


  end subroutine update_veloc_elastic_Newmark_backward



!------------------------------------------------------------------------------------------------
!
! acoustic domains
!
!------------------------------------------------------------------------------------------------

  subroutine update_veloc_acoustic_Newmark()

! updates velocity potential (corrector)

  use specfem_par

  implicit none

  !! DK DK this should be vectorized
  potential_dot_acoustic(:) = potential_dot_acoustic(:) + deltatover2 * potential_dot_dot_acoustic(:)

  ! update the potential field (use a new array here) for coupling terms
  if (SIMULATION_TYPE == 3) then
    potential_acoustic_adj_coupling(:) = potential_acoustic(:) + deltat * potential_dot_acoustic(:) + &
                                         deltatsquareover2 * potential_dot_dot_acoustic(:)
  endif

  end subroutine update_veloc_acoustic_Newmark


!------------------------------------------------------------------------------------------------

  subroutine update_veloc_acoustic_Newmark_backward()

! updates velocity potential (corrector)

  use specfem_par

  !! DK DK this should be vectorized
  b_potential_dot_acoustic(:) = b_potential_dot_acoustic(:) + b_deltatover2 * b_potential_dot_dot_acoustic(:)

  end subroutine update_veloc_acoustic_Newmark_backward


!------------------------------------------------------------------------------------------------
!
! poroelastic domains
!
!------------------------------------------------------------------------------------------------

  subroutine update_veloc_poroelastic_Newmark()

! updates velocity (corrector)

  use specfem_par

  implicit none

  ! solid
  velocs_poroelastic(:,:) = velocs_poroelastic(:,:) + deltatover2 * accels_poroelastic(:,:)
  ! fluid
  velocw_poroelastic(:,:) = velocw_poroelastic(:,:) + deltatover2 * accelw_poroelastic(:,:)

  end subroutine update_veloc_poroelastic_Newmark


!------------------------------------------------------------------------------------------------

  subroutine update_veloc_poroelastic_Newmark_backward()

! updates velocity (corrector)

  use specfem_par

  ! solid
  b_velocs_poroelastic(:,:) = b_velocs_poroelastic(:,:) + b_deltatover2 * b_accels_poroelastic(:,:)
  ! fluid
  b_velocw_poroelastic(:,:) = b_velocw_poroelastic(:,:) + b_deltatover2 * b_accelw_poroelastic(:,:)

  end subroutine update_veloc_poroelastic_Newmark_backward

