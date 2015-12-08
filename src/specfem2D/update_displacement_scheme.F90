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


  subroutine update_displacement_newmark_acoustic(deltat,deltatover2,deltatsquareover2,&
                                                  potential_dot_dot_acoustic,potential_dot_acoustic,&
                                                  potential_acoustic,potential_acoustic_old,&
                                                  PML_BOUNDARY_CONDITIONS)

  use specfem_par, only : nglob_acoustic,CUSTOM_REAL

  implicit none

  double precision :: deltat,deltatover2,deltatsquareover2
  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: potential_acoustic,potential_dot_acoustic,&
                                                       potential_dot_dot_acoustic,potential_acoustic_old
  logical :: PML_BOUNDARY_CONDITIONS

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

  subroutine update_displacement_newmark_elastic(deltat,deltatover2,deltatsquareover2,&
                                                 accel_elastic,veloc_elastic,&
                                                 displ_elastic,displ_elastic_old,&
                                                 PML_BOUNDARY_CONDITIONS)

  use specfem_par, only : nglob_elastic,ATTENUATION_VISCOELASTIC_SOLID,CUSTOM_REAL,TWO

  implicit none

  double precision :: deltat,deltatover2,deltatsquareover2
  real(kind=CUSTOM_REAL), dimension(3,nglob_elastic) :: accel_elastic,veloc_elastic, &
                                                        displ_elastic,displ_elastic_old

  logical :: PML_BOUNDARY_CONDITIONS

#ifdef FORCE_VECTORIZATION
  integer :: i
#endif

  ! attenuation/PML simulations
  if (PML_BOUNDARY_CONDITIONS .or. ATTENUATION_VISCOELASTIC_SOLID) then

    ! note: todo - there is an additional factor 1/TWO to the default deltasquareover2 for the acceleration term
    !       find explanations where?
#ifdef FORCE_VECTORIZATION
    do i = 1,3*nglob_elastic
      displ_elastic_old(i,1) = displ_elastic(i,1) + deltatsquareover2/TWO * accel_elastic(i,1)
    enddo
#else
    displ_elastic_old(:,:) = displ_elastic(:,:) + deltatsquareover2/TWO * accel_elastic(:,:)
#endif

  endif ! PM

#ifdef FORCE_VECTORIZATION
  !! DK DK this allows for full vectorization by using a trick to see the 2D array as a 1D array
  !! DK DK whose dimension is the product of the two dimensions, the second dimension being equal to 1
  do i = 1,3*nglob_elastic !! DK DK here change 3 to NDIM when/if we suppress the 2nd component of the arrays (the SH component)
  !do i = 1,NDIM*nglob_elastic  !! DK DK this should be the correct size in principle, but not here because of the SH component
    displ_elastic(i,1) = displ_elastic(i,1) + deltat * veloc_elastic(i,1) + deltatsquareover2 * accel_elastic(i,1)
    veloc_elastic(i,1) = veloc_elastic(i,1) + deltatover2 * accel_elastic(i,1)
    accel_elastic(i,1) = 0._CUSTOM_REAL
  enddo
#else
  displ_elastic(:,:) = displ_elastic(:,:) + deltat * veloc_elastic(:,:) + deltatsquareover2 * accel_elastic(:,:)
  veloc_elastic(:,:) = veloc_elastic(:,:) + deltatover2 * accel_elastic(:,:)
  accel_elastic(:,:) = 0._CUSTOM_REAL
#endif

  end subroutine update_displacement_newmark_elastic

!
!------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_newmark_poroelastic(deltat,deltatover2,deltatsquareover2,&
                                                     accels_poroelastic,velocs_poroelastic,&
                                                     displs_poroelastic,accelw_poroelastic,&
                                                     velocw_poroelastic,displw_poroelastic)

  use specfem_par, only : nglob_poroelastic,CUSTOM_REAL,NDIM,PML_BOUNDARY_CONDITIONS


  implicit none

  double precision :: deltat,deltatover2,deltatsquareover2
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic) :: accels_poroelastic,velocs_poroelastic,displs_poroelastic,&
                                                            accelw_poroelastic,velocw_poroelastic,displw_poroelastic

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

  use specfem_par, only : SIMULATION_TYPE,any_acoustic,any_elastic,any_poroelastic,&
    PML_BOUNDARY_CONDITIONS,myrank

  use specfem_par_gpu, only : Mesh_pointer,deltatf,deltatover2f,deltatsquareover2f,b_deltatf,b_deltatover2f,&
    b_deltatsquareover2f

  implicit none

  ! update displacement using finite-difference time scheme (Newmark)

  if (any_acoustic) then
    ! wavefields on GPU
    ! check
    if (SIMULATION_TYPE == 3) then
      if (PML_BOUNDARY_CONDITIONS) then
        call exit_MPI(myrank,'acoustic time marching scheme with PML_CONDITIONS on GPU not implemented yet...')
      endif
    endif

    ! updates acoustic potentials
    call it_update_displacement_ac_cuda(Mesh_pointer,deltatf,deltatsquareover2f,deltatover2f,&
                                        b_deltatf,b_deltatsquareover2f,b_deltatover2f)
  endif

  if (any_elastic) then
    ! wavefields on GPU
    ! check
    if (SIMULATION_TYPE == 3) then
      if (PML_BOUNDARY_CONDITIONS) then
        call exit_MPI(myrank,'elastic time marching scheme with PML_CONDITIONS on GPU not implemented yet...')
      endif
    endif

    ! updates elastic displacement and velocity
    ! Includes SIM_TYPE 1 & 3 (for noise tomography)
    call update_displacement_cuda(Mesh_pointer,deltatf,deltatsquareover2f,deltatover2f,&
                                  b_deltatf,b_deltatsquareover2f,b_deltatover2f)
  endif

  if (any_poroelastic) then
    call exit_MPI(myrank,'poroelastic time marching scheme on GPU not implemented yet...')
  endif

  end subroutine update_displacement_newmark_GPU

