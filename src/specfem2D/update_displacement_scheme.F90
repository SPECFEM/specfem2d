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
subroutine update_displacement_precondition_newmark_acoustic(deltat,deltatover2,deltatsquareover2,&
                                                             potential_dot_dot_acoustic,potential_dot_acoustic,&
                                                             potential_acoustic,potential_acoustic_old,&
                                                             PML_BOUNDARY_CONDITIONS)
  use specfem_par, only : nglob_acoustic
  implicit none
  include 'constants.h'

  double precision :: deltat,deltatover2,deltatsquareover2
  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: potential_acoustic,potential_dot_acoustic,&
                                                       potential_dot_dot_acoustic,potential_acoustic_old
  logical :: PML_BOUNDARY_CONDITIONS

#ifdef FORCE_VECTORIZATION
  integer :: i

  if( PML_BOUNDARY_CONDITIONS ) then
    do i = 1,nglob_acoustic
      potential_acoustic_old(i) = potential_acoustic(i) + deltatsquareover2/TWO*potential_dot_dot_acoustic(i)
    enddo
  endif

  do i = 1,nglob_acoustic
    potential_acoustic(i) = potential_acoustic(i) + deltat*potential_dot_acoustic(i) + &
                            deltatsquareover2*potential_dot_dot_acoustic(i)
    potential_dot_acoustic(i) = potential_dot_acoustic(i) + deltatover2*potential_dot_dot_acoustic(i)
    potential_dot_dot_acoustic(i) = ZERO
  enddo

#else

  if( PML_BOUNDARY_CONDITIONS ) then
    potential_acoustic_old = potential_acoustic + deltatsquareover2*potential_dot_dot_acoustic
  endif

  potential_acoustic = potential_acoustic + deltat*potential_dot_acoustic + deltatsquareover2*potential_dot_dot_acoustic
  potential_dot_acoustic = potential_dot_acoustic + deltatover2*potential_dot_dot_acoustic
  potential_dot_dot_acoustic = ZERO

#endif

end subroutine update_displacement_precondition_newmark_acoustic
!========================================================================

!========================================================================
subroutine update_displacement_precondition_newmark_elastic(deltat,deltatover2,deltatsquareover2,&
                                                            accel_elastic,veloc_elastic,&
                                                            displ_elastic,displ_elastic_old,&
                                                            PML_BOUNDARY_CONDITIONS)
  use specfem_par, only : nglob_elastic,ATTENUATION_VISCOELASTIC_SOLID
  implicit none
  include 'constants.h'

  double precision :: deltat,deltatover2,deltatsquareover2
  real(kind=CUSTOM_REAL), dimension(3,nglob_elastic) :: accel_elastic,veloc_elastic, &
                                                        displ_elastic,displ_elastic_old

  logical :: PML_BOUNDARY_CONDITIONS

#ifdef FORCE_VECTORIZATION
  integer :: i

  if( PML_BOUNDARY_CONDITIONS .or. ATTENUATION_VISCOELASTIC_SOLID ) then
    do i = 1,3*nglob_elastic
      displ_elastic_old(i,1) = displ_elastic(i,1) + deltatsquareover2/TWO * accel_elastic(i,1)
    enddo
  endif

  !! DK DK this allows for full vectorization by using a trick to see the 2D array as a 1D array
  !! DK DK whose dimension is the product of the two dimensions, the second dimension being equal to 1
  do i = 1,3*nglob_elastic !! DK DK here change 3 to NDIM when/if we suppress the 2nd component of the arrays (the SH component)
  !do i = 1,NDIM*nglob_elastic  !! DK DK this should be the correct size in principle, but not here because of the SH component
    displ_elastic(i,1) = displ_elastic(i,1) + deltat*veloc_elastic(i,1) + deltatsquareover2*accel_elastic(i,1)
    veloc_elastic(i,1) = veloc_elastic(i,1) + deltatover2*accel_elastic(i,1)
    accel_elastic(i,1) = ZERO
  enddo
#else

  if( PML_BOUNDARY_CONDITIONS .or. ATTENUATION_VISCOELASTIC_SOLID  ) then
    displ_elastic_old = displ_elastic + deltatsquareover2/TWO * accel_elastic
  endif

  displ_elastic = displ_elastic + deltat*veloc_elastic + deltatsquareover2*accel_elastic
  veloc_elastic = veloc_elastic + deltatover2*accel_elastic
  accel_elastic = ZERO

#endif

end subroutine update_displacement_precondition_newmark_elastic
!========================================================================

!========================================================================
subroutine update_displacement_precondition_newmark_poroelastic(deltat,deltatover2,deltatsquareover2,&
                                                                accels_poroelastic,velocs_poroelastic,&
                                                                displs_poroelastic,accelw_poroelastic,&
                                                                velocw_poroelastic,displw_poroelastic)

  use specfem_par, only : nglob_poroelastic
  implicit none
  include 'constants.h'

  double precision :: deltat,deltatover2,deltatsquareover2
  real(kind=CUSTOM_REAL), dimension(3,nglob_poroelastic) :: accels_poroelastic,velocs_poroelastic,displs_poroelastic,&
                                                            accelw_poroelastic,velocw_poroelastic,displw_poroelastic

  !PML did not implemented for poroelastic simulation

  !for the solid
  displs_poroelastic = displs_poroelastic + deltat*velocs_poroelastic + deltatsquareover2*accels_poroelastic
  velocs_poroelastic = velocs_poroelastic + deltatover2*accels_poroelastic
  accels_poroelastic = ZERO

  !for the fluid
  displw_poroelastic = displw_poroelastic + deltat*velocw_poroelastic + deltatsquareover2*accelw_poroelastic
  velocw_poroelastic = velocw_poroelastic + deltatover2*accelw_poroelastic
  accelw_poroelastic = ZERO

end subroutine update_displacement_precondition_newmark_poroelastic
!========================================================================


  subroutine update_displacement_precondition_newmark_GPU()

  use specfem_par, only : Mesh_pointer,deltatf,deltatover2f,deltatsquareover2f,b_deltatf,b_deltatover2f,&
                          b_deltatsquareover2f,SIMULATION_TYPE,&
                          any_acoustic,any_elastic,any_poroelastic,&
                          PML_BOUNDARY_CONDITIONS

  implicit none
  include 'constants.h'


! update displacement using finite-difference time scheme (Newmark)

  if(any_acoustic) then

    ! wavefields on GPU
    ! check
    if( SIMULATION_TYPE == 3 ) then
      if( PML_BOUNDARY_CONDITIONS )then
        call exit_MPI('acoustic time marching scheme with PML_CONDITIONS on GPU not implemented yet...')
      endif
    endif

    ! updates acoustic potentials
    call it_update_displacement_ac_cuda(Mesh_pointer,deltatf,deltatsquareover2f,deltatover2f,&
             b_deltatf,b_deltatsquareover2f,b_deltatover2f)

  endif

  if(any_elastic) then
    ! wavefields on GPU
    ! check
    if( SIMULATION_TYPE == 3 ) then
      if( PML_BOUNDARY_CONDITIONS )then

          call exit_MPI('elastic time marching scheme with PML_CONDITIONS on GPU not implemented yet...')

      endif
    endif
    ! updates elastic displacement and velocity
    ! Includes SIM_TYPE 1 & 3 (for noise tomography)
    call update_displacement_cuda(Mesh_pointer,deltatf,deltatsquareover2f,&
                                  deltatover2f,b_deltatf,b_deltatsquareover2f,b_deltatover2f)

  endif

  if(any_poroelastic) then
    call exit_MPI('poroelastic time marching scheme on GPU not implemented yet...')
  endif

  end subroutine update_displacement_precondition_newmark_GPU
