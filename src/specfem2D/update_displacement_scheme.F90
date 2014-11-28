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


  subroutine update_displacement_precondition_newmark()

  use specfem_par, only : deltat,deltatover2,deltatsquareover2,b_deltat,b_deltatover2,b_deltatsquareover2, &
                          myrank,time_stepping_scheme,SIMULATION_TYPE,&
                          nglob_acoustic,nglob_elastic,nglob_poroelastic,&
                          any_acoustic,any_elastic,any_poroelastic,&
                          potential_acoustic,potential_dot_acoustic,&
                          potential_dot_dot_acoustic,potential_acoustic_old,&
                          displ_elastic,veloc_elastic,accel_elastic,displ_elastic_old,&
                          displs_poroelastic,velocs_poroelastic,accels_poroelastic,&
                          displs_poroelastic_old,displw_poroelastic,velocw_poroelastic,&
                          accelw_poroelastic,&
                          b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic,&
                          b_potential_acoustic_old,&
                          b_displ_elastic,b_veloc_elastic,b_accel_elastic,b_displ_elastic_old,&
                          accel_elastic_adj_coupling,&
                          b_displs_poroelastic,b_velocs_poroelastic,b_accels_poroelastic,&
                          accels_poroelastic_adj_coupling,&
                          b_displw_poroelastic,b_velocw_poroelastic,b_accelw_poroelastic,&
                          accelw_poroelastic_adj_coupling,PML_BOUNDARY_CONDITIONS

  implicit none
  include 'constants.h'


! update displacement using finite-difference time scheme (Newmark)

    if(any_acoustic) then

      if(time_stepping_scheme==1)then
      ! Newmark time scheme
!! DK DK this should be vectorized
        potential_acoustic_old = potential_acoustic + deltatsquareover2*potential_dot_dot_acoustic
        potential_acoustic = potential_acoustic + &
                 deltat*potential_dot_acoustic + deltatsquareover2*potential_dot_dot_acoustic
        potential_dot_acoustic = potential_dot_acoustic + deltatover2*potential_dot_dot_acoustic

      endif
      potential_dot_dot_acoustic = ZERO

      if(SIMULATION_TYPE == 3) then ! Adjoint calculation
!! DK DK this should be vectorized
        b_potential_acoustic_old = b_potential_acoustic + b_deltatsquareover2*b_potential_dot_dot_acoustic
        b_potential_acoustic = b_potential_acoustic + b_deltat*b_potential_dot_acoustic + &
                               b_deltatsquareover2*b_potential_dot_dot_acoustic
        b_potential_dot_acoustic = b_potential_dot_acoustic &
                                  + b_deltatover2*b_potential_dot_dot_acoustic
        b_potential_dot_dot_acoustic = ZERO
      endif

    endif

    if(any_elastic) then

      if(time_stepping_scheme==1)then
#ifdef FORCE_VECTORIZATION
!! DK DK this allows for full vectorization by using a trick to see the 2D array as a 1D array
!! DK DK whose dimension is the product of the two dimensions, the second dimension being equal to 1
     do i = 1,3*nglob_elastic !! DK DK here change 3 to NDIM when/if we suppress the 2nd component of the arrays (the SH component)
!    do i = 1,NDIM*nglob_elastic  !! DK DK this should be the correct size in principle, but not here because of the SH component
      displ_elastic_old(i,1) = displ_elastic(i,1) + deltatsquareover2 * accel_elastic(i,1)
      displ_elastic(i,1) = displ_elastic(i,1) &
                    + deltat*veloc_elastic(i,1) &
                    + deltatsquareover2*accel_elastic(i,1)
      veloc_elastic(i,1) = veloc_elastic(i,1) + deltatover2*accel_elastic(i,1)
     enddo
#else
      displ_elastic_old = displ_elastic + deltatsquareover2 * accel_elastic
      displ_elastic = displ_elastic &
                    + deltat*veloc_elastic &
                    + deltatsquareover2*accel_elastic
      veloc_elastic = veloc_elastic + deltatover2*accel_elastic
#endif
      endif
      accel_elastic_adj_coupling = accel_elastic
      accel_elastic = ZERO

      if(SIMULATION_TYPE == 3) then ! Adjoint calculation
!! DK DK this should be fully vectorized
        b_displ_elastic_old = b_displ_elastic + b_deltatsquareover2 * b_accel_elastic
        b_displ_elastic = b_displ_elastic &
                        + b_deltat*b_veloc_elastic &
                        + b_deltatsquareover2*b_accel_elastic
        b_veloc_elastic = b_veloc_elastic + b_deltatover2*b_accel_elastic
        b_accel_elastic = ZERO
      endif
    endif

    if(any_poroelastic) then

      if(time_stepping_scheme==1)then
      !for the solid
      displs_poroelastic_old = displs_poroelastic + deltatover2 * accels_poroelastic
      displs_poroelastic = displs_poroelastic &
                         + deltat*velocs_poroelastic &
                         + deltatsquareover2*accels_poroelastic
      velocs_poroelastic = velocs_poroelastic + deltatover2*accels_poroelastic
      accels_poroelastic_adj_coupling = accels_poroelastic
      accels_poroelastic = ZERO
      !for the fluid
      displw_poroelastic = displw_poroelastic &
                         + deltat*velocw_poroelastic &
                         + deltatsquareover2*accelw_poroelastic
      velocw_poroelastic = velocw_poroelastic + deltatover2*accelw_poroelastic
      accelw_poroelastic_adj_coupling = accelw_poroelastic
      accelw_poroelastic = ZERO
      endif

      if(time_stepping_scheme==2)then
      !for the solid
      accels_poroelastic_adj_coupling = accels_poroelastic
      accels_poroelastic = ZERO
      !for the fluid
      accelw_poroelastic_adj_coupling = accelw_poroelastic
      accelw_poroelastic = ZERO
      endif

      if(time_stepping_scheme==3)then
      !for the solid
      accels_poroelastic_adj_coupling = accels_poroelastic
      accels_poroelastic = ZERO
      !for the fluid
      accelw_poroelastic_adj_coupling = accelw_poroelastic
      accelw_poroelastic = ZERO
      endif

      if(SIMULATION_TYPE == 3) then ! Adjoint calculation
        !for the solid
        b_displs_poroelastic = b_displs_poroelastic &
                             + b_deltat*b_velocs_poroelastic &
                             + b_deltatsquareover2*b_accels_poroelastic
        b_velocs_poroelastic = b_velocs_poroelastic + b_deltatover2*b_accels_poroelastic
        b_accels_poroelastic = ZERO
        !for the fluid
        b_displw_poroelastic = b_displw_poroelastic &
                             + b_deltat*b_velocw_poroelastic &
                             + b_deltatsquareover2*b_accelw_poroelastic
        b_velocw_poroelastic = b_velocw_poroelastic + b_deltatover2*b_accelw_poroelastic
        b_accelw_poroelastic = ZERO
      endif
    endif

  end subroutine update_displacement_precondition_newmark
