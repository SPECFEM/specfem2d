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
!=====================================================================

! for viscoelastic solver

  subroutine compute_add_sources_viscoelastic(accel_elastic,it,i_stage)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM

  use specfem_par, only: P_SV,ispec_is_elastic,nglob_elastic, &
                         NSOURCES,source_time_function, &
                         is_proc_source,ispec_selected_source,sourcearrays, &
                         ibool
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic) :: accel_elastic
  integer :: it, i_stage

  !local variable
  integer :: i_source,i,j,iglob,ispec

  ! --- add the source
  do i_source = 1,NSOURCES

    ! if this processor core carries the source
    if (is_proc_source(i_source) == 1) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! source element is elastic
      if (ispec_is_elastic(ispec)) then

        ! adds source term
        ! note: we use sourcearrays for both collocated forces and moment tensors
        !       (see setup in setup_source_interpolation() routine)
        if (P_SV) then
          ! P-SV calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              accel_elastic(1,iglob) = accel_elastic(1,iglob) + &
                                       sourcearrays(i_source,1,i,j) * source_time_function(i_source,it,i_stage)
              accel_elastic(2,iglob) = accel_elastic(2,iglob) + &
                                       sourcearrays(i_source,1,i,j) * source_time_function(i_source,it,i_stage)
            enddo
          enddo
        else
          ! SH (membrane) calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              accel_elastic(1,iglob) = accel_elastic(1,iglob) + &
                                       sourcearrays(i_source,1,i,j) * source_time_function(i_source,it,i_stage)

              ! daniel debug 
              if (iglob == 37905) &
              write(1234,*) it, dble(sourcearrays(i_source,1,i,j) * source_time_function(i_source,it,i_stage)), &
                            accel_elastic(1,iglob),source_time_function(i_source,it,i_stage),sourcearrays(i_source,1,i,j)


            enddo
          enddo
        endif

      endif ! source element is elastic
    endif ! if this processor core carries the source
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_viscoelastic

!
!=====================================================================
!

! for viscoelastic solver for adjoint propagation wave field

  subroutine compute_add_sources_viscoelastic_adjoint()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: myrank,P_SV,accel_elastic,ispec_is_elastic,NSTEP,it,&
                         nrec,which_proc_receiver,ispec_selected_rec,adj_sourcearrays,&
                         ibool
  implicit none

  !local variables
  integer :: irec_local,irec,i,j,iglob,ispec
  integer :: it_tmp

  ! time step index
  it_tmp = NSTEP - it + 1

  irec_local = 0
  do irec = 1,nrec
    !   add the source (only if this proc carries the source)
    if (myrank == which_proc_receiver(irec)) then
      irec_local = irec_local + 1

      ! element containing adjoint source
      ispec = ispec_selected_rec(irec)

      if (ispec_is_elastic(ispec)) then
        ! add source array
        if (P_SV) then
          ! P-SH waves
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              accel_elastic(1,iglob) = accel_elastic(1,iglob) + adj_sourcearrays(irec_local,it_tmp,1,i,j)
              accel_elastic(2,iglob) = accel_elastic(2,iglob) + adj_sourcearrays(irec_local,it_tmp,2,i,j)
            enddo
          enddo
        else
          ! SH (membrane) wavescompute_forces_v
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              accel_elastic(1,iglob) = accel_elastic(1,iglob) + adj_sourcearrays(irec_local,it_tmp,1,i,j)
            enddo
          enddo
        endif
      endif ! if element is elastic
    endif ! if this processor core carries the adjoint source and the source element is elastic

  enddo ! irec = 1,nrec

  end subroutine compute_add_sources_viscoelastic_adjoint

