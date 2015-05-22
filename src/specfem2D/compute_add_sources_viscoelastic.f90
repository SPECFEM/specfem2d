
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
!=====================================================================

! for viscoelastic solver

  subroutine compute_add_sources_viscoelastic(accel_elastic,it,i_stage)

  use specfem_par, only: p_sv,elastic,nglob_elastic,&
                         NSOURCES,source_type,anglesource,source_time_function,&
                         is_proc_source,ispec_selected_source,sourcearray,&
                         hxis_store,hgammas_store,ibool
  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL), dimension(3,nglob_elastic) :: accel_elastic
  integer :: it, i_stage

  !local variable
  integer :: i_source,i,j,iglob
  double precision :: hlagrange

  ! --- add the source
  do i_source=1,NSOURCES
    ! if this processor core carries the source and the source element is elastic
    if( is_proc_source(i_source) == 1 .and. elastic(ispec_selected_source(i_source)) ) then
      ! collocated force
      if( source_type(i_source) == 1 ) then
        if( p_sv ) then ! P-SV calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec_selected_source(i_source))
              hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
              accel_elastic(1,iglob) = accel_elastic(1,iglob) - &
                                       sin(anglesource(i_source))*source_time_function(i_source,it,i_stage)*hlagrange
              accel_elastic(3,iglob) = accel_elastic(3,iglob) + &
                                       cos(anglesource(i_source))*source_time_function(i_source,it,i_stage)*hlagrange
            enddo
          enddo
        else    ! SH (membrane) calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec_selected_source(i_source))
              hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
              accel_elastic(2,iglob) = accel_elastic(2,iglob) + source_time_function(i_source,it,i_stage)*hlagrange
            enddo
          enddo
        endif
      endif

      ! moment tensor
      if( source_type(i_source) == 2 ) then
        if( .not. p_sv )  call exit_MPI('cannot have moment tensor source in SH (membrane) waves calculation')
        ! add source array
        do j=1,NGLLZ;
          do i=1,NGLLX
            iglob = ibool(i,j,ispec_selected_source(i_source))
            accel_elastic(1,iglob) = accel_elastic(1,iglob) + &
                                     sourcearray(i_source,1,i,j) * source_time_function(i_source,it,i_stage)
            accel_elastic(3,iglob) = accel_elastic(3,iglob) + &
                                     sourcearray(i_source,2,i,j) * source_time_function(i_source,it,i_stage)
          enddo
        enddo
      endif !if( source_type(i_source) == 2)
    endif ! if this processor core carries the source and the source element is elastic
  enddo ! do i_source=1,NSOURCES

  end subroutine compute_add_sources_viscoelastic
!
!=====================================================================
! for viscoelastic solver for adjoint propagation wave field
  subroutine compute_add_sources_viscoelastic_adjoint()

  use specfem_par, only: myrank,p_sv,accel_elastic,elastic,NSTEP,it,&
                         nrec,which_proc_receiver,ispec_selected_rec,adj_sourcearrays,&
                         ibool
  implicit none
  include "constants.h"

  !local variables
  integer :: irec_local,irec,i,j,iglob

  irec_local = 0
  do irec = 1,nrec
    !   add the source (only if this proc carries the source)
    if( myrank == which_proc_receiver(irec) ) then
      irec_local = irec_local + 1
      if( elastic(ispec_selected_rec(irec)) ) then
        ! add source array
        do j=1,NGLLZ
          do i=1,NGLLX
            iglob = ibool(i,j,ispec_selected_rec(irec))
            if( p_sv ) then !P-SH waves
              accel_elastic(1,iglob) = accel_elastic(1,iglob) + adj_sourcearrays(irec_local,NSTEP-it+1,1,i,j)
              accel_elastic(3,iglob) = accel_elastic(3,iglob) + adj_sourcearrays(irec_local,NSTEP-it+1,3,i,j)
            else !SH (membrane) wavescompute_forces_v
              accel_elastic(2,iglob) = accel_elastic(2,iglob) + adj_sourcearrays(irec_local,NSTEP-it+1,2,i,j)
            endif
          enddo
        enddo
      endif ! if element is elastic
    endif ! if this processor core carries the adjoint source and the source element is elastic

  enddo ! irec = 1,nrec

  end subroutine compute_add_sources_viscoelastic_adjoint

