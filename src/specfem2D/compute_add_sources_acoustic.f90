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
!=====================================================================

! for acoustic solver

  subroutine compute_add_sources_acoustic(potential_dot_dot_acoustic,it,i_stage)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: ispec_is_acoustic,nglob_acoustic, &
                         NSOURCES,source_type,source_time_function,sourcearrays, &
                         islice_selected_source,ispec_selected_source,ibool,kappastore,myrank
  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic),intent(inout) :: potential_dot_dot_acoustic
  integer,intent(in) :: it,i_stage

  !local variables
  integer :: i_source,i,j,iglob,ispec
  real(kind=CUSTOM_REAL) :: stf_used

  do i_source = 1,NSOURCES
    ! if this processor core carries the source
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! source element is acoustic
      if (ispec_is_acoustic(ispec)) then

        ! source time function
        stf_used = source_time_function(i_source,it,i_stage)

        ! collocated force
        ! beware, for an acoustic medium, the source is pressure divided by Kappa of the fluid
        if (source_type(i_source) == 1) then
          ! forward wavefield
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + &
                                  real(sourcearrays(i_source,1,i,j) * stf_used / kappastore(i,j,ispec),kind=CUSTOM_REAL)
            enddo
          enddo

        else if (source_type(i_source) == 2) then
          ! moment tensor
          call exit_MPI(myrank,'cannot have moment tensor source in acoustic element')
        endif

      endif
    endif ! if this processor core carries the source and the source element is acoustic
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_acoustic

!
!=====================================================================
!

  subroutine compute_add_sources_acoustic_moving_source(potential_dot_dot_acoustic,it,i_stage)

! This subroutine is the same than the previous one but with a moving source

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,TINYVAL,SOURCE_IS_MOVING,IMAIN

  use specfem_par, only: ispec_is_acoustic,nglob_acoustic, &
                         NSOURCES,source_type,source_time_function, &
                         islice_selected_source,ispec_selected_source, &
                         hxis_store,hgammas_store,ibool,kappastore,myrank,deltat,t0,tshift_src, &
! These three lines are for a moving acoustic source
                         coord,nspec,nglob,xigll,zigll,z_source,NPROC,xi_source, &
                         gamma_source,coorg,knods,ngnod,npgeo,iglob_source,x_source,z_source, &
                         time_stepping_scheme, &
                         hxis,hpxis,hgammas,hpgammas
  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic),intent(inout) :: potential_dot_dot_acoustic
  integer,intent(in) :: it,i_stage

  !local variables
  integer :: i_source,i,j,iglob,ispec
  double precision :: hlagrange
  double precision :: xminSource,vSource,timeval,t_used

  ! checks if anything to do
  if (.not. SOURCE_IS_MOVING) return

  xminSource = -60.0d0 ! m
  vSource = 60.0d0 ! m/s

  if (time_stepping_scheme == 1) then
    ! Newmark
    timeval = (it-1)*deltat
  else
    call exit_MPI(myrank,'Not implemented!')
  endif

  ! moves and re-locates sources along x-axis
  do i_source = 1,NSOURCES
    if (abs(source_time_function(i_source,it,i_stage)) > TINYVAL) then
      t_used = (timeval-t0-tshift_src(i_source))

      x_source(i_source) = xminSource + vSource*t_used !timeval?

      ! collocated force source
      call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                         x_source(i_source),z_source(i_source), &
                         ispec_selected_source(i_source),islice_selected_source(i_source), &
                         NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo, &
                         iglob_source(i_source),.true.)

      ! define and store Lagrange interpolators (hxis,hpxis,hgammas,hpgammas) at all the sources
      !if (AXISYM) then
      !  if (is_on_the_axis(ispec_selected_source(i_source)) .and. myrank == islice_selected_source(i_source)) then
      !    call lagrange_any(xi_source(i_source),NGLJ,xiglj,hxis,hpxis)
      !    !do j = 1,NGLJ ! ABAB same result with that loop, this is good
      !    !  hxis(j) = hglj(j-1,xi_source(i),xiglj,NGLJ)
      !    !enddo
      !  else
      !    call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
      !  endif
      !else
        call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
      !endif
      call lagrange_any(gamma_source(i_source),NGLLZ,zigll,hgammas,hpgammas)

      ! stores Lagrangians for source
      hxis_store(i_source,:) = hxis(:)
      hgammas_store(i_source,:) = hgammas(:)

    endif
  enddo

  ! adds source contributions
  do i_source = 1,NSOURCES
    ! if this processor core carries the source and the source element is acoustic
    ! .and. acoustic(ispec_selected_source(i_source)) ??
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      if (ispec_is_acoustic(ispec)) then
        ! collocated force
        ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
        ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
        ! to add minus the source to Chi_dot_dot to get plus the source in pressure
        if (source_type(i_source) == 1) then
          ! forward wavefield
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)

              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + &
                      real(source_time_function(i_source,it,i_stage)*hlagrange / kappastore(i,j,ispec),kind=CUSTOM_REAL)
            enddo
          enddo
          ! moment tensor
          else if (source_type(i_source) == 2) then
            call exit_MPI(myrank,'cannot have moment tensor source in acoustic element')
        endif
      endif
    endif ! if this processor core carries the source and the source element is acoustic
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_acoustic_moving_source

!
!=====================================================================
!

! for acoustic solver for adjoint propagation wave field

  subroutine compute_add_sources_acoustic_adjoint()

  use constants, only: NGLLX,NGLLZ,CUSTOM_REAL

  use specfem_par, only: potential_dot_dot_acoustic,ispec_is_acoustic,NSTEP,it, &
                         nrecloc,ispec_selected_rec_loc, &
                         ibool,kappastore,source_adjoint,xir_store_loc,gammar_store_loc
  implicit none

  !local variables
  integer :: irec_local,i,j,iglob,ispec
  integer :: it_tmp

  ! time step index
  it_tmp = NSTEP - it + 1

  do irec_local = 1,nrecloc

    ! element containing adjoint source
    ispec = ispec_selected_rec_loc(irec_local)

    if (ispec_is_acoustic(ispec)) then
      ! add source array
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)

          !ZN Zinan Xie, Year 2012:
          !ZN be careful, the following line is newly added, thus when doing a comparison
          !ZN of the new code with the old code, you will have a big difference if you
          !ZN do not adjust the amplitude of the adjoint source
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
                                              real(xir_store_loc(irec_local,i)*gammar_store_loc(irec_local,j)* &
                                              source_adjoint(irec_local,it_tmp,1),kind=CUSTOM_REAL) &
                                              / kappastore(i,j,ispec)
        enddo
      enddo
    endif ! if element acoustic
  enddo ! irec_local = 1,nrecloc

  end subroutine compute_add_sources_acoustic_adjoint

