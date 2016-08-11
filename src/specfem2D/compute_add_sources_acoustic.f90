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
        ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
        ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
        ! to add minus the source to Chi_dot_dot to get plus the source in pressure
        if (source_type(i_source) == 1) then
          ! forward wavefield
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              ! old way: source without factor 1/kappa
              !potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
              !                                    sourcearrays(i_source,1,i,j) * stf_used

              !ZN becareful the following line is new added, thus when do comparison
              !ZN of the new code with the old code, you will have big difference if you
              !ZN do not tune the source
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
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
                         coord,nspec,nglob,xigll,zigll,z_source,NPROC,xi_source,& !These 3 lines are for moving src
                         gamma_source,coorg,knods,ngnod,npgeo,iglob_source,x_source,z_source, &
                         time_stepping_scheme, &
                         hxis,hpxis,hgammas,hpgammas !,AXISYM,xiglj,is_on_the_axis
  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic),intent(inout) :: potential_dot_dot_acoustic
  integer,intent(in) :: it,i_stage

  !local variables
  integer :: i_source,i,j,iglob,ispec
  double precision :: hlagrange
  double precision :: xminSource,vSource,timeval,t_used

  ! checks if anything to do
  if (.not. SOURCE_IS_MOVING) return

  xminSource = -15000.0d0 !m
  vSource = 1250.0d0 !1425.0d0 !1250.0 !m/s

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

!      if (mod(it,10) == 0) then
!          !  write(IMAIN,*) "myrank:",myrank
!          ! user output
!          if (myrank == islice_selected_source(i_source)) then
!            iglob = ibool(2,2,ispec_selected_source(i_source))
!            !write(IMAIN,*) 'xcoord: ',coord(1,iglob)
!            write(IMAIN,*) 'it?: ',it,'xcoord: ',coord(1,iglob)," iglob",iglob
!            !'source carried by proc',myrank,"  source x:",x_source(i_source)," ispec:",ispec_selected_source(i_source)

!            !call flush_IMAIN()
!          endif

!      endif
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

              !if (mod(it,10) == 0 .and. i == 2 .and. j == 2) write(IMAIN,*) 'it',it,'source carried by proc',myrank, &
              !"iglob",iglob !"  source x:",x_source(i_source)," xcoord:", coord(1,iglob)," ispec:",ispec_selected_source(i_source)
              hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)

              !ZN becareful the following line is new added, thus when do comparison
              !ZN of the new code with the old code, you will have big difference if you
              !ZN do not tune the source
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
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

  use constants, only: NGLLX,NGLLZ

  use specfem_par, only: myrank,potential_dot_dot_acoustic,ispec_is_acoustic,NSTEP,it, &
                         nrec,islice_selected_rec,ispec_selected_rec,adj_sourcearrays, &
                         ibool,kappastore
  implicit none

  !local variables
  integer :: irec_local,irec,i,j,iglob,ispec
  integer :: it_tmp

  ! time step index
  it_tmp = NSTEP - it + 1

  irec_local = 0
  do irec = 1,nrec
    ! add the source (only if this proc carries the source)
    if (myrank == islice_selected_rec(irec)) then
      irec_local = irec_local + 1

      ! element containing adjoint source
      ispec = ispec_selected_rec(irec)

      if (ispec_is_acoustic(ispec)) then
        ! add source array
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)

            !ZN becareful the following line is new added, thus when do comparison
            !ZN of the new code with the old code, you will have big difference if you
            !ZN do not tune the source
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + &
                                                adj_sourcearrays(irec_local,it_tmp,1,i,j) &
                                                / kappastore(i,j,ispec)
          enddo
        enddo
      endif ! if element acoustic
    endif ! if this processor core carries the adjoint source
  enddo ! irec = 1,nrec

  end subroutine compute_add_sources_acoustic_adjoint

