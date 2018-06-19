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

! for viscoelastic solver

  subroutine compute_add_sources_viscoelastic(accel_elastic,it,i_stage)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM

  use specfem_par, only: myrank,P_SV,ispec_is_elastic,nglob_elastic, &
                         NSOURCES,source_time_function, &
                         islice_selected_source,ispec_selected_source,sourcearrays, &
                         ibool
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic) :: accel_elastic
  integer :: it, i_stage

  !local variable
  integer :: i_source,i,j,iglob,ispec
  real(kind=CUSTOM_REAL) :: stf_used

  ! --- add the source
  do i_source = 1,NSOURCES

    ! if this processor core carries the source
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! source element is elastic
      if (ispec_is_elastic(ispec)) then

        ! source time function
        stf_used = source_time_function(i_source,it,i_stage)

        ! adds source term
        ! note: we use sourcearrays for both collocated forces and moment tensors
        !       (see setup in setup_source_interpolation() routine)
        if (P_SV) then
          ! P-SV calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              accel_elastic(1,iglob) = accel_elastic(1,iglob) + sourcearrays(i_source,1,i,j) * stf_used
              accel_elastic(2,iglob) = accel_elastic(2,iglob) + sourcearrays(i_source,2,i,j) * stf_used
            enddo
          enddo
        else
          ! SH (membrane) calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              accel_elastic(1,iglob) = accel_elastic(1,iglob) + sourcearrays(i_source,1,i,j) * stf_used

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

  subroutine compute_add_sources_viscoelastic_moving_source(accel_elastic,it,i_stage)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,SOURCE_IS_MOVING,TINYVAL,NGLJ,IMAIN

  use specfem_par, only: P_SV,ispec_is_elastic,nglob_elastic, &
                         NSOURCES,source_time_function, &
                         islice_selected_source,ispec_selected_source,sourcearrays, &
                         ibool,coord,nspec,nglob,xigll,zigll,z_source,NPROC, & !These 3 lines are for moving src
                         xi_source,gamma_source,coorg,knods,ngnod,npgeo,iglob_source,x_source,z_source,deltat,t0,myrank, &
                         time_stepping_scheme,hxis_store,hgammas_store,tshift_src,source_type,ispec_is_acoustic, &
                         hxis,hpxis,hgammas,hpgammas,anglesource,ispec_is_poroelastic,Mxx,Mxz,Mzz,gammax,gammaz,xix,xiz, &
                         AXISYM,xiglj,is_on_the_axis,initialfield
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic) :: accel_elastic
  integer :: it, i_stage

  !local variable
  integer :: i_source,i,j,iglob,ispec
  real(kind=CUSTOM_REAL) :: stf_used
  double precision :: hlagrange
  double precision :: xminSource,vSource,timeval,t_used
  ! single source array
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  ! checks if anything to do
  if (.not. SOURCE_IS_MOVING) return

  !xminSource = -5000.0d0 !m
  !vSource = 2150.0d0 !1425.0d0 !1250.0 !m/s
  xminSource = -60.0d0 !m
  vSource = 60.0d0 !m/s

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
      if (source_type(i_source) == 1) then
        call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                           x_source(i_source),z_source(i_source), &
                           ispec_selected_source(i_source),islice_selected_source(i_source), &
                           NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo, &
                           iglob_source(i_source),.true.)

      else if (source_type(i_source) == 2) then
        ! moment-tensor source
        call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                           x_source(i_source),z_source(i_source), &
                           ispec_selected_source(i_source),islice_selected_source(i_source), &
                           NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo, &
                           iglob_source(i_source),.false.)

      else if (.not. initialfield) then

        call exit_MPI(myrank,'incorrect source type')

      endif

      ispec = ispec_selected_source(i_source)
      ! source element is elastic
      if (ispec_is_elastic(ispec)) then

        ! Lagrange interpolators
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            call lagrange_any(xi_source(i_source),NGLJ,xiglj,hxis,hpxis)
            !do j = 1,NGLJ ! ABAB same result with that loop, this is good
            !  hxis(j) = hglj(j-1,xi_source(i_source),xiglj,NGLJ)
            !enddo
          else
            call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
          endif
        else
          call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
        endif
        call lagrange_any(gamma_source(i_source),NGLLZ,zigll,hgammas,hpgammas)

        if (mod(it,10000) == 0) then
            !  write(IMAIN,*) "myrank:",myrank
            ! user output
            if (myrank == islice_selected_source(i_source)) then
              iglob = ibool(2,2,ispec_selected_source(i_source))
              !write(IMAIN,*) 'xcoord: ',coord(1,iglob)
              write(IMAIN,*) 'Problem... it??: ',it,'xcoord: ',coord(1,iglob)," iglob",iglob
              !'source carried by proc',myrank,"  source x:",x_source(i_source)," ispec:",ispec_selected_source(i_source)

              !call flush_IMAIN()
            endif

        endif

        ! stores Lagrangians for source
        hxis_store(i_source,:) = hxis(:)
        hgammas_store(i_source,:) = hgammas(:)

        sourcearray(:,:,:) = 0._CUSTOM_REAL

        ! computes source arrays
        select case (source_type(i_source))
        case (1)
          ! collocated force source
          do j = 1,NGLLZ
            do i = 1,NGLLX
              hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)

              ! source element is acoustic
              if (ispec_is_acoustic(ispec)) then
                sourcearray(:,i,j) = real(hlagrange,kind=CUSTOM_REAL)
              endif

              ! source element is elastic
              if (ispec_is_elastic(ispec)) then
                if (P_SV) then
                  ! P_SV case
                  sourcearray(1,i,j) = real(- sin(anglesource(i_source)) * hlagrange,kind=CUSTOM_REAL)
                  sourcearray(2,i,j) = real(cos(anglesource(i_source)) * hlagrange,kind=CUSTOM_REAL)
                else
                  ! SH case (membrane)
                  sourcearray(:,i,j) = real(hlagrange,kind=CUSTOM_REAL)
                endif
              endif

              ! source element is poroelastic
              if (ispec_is_poroelastic(ispec)) then
                sourcearray(1,i,j) = real(- sin(anglesource(i_source)) * hlagrange,kind=CUSTOM_REAL)
                sourcearray(2,i,j) = real(cos(anglesource(i_source)) * hlagrange,kind=CUSTOM_REAL)
              endif

            enddo
          enddo

        case (2)
          ! moment-tensor source
          call compute_arrays_source(ispec,xi_source(i_source),gamma_source(i_source),sourcearray, &
                                     Mxx(i_source),Mzz(i_source),Mxz(i_source),xix,xiz,gammax,gammaz,xigll,zigll,nspec)
          ! checks source
          if (ispec_is_acoustic(ispec)) then
            call exit_MPI(myrank,'cannot have moment tensor source in acoustic element')
          endif

          ! checks wave type
          if (ispec_is_elastic(ispec)) then
            if (.not. P_SV ) call exit_MPI(myrank,'cannot have moment tensor source in SH (membrane) waves calculation')
          endif

        end select

        ! stores sourcearray for all sources
        sourcearrays(i_source,:,:,:) = sourcearray(:,:,:)

      endif
    endif
  enddo

  ! --- add the source
  do i_source = 1,NSOURCES

    ! if this processor core carries the source
    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! source element is elastic
      if (ispec_is_elastic(ispec)) then

        ! source time function
        stf_used = source_time_function(i_source,it,i_stage)

        ! adds source term
        ! note: we use sourcearrays for both collocated forces and moment tensors
        !       (see setup in setup_source_interpolation() routine)
        if (P_SV) then
          ! P-SV calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              accel_elastic(1,iglob) = accel_elastic(1,iglob) + &
                                       sourcearrays(i_source,1,i,j) * stf_used
              accel_elastic(2,iglob) = accel_elastic(2,iglob) + &
                                       sourcearrays(i_source,2,i,j) * stf_used
            enddo
          enddo
        else
          ! SH (membrane) calculation
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)

              accel_elastic(1,iglob) = accel_elastic(1,iglob) + &
                                       sourcearrays(i_source,1,i,j) * stf_used

              ! daniel debug source contribution
              !if (iglob == 37905) &
              !write(1234,*) it, dble(sourcearrays(i_source,1,i,j) * source_time_function(i_source,it,i_stage)), &
              !              accel_elastic(1,iglob),source_time_function(i_source,it,i_stage),sourcearrays(i_source,1,i,j)


            enddo
          enddo
        endif

      endif ! source element is elastic
    endif ! if this processor core carries the source
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_viscoelastic_moving_source

!
!=====================================================================
!

! for viscoelastic solver for adjoint propagation wave field

  subroutine compute_add_sources_viscoelastic_adjoint()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: P_SV,accel_elastic,ispec_is_elastic,NSTEP,it, &
                         nrecloc,ispec_selected_rec_loc,ibool, &
                         source_adjoint,xir_store_loc,gammar_store_loc
  implicit none

  !local variables
  integer :: irec_local,i,j,iglob,ispec
  integer :: it_tmp

  ! time step index
  it_tmp = NSTEP - it + 1

  do irec_local = 1,nrecloc

    ! element containing adjoint source
    ispec = ispec_selected_rec_loc(irec_local)

    if (ispec_is_elastic(ispec)) then
      ! add source array
      if (P_SV) then
        ! P-SV waves
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            accel_elastic(1,iglob) = accel_elastic(1,iglob) + real(xir_store_loc(irec_local,i)*gammar_store_loc(irec_local,j)* &
                                        source_adjoint(irec_local,it_tmp,1),kind=CUSTOM_REAL)
            accel_elastic(2,iglob) = accel_elastic(2,iglob) + real(xir_store_loc(irec_local,i)*gammar_store_loc(irec_local,j)* &
                                        source_adjoint(irec_local,it_tmp,2),kind=CUSTOM_REAL)
          enddo
        enddo
      else
        ! SH (membrane) wavescompute_forces_v
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            accel_elastic(1,iglob) = accel_elastic(1,iglob) +  real(xir_store_loc(irec_local,i)*gammar_store_loc(irec_local,j)* &
                                        source_adjoint(irec_local,it_tmp,1),kind=CUSTOM_REAL)
          enddo
        enddo
      endif
    endif ! if element is elastic

  enddo ! irec_local = 1,nrecloc

  end subroutine compute_add_sources_viscoelastic_adjoint

