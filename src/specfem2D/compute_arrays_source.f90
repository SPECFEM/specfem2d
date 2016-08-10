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

  subroutine compute_arrays_source(ispec_selected_source,xi_source,gamma_source,sourcearray, &
                                   Mxx,Mzz,Mxz,xix,xiz,gammax,gammaz,xigll,zigll,nspec)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM,ZERO

  use specfem_par, only: AXISYM,is_on_the_axis,xiglj

  implicit none

  integer,intent(in) :: ispec_selected_source
  integer,intent(in) :: nspec

  double precision,intent(in) :: xi_source,gamma_source
  double precision,intent(in) :: Mxx,Mzz,Mxz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  double precision :: xixd,xizd,gammaxd,gammazd

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll

! source arrays
  double precision, dimension(NGLLX,NGLLZ) :: G11,G13,G31,G33
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  integer :: k,m
  integer :: ir,iv

! calculate G_ij for general source location
! the source does not necessarily correspond to a Gauss-Lobatto point
  do m = 1,NGLLZ
      do k = 1,NGLLX

        xixd    = xix(k,m,ispec_selected_source)
        xizd    = xiz(k,m,ispec_selected_source)
        gammaxd = gammax(k,m,ispec_selected_source)
        gammazd = gammaz(k,m,ispec_selected_source)

        G11(k,m) = Mxx*xixd + Mxz*xizd
        G13(k,m) = Mxx*gammaxd + Mxz*gammazd
        G31(k,m) = Mxz*xixd + Mzz*xizd
        G33(k,m) = Mxz*gammaxd + Mzz*gammazd

      enddo
  enddo

! compute Lagrange polynomials at the source location

  if (AXISYM) then
    if (is_on_the_axis(ispec_selected_source)) then ! TODO verify if we have to add : .and. myrank == islice_selected_source(i)
      call lagrange_any(xi_source,NGLJ,xiglj,hxis,hpxis)
    else
      call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
    endif
  else
    call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
  endif
  call lagrange_any(gamma_source,NGLLZ,zigll,hgammas,hpgammas)


! calculate source array
  sourcearray(:,:,:) = ZERO

  do m = 1,NGLLZ
    do k = 1,NGLLX

      do iv = 1,NGLLZ
        do ir = 1,NGLLX

          sourcearray(1,k,m) = sourcearray(1,k,m) + &
            real(hxis(ir)*hgammas(iv)*(G11(ir,iv)*hpxis(k)*hgammas(m) + G13(ir,iv)*hxis(k)*hpgammas(m)),kind=CUSTOM_REAL)

          sourcearray(2,k,m) = sourcearray(2,k,m) + &
            real(hxis(ir)*hgammas(iv) *(G31(ir,iv)*hpxis(k)*hgammas(m) + G33(ir,iv)*hxis(k)*hpgammas(m)),kind=CUSTOM_REAL)

        enddo
      enddo

    enddo
  enddo

  end subroutine compute_arrays_source

!
!-----------------------------------------------------------------------------------------
!

  subroutine compute_arrays_adj_source(xi_rec,gamma_rec,irec_local,adj_source_file,adj_sourcearray)

! reads in adjoint source file and computes source arrays for receiver location

  use constants, only: IIN,MAX_STRING_LEN,NGLLX,NGLLZ,NGLJ,CUSTOM_REAL,NDIM

  use specfem_par, only: myrank,NSTEP,P_SV,AXISYM,is_on_the_axis, &
                        xigll,zigll,hxir,hpxir,hgammar,hpgammar,xiglj, &
                        ispec_selected_rec,seismotype,GPU_MODE

  use specfem_par_gpu, only: source_adjointe

  implicit none

  double precision,intent(in) :: xi_rec, gamma_rec
  integer,intent(in) :: irec_local

  character(len=MAX_STRING_LEN),intent(in) :: adj_source_file

  real(kind=CUSTOM_REAL),dimension(NSTEP,NDIM,NGLLX,NGLLZ),intent(out) :: adj_sourcearray

  ! local parameters
  integer :: icomp, itime, i, k
  integer :: ier
  double precision :: junk
  double precision,dimension(NSTEP,3) :: adj_src_s

  character(len=3) :: comp(3)
  character(len=MAX_STRING_LEN) :: filename

  ! initializes temporary array
  ! note: we need to explicitly initialize the full array,
  !       otherwise it can lead to a floating overflow problem (with intel compiler 15.x)
  !       when copying the temporary array into the actual storage array adj_sourcearray
  adj_src_s(:,:) = 0.d0

  if (seismotype == 1 .or. seismotype == 2 .or. seismotype == 3) then
    ! displacement/velocity/acceleration
    comp = (/"BXX","BXY","BXZ"/)

    ! reads corresponding components
    do icomp = 1,3
      ! skips unnecessary components
      if (P_SV) then
        ! P_SV-case: skips BXY component
        if (icomp == 2) cycle
      else
        ! SH-case: skips BXX and BXZ components
        if (icomp == 1 .or. icomp == 3) cycle
      endif

      ! reads in ascii adjoint source files **.adj
      filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'
      open(unit = IIN, file = trim(filename), iostat = ier)
      if (ier /= 0) then
        print *,'Error: could not find adjoint source file ',trim(filename)
        print *,'Please check if file exists...'
        call exit_MPI(myrank,'file '//trim(filename)//' does not exist')
      endif

      do itime = 1, NSTEP
        read(IIN,*) junk, adj_src_s(itime,icomp)
      enddo
      close(IIN)
    enddo

  else if (seismotype == 4) then
    ! pressure
    ! reads in ascii adjoint source files **.PRE.adj
    filename = 'SEM/'//trim(adj_source_file) // '.PRE.adj'
    open(unit = IIN, file = trim(filename), iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

    ! reads only single component
    do itime = 1, NSTEP
      read(IIN,*) junk, adj_src_s(itime,1)
    enddo
    close(IIN)

  else if (seismotype == 6) then
    ! potential
    ! reads in ascii adjoint source files **.POT.adj
    filename = 'SEM/'//trim(adj_source_file) // '.POT.adj'
    open(unit = IIN, file = trim(filename), iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

    ! reads only single component
    do itime = 1, NSTEP
      read(IIN,*) junk, adj_src_s(itime,1)
    enddo
    close(IIN)

  endif

  ! gpu arrays
  if (GPU_MODE) then
    if (P_SV) then
      ! P_SV-case
      source_adjointe(irec_local,:,1) = real(adj_src_s(:,1),kind=CUSTOM_REAL)
      source_adjointe(irec_local,:,2) = real(adj_src_s(:,3),kind=CUSTOM_REAL)
    else
      ! SH-case
      source_adjointe(irec_local,:,1) = real(adj_src_s(:,2),kind=CUSTOM_REAL)
    endif
  endif

  ! location interpolation
  if (AXISYM) then
    if (is_on_the_axis(ispec_selected_rec(irec_local))) then ! TODO verify irec_local...
      call lagrange_any(xi_rec,NGLJ,xiglj,hxir,hpxir)
    else
      call lagrange_any(xi_rec,NGLLX,xigll,hxir,hpxir)
    endif
  else
    call lagrange_any(xi_rec,NGLLX,xigll,hxir,hpxir)
  endif
  call lagrange_any(gamma_rec,NGLLZ,zigll,hgammar,hpgammar)

  ! fills into full adjoint source array (NSTEP,NDIM,NGLLX,NGLLZ)
  ! note: adj_sourcearray is defined as CUSTOM_REAL
  adj_sourcearray(:,:,:,:) = 0._CUSTOM_REAL
  do k = 1, NGLLZ
    do i = 1, NGLLX
      if (P_SV) then
        ! P_SV-case
        adj_sourcearray(:,1,i,k) = real(hxir(i) * hgammar(k) * adj_src_s(:,1),kind=CUSTOM_REAL)
        adj_sourcearray(:,2,i,k) = real(hxir(i) * hgammar(k) * adj_src_s(:,3),kind=CUSTOM_REAL)
      else
        ! SH-case
        adj_sourcearray(:,1,i,k) = real(hxir(i) * hgammar(k) * adj_src_s(:,2),kind=CUSTOM_REAL)
      endif
    enddo
  enddo

  end subroutine compute_arrays_adj_source

!
!-----------------------------------------------------------------------------------------
!

  subroutine compute_arrays_adj_source_SU(seismotype)

! reads in all adjoint sources in SU-file format and computes source arrays for receiver locations

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLZ,NGLJ,NDIM

  use specfem_par, only: AXISYM,xiglj,is_on_the_axis, &
                         myrank, NSTEP, nrec, &
                         xi_receiver, gamma_receiver, islice_selected_rec, &
                         xigll,zigll,hxir,hgammar,hpxir,hpgammar, &
                         adj_sourcearrays, &
                         GPU_MODE, ispec_selected_rec,P_SV

  use specfem_par_gpu, only: source_adjointe

  implicit none

  integer,intent(in) :: seismotype

  ! local parameters
  integer :: i, j, irec, irec_local, ier
  character(len=MAX_STRING_LEN) :: filename

  ! SU
  integer(kind=4) :: r4head(60)
  integer(kind=2) :: header2(2)

  real(kind=4),dimension(:,:),allocatable :: adj_src_s

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: adj_sourcearray

  ! opens adjoint source files in SU format
  if (seismotype == 4 .or. seismotype == 6) then
    ! pressure/potential type
    write(filename, "('./SEM/Up_file_single.su.adj')")
    open(111,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')
  else
    ! displacement/velocity/acceleration type
    if (P_SV) then
      ! P_SV-case
      write(filename, "('./SEM/Ux_file_single.su.adj')")
      open(111,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
      if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

      write(filename, "('./SEM/Uz_file_single.su.adj')")
      open(113,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
      if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')
    else
      ! SH-case
      write(filename, "('./SEM/Uy_file_single.su.adj')")
      open(112,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
      if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')
    endif

  endif

  ! allocates temporary array
  allocate(adj_src_s(NSTEP,NDIM))

  ! temporary array
  allocate(adj_sourcearray(NSTEP,NDIM,NGLLX,NGLLZ))

  irec_local = 0
  do irec = 1, nrec

    ! only process/slice holding receiver
    if (myrank == islice_selected_rec(irec)) then
      irec_local = irec_local + 1

      adj_sourcearray(:,:,:,:) = 0._CUSTOM_REAL
      adj_src_s(:,:) = 0.0

      if (seismotype == 4 .or. seismotype == 6) then
        ! pressure/potential
        ! single component
        read(111,rec=irec,iostat=ier) r4head, adj_src_s(:,1)
        if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' read error')
      else
        ! displacement/velocity/acceleration
        if (P_SV) then
          ! P_SV-case
          read(111,rec=irec,iostat=ier) r4head, adj_src_s(:,1)
          if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' read error')

          read(113,rec=irec,iostat=ier) r4head, adj_src_s(:,2)
          if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' read error')
        else
          ! SH-case
          read(112,rec=irec,iostat=ier) r4head, adj_src_s(:,1)
          if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' read error')
        endif
      endif

      header2 = int(r4head(29), kind=2)
      if (irec == 1) print *, r4head(1),r4head(19),r4head(20),r4head(21),r4head(22),header2(2)

      ! location interpolation
      if (AXISYM) then
        if (is_on_the_axis(ispec_selected_rec(irec))) then ! TODO verify ispec_selected_rec and not ispec_selected_source
          call lagrange_any(xi_receiver(irec),NGLJ,xiglj,hxir,hpxir)
        else
          call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
        endif
      else
        call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
      endif

      call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)

      ! gpu arrays
      if (GPU_MODE) then
        if (P_SV) then
          ! P_SV-case
          source_adjointe(irec_local,:,1) = real(adj_src_s(:,1),kind=CUSTOM_REAL)
          source_adjointe(irec_local,:,2) = real(adj_src_s(:,2),kind=CUSTOM_REAL)
        else
          ! SH-case
          source_adjointe(irec_local,:,1) = real(adj_src_s(:,1),kind=CUSTOM_REAL)
        endif
      endif

      ! adjoint source interpolated on element
      if (.not. GPU_MODE) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            if (P_SV) then
              ! P_SV-case
              adj_sourcearray(:,1,i,j) = real(hxir(i) * hgammar(j) * adj_src_s(:,1),kind=CUSTOM_REAL)
              adj_sourcearray(:,2,i,j) = real(hxir(i) * hgammar(j) * adj_src_s(:,2),kind=CUSTOM_REAL)
            else
              ! SH-case
              adj_sourcearray(:,1,i,j) = real(hxir(i) * hgammar(j) * adj_src_s(:,1),kind=CUSTOM_REAL)
            endif
          enddo
        enddo

        ! stores all local adjoint sources
        adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
      endif

    endif !  if (myrank == islice_selected_rec(irec))
  enddo ! irec

  ! closes files
  if (seismotype == 4) then
    close(111)
  else
    if (P_SV) then
      ! P_SV-case
      close(111)
      close(113)
    else
      ! SH-case
      close(112)
    endif
  endif

  ! frees memory
  deallocate(adj_src_s)
  deallocate(adj_sourcearray)

  end subroutine compute_arrays_adj_source_SU
