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
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  integer :: k,m

  double precision :: hlagrange
  double precision :: dsrc_dx,dsrc_dz
  double precision :: dxis_dx,dgammas_dx
  double precision :: dxis_dz,dgammas_dz

  ! debugging (original version)
  logical, parameter :: DEBUG = .false.
  double precision, dimension(NGLLX,NGLLZ) :: G11,G13,G31,G33
  integer :: ir,iv

! compute Lagrange polynomials at the source location
! the source does not necessarily correspond to a Gauss-Lobatto point
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

  dxis_dx = ZERO
  dxis_dz = ZERO
  dgammas_dx = ZERO
  dgammas_dz = ZERO

  do m = 1,NGLLZ
    do k = 1,NGLLX

        xixd    = xix(k,m,ispec_selected_source)
        xizd    = xiz(k,m,ispec_selected_source)
        gammaxd = gammax(k,m,ispec_selected_source)
        gammazd = gammaz(k,m,ispec_selected_source)

        hlagrange = hxis(k) * hgammas(m)

        dxis_dx = dxis_dx + hlagrange * xixd
        dxis_dz = dxis_dz + hlagrange * xizd
        dgammas_dx = dgammas_dx + hlagrange * gammaxd
        dgammas_dz = dgammas_dz + hlagrange * gammazd

    enddo
  enddo

! calculate source array
  sourcearray(:,:,:) = ZERO

  do m = 1,NGLLZ
    do k = 1,NGLLX

        dsrc_dx = (hpxis(k)*dxis_dx)*hgammas(m) + hxis(k)*(hpgammas(m)*dgammas_dx)
        dsrc_dz = (hpxis(k)*dxis_dz)*hgammas(m) + hxis(k)*(hpgammas(m)*dgammas_dz)

        ! formula: see notes in doc/notes_from_Youshan_Liu_The_point_moment_tensor_source_with_merged_loops_in_2D.pdf
        sourcearray(1,k,m) = sourcearray(1,k,m) + real(Mxx*dsrc_dx + Mxz*dsrc_dz,kind=CUSTOM_REAL)
        sourcearray(2,k,m) = sourcearray(2,k,m) + real(Mxz*dsrc_dx + Mzz*dsrc_dz,kind=CUSTOM_REAL)

    enddo
  enddo

  ! debugging
  if (DEBUG) then
    ! compares against original double-loop version
    !
    ! calculate G_ij for general source location
    ! the source does not necessarily correspond to a Gauss-Lobatto point
    do m = 1,NGLLZ
      do k = 1,NGLLX
        xixd    = xix(k,m,ispec_selected_source)
        xizd    = xiz(k,m,ispec_selected_source)
        gammaxd = gammax(k,m,ispec_selected_source)
        gammazd = gammaz(k,m,ispec_selected_source)
        G11(k,m) = Mxx * xixd    + Mxz * xizd
        G13(k,m) = Mxx * gammaxd + Mxz * gammazd
        G31(k,m) = Mxz * xixd    + Mzz * xizd
        G33(k,m) = Mxz * gammaxd + Mzz * gammazd
      enddo
    enddo
    sourcearray(:,:,:) = ZERO
    do m = 1,NGLLZ
      do k = 1,NGLLX
        do iv = 1,NGLLZ
          do ir = 1,NGLLX
            sourcearray(1,k,m) = sourcearray(1,k,m) + &
              real(hxis(ir)*hgammas(iv) * (G11(ir,iv)*hpxis(k)*hgammas(m) + G13(ir,iv)*hxis(k)*hpgammas(m)),kind=CUSTOM_REAL)
            sourcearray(2,k,m) = sourcearray(2,k,m) + &
              real(hxis(ir)*hgammas(iv) * (G31(ir,iv)*hpxis(k)*hgammas(m) + G33(ir,iv)*hxis(k)*hpgammas(m)),kind=CUSTOM_REAL)
          enddo
        enddo
      enddo
    enddo
  endif

  end subroutine compute_arrays_source

!
!-----------------------------------------------------------------------------------------
!

  subroutine read_adj_source(irec_local,seismotype_adj,adj_source_file)

! reads in adjoint source file

  use constants, only: IIN, MAX_STRING_LEN, CUSTOM_REAL, NDIM

  use specfem_par, only: myrank, NSTEP, DT, P_SV, source_adjoint, &
                         f0_source, vx_source, vz_source, SOURCE_IS_MOVING
  use specfem_par_movie, only: vpImax

  implicit none

  integer,intent(in) :: irec_local
  integer,intent(in) :: seismotype_adj
  character(len=MAX_STRING_LEN),intent(in) :: adj_source_file

  ! local parameters
  integer :: icomp, itime
  integer :: ier, irek, norder
  double precision :: junk, vx_source_max, vz_source_max, v_source_max
  double precision :: f1, f2
  real(kind=CUSTOM_REAL), dimension(NSTEP) :: temp
  real(kind=CUSTOM_REAL), dimension(NSTEP,3) :: adj_src_s

  character(len=3) :: comp(3)
  character(len=MAX_STRING_LEN) :: filename

  ! initializes temporary array
  ! note: we need to explicitly initialize the full array,
  !       otherwise it can lead to a floating overflow problem (with intel compiler 15.x)
  !       when copying the temporary array into the actual storage array adj_sourcearray
  adj_src_s(:,:) = 0.d0
  temp(:) = 0.d0

  select case(seismotype_adj)
  case (1,2,3)
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

  case (4)
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

  case (6)
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

  case default
    ! unknown, not implement yet
    print *,'Error: unknown adjoint source type ',seismotype_adj
    call exit_MPI(myrank,'Invalid adjoint source type')

  end select

  if (SOURCE_IS_MOVING) then
    ! In this case high frequency noise is produced. We have to filter it out before sending it at the receivers
    irek = 1
    f1 = 0.0d0 ! High pass filter.
    f2 = maxval(f0_source) * 5.0 ! Cutoff frequency
    ! Takes also Doppler shift into account:
    vx_source_max = maxval(vx_source)
    vz_source_max = maxval(vz_source)
    v_source_max = max(vx_source_max, vz_source_max)
    if ((vpImax > 0) .and. (vpImax < 15000.0)) then  ! Check that it has been calculated (value makes sense)
      f2 = f2 / (1.0d0 - v_source_max/vpImax)
    else
      ! It should have been calculated in check_grid.F90
      call exit_MPI(myrank, 'It looks like vpImax has not been calculated here yet...')
    endif
    norder = 4
    ! Filter adjoint signals :
    call exit_MPI(myrank, 'This part has to be tested. Make sure the following filtering is ' // &
                          'going well in single and double precicison. This is just a 5 min check, ' // &
                          'then you can remove this stop')
    if (P_SV) then
      call bwfilt(adj_src_s(:,1), temp, DT, NSTEP, irek, norder, f1, f2)
      adj_src_s(:,1) = temp
      call bwfilt(adj_src_s(:,3), temp, DT, NSTEP, irek, norder, f1, f2)
      adj_src_s(:,3) = temp
    else ! SH-case
      call bwfilt(adj_src_s(:,2), temp, DT, NSTEP, irek, norder, f1, f2)
      adj_src_s(:,2) = temp
    endif
  endif

  if (P_SV) then
    ! P_SV-case
    source_adjoint(irec_local,:,1) = real(adj_src_s(:,1),kind=CUSTOM_REAL)
    source_adjoint(irec_local,:,2) = real(adj_src_s(:,3),kind=CUSTOM_REAL)
  else
    ! SH-case
    source_adjoint(irec_local,:,1) = real(adj_src_s(:,2),kind=CUSTOM_REAL)
  endif

  end subroutine read_adj_source

!
!-----------------------------------------------------------------------------------------
!

  subroutine read_adj_source_SU(seismotype_adj)

! reads in all adjoint sources in SU-file format

  use constants, only: CUSTOM_REAL, MAX_STRING_LEN, NDIM

  use specfem_par, only: myrank, NSTEP, nrec, DT, f0_source, &
                         islice_selected_rec, P_SV, source_adjoint, &
                         SOURCE_IS_MOVING

  implicit none

  integer,intent(in) :: seismotype_adj

  ! local parameters
  integer :: irec, irec_local, ier
  integer :: irek, norder
  double precision :: f1, f2
  real(kind=CUSTOM_REAL), dimension(NSTEP) :: temp
  character(len=MAX_STRING_LEN) :: filename

  ! SU
  integer(kind=4) :: r4head(60)
  integer(kind=2) :: header2(2)

  real(kind=4),dimension(:,:),allocatable :: adj_src_s

  ! opens adjoint source files in SU format
  if (seismotype_adj == 4 .or. seismotype_adj == 6) then
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

  irec_local = 0
  do irec = 1, nrec

    ! only process/slice holding receiver
    if (myrank == islice_selected_rec(irec)) then
      irec_local = irec_local + 1

      adj_src_s(:,:) = 0.d0
      temp(:) = 0.d0

      if (seismotype_adj == 4 .or. seismotype_adj == 6) then
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

      if (SOURCE_IS_MOVING) then
        ! In this case high frequency noise is produced. We have to filter it out before sending it at the receivers
        irek = 1
        f1 = 0.d0
        f2 = maxval(f0_source) * 5.0
        norder = 4
        ! Filter adjoint signals :
        if (P_SV) then
          call bwfilt(real(adj_src_s(:,1),kind=CUSTOM_REAL), temp, DT, NSTEP, irek, norder, f1, f2)
          adj_src_s(:,1) = real(temp(:),kind=4)
          call bwfilt(real(adj_src_s(:,2),kind=CUSTOM_REAL), temp, DT, NSTEP, irek, norder, f1, f2)
          adj_src_s(:,2) = real(temp(:),kind=4)
        else ! SH-case
          call bwfilt(real(adj_src_s(:,1),kind=CUSTOM_REAL), temp, DT, NSTEP, irek, norder, f1, f2)
          adj_src_s(:,1) = real(temp(:),kind=4)
        endif
      endif

      header2 = int(r4head(29), kind=2)

      if (P_SV) then
        ! P_SV-case
        source_adjoint(irec_local,:,1) = real(adj_src_s(:,1),kind=CUSTOM_REAL)
        source_adjoint(irec_local,:,2) = real(adj_src_s(:,2),kind=CUSTOM_REAL)
      else
        ! SH-case
        source_adjoint(irec_local,:,1) = real(adj_src_s(:,1),kind=CUSTOM_REAL)
      endif

    endif !  if (myrank == islice_selected_rec(irec))
  enddo ! irec

  ! closes files
  if (seismotype_adj == 4) then
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

  end subroutine read_adj_source_SU

!
!-----------------------------------------------------------------------------------------
!

  subroutine bwfilt(x, y, dt, n, irek, norder, f1, f2)

    ! recursive filtering of data with butterworth filter
    ! x: input array
    ! y: output array
    ! dt: time increment
    ! n: number of data points

    ! irek=0: forward filtering only
    ! irek=1: forward and backward filtering

    ! norder: order of butterworth filter
    ! norder=0: only filtering, no determination of coefficients
    ! norder < 0: no starplots of transfer function and impulse response

    ! f1: low cutoff frequency (Hz)
    ! f1=0: low pass filter

    ! f2: high cutoff frequency (Hz)
    ! f2>0.5/dt: high pass filter

    use constants, only: CUSTOM_REAL

    implicit none

    real(kind=CUSTOM_REAL), dimension(1)::x,y
    real(kind=CUSTOM_REAL), dimension (10) ::  a, b1, b2
    double precision :: dt,f1,f2
    integer :: iunit, npoles,norder,irek,n,lx
    !real(kind(0d0)) :: x(n),y(n)

     iunit = 3

     if (norder /= 0) then
        npoles=iabs(norder)
        !determination of filter coefficients
        call bpcoeff(f1,f2,npoles, dt, a,b1, b2)
        if (norder >= 0) then
           !plot of transfer function and impuulse response
           lx = 100
           !filtering
        endif
     endif


     if (n /= 0) then
        call rekurs(x,y,n,a,b1,b2,npoles,irek)
     endif
     return
   end subroutine bwfilt

  !---------------------------------------------------------------


  subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)
    ! performs recursive filtering of data in array x of length ndat
    ! filtered output in y
    ! a, b1, b2 are the filtercoefficients previously determined in bwcoef
    ! npoles is the number of poles
    ! iflag=0: forward filtering only
    ! iflag /= 0: forward and backward filtering

    use constants, only: CUSTOM_REAL

    implicit none

    real(kind=CUSTOM_REAL), dimension(10) :: z,z1,z2 ,a,b1,b2
    real(kind=CUSTOM_REAL) ::  x1,x2
    integer :: ndat, npoles, iflag, n,i
    real(kind=CUSTOM_REAL) :: x(ndat), y(ndat)

    !forward

    x1 = 0.d0
    x2 = 0.d0

    do i = 1, npoles
       z1(i) = 0.d0
       z2(i) = 0.d0
    enddo

    do n = 1, ndat
       z(1) = a(1)*(x(n)-x2) -b1(1)*z1(1) -b2(1)*z2(1)
       do i = 2, npoles
          z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
       enddo
       x2=x1
       x1=x(n)
       do i = 1, npoles
          z2(i) =z1(i)
          z1(i) =z(i)
       enddo
       y(n) = z(npoles)
    enddo

    if (iflag == 0) then
       return
    endif

    !backward

    x1 =0.d0
    x2 =0.d0

    do i = 1, npoles
       z1(i) = 0.d0
       z2(i) = 0.d0
    enddo

    do n = ndat, 1, -1
       z(1) = a(1)*(y(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
       do i =2, npoles
          z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
       enddo
       x2=x1
       x1=y(n)
       do i = 1,npoles
          z2(i)=z1(i)
          z1(i)=z(i)
       enddo
       y(n) = z(npoles)
    enddo
    return
  end subroutine rekurs


  !---------------------------------------------------------------


  subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)
    !determines filtercoefficients for recursive bandpassfilter

    use constants, only: CUSTOM_REAL

    real(kind=CUSTOM_REAL),dimension(10) :: a,b1,b2
    complex(kind=CUSTOM_REAL) :: s(20), t1,t2,p
    real(kind=CUSTOM_REAL), parameter :: pi = 3.141592653589793d0
    double precision :: f1,f2,dt,d2,w0,w1,w2,ssum, sprod,fact1,fact2,fact3
    integer :: i,npol2,n,npoles

    if (npoles > 10) then
       stop ' npoles greater than 10: STOP '
    endif

    d2= 2.d0/dt
    w1=d2*tan(2.d0*pi*f1/d2)
    w2=d2*tan(2.d0*pi*f2/d2)
    w0=0.5*(w2-w1)

    i=1
    npol2=npoles/2+1
    do n =1,npoles
       p = cexp(cmplx(0.d0,dble(2*n-1+npoles)*pi/dble(2*npoles)))
       t1 = p*cmplx(w0,0.d0)
       t2 = sqrt(t1*t1-cmplx(w1*w2,0.d0))
       s(i)=t1+t2
       s(i+1)=t1-t2
       i=i+2
    enddo

    do n=1,npoles
       ssum=2*real(s(n))
       sprod=dble(s(n)*conjg(s(n)))
       fact1=d2*d2-d2*ssum+sprod
       fact2=2.d0*(sprod-d2*d2)
       fact3=d2*d2+d2*ssum+sprod
       a(n)=2.d0*d2*w0/fact1
       b1(n)=fact2/fact1
       b2(n)=fact3/fact1
    enddo
    return
  end subroutine bpcoeff
