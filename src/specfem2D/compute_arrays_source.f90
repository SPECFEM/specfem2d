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

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM,ZERO

  use specfem_par, only : AXISYM,is_on_the_axis,xiglj

  implicit none

  integer ispec_selected_source
  integer nspec

  double precision xi_source,gamma_source
  double precision Mxx,Mzz,Mxz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  double precision xixd,xizd,gammaxd,gammazd

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll

! source arrays
  double precision, dimension(NGLLX,NGLLZ) :: G11,G13,G31,G33
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  integer k,m
  integer ir,iv

! calculate G_ij for general source location
! the source does not necessarily correspond to a Gauss-Lobatto point
  do m = 1,NGLLZ
      do k = 1,NGLLX

        xixd    = xix(k,m,ispec_selected_source)
        xizd    = xiz(k,m,ispec_selected_source)
        gammaxd = gammax(k,m,ispec_selected_source)
        gammazd = gammaz(k,m,ispec_selected_source)

        G11(k,m) = Mxx*xixd+Mxz*xizd
        G13(k,m) = Mxx*gammaxd+Mxz*gammazd
        G31(k,m) = Mxz*xixd+Mzz*xizd
        G33(k,m) = Mxz*gammaxd+Mzz*gammazd

!!!!        G21(k,m) = Mxy*xixd+Myz*xizd
!!!!        G23(k,m) = Mxy*gammaxd+Myz*gammazd

      enddo
  enddo

! compute Lagrange polynomials at the source location

  if (AXISYM) then
    if (is_on_the_axis(ispec_selected_source)) then ! TODO verify if we have to add : .and. is_proc_source(i)
      call lagrange_any(xi_source,NGLJ,xiglj,hxis,hpxis)
    else
      call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
    endif
  else
    call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
  endif

  call lagrange_any(gamma_source,NGLLZ,zigll,hgammas,hpgammas)

! calculate source array
  do m = 1,NGLLZ
    do k = 1,NGLLX

      sourcearray(:,k,m) = ZERO

      do iv = 1,NGLLZ
        do ir = 1,NGLLX

          sourcearray(1,k,m) = sourcearray(1,k,m) + hxis(ir)*hgammas(iv) &
                                 *(G11(ir,iv)*hpxis(k)*hgammas(m) &
                                 +G13(ir,iv)*hxis(k)*hpgammas(m))

!        sourcearray(2,k,m) = sourcearray(2,k,m) + hxis(ir)*hgammas(iv) &
!                               *(G21(ir,iv)*hpxis(k)*hgammas(m) &
!                               +G23(ir,iv)*hxis(k)*hpgammas(m))

          sourcearray(2,k,m) = sourcearray(2,k,m) + hxis(ir)*hgammas(iv) &
                                 *(G31(ir,iv)*hpxis(k)*hgammas(m) &
                                 +G33(ir,iv)*hxis(k)*hpgammas(m))

        enddo
      enddo

    enddo
  enddo

  end subroutine compute_arrays_source

! ------------------------------------------------------------------------------------------------------


  subroutine compute_arrays_adj_source(xi_rec,gamma_rec,irec_local,adj_source_file,adj_sourcearray)

  use constants,only: IIN,MAX_STRING_LEN,NGLLX,NGLLZ,NGLJ,CUSTOM_REAL,NDIM

  use specfem_par,only: myrank,NSTEP,P_SV, &
                        AXISYM,is_on_the_axis, &
                        xigll,zigll,hxir,hpxir,hgammar,hpgammar,xiglj, &
                        ispec_selected_rec,seismotype

  use specfem_par_gpu,only: source_adjointe

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

    comp = (/"BXX","BXY","BXZ"/)

    ! reads all components
    do icomp = 1,3

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

    filename = 'SEM/'//trim(adj_source_file) // '.PRE.adj'
    open(unit = IIN, file = trim(filename), iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

    ! reads only single component
    do itime = 1, NSTEP
      read(IIN,*) junk, adj_src_s(itime,1)
    enddo
    close(IIN)

  else if (seismotype == 6) then

    filename = 'SEM/'//trim(adj_source_file) // '.POT.adj'
    open(unit = IIN, file = trim(filename), iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

    ! reads only single component
    do itime = 1, NSTEP
      read(IIN,*) junk, adj_src_s(itime,1)
    enddo
    close(IIN)

  endif

  if (P_SV) then
    ! P_SV-case
    source_adjointe(irec_local,:,1) = adj_src_s(:,1)
    source_adjointe(irec_local,:,2) = adj_src_s(:,3)
  else
    ! SH-case
    source_adjointe(irec_local,:,1) = adj_src_s(:,2)
  endif


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
