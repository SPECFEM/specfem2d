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

! subroutines to sort indexing arrays based on geometrical coordinates instead of based on topology (because that is much faster)

! subroutines to sort MPI buffers to assemble between chunks

  subroutine sort_array_coordinates(npointot,x,z,ibool,iglob,locval,ifseg, &
                                    nglob,ind,ninseg,iwork,work)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points
!
! returns: sorted indexing array (ibool),  reordering array (iglob) & number of global points (nglob)

  use constants, only: NGLLX,NGLLZ,NDIM,SMALLVALTOL

  implicit none

  integer, intent(in) :: npointot
  integer, intent(out) :: nglob

  integer, intent(inout) :: ibool(npointot)

  integer iglob(npointot),locval(npointot)
  integer ind(npointot),ninseg(npointot)
  logical ifseg(npointot)
  double precision, intent(in) :: x(npointot),z(npointot)
  integer iwork(npointot)
  double precision work(npointot)

  ! local parameters
  integer :: ipoin, i, j
  integer :: nseg, ioff, iseg, ig
  ! define a tolerance, normalized radius is 1., so let's use a small value
  double precision,parameter :: xtol = SMALLVALTOL

  ! establish initial pointers
  do ipoin = 1,npointot
    locval(ipoin) = ipoin
  enddo

  ifseg(:) = .false.

  nseg = 1
  ifseg(1) = .true.
  ninseg(1) = npointot

  do j = 1,NDIM

    ! sort within each segment
    ioff = 1
    do iseg = 1,nseg
      if (j == 1) then
        call rank_buffers(x(ioff),ind,ninseg(iseg))
      else if (j == 2) then
        call rank_buffers(z(ioff),ind,ninseg(iseg))
      endif

      call swap_all_buffers(ibool(ioff),locval(ioff), &
                  x(ioff),z(ioff),iwork,work,ind,ninseg(iseg))

      ioff = ioff + ninseg(iseg)
    enddo

    ! check for jumps in current coordinate
    ! define a tolerance, normalized radius is 1., so let's use a small value
    if (j == 1) then
      do i = 2,npointot
        if (dabs(x(i) - x(i-1)) > xtol) ifseg(i) = .true.
      enddo
    else if (j == 2) then
      do i = 2,npointot
        if (dabs(z(i) - z(i-1)) > xtol) ifseg(i) = .true.
      enddo
    endif

    ! count up number of different segments
    nseg = 0
    do i = 1,npointot
      if (ifseg(i)) then
        nseg = nseg + 1
        ninseg(nseg) = 1
      else
        ninseg(nseg) = ninseg(nseg) + 1
      endif
    enddo

  enddo

  ! assign global node numbers (now sorted lexicographically)
  ig = 0
  do i = 1,npointot
    ! eliminate the multiples by using a single (new) point number for all the points that have the same X Y Z after sorting
    if (ifseg(i)) ig = ig + 1
    iglob(locval(i)) = ig
  enddo

  nglob = ig

  end subroutine sort_array_coordinates

!
!--------------------
!


! sorting routine left here for inlining

! -------------------- library for sorting routine ------------------

! sorting routines put here in same file to allow for inlining

  subroutine rank_buffers(A,IND,N)
!
! Use Heap Sort (Numerical Recipes)
!
  implicit none

  integer n
  double precision A(n)
  integer IND(n)

  integer i,j,l,ir,indx
  double precision q

  do j = 1,n
    IND(j)=j
  enddo

  if (n == 1) return

  L=n/2+1
  ir=n
  100 continue
   if (l > 1) then
      l=l-1
      indx=IND(l)
      q=A(indx)
   ELSE
      indx=IND(ir)
      q=A(indx)
      IND(ir)=IND(1)
      ir=ir-1
      if (ir == 1) then
         IND(1)=indx
         return
      endif
   endif
   i=l
   j=l+l
  200    continue
   if (j <= ir) then
      if (j < ir) then
         if (A(IND(j)) < A(IND(j+1))) j=j+1
      endif
      if (q < A(IND(j))) then
         IND(i)=IND(j)
         i=j
         j=j+j
      ELSE
         j=ir+1
      endif
   goto 200
   endif
   IND(i)=indx
  goto 100
  end subroutine rank_buffers

! -------------------------------------------------------------------

  subroutine swap_all_buffers(IA,IB,A,B,IW,W,ind,n)
!
! swap arrays IA, IB, A and B according to addressing in array IND
!
  implicit none

  integer n

  integer IND(n)
  integer IA(n),IB(n),IW(n)
  double precision A(n),B(n),W(n)

  integer i

  do i = 1,n
    W(i)=A(i)
    IW(i)=IA(i)
  enddo

  do i = 1,n
    A(i)=W(ind(i))
    IA(i)=IW(ind(i))
  enddo

  do i = 1,n
    W(i)=B(i)
    IW(i)=IB(i)
  enddo

  do i = 1,n
    B(i)=W(ind(i))
    IB(i)=IW(ind(i))
  enddo

  end subroutine swap_all_buffers

