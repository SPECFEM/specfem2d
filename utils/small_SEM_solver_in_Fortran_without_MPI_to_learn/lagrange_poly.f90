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

  double precision function hgll(I,Z,ZGLL,NZ)

!-------------------------------------------------------------
!
!  Compute the value of the Lagrangian interpolant L through
!  the NZ Gauss-Lobatto Legendre points ZGLL at point Z
!
!-------------------------------------------------------------

  implicit none

  integer i,nz
  double precision z
  double precision ZGLL(0:nz-1)

  integer n
  double precision EPS,DZ,ALFAN
  double precision, external :: PNLEG,PNDLEG

  EPS = 1.d-5
  DZ = Z - ZGLL(I)
  if (abs(DZ) < EPS) then
   HGLL = 1.d0
   return
  endif
  N = NZ - 1
  ALFAN = dble(N)*(dble(N)+1.d0)
  HGLL = - (1.d0-Z*Z)*PNDLEG(Z,N)/ (ALFAN*PNLEG(ZGLL(I),N)*(Z-ZGLL(I)))

  end function hgll

!
!=====================================================================
!

  double precision function hglj(I,Z,ZGLJ,NZ)

!-------------------------------------------------------------
!
!  Compute the value of the Lagrangian interpolant L through
!  the NZ Gauss-Lobatto Jacobi points ZGLJ at point Z
!
!-------------------------------------------------------------

  implicit none

  integer i,nz
  double precision z
  double precision ZGLJ(0:nz-1)

  integer n
  double precision EPS,DZ,ALFAN
  double precision, external :: PNGLJ,PNDGLJ

  EPS = 1.d-5
  DZ = Z - ZGLJ(I)
  if (abs(DZ) < EPS) then
   HGLJ = 1.d0
   return
  endif
  N = NZ - 1
  ALFAN = dble(N)*(dble(N)+2.d0)
  HGLJ = - (1.d0-Z*Z)*PNDGLJ(Z,N)/ (ALFAN*PNGLJ(ZGLJ(I),N)*(Z-ZGLJ(I)))

  end function hglj

!
!=====================================================================
!

  subroutine lagrange_any(xi,NGLL,xigll,h,hprime)

! subroutine to compute the Lagrange interpolants based upon the GLL points
! and their first derivatives at any point xi in [-1,1]

  implicit none

  integer NGLL
  double precision xi,xigll(NGLL),h(NGLL),hprime(NGLL)

  integer dgr,i,j
  double precision prod1,prod2

  do dgr=1,NGLL

    prod1 = 1.0d0
    prod2 = 1.0d0
    do i=1,NGLL
      if (i /= dgr) then
        prod1 = prod1*(xi-xigll(i))
        prod2 = prod2*(xigll(dgr)-xigll(i))
      endif
    enddo
    h(dgr)=prod1/prod2

    hprime(dgr)=0.0d0
    do i=1,NGLL
      if (i /= dgr) then
        prod1=1.0d0
        do j=1,NGLL
          if (j /= dgr .and. j /= i) prod1 = prod1*(xi-xigll(j))
        enddo
        hprime(dgr) = hprime(dgr)+prod1
      endif
    enddo
    hprime(dgr) = hprime(dgr)/prod2

  enddo

  end subroutine lagrange_any

!
!=====================================================================
!

! subroutine to compute the derivative of the Lagrange interpolants
! at the GLL points at any given GLL point

  double precision function lagrange_deriv_GLL(I,j,ZGLL,NZ)

!------------------------------------------------------------------------
!
!     Compute the value of the derivative of the I-th
!     Lagrange interpolant through the
!     NZ Gauss-Lobatto Legendre points ZGLL at point ZGLL(j)
!
!------------------------------------------------------------------------

  implicit none

  integer i,j,nz
  double precision zgll(0:nz-1)

  integer degpoly

  double precision, external :: pnleg,pndleg

  degpoly = nz - 1
  if (i == 0 .and. j == 0) then
    lagrange_deriv_GLL = - dble(degpoly)*(dble(degpoly)+1.d0) / 4.d0
  else if (i == degpoly .and. j == degpoly) then
    lagrange_deriv_GLL = dble(degpoly)*(dble(degpoly)+1.d0) / 4.d0
  else if (i == j) then
    lagrange_deriv_GLL = 0.d0
  else
    lagrange_deriv_GLL = pnleg(zgll(j),degpoly) / &
      (pnleg(zgll(i),degpoly)*(zgll(j)-zgll(i))) &
      + (1.d0-zgll(j)*zgll(j))*pndleg(zgll(j),degpoly) / (dble(degpoly)* &
      (dble(degpoly)+1.d0)*pnleg(zgll(i),degpoly)*(zgll(j)-zgll(i))*(zgll(j)-zgll(i)))
  endif

  end function lagrange_deriv_GLL

!
!=======================================================================
!

! subroutine to compute the derivative of the interpolants of the GLJ
! quadrature at the GLJ points at any given GLJ point

  double precision function poly_deriv_GLJ(I,j,ZGLJ,NZ)

!------------------------------------------------------------------------
!
!     Compute the value of the derivative of the I-th
!     polynomial interpolant of the GLJ quadrature through the
!     NZ Gauss-Lobatto-Jacobi (0,1) points ZGLJ at point ZGLJ(j)
!
!------------------------------------------------------------------------

  implicit none

  integer i,j,nz
  double precision zglj(0:nz-1)

  integer degpoly

  double precision, external :: pnglj

  degpoly = nz - 1

  if (i == 0 .and. j == 0) then ! Case 1
    poly_deriv_GLJ = -dble(degpoly)*(dble(degpoly)+2.d0)/6.d0
  else if (i == 0 .and. 0 < j .and. j < degpoly) then ! Case 2
    poly_deriv_GLJ = 2.d0*(-1)**degpoly*pnglj(zglj(j),degpoly)/((1.d0+zglj(j))*(dble(degpoly)+1.d0))
  else if (i == 0 .and. j == degpoly) then ! Case 3
    poly_deriv_GLJ = (-1)**degpoly/(dble(degpoly)+1.d0)
  else if (0 < i .and. i < degpoly .and. j == 0) then ! Case 4
    poly_deriv_GLJ = (-1)**(degpoly+1)*(dble(degpoly)+1.d0)/(2.d0*pnglj(zglj(i),degpoly)*(1.d0+zglj(i)))
  else if (0 < i .and. i < degpoly .and. 0 < j .and. j < degpoly .and. i /= j) then ! Case 5
    poly_deriv_GLJ = 1.d0/(zglj(j)-zglj(i))*pnglj(zglj(j),degpoly)/pnglj(zglj(i),degpoly)
  else if (0 < i .and. i < degpoly .and. i == j) then  ! Case 6
    poly_deriv_GLJ = -1.d0/(2.d0*(1.d0+zglj(i)))
  else if (0 < i .and. i < degpoly .and. j == degpoly) then ! Case 7
    poly_deriv_GLJ = 1.d0/(pnglj(zglj(i),degpoly)*(1.d0-zglj(i)))
  else if (i == degpoly .and. j == 0) then ! Case 8
    poly_deriv_GLJ = (-1)**(degpoly+1)*(dble(degpoly)+1.d0)/4.d0
  else if (i == degpoly .and. 0 < j .and. j < degpoly) then ! Case 9
    poly_deriv_GLJ = -1.d0/(1.d0-zglj(j))*pnglj(zglj(j),degpoly)
  else if (i == degpoly .and. j == degpoly) then ! Case 10
    poly_deriv_GLJ = (dble(degpoly)*(dble(degpoly)+2.d0)-1.d0)/4.d0
  else
    stop 'Problem in poly_deriv_GLJ: in a perfect world this would NEVER appear'
  endif

  end function poly_deriv_GLJ

