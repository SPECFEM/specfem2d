!=====================================================================
!
!                 S p e c f e m  V e r s i o n  4 . 2
!                 -----------------------------------
!
!                         Dimitri Komatitsch
!    Department of Earth and Planetary Sciences - Harvard University
!                         Jean-Pierre Vilotte
!                 Departement de Sismologie - IPGP - Paris
!                           (c) June 1998
!
!=====================================================================

  double precision FUNCTION HDGLL (I,j,ZGLL,NZ)
!-------------------------------------------------------------
!
!     Compute the value of the derivative of the I-th
!     Lagrangian interpolant through the
!     NZ Gauss-Lobatto Legendre points ZGLL at point ZGLL(j).
!
!-------------------------------------------------------------

  implicit none

  integer i,j,nz
  double precision zgll(0:nz-1)

  integer idegpoly
  double precision rlegendre1,rlegendre2,rlegendre3

  double precision, external :: pnleg,pndleg

    idegpoly = nz - 1
  if ((i == 0).and.(j == 0)) then
          hdgll = - dble(idegpoly)*(dble(idegpoly)+1.d0)/4.d0
  else if ((i == idegpoly).and.(j == idegpoly)) then
          hdgll = dble(idegpoly)*(dble(idegpoly)+1.d0)/4.d0
  else if (i == j) then
          hdgll = 0.d0
  else
         rlegendre1 = pnleg(zgll(j),idegpoly)
         rlegendre2 = pndleg(zgll(j),idegpoly)
         rlegendre3 = pnleg(zgll(i),idegpoly)
  hdgll = rlegendre1 / (rlegendre3*(zgll(j)-zgll(i))) &
    + (1.d0-zgll(j)*zgll(j))*rlegendre2/(dble(idegpoly)* &
    (dble(idegpoly)+1.d0)*rlegendre3* &
    (zgll(j)-zgll(i))*(zgll(j)-zgll(i)))
  endif

  return
  end FUNCTION hdgll
