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

  double precision FUNCTION HGLL (I,Z,ZGLL,NZ)
!-------------------------------------------------------------
!
!     Compute the value of the Lagrangian interpolant L through
!     the NZ Gauss-Lobatto Legendre points ZGLL at the point Z.
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
  IF (dABS(DZ)  <  EPS) THEN
   HGLL = 1.d0
   RETURN
  ENDIF
  N = NZ - 1
  ALFAN = dble(N)*(dble(N)+1.d0)
  HGLL = - (1.d0-Z*Z)*PNDLEG(Z,N)/ (ALFAN*PNLEG(ZGLL(I),N)*(Z-ZGLL(I)))
  RETURN
  end function hgll
