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

  double precision FUNCTION PNDLEG (Z,N)
!-------------------------------------------------------------
!
!     Compute the derivative of the Nth order Legendre polynomial at Z.
!     Based on the recursion formula for the Legendre polynomials.
!
!-------------------------------------------------------------
  implicit none

  double precision z
  integer n

  double precision P1,P2,P1D,P2D,P3D,FK,P3
  integer k

  P1   = 1.d0
  P2   = Z
  P1D  = 0.d0
  P2D  = 1.d0
  P3D  = 1.d0
  DO 10 K = 1, N-1
   FK  = dble(K)
   P3  = ((2.d0*FK+1.d0)*Z*P2 - FK*P1)/(FK+1.d0)
   P3D = ((2.d0*FK+1.d0)*P2 + (2.d0*FK+1.d0)*Z*P2D - FK*P1D) &
                          /(FK+1.d0)
   P1  = P2
   P2  = P3
   P1D = P2D
   P2D = P3D
 10   CONTINUE
  PNDLEG = P3D
  RETURN
  end function pndleg
