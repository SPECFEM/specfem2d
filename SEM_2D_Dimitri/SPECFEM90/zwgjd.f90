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

  subroutine zwgjd(z,w,np,alpha,beta)
!
!=======================================================================
!
!     Z w g j d : Generate np Gauss-Jacobi points and weights
!                 associated with Jacobi polynomial of degree n = np-1
!
!=======================================================================
!
!     Note : Coefficients alpha and beta must be greater than -1.
!     ----
!
!=======================================================================
!
  implicit none

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0

  integer np
  double precision z(np),w(np)
  double precision alpha,beta

  integer n,np1,np2,i
  double precision p,pd,pm1,pdm1,pm2,pdm2
  double precision apb,dnp1,dnp2,fac1,fac2,fac3,fnorm,rcoef
  double precision, external :: gammaf,pnormj
!
!-----------------------------------------------------------------------
!
  pd = zero
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  n    = np-1
  apb  = alpha+beta
  p    = zero
  pdm1 = zero

  if (np <= 0) stop 'Minimum number of Gauss points is 1'

  if ((alpha <= -one).or.(beta <= -one)) &
    stop 'Alpha and Beta must be greater than -1'

  if (np == 1) then
   z(1) = (beta-alpha)/(apb+two)
   w(1) = gammaf(alpha+one)*gammaf(beta+one)/gammaf(apb+two) * two**(apb+one)
   return
  endif

  call jacg (z,np,alpha,beta)

  np1   = n+1
  np2   = n+2
  dnp1  = dble(np1)
  dnp2  = dble(np2)
  fac1  = dnp1+alpha+beta+one
  fac2  = fac1+dnp1
  fac3  = fac2+one
  fnorm = pnormj(np1,alpha,beta)
  rcoef = (fnorm*fac2*fac3)/(two*fac1*dnp2)
  do i=1,np
    call jacobf (p,pd,pm1,pdm1,pm2,pdm2,np2,alpha,beta,z(i))
    w(i) = -rcoef/(p*pdm1)
  enddo

  return
  end subroutine zwgjd
