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

  subroutine zwgljd (z,w,np,alpha,beta)
!
!=======================================================================
!
!     Z w g l j d : Generate np Gauss-Lobatto-Jacobi points and the
!     -----------   weights associated with Jacobi polynomials of degree
!                   n = np-1.
!
!=======================================================================
!
!     Note : alpha and beta coefficients must be greater than -1.
!     ----
!            Legendre polynomials are special case of Jacobi polynomials
!            just by setting alpha and beta to 0.
!
!=======================================================================
!
  implicit none

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0

  integer np
  double precision alpha,beta
  double precision z(np), w(np)

  integer n,nm1,i
  double precision p,pd,pm1,pdm1,pm2,pdm2
  double precision alpg,betg
  double precision, external :: endw1,endw2
!
!-----------------------------------------------------------------------
!
  p = zero
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  n   = np-1
  nm1 = n-1
  pd  = zero

  if (np <= 1) stop 'Minimum number of Gauss-Lobatto points is 2'

  if ((alpha <= -one).or.(beta <= -one)) &
    stop 'Alpha and Beta must be greater than -1'

  if (nm1 > 0) then
   alpg  = alpha+one
   betg  = beta+one
   call zwgjd (z(2),w(2),nm1,alpg,betg)
  endif
  z(1)  = - one
  z(np) =  one
  do  110 i=2,np-1
   w(i) = w(i)/(one-z(i)**2)
  110 continue
  call jacobf (p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(1))
  w(1)  = endw1 (n,alpha,beta)/(two*pd)
  call jacobf (p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(np))
  w(np) = endw2 (n,alpha,beta)/(two*pd)

  return
  end subroutine zwgljd
