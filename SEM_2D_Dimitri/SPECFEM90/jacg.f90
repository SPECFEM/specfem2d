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

  subroutine jacg (xjac,np,alpha,beta)
!
!=======================================================================
!
!     J a c g : Compute np Gauss points, which are the zeros of the
!               Jacobi polynomial with parameter alpha and beta.
!
!=======================================================================
!
!     Note :
!     ----
!          .Alpha and Beta determines the specific type of gauss points.
!                  .alpha = beta =  0.0  ->  Legendre points
!                  .alpha = beta = -0.5  ->  Chebyshev points
!
!=======================================================================
!
  implicit none

  integer np
  double precision alpha,beta
  double precision xjac(np)

  integer k,j,i,jmin,jm,n
  double precision xlast,dth,x,x1,x2,recsum,delx,xmin,swap
  double precision p,pd,pm1,pdm1,pm2,pdm2

  integer, parameter :: kstop = 10
  double precision, parameter :: zero = 0.d0, eps = 1.0d-12

!
!-----------------------------------------------------------------------
!
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  xlast = 0.d0
  n   = np-1
  dth = 4.d0*datan(1.d0)/(2.d0*dble(n)+2.d0)
  p = 0.d0
  pd = 0.d0
  jmin = 0
  do 40 j=1,np
   if (j == 1) then
      x = dcos((2.d0*(dble(j)-1.d0)+1.d0)*dth)
   else
      x1 = dcos((2.d0*(dble(j)-1.d0)+1.d0)*dth)
      x2 = xlast
      x  = (x1+x2)/2.d0
   endif
   do 30 k=1,kstop
      call jacobf (p,pd,pm1,pdm1,pm2,pdm2,np,alpha,beta,x)
      recsum = 0.d0
      jm = j-1
      do 29 i=1,jm
         recsum = recsum+1.d0/(x-xjac(np-i+1))
 29         continue
      delx = -p/(pd-recsum*p)
      x    = x+delx
      if (abs(delx) < eps) goto 31
 30      continue
 31      continue
   xjac(np-j+1) = x
   xlast        = x
 40   continue
  do 200 i=1,np
   xmin = 2.d0
   do 100 j=i,np
      if (xjac(j) < xmin) then
         xmin = xjac(j)
         jmin = j
      endif
 100     continue
   if (jmin /= i) then
      swap = xjac(i)
      xjac(i) = xjac(jmin)
      xjac(jmin) = swap
   endif
 200  continue
  return
  end subroutine jacg
