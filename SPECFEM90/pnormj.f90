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

  double precision function pnormj (n,alpha,beta)
!
!=======================================================================
!
!     P n o r m j
!     -----------
!
!=======================================================================
!
  implicit none

  double precision alpha,beta
  integer n

  double precision one,two,dn,const,prod,dindx,frac
  double precision, external :: gammaf
  integer i

  one   = 1.d0
  two   = 2.d0
  dn    = dble(n)
  const = alpha+beta+one
  if (n <= 1) then
   prod   = gammaf(dn+alpha)*gammaf(dn+beta)
   prod   = prod/(gammaf(dn)*gammaf(dn+alpha+beta))
   pnormj = prod * two**const/(two*dn+const)
   return
  endif
  prod  = gammaf(alpha+one)*gammaf(beta+one)
  prod  = prod/(two*(one+const)*gammaf(const+one))
  prod  = prod*(one+alpha)*(two+alpha)
  prod  = prod*(one+beta)*(two+beta)
  do 100 i=3,n
   dindx = dble(i)
   frac  = (dindx+alpha)*(dindx+beta)/(dindx*(dindx+alpha+beta))
   prod  = prod*frac
 100  continue
  pnormj = prod * two**const/(two*dn+const)

  return
  end function pnormj
