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

  double precision function endw1 (n,alpha,beta)
!
!=======================================================================
!
!     E n d w 1 :
!     ---------
!
!=======================================================================
!
  implicit none

  integer n
  double precision alpha,beta

  double precision, parameter :: zero=0.d0,one=1.d0,two=2.d0, &
        three=3.d0,four=4.d0

  double precision apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3

  double precision, external :: gammaf

  integer i
!
!-----------------------------------------------------------------------
!
  f3 = zero
  apb   = alpha+beta
  if (n == 0) then
   endw1 = zero
   return
  endif
  f1   = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
  f1   = f1*(apb+two)*two**(apb+two)/two
  if (n == 1) then
   endw1 = f1
   return
  endif
  fint1 = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
  fint1 = fint1*two**(apb+two)
  fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
  fint2 = fint2*two**(apb+three)
  f2    = (-two*(beta+two)*fint1 + (apb+four)*fint2) * (apb+three)/four
  if (n == 2) then
   endw1 = f2
   return
  endif
  do i=3,n
   di   = dble(i-1)
   abn  = alpha+beta+di
   abnn = abn+di
   a1   = -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
   a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
   a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
   f3   =  -(a2*f2+a1*f1)/a3
   f1   = f2
   f2   = f3
  enddo
  endw1  = f3

  return
  end function endw1
