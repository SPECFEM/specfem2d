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

  subroutine calcforce(F,ndime,gltfu,nltfl,t)

! calcul de la force source en temps

  implicit none

  integer ndime,nltfl
  double precision t
  double precision gltfu(20,nltfl)
  double precision F(ndime,nltfl)

  integer n,isource
  double precision funct,angle
  double precision, external :: ricker,dirac

  do n=1,nltfl

! determiner type de source
  isource = nint(gltfu(1,n))

! la source est une force colloquee
  if(nint(gltfu(2,n)) == 1) then

! introduire source suivant son type
  if(isource == 6) then
      funct = ricker(t,n,gltfu,nltfl)
  else if(isource == 7) then
      funct = dirac(t,n,gltfu,nltfl)
  else
      funct = 0.d0
  endif

  angle = gltfu(8,n)
  F(1,n) = - dsin(angle) * funct
  F(2,n) = + dcos(angle) * funct

  endif

  enddo

  return
  end subroutine calcforce
