
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.0
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) May 2004
!
!========================================================================

  subroutine qsumspec(hprime,hTprime, &
         a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13x,a13z, &
         ibool,displ,veloc,accel,Uxnewloc,Uznewloc, &
         rmass,npoin,nspec,gltfu,initialfield, &
         is_bordabs,nelemabs,anyabs,time)

  implicit none

  include "constants.h"

  integer npoin,nspec,nelemabs
  logical anyabs

  double precision hprime(NGLLX,NGLLX),hTprime(NGLLX,NGLLX)
  double precision a1(NGLLX,NGLLX,nspec),a2(NGLLX,NGLLX,nspec), &
     a3(NGLLX,NGLLX,nspec),a4(NGLLX,NGLLX,nspec),a5(NGLLX,NGLLX,nspec), &
     a6(NGLLX,NGLLX,nspec),a7(NGLLX,NGLLX,nspec), &
     a8(NGLLX,NGLLX,nspec),a9(NGLLX,NGLLX,nspec),a10(NGLLX,NGLLX,nspec)
  double precision a11(NGLLX,NGLLX),a12(NGLLX,NGLLX)
  double precision a13x(NGLLX,NGLLX,nelemabs)
  double precision a13z(NGLLX,NGLLX,nelemabs)
  double precision Uxnewloc(NGLLX,NGLLX,nspec)
  double precision Uznewloc(NGLLX,NGLLX,nspec)

  integer is_bordabs(nspec)

! local arrays
  double precision Uxoldloc(NGLLX,NGLLX)
  double precision Uzoldloc(NGLLX,NGLLX)
  double precision t1(NGLLX,NGLLX)
  double precision t2(NGLLX,NGLLX)
  double precision t3(NGLLX,NGLLX)
  double precision t4(NGLLX,NGLLX)

  double precision dUx_dxi,dUz_dxi,dUx_deta,dUz_deta
  double precision hprimex,hTprimex,hprimez,hTprimez

  integer ibool(NGLLX,NGLLX,nspec)

  double precision rmass(npoin)
  double precision displ(NDIME,npoin),veloc(NDIME,npoin),accel(NDIME,npoin)

  double precision gltfu(20)

  integer i,j,k,l,ielems,iglobsource,iglobnum,numer_abs
  double precision ricker,time
  logical initialfield

  double precision f0,t0,factor,a,angleforce

! main loop on all the spectral elements
  do k=1,nspec

! map the global displacement field to the local mesh
    do j=1,NGLLX
      do i=1,NGLLX
        iglobnum = ibool(i,j,k)
        Uxoldloc(i,j) = displ(1,iglobnum)
        Uzoldloc(i,j) = displ(2,iglobnum)
      enddo
    enddo

    do j=1,NGLLX
    do i=1,NGLLX

! compute the gradient of the displacement field (matrix products)
    dUx_dxi  = zero
    dUz_dxi  = zero
    dUx_deta = zero
    dUz_deta = zero

      do l=1,NGLLX

        hTprimex = hTprime(i,l)
        hprimez = hprime(l,j)

        dUx_dxi = dUx_dxi + hTprimex*Uxoldloc(l,j)
        dUz_dxi = dUz_dxi + hTprimex*Uzoldloc(l,j)
        dUx_deta = dUx_deta + Uxoldloc(i,l)*hprimez
        dUz_deta = dUz_deta + Uzoldloc(i,l)*hprimez

      enddo

! compute the local arrays using the components of the stiffness matrix
  t1(i,j) = a1(i,j,k)*dUx_dxi + a2(i,j,k)*dUx_deta + &
                     a3(i,j,k)*dUz_dxi + a4(i,j,k)*dUz_deta
  t2(i,j) = a2(i,j,k)*dUx_dxi + a6(i,j,k)*dUx_deta + &
                     a7(i,j,k)*dUz_dxi + a8(i,j,k)*dUz_deta
  t3(i,j)=  a3(i,j,k)*dUx_dxi + a7(i,j,k)*dUx_deta + &
                     a9(i,j,k)*dUz_dxi + a10(i,j,k)*dUz_deta
  t4(i,j)=  a4(i,j,k)*dUx_dxi + a8(i,j,k)*dUx_deta + &
                     a10(i,j,k)*dUz_dxi + a5(i,j,k)*dUz_deta

      enddo
    enddo

! compute the local forces (sum of two matrix products)
    do j=1,NGLLX
      do i=1,NGLLX

        Uxnewloc(i,j,k) = zero
        Uznewloc(i,j,k) = zero

        do l=1,NGLLX
          hprimex = hprime(i,l)
          hTprimez = hTprime(l,j)

          Uxnewloc(i,j,k) = Uxnewloc(i,j,k) + hprimex*t1(l,j) + t2(i,l)*hTprimez
          Uznewloc(i,j,k) = Uznewloc(i,j,k) + hprimex*t3(l,j) + t4(i,l)*hTprimez

        enddo

      enddo
    enddo

! conditions absorbantes nouvelle formulation
! pas de dependance par l'adressage indirect
! car chaque element absorbant est mappe sur un element spectral different
  if(anyabs) then
    numer_abs = is_bordabs(k)
    if(numer_abs .gt. 0) then
      do j=1,NGLLX
      do i=1,NGLLX
        if(a13x(i,j,numer_abs) .ne. zero) then
          iglobnum = ibool(i,j,k)
          Uxnewloc(i,j,k) = Uxnewloc(i,j,k) - a13x(i,j,numer_abs)*veloc(1,iglobnum)
          Uznewloc(i,j,k) = Uznewloc(i,j,k) - a13z(i,j,numer_abs)*veloc(2,iglobnum)
        endif
      enddo
      enddo
    endif
  endif

! assemblage des contributions des differents elements
    do j=1,NGLLX
      do i=1,NGLLX
        iglobnum = ibool(i,j,k)
        accel(1,iglobnum) = accel(1,iglobnum) + Uxnewloc(i,j,k)
        accel(2,iglobnum) = accel(2,iglobnum) + Uznewloc(i,j,k)
      enddo
    enddo

  enddo

! --- add the source
  if(.not. initialfield) then

  f0 = gltfu(5)
  t0 = gltfu(6)
  factor = gltfu(7)
  angleforce = gltfu(8)

! Ricker wavelet for the source time function
  a = pi*pi*f0*f0
  ricker = - factor * (1.d0-2.d0*a*(time-t0)**2)*exp(-a*(time-t0)**2)

! --- collocated force
  if(nint(gltfu(2))  ==  1) then
    iglobsource = nint(gltfu(9))
    accel(1,iglobsource) = accel(1,iglobsource) - dsin(angleforce)*ricker
    accel(2,iglobsource) = accel(2,iglobsource) + dcos(angleforce)*ricker

!---- explosion
  else if(nint(gltfu(2))  ==  2) then
!   recuperer numero d'element de la source
    ielems = nint(gltfu(12))
    do i=1,NGLLX
      do j=1,NGLLX
        iglobnum = ibool(i,j,ielems)
        accel(1,iglobnum) = accel(1,iglobnum) + a11(i,j)*ricker
        accel(2,iglobnum) = accel(2,iglobnum) + a12(i,j)*ricker
      enddo
    enddo
  endif

  else
    stop 'wrong source type'
  endif

! --- multiplier par l'inverse de la matrice de masse
  accel(1,:) =  accel(1,:)*rmass(:)
  accel(2,:) =  accel(2,:)*rmass(:)

  end subroutine qsumspec

