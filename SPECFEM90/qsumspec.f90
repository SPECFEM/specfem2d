
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

  subroutine qsumspec(hprime,hTprime, &
         a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13x,a13z,force, &
         ibool,displ,veloc,accel,Uxnewloc,Uznewloc, &
         rmass,nxgll,npoin,ndime,nspec,gltfu,nltfl,initialfield, &
         numabs,is_bordabs,nelemabs,anyabs)

  use timeparams

  implicit none

  integer nxgll,npoin,ndime,nspec,nltfl,nelemabs
  logical anyabs

  double precision hprime(nxgll,nxgll),hTprime(nxgll,nxgll)
  double precision a1(nxgll,nxgll,nspec),a2(nxgll,nxgll,nspec), &
     a3(nxgll,nxgll,nspec),a4(nxgll,nxgll,nspec),a5(nxgll,nxgll,nspec), &
         a6(nxgll,nxgll,nspec),a7(nxgll,nxgll,nspec), &
     a8(nxgll,nxgll,nspec),a9(nxgll,nxgll,nspec),a10(nxgll,nxgll,nspec)
  double precision a11(nxgll,nxgll,nltfl),a12(nxgll,nxgll,nltfl)
  double precision a13x(nxgll,nxgll,nelemabs)
  double precision a13z(nxgll,nxgll,nelemabs)
  double precision Uxnewloc(nxgll,nxgll,nspec)
  double precision Uznewloc(nxgll,nxgll,nspec)

  integer numabs(nelemabs)
  integer is_bordabs(nspec)

! petits tableaux locaux (could be suppressed if needed)
! maxnxgll est la valeur maximale possible du degre polynomial (10 par exemple)
  integer, parameter :: maxnxgll = 10
  double precision Uxoldloc(maxnxgll,maxnxgll)
  double precision Uzoldloc(maxnxgll,maxnxgll)
  double precision t1(maxnxgll,maxnxgll)
  double precision t2(maxnxgll,maxnxgll)
  double precision t3(maxnxgll,maxnxgll)
  double precision t4(maxnxgll,maxnxgll)

  double precision dUx_dxi,dUz_dxi,dUx_deta,dUz_deta
  double precision hprimex,hTprimex,hprimez,hTprimez

  integer ibool(nxgll,nxgll,nspec)

  double precision rmass(npoin)
  double precision force(ndime,nltfl)
  double precision displ(ndime,npoin),veloc(ndime,npoin),accel(ndime,npoin)

  double precision gltfu(20,nltfl)

  double precision, external :: dirac,ricker

  integer i,j,k,l,n,isource,ielems,iglobsource,iglobnum,ip,numer_abs
  double precision sig
  logical initialfield

  double precision, parameter :: zero=0.d0

! main loop on all the spectral elements
  do k=1,nspec

! map the global displacement field to the local mesh
    do j=1,nxgll
    do i=1,nxgll
                iglobnum = ibool(i,j,k)
                Uxoldloc(i,j) = displ(1,iglobnum)
                Uzoldloc(i,j) = displ(2,iglobnum)
          enddo
    enddo

    do j=1,nxgll
    do i=1,nxgll

! compute the gradient of the displacement field (matrix products)
    dUx_dxi  = zero
    dUz_dxi  = zero
    dUx_deta = zero
    dUz_deta = zero

          do l=1,nxgll

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
    do j=1,nxgll
    do i=1,nxgll
          Uxnewloc(i,j,k) = zero
          Uznewloc(i,j,k) = zero

          do l=1,nxgll

                hprimex = hprime(i,l)
                hTprimez = hTprime(l,j)

                Uxnewloc(i,j,k) = Uxnewloc(i,j,k) + &
                          hprimex*t1(l,j) + t2(i,l)*hTprimez
                Uznewloc(i,j,k) = Uznewloc(i,j,k) + &
                          hprimex*t3(l,j) + t4(i,l)*hTprimez

          enddo

    enddo
    enddo

! conditions absorbantes nouvelle formulation
! pas de dependance par l'adressage indirect
! car chaque element absorbant est mappe sur un element spectral different
  if(anyabs) then
    numer_abs = is_bordabs(k)
    if(numer_abs .gt. 0) then
        do j=1,nxgll
        do i=1,nxgll
        if(a13x(i,j,numer_abs) .ne. zero) then
              iglobnum = ibool(i,j,k)
              Uxnewloc(i,j,k) = Uxnewloc(i,j,k) - &
                   a13x(i,j,numer_abs)*veloc(1,iglobnum)
              Uznewloc(i,j,k) = Uznewloc(i,j,k) - &
                   a13z(i,j,numer_abs)*veloc(2,iglobnum)
          endif
        enddo
        enddo
    endif
  endif

! assemblage des contributions des differents elements
    do j=1,nxgll
    do i=1,nxgll
            iglobnum = ibool(i,j,k)
            accel(1,iglobnum) = accel(1,iglobnum) + Uxnewloc(i,j,k)
            accel(2,iglobnum) = accel(2,iglobnum) + Uznewloc(i,j,k)
    enddo
    enddo

  enddo

! --- ajouter sources forces colloquees

  if(.not. initialfield) then
    do n=1,nltfl
      iglobsource = nint(gltfu(9,n))
      accel(:,iglobsource) = accel(:,iglobsource) + force(:,n)
    enddo
  endif

!---- ajouter sources explosives

  if(.not. initialfield) then

  do n=1,nltfl

! seulement si source explosive
  if(nint(gltfu(2,n))  ==  2) then

! determiner type de source en temps
  isource = nint(gltfu(1,n))

! introduire source suivant son type
  if(isource  ==  6) then
    sig = ricker(time,n,gltfu,nltfl)
  else if(isource  ==  7) then
    sig = dirac(time,n,gltfu,nltfl)
  else
    sig = zero
  endif

! recuperer numero d'element de la source
  ielems = nint(gltfu(12,n))

 do i=1,nxgll
   do j=1,nxgll
     iglobnum = ibool(i,j,ielems)
     accel(1,iglobnum) = accel(1,iglobnum) + a11(i,j,n)*sig
     accel(2,iglobnum) = accel(2,iglobnum) + a12(i,j,n)*sig
   enddo
 enddo

  endif

  enddo

  endif

! --- multiplier par l'inverse de la matrice de masse

  accel(1,:) =  accel(1,:)*rmass(:)
  accel(2,:) =  accel(2,:)*rmass(:)

  return
  end subroutine qsumspec
