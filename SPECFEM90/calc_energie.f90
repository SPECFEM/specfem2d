
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

  subroutine calc_energie(hprime,hTprime,ibool,displ,veloc, &
       Uxnewloc,Uznewloc,kmato,dvolu,xjaci,density,elastcoef,wx,wy, &
       nxgll,npoin,ndime,nspec,numat)

  use timeparams
  use energie

  implicit none

  double precision, parameter :: zero = 0.d0, two = 2.d0

  integer nxgll,nspec,ndime,npoin,numat

  double precision hprime(nxgll,nxgll),hTprime(nxgll,nxgll)
  double precision Uxnewloc(nxgll,nxgll,nspec)
  double precision Uznewloc(nxgll,nxgll,nspec)

  double precision dUx_dxi,dUz_dxi,dUx_deta,dUz_deta
  double precision hTprimex,hprimez

  double precision density(numat),elastcoef(4,numat)
  double precision dvolu(nspec,nxgll,nxgll)
  double precision xjaci(nspec,ndime,ndime,nxgll,nxgll)
  double precision wx(nxgll),wy(nxgll)

  integer ibool(nxgll,nxgll,nspec)
  integer kmato(nspec)

  double precision displ(ndime,npoin),veloc(ndime,npoin)

  integer i,j,k,l,iglobnum,material
  double precision energie_pot,energie_cin
  double precision dxux,dzux,dxuz,dzuz
  double precision rKmod,rlamda,rmu,xix,xiz,etax,etaz,denst,rjacob

! map the global displacement field to the local mesh
!$PAR DOALL
!$PAR&      READONLY(ibool,displ)
  do k=1,nspec
    do j=1,nxgll
    do i=1,nxgll
                iglobnum = ibool(i,j,k)
                Uxnewloc(i,j,k) = displ(1,iglobnum)
                Uznewloc(i,j,k) = displ(2,iglobnum)
          enddo
    enddo
  enddo

  energie_pot = zero
  energie_cin = zero

! this loop is simply a reduction
! on the two scalar variables "energie_cin" and "energie_pot"
!$PAR DOALL_REDUCTION
  do k=1,nspec

! get the elastic parameters
 material = kmato(k)

 rlamda = elastcoef(1,material)
 rmu    = elastcoef(2,material)
 rKmod  = elastcoef(3,material)
 denst  = density(material)

    do j=1,nxgll
    do i=1,nxgll

! compute the gradient of the displacement field (matrix products)
    dUx_dxi = zero
    dUz_dxi = zero
    dUx_deta = zero
    dUz_deta = zero

          do l=1,nxgll

                hTprimex = hTprime(i,l)
                hprimez = hprime(l,j)

                dUx_dxi = dUx_dxi + hTprimex*Uxnewloc(l,j,k)
                dUz_dxi = dUz_dxi + hTprimex*Uznewloc(l,j,k)
                dUx_deta = dUx_deta + Uxnewloc(i,l,k)*hprimez
                dUz_deta = dUz_deta + Uznewloc(i,l,k)*hprimez

          enddo

! apply the chain rule to get this gradient in the physical domain
  xix = xjaci(k,1,1,i,j)
  xiz = xjaci(k,1,2,i,j)
  etax = xjaci(k,2,1,i,j)
  etaz = xjaci(k,2,2,i,j)
  rjacob = dvolu(k,i,j)

  dxux = dUx_dxi*xix + dUx_deta*etax
  dzux = dUx_dxi*xiz + dUx_deta*etaz

  dxuz = dUz_dxi*xix + dUz_deta*etax
  dzuz = dUz_dxi*xiz + dUz_deta*etaz

  iglobnum = ibool(i,j,k)

! calcul de l'energie cinetique
  energie_cin = energie_cin + &
          denst*(veloc(1,iglobnum)**2 + veloc(2,iglobnum)**2) &
          *wx(i)*wy(j)*rjacob

! calcul de l'energie potentielle elastique
  energie_pot = energie_pot + &
          (rKmod*dxux**2 + rKmod*dzuz**2 + two*rlamda*dxux*dzuz + &
          rmu*(dzux + dxuz)**2)*wx(i)*wy(j)*rjacob

          enddo
    enddo
  enddo

! do not forget to divide by two at the end
  energie_cin = energie_cin / two
  energie_pot = energie_pot / two

! on sauvegarde aussi l'energie totale qui doit etre constante
! au cours du temps (une fois que la source a fini d'agir)
! en l'absence de bords absorbants
! et decroitre au cours du temps en presence de bords absorbants
  write(ienergy,*) sngl(time),sngl(energie_cin),sngl(energie_pot), &
                sngl(energie_cin + energie_pot)

  return
  end subroutine calc_energie
