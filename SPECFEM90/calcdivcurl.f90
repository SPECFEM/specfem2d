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

  subroutine calcdivcurl(displ,div,curl,hprime,hTprime,ibool, &
        Uxloc,Uzloc,dUx_dxi,dUz_dxi,dUx_deta,dUz_deta,xjaci)
!
!=======================================================================
!
!     "c a l c d i v c u r l" : Compute the divergence and the curl
!                               of the displacement field
!
!=======================================================================
!
  use mesh01
  use spela202

  implicit none

  double precision Uxloc(nxgll,nxgll,nspec)
  double precision Uzloc(nxgll,nxgll,nspec)
  double precision hprime(nxgll,nxgll)
  double precision hTprime(nxgll,nxgll)
  double precision xjaci(nspec,ndime,ndime,nxgll,nxgll)
  double precision dUx_dxi(nxgll,nxgll,nspec)
  double precision dUz_dxi(nxgll,nxgll,nspec)
  double precision dUx_deta(nxgll,nxgll,nspec)
  double precision dUz_deta(nxgll,nxgll,nspec)
  double precision displ(ndime,npoin)
  double precision div(npoin)
  double precision curl(npoin)
  integer ibool(nxgll,nxgll,nspec)

  integer i,j,k,l,iglobnum
  double precision xix,xiz,etax,etaz


! definir div et curl

  do i=1,nxgll
      do j=1,nxgll
            do k=1,nspec
              iglobnum = ibool(i,j,k)
              Uxloc(i,j,k) = displ(1,iglobnum)
              Uzloc(i,j,k) = displ(2,iglobnum)
            enddo
      enddo
  enddo

  do k=1,nspec
   do i=1,nxgll
      do j=1,nxgll
      dUx_dxi(i,j,k) = 0.d0
      dUz_dxi(i,j,k) = 0.d0
      dUx_deta(i,j,k) = 0.d0
      dUz_deta(i,j,k) = 0.d0
            do l=1,nxgll

      dUx_dxi(i,j,k) = dUx_dxi(i,j,k) + hTprime(i,l)*Uxloc(l,j,k)
      dUz_dxi(i,j,k) = dUz_dxi(i,j,k) + hTprime(i,l)*Uzloc(l,j,k)
      dUx_deta(i,j,k) = dUx_deta(i,j,k) + Uxloc(i,l,k)*hprime(l,j)
      dUz_deta(i,j,k) = dUz_deta(i,j,k) + Uzloc(i,l,k)*hprime(l,j)

               enddo
            enddo
      enddo
  enddo

  do k=1,nspec
   do i=1,nxgll
      do j=1,nxgll

  xix = xjaci(k,1,1,i,j)
  xiz = xjaci(k,1,2,i,j)
  etax = xjaci(k,2,1,i,j)
  etaz = xjaci(k,2,2,i,j)

  iglobnum = ibool(i,j,k)

  div(iglobnum) = dUx_dxi(i,j,k)*xix + dUx_deta(i,j,k)*etax + &
        dUz_dxi(i,j,k)*xiz + dUz_deta(i,j,k)*etaz
  curl(iglobnum) = dUx_dxi(i,j,k)*xiz + dUx_deta(i,j,k)*etaz - &
        dUz_dxi(i,j,k)*xix - dUz_deta(i,j,k)*etax

      enddo
      enddo
  enddo

  return
  end subroutine calcdivcurl
