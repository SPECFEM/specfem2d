
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

  subroutine compute_gradient_attenuation(displ,duxdxl,duzdxl,duxdzl,duzdzl,xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,NSPEC,npoin)

! compute Grad(displ) for attenuation

  implicit none

  include "constants.h"

  integer NSPEC,npoin

  integer, dimension(NGLLX,NGLLZ,NSPEC) :: ibool

  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: duxdxl,duzdxl,duxdzl,duzdzl,xix,xiz,gammax,gammaz

  double precision, dimension(NDIME,npoin) :: displ

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz

! local variables
  integer i,j,k,ispec,iglob

! space derivatives
  double precision tempx1l,tempx2l,tempz1l,tempz2l
  double precision hp1,hp2

! jacobian
  double precision xixl,xizl,gammaxl,gammazl

! loop over spectral elements
  do ispec = 1,NSPEC

! double loop over GLL to compute and store gradients
    do j = 1,NGLLZ
      do i = 1,NGLLX

! derivative along x
          tempx1l = ZERO
          tempz1l = ZERO
          do k = 1,NGLLX
            hp1 = hprime_xx(k,i)
            iglob = ibool(k,j,ispec)
            tempx1l = tempx1l + displ(1,iglob)*hp1
            tempz1l = tempz1l + displ(2,iglob)*hp1
          enddo

! derivative along z
          tempx2l = ZERO
          tempz2l = ZERO
          do k = 1,NGLLZ
            hp2 = hprime_zz(k,j)
            iglob = ibool(i,k,ispec)
            tempx2l = tempx2l + displ(1,iglob)*hp2
            tempz2l = tempz2l + displ(2,iglob)*hp2
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

! derivatives of displacement
          duxdxl(i,j,ispec) = tempx1l*xixl + tempx2l*gammaxl
          duxdzl(i,j,ispec) = tempx1l*xizl + tempx2l*gammazl

          duzdxl(i,j,ispec) = tempz1l*xixl + tempz2l*gammaxl
          duzdzl(i,j,ispec) = tempz1l*xizl + tempz2l*gammazl

      enddo
    enddo
  enddo

  end subroutine compute_gradient_attenuation

