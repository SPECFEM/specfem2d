
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) January 2005
!
!========================================================================

  subroutine compute_gradient_fluid(potential,veloc_field_postscript, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,NSPEC,npoin)

! compute Grad(potential) in fluid medium

  implicit none

  include "constants.h"

  integer NSPEC,npoin

  integer, dimension(NGLLX,NGLLZ,NSPEC) :: ibool

  double precision, dimension(NGLLX,NGLLZ,NSPEC) :: xix,xiz,gammax,gammaz

! for compatibility with elastic arrays, potential is declared as a vector but correctly used below as a scalar
  double precision, dimension(NDIM,npoin) :: potential,veloc_field_postscript

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz

! local variables
  integer i,j,k,ispec,iglob

! space derivatives
  double precision tempx1l,tempx2l
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
          do k = 1,NGLLX
            hp1 = hprime_xx(k,i)
            iglob = ibool(k,j,ispec)
            tempx1l = tempx1l + potential(1,iglob)*hp1
          enddo

! derivative along z
          tempx2l = ZERO
          do k = 1,NGLLZ
            hp2 = hprime_zz(k,j)
            iglob = ibool(i,k,ispec)
            tempx2l = tempx2l + potential(1,iglob)*hp2
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

! derivatives of velocity potential
          iglob = ibool(i,j,ispec)
          veloc_field_postscript(1,iglob) = tempx1l*xixl + tempx2l*gammaxl
          veloc_field_postscript(2,iglob) = tempx1l*xizl + tempx2l*gammazl

      enddo
    enddo
  enddo

  end subroutine compute_gradient_fluid

