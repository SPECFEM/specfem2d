
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

  subroutine define_derivative_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz)

  implicit none

  include "constants.h"

! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll

! weights
  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLZ) :: wzgll

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz

! function for calculating derivatives of Lagrange polynomials
  double precision, external :: lagrange_deriv_GLL

  integer i1,i2,k1,k2

! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly zero
  if(mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = ZERO
  if(mod(NGLLZ,2) /= 0) zigll((NGLLZ-1)/2+1) = ZERO

! calculate derivatives of the Lagrange polynomials
! and precalculate some products in double precision
! hprime(i,j) = h'_i(xigll_j) by definition of the derivative matrix
  do i1=1,NGLLX
    do i2=1,NGLLX
      hprime_xx(i1,i2) = lagrange_deriv_GLL(i1-1,i2-1,xigll,NGLLX)
    enddo
  enddo

  do k1=1,NGLLZ
    do k2=1,NGLLZ
      hprime_zz(k1,k2) = lagrange_deriv_GLL(k1-1,k2-1,zigll,NGLLZ)
    enddo
  enddo

  end subroutine define_derivative_matrices

