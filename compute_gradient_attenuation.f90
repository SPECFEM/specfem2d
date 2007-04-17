
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!                         Dimitri Komatitsch
!                     University of Pau, France
!
!                          (c) April 2007
!
!========================================================================

  subroutine compute_gradient_attenuation(displ_elastic,dux_dxl,duz_dxl,dux_dzl,duz_dzl, &
         xix,xiz,gammax,gammaz,ibool,elastic,hprime_xx,hprime_zz,nspec,npoin)

! compute Grad(displ_elastic) for attenuation

  implicit none

  include "constants.h"

  integer :: nspec,npoin

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

  logical, dimension(nspec) :: elastic

  double precision, dimension(NGLLX,NGLLZ,nspec) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl,xix,xiz,gammax,gammaz

  double precision, dimension(NDIM,npoin) :: displ_elastic

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz

! local variables
  integer :: i,j,k,ispec

! spatial derivatives
  double precision :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma

! jacobian
  double precision :: xixl,xizl,gammaxl,gammazl

! loop over spectral elements
  do ispec = 1,nspec

!---
!--- elastic spectral element
!---
    if(elastic(ispec)) then

! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX

! derivative along x and along z
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(k,i)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(k,i)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(k,j)
            duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(k,j)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

! derivatives of displacement
          dux_dxl(i,j,ispec) = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl(i,j,ispec) = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl(i,j,ispec) = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl(i,j,ispec) = duz_dxi*xizl + duz_dgamma*gammazl

        enddo
      enddo

    endif

  enddo

  end subroutine compute_gradient_attenuation

