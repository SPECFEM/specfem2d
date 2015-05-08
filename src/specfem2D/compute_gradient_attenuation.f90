
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

  subroutine compute_gradient_attenuation(displ_elastic,dux_dxl,duz_dxl,dux_dzl,duz_dzl, &
         xix,xiz,gammax,gammaz,ibool,elastic,hprime_xx,hprime_zz,nspec,nglob)

  use specfem_par, only: AXISYM,is_on_the_axis,hprimeBar_xx

! compute Grad(displ_elastic) for attenuation

  implicit none

  include "constants.h"

  integer :: nspec,nglob

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

  logical, dimension(nspec) :: elastic

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec)  :: xix,xiz,gammax,gammaz

  real(kind=CUSTOM_REAL), dimension(3,nglob) :: displ_elastic

! array with derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx!,hprimeBar_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

! local variables
  integer :: i,j,k,ispec

! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma

! jacobian
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl

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

          if(AXISYM) then
            if (is_on_the_axis(ispec)) then
              do k = 1,NGLJ
                dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
                duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)
              enddo
            else
              do k = 1,NGLJ
                dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
                duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
                dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
                duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)
              enddo
            endif
          else
            do k = 1,NGLLX
              dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
              duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
              dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
              duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)
            enddo
          endif

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

! derivatives of displacement
          dux_dxl(i,j,ispec) = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl(i,j,ispec) = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl(i,j,ispec) = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl(i,j,ispec) = duz_dxi*xizl + duz_dgamma*gammazl

          if (AXISYM .and. is_on_the_axis(ispec) .and. i == 1) then ! d_uz/dr=0 on the axis
            duz_dxl(i,j,ispec) = 0.d0
          endif

        enddo
      enddo

    endif

  enddo

  end subroutine compute_gradient_attenuation

