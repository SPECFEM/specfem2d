!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================


  subroutine setup_GLL_points()

! computes shape functions and their derivatives for SEM grid

  use specfem_par

  implicit none

  ! local parameters
  integer :: i,j,ier

  ! set up Gauss-Lobatto-Legendre points, weights and also derivation matrices (hprime_**,..)
  call define_derivation_matrices()

  if (AXISYM) then
    ! set up Gauss-Lobatto-Jacobi points, weights and also derivation matrices
    call define_GLJ_derivation_matrix()
  endif

  ! shape arrays
  allocate(shape2D(ngnod,NGLLX,NGLLZ),dershape2D(NDIM,ngnod,NGLLX,NGLLZ),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating shape arrays')

  do j = 1,NGLLZ
    do i = 1,NGLLX
      call define_shape_functions(shape2D(:,i,j),dershape2D(:,:,i,j),xigll(i),zigll(j),ngnod)
    enddo
  enddo

  ! allocate 1-D Lagrange interpolators and derivatives
  ! for source and receivers
  allocate(hxir(NGLLX), &
           hxis(NGLLX), &
           hpxir(NGLLX), &
           hpxis(NGLLX), &
           hgammar(NGLLZ), &
           hgammas(NGLLZ), &
           hpgammar(NGLLZ), &
           hpgammas(NGLLZ),stat=ier)
  if (ier /= 0) call stop_the_code('error allocating arrays for interpolators')

  end subroutine setup_GLL_points

