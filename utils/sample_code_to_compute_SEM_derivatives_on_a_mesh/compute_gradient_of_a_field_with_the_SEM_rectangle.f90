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

! Dimitri Komatitsch, CNRS, Marseille, August 2014

! this is a sample program to illustrate how to compute the gradient (in Cartesian coordinates)
! of a 2D displacement field using the spectral-element method, and check it at a given point in space;
! we define the mesh as a rectangle

  program compute_gradient

  implicit none

  include "constants.h"

! physical coordinates of the mesh
  double precision, parameter :: xmin =  +12.d0, xmax = +23.d0
  double precision, parameter :: zmin = -27.d0, zmax = +34.d0

! spectral elements
  integer, parameter :: NEX_XI = 280, NEX_ETA = 400

! physical location of the receiver at which we compute the gradient of the displacement vector below
  double precision, parameter :: x_receiver = 14.378d0, z_receiver = + 2.674d0

! "size" of the whole model and of each spectral element in the topologically-regular grid
  double precision, parameter :: sizex = xmax - xmin
  double precision, parameter :: sizez = zmax - zmin
  double precision, parameter :: deltax = sizex / NEX_XI
  double precision, parameter :: deltaz = sizez / NEX_ETA

! number of spectral elements
  integer, parameter :: nspec = NEX_XI * NEX_ETA

! intermediate variables needed for the calculation of NGLOB below
  integer, parameter :: number_of_unique_points_along_x = (NEX_XI * (NGLLX-1) + 1)
  integer, parameter :: number_of_unique_points_along_z = (NEX_ETA * (NGLLZ-1) + 1)

! number of unique points in the global vector of unknowns
  integer, parameter :: nglob = number_of_unique_points_along_x * number_of_unique_points_along_z

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

! array containing the geometrical coordinates of the GLL points of the whole mesh
  double precision, dimension(NDIM,nglob) :: coord

! geometrical coordinates of the anchor points of the elements of the mesh
  double precision, dimension(NDIM,NGNOD,NSPEC) :: coord_of_anchor_points

! the displacement vector (the goal of the program below is to compute its derivative)
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: displacement_vector

! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll

! weights
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  real(kind=CUSTOM_REAL), dimension(NGLLZ) :: wzgll

! array with derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: dux_dx,duz_dx,dux_dz,duz_dz
  real(kind=CUSTOM_REAL) :: dux_dx_interpolated,duz_dx_interpolated,dux_dz_interpolated,duz_dz_interpolated
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz,jacobian

! local variables
  integer :: i,j,k,ispec,iglob,ia
  integer :: ixe,ize,ix_in_global_grid,iz_in_global_grid
  integer :: ispec_selected_receiver, ix_initial_guess, iz_initial_guess, iter_loop

! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma

! Lagrange interpolants at the receiver
  double precision, dimension(NGLLX) :: hxir,hpxir
  double precision, dimension(NGLLZ) :: hgammar,hpgammar

! Jacobian matrix
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl,weight,area

  double precision :: x,z
  double precision :: xi,gamma
  double precision :: xi_in_global_grid,gamma_in_global_grid
  double precision :: dx,dz,dxi,dgamma,xi_receiver,gamma_receiver,hlagrange,final_distance
  double precision :: dist,distmin

! define the Gauss-Lobatto-Legendre (GLL) integration points, integration weights, and derivative matrix
  call define_derivation_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz)

! build the local-to-global and global-to-local topological mapping between
! the GLL points of each spectral element and the global unique points of the vector for which we want to compute the derivative
! loop on all the spectral elements of the grid to create

  ispec = 0
  ibool(:,:,:) = 0

  ! loop on all the spectral elements to create in the mesh
  do ize = 1,NEX_ETA
    do ixe = 1,NEX_XI

      ispec = ispec + 1

     ! loop on all the GLL points to create inside a given spectral element of the mesh
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ix_in_global_grid = (ixe-1)*(NGLLX-1) + i
          iz_in_global_grid = (ize-1)*(NGLLZ-1) + j
          ibool(i,j,ispec) = (iz_in_global_grid - 1)*number_of_unique_points_along_x + ix_in_global_grid
        enddo
      enddo

! build the grid of control nodes that define the geometry of the anchor points of the mesh elements
! in order to be able to honor the curvature of the Earth we use elements with 9 anchor points rather than 4
! (i.e. the geometry is defined with lower-order functions than the unknowns
!  of the problem; see for instance Chapter 16 of the finite-element course of
!  Carlos Felippa in Colorado for a discussion about this).
!  The 9 control nodes are defined as follows:
!
!                               4 . . . . 7 . . . . 3
!                               .                   .
!                               .         gamma     .
!                               .                   .
!                               8         9  xi     6
!                               .                   .
!                               .                   .
!                               .                   .
!                               1 . . . . 5 . . . . 2
!
!                           Local coordinate system : (xi,gamma)

! point 1
    ia = 1
    xi_in_global_grid = ((ixe-1)*deltax) / sizex
    gamma_in_global_grid = ((ize-1)*deltaz) / sizez
    coord_of_anchor_points(1,ia,ispec) = xmin + xi_in_global_grid * sizex
    coord_of_anchor_points(2,ia,ispec) = zmin + gamma_in_global_grid * sizez

! point 2
    ia = 2
    xi_in_global_grid = ((ixe-1)*deltax + deltax) / sizex
    gamma_in_global_grid = ((ize-1)*deltaz) / sizez
    coord_of_anchor_points(1,ia,ispec) = xmin + xi_in_global_grid * sizex
    coord_of_anchor_points(2,ia,ispec) = zmin + gamma_in_global_grid * sizez

! point 3
    ia = 3
    xi_in_global_grid = ((ixe-1)*deltax + deltax) / sizex
    gamma_in_global_grid = ((ize-1)*deltaz + deltaz) / sizez
    coord_of_anchor_points(1,ia,ispec) = xmin + xi_in_global_grid * sizex
    coord_of_anchor_points(2,ia,ispec) = zmin + gamma_in_global_grid * sizez

! point 4
    ia = 4
    xi_in_global_grid = ((ixe-1)*deltax) / sizex
    gamma_in_global_grid = ((ize-1)*deltaz + deltaz) / sizez
    coord_of_anchor_points(1,ia,ispec) = xmin + xi_in_global_grid * sizex
    coord_of_anchor_points(2,ia,ispec) = zmin + gamma_in_global_grid * sizez

! point 5
    ia = 5
    xi_in_global_grid = ((ixe-1)*deltax + deltax / 2) / sizex
    gamma_in_global_grid = ((ize-1)*deltaz) / sizez
    coord_of_anchor_points(1,ia,ispec) = xmin + xi_in_global_grid * sizex
    coord_of_anchor_points(2,ia,ispec) = zmin + gamma_in_global_grid * sizez

! point 6
    ia = 6
    xi_in_global_grid = ((ixe-1)*deltax + deltax) / sizex
    gamma_in_global_grid = ((ize-1)*deltaz + deltaz / 2) / sizez
    coord_of_anchor_points(1,ia,ispec) = xmin + xi_in_global_grid * sizex
    coord_of_anchor_points(2,ia,ispec) = zmin + gamma_in_global_grid * sizez

! point 7
    ia = 7
    xi_in_global_grid = ((ixe-1)*deltax + deltax / 2) / sizex
    gamma_in_global_grid = ((ize-1)*deltaz + deltaz) / sizez
    coord_of_anchor_points(1,ia,ispec) = xmin + xi_in_global_grid * sizex
    coord_of_anchor_points(2,ia,ispec) = zmin + gamma_in_global_grid * sizez

! point 8
    ia = 8
    xi_in_global_grid = ((ixe-1)*deltax) / sizex
    gamma_in_global_grid = ((ize-1)*deltaz + deltaz / 2) / sizez
    coord_of_anchor_points(1,ia,ispec) = xmin + xi_in_global_grid * sizex
    coord_of_anchor_points(2,ia,ispec) = zmin + gamma_in_global_grid * sizez

! point 9
    ia = 9
    xi_in_global_grid = ((ixe-1)*deltax + deltax / 2) / sizex
    gamma_in_global_grid = ((ize-1)*deltaz + deltaz / 2) / sizez
    coord_of_anchor_points(1,ia,ispec) = xmin + xi_in_global_grid * sizex
    coord_of_anchor_points(2,ia,ispec) = zmin + gamma_in_global_grid * sizez

    enddo
  enddo

  if (ispec /= NSPEC) stop 'the total number of spectral elements created is not correct'

! check that the numbering created is correct
  if (minval(ibool) /= 1 .or. maxval(ibool) /= NGLOB) stop 'the grid numbering created is not correct'

! print *,'x of anchor points min, max = ',minval(coord_of_anchor_points(1,:,:)),maxval(coord_of_anchor_points(1,:,:))
! print *,'z of anchor points min, max = ',minval(coord_of_anchor_points(2,:,:)),maxval(coord_of_anchor_points(2,:,:))
! print *

! define the geometrical coordinates of the points of the global grid
! and compute the 2D Jacobian at a given point in a 4-control-node or 9-control-node geometrical element
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX

          ! this is the local GLL point inside the current spectral element
          xi = xigll(i)
          gamma = zigll(j)

          ! compute the Jacobian matrix and the physical coordinates of this mesh point
          call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl,jacobianl, &
                    coord_of_anchor_points,ispec,nspec)

          coord(1,ibool(i,j,ispec)) = x
          coord(2,ibool(i,j,ispec)) = z

          xix(i,j,ispec) = xixl
          xiz(i,j,ispec) = xizl
          gammax(i,j,ispec) = gammaxl
          gammaz(i,j,ispec) = gammazl
          jacobian(i,j,ispec) = jacobianl

        enddo
      enddo
    enddo

  print *
  print *,'x of the mesh min, max = ',minval(coord(1,:)),maxval(coord(1,:))
  print *,'z of the mesh min, max = ',minval(coord(2,:)),maxval(coord(2,:))

! as a verification that the Jacobian matrix is fine, compute the area of the mesh
! by performing a numerical integration of the Jacobian and compare it to the exact area
! for the simple geometry that we use
  area = ZERO
  do ispec = 1,nspec
    do i=1,NGLLX
      do j=1,NGLLZ
        weight=wxgll(i)*wzgll(j)
        area = area + jacobian(i,j,ispec)*weight
      enddo
    enddo
  enddo
  print *
  print *,'calculated area of the mesh = ',sngl(area)
  print *,'exact area of the mesh = ',sngl(sizex * sizez)

! assign a test analytical displacement field, so that we can then check if the derivative below is fine
  do iglob = 1,nglob
    x = coord(1,iglob)
    z = coord(2,iglob)
    displacement_vector(1,iglob) = cos(x) + cos(z)
    displacement_vector(2,iglob) = cos(2*x) + cos(2*z)
  enddo

!---
!--- this holds for elastic spectral elements (for acoustic elements in the outer core the formulation would be with a scalar)
!---

! loop over all the spectral elements of the mesh
  do ispec = 1,nspec

! double loop over the GLL points of each spectral element to compute the gradient
      do j = 1,NGLLZ
        do i = 1,NGLLX

! derivative along the local (xi,gamma) coordinate system of a given spectral element
          dux_dxi = 0._CUSTOM_REAL
          duz_dxi = 0._CUSTOM_REAL

          dux_dgamma = 0._CUSTOM_REAL
          duz_dgamma = 0._CUSTOM_REAL

! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
! derivatives in the (xi,gamma) local coordinate system
! get the displacement field from the global vector of unknowns
! component 1 is Ux and component 2 is Uz
            dux_dxi = dux_dxi + displacement_vector(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displacement_vector(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displacement_vector(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displacement_vector(2,ibool(i,k,ispec))*hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

! use the chain rule and the Jacobian matrix to compute the derivatives of displacement in the Cartesian frame
! the notation we use here is d(Ux) / dz = dux_dz, d(Uz) / dx = duz_dx etc.
          dux_dx(i,j,ispec) = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dz(i,j,ispec) = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dx(i,j,ispec) = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dz(i,j,ispec) = duz_dxi*xizl + duz_dgamma*gammazl

        enddo
      enddo

  enddo

! find where the receiver is located in the mesh.
! in general it will NOT fall exactly on a GLL point and thus we need to solve a small nonlinear system
! in order to determine at which coordinates (xi,gamma) it is located between points

! set distance to huge initial value
  distmin = HUGEVAL

  do ispec = 1,nspec

! loop only on points inside the element
! exclude edges to ensure this point is not shared with other elements
     do j = 2,NGLLZ-1
        do i = 2,NGLLX-1

           iglob = ibool(i,j,ispec)
           dist = sqrt((x_receiver-dble(coord(1,iglob)))**2 + (z_receiver-dble(coord(2,iglob)))**2)

!          keep this point if it is closer to the receiver
           if (dist < distmin) then
              distmin = dist
              ispec_selected_receiver = ispec
              ix_initial_guess = i
              iz_initial_guess = j
           endif

        enddo
     enddo

! end of loop on all the spectral elements
  enddo

! now find the best (xi,gamma) for the receiver

! use initial guess in xi and gamma
  xi = xigll(ix_initial_guess)
  gamma = zigll(iz_initial_guess)

! iterate to solve the non linear system
  do iter_loop = 1,NUM_ITER

! recompute jacobian for the new point
    call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl,jacobianl, &
                    coord_of_anchor_points,ispec_selected_receiver,nspec)

! compute distance to target location
    dx = - (x - x_receiver)
    dz = - (z - z_receiver)

! compute increments
    dxi  = xixl*dx + xizl*dz
    dgamma = gammaxl*dx + gammazl*dz

! update values
    xi = xi + dxi
    gamma = gamma + dgamma

! impose that we stay in that element
! (useful if user gives a receiver outside the mesh for instance)
! we can go slightly outside the [1,1] segment since with finite elements
! the polynomial solution is defined everywhere
! this can be useful for convergence of itertive scheme with distorted elements
    if (xi > 1.10d0) xi = 1.10d0
    if (xi < -1.10d0) xi = -1.10d0
    if (gamma > 1.10d0) gamma = 1.10d0
    if (gamma < -1.10d0) gamma = -1.10d0

! end of nonlinear iterations
  enddo

! compute final coordinates of point found
    call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl,jacobianl, &
                    coord_of_anchor_points,ispec_selected_receiver,nspec)

! store xi,gamma of point found
  xi_receiver = xi
  gamma_receiver = gamma

! compute final distance between asked and found
  final_distance = sqrt((x_receiver-x)**2 + (z_receiver-z)**2)

  print *
  print *,'Location of the receiver:'

  if (final_distance == HUGEVAL) stop 'error locating the receiver'

  print *,'            original x: ',sngl(x_receiver)
  print *,'            original z: ',sngl(z_receiver)
  print *,'closest estimate found: ',sngl(final_distance),' m away'
  print *,' in element ',ispec_selected_receiver
  print *,' at xi,gamma coordinates = ',xi_receiver,gamma_receiver
  print *

! now that the receiver is located in the mesh, we need to compute the Lagrange interpolants
! (i.e. the interpolation matrix) at its (xi,gamma) location between grid points
! in order to be able to evaluate the derivative computed above there
! (for now it has been computed at the GLL points only)
  call lagrange_any(xi_receiver,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(gamma_receiver,NGLLZ,zigll,hgammar,hpgammar)

  dux_dx_interpolated = 0._CUSTOM_REAL
  dux_dz_interpolated = 0._CUSTOM_REAL
  duz_dx_interpolated = 0._CUSTOM_REAL
  duz_dz_interpolated = 0._CUSTOM_REAL

  do j = 1,NGLLZ
    do i = 1,NGLLX

      hlagrange = hxir(i)*hgammar(j)

      ! compute interpolated field
      dux_dx_interpolated = dux_dx_interpolated + dux_dx(i,j,ispec_selected_receiver)*hlagrange
      dux_dz_interpolated = dux_dz_interpolated + dux_dz(i,j,ispec_selected_receiver)*hlagrange
      duz_dx_interpolated = duz_dx_interpolated + duz_dx(i,j,ispec_selected_receiver)*hlagrange
      duz_dz_interpolated = duz_dz_interpolated + duz_dz(i,j,ispec_selected_receiver)*hlagrange

    enddo
  enddo

! compare the computed derivative to the analytical one for validation,
! since we set the vector field above to an analytical value in order to be able to compare to the exact derivative here
  print *,'dUx / dx numerical, exact = ',dux_dx_interpolated, -sin(x_receiver)
  print *,'dUx / dz numerical, exact = ',dux_dz_interpolated, -sin(z_receiver)
  print *,'dUz / dx numerical, exact = ',duz_dx_interpolated, -2*sin(2*x_receiver)
  print *,'dUz / dz numerical, exact = ',duz_dz_interpolated, -2*sin(2*z_receiver)
  print *

  end program compute_gradient

