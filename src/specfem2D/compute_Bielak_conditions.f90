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

! compute analytical initial plane wave for Bielak's conditions

  subroutine compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                                       x0_source, z0_source, A_plane, B_plane, C_plane, anglesource, anglesource_refl, &
                                       c_inc, c_refl, time_offset,f0)

  use constants, only: NDIM

  implicit none

  integer, intent(in) :: iglob,nglob,it

  double precision, dimension(NDIM,nglob), intent(in) :: coord
  double precision, intent(in) :: deltat

  double precision, intent(out) :: dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert
  double precision, intent(in) :: x0_source, z0_source

  double precision, dimension(NDIM),intent(in) :: A_plane, B_plane, C_plane
  double precision, intent(in) :: anglesource, anglesource_refl

  double precision :: c_inc, c_refl, time_offset, f0

  ! local parameters
  double precision :: time_veloc,time_traction,t,x,z
  double precision, external :: ricker_Bielak_veloc

! get the coordinates of the mesh point
  x = coord(1,iglob) - x0_source
  z = z0_source - coord(2,iglob)

! times for velocity and traction are staggered i.e. separated by deltat/2.d0
  time_veloc = (it-1)*deltat + deltat/2.d0 + time_offset
  time_traction = time_veloc + deltat/2.d0

  t = time_traction

!!$!SV30
!!$
!!$!analytical expression of the displacement for a SV 30 degrees and 0.3333 poisson ratio
!!$!  Ux = sqrt(3.d0)/2.d0 * rickertest(t - x/2.d0 + (9 - z) * sqrt(3.d0)/2.d0) &
!!$!     + sqrt(3.d0)/2.d0 * rickertest(t - x/2.d0 - (9 - z) * sqrt(3.d0)/2.d0) &
!!$!     + sqrt(3.d0) * rickertest(t - x/2.d0)
!!$!  Uz = - HALF * rickertest(t - x/2.d0 + (9 - z) * sqrt(3.d0)/2.d0) &
!!$!       + HALF * rickertest(t - x/2.d0 - (9 - z) * sqrt(3.d0)/2.d0)
!!$
!!$
!!$! derivatives of analytical expression of horizontal and vertical displacements,
!!$! computed using the "Mathematica" script in UTILS/deriv_ricker_spatial.m
!!$  dxUx = (sqrt(3.d0)*a*((-8*t + 4*x)*exp(-a*(t - x/2.d0)**2) + &
!!$      ((2*t - x)*(-2 + a*(-2*t + x)**2))*exp(-a*(t - x/2.d0)**2) + &
!!$      (2*(-2*t + x - sqrt(3.d0)*(-9 + z)))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
!!$      ((1 - (a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/2.d0)* &
!!$         (-2*t + x - sqrt(3.d0)*(-9 + z)))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
!!$      (2*(-2*t + x + sqrt(3.d0)*(-9 + z)))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
!!$      ((1 - (a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/2.d0)* &
!!$         (-2*t + x + sqrt(3.d0)*(-9 + z)))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0)))/4.d0
!!$
!!$  dzUx = (3*a*(((t + (-x + sqrt(3.d0)*(-9 + z))/2.d0)* &
!!$         (1 - (a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/2.d0))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) - &
!!$      ((1 - (a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/2.d0)* &
!!$         (t - x/2.d0 - (sqrt(3.d0)*(-9 + z))/2.d0))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
!!$      (2*t - x + sqrt(3.d0)*(-9 + z))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
!!$      (-2*t + x + sqrt(3.d0)*(-9 + z))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0)))/2.d0
!!$
!!$  dxUz = (a*((2*t - x - sqrt(3.d0)*(-9 + z))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
!!$      (-2*t + x - sqrt(3.d0)*(-9 + z))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
!!$      ((1 - (a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/2.d0)* &
!!$         (-2*t + x - sqrt(3.d0)*(-9 + z)))/2.d0*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) - &
!!$      ((1 - (a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/2.d0)* &
!!$         (-2*t + x + sqrt(3.d0)*(-9 + z)))/2.d0*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0)))/2.d0
!!$
!!$  dzUz = (sqrt(3.d0)*a*(((t + (-x + sqrt(3.d0)*(-9 + z))/2.d0)* &
!!$         (1 - (a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/2.d0))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
!!$      (2*t - x - sqrt(3.d0)*(-9 + z))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
!!$      ((1 - (a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/2.d0)* &
!!$         (t - x/2.d0 - (sqrt(3.d0)*(-9 + z))/2.d0))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
!!$      (2*t - x + sqrt(3.d0)*(-9 + z))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0)))/2.d0

! to compute the derivative of the displacement, we take the velocity ricker expression and we multiply by
! the derivative of the interior argument of ricker_Bielak_veloc

  dxUx = A_plane(1) * ricker_Bielak_veloc(t - sin(anglesource)*x/c_inc + cos(anglesource)*z/c_inc,f0) * (-sin(anglesource)/c_inc)&
       + B_plane(1) * ricker_Bielak_veloc(t - sin(anglesource)*x/c_inc - cos(anglesource)*z/c_inc,f0) * (-sin(anglesource)/c_inc)&
       + C_plane(1) * ricker_Bielak_veloc(t - sin(anglesource_refl)*x/c_refl - cos(anglesource_refl)*z/c_refl,f0)&
       * (-sin(anglesource_refl)/c_refl)

  dzUx = A_plane(1) * ricker_Bielak_veloc(t - sin(anglesource)*x/c_inc + cos(anglesource)*z/c_inc,f0) * (-cos(anglesource)/c_inc)&
       + B_plane(1) * ricker_Bielak_veloc(t - sin(anglesource)*x/c_inc - cos(anglesource)*z/c_inc,f0) * (cos(anglesource)/c_inc)&
       + C_plane(1) * ricker_Bielak_veloc(t - sin(anglesource_refl)*x/c_refl - cos(anglesource_refl)*z/c_refl,f0)&
       * (cos(anglesource_refl)/c_refl)

  dxUz = A_plane(2) * ricker_Bielak_veloc(t - sin(anglesource)*x/c_inc + cos(anglesource)*z/c_inc,f0) * (-sin(anglesource)/c_inc)&
       + B_plane(2) * ricker_Bielak_veloc(t - sin(anglesource)*x/c_inc - cos(anglesource)*z/c_inc,f0) * (-sin(anglesource)/c_inc)&
       + C_plane(2) * ricker_Bielak_veloc(t - sin(anglesource_refl)*x/c_refl - cos(anglesource_refl)*z/c_refl,f0)&
       * (-sin(anglesource_refl)/c_refl)

  dzUz = A_plane(2) * ricker_Bielak_veloc(t - sin(anglesource)*x/c_inc + cos(anglesource)*z/c_inc,f0) * (-cos(anglesource)/c_inc)&
       + B_plane(2) * ricker_Bielak_veloc(t - sin(anglesource)*x/c_inc - cos(anglesource)*z/c_inc,f0) * (cos(anglesource)/c_inc)&
       + C_plane(2) * ricker_Bielak_veloc(t - sin(anglesource_refl)*x/c_refl - cos(anglesource_refl)*z/c_refl,f0)&
       * (cos(anglesource_refl)/c_refl)

  t = time_veloc

!!$!SV30
!!$! analytical expression of the two components of the velocity vector
!!$      veloc_horiz = (sqrt(3.d0)/2.d0) * ricker_Bielak_veloc(t - x/2.d0 + (9 - z) * (sqrt(3.d0)/2.d0)) &
!!$        + (sqrt(3.d0)/2.d0) * ricker_Bielak_veloc(t - x/2.d0 - (9 - z) * (sqrt(3.d0)/2.d0)) &
!!$        + sqrt(3.d0) * ricker_Bielak_veloc(t - x/2.d0)
!!$      veloc_vert = - HALF * ricker_Bielak_veloc(t - x/2.d0 + (9 - z) * (sqrt(3.d0)/2.d0)) &
!!$        + HALF * ricker_Bielak_veloc(t - x/2.d0 - (9 - z) * (sqrt(3.d0)/2.d0))

  veloc_horiz = A_plane(1) * ricker_Bielak_veloc(t - sin(anglesource)*x/c_inc + cos(anglesource)*z/c_inc,f0) &
       + B_plane(1) * ricker_Bielak_veloc(t - sin(anglesource)*x/c_inc - cos(anglesource)*z/c_inc,f0) &
       + C_plane(1) * ricker_Bielak_veloc(t - sin(anglesource_refl)*x/c_refl - cos(anglesource_refl)*z/c_refl,f0)
  veloc_vert = A_plane(2) * ricker_Bielak_veloc(t - sin(anglesource)*x/c_inc + cos(anglesource)*z/c_inc,f0) &
       + B_plane(2) * ricker_Bielak_veloc(t - sin(anglesource)*x/c_inc - cos(anglesource)*z/c_inc,f0) &
       + C_plane(2) * ricker_Bielak_veloc(t - sin(anglesource_refl)*x/c_refl - cos(anglesource_refl)*z/c_refl,f0)

  end subroutine compute_Bielak_conditions

! ********

! compute time variation of the source for analytical initial plane wave
  double precision function ricker_Bielak_integrale_displ(t,f0)

  use constants, only: PI

  implicit none

  double precision :: t,f0,a

  a = pi*pi*f0*f0

! integral of a Ricker, i.e. first derivative of a Gaussian
  ricker_Bielak_integrale_displ = t*exp(-a*t**2)

  end function ricker_Bielak_integrale_displ

! ********

! compute time variation of the source for analytical initial plane wave
  double precision function ricker_Bielak_displ(t,f0)

  use constants, only: PI

  implicit none

  double precision :: t,f0,a

  a = pi*pi*f0*f0

! Ricker, i.e. second derivative of a Gaussian
  ricker_Bielak_displ = (1 - 2*a*t**2)*exp(-a*t**2)

  end function ricker_Bielak_displ

! *******

! compute time variation of the source for analytical initial plane wave
  double precision function ricker_Bielak_veloc(t,f0)

  use constants, only: PI

  implicit none

  double precision :: t,f0,a

  a = pi*pi*f0*f0

! first time derivative of a Ricker, i.e. third derivative of a Gaussian
  ricker_Bielak_veloc = - 2*a*t*(3 - 2*a*t**2)*exp(-a*t**2)

  end function ricker_Bielak_veloc

! *******

! compute time variation of the source for analytical initial plane wave
  double precision function ricker_Bielak_accel(t,f0)

  use constants, only: PI

  implicit none

  double precision :: t,f0,a

  a = pi*pi*f0*f0

! second time derivative of a Ricker, i.e. fourth derivative of a Gaussian
  ricker_Bielak_accel = - 2*a*(3 - 12*a*t**2 + 4*a**2*t**4)* exp(-a*t**2)

  end function ricker_Bielak_accel

! *******

! compute time variation of the source for analytical initial plane wave
  double precision function ricker_Bielak_third_derivative(t,f0)

  use constants, only: PI

  implicit none

  double precision :: t,f0,a

  a = pi*pi*f0*f0

! third time derivative of a Ricker, i.e. fifth derivative of a Gaussian
  ricker_Bielak_third_derivative = 4*a**2*t*exp(-a*t**2)*(4*a**2*t**4 - 20*a*t**2 + 15)

  end function ricker_Bielak_third_derivative

! *******

! compute time variation of the source for analytical initial plane wave
  double precision function ricker_Bielak_fourth_derivative(t,f0)

  use constants, only: PI

  implicit none

  double precision :: t,f0,a

  a = pi*pi*f0*f0

! fourth time derivative of a Ricker, i.e. sixth derivative of a Gaussian
  ricker_Bielak_fourth_derivative = -4*a**2*exp(-a*t**2)*(8*a**3*t**6 - 60*a**2*t**4 + 90*a*t**2 - 15)

  end function ricker_Bielak_fourth_derivative

! *******

! compute time variation of the source for analytical initial plane wave
  double precision function ricker_Bielak_fifth_derivative(t,f0)

  use constants, only: PI

  implicit none

  double precision :: t,f0,a

  a = pi*pi*f0*f0

! fifth time derivative of a Ricker, i.e. seventh derivative of a Gaussian
  ricker_Bielak_fifth_derivative = 8*a**3*t*exp(-a*t**2)*(8*a**3*t**6 - 84*a**2*t**4 + 210*a*t**2 - 105)

  end function ricker_Bielak_fifth_derivative

! *******

! compute time variation of the source for analytical initial plane wave
  double precision function ricker_Bielak_sixth_derivative(t,f0)

  use constants, only: PI

  implicit none

  double precision :: t,f0,a

  a = pi*pi*f0*f0

! sixth time derivative of a Ricker, i.e. eighth derivative of a Gaussian
  ricker_Bielak_sixth_derivative = -8*a**3*exp(-a*t**2)*(16*a**4*t**8 - 224*a**3*t**6 + 840*a**2*t**4 - 840*a*t**2 + 105)

  end function ricker_Bielak_sixth_derivative

