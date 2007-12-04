
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour and CNRS, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
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

! compute analytical initial plane wave for Bielak's conditions

  subroutine compute_Bielak_conditions(coord,iglob,npoin,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert)

  implicit none

  include "constants.h"

  integer, intent(in) :: iglob,npoin,it

  double precision, intent(in) :: deltat

  double precision, intent(out) :: dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert

  double precision, dimension(NDIM,npoin), intent(in) :: coord

  double precision :: time_veloc,time_traction,t,x,z

  double precision, external :: ricker_Bielak_veloc

! get the coordinates of the mesh point
      x = coord(1,iglob)
      z = coord(2,iglob)

! times for velocity and traction are staggered i.e. separated by deltat/2.d0
!! DK DK YYYYYYYYYYYYYYYY UGLY il faudra peut-etre ajouter ou enlever deltat ici
      time_veloc = (it-1)*deltat + deltat/2.d0 + time_offset
      time_traction = time_veloc + deltat/2.d0

      t = time_traction

! derivatives of analytical expression of horizontal and vertical displacements,
! computed using the "Mathematica" script in UTILS/deriv_ricker_spatial.m
  dxUx = (sqrt(3.d0)*a*((-8*t + 4*x)*exp(-a*(t - x/2.d0)**2) + &
      ((2*t - x)*(-2 + a*(-2*t + x)**2))*exp(-a*(t - x/2.d0)**2) + &
      (2*(-2*t + x - sqrt(3.d0)*(-9 + z)))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
      ((1 - (a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/2.d0)* &
         (-2*t + x - sqrt(3.d0)*(-9 + z)))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
      (2*(-2*t + x + sqrt(3.d0)*(-9 + z)))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
      ((1 - (a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/2.d0)* &
         (-2*t + x + sqrt(3.d0)*(-9 + z)))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0)))/4.d0

  dzUx = (3*a*(((t + (-x + sqrt(3.d0)*(-9 + z))/2.d0)* &
         (1 - (a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/2.d0))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) - &
      ((1 - (a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/2.d0)* &
         (t - x/2.d0 - (sqrt(3.d0)*(-9 + z))/2.d0))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
      (2*t - x + sqrt(3.d0)*(-9 + z))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
      (-2*t + x + sqrt(3.d0)*(-9 + z))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0)))/2.d0

  dxUz = (a*((2*t - x - sqrt(3.d0)*(-9 + z))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
      (-2*t + x - sqrt(3.d0)*(-9 + z))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
      ((1 - (a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/2.d0)* &
         (-2*t + x - sqrt(3.d0)*(-9 + z)))/2.d0*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) - &
      ((1 - (a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/2.d0)* &
         (-2*t + x + sqrt(3.d0)*(-9 + z)))/2.d0*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0)))/2.d0

  dzUz = (sqrt(3.d0)*a*(((t + (-x + sqrt(3.d0)*(-9 + z))/2.d0)* &
         (1 - (a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/2.d0))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
      (2*t - x - sqrt(3.d0)*(-9 + z))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
      ((1 - (a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/2.d0)* &
         (t - x/2.d0 - (sqrt(3.d0)*(-9 + z))/2.d0))*exp(-(a*(-2*t + x + sqrt(3.d0)*(-9 + z))**2)/4.d0) + &
      (2*t - x + sqrt(3.d0)*(-9 + z))*exp(-(a*(2*t - x + sqrt(3.d0)*(-9 + z))**2)/4.d0)))/2.d0

      t = time_veloc

! analytical expression of the two components of the velocity vector
      veloc_horiz = rac3sur2 * ricker_Bielak_veloc(t - x/2.d0 + (9 - z) * rac3sur2) &
        + rac3sur2 * ricker_Bielak_veloc(t - x/2.d0 - (9 - z) * rac3sur2) &
        + rac3 * ricker_Bielak_veloc(t - x/2.d0)
      veloc_vert = - HALF * ricker_Bielak_veloc(t - x/2.d0 + (9 - z) * rac3sur2) &
        + HALF * ricker_Bielak_veloc(t - x/2.d0 - (9 - z) * rac3sur2)

  end subroutine compute_Bielak_conditions

! ********

! compute time variation of the source for analytical initial plane wave
  double precision function ricker_Bielak_displ(t)

  implicit none

  include "constants.h"

  double precision :: t

! Ricker
  ricker_Bielak_displ = - (1 - 2*a*t**2)*exp(-a*t**2)

  end function ricker_Bielak_displ

! *******

! compute time variation of the source for analytical initial plane wave
  double precision function ricker_Bielak_veloc(t)

  implicit none

  include "constants.h"

  double precision :: t

! first time derivative of a Ricker
  ricker_Bielak_veloc = 2*a*t*(3 - 2*a*t**2)*exp(-a*t**2)

  end function ricker_Bielak_veloc

! *******

! compute time variation of the source for analytical initial plane wave
  double precision function ricker_Bielak_accel(t)

  implicit none

  include "constants.h"

  double precision :: t

! second time derivative of a Ricker
  ricker_Bielak_accel = 2*a*(3 - 12*a*t**2 + 4*a**2*t**4)* exp(-a*t**2)

  end function ricker_Bielak_accel

