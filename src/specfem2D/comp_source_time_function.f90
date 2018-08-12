!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%% calls with hdur as an argument are below %%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  double precision function comp_source_time_function_heaviside_hdur(t,hdur)

  implicit none

  double precision :: t,hdur

  double precision, external :: netlib_specfun_erf

  ! compared with calling these same functions below with f0,
  ! one has the relationship hdur = 1 / (PI * f0), or equivalently f0 = 1 / (PI * hdur)

  ! quasi Heaviside, small Gaussian moment-rate tensor with hdur
  comp_source_time_function_heaviside_hdur = 0.5d0*(1.0d0 + netlib_specfun_erf(t/hdur))

  end function comp_source_time_function_heaviside_hdur

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%% calls with f0 as an argument are below %%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  double precision function comp_source_time_function_Gaussian(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0
  double precision :: a

  ! Gaussian wavelet i.e. second integral of a Ricker wavelet
  a = PI**2 * f0**2
  comp_source_time_function_Gaussian = - exp(-a * t**2) / (2.d0 * a)

  end function comp_source_time_function_Gaussian

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_dGaussian(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0
  double precision :: a

  ! first integral of a Ricker wavelet
  a = PI**2 * f0**2
  comp_source_time_function_dGaussian = t * exp(-a * t**2)

  end function comp_source_time_function_dGaussian

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_d2Gaussian(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0
  double precision :: a

  ! Ricker wavelet (second derivative of a Gaussian)
  a = PI**2 * f0**2
  comp_source_time_function_d2Gaussian = (1.d0 - 2.d0 * a * t**2) * exp(-a * t**2)

  end function comp_source_time_function_d2Gaussian

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_d3Gaussian(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0
  double precision :: a

  ! first derivative of a Ricker wavelet (third derivative of a Gaussian)
  a = PI**2 * f0**2
  comp_source_time_function_d3Gaussian = 2.d0 * a * t * exp(-a * t**2) * (- 3.d0 + 2.d0 * a * t**2)

  end function comp_source_time_function_d3Gaussian

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_d4Gaussian(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0

  ! local variables
  double precision :: a

  ! second derivative of a Ricker wavelet
  a = PI**2 * f0**2
  comp_source_time_function_d4Gaussian = - 2.d0 * a * (3.d0 - 12.d0 * a * t*t + 4.d0 * a**2 * t*t*t*t) * exp( -a * t*t )

  end function comp_source_time_function_d4Gaussian

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_Ricker(t,f0)

! Ricker wavelet (second derivative of a Gaussian)

  implicit none

  double precision, intent(in) :: t,f0

  double precision, external :: comp_source_time_function_d2Gaussian

  ! Ricker wavelet
  comp_source_time_function_Ricker = comp_source_time_function_d2Gaussian(t,f0)

  !! another source time function that is improperly called 'Ricker' in some old papers,
  !! e.g., 'Finite-Frequency Kernels Based on Adjoint Methods' by Liu & Tromp, BSSA (2006), is:
  ! comp_source_time_function_Ricker = -2.d0*PI*PI*f0*f0*f0*t * exp(-PI*PI*f0*f0*t*t)

  end function comp_source_time_function_Ricker

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_dRicker(t,f0)

  implicit none

  double precision, intent(in) :: t,f0

  double precision, external :: comp_source_time_function_d3Gaussian

  ! first derivative of a Ricker wavelet
  comp_source_time_function_dRicker = comp_source_time_function_d3Gaussian(t,f0)

  end function comp_source_time_function_dRicker

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_d2Ricker(t,f0)

  implicit none

  double precision, intent(in) :: t,f0

  double precision, external :: comp_source_time_function_d4Gaussian

  ! second derivative of a Ricker wavelet
  comp_source_time_function_d2Ricker = comp_source_time_function_d4Gaussian(t,f0)

  end function comp_source_time_function_d2Ricker

!
!------------------------------------------------------------
!

  double precision function sinc(a)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: a

  if (abs(a) < 1.0d-10) then
    sinc = 1.0d0
  else
    sinc = sin(a)/a
  endif

  end function sinc

!
!------------------------------------------------------------
!

  double precision function cos_taper(a,hdur)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: a,hdur

  double precision :: b

  b = abs(a)
  cos_taper = 0.0

  if (b <= hdur) then
    cos_taper = 1.0;
  else if (b > hdur .and. b < 2.0 * hdur) then
    cos_taper = cos( PI*0.5*(b-hdur)/hdur )
  endif

  end function cos_taper

!
!------------------------------------------------------------
!

  double precision function marmousi_ormsby_wavelet(a)

  use constants, only: PI

  implicit none

  double precision :: sinc

  double precision, intent(in) :: a

  double precision :: f1,f2,f3,f4,b,c,tmp

  ! 5-10-60-80 Hz Ormsby Wavelet for Marmousi2 Model (Gray S. Martin, 2006)
  ! Please find the Ormsby Wavelet here http://subsurfwiki.org/wiki/Ormsby_filter

  f1 = 5.0d0  ! low-cut frequency
  f2 = 10.0d0 ! low-pass frequency
  f3 = 60.0d0 ! high-pass frequency
  f4 = 80.0d0 ! high-cut frequency

  b = sinc( f1 * a )
  tmp = PI * f1
  tmp = tmp * tmp
  tmp = tmp / PI /( f2 - f1 )
  c = b * b * tmp

  b = sinc( f2 * a )
  tmp = PI * f2
  tmp = tmp * tmp
  tmp = tmp / PI /( f2 - f1 )
  c = c - b * b * tmp

  b = sinc( f3 * a )
  tmp = PI * f3
  tmp = tmp * tmp
  tmp = tmp / PI /( f4 - f3 )
  c = c - b * b * tmp

  b = sinc( f4 * a )
  tmp = PI * f4
  tmp = tmp * tmp
  tmp = tmp / PI /( f4 - f3 )
  c = c + b * b * tmp

  marmousi_ormsby_wavelet = c

  end function marmousi_ormsby_wavelet

