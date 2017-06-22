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

  double precision function comp_source_time_function_heavi(t,hdur)

  implicit none

  double precision :: t,hdur

  double precision, external :: netlib_specfun_erf

  ! quasi Heaviside, small Gaussian moment-rate tensor with hdur
  comp_source_time_function_heavi = 0.5d0*(1.0d0 + netlib_specfun_erf(t/hdur))

  end function comp_source_time_function_heavi


!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_gauss(t,hdur)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,hdur
  double precision :: hdur_decay,a

  ! note: hdur given is hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
  !           and SOURCE_DECAY_MIMIC_TRIANGLE ~ 1.68
  hdur_decay = hdur

  ! this here uses a stronger Gaussian decay rate (empirical value) to avoid non-zero onset times;
  ! however, it should mimik a triangle source time function...
  !hdur_decay = hdur  / SOURCE_DECAY_STRONG

  ! note: a nonzero time to start the simulation with would lead to more high-frequency noise
  !          due to the (spatial) discretization of the point source on the mesh

  ! Gaussian wavelet
  a = 1.d0 / (hdur_decay**2)
  comp_source_time_function_gauss = exp(-a * t**2) / (sqrt(PI)*hdur_decay)

  end function comp_source_time_function_gauss

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_dgau(t,hdur)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,hdur
  double precision :: hdur_decay,a

  ! note: hdur given is hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
  !           and SOURCE_DECAY_MIMIC_TRIANGLE ~ 1.68
  hdur_decay = hdur

  ! this here uses a stronger Gaussian decay rate (empirical value) to avoid non-zero onset times;
  ! however, it should mimik a triangle source time function...
  !hdur_decay = hdur  / SOURCE_DECAY_STRONG

  ! note: a nonzero time to start the simulation with would lead to more high-frequency noise
  !          due to the (spatial) discretization of the point source on the mesh

  ! first derivative of a Gaussian wavelet
  a = 1.d0 / (hdur_decay**2)
  comp_source_time_function_dgau = - 2.d0 * a * t * exp(-a * t**2) / (sqrt(PI)*hdur_decay)

  end function comp_source_time_function_dgau

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_d2gau(t,hdur)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,hdur
  double precision :: hdur_decay,a

  ! note: hdur given is hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
  !           and SOURCE_DECAY_MIMIC_TRIANGLE ~ 1.68
  hdur_decay = hdur

  ! this here uses a stronger Gaussian decay rate (empirical value) to avoid non-zero onset times;
  ! however, it should mimik a triangle source time function...
  !hdur_decay = hdur  / SOURCE_DECAY_STRONG

  ! note: a nonzero time to start the simulation with would lead to more high-frequency noise
  !          due to the (spatial) discretization of the point source on the mesh

  ! second derivative of a Gaussian wavelet
  a = 1.d0 / (hdur_decay**2)
  comp_source_time_function_d2gau = 2.d0 * a * (-1.d0 + 2.d0 * a * t**2) * exp(-a * t**2) / (sqrt(PI)*hdur_decay)

  end function comp_source_time_function_d2gau


!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_gaussB(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0
  double precision :: a

  ! Gaussian wavelet
  a = PI**2 * f0**2

  ! modified Gaussian
  comp_source_time_function_gaussB = exp(-a * t**2)

  end function comp_source_time_function_gaussB

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_dgaussB(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0
  double precision :: a

  ! Gaussian wavelet
  a = PI**2 * f0**2

  ! modified Gaussian
  comp_source_time_function_dgaussB = - 2.d0 * a * t * exp(-a * t**2)

  end function comp_source_time_function_dgaussB

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_d2gaussB(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0
  double precision :: a

  ! Gaussian wavelet
  a = PI**2 * f0**2

  ! modified Gaussian
  comp_source_time_function_d2gaussB = 2.0d0 * a * (2.0d0 * a * t**2 - 1.0d0) * exp(-a * t**2)

  end function comp_source_time_function_d2gaussB


!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_d3gaussB(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0
  double precision :: a

  ! Gaussian wavelet
  a = PI**2 * f0**2

  ! modified Gaussian
  comp_source_time_function_d3gaussB = 4.0d0 * a**2 * t * (3.0d0 - 2.d0 * a * t**2) * exp(-a * t**2)

  end function comp_source_time_function_d3gaussB



!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_rickr(t,f0)

!  Ricker (second derivative of a Gaussian)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0

  ! local variables
  double precision :: a

  ! Ricker wavelet
  a = PI**2 * f0**2
  comp_source_time_function_rickr = (1.d0 - 2.d0 * a * t*t) * exp( -a * t*t )

  !!! another source time function they have called 'ricker' in some old papers,
  !!! e.g., 'Finite-Frequency Kernels Based on Adjoint Methods' by Liu & Tromp, BSSA (2006)
  !!! in order to benchmark those simulations, the following formula is needed.
  ! comp_source_time_function_rickr = -2.d0*PI*PI*f0*f0*f0*t * exp(-PI*PI*f0*f0*t*t)

  end function comp_source_time_function_rickr

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_drck(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0

  ! local variables
  double precision :: a

  ! first derivative of a Ricker wavelet
  a = PI**2 * f0**2
  comp_source_time_function_drck = 2.d0 * a * t * (-3.d0 + 2.d0 * a * t*t) * exp( -a * t*t )

  end function comp_source_time_function_drck

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_d2rck(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0

  ! local variables
  double precision :: a

  ! second derivative of a Ricker wavelet
  a = PI**2 * f0**2
  comp_source_time_function_d2rck = - 2.d0 * a * (3.d0 - 12.d0 * a * t*t + 4.d0 * a**2 * t*t*t*t) * exp( -a * t*t )

  end function comp_source_time_function_d2rck

!-------------------------------------------------------------------------------------------------

  double precision function sinc(a)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: a

  ! sinc function defined here

  if (abs(a) < 1.0d-10) then
    sinc = 1.0d0
  else
    sinc = sin(a)/a
  endif

  end function sinc

!-------------------------------------------------------------------------------------------------

  double precision function cos_taper(a,hdur)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: a,hdur

  double precision :: b

  ! cos_taper function defined here

  b = abs(a)
  cos_taper = 0.0

  if (b <= hdur) then
    cos_taper = 1.0;
  else if (b > hdur .and. b < 2.0 * hdur) then
    cos_taper = cos( PI*0.5*(b-hdur)/hdur )
  endif

  end function cos_taper

!-----------------------------------------------------------------------------------------------------

  double precision function marmousi_ormsby_wavelet(a)

  use constants, only: PI

  implicit none

  double precision :: sinc

  double precision, intent(in) :: a

  double precision :: f1,f2,f3,f4,b,c,tmp

  ! 5-10-60-80 Hz Ormsby Wavelet for Marmousi2 Model (Gray S. Martin, 2006)
  ! Please find the Ormsby Wavelet here http://subsurfwiki.org/wiki/Ormsby_filter

  f1 = 5.0d0;  ! low-cut frequency
  f2 = 10.0d0; ! low-pass frequency
  f3 = 60.0d0; ! high-pass frequency
  f4 = 80.0d0; ! high-cut frequency

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

!-------------------------------------------------------------------------------------------------


