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

  double precision, intent(in) :: t,hdur

  double precision, external :: netlib_specfun_erf

  ! compared with calling these same functions below with f0,
  ! one has the relationship hdur = 1 / (PI * f0), or equivalently f0 = 1 / (PI * hdur)

  ! quasi Heaviside, small Gaussian moment-rate tensor with hdur
  comp_source_time_function_heaviside_hdur = 0.5d0 * (1.0d0 + netlib_specfun_erf(t/hdur))

  end function comp_source_time_function_heaviside_hdur

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_Gaussian_norm(t,hdur)

! normalized Gaussian source time function
!
! note: here, "normalized" means that the integral value of the source time function becomes 1.
!       this function is also used as default in the SPECFEM3D versions.
!
!       this Gaussian function is using a different normalization factor compared to the comp_source_time_function_Gaussian()
!       below that implements a Gaussian formulation derived from the second integral of the Ricker wavelet
!       (given by comp_source_time_function_Ricker() or comp_source_time_function_d2Gaussian()).
!
!       thus, amplitudes and hdur are different here as compared to comp_source_time_function_Gaussian() values.

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,hdur
  double precision :: hdur_decay,a

  ! note: hdur given is hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
  !           and SOURCE_DECAY_MIMIC_TRIANGLE ~ 1.68
  hdur_decay = hdur

  ! Gaussian wavelet
  a = 1.d0 / (hdur_decay**2)
  comp_source_time_function_Gaussian_norm = exp(-a * t**2) / (sqrt(PI) * hdur_decay)

  end function comp_source_time_function_Gaussian_norm

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_d2Gaussian_norm(t,hdur)

! second derivative of the normalized Gaussian function above

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,hdur
  double precision :: hdur_decay,a

  ! note: hdur given is hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
  !           and SOURCE_DECAY_MIMIC_TRIANGLE ~ 1.68
  hdur_decay = hdur

  ! second derivative of a Gaussian wavelet
  a = 1.d0 / (hdur_decay**2)
  comp_source_time_function_d2Gaussian_norm = 2.d0 * a * (-1.d0 + 2.d0 * a * t**2) * exp(-a * t**2) / (sqrt(PI) * hdur_decay)

  end function comp_source_time_function_d2Gaussian_norm



! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%% calls with f0 as an argument are below %%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  double precision function comp_source_time_function_Gaussian(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0
  double precision :: a

  ! Gaussian wavelet i.e. second integral of a Ricker wavelet
  ! note: The Gaussian here is defined as the second integral of the Ricker wavelet defined below by
  !       routine comp_source_time_function_d2Gaussian(t,f0).
  !       Integrating that Ricker function twice will lead to a minus sign in front of the exponential term here.
  !       Thus, the Gaussian will be inverted, going from zero to
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

  ! first integral of a Ricker wavelet (or equivalent first derivative of a Gaussian)
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
  comp_source_time_function_d3Gaussian = 2.d0 * a * t * (- 3.d0 + 2.d0 * a * t**2) * exp(-a * t**2)

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
  comp_source_time_function_d4Gaussian = - 2.d0 * a * (3.d0 - 12.d0 * a * t*t + 4.d0 * a**2 * t*t*t*t) * exp(-a * t**2)

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


!------------------------------------------------------------
!
! specialized source time functions
!
!------------------------------------------------------------



  double precision function comp_source_time_function_ocean_I(t,f0)

! ocean acoustics type I

  use constants, only: PI,QUARTER,HALF,ONE,TWO
  implicit none

  double precision, intent(in) :: t,f0
  double precision :: Tc,omega_coa,omegat

  Tc = 4.d0 / f0

  omega_coa = TWO * PI * f0
  omegat = omega_coa * t

  if (t > 0.d0 .and. t < Tc) then
    ! source time function from Computational Ocean Acoustics
    comp_source_time_function_ocean_I = HALF * sin(omegat) * (ONE - cos(QUARTER * omegat))

    ! alternative
    !comp_source_time_function_ocean_I = HALF / omega_coa / omega_coa * &
    !      ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) - 8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) )
  else
    comp_source_time_function_ocean_I = 0.d0
  endif

  end function comp_source_time_function_ocean_I

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_ocean_II(t,f0)

! ocean acoustics type II

  use constants, only: PI,ZERO,QUARTER,HALF,ONE,TWO
  implicit none

  double precision, intent(in) :: t,f0
  double precision :: Tc,omega_coa,omegat

  Tc = 4.d0 / f0

  omega_coa = TWO * PI * f0
  omegat = omega_coa * t

  if (t > 0.d0 .and. t < Tc) then
    ! source time function from Computational Ocean Acoustics
    comp_source_time_function_ocean_II = HALF / omega_coa / omega_coa * ( - sin(omegat) + 8.d0 / 9.d0 * sin(3.d0 / 4.d0 * omegat) &
                                           + 8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) - 1.d0 / 15.d0 * omegat )

    ! alternative
    !comp_source_time_function_ocean_II = HALF / omega_coa / omega_coa * ( sin(omegat) - 8.d0 / 9.d0 * sin(3.d0/ 4.d0 * omegat) &
    !                       - 8.d0 / 25.d0 * sin(5.d0 / 4.d0 * omegat) - 1.d0 / 15.d0 * t + 1.d0 / 15.d0 * Tc )
  else if (t > 0.d0) then
    comp_source_time_function_ocean_II = - HALF / omega_coa / 15.d0 * Tc
  else
    comp_source_time_function_ocean_II = ZERO
  endif

  ! alternative - For source 1 OASES
  !Tc = 1.d0 / f0
  !
  !if (t > 0.d0 .and. t < Tc) then
  !  comp_source_time_function_ocean_II = 0.75d0 - cos(omegat) + 0.25d0*cos(TWO * omegat)
  !else
  !  comp_source_time_function_ocean_II = ZERO
  !endif

  end function comp_source_time_function_ocean_II


!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_d2ocean_II(t,f0)

! second derivative of ocean acoustics type II

  use constants, only: PI,ZERO,QUARTER,HALF,ONE,TWO
  implicit none

  double precision, intent(in) :: t,f0
  double precision :: Tc,omega_coa,omegat

  Tc = 4.d0 / f0

  omega_coa = TWO * PI * f0
  omegat = omega_coa * t

  if (t > 0.d0 .and. t < Tc) then
    ! second derivate of source time function from Computational Ocean Acoustics
    comp_source_time_function_d2ocean_II = HALF * (ONE - cos(omegat / 4.0d0)) * sin(omegat)
  else
    comp_source_time_function_d2ocean_II = ZERO
  endif

  end function comp_source_time_function_d2ocean_II

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_burst(t,f0,burst_band_width)

! burst type source time function

  use constants, only: PI,ZERO,ONE,TWO

  implicit none

  double precision,intent(in) :: t,f0,burst_band_width
  double precision :: Tc,Nc,omega_coa,omegat

  Nc = TWO * f0 / burst_band_width
  Tc = Nc / f0

  omega_coa = TWO * PI * f0
  omegat = omega_coa * t

  if (t > 0.d0 .and. t < Tc) then
    comp_source_time_function_burst = 0.5d0 * (ONE - cos(omegat / Nc)) * sin(omegat)
  else
    comp_source_time_function_burst = ZERO
  endif

  end function comp_source_time_function_burst

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_d2burst(t,f0,burst_band_width)

! second derivative of burst type source time function

  use constants, only: PI,ZERO,ONE,TWO

  implicit none

  double precision,intent(in) :: t,f0,burst_band_width
  double precision :: Tc,Nc,omega_coa,omegat

  Nc = TWO * f0 / burst_band_width
  Tc = Nc / f0

  omega_coa = TWO * PI * f0
  omegat = omega_coa * t

  if (t > 0.d0 .and. t < Tc) then ! t_used > 0 t_used < Nc/f0) then
    comp_source_time_function_d2burst = 0.5d0 * (omega_coa)**2 * sin(omegat) * cos(omegat / Nc) / Nc**2 &
                                         - 0.5d0 * (omega_coa)**2 * sin(omegat) * (ONE - cos(omegat / Nc)) &
                                         + (omega_coa)**2 * cos(omegat) * sin(omegat / Nc) / Nc
  else
    comp_source_time_function_d2burst = ZERO
  endif

  ! alternative - Integral of burst
  !if (t > 0.d0 .and. t < Tc) then
  !  comp_source_time_function_d2burst = - ( Nc*( (Nc+1.0d0)*cos((omega_coa*(Nc-1.0d0)*t_used)/Nc) + &
  !                                          (Nc-1.0d0)*cos((omega_coa*(Nc+1.0d0)*t_used)/Nc)) - &
  !                                           TWO*(Nc**2-1.0d0)*cos(omega_coa*t_used) ) / (8.0d0*PI*f0*(Nc-1)*(Nc+1))
  !else
  !  comp_source_time_function_d2burst = ZERO
  !endif

  ! Double integral of burst
  !if (t > 0.d0 .and. t < Tc) then
  !  comp_source_time_function_d2burst = - ( -sin(TWO*f0*Pi*t_used)/(8.0d0*f0**TWO*Pi**2) + &
  !                                          (Nc**2*sin((TWO*f0*(Nc-1)*PI*t_used)/Nc))/(16.0d0*f0**2*(Nc-1)**2*Pi**2) + &
  !                                          (Nc**2*sin((TWO*f0*(Nc+1)*PI*t_used)/Nc))/(16.0d0*f0**2*(Nc+1)**2*Pi**2) )
  !else
  !  comp_source_time_function_d2burst = ZERO
  !endif

  end function comp_source_time_function_d2burst

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
  cos_taper = 0.d0

  if (b <= hdur) then
    cos_taper = 1.d0
  else if (b > hdur .and. b < 2.d0 * hdur) then
    cos_taper = cos(PI * 0.5d0 * (b-hdur)/hdur)
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


!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_marmousi(t,f0)

! source time function with Marmousi Ormsby wavelet

  use constants, only: PI

  implicit none

  double precision,intent(in) :: t,f0
  double precision :: hdur
  double precision, external :: cos_taper, marmousi_ormsby_wavelet

  hdur = 1.0 / f0  ! hdur = 1.0 / 35.0

  comp_source_time_function_marmousi = cos_taper(t,hdur) * marmousi_ormsby_wavelet(PI * t)

  end function comp_source_time_function_marmousi
!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_mono(t,f0)

  ! monochromatic source time function

  use constants, only: PI,TAPER_MONOCHROMATIC_SOURCE

  implicit none

  double precision,intent(in) :: t,f0
  double precision :: tt,omega_coa,taper_val
  integer :: taper

  omega_coa = 2.d0 * PI * f0
  tt = omega_coa * t

  taper = ceiling(TAPER_MONOCHROMATIC_SOURCE * f0)
  if (t < taper / f0) then
    taper_val = 0.5d0 - 0.5d0 * cos(tt / taper / 2.d0)
  else
    taper_val = 1.d0
  endif

  comp_source_time_function_mono = sin(tt) * taper_val

  end function comp_source_time_function_mono

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_d2mono(t,f0)

  ! second derivative of monochromatic source time function

  use constants, only: PI,TAPER_MONOCHROMATIC_SOURCE

  implicit none

  double precision,intent(in) :: t,f0
  double precision :: tt,omega_coa,taper_val
  integer :: taper

  omega_coa = 2.d0 * PI * f0
  tt = omega_coa * t

  taper = ceiling(TAPER_MONOCHROMATIC_SOURCE * f0)
  if (t < taper / f0) then
    taper_val = 0.5d0 - 0.5d0 * cos(tt / taper / 2.d0)
  else
    taper_val = 1.d0
  endif

  comp_source_time_function_d2mono = - omega_coa*omega_coa * sin(tt) * taper_val

  end function comp_source_time_function_d2mono

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_ext(it_index,isource)

! reads from external source time function file

  use specfem_par, only: myrank, NSTEP, name_of_source_file

  implicit none

  integer,intent(in) :: it_index,isource

  double precision :: dummy_t,stf_val
  integer :: file_unit,ier
  character(len=250) :: error_msg
  character(len=150), parameter :: error_msg1 = 'Error opening the file that contains the external source: '
  character(len=256) :: line

  ! file unit for external source time function files
  file_unit = 800 + isource

  ! opens external file to read in source time function
  if (it_index == 1) then
    ! reads in from external source time function file
    open(unit=file_unit,file=trim(name_of_source_file(isource)),status='old',action='read',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening source time function file: ',trim(name_of_source_file(isource))
      error_msg = trim(error_msg1)//trim(name_of_source_file(isource))
      call exit_MPI(myrank,error_msg)
    endif
  endif

  ! reads in 2-column file values (time value in first column will be ignored)
  ! format: #time #stf-value
  !read(file_unit,*,iostat=ier) dummy_t, stf_val
  !
  ! reads in either 2-column or single-column file values
  ! 2-column format     : #time #stf-value
  ! single-column format: #stf-value
  !
  ! reads line value(s) as a string
  read(file_unit,'(A)',iostat=ier) line
  if (ier == 0) then
    ! suppress leading white spaces, if any
    line = adjustl(line)

    ! check if file contains single-column or 2-columns values
    ! try first to read in 2-column line
    read(line,*,iostat=ier) dummy_t, stf_val
    if (ier /= 0) then
      ! try single-column line
      read(line,*,iostat=ier) stf_val
    endif
  endif

  ! check if reading okay
  if (ier /= 0) then
    print *,'Error reading source time function file: ',trim(name_of_source_file(isource)),' at line ',it_index
    print *,'Please make sure the file contains the same number of lines as the number of timesteps NSTEP ',NSTEP
    call exit_MPI(myrank,'Error reading source time function file')
  endif

  ! closes external file
  if (it_index == NSTEP) close(file_unit)

  comp_source_time_function_ext = stf_val

  end function comp_source_time_function_ext

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_Brune(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0

  ! local variables
  double precision :: omega,omegat,stf_val

  ! Brune source-time function
  ! Moment function
  !if (t < 0.d0) then
  !  stf_val = 0.d0
  !else
  !  omegat = 2.d0*PI*f0*t
  !  stf_val = 1.d0 - exp( -omegat ) * (1.0d0+omegat)
  !endif

  ! Moment rate function
  omega = 2.d0 * PI * f0
  omegat = omega * t

  if (t < 0.d0) then
    stf_val = 0.d0
  else
    stf_val = omega * omegat * exp(-omegat)
  endif

  comp_source_time_function_Brune = stf_val

  end function comp_source_time_function_Brune

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_Smooth_Brune(t,f0)

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0

  ! local variables
  double precision,parameter :: tau0 = 2.31d0
  double precision :: omega,omegat,stf_val

  ! Brune source-time function
  ! Moment function
  !omegat = 2.d0*PI*f0*t
  !if (t < 0.d0) then
  !  stf_val = 0.d0
  !else if (omegat >= 0.d0 .and. omegat < tau0) then
  !  stf_val = 1.d0 - exp(-omegat)*( 1.0d0 + omegat +  &
  !            0.5d0*omegat**2 - (1.5d0*omegat**3)/tau0 + (1.5d0*omegat**4)/(tau0**2) - (0.5d0*omegat**5)/(tau0**3) )
  !else ! (omegat > tau0) then
  !  stf_val = 1.d0 - exp( -omegat ) * (1.0d0+omegat)
  !endif

  ! Moment rate function
  omega = 2.d0 * PI * f0
  omegat = omega * t

  if (t < 0.d0) then
    stf_val = 0.d0
  else if (omegat >= 0.d0 .and. omegat < tau0) then
    ! 0 <= omega * t < tau0
    stf_val = ( 0.5d0 * omega * (omegat**2) * exp(-omegat)/tau0**3 ) &
              * ( tau0**3 - 3.d0 * (tau0**2) * (omegat-3.d0) + &
                  3.d0 * tau0 * omegat * (omegat-4.d0) - (omegat**2) * (omegat-5.d0) )
  else
    ! omega * t >= tau0
    stf_val = omega * omegat * exp(-omegat)
  endif

  comp_source_time_function_Smooth_Brune = stf_val

  end function comp_source_time_function_Smooth_Brune

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_Yoffe(t,f0,burst_band_width)

! Regularized Yoffe function

  use constants, only: PI

  implicit none

  double precision, intent(in) :: t,f0,burst_band_width

  ! local variables
  double precision :: T_acc,T_eff,tauR,tauS
  double precision :: K_val,stf_val

  ! regularized Yoffe function defined in Appendix (A13) - (A20) of
  ! Tinti et al. (2005),
  ! A kinematic source-time function compatible with earthquake dynamics,
  ! BSSA, 95 (4), 1211-1223. https://doi.org/10.1785/0120040177

  ! fixed parameter example:
  ! acceleration time T_acc (time to peak slip velocity)
  !    T_acc = 0.2
  ! effective duration time T_eff
  !    T_eff = 0.9
  ! rise times tau
  !    tauS = T_acc / 1.27d0         ! tauS - half-duration of the triangular function used for regularizing the Yoffe function
  !    tauR = T_eff - 2.d0 * tauS    ! tauR - Yoffe rise-time

  ! to avoid adding new parameter lines to the DATA/SOURCE files,
  ! here, we re-interprete the source parameters for
  !   frequency        -> T_acc == 1/f0   as slip acceleration duration
  ! and
  !   burst_band_width -> T_eff == 1/bbw  as effective final duration
  T_acc = 1.d0 / f0
  T_eff = 1.d0 / burst_band_width

  ! computes related rise times
  tauS = T_acc / 1.27d0              ! uses factor 1.27 from paper
  tauR = T_eff - 2.d0 * tauS         ! Yoffe rise time

  ! imposes tauS >= 0.0
  if (tauS < 0.d0) tauS = 0.d0

  ! imposes tauR >= tauS
  if (tauR < tauS) tauR = tauS

  ! check rise times
  if (tauR == 0.d0 .or. tauS == 0.d0) then
    comp_source_time_function_Yoffe = 0.d0
    return
  endif

  ! defined only for times t > 0
  if (t <= 0.d0) then
    comp_source_time_function_Yoffe = 0.d0
    return
  endif

  ! constant
  K_val = 2.d0 / (PI * tauR * tauS*tauS)

  ! analytical expressions of regularized Yoffe source time function values
  ! (integrations of triangular function convolved with Yoffe function)
  if (tauR > 2.d0 * tauS) then
    ! (A13)
    if (t <= tauS) then
      stf_val = K_val * (C1() + C2())
    else if (t <= 2.0 * tauS) then
      stf_val = K_val * (C1() - C2() + C3())
    else if (t < tauR) then
      stf_val = K_val * (C1() + C3() + C4())
    else if (t < tauR + tauS) then
      stf_val = K_val * (C3() + C4() + C5())
    else if (t < tauR + 2.0 * tauS) then
      stf_val = K_val * (C4() + C6())
    else
      stf_val = 0.d0
    endif
  else
    ! (A14)
    ! using correction as in SeisSol implementation:
    ! integration boundaries have been fixed here,
    !       from
    !         tauS < t < tauR        instead of   tauS <= t < tauS  (A14, 3. interval case)
    !       from
    !         tauR <= t <= 2 tauS    instead of   tauS <= t < 2 tauR (A14, 4. interval case)
    if (t <= tauS) then
      stf_val = K_val * (C1() + C2())
    else if (t < tauR) then
      stf_val = K_val * (C1() - C2() + C3())
    else if (t <= 2.0 * tauS) then
      stf_val = K_val * (C5() + C3() - C2())
    else if (t < tauR + tauS) then
      stf_val = K_val * (C3() + C4() + C5())
    else if (t < tauR + 2.0 * tauS) then
      stf_val = K_val * (C4() + C6())
    else
      stf_val = 0.d0
    endif
  endif

  ! returns source time function value
  comp_source_time_function_Yoffe = stf_val

contains

  ! C factor functions
  double precision function C1()
    implicit none
    C1 = (0.5d0 * t + 0.25d0 * tauR) * sqrt(t * (tauR - t)) &
         + (t * tauR - tauR*tauR) * asin(sqrt(t / tauR)) &
         - 0.75d0 * tauR*tauR * atan(sqrt((tauR - t) / t))
  end function C1

  double precision function C2()
    implicit none
    C2 = 0.375d0 * PI * tauR*tauR
  end function C2

  double precision function C3()
    implicit none
    C3 = (tauS - t - 0.5d0 * tauR) * sqrt((t - tauS) * (tauR - t + tauS)) &
         + tauR * (2.d0 * tauR - 2.d0 * t + 2.d0 * tauS) * asin(sqrt((t - tauS) / tauR)) &
         + 1.5d0 * tauR*tauR * atan(sqrt((tauR - t + tauS) / (t - tauS)))
  end function C3

  double precision function C4()
    implicit none
    ! applies correction from:
    !   Bizzarri (2012),
    !   Analytical representation of the fault slip velocity from spontaneous dynamic earthquake models,
    !   JGR, 117, doi:10.1029/2011JB009097
    ! see C4 in eq (3).
    ! 2 typos fixed in the second term compared with Tinti et al. 2005, using
    !   .. - tauR * (tauR - t + 2.d0 * tauS) ..
    ! instead of original formula in (A18)
    !   .. - tauR * (tauR + t - 2.d0 * tauS) ..
    ! this is the same correction as in SeisSol implementation:
    !   https://github.com/SeisSol/SeisSol/blob/master/src/Numerical/RegularizedYoffe.h
    C4 = (-tauS + 0.5d0 * t + 0.25d0 * tauR) * sqrt((t - 2.d0 * tauS) * (tauR - t + 2.d0 * tauS)) &
         - tauR * (tauR - t + 2.d0 * tauS) * asin(sqrt((t - 2.d0 * tauS) / tauR)) &
         - 0.75d0 * tauR*tauR * atan(sqrt((tauR - t + 2.d0 * tauS) / (t - 2.d0 * tauS)))
  end function C4

  double precision function C5()
    implicit none
    C5 = 0.5d0 * PI * tauR * (t - tauR)
  end function C5

  double precision function C6()
    implicit none
    C6 = 0.5d0 * PI * tauR * (2.d0 * tauS - t + tauR)
  end function C6

  end function comp_source_time_function_Yoffe

!
!------------------------------------------------------------
!

  double precision function comp_source_time_function_Yoffe_integrated(t,f0,burst_band_width)

! integrated regularized Yoffe function using Simpson's rule

  use constants, only: PI

  use shared_parameters, only: DT

  implicit none

  double precision, intent(in) :: t,f0,burst_band_width

  ! local variables
  double precision :: T_acc,T_eff,tauR,tauS
  double precision :: integrated_val
  integer :: n_steps, i
  double precision :: dt_step, t_current, stf_val

  ! external function
  double precision, external :: comp_source_time_function_Yoffe

! note: The regularized Yoffe function is a moment-rate (or slip rate or slip velocity) function.
!       Here, we provide the moment function that is the integral from 0 to t of the Yoffe function.
!       This is thus the moment (or cummulative slip) function that applies to CMT sources.
!       Integration is done numerically using Simpson's 1/3 rule.

  ! defined only for times t > 0
  if (t <= 0.d0) then
    comp_source_time_function_Yoffe_integrated = 0.d0
    return
  endif

  ! regularized Yoffe function parameters
  T_acc = 1.d0 / f0
  T_eff = 1.d0 / burst_band_width

  tauS = T_acc / 1.27d0
  tauR = T_eff - 2.d0 * tauS

  ! imposes tauS >= 0.0
  if (tauS < 0.d0) tauS = 0.d0

  ! imposes tauR >= tauS
  if (tauR < tauS) tauR = tauS

  ! check rise times
  if (tauR == 0.d0 .or. tauS == 0.d0) then
    comp_source_time_function_Yoffe_integrated = 0.d0
    return
  endif

  ! determine integration parameters (must be even for Simpson's rule)
  n_steps = max(100, int(t / (tauS * 0.01d0)))
  n_steps = max(n_steps, int(t / DT))         ! to double-check enough fine stepping in case tauS is large

  ! ensure even number
  if (mod(n_steps, 2) == 1) n_steps = n_steps + 1
  dt_step = t / dble(n_steps)

  ! Simpson's rule integration
  integrated_val = 0.d0

  ! first point
  stf_val = 0.d0  ! or comp_source_time_function_Yoffe(0.d0, f0, burst_band_width)
                  !    integrated_val = integrated_val + stf_val

  ! intermediate points
  do i = 1, n_steps - 1
    t_current = dble(i) * dt_step
    stf_val = comp_source_time_function_Yoffe(t_current, f0, burst_band_width)
    if (mod(i, 2) == 1) then
      ! odd indices
      integrated_val = integrated_val + 4.d0 * stf_val
    else
      ! even indices
      integrated_val = integrated_val + 2.d0 * stf_val
    endif
  enddo

  ! last point
  stf_val = comp_source_time_function_Yoffe(t, f0, burst_band_width)
  integrated_val = integrated_val + stf_val

  ! apply Simpson's rule factor
  integrated_val = integrated_val * dt_step / 3.d0

  ! returns integrated source time function value
  comp_source_time_function_Yoffe_integrated = integrated_val

  end function comp_source_time_function_Yoffe_integrated
