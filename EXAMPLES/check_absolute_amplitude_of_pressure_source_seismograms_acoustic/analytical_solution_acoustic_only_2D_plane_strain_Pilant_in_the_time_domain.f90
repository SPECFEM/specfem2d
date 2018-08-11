
  program analytical_solution

! This program implements the analytical solution for pressure in a 2D plane-strain acoustic medium
! with a pressure source located in (0,0),
! from Appendix B of Carcione et al., Wave propagation simulation in a linear viscoacoustic medium,
! Geophysical Journal, vol. 93, p. 393-407 (1988),
! modified to compute the solution in the time domain instead of in the frequency domain
! (and thus for the acoustic case only, instead of for the viscoacoustic case).
! The amplitude of the source is called F and is defined below.

!! DK DK program written by Dimitri Komatitsch, CNRS, Marseille, France, May 2018.

  implicit none

! for time-domain calculation

! number of time steps used to discretize the time interval
! use a very high value here, because the convolution algorithm that I implemented is very basic
! and thus needs a tiny time interval in order to be accurate
  integer, parameter :: nt = 1000000

  double precision, parameter :: Tmax = 1.d0 ! seconds
  integer :: it
  integer :: j
  double precision :: tau_j
  double precision :: ricker,convolution_value,a,r
  double precision, dimension(nt) :: Green

  double precision, parameter :: pi = 3.141592653589793d0

! unrelaxed (f = +infinity) values
! this value for the unrelaxed state is computed from the relaxed state value (Vp = 3000)
! given in Carcione et al. 1988 GJI vol 95 p 604 Table 1
  double precision, parameter :: Vp = 3297.849d0

! amplitude of the force source
  double precision, parameter :: F = 1.d0

! definition position recepteur Carcione
  double precision x1,x2

! Definition source Dimitri
  double precision, parameter :: f0 = 18.d0
  double precision, parameter :: t0 = 1.2d0 / f0

  double precision :: deltat,time

! external functions
  double precision, external :: G1

! ********** end of variable declarations ************

!! duration of a time step
  deltat = Tmax / (nt-1)

! position of the receiver
  x1 = +500.d0
  x2 = +500.d0

! x1 = +1.d0
! x2 = +1.d0

! x1 = +0.1d0
! x2 = +0.1d0

  print *,'Force source located at the origin (0,0)'
  print *,'Receiver located in (x,z) = ',x1,x2

! ****************
! Compute pressure
! ****************

! save time result for pressure
  open(unit=11,file='pressure_time_analytical_solution_acoustic_time_domain.dat',status='unknown')

! source-receiver distance
  r = dsqrt(x1**2 + x2**2)

! store the Green function
  do it = 1,nt
! subtract t0 to be consistent with the SPECFEM2D code
    time = dble(it-1)*deltat - t0
    Green(it) = F * G1(r,time,Vp)
  enddo

!! DK DK to compare to our finite-difference codes from SEISMIC_CPML or SOUNDVIEW,
!! DK DK we divide the source by 4 * PI * cp^2 to get the right amplitude (our convention being to use a source of amplitude 1,
!! DK DK while the convention used by Carcione in his 1988 paper is to use a source of amplitude 4 * PI * cp^2
  Green(:) = Green(:) / (4.d0 * PI * Vp**2)

! to avoid writing a huge file, since we purposely used a huge number of time steps, write only every 100 time steps
  do it = 1,nt,100

! subtract t0 to be consistent with the SPECFEM2D code
    time = dble(it-1)*deltat - t0

    if (mod(it-1,5000) == 0) print *,'For pressure, computing ',it,' out of ',nt

    convolution_value = 0.d0

    do j = it-nt,it-1

      tau_j = dble(j)*deltat

! convolve with a Ricker wavelet
      a = PI**2 * f0**2
      ricker = (1.d0 - 2.d0 * a * tau_j*tau_j) * exp( -a * tau_j*tau_j )

      convolution_value = convolution_value + Green(it-j)*ricker*deltat

    enddo

    if (time >= 0.d0) write(11,*) sngl(time),sngl(convolution_value)
!   write(11,*) sngl(time),sngl(convolution_value)

  enddo

  close(11)

  end

! -----------

  double precision function G1(r,t,v1)

  implicit none

  double precision r,t
  double precision v1

  integer heaviside
  external heaviside

  double precision :: tau1
  integer :: heaviside_tau1

  double precision, parameter :: pi = 3.141592653589793d0

! equation (B2) of Carcione et al., Wave propagation simulation in a linear viscoacoustic medium,
! Geophysical Journal, vol. 93, p. 393-407 (1988)

  tau1 = r / v1

  heaviside_tau1 = heaviside(t - tau1)

  G1 = 0.d0
  if (heaviside_tau1 == 1) G1 = G1 + 2.d0 / sqrt(t**2 - tau1**2)

  end

! -----------

  integer function heaviside(t)

  implicit none

  double precision t

  if (t > 0.d0) then
    heaviside = 1
  else
    heaviside = 0
  endif

  end

