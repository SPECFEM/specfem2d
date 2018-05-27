
  program analytical_solution

! This program implements the analytical solution for the displacement vector in a 2D plane-strain elastic medium
! with a vertical force source located in (0,0),
! from Appendix B of Carcione et al., Wave propagation simulation in a linear viscoelastic medium, GJI, vol. 95, p. 597-611 (1988),
! modified to compute the solution in the time domain instead of in the frequency domain
! (and thus for the elastic case only, instead of for the viscoelastic case).
! The amplitude of the force is called F and is defined below.

!! DK DK program written by Dimitri Komatitsch, CNRS, Marseille, France, May 2018.

  implicit none

! to see how small the contribution of the near-field term is,
! here the user can ask not to include it, to then compare with the full result obtained with this flag set to false
  logical, parameter :: DO_NOT_COMPUTE_THE_NEAR_FIELD = .false.

! for time-domain calculation

! number of time steps used to discretize the time interval
! use a very high value here, because the convolution algorithm that I implemented is very basic
! and thus needs a tiny time interval in order to be accurate
  integer, parameter :: nt = 100000

  double precision, parameter :: Tmax = 1.d0 ! seconds
  integer :: it
  integer :: j
  double precision :: tau_j
  double precision :: ricker,convolution_value,a
  double precision, dimension(nt) :: Green

  double precision, parameter :: pi = 3.141592653589793d0

! density of the medium
  double precision, parameter :: rho = 2000.d0

! unrelaxed (f = +infinity) values
! these values for the unrelaxed state are computed from the relaxed state values (Vp = 3000, Vs = 2000, rho = 2000)
! given in Carcione et al. 1988 GJI vol 95 p 604 Table 1
  double precision, parameter :: Vp = 3297.849d0
  double precision, parameter :: Vs = 2222.536d0

! amplitude of the force source
  double precision, parameter :: F = 1.d0

! definition position recepteur Carcione
  double precision x1,x2

! Definition source Dimitri
  double precision, parameter :: f0 = 18.d0
  double precision, parameter :: t0 = 1.2d0 / f0

  double precision :: deltat,time

! external functions
  double precision, external :: u1,u2

! ********** end of variable declarations ************

!! duration of a time step
  deltat = Tmax / (nt-1)

! position of the receiver
  x1 = +500.
  x2 = +500.

  print *,'Force source located at the origin (0,0)'
  print *,'Receiver located in (x,z) = ',x1,x2

  if (DO_NOT_COMPUTE_THE_NEAR_FIELD) then
    print *,'BEWARE: computing the far-field solution only, rather than the full Green function'
  else
    print *,'Computing the full solution, including the near-field term of the Green function'
  endif

! **********
! Compute Ux
! **********

! save time result for Ux
  if (DO_NOT_COMPUTE_THE_NEAR_FIELD) then
    open(unit=11,file='Ux_time_analytical_solution_elastic_without_near_field_time_domain.dat',status='unknown')
  else
    open(unit=11,file='Ux_time_analytical_solution_elastic_time_domain.dat',status='unknown')
  endif

! store the Green function
  do it = 1,nt
! subtract t0 to be consistent with the SPECFEM2D code
    time = dble(it-1)*deltat - t0
    Green(it) = u1(time,Vp,Vs,x1,x2,rho,F,DO_NOT_COMPUTE_THE_NEAR_FIELD)
  enddo

! to avoid writing a huge file, since we purposely used a huge number of time steps, write only every 10 time steps
  do it = 1,nt,10

! subtract t0 to be consistent with the SPECFEM2D code
    time = dble(it-1)*deltat - t0

    if (mod(it-1,5000) == 0) print *,'For Ux, computing ',it,' out of ',nt

    convolution_value = 0.d0

    do j = it-nt,it-1

      tau_j = dble(j)*deltat

! convolve with a Ricker wavelet
      a = PI**2 * f0**2
      ricker = (1.d0 - 2.d0 * a * tau_j*tau_j) * exp( -a * tau_j*tau_j )

      convolution_value = convolution_value + Green(it-j)*ricker*deltat

    enddo

! to avoid writing a huge file, since we purposely used a huge number of time steps, write only every 10 time steps
    if (time >= 0.d0) write(11,*) sngl(time),sngl(convolution_value)

  enddo

  close(11)

! **********
! Compute Uz
! **********

! save time result for Uz
  if (DO_NOT_COMPUTE_THE_NEAR_FIELD) then
    open(unit=11,file='Uz_time_analytical_solution_elastic_without_near_field_time_domain.dat',status='unknown')
  else
    open(unit=11,file='Uz_time_analytical_solution_elastic_time_domain.dat',status='unknown')
  endif

! store the Green function
  do it = 1,nt
! subtract t0 to be consistent with the SPECFEM2D code
    time = dble(it-1)*deltat - t0
    Green(it) = u2(time,Vp,Vs,x1,x2,rho,F,DO_NOT_COMPUTE_THE_NEAR_FIELD)
  enddo

! to avoid writing a huge file, since we purposely used a huge number of time steps, write only every 10 time steps
  do it = 1,nt,10

! subtract t0 to be consistent with the SPECFEM2D code
    time = dble(it-1)*deltat - t0

    if (mod(it-1,5000) == 0) print *,'For Uz, computing ',it,' out of ',nt

    convolution_value = 0.d0

    do j = it-nt,it-1

      tau_j = dble(j)*deltat

! convolve with a Ricker wavelet
      a = PI**2 * f0**2
      ricker = (1.d0 - 2.d0 * a * tau_j*tau_j) * exp( -a * tau_j*tau_j )

      convolution_value = convolution_value + Green(it-j)*ricker*deltat

    enddo

    if (time >= 0.d0) write(11,*) sngl(time),sngl(convolution_value)

  enddo

  close(11)

  end

! -----------

  double precision function u1(t,v1,v2,x1,x2,rho,F,DO_NOT_COMPUTE_THE_NEAR_FIELD)

  implicit none

  double precision t
  double precision v1,v2

  logical :: DO_NOT_COMPUTE_THE_NEAR_FIELD

  double precision G1,G2
  external G1,G2

  double precision, parameter :: pi = 3.141592653589793d0

! amplitude of the force
  double precision F

  double precision x1,x2,r,rho

! source-receiver distance
  r = dsqrt(x1**2 + x2**2)

  u1 = F * x1 * x2 * (G1(r,t,v1,v2,DO_NOT_COMPUTE_THE_NEAR_FIELD) + &
                      G2(r,t,v1,v2,DO_NOT_COMPUTE_THE_NEAR_FIELD)) / (2.d0 * pi * rho * r**2)

  end

! -----------

  double precision function u2(t,v1,v2,x1,x2,rho,F,DO_NOT_COMPUTE_THE_NEAR_FIELD)

  implicit none

  double precision t
  double precision v1,v2

  logical :: DO_NOT_COMPUTE_THE_NEAR_FIELD

  double precision G1,G2
  external G1,G2

  double precision, parameter :: pi = 3.141592653589793d0

! amplitude of the force
  double precision F

  double precision x1,x2,r,rho

! source-receiver distance
  r = dsqrt(x1**2 + x2**2)

  u2 = F * (x2*x2*G1(r,t,v1,v2,DO_NOT_COMPUTE_THE_NEAR_FIELD) - &
            x1*x1*G2(r,t,v1,v2,DO_NOT_COMPUTE_THE_NEAR_FIELD)) / (2.d0 * pi * rho * r**2)

  end

! -----------

  double precision function G1(r,t,v1,v2,DO_NOT_COMPUTE_THE_NEAR_FIELD)

  implicit none

  double precision r,t
  double precision v1,v2

  logical :: DO_NOT_COMPUTE_THE_NEAR_FIELD

  integer heaviside
  external heaviside

  double precision :: tau1,tau2
  integer :: heaviside_tau1,heaviside_tau2

  double precision, parameter :: pi = 3.141592653589793d0

! from equation (B2a) of Carcione et al., Wave propagation simulation in a linear viscoelastic medium,
! Geophysical Journal, vol. 95, p. 597-611 (1988)

  tau1 = r / v1
  tau2 = r / v2

  heaviside_tau1 = heaviside(t - tau1)
  heaviside_tau2 = heaviside(t - tau2)

  G1 = 0.d0
  if (heaviside_tau1 == 1) G1 = G1 + 1.d0 / (v1**2 * sqrt(t**2 - tau1**2))
  if (.not. DO_NOT_COMPUTE_THE_NEAR_FIELD) then
    if (heaviside_tau1 == 1) G1 = G1 + sqrt(t**2 - tau1**2) / (r**2)
    if (heaviside_tau2 == 1) G1 = G1 - sqrt(t**2 - tau2**2) / (r**2)
  endif

  end

! -----------

  double precision function G2(r,t,v1,v2,DO_NOT_COMPUTE_THE_NEAR_FIELD)

  implicit none

  double precision r,t
  double precision v1,v2

  logical :: DO_NOT_COMPUTE_THE_NEAR_FIELD

  integer heaviside
  external heaviside

  double precision :: tau1,tau2
  integer :: heaviside_tau1,heaviside_tau2

  double precision, parameter :: pi = 3.141592653589793d0

! from equation (B2a) of Carcione et al., Wave propagation simulation in a linear viscoelastic medium,
! Geophysical Journal, vol. 95, p. 597-611 (1988)

  tau1 = r / v1
  tau2 = r / v2

  heaviside_tau1 = heaviside(t - tau1)
  heaviside_tau2 = heaviside(t - tau2)

  G2 = 0.d0
  if (heaviside_tau2 == 1) G2 = G2 - 1.d0 / (v2**2 * sqrt(t**2 - tau2**2))
  if (.not. DO_NOT_COMPUTE_THE_NEAR_FIELD) then
    if (heaviside_tau1 == 1) G2 = G2 + sqrt(t**2 - tau1**2) / (r**2)
    if (heaviside_tau2 == 1) G2 = G2 - sqrt(t**2 - tau2**2) / (r**2)
  endif

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

