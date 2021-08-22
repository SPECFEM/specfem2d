
  program analytical_solution

! This program implements the analytical solution for the displacement vector in a 2D plane-strain viscoelastic medium
! with a vertical force source located in (0,0),
! from Appendix B of Carcione et al., Wave propagation simulation in a linear viscoelastic medium, GJI, vol. 95, p. 597-611 (1988)
! (note that that Appendix contains two typos, fixed in this code; I added two comments below to mention them).
! The amplitude of the force is called F and is defined below.

  implicit none

!! DK DK May 2018: the missing 1/L factor in older Carcione papers
!! DK DK May 2018: has been added to this code by Quentin Brissaud and by Etienne Bachmann
!! DK DK for the viscoacoustic code in directory EXAMPLES/attenuation/viscoacoustic,
!! DK DK it would be very easy to copy the changes from there to this viscoelastic version;
!! DK DK but then all the values of the tau_epsilon in the code below would need to change.

!! DK DK Dimitri Komatitsch, CNRS Marseille, France, April 2017: added the elastic reference calculation.

! compute the non-viscoacoustic case as a reference if needed, i.e. turn attenuation off
  logical, parameter :: TURN_ATTENUATION_OFF = .true. ! .false.

! to see how small the contribution of the near-field term is,
! here the user can ask not to include it, to then compare with the full result obtained with this flag set to false
  logical, parameter :: DO_NOT_COMPUTE_THE_NEAR_FIELD = .false.

  integer, parameter :: iratio = 32

  integer, parameter :: nfreq = 524288
  integer, parameter :: nt = iratio * nfreq

  double precision, parameter :: freqmax = 80.d0
!! DK DK to print the velocity if we want to display the curve of how velocity varies with frequency
!! DK DK for instance to compute the unrelaxed velocity in the Zener model
! double precision, parameter :: freqmax = 20000.d0

  double precision, parameter :: freqseuil = 0.00005d0

  double precision, parameter :: pi = 3.141592653589793d0

! for the solution in time domain
  integer it,i
  real wsave(4*nt+15)
  complex c(nt)

!! DK DK for my slow inverse Discrete Fourier Transform using a double loop
  complex :: input(nt), i_imaginary_constant
  integer :: j,m

! density of the medium
  double precision, parameter :: rho = 2000.d0

! unrelaxed (f = +infinity) values
! these values for the unrelaxed state are computed from the relaxed state values (Vp = 3000, Vs = 2000, rho = 2000)
! given in Carcione et al. 1988 GJI vol 95 p 604 Table 1
  double precision, parameter :: Vp = 3297.849d0
  double precision, parameter :: Vs = 2222.536d0

! unrelaxed (f = +infinity) values, i.e. using the fastest Vp and Vs velocities
  double precision, parameter :: M2_unrelaxed = Vs**2 * 2.d0 * rho
  double precision, parameter :: M1_unrelaxed = 2.d0 * Vp**2 * rho - M2_unrelaxed

! amplitude of the force source
  double precision, parameter :: F = 1.d0

! definition position recepteur Carcione
  double precision x1,x2

! Definition source Dimitri
  double precision, parameter :: f0 = 18.d0
  double precision, parameter :: t0 = 1.2d0 / f0

! Definition source Carcione
! double precision f0,t0,eta,epsil
! parameter(f0 = 50.d0)
! parameter(t0 = 0.075d0)
! parameter(epsil = 1.d0)
! parameter(eta = 0.5d0)

! number of Zener standard linear solids in parallel
  integer, parameter :: Lnu = 3

! DK DK I implemented a very simple and slow inverse Discrete Fourier Transform
! DK DK at some point, for verification, using a double loop. I keep it just in case.
! DK DK For large number of points it is extremely slow because of the double loop.
! DK DK Thus there is no reason to turn this flag on.
  logical, parameter :: USE_SLOW_FOURIER_TRANSFORM = .false.

!! DK DK March 2018: this missing 1/L factor has been added to this code by Quentin Brissaud
!! DK DK for the viscoacoustic code in directory EXAMPLES/attenuation/viscoacoustic,
!! DK DK it would be very easy to copy the changes from there to this viscoelastic version;
!! DK DK but then all the values of the tau_epsilon below would need to change.

 double precision, dimension(Lnu) :: tau_sigma_nu1,tau_sigma_nu2,tau_epsilon_nu1,tau_epsilon_nu2

  integer :: ifreq
  double precision :: deltafreq,freq,omega,omega0,deltat,time,a
  double complex :: comparg

! Fourier transform of the Ricker wavelet source
  double complex fomega(0:nfreq)

! real and imaginary parts
  double precision ra(0:nfreq),rb(0:nfreq)

! spectral amplitude
  double precision ampli(0:nfreq)

! analytical solution for the two components
  double complex phi1(-nfreq:nfreq)
  double complex phi2(-nfreq:nfreq)

! external functions
  double complex, external :: u1,u2

! modules elastiques
  double complex :: M1C, M2C, E, V1, V2, temp

! ********** end of variable declarations ************

! classical least-squares constants
  tau_epsilon_nu1 =  (/ 0.109527114743452     ,  1.070028707488438E-002,  1.132519034287800E-003/)
  tau_sigma_nu1 = (/  8.841941282883074E-002 , 8.841941282883075E-003,  8.841941282883074E-004/)
  tau_epsilon_nu2 = (/  0.112028084581976    ,   1.093882462934487E-002,  1.167173427475064E-003/)
  tau_sigma_nu2 = (/  8.841941282883074E-002,  8.841941282883075E-003,  8.841941282883074E-004/)

! position of the receiver
  x1 = +500.
  x2 = +500.

  print *,'Force source located at the origin (0,0)'
  print *,'Receiver located in (x,z) = ',x1,x2

  if (TURN_ATTENUATION_OFF) then
    print *,'BEWARE: computing the elastic reference solution (i.e., without attenuation) instead of the viscoelastic solution'
  else
    print *,'Computing the viscoelastic solution'
  endif

  if (DO_NOT_COMPUTE_THE_NEAR_FIELD) then
    print *,'BEWARE: computing the far-field solution only, rather than the full Green function'
  else
    print *,'Computing the full solution, including the near-field term of the Green function'
  endif

! step in frequency
  deltafreq = freqmax / dble(nfreq)

! define parameters for the Ricker source
  omega0 = 2.d0 * pi * f0
  a = pi**2 * f0**2

  deltat = 1.d0 / (freqmax*dble(iratio))

! define the spectrum of the source
  do ifreq=0,nfreq
      freq = deltafreq * dble(ifreq)
      omega = 2.d0 * pi * freq

! typo in equation (B7) of Carcione et al., Wave propagation simulation in a linear viscoelastic medium,
! Geophysical Journal, vol. 95, p. 597-611 (1988), the exponential should be of -i omega t0,
! fixed here by adding the minus sign
      comparg = dcmplx(0.d0,-omega*t0)

! definir le spectre du Ricker de Carcione avec cos()
! equation (B7) of Carcione et al., Wave propagation simulation in a linear viscoelastic medium,
! Geophysical Journal, vol. 95, p. 597-611 (1988)
!     fomega(ifreq) = pi * dsqrt(pi/eta) * (1.d0/omega0) * cdexp(comparg) * ( dexp(- (pi*pi/eta) * (epsil/2 - omega/omega0)**2) &
!         + dexp(- (pi*pi/eta) * (epsil/2 + omega/omega0)**2) )

! definir le spectre d'un Ricker classique (centre en t0)
      fomega(ifreq) = dsqrt(pi) * cdexp(comparg) * omega**2 * dexp(-omega**2/(4.d0*a)) / (2.d0 * dsqrt(a**3))

      ra(ifreq) = dreal(fomega(ifreq))
      rb(ifreq) = dimag(fomega(ifreq))
! prendre le module de l'amplitude spectrale
      ampli(ifreq) = dsqrt(ra(ifreq)**2 + rb(ifreq)**2)
  enddo

! sauvegarde du spectre d'amplitude de la source en Hz au format Gnuplot
  open(unit=10,file='spectrum_of_the_source_used.gnu',status='unknown')
  do ifreq = 0,nfreq
    freq = deltafreq * dble(ifreq)
    write(10,*) sngl(freq),sngl(ampli(ifreq))
  enddo
  close(10)

! ************** calcul solution analytique ****************

! d'apres Carcione GJI vol 95 p 611 (1988)
  do ifreq=0,nfreq
      freq = deltafreq * dble(ifreq)
      omega = 2.d0 * pi * freq

! critere ad-hoc pour eviter singularite en zero
  if (freq < freqseuil) omega = 2.d0 * pi * freqseuil

! use standard infinite frequency (unrelaxed) reference,
! in which waves slow down when attenuation is turned on.
  temp = dcmplx(0.d0,0.d0)
  do i=1,Lnu
    temp = temp + dcmplx(1.d0,omega*tau_epsilon_nu1(i)) / dcmplx(1.d0,omega*tau_sigma_nu1(i))
  enddo

  M1C = (M1_unrelaxed /(sum(tau_epsilon_nu1(:)/tau_sigma_nu1(:)))) * temp

  temp = dcmplx(0.d0,0.d0)
  do i=1,Lnu
    temp = temp + dcmplx(1.d0,omega*tau_epsilon_nu2(i)) / dcmplx(1.d0,omega*tau_sigma_nu2(i))
  enddo

  M2C = (M2_unrelaxed /(sum(tau_epsilon_nu2(:)/tau_sigma_nu2(:)))) * temp

  if (TURN_ATTENUATION_OFF) then
! from Etienne Bachmann, May 2018: pour calculer la solution sans attenuation, il faut donner le Mu_unrelaxed et pas le Mu_relaxed.
! En effet, pour comparer avec SPECFEM, il faut simplement partir de la bonne reference.
! SPECFEM est defini en unrelaxed et les constantes unrelaxed dans Carcione matchent parfaitement les Vp et Vs definis dans SPECFEM.
    M1C = M1_unrelaxed
    M2C = M2_unrelaxed
  endif

  E = (M1C + M2C) / 2
  V1 = cdsqrt(E / rho)  !! DK DK this is Vp
!! DK DK print the velocity if we want to display the curve of how velocity varies with frequency
!! DK DK for instance to compute the unrelaxed velocity in the Zener model
! print *,freq,dsqrt(real(V1)**2 + imag(V1)**2)
  V2 = cdsqrt(M2C / (2.d0 * rho))  !! DK DK this is Vs
!! DK DK print the velocity if we want to display the curve of how velocity varies with frequency
!! DK DK for instance to compute the unrelaxed velocity in the Zener model
! print *,freq,dsqrt(real(V2)**2 + imag(V2)**2)

! calcul de la solution analytique en frequence
  phi1(ifreq) = u1(omega,V1,V2,x1,x2,rho,F,DO_NOT_COMPUTE_THE_NEAR_FIELD) * fomega(ifreq)
  phi2(ifreq) = u2(omega,V1,V2,x1,x2,rho,F,DO_NOT_COMPUTE_THE_NEAR_FIELD) * fomega(ifreq)

  enddo

! take the conjugate value for negative frequencies
  do ifreq=-nfreq,-1
      phi1(ifreq) = dconjg(phi1(-ifreq))
      phi2(ifreq) = dconjg(phi2(-ifreq))
  enddo

! save the result in the frequency domain
! open(unit=11,file='cmplx_phi',status='unknown')
! do ifreq=-nfreq,nfreq
!     freq = deltafreq * dble(ifreq)
!     write(11,*) sngl(freq),sngl(dreal(phi1(ifreq))),sngl(dimag(phi1(ifreq))),sngl(dreal(phi2(ifreq))),sngl(dimag(phi2(ifreq)))
! enddo
! close(11)

! ***************************************************************************
! Calculation of the time domain solution (using routine "cfftb" from Netlib)
! ***************************************************************************

! **********
! Compute Ux
! **********

! initialize FFT arrays
  call cffti(nt,wsave)

! clear array of Fourier coefficients
  do it = 1,nt
    c(it) = cmplx(0.,0.)
  enddo

! use the Fourier values for Ux
  c(1) = cmplx(phi1(0))
  do ifreq=1,nfreq-2
      c(ifreq+1) = cmplx(phi1(ifreq))
      c(nt+1-ifreq) = conjg(cmplx(phi1(ifreq)))
  enddo

! perform the inverse FFT for Ux
  if (.not. USE_SLOW_FOURIER_TRANSFORM) then
    call cfftb(nt,c,wsave)
  else
! DK DK I implemented a very simple and slow inverse Discrete Fourier Transform here
! DK DK at some point, for verification, using a double loop. I keep it just in case.
! DK DK For large number of points it is extremely slow because of the double loop.
    input(:) = c(:)
!   imaginary constant "i"
    i_imaginary_constant = (0.,1.)
    do it = 1,nt
      if (mod(it,1000) == 0) print *,'FFT inverse it = ',it,' out of ',nt
      j = it
      c(j) = cmplx(0.,0.)
      do m = 1,nt
        c(j) = c(j) + input(m) * exp(2.d0 * PI * i_imaginary_constant * dble((m-1) * (j-1)) / nt)
      enddo
    enddo
  endif

! in the inverse Discrete Fourier transform one needs to divide by N, the number of samples (number of time steps here)
  c(:) = c(:) / nt

! value of a time step
  deltat = 1.d0 / (freqmax*dble(iratio))

! to get the amplitude right, we need to divide by the time step
  c(:) = c(:) / deltat

! save time result inverse FFT for Ux

  if (TURN_ATTENUATION_OFF) then
    if (DO_NOT_COMPUTE_THE_NEAR_FIELD) then
      open(unit=11,file='Ux_time_analytical_solution_elastic_without_near_field.dat',status='unknown')
    else
      open(unit=11,file='Ux_time_analytical_solution_elastic.dat',status='unknown')
    endif
  else
    if (DO_NOT_COMPUTE_THE_NEAR_FIELD) then
      open(unit=11,file='Ux_time_analytical_solution_viscoelastic_without_near_field.dat',status='unknown')
    else
      open(unit=11,file='Ux_time_analytical_solution_viscoelastic.dat',status='unknown')
    endif
  endif
  do it=1,nt
! DK DK Dec 2011: subtract t0 to be consistent with the SPECFEM2D code
        time = dble(it-1)*deltat - t0
! the seismograms are very long due to the very large number of FFT points used,
! thus keeping the useful part of the signal only (the first six seconds of the seismogram)
        if (time >= 0.d0 .and. time <= 6.d0) write(11,*) sngl(time),real(c(it))
  enddo
  close(11)

! **********
! Compute Uz
! **********

! clear array of Fourier coefficients
  do it = 1,nt
    c(it) = cmplx(0.,0.)
  enddo

! use the Fourier values for Uz
  c(1) = cmplx(phi2(0))
  do ifreq=1,nfreq-2
      c(ifreq+1) = cmplx(phi2(ifreq))
      c(nt+1-ifreq) = conjg(cmplx(phi2(ifreq)))
  enddo

! perform the inverse FFT for Ux
  if (.not. USE_SLOW_FOURIER_TRANSFORM) then
    call cfftb(nt,c,wsave)
  else
! DK DK I implemented a very simple and slow inverse Discrete Fourier Transform here
! DK DK at some point, for verification, using a double loop. I keep it just in case.
! DK DK For large number of points it is extremely slow because of the double loop.
    input(:) = c(:)
!   imaginary constant "i"
    i_imaginary_constant = (0.,1.)
    do it = 1,nt
      if (mod(it,1000) == 0) print *,'FFT inverse it = ',it,' out of ',nt
      j = it
      c(j) = cmplx(0.,0.)
      do m = 1,nt
        c(j) = c(j) + input(m) * exp(2.d0 * PI * i_imaginary_constant * dble((m-1) * (j-1)) / nt)
      enddo
    enddo
  endif

! in the inverse Discrete Fourier transform one needs to divide by N, the number of samples (number of time steps here)
  c(:) = c(:) / nt

! value of a time step
  deltat = 1.d0 / (freqmax*dble(iratio))

! to get the amplitude right, we need to divide by the time step
  c(:) = c(:) / deltat

! save time result inverse FFT for Uz
  if (TURN_ATTENUATION_OFF) then
    if (DO_NOT_COMPUTE_THE_NEAR_FIELD) then
      open(unit=11,file='Uz_time_analytical_solution_elastic_without_near_field.dat',status='unknown')
    else
      open(unit=11,file='Uz_time_analytical_solution_elastic.dat',status='unknown')
    endif
  else
    if (DO_NOT_COMPUTE_THE_NEAR_FIELD) then
      open(unit=11,file='Uz_time_analytical_solution_viscoelastic_without_near_field.dat',status='unknown')
    else
      open(unit=11,file='Uz_time_analytical_solution_viscoelastic.dat',status='unknown')
    endif
  endif
  do it=1,nt
! DK DK Dec 2011: subtract t0 to be consistent with the SPECFEM2D code
        time = dble(it-1)*deltat - t0
! the seismograms are very long due to the very large number of FFT points used,
! thus keeping the useful part of the signal only (the first six seconds of the seismogram)
        if (time >= 0.d0 .and. time <= 6.d0) write(11,*) sngl(time),real(c(it))
  enddo
  close(11)

  end

! -----------

  double complex function u1(omega,v1,v2,x1,x2,rho,F,DO_NOT_COMPUTE_THE_NEAR_FIELD)

  implicit none

  double precision omega
  double complex v1,v2

  logical :: DO_NOT_COMPUTE_THE_NEAR_FIELD

  double complex G1,G2
  external G1,G2

  double precision, parameter :: pi = 3.141592653589793d0

! amplitude of the force
  double precision F

  double precision x1,x2,r,rho

! source-receiver distance
  r = dsqrt(x1**2 + x2**2)

  u1 = F * x1 * x2 * (G1(r,omega,v1,v2,DO_NOT_COMPUTE_THE_NEAR_FIELD) + &
                      G2(r,omega,v1,v2,DO_NOT_COMPUTE_THE_NEAR_FIELD)) / (2.d0 * pi * rho * r**2)

  end

! -----------

  double complex function u2(omega,v1,v2,x1,x2,rho,F,DO_NOT_COMPUTE_THE_NEAR_FIELD)

  implicit none

  double precision omega
  double complex v1,v2

  logical :: DO_NOT_COMPUTE_THE_NEAR_FIELD

  double complex G1,G2
  external G1,G2

  double precision, parameter :: pi = 3.141592653589793d0

! amplitude of the force
  double precision F

  double precision x1,x2,r,rho

! source-receiver distance
  r = dsqrt(x1**2 + x2**2)

  u2 = F * (x2*x2*G1(r,omega,v1,v2,DO_NOT_COMPUTE_THE_NEAR_FIELD) - &
            x1*x1*G2(r,omega,v1,v2,DO_NOT_COMPUTE_THE_NEAR_FIELD)) / (2.d0 * pi * rho * r**2)

  end

! -----------

  double complex function G1(r,omega,v1,v2,DO_NOT_COMPUTE_THE_NEAR_FIELD)

  implicit none

  double precision r,omega
  double complex v1,v2

  logical :: DO_NOT_COMPUTE_THE_NEAR_FIELD

  double complex hankel0,hankel1
  external hankel0,hankel1

  double precision, parameter :: pi = 3.141592653589793d0

! typo in equations (B4a) and (B4b) of Carcione et al., Wave propagation simulation in a linear viscoelastic medium,
! Geophysical Journal, vol. 95, p. 597-611 (1988), fixed here: omega/(r*v) -> omega*r/v

  if (DO_NOT_COMPUTE_THE_NEAR_FIELD) then
   G1 = (hankel0(omega*r/v1)/(v1**2)) * dcmplx(0.d0,-pi/2.d0)
  else
   G1 = (hankel0(omega*r/v1)/(v1**2) + hankel1(omega*r/v2)/(omega*r*v2) - hankel1(omega*r/v1)/(omega*r*v1)) * dcmplx(0.d0,-pi/2.d0)
  endif

  end

! -----------

  double complex function G2(r,omega,v1,v2,DO_NOT_COMPUTE_THE_NEAR_FIELD)

  implicit none

  double precision r,omega
  double complex v1,v2

  logical :: DO_NOT_COMPUTE_THE_NEAR_FIELD

  double complex hankel0,hankel1
  external hankel0,hankel1

  double precision, parameter :: pi = 3.141592653589793d0

! typo in equations (B4a) and (B4b) of Carcione et al., Wave propagation simulation in a linear viscoelastic medium,
! Geophysical Journal, vol. 95, p. 597-611 (1988), fixed here: omega/(r*v) -> omega*r/v

  if (DO_NOT_COMPUTE_THE_NEAR_FIELD) then
   G2 = (hankel0(omega*r/v2)/(v2**2)) * dcmplx(0.d0,+pi/2.d0)
  else
   G2 = (hankel0(omega*r/v2)/(v2**2) - hankel1(omega*r/v2)/(omega*r*v2) + hankel1(omega*r/v1)/(omega*r*v1)) * dcmplx(0.d0,+pi/2.d0)
  endif

  end

! -----------

  double complex function hankel0(z)

  implicit none

  double complex z

! on utilise la routine NAG appelee S17DLE (simple precision)

  integer ifail,nz
  complex result

  ifail = -1
  call S17DLE(2,0.0,cmplx(z),1,'U',result,nz,ifail)
  if (ifail /= 0) stop 'S17DLE failed in hankel0'
  if (nz > 0) print *,nz,' termes mis a zero par underflow'

  hankel0 = dcmplx(result)

  end

! -----------

  double complex function hankel1(z)

  implicit none

  double complex z

! on utilise la routine NAG appelee S17DLE (simple precision)

  integer ifail,nz
  complex result

  ifail = -1
  call S17DLE(2,1.0,cmplx(z),1,'U',result,nz,ifail)
  if (ifail /= 0) stop 'S17DLE failed in hankel1'
  if (nz > 0) print *,nz,' termes mis a zero par underflow'

  hankel1 = dcmplx(result)

  end

! ***************** routine de FFT pour signal en temps ****************

! FFT routine taken from Netlib

  subroutine CFFTB (N,C,WSAVE)
  DIMENSION       C(1)       ,WSAVE(1)
  if (N == 1) return
  IW1 = N+N+1
  IW2 = IW1+N+N
  CALL CFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
  return
  END
  subroutine CFFTB1 (N,C,CH,WA,IFAC)
  DIMENSION       CH(1)      ,C(1)       ,WA(1)      ,IFAC(1)
  NF = IFAC(2)
  NA = 0
  L1 = 1
  IW = 1
  DO 116 K1=1,NF
   IP = IFAC(K1+2)
   L2 = IP*L1
   IDO = N/L2
   IDOT = IDO+IDO
   IDL1 = IDOT*L1
   if (IP /= 4) goto 103
   IX2 = IW+IDOT
   IX3 = IX2+IDOT
   if (NA /= 0) goto 101
   CALL PASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
   goto 102
  101    CALL PASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
   goto 115
  103    if (IP /= 2) goto 106
   if (NA /= 0) goto 104
   CALL PASSB2 (IDOT,L1,C,CH,WA(IW))
   goto 105
  104    CALL PASSB2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
   goto 115
  106    if (IP /= 3) goto 109
   IX2 = IW+IDOT
   if (NA /= 0) goto 107
   CALL PASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
   goto 108
  107    CALL PASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
   goto 115
  109    if (IP /= 5) goto 112
   IX2 = IW+IDOT
   IX3 = IX2+IDOT
   IX4 = IX3+IDOT
   if (NA /= 0) goto 110
   CALL PASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
   goto 111
  110    CALL PASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
   goto 115
  112    if (NA /= 0) goto 113
   CALL PASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
   goto 114
  113    CALL PASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    if (NAC /= 0) NA = 1-NA
  115    L1 = L2
   IW = IW+(IP-1)*IDOT
  116 continue
  if (NA == 0) return
  N2 = N+N
  DO 117 I=1,N2
   C(I) = CH(I)
  117 continue
  return
  END
  subroutine PASSB (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
  DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1), &
                  C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP), &
                  CH2(IDL1,IP)
  IDOT = IDO/2
  NT = IP*IDL1
  IPP2 = IP+2
  IPPH = (IP+1)/2
  IDP = IP*IDO
!
  if (IDO < L1) goto 106
  DO 103 J=2,IPPH
   JC = IPP2-J
   DO 102 K=1,L1
      DO 101 I=1,IDO
         CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
         CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       continue
  102    continue
  103 continue
  DO 105 K=1,L1
   DO 104 I=1,IDO
      CH(I,K,1) = CC(I,1,K)
  104    continue
  105 continue
  goto 112
  106 DO 109 J=2,IPPH
   JC = IPP2-J
   DO 108 I=1,IDO
      DO 107 K=1,L1
         CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
         CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       continue
  108    continue
  109 continue
  DO 111 I=1,IDO
   DO 110 K=1,L1
      CH(I,K,1) = CC(I,1,K)
  110    continue
  111 continue
  112 IDL = 2-IDO
  INC = 0
  DO 116 L=2,IPPH
   LC = IPP2-L
   IDL = IDL+IDO
   DO 113 IK=1,IDL1
      C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
      C2(IK,LC) = WA(IDL)*CH2(IK,IP)
  113    continue
   IDLJ = IDL
   INC = INC+IDO
   DO 115 J=3,IPPH
      JC = IPP2-J
      IDLJ = IDLJ+INC
      if (IDLJ > IDP) IDLJ = IDLJ-IDP
      WAR = WA(IDLJ-1)
      WAI = WA(IDLJ)
      DO 114 IK=1,IDL1
         C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
         C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
  114       continue
  115    continue
  116 continue
  DO 118 J=2,IPPH
   DO 117 IK=1,IDL1
      CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    continue
  118 continue
  DO 120 J=2,IPPH
   JC = IPP2-J
   DO 119 IK=2,IDL1,2
      CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
      CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
      CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
      CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    continue
  120 continue
  NAC = 1
  if (IDO == 2) return
  NAC = 0
  DO 121 IK=1,IDL1
   C2(IK,1) = CH2(IK,1)
  121 continue
  DO 123 J=2,IP
   DO 122 K=1,L1
      C1(1,K,J) = CH(1,K,J)
      C1(2,K,J) = CH(2,K,J)
  122    continue
  123 continue
  if (IDOT > L1) goto 127
  IDIJ = 0
  DO 126 J=2,IP
   IDIJ = IDIJ+2
   DO 125 I=4,IDO,2
      IDIJ = IDIJ+2
      DO 124 K=1,L1
         C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
         C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  124       continue
  125    continue
  126 continue
  return
  127 IDJ = 2-IDO
  DO 130 J=2,IP
   IDJ = IDJ+IDO
   DO 129 K=1,L1
      IDIJ = IDJ
      DO 128 I=4,IDO,2
         IDIJ = IDIJ+2
         C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
         C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  128       continue
  129    continue
  130 continue
  return
  END
  subroutine PASSB2 (IDO,L1,CC,CH,WA1)
  DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2), &
                  WA1(1)
  if (IDO > 2) goto 102
  DO 101 K=1,L1
   CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
   CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
   CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
   CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 continue
  return
  102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
      TR2 = CC(I-1,1,K)-CC(I-1,2,K)
      CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
      TI2 = CC(I,1,K)-CC(I,2,K)
      CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
      CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  103    continue
  104 continue
  return
  END
  subroutine PASSB3 (IDO,L1,CC,CH,WA1,WA2)
  DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3), &
                  WA1(1)     ,WA2(1)
  DATA TAUR,TAUI /-.5,.866025403784439/
  if (IDO /= 2) goto 102
  DO 101 K=1,L1
   TR2 = CC(1,2,K)+CC(1,3,K)
   CR2 = CC(1,1,K)+TAUR*TR2
   CH(1,K,1) = CC(1,1,K)+TR2
   TI2 = CC(2,2,K)+CC(2,3,K)
   CI2 = CC(2,1,K)+TAUR*TI2
   CH(2,K,1) = CC(2,1,K)+TI2
   CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
   CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
   CH(1,K,2) = CR2-CI3
   CH(1,K,3) = CR2+CI3
   CH(2,K,2) = CI2+CR3
   CH(2,K,3) = CI2-CR3
  101 continue
  return
  102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      TR2 = CC(I-1,2,K)+CC(I-1,3,K)
      CR2 = CC(I-1,1,K)+TAUR*TR2
      CH(I-1,K,1) = CC(I-1,1,K)+TR2
      TI2 = CC(I,2,K)+CC(I,3,K)
      CI2 = CC(I,1,K)+TAUR*TI2
      CH(I,K,1) = CC(I,1,K)+TI2
      CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
      CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
      DR2 = CR2-CI3
      DR3 = CR2+CI3
      DI2 = CI2+CR3
      DI3 = CI2-CR3
      CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
      CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
      CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
      CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
  103    continue
  104 continue
  return
  END
  subroutine PASSB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
  DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4), &
                  WA1(1)     ,WA2(1)     ,WA3(1)
  if (IDO /= 2) goto 102
  DO 101 K=1,L1
   TI1 = CC(2,1,K)-CC(2,3,K)
   TI2 = CC(2,1,K)+CC(2,3,K)
   TR4 = CC(2,4,K)-CC(2,2,K)
   TI3 = CC(2,2,K)+CC(2,4,K)
   TR1 = CC(1,1,K)-CC(1,3,K)
   TR2 = CC(1,1,K)+CC(1,3,K)
   TI4 = CC(1,2,K)-CC(1,4,K)
   TR3 = CC(1,2,K)+CC(1,4,K)
   CH(1,K,1) = TR2+TR3
   CH(1,K,3) = TR2-TR3
   CH(2,K,1) = TI2+TI3
   CH(2,K,3) = TI2-TI3
   CH(1,K,2) = TR1+TR4
   CH(1,K,4) = TR1-TR4
   CH(2,K,2) = TI1+TI4
   CH(2,K,4) = TI1-TI4
  101 continue
  return
  102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      TI1 = CC(I,1,K)-CC(I,3,K)
      TI2 = CC(I,1,K)+CC(I,3,K)
      TI3 = CC(I,2,K)+CC(I,4,K)
      TR4 = CC(I,4,K)-CC(I,2,K)
      TR1 = CC(I-1,1,K)-CC(I-1,3,K)
      TR2 = CC(I-1,1,K)+CC(I-1,3,K)
      TI4 = CC(I-1,2,K)-CC(I-1,4,K)
      TR3 = CC(I-1,2,K)+CC(I-1,4,K)
      CH(I-1,K,1) = TR2+TR3
      CR3 = TR2-TR3
      CH(I,K,1) = TI2+TI3
      CI3 = TI2-TI3
      CR2 = TR1+TR4
      CR4 = TR1-TR4
      CI2 = TI1+TI4
      CI4 = TI1-TI4
      CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
      CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
      CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
      CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
      CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
      CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  103    continue
  104 continue
  return
  END
  subroutine PASSB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
  DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5), &
                  WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)
  DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154, &
  -.809016994374947,.587785252292473/
  if (IDO /= 2) goto 102
  DO 101 K=1,L1
   TI5 = CC(2,2,K)-CC(2,5,K)
   TI2 = CC(2,2,K)+CC(2,5,K)
   TI4 = CC(2,3,K)-CC(2,4,K)
   TI3 = CC(2,3,K)+CC(2,4,K)
   TR5 = CC(1,2,K)-CC(1,5,K)
   TR2 = CC(1,2,K)+CC(1,5,K)
   TR4 = CC(1,3,K)-CC(1,4,K)
   TR3 = CC(1,3,K)+CC(1,4,K)
   CH(1,K,1) = CC(1,1,K)+TR2+TR3
   CH(2,K,1) = CC(2,1,K)+TI2+TI3
   CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
   CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
   CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
   CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
   CR5 = TI11*TR5+TI12*TR4
   CI5 = TI11*TI5+TI12*TI4
   CR4 = TI12*TR5-TI11*TR4
   CI4 = TI12*TI5-TI11*TI4
   CH(1,K,2) = CR2-CI5
   CH(1,K,5) = CR2+CI5
   CH(2,K,2) = CI2+CR5
   CH(2,K,3) = CI3+CR4
   CH(1,K,3) = CR3-CI4
   CH(1,K,4) = CR3+CI4
   CH(2,K,4) = CI3-CR4
   CH(2,K,5) = CI2-CR5
  101 continue
  return
  102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      TI5 = CC(I,2,K)-CC(I,5,K)
      TI2 = CC(I,2,K)+CC(I,5,K)
      TI4 = CC(I,3,K)-CC(I,4,K)
      TI3 = CC(I,3,K)+CC(I,4,K)
      TR5 = CC(I-1,2,K)-CC(I-1,5,K)
      TR2 = CC(I-1,2,K)+CC(I-1,5,K)
      TR4 = CC(I-1,3,K)-CC(I-1,4,K)
      TR3 = CC(I-1,3,K)+CC(I-1,4,K)
      CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
      CH(I,K,1) = CC(I,1,K)+TI2+TI3
      CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
      CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
      CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
      CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
      CR5 = TI11*TR5+TI12*TR4
      CI5 = TI11*TI5+TI12*TI4
      CR4 = TI12*TR5-TI11*TR4
      CI4 = TI12*TI5-TI11*TI4
      DR3 = CR3-CI4
      DR4 = CR3+CI4
      DI3 = CI3+CR4
      DI4 = CI3-CR4
      DR5 = CR2+CI5
      DR2 = CR2-CI5
      DI5 = CI2-CR5
      DI2 = CI2+CR5
      CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
      CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
      CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
      CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
      CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
      CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
      CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
      CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  103    continue
  104 continue
  return
  END



  subroutine CFFTI (N,WSAVE)
  DIMENSION       WSAVE(1)
  if (N == 1) return
  IW1 = N+N+1
  IW2 = IW1+N+N
  CALL CFFTI1 (N,WSAVE(IW1),WSAVE(IW2))
  return
  END
  subroutine CFFTI1 (N,WA,IFAC)
  DIMENSION       WA(1)      ,IFAC(1)    ,NTRYH(4)
  DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
  NL = N
  NF = 0
  J = 0
  101 J = J+1
  if (J-4) 102,102,103
  102 NTRY = NTRYH(J)
  goto 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
  NR = NL-NTRY*NQ
  if (NR) 101,105,101
  105 NF = NF+1
  IFAC(NF+2) = NTRY
  NL = NQ
  if (NTRY /= 2) goto 107
  if (NF == 1) goto 107
  DO 106 I=2,NF
   IB = NF-I+2
   IFAC(IB+2) = IFAC(IB+1)
  106 continue
  IFAC(3) = 2
  107 if (NL /= 1) goto 104
  IFAC(1) = N
  IFAC(2) = NF
  TPI = 6.28318530717959
  ARGH = TPI/FLOAT(N)
  I = 2
  L1 = 1
  DO 110 K1=1,NF
   IP = IFAC(K1+2)
   LD = 0
   L2 = L1*IP
   IDO = N/L2
   IDOT = IDO+IDO+2
   IPM = IP-1
   DO 109 J=1,IPM
      I1 = I
      WA(I-1) = 1.
      WA(I) = 0.
      LD = LD+L1
      FI = 0.
      ARGLD = FLOAT(LD)*ARGH
      DO 108 II=4,IDOT,2
         I = I+2
         FI = FI+1.
         ARG = FI*ARGLD
         WA(I-1) = COS(ARG)
         WA(I) = SIN(ARG)
  108       continue
      if (IP <= 5) goto 109
      WA(I1-1) = WA(I-1)
      WA(I1) = WA(I)
  109    continue
   L1 = L2
  110 continue
  return
  END





  subroutine CFFTF (N,C,WSAVE)
  DIMENSION       C(1)       ,WSAVE(1)
  if (N == 1) return
  IW1 = N+N+1
  IW2 = IW1+N+N
  CALL CFFTF1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
  return
  END
  subroutine CFFTF1 (N,C,CH,WA,IFAC)
  DIMENSION       CH(1)      ,C(1)       ,WA(1)      ,IFAC(1)
  NF = IFAC(2)
  NA = 0
  L1 = 1
  IW = 1
  DO 116 K1=1,NF
   IP = IFAC(K1+2)
   L2 = IP*L1
   IDO = N/L2
   IDOT = IDO+IDO
   IDL1 = IDOT*L1
   if (IP /= 4) goto 103
   IX2 = IW+IDOT
   IX3 = IX2+IDOT
   if (NA /= 0) goto 101
   CALL PASSF4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
   goto 102
  101    CALL PASSF4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
   goto 115
  103    if (IP /= 2) goto 106
   if (NA /= 0) goto 104
   CALL PASSF2 (IDOT,L1,C,CH,WA(IW))
   goto 105
  104    CALL PASSF2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
   goto 115
  106    if (IP /= 3) goto 109
   IX2 = IW+IDOT
   if (NA /= 0) goto 107
   CALL PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
   goto 108
  107    CALL PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
   goto 115
  109    if (IP /= 5) goto 112
   IX2 = IW+IDOT
   IX3 = IX2+IDOT
   IX4 = IX3+IDOT
   if (NA /= 0) goto 110
   CALL PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
   goto 111
  110    CALL PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
   goto 115
  112    if (NA /= 0) goto 113
   CALL PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
   goto 114
  113    CALL PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    if (NAC /= 0) NA = 1-NA
  115    L1 = L2
   IW = IW+(IP-1)*IDOT
  116 continue
  if (NA == 0) return
  N2 = N+N
  DO 117 I=1,N2
   C(I) = CH(I)
  117 continue
  return
  END
  subroutine PASSF (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
  DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1), &
                  C1(IDO,L1,IP)          ,WA(1)      ,C2(IDL1,IP), &
                  CH2(IDL1,IP)
  IDOT = IDO/2
  NT = IP*IDL1
  IPP2 = IP+2
  IPPH = (IP+1)/2
  IDP = IP*IDO
!
  if (IDO < L1) goto 106
  DO 103 J=2,IPPH
   JC = IPP2-J
   DO 102 K=1,L1
      DO 101 I=1,IDO
         CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
         CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       continue
  102    continue
  103 continue
  DO 105 K=1,L1
   DO 104 I=1,IDO
      CH(I,K,1) = CC(I,1,K)
  104    continue
  105 continue
  goto 112
  106 DO 109 J=2,IPPH
   JC = IPP2-J
   DO 108 I=1,IDO
      DO 107 K=1,L1
         CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
         CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       continue
  108    continue
  109 continue
  DO 111 I=1,IDO
   DO 110 K=1,L1
      CH(I,K,1) = CC(I,1,K)
  110    continue
  111 continue
  112 IDL = 2-IDO
  INC = 0
  DO 116 L=2,IPPH
   LC = IPP2-L
   IDL = IDL+IDO
   DO 113 IK=1,IDL1
      C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
      C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
  113    continue
   IDLJ = IDL
   INC = INC+IDO
   DO 115 J=3,IPPH
      JC = IPP2-J
      IDLJ = IDLJ+INC
      if (IDLJ > IDP) IDLJ = IDLJ-IDP
      WAR = WA(IDLJ-1)
      WAI = WA(IDLJ)
      DO 114 IK=1,IDL1
         C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
         C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
  114       continue
  115    continue
  116 continue
  DO 118 J=2,IPPH
   DO 117 IK=1,IDL1
      CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    continue
  118 continue
  DO 120 J=2,IPPH
   JC = IPP2-J
   DO 119 IK=2,IDL1,2
      CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
      CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
      CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
      CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    continue
  120 continue
  NAC = 1
  if (IDO == 2) return
  NAC = 0
  DO 121 IK=1,IDL1
   C2(IK,1) = CH2(IK,1)
  121 continue
  DO 123 J=2,IP
   DO 122 K=1,L1
      C1(1,K,J) = CH(1,K,J)
      C1(2,K,J) = CH(2,K,J)
  122    continue
  123 continue
  if (IDOT > L1) goto 127
  IDIJ = 0
  DO 126 J=2,IP
   IDIJ = IDIJ+2
   DO 125 I=4,IDO,2
      IDIJ = IDIJ+2
      DO 124 K=1,L1
         C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
         C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  124       continue
  125    continue
  126 continue
  return
  127 IDJ = 2-IDO
  DO 130 J=2,IP
   IDJ = IDJ+IDO
   DO 129 K=1,L1
      IDIJ = IDJ
      DO 128 I=4,IDO,2
         IDIJ = IDIJ+2
         C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
         C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  128       continue
  129    continue
  130 continue
  return
  END
  subroutine PASSF2 (IDO,L1,CC,CH,WA1)
  DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2), &
                  WA1(1)
  if (IDO > 2) goto 102
  DO 101 K=1,L1
   CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
   CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
   CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
   CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 continue
  return
  102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
      TR2 = CC(I-1,1,K)-CC(I-1,2,K)
      CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
      TI2 = CC(I,1,K)-CC(I,2,K)
      CH(I,K,2) = WA1(I-1)*TI2-WA1(I)*TR2
      CH(I-1,K,2) = WA1(I-1)*TR2+WA1(I)*TI2
  103    continue
  104 continue
  return
  END
  subroutine PASSF3 (IDO,L1,CC,CH,WA1,WA2)
  DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3), &
                  WA1(1)     ,WA2(1)
  DATA TAUR,TAUI /-.5,-.866025403784439/
  if (IDO /= 2) goto 102
  DO 101 K=1,L1
   TR2 = CC(1,2,K)+CC(1,3,K)
   CR2 = CC(1,1,K)+TAUR*TR2
   CH(1,K,1) = CC(1,1,K)+TR2
   TI2 = CC(2,2,K)+CC(2,3,K)
   CI2 = CC(2,1,K)+TAUR*TI2
   CH(2,K,1) = CC(2,1,K)+TI2
   CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
   CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
   CH(1,K,2) = CR2-CI3
   CH(1,K,3) = CR2+CI3
   CH(2,K,2) = CI2+CR3
   CH(2,K,3) = CI2-CR3
  101 continue
  return
  102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      TR2 = CC(I-1,2,K)+CC(I-1,3,K)
      CR2 = CC(I-1,1,K)+TAUR*TR2
      CH(I-1,K,1) = CC(I-1,1,K)+TR2
      TI2 = CC(I,2,K)+CC(I,3,K)
      CI2 = CC(I,1,K)+TAUR*TI2
      CH(I,K,1) = CC(I,1,K)+TI2
      CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
      CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
      DR2 = CR2-CI3
      DR3 = CR2+CI3
      DI2 = CI2+CR3
      DI3 = CI2-CR3
      CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
      CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
      CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
      CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
  103    continue
  104 continue
  return
  END
  subroutine PASSF4 (IDO,L1,CC,CH,WA1,WA2,WA3)
  DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4), &
                  WA1(1)     ,WA2(1)     ,WA3(1)
  if (IDO /= 2) goto 102
  DO 101 K=1,L1
   TI1 = CC(2,1,K)-CC(2,3,K)
   TI2 = CC(2,1,K)+CC(2,3,K)
   TR4 = CC(2,2,K)-CC(2,4,K)
   TI3 = CC(2,2,K)+CC(2,4,K)
   TR1 = CC(1,1,K)-CC(1,3,K)
   TR2 = CC(1,1,K)+CC(1,3,K)
   TI4 = CC(1,4,K)-CC(1,2,K)
   TR3 = CC(1,2,K)+CC(1,4,K)
   CH(1,K,1) = TR2+TR3
   CH(1,K,3) = TR2-TR3
   CH(2,K,1) = TI2+TI3
   CH(2,K,3) = TI2-TI3
   CH(1,K,2) = TR1+TR4
   CH(1,K,4) = TR1-TR4
   CH(2,K,2) = TI1+TI4
   CH(2,K,4) = TI1-TI4
  101 continue
  return
  102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      TI1 = CC(I,1,K)-CC(I,3,K)
      TI2 = CC(I,1,K)+CC(I,3,K)
      TI3 = CC(I,2,K)+CC(I,4,K)
      TR4 = CC(I,2,K)-CC(I,4,K)
      TR1 = CC(I-1,1,K)-CC(I-1,3,K)
      TR2 = CC(I-1,1,K)+CC(I-1,3,K)
      TI4 = CC(I-1,4,K)-CC(I-1,2,K)
      TR3 = CC(I-1,2,K)+CC(I-1,4,K)
      CH(I-1,K,1) = TR2+TR3
      CR3 = TR2-TR3
      CH(I,K,1) = TI2+TI3
      CI3 = TI2-TI3
      CR2 = TR1+TR4
      CR4 = TR1-TR4
      CI2 = TI1+TI4
      CI4 = TI1-TI4
      CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
      CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
      CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
      CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
      CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
      CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
  103    continue
  104 continue
  return
  END
  subroutine PASSF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
  DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5), &
                  WA1(1)     ,WA2(1)     ,WA3(1)     ,WA4(1)
  DATA TR11,TI11,TR12,TI12 /.309016994374947,-.951056516295154, &
  -.809016994374947,-.587785252292473/
  if (IDO /= 2) goto 102
  DO 101 K=1,L1
   TI5 = CC(2,2,K)-CC(2,5,K)
   TI2 = CC(2,2,K)+CC(2,5,K)
   TI4 = CC(2,3,K)-CC(2,4,K)
   TI3 = CC(2,3,K)+CC(2,4,K)
   TR5 = CC(1,2,K)-CC(1,5,K)
   TR2 = CC(1,2,K)+CC(1,5,K)
   TR4 = CC(1,3,K)-CC(1,4,K)
   TR3 = CC(1,3,K)+CC(1,4,K)
   CH(1,K,1) = CC(1,1,K)+TR2+TR3
   CH(2,K,1) = CC(2,1,K)+TI2+TI3
   CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
   CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
   CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
   CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
   CR5 = TI11*TR5+TI12*TR4
   CI5 = TI11*TI5+TI12*TI4
   CR4 = TI12*TR5-TI11*TR4
   CI4 = TI12*TI5-TI11*TI4
   CH(1,K,2) = CR2-CI5
   CH(1,K,5) = CR2+CI5
   CH(2,K,2) = CI2+CR5
   CH(2,K,3) = CI3+CR4
   CH(1,K,3) = CR3-CI4
   CH(1,K,4) = CR3+CI4
   CH(2,K,4) = CI3-CR4
   CH(2,K,5) = CI2-CR5
  101 continue
  return
  102 DO 104 K=1,L1
   DO 103 I=2,IDO,2
      TI5 = CC(I,2,K)-CC(I,5,K)
      TI2 = CC(I,2,K)+CC(I,5,K)
      TI4 = CC(I,3,K)-CC(I,4,K)
      TI3 = CC(I,3,K)+CC(I,4,K)
      TR5 = CC(I-1,2,K)-CC(I-1,5,K)
      TR2 = CC(I-1,2,K)+CC(I-1,5,K)
      TR4 = CC(I-1,3,K)-CC(I-1,4,K)
      TR3 = CC(I-1,3,K)+CC(I-1,4,K)
      CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
      CH(I,K,1) = CC(I,1,K)+TI2+TI3
      CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
      CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
      CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
      CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
      CR5 = TI11*TR5+TI12*TR4
      CI5 = TI11*TI5+TI12*TI4
      CR4 = TI12*TR5-TI11*TR4
      CI4 = TI12*TI5-TI11*TI4
      DR3 = CR3-CI4
      DR4 = CR3+CI4
      DI3 = CI3+CR4
      DI4 = CI3-CR4
      DR5 = CR2+CI5
      DR2 = CR2-CI5
      DI5 = CI2-CR5
      DI2 = CI2+CR5
      CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
      CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
      CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
      CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
      CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
      CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
      CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
      CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
  103    continue
  104 continue
  return
  END

! !!!!!!!! DK DK NAG routines included below

! DK DK march99 : routines recuperees sur le Cray (simple precision)

  subroutine ABZP01
!     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
!
!     Terminates execution when a hard failure occurs.
!
!     ******************** IMPLEMENTATION NOTE ********************
!     The following STOP statement may be replaced by a call to an
!     implementation-dependent routine to display a message and/or
!     to abort the program.
!     *************************************************************
!     .. Executable Statements ..
  STOP
  END

  subroutine DCYS18(Z,FNU,KODE,MR,N,Y,NZ,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-785 (DEC 1989).
!
!     Original name: CUNK2
!
!     DCYS18 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
!     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
!     UNIFORM ASYMPTOTIC EXPANSIONS FOR H(KIND,FNU,ZN) AND J(FNU,ZN)
!     WHERE ZN IS IN THE RIGHT HALF PLANE, KIND=(3-MR)/2, MR=+1 OR
!     -1. HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT
!     HALF PLANE OR ZR=-Z IF Z IS IN THE LEFT HALF PLANE. MR INDIC-
!     ATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
!     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, TOL
  INTEGER           KODE, MR, N, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           AI, ARGD, ASUMD, BSUMD, C1, C2, CFN, CI, CK, &
                    CONE, CR1, CR2, CRSC, CS, CSCL, CSGN, CSPN, &
                    CZERO, DAI, PHID, RZ, S1, S2, ZB, ZETA1D, &
                    ZETA2D, ZN, ZR
  REAL              AARG, AIC, ANG, APHI, ASC, ASCLE, C2I, C2M, C2R, &
                    CAR, CPN, FMR, FN, FNF, HPI, PI, RS1, SAR, SGN, &
                    SPN, X, YY
  INTEGER           I, IB, IC, IDUM, IFLAG, IFN, IL, IN, INU, IPARD, &
                    IUF, J, K, KDFLG, KFLAG, KK, NAI, NDAI, NW
!     .. Local Arrays ..
  COMPLEX           ARG(2), ASUM(2), BSUM(2), CIP(4), CSR(3), &
                    CSS(3), CY(2), PHI(2), ZETA1(2), ZETA2(2)
  REAL              BRY(3)
!     .. External functions ..
  REAL              X02AME, X02ALE
  EXTERNAL          X02AME, X02ALE
!     .. External subroutines ..
  EXTERNAL          DEUS17, S17DGE, DGSS17, DGVS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, CONJG, COS, EXP, INT, LOG, &
                    MAX, MOD, REAL, SIGN, SIN
!     .. Data statements ..
  DATA              CZERO, CONE, CI, CR1, CR2/(0.0E0,0.0E0), &
                    (1.0E0,0.0E0), (0.0E0,1.0E0), &
                    (1.0E0,1.73205080756887729E0), &
                    (-0.5E0,-8.66025403784438647E-01)/
  DATA              HPI, PI, AIC/1.57079632679489662E+00, &
                    3.14159265358979324E+00, &
                    1.26551212348464539E+00/
  DATA              CIP(1), CIP(2), CIP(3), CIP(4)/(1.0E0,0.0E0), &
                    (0.0E0,-1.0E0), (-1.0E0,0.0E0), (0.0E0,1.0E0)/
!     .. Executable Statements ..
!
  KDFLG = 1
  NZ = 0
!     ------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
!     THE UNDERFLOW LIMIT
!     ------------------------------------------------------------------
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CRSC = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CRSC
  CSR(1) = CRSC
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = (1.0E+3*X02AME())/TOL
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = X02ALE()
  X = REAL(Z)
  ZR = Z
  if (X < 0.0E0) ZR = -Z
  YY = AIMAG(ZR)
  ZN = -ZR*CI
  ZB = ZR
  INU = INT(FNU)
  FNF = FNU - INU
  ANG = -HPI*FNF
  CAR = COS(ANG)
  SAR = SIN(ANG)
  CPN = -HPI*CAR
  SPN = -HPI*SAR
  C2 = CMPLX(-SPN,CPN)
  KK = MOD(INU,4) + 1
  CS = CR1*C2*CIP(KK)
  if (YY <= 0.0E0) then
   ZN = CONJG(-ZN)
   ZB = CONJG(ZB)
  endif
!     ------------------------------------------------------------------
!     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
!     QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE COMPUTED BY
!     CONJUGATION SINCE THE K function IS REAL ON THE POSITIVE REAL AXIS
!     ------------------------------------------------------------------
  J = 2
  DO 40 I = 1, N
!        ---------------------------------------------------------------
!        J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
!        ---------------------------------------------------------------
   J = 3 - J
   FN = FNU + I - 1
   CALL DEUS17(ZN,FN,0,TOL,PHI(J),ARG(J),ZETA1(J),ZETA2(J),ASUM(J) &
                 ,BSUM(J),ELIM)
   if (KODE == 1) then
      S1 = ZETA1(J) - ZETA2(J)
   ELSE
      CFN = CMPLX(FN,0.0E0)
      S1 = ZETA1(J) - CFN*(CFN/(ZB+ZETA2(J)))
   endif
!        ---------------------------------------------------------------
!        TEST FOR UNDERFLOW AND OVERFLOW
!        ---------------------------------------------------------------
   RS1 = REAL(S1)
   if (ABS(RS1) <= ELIM) then
      if (KDFLG == 1) KFLAG = 2
      if (ABS(RS1) >= ALIM) then
!              ---------------------------------------------------------
!              REFINE  TEST AND SCALE
!              ---------------------------------------------------------
         APHI = ABS(PHI(J))
         AARG = ABS(ARG(J))
         RS1 = RS1 + LOG(APHI) - 0.25E0*LOG(AARG) - AIC
         if (ABS(RS1) > ELIM) then
            goto 20
         ELSE
            if (KDFLG == 1) KFLAG = 1
            if (RS1 >= 0.0E0) then
               if (KDFLG == 1) KFLAG = 3
            endif
         endif
      endif
!           ------------------------------------------------------------
!           SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!           EXPONENT EXTREMES
!           ------------------------------------------------------------
      C2 = ARG(J)*CR2
      IDUM = 1
!           S17DGE assumed not to fail, therefore IDUM set to one.
      CALL S17DGE('F',C2,'S',AI,NAI,IDUM)
      IDUM = 1
      CALL S17DGE('D',C2,'S',DAI,NDAI,IDUM)
      S2 = CS*PHI(J)*(AI*ASUM(J)+CR2*DAI*BSUM(J))
      C2R = REAL(S1)
      C2I = AIMAG(S1)
      C2M = EXP(C2R)*REAL(CSS(KFLAG))
      S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
      S2 = S2*S1
      if (KFLAG == 1) then
         CALL DGVS17(S2,NW,BRY(1),TOL)
         if (NW /= 0) goto 20
      endif
      if (YY <= 0.0E0) S2 = CONJG(S2)
      CY(KDFLG) = S2
      Y(I) = S2*CSR(KFLAG)
      CS = -CI*CS
      if (KDFLG == 2) then
         goto 60
      ELSE
         KDFLG = 2
         goto 40
      endif
   endif
   20    if (RS1 > 0.0E0) then
      goto 280
!           ------------------------------------------------------------
!           FOR X < 0.0, THE I function TO BE ADDED WILL OVERFLOW
!           ------------------------------------------------------------
   else if (X < 0.0E0) then
      goto 280
   ELSE
      KDFLG = 1
      Y(I) = CZERO
      CS = -CI*CS
      NZ = NZ + 1
      if (I /= 1) then
         if (Y(I-1) /= CZERO) then
            Y(I-1) = CZERO
            NZ = NZ + 1
         endif
      endif
   endif
   40 continue
  I = N
   60 RZ = CMPLX(2.0E0,0.0E0)/ZR
  CK = CMPLX(FN,0.0E0)*RZ
  IB = I + 1
  if (N >= IB) then
!        ---------------------------------------------------------------
!        TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO
!        ZERO ON UNDERFLOW
!        ---------------------------------------------------------------
   FN = FNU + N - 1
   IPARD = 1
   if (MR /= 0) IPARD = 0
   CALL DEUS17(ZN,FN,IPARD,TOL,PHID,ARGD,ZETA1D,ZETA2D,ASUMD, &
                 BSUMD,ELIM)
   if (KODE == 1) then
      S1 = ZETA1D - ZETA2D
   ELSE
      CFN = CMPLX(FN,0.0E0)
      S1 = ZETA1D - CFN*(CFN/(ZB+ZETA2D))
   endif
   RS1 = REAL(S1)
   if (ABS(RS1) <= ELIM) then
      if (ABS(RS1) >= ALIM) then
!              ---------------------------------------------------------
!              REFINE ESTIMATE AND TEST
!              ---------------------------------------------------------
         APHI = ABS(PHID)
         AARG = ABS(ARGD)
         RS1 = RS1 + LOG(APHI) - 0.25E0*LOG(AARG) - AIC
         if (ABS(RS1) >= ELIM) goto 100
      endif
!           ------------------------------------------------------------
!           SCALED FORWARD RECURRENCE FOR REMAINDER OF THE SEQUENCE
!           ------------------------------------------------------------
      S1 = CY(1)
      S2 = CY(2)
      C1 = CSR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 80 I = IB, N
         C2 = S2
         S2 = CK*S2 + S1
         S1 = C2
         CK = CK + RZ
         C2 = S2*C1
         Y(I) = C2
         if (KFLAG < 3) then
            C2R = REAL(C2)
            C2I = AIMAG(C2)
            C2R = ABS(C2R)
            C2I = ABS(C2I)
            C2M = MAX(C2R,C2I)
            if (C2M > ASCLE) then
               KFLAG = KFLAG + 1
               ASCLE = BRY(KFLAG)
               S1 = S1*C1
               S2 = C2
               S1 = S1*CSS(KFLAG)
               S2 = S2*CSS(KFLAG)
               C1 = CSR(KFLAG)
            endif
         endif
   80       continue
      goto 140
   endif
  100    if (RS1 > 0.0E0) then
      goto 280
!           ------------------------------------------------------------
!           FOR X < 0.0, THE I function TO BE ADDED WILL OVERFLOW
!           ------------------------------------------------------------
   else if (X < 0.0E0) then
      goto 280
   ELSE
      NZ = N
      DO 120 I = 1, N
         Y(I) = CZERO
  120       continue
      return
   endif
  endif
  140 if (MR == 0) then
   return
  ELSE
!        ---------------------------------------------------------------
!        ANALYTIC CONTINUATION FOR RE(Z) < 0.0E0
!        ---------------------------------------------------------------
   NZ = 0
   FMR = MR
   SGN = -SIGN(PI,FMR)
!        ---------------------------------------------------------------
!        CSPN AND CSGN ARE COEFF OF K AND I functionS RESP.
!        ---------------------------------------------------------------
   CSGN = CMPLX(0.0E0,SGN)
   if (YY <= 0.0E0) CSGN = CONJG(CSGN)
   IFN = INU + N - 1
   ANG = FNF*SGN
   CPN = COS(ANG)
   SPN = SIN(ANG)
   CSPN = CMPLX(CPN,SPN)
   if (MOD(IFN,2) == 1) CSPN = -CSPN
!        ---------------------------------------------------------------
!        CS=COEFF OF THE J function TO GET THE I function. I(FNU,Z) IS
!        COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE
!        FIRST QUADRANT. FOURTH QUADRANT VALUES (YY <= 0.0E0) ARE
!        COMPUTED BY CONJUGATION SINCE THE I function IS REAL ON THE
!        POSITIVE REAL AXIS
!        ---------------------------------------------------------------
   CS = CMPLX(CAR,-SAR)*CSGN
   IN = MOD(IFN,4) + 1
   C2 = CIP(IN)
   CS = CS*CONJG(C2)
   ASC = BRY(1)
   KK = N
   KDFLG = 1
   IB = IB - 1
   IC = IB - 1
   IUF = 0
   DO 220 K = 1, N
!           ------------------------------------------------------------
!           LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
!           function ABOVE
!           ------------------------------------------------------------
      FN = FNU + KK - 1
      if (N > 2) then
         if ((KK == N) .and. (IB < N)) then
            goto 160
         else if ((KK /= IB) .and. (KK /= IC)) then
            CALL DEUS17(ZN,FN,0,TOL,PHID,ARGD,ZETA1D,ZETA2D,ASUMD, &
                          BSUMD,ELIM)
            goto 160
         endif
      endif
      PHID = PHI(J)
      ARGD = ARG(J)
      ZETA1D = ZETA1(J)
      ZETA2D = ZETA2(J)
      ASUMD = ASUM(J)
      BSUMD = BSUM(J)
      J = 3 - J
  160       if (KODE == 1) then
         S1 = -ZETA1D + ZETA2D
      ELSE
         CFN = CMPLX(FN,0.0E0)
         S1 = -ZETA1D + CFN*(CFN/(ZB+ZETA2D))
      endif
!           ------------------------------------------------------------
!           TEST FOR UNDERFLOW AND OVERFLOW
!           ------------------------------------------------------------
      RS1 = REAL(S1)
      if (ABS(RS1) <= ELIM) then
         if (KDFLG == 1) IFLAG = 2
         if (ABS(RS1) >= ALIM) then
!                 ------------------------------------------------------
!                 REFINE  TEST AND SCALE
!                 ------------------------------------------------------
            APHI = ABS(PHID)
            AARG = ABS(ARGD)
            RS1 = RS1 + LOG(APHI) - 0.25E0*LOG(AARG) - AIC
            if (ABS(RS1) > ELIM) then
               goto 180
            ELSE
               if (KDFLG == 1) IFLAG = 1
               if (RS1 >= 0.0E0) then
                  if (KDFLG == 1) IFLAG = 3
               endif
            endif
         endif
         IDUM = 1
!              S17DGE assumed not to fail, therefore IDUM set to one.
         CALL S17DGE('F',ARGD,'S',AI,NAI,IDUM)
         IDUM = 1
         CALL S17DGE('D',ARGD,'S',DAI,NDAI,IDUM)
         S2 = CS*PHID*(AI*ASUMD+DAI*BSUMD)
         C2R = REAL(S1)
         C2I = AIMAG(S1)
         C2M = EXP(C2R)*REAL(CSS(IFLAG))
         S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
         S2 = S2*S1
         if (IFLAG == 1) then
            CALL DGVS17(S2,NW,BRY(1),TOL)
            if (NW /= 0) S2 = CMPLX(0.0E0,0.0E0)
         endif
         goto 200
      endif
  180       if (RS1 > 0.0E0) then
         goto 280
      ELSE
         S2 = CZERO
      endif
  200       if (YY <= 0.0E0) S2 = CONJG(S2)
      CY(KDFLG) = S2
      C2 = S2
      S2 = S2*CSR(IFLAG)
!           ------------------------------------------------------------
!           ADD I AND K functionS, K SEQUENCE IN Y(I), I=1,N
!           ------------------------------------------------------------
      S1 = Y(KK)
      if (KODE /= 1) then
         CALL DGSS17(ZR,S1,S2,NW,ASC,ALIM,IUF)
         NZ = NZ + NW
      endif
      Y(KK) = S1*CSPN + S2
      KK = KK - 1
      CSPN = -CSPN
      CS = -CS*CI
      if (C2 == CZERO) then
         KDFLG = 1
      else if (KDFLG == 2) then
         goto 240
      ELSE
         KDFLG = 2
      endif
  220    continue
   K = N
  240    IL = N - K
   if (IL /= 0) then
!           ------------------------------------------------------------
!           RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
!           K functionS, SCALING THE I SEQUENCE DURING RECURRENCE TO
!           KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT
!           EXTREMES.
!           ------------------------------------------------------------
      S1 = CY(1)
      S2 = CY(2)
      CS = CSR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = INU + IL
      DO 260 I = 1, IL
         C2 = S2
         S2 = S1 + CMPLX(FN+FNF,0.0E0)*RZ*S2
         S1 = C2
         FN = FN - 1.0E0
         C2 = S2*CS
         CK = C2
         C1 = Y(KK)
         if (KODE /= 1) then
            CALL DGSS17(ZR,C1,C2,NW,ASC,ALIM,IUF)
            NZ = NZ + NW
         endif
         Y(KK) = C1*CSPN + C2
         KK = KK - 1
         CSPN = -CSPN
         if (IFLAG < 3) then
            C2R = REAL(CK)
            C2I = AIMAG(CK)
            C2R = ABS(C2R)
            C2I = ABS(C2I)
            C2M = MAX(C2R,C2I)
            if (C2M > ASCLE) then
               IFLAG = IFLAG + 1
               ASCLE = BRY(IFLAG)
               S1 = S1*CS
               S2 = CK
               S1 = S1*CSS(IFLAG)
               S2 = S2*CSS(IFLAG)
               CS = CSR(IFLAG)
            endif
         endif
  260       continue
   endif
   return
  endif
  280 NZ = -1
  return
  END
  subroutine DCZS18(Z,FNU,KODE,MR,N,Y,NZ,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-786 (DEC 1989).
!
!     Original name: CUNK1
!
!     DCZS18 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
!     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
!     UNIFORM ASYMPTOTIC EXPANSION.
!     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
!     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, TOL
  INTEGER           KODE, MR, N, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           C1, C2, CFN, CK, CONE, CRSC, CS, CSCL, CSGN, &
                    CSPN, CZERO, PHID, RZ, S1, S2, SUMD, ZETA1D, &
                    ZETA2D, ZR
  REAL              ANG, APHI, ASC, ASCLE, C2I, C2M, C2R, CPN, FMR, &
                    FN, FNF, PI, RS1, SGN, SPN, X
  INTEGER           I, IB, IC, IFLAG, IFN, IL, INITD, INU, IPARD, &
                    IUF, J, K, KDFLG, KFLAG, KK, M, NW
!     .. Local Arrays ..
  COMPLEX           CSR(3), CSS(3), CWRK(16,3), CY(2), PHI(2), &
                    SUM(2), ZETA1(2), ZETA2(2)
  REAL              BRY(3)
  INTEGER           INIT(2)
!     .. External functions ..
  REAL              X02AME, X02ALE
  EXTERNAL          X02AME, X02ALE
!     .. External subroutines ..
  EXTERNAL          DEWS17, DGSS17, DGVS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, COS, EXP, INT, LOG, MAX, MOD, &
                    REAL, SIGN, SIN
!     .. Data statements ..
  DATA              CZERO, CONE/(0.0E0,0.0E0), (1.0E0,0.0E0)/
  DATA              PI/3.14159265358979324E0/
!     .. Executable Statements ..
!
  KDFLG = 1
  NZ = 0
!     ------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
!     THE UNDERFLOW LIMIT
!     ------------------------------------------------------------------
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CRSC = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CRSC
  CSR(1) = CRSC
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = (1.0E+3*X02AME())/TOL
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = X02ALE()
  X = REAL(Z)
  ZR = Z
  if (X < 0.0E0) ZR = -Z
  J = 2
  DO 40 I = 1, N
!        ---------------------------------------------------------------
!        J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
!        ---------------------------------------------------------------
   J = 3 - J
   FN = FNU + I - 1
   INIT(J) = 0
   CALL DEWS17(ZR,FN,2,0,TOL,INIT(J),PHI(J),ZETA1(J),ZETA2(J), &
                 SUM(J),CWRK(1,J),ELIM)
   if (KODE == 1) then
      S1 = ZETA1(J) - ZETA2(J)
   ELSE
      CFN = CMPLX(FN,0.0E0)
      S1 = ZETA1(J) - CFN*(CFN/(ZR+ZETA2(J)))
   endif
!        ---------------------------------------------------------------
!        TEST FOR UNDERFLOW AND OVERFLOW
!        ---------------------------------------------------------------
   RS1 = REAL(S1)
   if (ABS(RS1) <= ELIM) then
      if (KDFLG == 1) KFLAG = 2
      if (ABS(RS1) >= ALIM) then
!              ---------------------------------------------------------
!              REFINE  TEST AND SCALE
!              ---------------------------------------------------------
         APHI = ABS(PHI(J))
         RS1 = RS1 + LOG(APHI)
         if (ABS(RS1) > ELIM) then
            goto 20
         ELSE
            if (KDFLG == 1) KFLAG = 1
            if (RS1 >= 0.0E0) then
               if (KDFLG == 1) KFLAG = 3
            endif
         endif
      endif
!           ------------------------------------------------------------
!           SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!           EXPONENT EXTREMES
!           ------------------------------------------------------------
      S2 = PHI(J)*SUM(J)
      C2R = REAL(S1)
      C2I = AIMAG(S1)
      C2M = EXP(C2R)*REAL(CSS(KFLAG))
      S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
      S2 = S2*S1
      if (KFLAG == 1) then
         CALL DGVS17(S2,NW,BRY(1),TOL)
         if (NW /= 0) goto 20
      endif
      CY(KDFLG) = S2
      Y(I) = S2*CSR(KFLAG)
      if (KDFLG == 2) then
         goto 60
      ELSE
         KDFLG = 2
         goto 40
      endif
   endif
   20    if (RS1 > 0.0E0) then
      goto 280
!           ------------------------------------------------------------
!           FOR X < 0.0, THE I function TO BE ADDED WILL OVERFLOW
!           ------------------------------------------------------------
   else if (X < 0.0E0) then
      goto 280
   ELSE
      KDFLG = 1
      Y(I) = CZERO
      NZ = NZ + 1
      if (I /= 1) then
         if (Y(I-1) /= CZERO) then
            Y(I-1) = CZERO
            NZ = NZ + 1
         endif
      endif
   endif
   40 continue
  I = N
   60 RZ = CMPLX(2.0E0,0.0E0)/ZR
  CK = CMPLX(FN,0.0E0)*RZ
  IB = I + 1
  if (N >= IB) then
!        ---------------------------------------------------------------
!        TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO
!        ZERO ON UNDERFLOW
!        ---------------------------------------------------------------
   FN = FNU + N - 1
   IPARD = 1
   if (MR /= 0) IPARD = 0
   INITD = 0
   CALL DEWS17(ZR,FN,2,IPARD,TOL,INITD,PHID,ZETA1D,ZETA2D,SUMD, &
                 CWRK(1,3),ELIM)
   if (KODE == 1) then
      S1 = ZETA1D - ZETA2D
   ELSE
      CFN = CMPLX(FN,0.0E0)
      S1 = ZETA1D - CFN*(CFN/(ZR+ZETA2D))
   endif
   RS1 = REAL(S1)
   if (ABS(RS1) <= ELIM) then
      if (ABS(RS1) >= ALIM) then
!              ---------------------------------------------------------
!              REFINE ESTIMATE AND TEST
!              ---------------------------------------------------------
         APHI = ABS(PHID)
         RS1 = RS1 + LOG(APHI)
         if (ABS(RS1) >= ELIM) goto 100
      endif
!           ------------------------------------------------------------
!           RECUR FORWARD FOR REMAINDER OF THE SEQUENCE
!           ------------------------------------------------------------
      S1 = CY(1)
      S2 = CY(2)
      C1 = CSR(KFLAG)
      ASCLE = BRY(KFLAG)
      DO 80 I = IB, N
         C2 = S2
         S2 = CK*S2 + S1
         S1 = C2
         CK = CK + RZ
         C2 = S2*C1
         Y(I) = C2
         if (KFLAG < 3) then
            C2R = REAL(C2)
            C2I = AIMAG(C2)
            C2R = ABS(C2R)
            C2I = ABS(C2I)
            C2M = MAX(C2R,C2I)
            if (C2M > ASCLE) then
               KFLAG = KFLAG + 1
               ASCLE = BRY(KFLAG)
               S1 = S1*C1
               S2 = C2
               S1 = S1*CSS(KFLAG)
               S2 = S2*CSS(KFLAG)
               C1 = CSR(KFLAG)
            endif
         endif
   80       continue
      goto 140
   endif
  100    if (RS1 > 0.0E0) then
      goto 280
!           ------------------------------------------------------------
!           FOR X < 0.0, THE I function TO BE ADDED WILL OVERFLOW
!           ------------------------------------------------------------
   else if (X < 0.0E0) then
      goto 280
   ELSE
      NZ = N
      DO 120 I = 1, N
         Y(I) = CZERO
  120       continue
      return
   endif
  endif
  140 if (MR == 0) then
   return
  ELSE
!        ---------------------------------------------------------------
!        ANALYTIC CONTINUATION FOR RE(Z) < 0.0E0
!        ---------------------------------------------------------------
   NZ = 0
   FMR = MR
   SGN = -SIGN(PI,FMR)
!        ---------------------------------------------------------------
!        CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
!        ---------------------------------------------------------------
   CSGN = CMPLX(0.0E0,SGN)
   INU = INT(FNU)
   FNF = FNU - INU
   IFN = INU + N - 1
   ANG = FNF*SGN
   CPN = COS(ANG)
   SPN = SIN(ANG)
   CSPN = CMPLX(CPN,SPN)
   if (MOD(IFN,2) == 1) CSPN = -CSPN
   ASC = BRY(1)
   KK = N
   IUF = 0
   KDFLG = 1
   IB = IB - 1
   IC = IB - 1
   DO 220 K = 1, N
      FN = FNU + KK - 1
!           ------------------------------------------------------------
!           LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
!           function ABOVE
!           ------------------------------------------------------------
      M = 3
      if (N > 2) then
         if ((KK == N) .and. (IB < N)) then
            goto 160
         else if ((KK /= IB) .and. (KK /= IC)) then
            INITD = 0
            goto 160
         endif
      endif
      INITD = INIT(J)
      PHID = PHI(J)
      ZETA1D = ZETA1(J)
      ZETA2D = ZETA2(J)
      SUMD = SUM(J)
      M = J
      J = 3 - J
  160       CALL DEWS17(ZR,FN,1,0,TOL,INITD,PHID,ZETA1D,ZETA2D,SUMD, &
                    CWRK(1,M),ELIM)
      if (KODE == 1) then
         S1 = -ZETA1D + ZETA2D
      ELSE
         CFN = CMPLX(FN,0.0E0)
         S1 = -ZETA1D + CFN*(CFN/(ZR+ZETA2D))
      endif
!           ------------------------------------------------------------
!           TEST FOR UNDERFLOW AND OVERFLOW
!           ------------------------------------------------------------
      RS1 = REAL(S1)
      if (ABS(RS1) <= ELIM) then
         if (KDFLG == 1) IFLAG = 2
         if (ABS(RS1) >= ALIM) then
!                 ------------------------------------------------------
!                 REFINE  TEST AND SCALE
!                 ------------------------------------------------------
            APHI = ABS(PHID)
            RS1 = RS1 + LOG(APHI)
            if (ABS(RS1) > ELIM) then
               goto 180
            ELSE
               if (KDFLG == 1) IFLAG = 1
               if (RS1 >= 0.0E0) then
                  if (KDFLG == 1) IFLAG = 3
               endif
            endif
         endif
         S2 = CSGN*PHID*SUMD
         C2R = REAL(S1)
         C2I = AIMAG(S1)
         C2M = EXP(C2R)*REAL(CSS(IFLAG))
         S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
         S2 = S2*S1
         if (IFLAG == 1) then
            CALL DGVS17(S2,NW,BRY(1),TOL)
            if (NW /= 0) S2 = CMPLX(0.0E0,0.0E0)
         endif
         goto 200
      endif
  180       if (RS1 > 0.0E0) then
         goto 280
      ELSE
         S2 = CZERO
      endif
  200       CY(KDFLG) = S2
      C2 = S2
      S2 = S2*CSR(IFLAG)
!           ------------------------------------------------------------
!           ADD I AND K functionS, K SEQUENCE IN Y(I), I=1,N
!           ------------------------------------------------------------
      S1 = Y(KK)
      if (KODE /= 1) then
         CALL DGSS17(ZR,S1,S2,NW,ASC,ALIM,IUF)
         NZ = NZ + NW
      endif
      Y(KK) = S1*CSPN + S2
      KK = KK - 1
      CSPN = -CSPN
      if (C2 == CZERO) then
         KDFLG = 1
      else if (KDFLG == 2) then
         goto 240
      ELSE
         KDFLG = 2
      endif
  220    continue
   K = N
  240    IL = N - K
   if (IL /= 0) then
!           ------------------------------------------------------------
!           RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
!           K functionS, SCALING THE I SEQUENCE DURING RECURRENCE TO
!           KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT
!           EXTREMES.
!           ------------------------------------------------------------
      S1 = CY(1)
      S2 = CY(2)
      CS = CSR(IFLAG)
      ASCLE = BRY(IFLAG)
      FN = INU + IL
      DO 260 I = 1, IL
         C2 = S2
         S2 = S1 + CMPLX(FN+FNF,0.0E0)*RZ*S2
         S1 = C2
         FN = FN - 1.0E0
         C2 = S2*CS
         CK = C2
         C1 = Y(KK)
         if (KODE /= 1) then
            CALL DGSS17(ZR,C1,C2,NW,ASC,ALIM,IUF)
            NZ = NZ + NW
         endif
         Y(KK) = C1*CSPN + C2
         KK = KK - 1
         CSPN = -CSPN
         if (IFLAG < 3) then
            C2R = REAL(CK)
            C2I = AIMAG(CK)
            C2R = ABS(C2R)
            C2I = ABS(C2I)
            C2M = MAX(C2R,C2I)
            if (C2M > ASCLE) then
               IFLAG = IFLAG + 1
               ASCLE = BRY(IFLAG)
               S1 = S1*CS
               S2 = CK
               S1 = S1*CSS(IFLAG)
               S2 = S2*CSS(IFLAG)
               CS = CSR(IFLAG)
            endif
         endif
  260       continue
   endif
   return
  endif
  280 NZ = -1
  return
  END
  subroutine DERS17(Z,FNU,N,CY,TOL)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-761 (DEC 1989).
!
!     Original name: CRATI
!
!     DERS17 COMPUTES RATIOS OF I BESSEL functionS BY BACKWARD
!     RECURRENCE.  THE STARTING INDEX IS DETERMINED BY FORWARD
!     RECURRENCE AS DESCRIBED IN J. RES. OF NAT. BUR. OF STANDARDS-B,
!     MATHEMATICAL SCIENCES, VOL 77B, P111-114, SEPTEMBER, 1973,
!     BESSEL functionS I AND J OF COMPLEX ARGUMENT AND INTEGER ORDER,
!     BY D. J. SOOKNE.
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              FNU, TOL
  INTEGER           N
!     .. Array Arguments ..
  COMPLEX           CY(N)
!     .. Local Scalars ..
  COMPLEX           CDFNU, CONE, CZERO, P1, P2, PT, RZ, T1
  REAL              AK, AMAGZ, AP1, AP2, ARG, AZ, DFNU, FDNU, FLAM, &
                    FNUP, RAP1, RHO, TEST, TEST1
  INTEGER           I, ID, IDNU, INU, ITIME, K, KK, MAGZ
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, INT, MAX, MIN, REAL, SQRT
!     .. Data statements ..
  DATA              CZERO, CONE/(0.0E0,0.0E0), (1.0E0,0.0E0)/
!     .. Executable Statements ..
!
  AZ = ABS(Z)
  INU = INT(FNU)
  IDNU = INU + N - 1
  FDNU = IDNU
  MAGZ = INT(AZ)
  AMAGZ = MAGZ + 1
  FNUP = MAX(AMAGZ,FDNU)
  ID = IDNU - MAGZ - 1
  ITIME = 1
  K = 1
  RZ = (CONE+CONE)/Z
  T1 = CMPLX(FNUP,0.0E0)*RZ
  P2 = -T1
  P1 = CONE
  T1 = T1 + RZ
  if (ID > 0) ID = 0
  AP2 = ABS(P2)
  AP1 = ABS(P1)
!     ------------------------------------------------------------------
!     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNX
!     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT
!     P2 VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR
!     PREMATURELY.
!     ------------------------------------------------------------------
  ARG = (AP2+AP2)/(AP1*TOL)
  TEST1 = SQRT(ARG)
  TEST = TEST1
  RAP1 = 1.0E0/AP1
  P1 = P1*CMPLX(RAP1,0.0E0)
  P2 = P2*CMPLX(RAP1,0.0E0)
  AP2 = AP2*RAP1
   20 continue
  K = K + 1
  AP1 = AP2
  PT = P2
  P2 = P1 - T1*P2
  P1 = PT
  T1 = T1 + RZ
  AP2 = ABS(P2)
  if (AP1 <= TEST) then
   goto 20
  else if (ITIME /= 2) then
   AK = ABS(T1)*0.5E0
   FLAM = AK + SQRT(AK*AK-1.0E0)
   RHO = MIN(AP2/AP1,FLAM)
   TEST = TEST1*SQRT(RHO/(RHO*RHO-1.0E0))
   ITIME = 2
   goto 20
  endif
  KK = K + 1 - ID
  AK = KK
  DFNU = FNU + N - 1
  CDFNU = CMPLX(DFNU,0.0E0)
  T1 = CMPLX(AK,0.0E0)
  P1 = CMPLX(1.0E0/AP2,0.0E0)
  P2 = CZERO
  DO 40 I = 1, KK
   PT = P1
   P1 = RZ*(CDFNU+T1)*P1 + P2
   P2 = PT
   T1 = T1 - CONE
   40 continue
  if (REAL(P1) == 0.0E0 .and. AIMAG(P1) == 0.0E0) P1 = CMPLX(TOL, &
      TOL)
  CY(N) = P2/P1
  if (N /= 1) then
   K = N - 1
   AK = K
   T1 = CMPLX(AK,0.0E0)
   CDFNU = CMPLX(FNU,0.0E0)*RZ
   DO 60 I = 2, N
      PT = CDFNU + T1*RZ + CY(K+1)
      if (REAL(PT) == 0.0E0 .and. AIMAG(PT) == 0.0E0) &
            PT = CMPLX(TOL,TOL)
      CY(K) = CONE/PT
      T1 = T1 - CONE
      K = K - 1
   60    continue
  endif
  return
  END
  subroutine DESS17(ZR,FNU,KODE,N,Y,NZ,CW,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-762 (DEC 1989).
!
!     Original name: CWRSK
!
!     DESS17 COMPUTES THE I BESSEL function FOR RE(Z) >= 0.0 BY
!     NORMALIZING THE I function RATIOS FROM DERS17 BY THE WRONSKIAN
!
!     .. Scalar Arguments ..
  COMPLEX           ZR
  REAL              ALIM, ELIM, FNU, TOL
  INTEGER           KODE, N, NZ
!     .. Array Arguments ..
  COMPLEX           CW(2), Y(N)
!     .. Local Scalars ..
  COMPLEX           C1, C2, CINU, CSCL, CT, RCT, ST
  REAL              ACT, ACW, ASCLE, S1, S2, YY
  INTEGER           I, NW
!     .. External functions ..
  REAL              X02AME
  EXTERNAL          X02AME
!     .. External subroutines ..
  EXTERNAL          DERS17, DGXS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, CONJG, COS, SIN
!     .. Executable Statements ..
!     ------------------------------------------------------------------
!     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
!     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM DERS17 NORMALIZED BY THE
!     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM DGXS17.
!     ------------------------------------------------------------------
  NZ = 0
  CALL DGXS17(ZR,FNU,KODE,2,CW,NW,TOL,ELIM,ALIM)
  if (NW /= 0) then
   NZ = -1
   if (NW == (-2)) NZ = -2
   if (NW == (-3)) NZ = -3
  ELSE
   CALL DERS17(ZR,FNU,N,Y,TOL)
!        ---------------------------------------------------------------
!        RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
!        R(FNU+J-1,Z)=Y(J),  J=1,...,N
!        ---------------------------------------------------------------
   CINU = CMPLX(1.0E0,0.0E0)
   if (KODE /= 1) then
      YY = AIMAG(ZR)
      S1 = COS(YY)
      S2 = SIN(YY)
      CINU = CMPLX(S1,S2)
   endif
!        ---------------------------------------------------------------
!        ON LOW EXPONENT MACHINES THE K functionS CAN BE CLOSE TO BOTH
!        THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
!        SCALED TO PREVENT OVER OR UNDERFLOW. DEVS17 HAS DETERMINED THAT
!        THE RESULT IS ON SCALE.
!        ---------------------------------------------------------------
   ACW = ABS(CW(2))
   ASCLE = (1.0E+3*X02AME())/TOL
   CSCL = CMPLX(1.0E0,0.0E0)
   if (ACW > ASCLE) then
      ASCLE = 1.0E0/ASCLE
      if (ACW >= ASCLE) CSCL = CMPLX(TOL,0.0E0)
   ELSE
      CSCL = CMPLX(1.0E0/TOL,0.0E0)
   endif
   C1 = CW(1)*CSCL
   C2 = CW(2)*CSCL
   ST = Y(1)
!        ---------------------------------------------------------------
!        CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0E0/CABS(CT) PREVENTS
!        UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
!        ---------------------------------------------------------------
   CT = ZR*(C2+ST*C1)
   ACT = ABS(CT)
   RCT = CMPLX(1.0E0/ACT,0.0E0)
   CT = CONJG(CT)*RCT
   CINU = CINU*RCT*CT
   Y(1) = CINU*CSCL
   if (N /= 1) then
      DO 20 I = 2, N
         CINU = ST*CINU
         ST = Y(I)
         Y(I) = CINU*CSCL
   20       continue
   endif
  endif
  return
  END
  subroutine DETS17(Z,FNU,KODE,N,Y,NZ,NLAST,FNUL,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-763 (DEC 1989).
!
!     Original name: CUNI2
!
!     DETS17 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
!     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
!     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.
!
!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST /= 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
!     Y(I)=CZERO FOR I=NLAST+1,N
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, FNUL, TOL
  INTEGER           KODE, N, NLAST, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           AI, ARG, ASUM, BSUM, C1, C2, CFN, CI, CID, CONE, &
                    CRSC, CSCL, CZERO, DAI, PHI, RZ, S1, S2, ZB, &
                    ZETA1, ZETA2, ZN
  REAL              AARG, AIC, ANG, APHI, ASCLE, AY, C2I, C2M, C2R, &
                    CAR, FN, HPI, RS1, SAR, YY
  INTEGER           I, IDUM, IFLAG, IN, INU, J, K, NAI, ND, NDAI, &
                    NN, NUF, NW
!     .. Local Arrays ..
  COMPLEX           CIP(4), CSR(3), CSS(3), CY(2)
  REAL              BRY(3)
!     .. External functions ..
  REAL              X02AME, X02ALE
  EXTERNAL          X02AME, X02ALE
!     .. External subroutines ..
  EXTERNAL          DEUS17, DEVS17, S17DGE, DGVS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, CONJG, COS, EXP, INT, LOG, &
                    MAX, MIN, MOD, REAL, SIN
!     .. Data statements ..
  DATA              CZERO, CONE, CI/(0.0E0,0.0E0), (1.0E0,0.0E0), &
                    (0.0E0,1.0E0)/
  DATA              CIP(1), CIP(2), CIP(3), CIP(4)/(1.0E0,0.0E0), &
                    (0.0E0,1.0E0), (-1.0E0,0.0E0), (0.0E0,-1.0E0)/
  DATA              HPI, AIC/1.57079632679489662E+00, &
                    1.265512123484645396E+00/
!     .. Executable Statements ..
!
  NZ = 0
  ND = N
  NLAST = 0
!     ------------------------------------------------------------------
!     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
!     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
!     EXP(ALIM)=EXP(ELIM)*TOL
!     ------------------------------------------------------------------
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CRSC = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CRSC
  CSR(1) = CRSC
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = (1.0E+3*X02AME())/TOL
  YY = AIMAG(Z)
!     ------------------------------------------------------------------
!     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
!     ------------------------------------------------------------------
  ZN = -Z*CI
  ZB = Z
  CID = -CI
  INU = INT(FNU)
  ANG = HPI*(FNU-INU)
  CAR = COS(ANG)
  SAR = SIN(ANG)
  C2 = CMPLX(CAR,SAR)
  IN = INU + N - 1
  IN = MOD(IN,4)
  C2 = C2*CIP(IN+1)
  if (YY <= 0.0E0) then
   ZN = CONJG(-ZN)
   ZB = CONJG(ZB)
   CID = -CID
   C2 = CONJG(C2)
  endif
!     ------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!     ------------------------------------------------------------------
  FN = MAX(FNU,1.0E0)
  CALL DEUS17(ZN,FN,1,TOL,PHI,ARG,ZETA1,ZETA2,ASUM,BSUM,ELIM)
  if (KODE == 1) then
   S1 = -ZETA1 + ZETA2
  ELSE
   CFN = CMPLX(FNU,0.0E0)
   S1 = -ZETA1 + CFN*(CFN/(ZB+ZETA2))
  endif
  RS1 = REAL(S1)
  if (ABS(RS1) <= ELIM) then
   20    continue
   NN = MIN(2,ND)
   DO 40 I = 1, NN
      FN = FNU + ND - I
      CALL DEUS17(ZN,FN,0,TOL,PHI,ARG,ZETA1,ZETA2,ASUM,BSUM,ELIM)
      if (KODE == 1) then
         S1 = -ZETA1 + ZETA2
      ELSE
         CFN = CMPLX(FN,0.0E0)
         AY = ABS(YY)
         S1 = -ZETA1 + CFN*(CFN/(ZB+ZETA2)) + CMPLX(0.0E0,AY)
      endif
!           ------------------------------------------------------------
!           TEST FOR UNDERFLOW AND OVERFLOW
!           ------------------------------------------------------------
      RS1 = REAL(S1)
      if (ABS(RS1) > ELIM) then
         goto 60
      ELSE
         if (I == 1) IFLAG = 2
         if (ABS(RS1) >= ALIM) then
!                 ------------------------------------------------------
!                 REFINE  TEST AND SCALE
!                 ------------------------------------------------------
!                 ------------------------------------------------------
            APHI = ABS(PHI)
            AARG = ABS(ARG)
            RS1 = RS1 + LOG(APHI) - 0.25E0*LOG(AARG) - AIC
            if (ABS(RS1) > ELIM) then
               goto 60
            ELSE
               if (I == 1) IFLAG = 1
               if (RS1 >= 0.0E0) then
                  if (I == 1) IFLAG = 3
               endif
            endif
         endif
!              ---------------------------------------------------------
!              SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!              EXPONENT EXTREMES
!              ---------------------------------------------------------
         IDUM = 1
!              S17DGE assumed not to fail, therefore IDUM set to one.
         CALL S17DGE('F',ARG,'S',AI,NAI,IDUM)
         IDUM = 1
         CALL S17DGE('D',ARG,'S',DAI,NDAI,IDUM)
         S2 = PHI*(AI*ASUM+DAI*BSUM)
         C2R = REAL(S1)
         C2I = AIMAG(S1)
         C2M = EXP(C2R)*REAL(CSS(IFLAG))
         S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
         S2 = S2*S1
         if (IFLAG == 1) then
            CALL DGVS17(S2,NW,BRY(1),TOL)
            if (NW /= 0) goto 60
         endif
         if (YY <= 0.0E0) S2 = CONJG(S2)
         J = ND - I + 1
         S2 = S2*C2
         CY(I) = S2
         Y(J) = S2*CSR(IFLAG)
         C2 = C2*CID
      endif
   40    continue
   goto 80
   60    if (RS1 > 0.0E0) then
      goto 160
   ELSE
!           ------------------------------------------------------------
!           SET UNDERFLOW AND UPDATE PARAMETERS
!           ------------------------------------------------------------
      Y(ND) = CZERO
      NZ = NZ + 1
      ND = ND - 1
      if (ND == 0) then
         return
      ELSE
         CALL DEVS17(Z,FNU,KODE,1,ND,Y,NUF,TOL,ELIM,ALIM)
         if (NUF < 0) then
            goto 160
         ELSE
            ND = ND - NUF
            NZ = NZ + NUF
            if (ND == 0) then
               return
            ELSE
               FN = FNU + ND - 1
               if (FN < FNUL) then
                  goto 120
               ELSE
!                        FN = AIMAG(CID)
!                        J = NUF + 1
!                        K = MOD(J,4) + 1
!                        S1 = CIP(K)
!                        if (FN < 0.0E0) S1 = CONJG(S1)
!                        C2 = C2*S1
!                   The above 6 lines were replaced by the 5 below
!                   to fix a bug discovered during implementation
!                   on a Multics machine, whereby some results
!                   were returned wrongly scaled by sqrt(-1.0). MWP.
                  C2 = CMPLX(CAR,SAR)
                  IN = INU + ND - 1
                  IN = MOD(IN,4) + 1
                  C2 = C2*CIP(IN)
                  if (YY <= 0.0E0) C2 = CONJG(C2)
                  goto 20
               endif
            endif
         endif
      endif
   endif
   80    if (ND > 2) then
      RZ = CMPLX(2.0E0,0.0E0)/Z
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = X02ALE()
      S1 = CY(1)
      S2 = CY(2)
      C1 = CSR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = K
      DO 100 I = 3, ND
         C2 = S2
         S2 = S1 + CMPLX(FNU+FN,0.0E0)*RZ*S2
         S1 = C2
         C2 = S2*C1
         Y(K) = C2
         K = K - 1
         FN = FN - 1.0E0
         if (IFLAG < 3) then
            C2R = REAL(C2)
            C2I = AIMAG(C2)
            C2R = ABS(C2R)
            C2I = ABS(C2I)
            C2M = MAX(C2R,C2I)
            if (C2M > ASCLE) then
               IFLAG = IFLAG + 1
               ASCLE = BRY(IFLAG)
               S1 = S1*C1
               S2 = C2
               S1 = S1*CSS(IFLAG)
               S2 = S2*CSS(IFLAG)
               C1 = CSR(IFLAG)
            endif
         endif
  100       continue
   endif
   return
  120    NLAST = ND
   return
  else if (RS1 <= 0.0E0) then
   NZ = N
   DO 140 I = 1, N
      Y(I) = CZERO
  140    continue
   return
  endif
  160 NZ = -1
  return
  END
  subroutine DEUS17(Z,FNU,IPMTR,TOL,PHI,ARG,ZETA1,ZETA2,ASUM,BSUM, &
                    ELIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-764 (DEC 1989).
!
!     Original name: CUNHJ
!
!     REFERENCES
!         HANDBOOK OF MATHEMATICAL functionS BY M. ABRAMOWITZ AND I.A.
!         STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.
!
!         ASYMPTOTICS AND SPECIAL functionS BY F.W.J. OLVER, ACADEMIC
!         PRESS, N.Y., 1974, PAGE 420
!
!     ABSTRACT
!         DEUS17 COMPUTES PARAMETERS FOR BESSEL functionS C(FNU,Z) =
!         J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
!         BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION
!
!         C(FNU,Z)=C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )
!
!         FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
!         AN AIRY function AND DAIRY IS ITS DERIVATIVE.
!
!               (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,
!
!         ZETA1=0.5*FNU*CLOG((1+W)/(1-W)), ZETA2=FNU*W FOR SCALING
!         PURPOSES IN AIRY functionS FROM S17DGE OR S17DHE.
!
!         MCONJ=SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
!         MUST BE SPECIFIED. IPMTR=0 returnS ALL PARAMETERS. IPMTR=
!         1 COMPUTES ALL EXCEPT ASUM AND BSUM.
!
!     .. Scalar Arguments ..
  COMPLEX           ARG, ASUM, BSUM, PHI, Z, ZETA1, ZETA2
  REAL              ELIM, FNU, TOL
  INTEGER           IPMTR
!     .. Local Scalars ..
  COMPLEX           CFNU, CONE, CZERO, PRZTH, PTFN, RFN13, RTZTA, &
                    RZTH, SUMA, SUMB, T2, TFN, W, W2, ZA, ZB, ZC, &
                    ZETA, ZTH
  REAL              ANG, ASUMI, ASUMR, ATOL, AW2, AZTH, BSUMI, &
                    BSUMR, BTOL, EX1, EX2, FN13, FN23, HPI, PI, PP, &
                    RFNU, RFNU2, TEST, THPI, TSTI, TSTR, WI, WR, &
                    ZCI, ZCR, ZETAI, ZETAR, ZTHI, ZTHR
  INTEGER           IAS, IBS, IS, J, JR, JU, K, KMAX, KP1, KS, L, &
                    L1, L2, LR, LRP1, M
!     .. Local Arrays ..
  COMPLEX           CR(14), DR(14), P(30), UP(14)
  REAL              ALFA(180), AP(30), AR(14), BETA(210), BR(14), &
                    C(105), GAMA(30)
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, ATAN, CMPLX, COS, EXP, LOG, REAL, &
                    SIN, SQRT
!     .. Data statements ..
  DATA              AR(1), AR(2), AR(3), AR(4), AR(5), AR(6), AR(7), &
                    AR(8), AR(9), AR(10), AR(11), AR(12), AR(13), &
                    AR(14)/1.00000000000000000E+00, &
                    1.04166666666666667E-01, &
                    8.35503472222222222E-02, &
                    1.28226574556327160E-01, &
                    2.91849026464140464E-01, &
                    8.81627267443757652E-01, &
                    3.32140828186276754E+00, &
                    1.49957629868625547E+01, &
                    7.89230130115865181E+01, &
                    4.74451538868264323E+02, &
                    3.20749009089066193E+03, &
                    2.40865496408740049E+04, &
                    1.98923119169509794E+05, &
                    1.79190200777534383E+06/
  DATA              BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), &
                    BR(8), BR(9), BR(10), BR(11), BR(12), BR(13), &
                    BR(14)/1.00000000000000000E+00, &
                    -1.45833333333333333E-01, &
                    -9.87413194444444444E-02, &
                    -1.43312053915895062E-01, &
                    -3.17227202678413548E-01, &
                    -9.42429147957120249E-01, &
                    -3.51120304082635426E+00, &
                    -1.57272636203680451E+01, &
                    -8.22814390971859444E+01, &
                    -4.92355370523670524E+02, &
                    -3.31621856854797251E+03, &
                    -2.48276742452085896E+04, &
                    -2.04526587315129788E+05, &
                    -1.83844491706820990E+06/
  DATA              C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), &
                    C(9), C(10), C(11), C(12), C(13), C(14), C(15), &
                    C(16)/1.00000000000000000E+00, &
                    -2.08333333333333333E-01, &
                    1.25000000000000000E-01, &
                    3.34201388888888889E-01, &
                    -4.01041666666666667E-01, &
                    7.03125000000000000E-02, &
                    -1.02581259645061728E+00, &
                    1.84646267361111111E+00, &
                    -8.91210937500000000E-01, &
                    7.32421875000000000E-02, &
                    4.66958442342624743E+00, &
                    -1.12070026162229938E+01, &
                    8.78912353515625000E+00, &
                    -2.36408691406250000E+00, &
                    1.12152099609375000E-01, &
                    -2.82120725582002449E+01/
  DATA              C(17), C(18), C(19), C(20), C(21), C(22), C(23), &
                    C(24)/8.46362176746007346E+01, &
                    -9.18182415432400174E+01, &
                    4.25349987453884549E+01, &
                    -7.36879435947963170E+00, &
                    2.27108001708984375E-01, &
                    2.12570130039217123E+02, &
                    -7.65252468141181642E+02, &
                    1.05999045252799988E+03/
  DATA              C(25), C(26), C(27), C(28), C(29), C(30), C(31), &
                    C(32), C(33), C(34), C(35), C(36), C(37), C(38), &
                    C(39), C(40)/-6.99579627376132541E+02, &
                    2.18190511744211590E+02, &
                    -2.64914304869515555E+01, &
                    5.72501420974731445E-01, &
                    -1.91945766231840700E+03, &
                    8.06172218173730938E+03, &
                    -1.35865500064341374E+04, &
                    1.16553933368645332E+04, &
                    -5.30564697861340311E+03, &
                    1.20090291321635246E+03, &
                    -1.08090919788394656E+02, &
                    1.72772750258445740E+00, &
                    2.02042913309661486E+04, &
                    -9.69805983886375135E+04, &
                    1.92547001232531532E+05, &
                    -2.03400177280415534E+05/
  DATA              C(41), C(42), C(43), C(44), C(45), C(46), C(47), &
                    C(48)/1.22200464983017460E+05, &
                    -4.11926549688975513E+04, &
                    7.10951430248936372E+03, &
                    -4.93915304773088012E+02, &
                    6.07404200127348304E+00, &
                    -2.42919187900551333E+05, &
                    1.31176361466297720E+06, &
                    -2.99801591853810675E+06/
  DATA              C(49), C(50), C(51), C(52), C(53), C(54), C(55), &
                    C(56), C(57), C(58), C(59), C(60), C(61), C(62), &
                    C(63), C(64)/3.76327129765640400E+06, &
                    -2.81356322658653411E+06, &
                    1.26836527332162478E+06, &
                    -3.31645172484563578E+05, &
                    4.52187689813627263E+04, &
                    -2.49983048181120962E+03, &
                    2.43805296995560639E+01, &
                    3.28446985307203782E+06, &
                    -1.97068191184322269E+07, &
                    5.09526024926646422E+07, &
                    -7.41051482115326577E+07, &
                    6.63445122747290267E+07, &
                    -3.75671766607633513E+07, &
                    1.32887671664218183E+07, &
                    -2.78561812808645469E+06, &
                    3.08186404612662398E+05/
  DATA              C(65), C(66), C(67), C(68), C(69), C(70), C(71), &
                    C(72)/-1.38860897537170405E+04, &
                    1.10017140269246738E+02, &
                    -4.93292536645099620E+07, &
                    3.25573074185765749E+08, &
                    -9.39462359681578403E+08, &
                    1.55359689957058006E+09, &
                    -1.62108055210833708E+09, &
                    1.10684281682301447E+09/
  DATA              C(73), C(74), C(75), C(76), C(77), C(78), C(79), &
                    C(80), C(81), C(82), C(83), C(84), C(85), C(86), &
                    C(87), C(88)/-4.95889784275030309E+08, &
                    1.42062907797533095E+08, &
                    -2.44740627257387285E+07, &
                    2.24376817792244943E+06, &
                    -8.40054336030240853E+04, &
                    5.51335896122020586E+02, &
                    8.14789096118312115E+08, &
                    -5.86648149205184723E+09, &
                    1.86882075092958249E+10, &
                    -3.46320433881587779E+10, &
                    4.12801855797539740E+10, &
                    -3.30265997498007231E+10, &
                    1.79542137311556001E+10, &
                    -6.56329379261928433E+09, &
                    1.55927986487925751E+09, &
                    -2.25105661889415278E+08/
  DATA              C(89), C(90), C(91), C(92), C(93), C(94), C(95), &
                    C(96)/1.73951075539781645E+07, &
                    -5.49842327572288687E+05, &
                    3.03809051092238427E+03, &
                    -1.46792612476956167E+10, &
                    1.14498237732025810E+11, &
                    -3.99096175224466498E+11, &
                    8.19218669548577329E+11, &
                    -1.09837515608122331E+12/
  DATA              C(97), C(98), C(99), C(100), C(101), C(102), &
                    C(103), C(104), C(105)/1.00815810686538209E+12, &
                    -6.45364869245376503E+11, &
                    2.87900649906150589E+11, &
                    -8.78670721780232657E+10, &
                    1.76347306068349694E+10, &
                    -2.16716498322379509E+09, &
                    1.43157876718888981E+08, &
                    -3.87183344257261262E+06, &
                    1.82577554742931747E+04/
  DATA              ALFA(1), ALFA(2), ALFA(3), ALFA(4), ALFA(5), &
                    ALFA(6), ALFA(7), ALFA(8), ALFA(9), ALFA(10), &
                    ALFA(11), ALFA(12), ALFA(13), &
                    ALFA(14)/-4.44444444444444444E-03, &
                    -9.22077922077922078E-04, &
                    -8.84892884892884893E-05, &
                    1.65927687832449737E-04, &
                    2.46691372741792910E-04, &
                    2.65995589346254780E-04, &
                    2.61824297061500945E-04, &
                    2.48730437344655609E-04, &
                    2.32721040083232098E-04, &
                    2.16362485712365082E-04, &
                    2.00738858762752355E-04, &
                    1.86267636637545172E-04, &
                    1.73060775917876493E-04, &
                    1.61091705929015752E-04/
  DATA              ALFA(15), ALFA(16), ALFA(17), ALFA(18), &
                    ALFA(19), ALFA(20), ALFA(21), &
                    ALFA(22)/1.50274774160908134E-04, &
                    1.40503497391269794E-04, &
                    1.31668816545922806E-04, &
                    1.23667445598253261E-04, &
                    1.16405271474737902E-04, &
                    1.09798298372713369E-04, &
                    1.03772410422992823E-04, &
                    9.82626078369363448E-05/
  DATA              ALFA(23), ALFA(24), ALFA(25), ALFA(26), &
                    ALFA(27), ALFA(28), ALFA(29), ALFA(30), &
                    ALFA(31), ALFA(32), ALFA(33), ALFA(34), &
                    ALFA(35), ALFA(36)/9.32120517249503256E-05, &
                    8.85710852478711718E-05, &
                    8.42963105715700223E-05, &
                    8.03497548407791151E-05, &
                    7.66981345359207388E-05, &
                    7.33122157481777809E-05, &
                    7.01662625163141333E-05, &
                    6.72375633790160292E-05, &
                    6.93735541354588974E-04, &
                    2.32241745182921654E-04, &
                    -1.41986273556691197E-05, &
                    -1.16444931672048640E-04, &
                    -1.50803558053048762E-04, &
                    -1.55121924918096223E-04/
  DATA              ALFA(37), ALFA(38), ALFA(39), ALFA(40), &
                    ALFA(41), ALFA(42), ALFA(43), &
                    ALFA(44)/-1.46809756646465549E-04, &
                    -1.33815503867491367E-04, &
                    -1.19744975684254051E-04, &
                    -1.06184319207974020E-04, &
                    -9.37699549891194492E-05, &
                    -8.26923045588193274E-05, &
                    -7.29374348155221211E-05, &
                    -6.44042357721016283E-05/
  DATA              ALFA(45), ALFA(46), ALFA(47), ALFA(48), &
                    ALFA(49), ALFA(50), ALFA(51), ALFA(52), &
                    ALFA(53), ALFA(54), ALFA(55), ALFA(56), &
                    ALFA(57), ALFA(58)/-5.69611566009369048E-05, &
                    -5.04731044303561628E-05, &
                    -4.48134868008882786E-05, &
                    -3.98688727717598864E-05, &
                    -3.55400532972042498E-05, &
                    -3.17414256609022480E-05, &
                    -2.83996793904174811E-05, &
                    -2.54522720634870566E-05, &
                    -2.28459297164724555E-05, &
                    -2.05352753106480604E-05, &
                    -1.84816217627666085E-05, &
                    -1.66519330021393806E-05, &
                    -1.50179412980119482E-05, &
                    -1.35554031379040526E-05/
  DATA              ALFA(59), ALFA(60), ALFA(61), ALFA(62), &
                    ALFA(63), ALFA(64), ALFA(65), &
                    ALFA(66)/-1.22434746473858131E-05, &
                    -1.10641884811308169E-05, &
                    -3.54211971457743841E-04, &
                    -1.56161263945159416E-04, &
                    3.04465503594936410E-05, &
                    1.30198655773242693E-04, &
                    1.67471106699712269E-04, &
                    1.70222587683592569E-04/
  DATA              ALFA(67), ALFA(68), ALFA(69), ALFA(70), &
                    ALFA(71), ALFA(72), ALFA(73), ALFA(74), &
                    ALFA(75), ALFA(76), ALFA(77), ALFA(78), &
                    ALFA(79), ALFA(80)/1.56501427608594704E-04, &
                    1.36339170977445120E-04, &
                    1.14886692029825128E-04, &
                    9.45869093034688111E-05, &
                    7.64498419250898258E-05, &
                    6.07570334965197354E-05, &
                    4.74394299290508799E-05, &
                    3.62757512005344297E-05, &
                    2.69939714979224901E-05, &
                    1.93210938247939253E-05, &
                    1.30056674793963203E-05, &
                    7.82620866744496661E-06, &
                    3.59257485819351583E-06, &
                    1.44040049814251817E-07/
  DATA              ALFA(81), ALFA(82), ALFA(83), ALFA(84), &
                    ALFA(85), ALFA(86), ALFA(87), &
                    ALFA(88)/-2.65396769697939116E-06, &
                    -4.91346867098485910E-06, &
                    -6.72739296091248287E-06, &
                    -8.17269379678657923E-06, &
                    -9.31304715093561232E-06, &
                    -1.02011418798016441E-05, &
                    -1.08805962510592880E-05, &
                    -1.13875481509603555E-05/
  DATA              ALFA(89), ALFA(90), ALFA(91), ALFA(92), &
                    ALFA(93), ALFA(94), ALFA(95), ALFA(96), &
                    ALFA(97), ALFA(98), ALFA(99), ALFA(100), &
                    ALFA(101), ALFA(102)/-1.17519675674556414E-05, &
                    -1.19987364870944141E-05, &
                    3.78194199201772914E-04, &
                    2.02471952761816167E-04, &
                    -6.37938506318862408E-05, &
                    -2.38598230603005903E-04, &
                    -3.10916256027361568E-04, &
                    -3.13680115247576316E-04, &
                    -2.78950273791323387E-04, &
                    -2.28564082619141374E-04, &
                    -1.75245280340846749E-04, &
                    -1.25544063060690348E-04, &
                    -8.22982872820208365E-05, &
                    -4.62860730588116458E-05/
  DATA              ALFA(103), ALFA(104), ALFA(105), ALFA(106), &
                    ALFA(107), ALFA(108), ALFA(109), &
                    ALFA(110)/-1.72334302366962267E-05, &
                    5.60690482304602267E-06, &
                    2.31395443148286800E-05, &
                    3.62642745856793957E-05, &
                    4.58006124490188752E-05, &
                    5.24595294959114050E-05, &
                    5.68396208545815266E-05, &
                    5.94349820393104052E-05/
  DATA              ALFA(111), ALFA(112), ALFA(113), ALFA(114), &
                    ALFA(115), ALFA(116), ALFA(117), ALFA(118), &
                    ALFA(119), ALFA(120), ALFA(121), &
                    ALFA(122)/6.06478527578421742E-05, &
                    6.08023907788436497E-05, &
                    6.01577894539460388E-05, &
                    5.89199657344698500E-05, &
                    5.72515823777593053E-05, &
                    5.52804375585852577E-05, &
                    5.31063773802880170E-05, &
                    5.08069302012325706E-05, &
                    4.84418647620094842E-05, &
                    4.60568581607475370E-05, &
                    -6.91141397288294174E-04, &
                    -4.29976633058871912E-04/
  DATA              ALFA(123), ALFA(124), ALFA(125), ALFA(126), &
                    ALFA(127), ALFA(128), ALFA(129), &
                    ALFA(130)/1.83067735980039018E-04, &
                    6.60088147542014144E-04, &
                    8.75964969951185931E-04, &
                    8.77335235958235514E-04, &
                    7.49369585378990637E-04, &
                    5.63832329756980918E-04, &
                    3.68059319971443156E-04, &
                    1.88464535514455599E-04/
  DATA              ALFA(131), ALFA(132), ALFA(133), ALFA(134), &
                    ALFA(135), ALFA(136), ALFA(137), ALFA(138), &
                    ALFA(139), ALFA(140), ALFA(141), &
                    ALFA(142)/3.70663057664904149E-05, &
                    -8.28520220232137023E-05, &
                    -1.72751952869172998E-04, &
                    -2.36314873605872983E-04, &
                    -2.77966150694906658E-04, &
                    -3.02079514155456919E-04, &
                    -3.12594712643820127E-04, &
                    -3.12872558758067163E-04, &
                    -3.05678038466324377E-04, &
                    -2.93226470614557331E-04, &
                    -2.77255655582934777E-04, &
                    -2.59103928467031709E-04/
  DATA              ALFA(143), ALFA(144), ALFA(145), ALFA(146), &
                    ALFA(147), ALFA(148), ALFA(149), &
                    ALFA(150)/-2.39784014396480342E-04, &
                    -2.20048260045422848E-04, &
                    -2.00443911094971498E-04, &
                    -1.81358692210970687E-04, &
                    -1.63057674478657464E-04, &
                    -1.45712672175205844E-04, &
                    -1.29425421983924587E-04, &
                    -1.14245691942445952E-04/
  DATA              ALFA(151), ALFA(152), ALFA(153), ALFA(154), &
                    ALFA(155), ALFA(156), ALFA(157), ALFA(158), &
                    ALFA(159), ALFA(160), ALFA(161), &
                    ALFA(162)/1.92821964248775885E-03, &
                    1.35592576302022234E-03, &
                    -7.17858090421302995E-04, &
                    -2.58084802575270346E-03, &
                    -3.49271130826168475E-03, &
                    -3.46986299340960628E-03, &
                    -2.82285233351310182E-03, &
                    -1.88103076404891354E-03, &
                    -8.89531718383947600E-04, &
                    3.87912102631035228E-06, &
                    7.28688540119691412E-04, &
                    1.26566373053457758E-03/
  DATA              ALFA(163), ALFA(164), ALFA(165), ALFA(166), &
                    ALFA(167), ALFA(168), ALFA(169), &
                    ALFA(170)/1.62518158372674427E-03, &
                    1.83203153216373172E-03, &
                    1.91588388990527909E-03, &
                    1.90588846755546138E-03, &
                    1.82798982421825727E-03, &
                    1.70389506421121530E-03, &
                    1.55097127171097686E-03, &
                    1.38261421852276159E-03/
  DATA              ALFA(171), ALFA(172), ALFA(173), ALFA(174), &
                    ALFA(175), ALFA(176), ALFA(177), ALFA(178), &
                    ALFA(179), ALFA(180)/1.20881424230064774E-03, &
                    1.03676532638344962E-03, &
                    8.71437918068619115E-04, &
                    7.16080155297701002E-04, &
                    5.72637002558129372E-04, &
                    4.42089819465802277E-04, &
                    3.24724948503090564E-04, &
                    2.20342042730246599E-04, &
                    1.28412898401353882E-04, &
                    4.82005924552095464E-05/
  DATA              BETA(1), BETA(2), BETA(3), BETA(4), BETA(5), &
                    BETA(6), BETA(7), BETA(8), BETA(9), BETA(10), &
                    BETA(11), BETA(12), BETA(13), &
                    BETA(14)/1.79988721413553309E-02, &
                    5.59964911064388073E-03, &
                    2.88501402231132779E-03, &
                    1.80096606761053941E-03, &
                    1.24753110589199202E-03, &
                    9.22878876572938311E-04, &
                    7.14430421727287357E-04, &
                    5.71787281789704872E-04, &
                    4.69431007606481533E-04, &
                    3.93232835462916638E-04, &
                    3.34818889318297664E-04, &
                    2.88952148495751517E-04, &
                    2.52211615549573284E-04, &
                    2.22280580798883327E-04/
  DATA              BETA(15), BETA(16), BETA(17), BETA(18), &
                    BETA(19), BETA(20), BETA(21), &
                    BETA(22)/1.97541838033062524E-04, &
                    1.76836855019718004E-04, &
                    1.59316899661821081E-04, &
                    1.44347930197333986E-04, &
                    1.31448068119965379E-04, &
                    1.20245444949302884E-04, &
                    1.10449144504599392E-04, &
                    1.01828770740567258E-04/
  DATA              BETA(23), BETA(24), BETA(25), BETA(26), &
                    BETA(27), BETA(28), BETA(29), BETA(30), &
                    BETA(31), BETA(32), BETA(33), BETA(34), &
                    BETA(35), BETA(36)/9.41998224204237509E-05, &
                    8.74130545753834437E-05, &
                    8.13466262162801467E-05, &
                    7.59002269646219339E-05, &
                    7.09906300634153481E-05, &
                    6.65482874842468183E-05, &
                    6.25146958969275078E-05, &
                    5.88403394426251749E-05, &
                    -1.49282953213429172E-03, &
                    -8.78204709546389328E-04, &
                    -5.02916549572034614E-04, &
                    -2.94822138512746025E-04, &
                    -1.75463996970782828E-04, &
                    -1.04008550460816434E-04/
  DATA              BETA(37), BETA(38), BETA(39), BETA(40), &
                    BETA(41), BETA(42), BETA(43), &
                    BETA(44)/-5.96141953046457895E-05, &
                    -3.12038929076098340E-05, &
                    -1.26089735980230047E-05, &
                    -2.42892608575730389E-07, &
                    8.05996165414273571E-06, &
                    1.36507009262147391E-05, &
                    1.73964125472926261E-05, &
                    1.98672978842133780E-05/
  DATA              BETA(45), BETA(46), BETA(47), BETA(48), &
                    BETA(49), BETA(50), BETA(51), BETA(52), &
                    BETA(53), BETA(54), BETA(55), BETA(56), &
                    BETA(57), BETA(58)/2.14463263790822639E-05, &
                    2.23954659232456514E-05, &
                    2.28967783814712629E-05, &
                    2.30785389811177817E-05, &
                    2.30321976080909144E-05, &
                    2.28236073720348722E-05, &
                    2.25005881105292418E-05, &
                    2.20981015361991429E-05, &
                    2.16418427448103905E-05, &
                    2.11507649256220843E-05, &
                    2.06388749782170737E-05, &
                    2.01165241997081666E-05, &
                    1.95913450141179244E-05, &
                    1.90689367910436740E-05/
  DATA              BETA(59), BETA(60), BETA(61), BETA(62), &
                    BETA(63), BETA(64), BETA(65), &
                    BETA(66)/1.85533719641636667E-05, &
                    1.80475722259674218E-05, &
                    5.52213076721292790E-04, &
                    4.47932581552384646E-04, &
                    2.79520653992020589E-04, &
                    1.52468156198446602E-04, &
                    6.93271105657043598E-05, &
                    1.76258683069991397E-05/
  DATA              BETA(67), BETA(68), BETA(69), BETA(70), &
                    BETA(71), BETA(72), BETA(73), BETA(74), &
                    BETA(75), BETA(76), BETA(77), BETA(78), &
                    BETA(79), BETA(80)/-1.35744996343269136E-05, &
                    -3.17972413350427135E-05, &
                    -4.18861861696693365E-05, &
                    -4.69004889379141029E-05, &
                    -4.87665447413787352E-05, &
                    -4.87010031186735069E-05, &
                    -4.74755620890086638E-05, &
                    -4.55813058138628452E-05, &
                    -4.33309644511266036E-05, &
                    -4.09230193157750364E-05, &
                    -3.84822638603221274E-05, &
                    -3.60857167535410501E-05, &
                    -3.37793306123367417E-05, &
                    -3.15888560772109621E-05/
  DATA              BETA(81), BETA(82), BETA(83), BETA(84), &
                    BETA(85), BETA(86), BETA(87), &
                    BETA(88)/-2.95269561750807315E-05, &
                    -2.75978914828335759E-05, &
                    -2.58006174666883713E-05, &
                    -2.41308356761280200E-05, &
                    -2.25823509518346033E-05, &
                    -2.11479656768912971E-05, &
                    -1.98200638885294927E-05, &
                    -1.85909870801065077E-05/
  DATA              BETA(89), BETA(90), BETA(91), BETA(92), &
                    BETA(93), BETA(94), BETA(95), BETA(96), &
                    BETA(97), BETA(98), BETA(99), BETA(100), &
                    BETA(101), BETA(102)/-1.74532699844210224E-05, &
                    -1.63997823854497997E-05, &
                    -4.74617796559959808E-04, &
                    -4.77864567147321487E-04, &
                    -3.20390228067037603E-04, &
                    -1.61105016119962282E-04, &
                    -4.25778101285435204E-05, &
                    3.44571294294967503E-05, &
                    7.97092684075674924E-05, &
                    1.03138236708272200E-04, &
                    1.12466775262204158E-04, &
                    1.13103642108481389E-04, &
                    1.08651634848774268E-04, &
                    1.01437951597661973E-04/
  DATA              BETA(103), BETA(104), BETA(105), BETA(106), &
                    BETA(107), BETA(108), BETA(109), &
                    BETA(110)/9.29298396593363896E-05, &
                    8.40293133016089978E-05, &
                    7.52727991349134062E-05, &
                    6.69632521975730872E-05, &
                    5.92564547323194704E-05, &
                    5.22169308826975567E-05, &
                    4.58539485165360646E-05, &
                    4.01445513891486808E-05/
  DATA              BETA(111), BETA(112), BETA(113), BETA(114), &
                    BETA(115), BETA(116), BETA(117), BETA(118), &
                    BETA(119), BETA(120), BETA(121), &
                    BETA(122)/3.50481730031328081E-05, &
                    3.05157995034346659E-05, &
                    2.64956119950516039E-05, &
                    2.29363633690998152E-05, &
                    1.97893056664021636E-05, &
                    1.70091984636412623E-05, &
                    1.45547428261524004E-05, &
                    1.23886640995878413E-05, &
                    1.04775876076583236E-05, &
                    8.79179954978479373E-06, &
                    7.36465810572578444E-04, &
                    8.72790805146193976E-04/
  DATA              BETA(123), BETA(124), BETA(125), BETA(126), &
                    BETA(127), BETA(128), BETA(129), &
                    BETA(130)/6.22614862573135066E-04, &
                    2.85998154194304147E-04, &
                    3.84737672879366102E-06, &
                    -1.87906003636971558E-04, &
                    -2.97603646594554535E-04, &
                    -3.45998126832656348E-04, &
                    -3.53382470916037712E-04, &
                    -3.35715635775048757E-04/
  DATA              BETA(131), BETA(132), BETA(133), BETA(134), &
                    BETA(135), BETA(136), BETA(137), BETA(138), &
                    BETA(139), BETA(140), BETA(141), &
                    BETA(142)/-3.04321124789039809E-04, &
                    -2.66722723047612821E-04, &
                    -2.27654214122819527E-04, &
                    -1.89922611854562356E-04, &
                    -1.55058918599093870E-04, &
                    -1.23778240761873630E-04, &
                    -9.62926147717644187E-05, &
                    -7.25178327714425337E-05, &
                    -5.22070028895633801E-05, &
                    -3.50347750511900522E-05, &
                    -2.06489761035551757E-05, &
                    -8.70106096849767054E-06/
  DATA              BETA(143), BETA(144), BETA(145), BETA(146), &
                    BETA(147), BETA(148), BETA(149), &
                    BETA(150)/1.13698686675100290E-06, &
                    9.16426474122778849E-06, &
                    1.56477785428872620E-05, &
                    2.08223629482466847E-05, &
                    2.48923381004595156E-05, &
                    2.80340509574146325E-05, &
                    3.03987774629861915E-05, &
                    3.21156731406700616E-05/
  DATA              BETA(151), BETA(152), BETA(153), BETA(154), &
                    BETA(155), BETA(156), BETA(157), BETA(158), &
                    BETA(159), BETA(160), BETA(161), &
                    BETA(162)/-1.80182191963885708E-03, &
                    -2.43402962938042533E-03, &
                    -1.83422663549856802E-03, &
                    -7.62204596354009765E-04, &
                    2.39079475256927218E-04, &
                    9.49266117176881141E-04, &
                    1.34467449701540359E-03, &
                    1.48457495259449178E-03, &
                    1.44732339830617591E-03, &
                    1.30268261285657186E-03, &
                    1.10351597375642682E-03, &
                    8.86047440419791759E-04/
  DATA              BETA(163), BETA(164), BETA(165), BETA(166), &
                    BETA(167), BETA(168), BETA(169), &
                    BETA(170)/6.73073208165665473E-04, &
                    4.77603872856582378E-04, &
                    3.05991926358789362E-04, &
                    1.60315694594721630E-04, &
                    4.00749555270613286E-05, &
                    -5.66607461635251611E-05, &
                    -1.32506186772982638E-04, &
                    -1.90296187989614057E-04/
  DATA              BETA(171), BETA(172), BETA(173), BETA(174), &
                    BETA(175), BETA(176), BETA(177), BETA(178), &
                    BETA(179), BETA(180), BETA(181), &
                    BETA(182)/-2.32811450376937408E-04, &
                    -2.62628811464668841E-04, &
                    -2.82050469867598672E-04, &
                    -2.93081563192861167E-04, &
                    -2.97435962176316616E-04, &
                    -2.96557334239348078E-04, &
                    -2.91647363312090861E-04, &
                    -2.83696203837734166E-04, &
                    -2.73512317095673346E-04, &
                    -2.61750155806768580E-04, &
                    6.38585891212050914E-03, &
                    9.62374215806377941E-03/
  DATA              BETA(183), BETA(184), BETA(185), BETA(186), &
                    BETA(187), BETA(188), BETA(189), &
                    BETA(190)/7.61878061207001043E-03, &
                    2.83219055545628054E-03, &
                    -2.09841352012720090E-03, &
                    -5.73826764216626498E-03, &
                    -7.70804244495414620E-03, &
                    -8.21011692264844401E-03, &
                    -7.65824520346905413E-03, &
                    -6.47209729391045177E-03/
  DATA              BETA(191), BETA(192), BETA(193), BETA(194), &
                    BETA(195), BETA(196), BETA(197), BETA(198), &
                    BETA(199), BETA(200), BETA(201), &
                    BETA(202)/-4.99132412004966473E-03, &
                    -3.45612289713133280E-03, &
                    -2.01785580014170775E-03, &
                    -7.59430686781961401E-04, &
                    2.84173631523859138E-04, &
                    1.10891667586337403E-03, &
                    1.72901493872728771E-03, &
                    2.16812590802684701E-03, &
                    2.45357710494539735E-03, &
                    2.61281821058334862E-03, &
                    2.67141039656276912E-03, &
                    2.65203073395980430E-03/
  DATA              BETA(203), BETA(204), BETA(205), BETA(206), &
                    BETA(207), BETA(208), BETA(209), &
                    BETA(210)/2.57411652877287315E-03, &
                    2.45389126236094427E-03, &
                    2.30460058071795494E-03, &
                    2.13684837686712662E-03, &
                    1.95896528478870911E-03, &
                    1.77737008679454412E-03, &
                    1.59690280765839059E-03, &
                    1.42111975664438546E-03/
  DATA              GAMA(1), GAMA(2), GAMA(3), GAMA(4), GAMA(5), &
                    GAMA(6), GAMA(7), GAMA(8), GAMA(9), GAMA(10), &
                    GAMA(11), GAMA(12), GAMA(13), &
                    GAMA(14)/6.29960524947436582E-01, &
                    2.51984209978974633E-01, &
                    1.54790300415655846E-01, &
                    1.10713062416159013E-01, &
                    8.57309395527394825E-02, &
                    6.97161316958684292E-02, &
                    5.86085671893713576E-02, &
                    5.04698873536310685E-02, &
                    4.42600580689154809E-02, &
                    3.93720661543509966E-02, &
                    3.54283195924455368E-02, &
                    3.21818857502098231E-02, &
                    2.94646240791157679E-02, &
                    2.71581677112934479E-02/
  DATA              GAMA(15), GAMA(16), GAMA(17), GAMA(18), &
                    GAMA(19), GAMA(20), GAMA(21), &
                    GAMA(22)/2.51768272973861779E-02, &
                    2.34570755306078891E-02, &
                    2.19508390134907203E-02, &
                    2.06210828235646240E-02, &
                    1.94388240897880846E-02, &
                    1.83810633800683158E-02, &
                    1.74293213231963172E-02, &
                    1.65685837786612353E-02/
  DATA              GAMA(23), GAMA(24), GAMA(25), GAMA(26), &
                    GAMA(27), GAMA(28), GAMA(29), &
                    GAMA(30)/1.57865285987918445E-02, &
                    1.50729501494095594E-02, &
                    1.44193250839954639E-02, &
                    1.38184805735341786E-02, &
                    1.32643378994276568E-02, &
                    1.27517121970498651E-02, &
                    1.22761545318762767E-02, &
                    1.18338262398482403E-02/
  DATA              EX1, EX2, HPI, PI, THPI/3.33333333333333333E-01, &
                    6.66666666666666667E-01, &
                    1.57079632679489662E+00, &
                    3.14159265358979324E+00, &
                    4.71238898038468986E+00/
  DATA              CZERO, CONE/(0.0E0,0.0E0), (1.0E0,0.0E0)/
!     .. Executable Statements ..
!
  RFNU = 1.0E0/FNU
  TSTR = REAL(Z)
  TSTI = AIMAG(Z)
  TEST = FNU*EXP(-ELIM)
  if (ABS(TSTR) < TEST) TSTR = 0.0E0
  if (ABS(TSTI) < TEST) TSTI = 0.0E0
  if (TSTR == 0.0E0 .and. TSTI == 0.0E0) then
   ZETA1 = CMPLX(ELIM+ELIM+FNU,0.0E0)
   ZETA2 = CMPLX(FNU,0.0E0)
   PHI = CONE
   ARG = CONE
   return
  endif
  ZB = CMPLX(TSTR,TSTI)*CMPLX(RFNU,0.0E0)
  RFNU2 = RFNU*RFNU
!     ------------------------------------------------------------------
!     COMPUTE IN THE FOURTH QUADRANT
!     ------------------------------------------------------------------
  FN13 = FNU**EX1
  FN23 = FN13*FN13
  RFN13 = CMPLX(1.0E0/FN13,0.0E0)
  W2 = CONE - ZB*ZB
  AW2 = ABS(W2)
  if (AW2 > 0.25E0) then
!        ---------------------------------------------------------------
!        CABS(W2)>0.25E0
!        ---------------------------------------------------------------
   W = SQRT(W2)
   WR = REAL(W)
   WI = AIMAG(W)
   if (WR < 0.0E0) WR = 0.0E0
   if (WI < 0.0E0) WI = 0.0E0
   W = CMPLX(WR,WI)
   ZA = (CONE+W)/ZB
   ZC = LOG(ZA)
   ZCR = REAL(ZC)
   ZCI = AIMAG(ZC)
   if (ZCI < 0.0E0) ZCI = 0.0E0
   if (ZCI > HPI) ZCI = HPI
   if (ZCR < 0.0E0) ZCR = 0.0E0
   ZC = CMPLX(ZCR,ZCI)
   ZTH = (ZC-W)*CMPLX(1.5E0,0.0E0)
   CFNU = CMPLX(FNU,0.0E0)
   ZETA1 = ZC*CFNU
   ZETA2 = W*CFNU
   AZTH = ABS(ZTH)
   ZTHR = REAL(ZTH)
   ZTHI = AIMAG(ZTH)
   ANG = THPI
   if (ZTHR < 0.0E0 .or. ZTHI >= 0.0E0) then
      ANG = HPI
      if (ZTHR /= 0.0E0) then
         ANG = ATAN(ZTHI/ZTHR)
         if (ZTHR < 0.0E0) ANG = ANG + PI
      endif
   endif
   PP = AZTH**EX2
   ANG = ANG*EX2
   ZETAR = PP*COS(ANG)
   ZETAI = PP*SIN(ANG)
   if (ZETAI < 0.0E0) ZETAI = 0.0E0
   ZETA = CMPLX(ZETAR,ZETAI)
   ARG = ZETA*CMPLX(FN23,0.0E0)
   RTZTA = ZTH/ZETA
   ZA = RTZTA/W
   PHI = SQRT(ZA+ZA)*RFN13
   if (IPMTR /= 1) then
      TFN = CMPLX(RFNU,0.0E0)/W
      RZTH = CMPLX(RFNU,0.0E0)/ZTH
      ZC = RZTH*CMPLX(AR(2),0.0E0)
      T2 = CONE/W2
      UP(2) = (T2*CMPLX(C(2),0.0E0)+CMPLX(C(3),0.0E0))*TFN
      BSUM = UP(2) + ZC
      ASUM = CZERO
      if (RFNU >= TOL) then
         PRZTH = RZTH
         PTFN = TFN
         UP(1) = CONE
         PP = 1.0E0
         BSUMR = REAL(BSUM)
         BSUMI = AIMAG(BSUM)
         BTOL = TOL*(ABS(BSUMR)+ABS(BSUMI))
         KS = 0
         KP1 = 2
         L = 3
         IAS = 0
         IBS = 0
         DO 100 LR = 2, 12, 2
            LRP1 = LR + 1
!                 ------------------------------------------------------
!                 COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE
!                 TERMS IN NEXT SUMA AND SUMB
!                 ------------------------------------------------------
            DO 40 K = LR, LRP1
               KS = KS + 1
               KP1 = KP1 + 1
               L = L + 1
               ZA = CMPLX(C(L),0.0E0)
               DO 20 J = 2, KP1
                  L = L + 1
                  ZA = ZA*T2 + CMPLX(C(L),0.0E0)
   20                continue
               PTFN = PTFN*TFN
               UP(KP1) = PTFN*ZA
               CR(KS) = PRZTH*CMPLX(BR(KS+1),0.0E0)
               PRZTH = PRZTH*RZTH
               DR(KS) = PRZTH*CMPLX(AR(KS+2),0.0E0)
   40             continue
            PP = PP*RFNU2
            if (IAS /= 1) then
               SUMA = UP(LRP1)
               JU = LRP1
               DO 60 JR = 1, LR
                  JU = JU - 1
                  SUMA = SUMA + CR(JR)*UP(JU)
   60                continue
               ASUM = ASUM + SUMA
               ASUMR = REAL(ASUM)
               ASUMI = AIMAG(ASUM)
               TEST = ABS(ASUMR) + ABS(ASUMI)
               if (PP < TOL .and. TEST < TOL) IAS = 1
            endif
            if (IBS /= 1) then
               SUMB = UP(LR+2) + UP(LRP1)*ZC
               JU = LRP1
               DO 80 JR = 1, LR
                  JU = JU - 1
                  SUMB = SUMB + DR(JR)*UP(JU)
   80                continue
               BSUM = BSUM + SUMB
               BSUMR = REAL(BSUM)
               BSUMI = AIMAG(BSUM)
               TEST = ABS(BSUMR) + ABS(BSUMI)
               if (PP < BTOL .and. TEST < TOL) IBS = 1
            endif
            if (IAS == 1 .and. IBS == 1) goto 120
  100          continue
      endif
  120       ASUM = ASUM + CONE
      BSUM = -BSUM*RFN13/RTZTA
   endif
  ELSE
!        ---------------------------------------------------------------
!        POWER SERIES FOR CABS(W2) <= 0.25E0
!        ---------------------------------------------------------------
   K = 1
   P(1) = CONE
   SUMA = CMPLX(GAMA(1),0.0E0)
   AP(1) = 1.0E0
   if (AW2 >= TOL) then
      DO 140 K = 2, 30
         P(K) = P(K-1)*W2
         SUMA = SUMA + P(K)*CMPLX(GAMA(K),0.0E0)
         AP(K) = AP(K-1)*AW2
         if (AP(K) < TOL) goto 160
  140       continue
      K = 30
   endif
  160    KMAX = K
   ZETA = W2*SUMA
   ARG = ZETA*CMPLX(FN23,0.0E0)
   ZA = SQRT(SUMA)
   ZETA2 = SQRT(W2)*CMPLX(FNU,0.0E0)
   ZETA1 = ZETA2*(CONE+ZETA*ZA*CMPLX(EX2,0.0E0))
   ZA = ZA + ZA
   PHI = SQRT(ZA)*RFN13
   if (IPMTR /= 1) then
!           ------------------------------------------------------------
!           SUM SERIES FOR ASUM AND BSUM
!           ------------------------------------------------------------
      SUMB = CZERO
      DO 180 K = 1, KMAX
         SUMB = SUMB + P(K)*CMPLX(BETA(K),0.0E0)
  180       continue
      ASUM = CZERO
      BSUM = SUMB
      L1 = 0
      L2 = 30
      BTOL = TOL*ABS(BSUM)
      ATOL = TOL
      PP = 1.0E0
      IAS = 0
      IBS = 0
      if (RFNU2 >= TOL) then
         DO 280 IS = 2, 7
            ATOL = ATOL/RFNU2
            PP = PP*RFNU2
            if (IAS /= 1) then
               SUMA = CZERO
               DO 200 K = 1, KMAX
                  M = L1 + K
                  SUMA = SUMA + P(K)*CMPLX(ALFA(M),0.0E0)
                  if (AP(K) < ATOL) goto 220
  200                continue
  220                ASUM = ASUM + SUMA*CMPLX(PP,0.0E0)
               if (PP < TOL) IAS = 1
            endif
            if (IBS /= 1) then
               SUMB = CZERO
               DO 240 K = 1, KMAX
                  M = L2 + K
                  SUMB = SUMB + P(K)*CMPLX(BETA(M),0.0E0)
                  if (AP(K) < ATOL) goto 260
  240                continue
  260                BSUM = BSUM + SUMB*CMPLX(PP,0.0E0)
               if (PP < BTOL) IBS = 1
            endif
            if (IAS == 1 .and. IBS == 1) then
               goto 300
            ELSE
               L1 = L1 + 30
               L2 = L2 + 30
            endif
  280          continue
      endif
  300       ASUM = ASUM + CONE
      PP = RFNU*REAL(RFN13)
      BSUM = BSUM*CMPLX(PP,0.0E0)
   endif
  endif
  return
  END
  subroutine DEVS17(Z,FNU,KODE,IKFLG,N,Y,NUF,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-765 (DEC 1989).
!
!     Original name: CUOIK
!
!     DEVS17 COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
!     EXPANSIONS FOR THE I AND K functionS AND COMPARES THEM
!     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
!     WHERE ALIM < ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
!     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
!     THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER
!     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
!     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
!     EXP(-ELIM)/TOL
!
!     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
!          =2 MEANS THE K SEQUENCE IS TESTED
!     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
!         =-1 MEANS AN OVERFLOW WOULD OCCUR
!     IKFLG=1 AND NUF>0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
!             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
!     IKFLG=2 AND NUF==N MEANS ALL Y VALUES WERE SET TO ZERO
!     IKFLG=2 AND 0 < NUF < N NOT CONSIDERED. Y MUST BE SET BY
!             ANOTHER ROUTINE
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, TOL
  INTEGER           IKFLG, KODE, N, NUF
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           ARG, ASUM, BSUM, CZ, CZERO, PHI, SUM, ZB, ZETA1, &
                    ZETA2, ZN, ZR
  REAL              AARG, AIC, APHI, ASCLE, AX, AY, FNN, GNN, GNU, &
                    RCZ, X, YY
  INTEGER           I, IFORM, INIT, NN, NW
!     .. Local Arrays ..
  COMPLEX           CWRK(16)
!     .. External functions ..
  REAL              X02AME
  EXTERNAL          X02AME
!     .. External subroutines ..
  EXTERNAL          DEUS17, DEWS17, DGVS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, CONJG, COS, EXP, LOG, MAX, &
                    REAL, SIN
!     .. Data statements ..
  DATA              CZERO/(0.0E0,0.0E0)/
  DATA              AIC/1.265512123484645396E+00/
!     .. Executable Statements ..
!
  NUF = 0
  NN = N
  X = REAL(Z)
  ZR = Z
  if (X < 0.0E0) ZR = -Z
  ZB = ZR
  YY = AIMAG(ZR)
  AX = ABS(X)*1.7321E0
  AY = ABS(YY)
  IFORM = 1
  if (AY > AX) IFORM = 2
  GNU = MAX(FNU,1.0E0)
  if (IKFLG /= 1) then
   FNN = NN
   GNN = FNU + FNN - 1.0E0
   GNU = MAX(GNN,FNN)
  endif
!     ------------------------------------------------------------------
!     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
!     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
!     THE SIGN OF THE IMAGINARY PART CORRECT.
!     ------------------------------------------------------------------
  if (IFORM == 2) then
   ZN = -ZR*CMPLX(0.0E0,1.0E0)
   if (YY <= 0.0E0) ZN = CONJG(-ZN)
   CALL DEUS17(ZN,GNU,1,TOL,PHI,ARG,ZETA1,ZETA2,ASUM,BSUM,ELIM)
   CZ = -ZETA1 + ZETA2
   AARG = ABS(ARG)
  ELSE
   INIT = 0
   CALL DEWS17(ZR,GNU,IKFLG,1,TOL,INIT,PHI,ZETA1,ZETA2,SUM,CWRK, &
                 ELIM)
   CZ = -ZETA1 + ZETA2
  endif
  if (KODE == 2) CZ = CZ - ZB
  if (IKFLG == 2) CZ = -CZ
  APHI = ABS(PHI)
  RCZ = REAL(CZ)
!     ------------------------------------------------------------------
!     OVERFLOW TEST
!     ------------------------------------------------------------------
  if (RCZ <= ELIM) then
   if (RCZ < ALIM) then
!           ------------------------------------------------------------
!           UNDERFLOW TEST
!           ------------------------------------------------------------
      if (RCZ >= (-ELIM)) then
         if (RCZ > (-ALIM)) then
            goto 40
         ELSE
            RCZ = RCZ + LOG(APHI)
            if (IFORM == 2) RCZ = RCZ - 0.25E0*LOG(AARG) - AIC
            if (RCZ > (-ELIM)) then
               ASCLE = (1.0E+3*X02AME())/TOL
               CZ = CZ + LOG(PHI)
               if (IFORM /= 1) CZ = CZ - CMPLX(0.25E0,0.0E0) &
                                      *LOG(ARG) - CMPLX(AIC,0.0E0)
               AX = EXP(RCZ)/TOL
               AY = AIMAG(CZ)
               CZ = CMPLX(AX,0.0E0)*CMPLX(COS(AY),SIN(AY))
               CALL DGVS17(CZ,NW,ASCLE,TOL)
               if (NW /= 1) goto 40
            endif
         endif
      endif
      DO 20 I = 1, NN
         Y(I) = CZERO
   20       continue
      NUF = NN
      return
   ELSE
      RCZ = RCZ + LOG(APHI)
      if (IFORM == 2) RCZ = RCZ - 0.25E0*LOG(AARG) - AIC
      if (RCZ > ELIM) goto 80
   endif
   40    if (IKFLG /= 2) then
      if (N /= 1) then
   60          continue
!              ---------------------------------------------------------
!              SET UNDERFLOWS ON I SEQUENCE
!              ---------------------------------------------------------
         GNU = FNU + NN - 1
         if (IFORM == 2) then
            CALL DEUS17(ZN,GNU,1,TOL,PHI,ARG,ZETA1,ZETA2,ASUM, &
                          BSUM,ELIM)
            CZ = -ZETA1 + ZETA2
            AARG = ABS(ARG)
         ELSE
            INIT = 0
            CALL DEWS17(ZR,GNU,IKFLG,1,TOL,INIT,PHI,ZETA1,ZETA2, &
                          SUM,CWRK,ELIM)
            CZ = -ZETA1 + ZETA2
         endif
         if (KODE == 2) CZ = CZ - ZB
         APHI = ABS(PHI)
         RCZ = REAL(CZ)
         if (RCZ >= (-ELIM)) then
            if (RCZ > (-ALIM)) then
               return
            ELSE
               RCZ = RCZ + LOG(APHI)
               if (IFORM == 2) RCZ = RCZ - 0.25E0*LOG(AARG) - AIC
               if (RCZ > (-ELIM)) then
                  ASCLE = (1.0E+3*X02AME())/TOL
                  CZ = CZ + LOG(PHI)
                  if (IFORM /= 1) CZ = CZ - CMPLX(0.25E0,0.0E0) &
                                         *LOG(ARG) - CMPLX(AIC, &
                                         0.0E0)
                  AX = EXP(RCZ)/TOL
                  AY = AIMAG(CZ)
                  CZ = CMPLX(AX,0.0E0)*CMPLX(COS(AY),SIN(AY))
                  CALL DGVS17(CZ,NW,ASCLE,TOL)
                  if (NW /= 1) return
               endif
            endif
         endif
         Y(NN) = CZERO
         NN = NN - 1
         NUF = NUF + 1
         if (NN /= 0) goto 60
      endif
   endif
   return
  endif
   80 NUF = -1
  return
  END
  subroutine DEWS17(ZR,FNU,IKFLG,IPMTR,TOL,INIT,PHI,ZETA1,ZETA2,SUM, &
                    CWRK,ELIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-766 (DEC 1989).
!
!     Original name: CUNIK
!
!        DEWS17 COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC
!        EXPANSIONS OF THE I AND K functionS ON IKFLG= 1 OR 2
!        RESPECTIVELY BY
!
!        W(FNU,ZR) = PHI*EXP(ZETA)*SUM
!
!        WHERE       ZETA=-ZETA1 + ZETA2       OR
!                          ZETA1 - ZETA2
!
!        THE FIRST CALL MUST HAVE INIT=0. SUBSEQUENT CALLS WITH THE
!        SAME ZR AND FNU WILL return THE I OR K function ON IKFLG=
!        1 OR 2 WITH NO CHANGE IN INIT. CWRK IS A COMPLEX WORK
!        ARRAY. IPMTR=0 COMPUTES ALL PARAMETERS. IPMTR=1 COMPUTES PHI,
!        ZETA1,ZETA2.
!
!     .. Scalar Arguments ..
  COMPLEX           PHI, SUM, ZETA1, ZETA2, ZR
  REAL              ELIM, FNU, TOL
  INTEGER           IKFLG, INIT, IPMTR
!     .. Array Arguments ..
  COMPLEX           CWRK(16)
!     .. Local Scalars ..
  COMPLEX           CFN, CONE, CRFN, CZERO, S, SR, T, T2, ZN
  REAL              AC, RFN, TEST, TSTI, TSTR
  INTEGER           I, J, K, L
!     .. Local Arrays ..
  COMPLEX           CON(2)
  REAL              C(120)
!bc
!     .. external functions ..
  real              x02ane
  external          x02ane
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, EXP, LOG, REAL, SQRT
!     .. Data statements ..
  DATA              CZERO, CONE/(0.0E0,0.0E0), (1.0E0,0.0E0)/
  DATA              CON(1), CON(2)/(3.98942280401432678E-01,0.0E0), &
                    (1.25331413731550025E+00,0.0E0)/
  DATA              C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), &
                    C(9), C(10), C(11), C(12), C(13), C(14), C(15), &
                    C(16)/1.00000000000000000E+00, &
                    -2.08333333333333333E-01, &
                    1.25000000000000000E-01, &
                    3.34201388888888889E-01, &
                    -4.01041666666666667E-01, &
                    7.03125000000000000E-02, &
                    -1.02581259645061728E+00, &
                    1.84646267361111111E+00, &
                    -8.91210937500000000E-01, &
                    7.32421875000000000E-02, &
                    4.66958442342624743E+00, &
                    -1.12070026162229938E+01, &
                    8.78912353515625000E+00, &
                    -2.36408691406250000E+00, &
                    1.12152099609375000E-01, &
                    -2.82120725582002449E+01/
  DATA              C(17), C(18), C(19), C(20), C(21), C(22), C(23), &
                    C(24)/8.46362176746007346E+01, &
                    -9.18182415432400174E+01, &
                    4.25349987453884549E+01, &
                    -7.36879435947963170E+00, &
                    2.27108001708984375E-01, &
                    2.12570130039217123E+02, &
                    -7.65252468141181642E+02, &
                    1.05999045252799988E+03/
  DATA              C(25), C(26), C(27), C(28), C(29), C(30), C(31), &
                    C(32), C(33), C(34), C(35), C(36), C(37), C(38), &
                    C(39), C(40)/-6.99579627376132541E+02, &
                    2.18190511744211590E+02, &
                    -2.64914304869515555E+01, &
                    5.72501420974731445E-01, &
                    -1.91945766231840700E+03, &
                    8.06172218173730938E+03, &
                    -1.35865500064341374E+04, &
                    1.16553933368645332E+04, &
                    -5.30564697861340311E+03, &
                    1.20090291321635246E+03, &
                    -1.08090919788394656E+02, &
                    1.72772750258445740E+00, &
                    2.02042913309661486E+04, &
                    -9.69805983886375135E+04, &
                    1.92547001232531532E+05, &
                    -2.03400177280415534E+05/
  DATA              C(41), C(42), C(43), C(44), C(45), C(46), C(47), &
                    C(48)/1.22200464983017460E+05, &
                    -4.11926549688975513E+04, &
                    7.10951430248936372E+03, &
                    -4.93915304773088012E+02, &
                    6.07404200127348304E+00, &
                    -2.42919187900551333E+05, &
                    1.31176361466297720E+06, &
                    -2.99801591853810675E+06/
  DATA              C(49), C(50), C(51), C(52), C(53), C(54), C(55), &
                    C(56), C(57), C(58), C(59), C(60), C(61), C(62), &
                    C(63), C(64)/3.76327129765640400E+06, &
                    -2.81356322658653411E+06, &
                    1.26836527332162478E+06, &
                    -3.31645172484563578E+05, &
                    4.52187689813627263E+04, &
                    -2.49983048181120962E+03, &
                    2.43805296995560639E+01, &
                    3.28446985307203782E+06, &
                    -1.97068191184322269E+07, &
                    5.09526024926646422E+07, &
                    -7.41051482115326577E+07, &
                    6.63445122747290267E+07, &
                    -3.75671766607633513E+07, &
                    1.32887671664218183E+07, &
                    -2.78561812808645469E+06, &
                    3.08186404612662398E+05/
  DATA              C(65), C(66), C(67), C(68), C(69), C(70), C(71), &
                    C(72)/-1.38860897537170405E+04, &
                    1.10017140269246738E+02, &
                    -4.93292536645099620E+07, &
                    3.25573074185765749E+08, &
                    -9.39462359681578403E+08, &
                    1.55359689957058006E+09, &
                    -1.62108055210833708E+09, &
                    1.10684281682301447E+09/
  DATA              C(73), C(74), C(75), C(76), C(77), C(78), C(79), &
                    C(80), C(81), C(82), C(83), C(84), C(85), C(86), &
                    C(87), C(88)/-4.95889784275030309E+08, &
                    1.42062907797533095E+08, &
                    -2.44740627257387285E+07, &
                    2.24376817792244943E+06, &
                    -8.40054336030240853E+04, &
                    5.51335896122020586E+02, &
                    8.14789096118312115E+08, &
                    -5.86648149205184723E+09, &
                    1.86882075092958249E+10, &
                    -3.46320433881587779E+10, &
                    4.12801855797539740E+10, &
                    -3.30265997498007231E+10, &
                    1.79542137311556001E+10, &
                    -6.56329379261928433E+09, &
                    1.55927986487925751E+09, &
                    -2.25105661889415278E+08/
  DATA              C(89), C(90), C(91), C(92), C(93), C(94), C(95), &
                    C(96)/1.73951075539781645E+07, &
                    -5.49842327572288687E+05, &
                    3.03809051092238427E+03, &
                    -1.46792612476956167E+10, &
                    1.14498237732025810E+11, &
                    -3.99096175224466498E+11, &
                    8.19218669548577329E+11, &
                    -1.09837515608122331E+12/
  DATA              C(97), C(98), C(99), C(100), C(101), C(102), &
                    C(103), C(104), C(105), C(106), C(107), C(108), &
                    C(109), C(110)/1.00815810686538209E+12, &
                    -6.45364869245376503E+11, &
                    2.87900649906150589E+11, &
                    -8.78670721780232657E+10, &
                    1.76347306068349694E+10, &
                    -2.16716498322379509E+09, &
                    1.43157876718888981E+08, &
                    -3.87183344257261262E+06, &
                    1.82577554742931747E+04, &
                    2.86464035717679043E+11, &
                    -2.40629790002850396E+12, &
                    9.10934118523989896E+12, &
                    -2.05168994109344374E+13, &
                    3.05651255199353206E+13/
  DATA              C(111), C(112), C(113), C(114), C(115), C(116), &
                    C(117), C(118), C(119), &
                    C(120)/-3.16670885847851584E+13, &
                    2.33483640445818409E+13, &
                    -1.23204913055982872E+13, &
                    4.61272578084913197E+12, &
                    -1.19655288019618160E+12, &
                    2.05914503232410016E+11, &
                    -2.18229277575292237E+10, &
                    1.24700929351271032E+09, &
                    -2.91883881222208134E+07, &
                    1.18838426256783253E+05/
!     .. Executable Statements ..
!
  if (INIT == 0) then
!        ---------------------------------------------------------------
!        INITIALIZE ALL VARIABLES
!        ---------------------------------------------------------------
   RFN = 1.0E0/FNU
   CRFN = CMPLX(RFN,0.0E0)
   TSTR = REAL(ZR)
   TSTI = AIMAG(ZR)
   TEST = FNU*EXP(-ELIM)
   if (ABS(TSTR) < TEST) TSTR = 0.0E0
   if (ABS(TSTI) < TEST) TSTI = 0.0E0
!bc         if (TSTR==0.0E0 .and. TSTI==0.0E0) then
   if (abs(tstr) <= x02ane() .and. abs(tsti) <= x02ane()) then
      ZETA1 = CMPLX(ELIM+ELIM+FNU,0.0E0)
      ZETA2 = CMPLX(FNU,0.0E0)
      PHI = CONE
      return
   endif
   T = CMPLX(TSTR,TSTI)*CRFN
   S = CONE + T*T
   SR = SQRT(S)
   CFN = CMPLX(FNU,0.0E0)
   ZN = (CONE+SR)/T
   ZETA1 = CFN*LOG(ZN)
   ZETA2 = CFN*SR
   T = CONE/SR
   SR = T*CRFN
   CWRK(16) = SQRT(SR)
   PHI = CWRK(16)*CON(IKFLG)
   if (IPMTR /= 0) then
      return
   ELSE
      T2 = CONE/S
      CWRK(1) = CONE
      CRFN = CONE
      AC = 1.0E0
      L = 1
      DO 40 K = 2, 15
         S = CZERO
         DO 20 J = 1, K
            L = L + 1
            S = S*T2 + CMPLX(C(L),0.0E0)
   20          continue
         CRFN = CRFN*SR
         CWRK(K) = CRFN*S
         AC = AC*RFN
         TSTR = REAL(CWRK(K))
         TSTI = AIMAG(CWRK(K))
         TEST = ABS(TSTR) + ABS(TSTI)
         if (AC < TOL .and. TEST < TOL) goto 60
   40       continue
      K = 15
   60       INIT = K
   endif
  endif
  if (IKFLG == 2) then
!        ---------------------------------------------------------------
!        COMPUTE SUM FOR THE K function
!        ---------------------------------------------------------------
   S = CZERO
   T = CONE
   DO 80 I = 1, INIT
      S = S + T*CWRK(I)
      T = -T
   80    continue
   SUM = S
   PHI = CWRK(16)*CON(2)
  ELSE
!        ---------------------------------------------------------------
!        COMPUTE SUM FOR THE I function
!        ---------------------------------------------------------------
   S = CZERO
   DO 100 I = 1, INIT
      S = S + CWRK(I)
  100    continue
   SUM = S
   PHI = CWRK(16)*CON(1)
  endif
  return
  END
  subroutine DEXS17(Z,FNU,KODE,N,Y,NZ,NLAST,FNUL,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-767 (DEC 1989).
!
!     Original name: CUNI1
!
!     DEXS17 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
!     EXPANSION FOR I(FNU,Z) IN -PI/3 <= ARG Z <= PI/3.
!
!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST /= 0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
!     Y(I)=CZERO FOR I=NLAST+1,N
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, FNUL, TOL
  INTEGER           KODE, N, NLAST, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           C1, C2, CFN, CONE, CRSC, CSCL, CZERO, PHI, RZ, &
                    S1, S2, SUM, ZETA1, ZETA2
  REAL              APHI, ASCLE, C2I, C2M, C2R, FN, RS1, YY
  INTEGER           I, IFLAG, INIT, K, M, ND, NN, NUF, NW
!     .. Local Arrays ..
  COMPLEX           CSR(3), CSS(3), CWRK(16), CY(2)
  REAL              BRY(3)
!     .. External functions ..
  REAL              X02AME, X02ALE
  EXTERNAL          X02AME, X02ALE
!     .. External subroutines ..
  EXTERNAL          DEVS17, DEWS17, DGVS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, COS, EXP, LOG, MAX, MIN, &
                    REAL, SIN
!     .. Data statements ..
  DATA              CZERO, CONE/(0.0E0,0.0E0), (1.0E0,0.0E0)/
!     .. Executable Statements ..
!
  NZ = 0
  ND = N
  NLAST = 0
!     ------------------------------------------------------------------
!     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
!     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
!     EXP(ALIM)=EXP(ELIM)*TOL
!     ------------------------------------------------------------------
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CRSC = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CRSC
  CSR(1) = CRSC
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = (1.0E+3*X02AME())/TOL
!     ------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!     ------------------------------------------------------------------
  FN = MAX(FNU,1.0E0)
  INIT = 0
  CALL DEWS17(Z,FN,1,1,TOL,INIT,PHI,ZETA1,ZETA2,SUM,CWRK,ELIM)
  if (KODE == 1) then
   S1 = -ZETA1 + ZETA2
  ELSE
   CFN = CMPLX(FN,0.0E0)
   S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2))
  endif
  RS1 = REAL(S1)
  if (ABS(RS1) <= ELIM) then
   20    continue
   NN = MIN(2,ND)
   DO 40 I = 1, NN
      FN = FNU + ND - I
      INIT = 0
      CALL DEWS17(Z,FN,1,0,TOL,INIT,PHI,ZETA1,ZETA2,SUM,CWRK,ELIM)
      if (KODE == 1) then
         S1 = -ZETA1 + ZETA2
      ELSE
         CFN = CMPLX(FN,0.0E0)
         YY = AIMAG(Z)
         S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2)) + CMPLX(0.0E0,YY)
      endif
!           ------------------------------------------------------------
!           TEST FOR UNDERFLOW AND OVERFLOW
!           ------------------------------------------------------------
      RS1 = REAL(S1)
      if (ABS(RS1) > ELIM) then
         goto 60
      ELSE
         if (I == 1) IFLAG = 2
         if (ABS(RS1) >= ALIM) then
!                 ------------------------------------------------------
!                 REFINE  TEST AND SCALE
!                 ------------------------------------------------------
            APHI = ABS(PHI)
            RS1 = RS1 + LOG(APHI)
            if (ABS(RS1) > ELIM) then
               goto 60
            ELSE
               if (I == 1) IFLAG = 1
               if (RS1 >= 0.0E0) then
                  if (I == 1) IFLAG = 3
               endif
            endif
         endif
!              ---------------------------------------------------------
!              SCALE S1 IF CABS(S1) < ASCLE
!              ---------------------------------------------------------
         S2 = PHI*SUM
         C2R = REAL(S1)
         C2I = AIMAG(S1)
         C2M = EXP(C2R)*REAL(CSS(IFLAG))
         S1 = CMPLX(C2M,0.0E0)*CMPLX(COS(C2I),SIN(C2I))
         S2 = S2*S1
         if (IFLAG == 1) then
            CALL DGVS17(S2,NW,BRY(1),TOL)
            if (NW /= 0) goto 60
         endif
         M = ND - I + 1
         CY(I) = S2
         Y(M) = S2*CSR(IFLAG)
      endif
   40    continue
   goto 80
!        ---------------------------------------------------------------
!        SET UNDERFLOW AND UPDATE PARAMETERS
!        ---------------------------------------------------------------
   60    continue
   if (RS1 > 0.0E0) then
      goto 160
   ELSE
      Y(ND) = CZERO
      NZ = NZ + 1
      ND = ND - 1
      if (ND == 0) then
         return
      ELSE
         CALL DEVS17(Z,FNU,KODE,1,ND,Y,NUF,TOL,ELIM,ALIM)
         if (NUF < 0) then
            goto 160
         ELSE
            ND = ND - NUF
            NZ = NZ + NUF
            if (ND == 0) then
               return
            ELSE
               FN = FNU + ND - 1
               if (FN >= FNUL) then
                  goto 20
               ELSE
                  goto 120
               endif
            endif
         endif
      endif
   endif
   80    if (ND > 2) then
      RZ = CMPLX(2.0E0,0.0E0)/Z
      BRY(2) = 1.0E0/BRY(1)
      BRY(3) = X02ALE()
      S1 = CY(1)
      S2 = CY(2)
      C1 = CSR(IFLAG)
      ASCLE = BRY(IFLAG)
      K = ND - 2
      FN = K
      DO 100 I = 3, ND
         C2 = S2
         S2 = S1 + CMPLX(FNU+FN,0.0E0)*RZ*S2
         S1 = C2
         C2 = S2*C1
         Y(K) = C2
         K = K - 1
         FN = FN - 1.0E0
         if (IFLAG < 3) then
            C2R = REAL(C2)
            C2I = AIMAG(C2)
            C2R = ABS(C2R)
            C2I = ABS(C2I)
            C2M = MAX(C2R,C2I)
            if (C2M > ASCLE) then
               IFLAG = IFLAG + 1
               ASCLE = BRY(IFLAG)
               S1 = S1*C1
               S2 = C2
               S1 = S1*CSS(IFLAG)
               S2 = S2*CSS(IFLAG)
               C1 = CSR(IFLAG)
            endif
         endif
  100       continue
   endif
   return
  120    NLAST = ND
   return
  else if (RS1 <= 0.0E0) then
   NZ = N
   DO 140 I = 1, N
      Y(I) = CZERO
  140    continue
   return
  endif
  160 NZ = -1
  return
  END
  subroutine DEYS17(Z,FNU,KODE,N,Y,NZ,NUI,NLAST,FNUL,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-768 (DEC 1989).
!
!     Original name: CBUNI
!
!     DEYS17 COMPUTES THE I BESSEL function FOR LARGE CABS(Z)>
!     FNUL AND FNU+N-1 < FNUL. THE ORDER IS INCREASED FROM
!     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
!     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
!     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, FNUL, TOL
  INTEGER           KODE, N, NLAST, NUI, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           CSCL, CSCR, RZ, S1, S2, ST
  REAL              ASCLE, AX, AY, DFNU, FNUI, GNU, STI, STM, STR, &
                    XX, YY
  INTEGER           I, IFLAG, IFORM, K, NL, NW
!     .. Local Arrays ..
  COMPLEX           CY(2)
  REAL              BRY(3)
!     .. External functions ..
  REAL              X02AME
  EXTERNAL          X02AME
!     .. External subroutines ..
  EXTERNAL          DETS17, DEXS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, MAX, REAL
!     .. Executable Statements ..
!
  NZ = 0
  XX = REAL(Z)
  YY = AIMAG(Z)
  AX = ABS(XX)*1.7321E0
  AY = ABS(YY)
  IFORM = 1
  if (AY > AX) IFORM = 2
  if (NUI == 0) then
   if (IFORM == 2) then
!           ------------------------------------------------------------
!           ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!           APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
!           AND HPI=PI/2
!           ------------------------------------------------------------
      CALL DETS17(Z,FNU,KODE,N,Y,NW,NLAST,FNUL,TOL,ELIM,ALIM)
   ELSE
!           ------------------------------------------------------------
!           ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
!           -PI/3 <= ARG(Z) <= PI/3
!           ------------------------------------------------------------
      CALL DEXS17(Z,FNU,KODE,N,Y,NW,NLAST,FNUL,TOL,ELIM,ALIM)
   endif
   if (NW >= 0) then
      NZ = NW
      return
   endif
  ELSE
   FNUI = NUI
   DFNU = FNU + N - 1
   GNU = DFNU + FNUI
   if (IFORM == 2) then
!           ------------------------------------------------------------
!           ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!           APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
!           AND HPI=PI/2
!           ------------------------------------------------------------
      CALL DETS17(Z,GNU,KODE,2,CY,NW,NLAST,FNUL,TOL,ELIM,ALIM)
   ELSE
!           ------------------------------------------------------------
!           ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
!           -PI/3 <= ARG(Z) <= PI/3
!           ------------------------------------------------------------
      CALL DEXS17(Z,GNU,KODE,2,CY,NW,NLAST,FNUL,TOL,ELIM,ALIM)
   endif
   if (NW >= 0) then
      if (NW /= 0) then
         NLAST = N
      ELSE
         AY = ABS(CY(1))
!              ---------------------------------------------------------
!              SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER
!              USED
!              ---------------------------------------------------------
         BRY(1) = (1.0E+3*X02AME())/TOL
         BRY(2) = 1.0E0/BRY(1)
         BRY(3) = BRY(2)
         IFLAG = 2
         ASCLE = BRY(2)
         AX = 1.0E0
         CSCL = CMPLX(AX,0.0E0)
         if (AY <= BRY(1)) then
            IFLAG = 1
            ASCLE = BRY(1)
            AX = 1.0E0/TOL
            CSCL = CMPLX(AX,0.0E0)
         else if (AY >= BRY(2)) then
            IFLAG = 3
            ASCLE = BRY(3)
            AX = TOL
            CSCL = CMPLX(AX,0.0E0)
         endif
         AY = 1.0E0/AX
         CSCR = CMPLX(AY,0.0E0)
         S1 = CY(2)*CSCL
         S2 = CY(1)*CSCL
         RZ = CMPLX(2.0E0,0.0E0)/Z
         DO 20 I = 1, NUI
            ST = S2
            S2 = CMPLX(DFNU+FNUI,0.0E0)*RZ*S2 + S1
            S1 = ST
            FNUI = FNUI - 1.0E0
            if (IFLAG < 3) then
               ST = S2*CSCR
               STR = REAL(ST)
               STI = AIMAG(ST)
               STR = ABS(STR)
               STI = ABS(STI)
               STM = MAX(STR,STI)
               if (STM > ASCLE) then
                  IFLAG = IFLAG + 1
                  ASCLE = BRY(IFLAG)
                  S1 = S1*CSCR
                  S2 = ST
                  AX = AX*TOL
                  AY = 1.0E0/AX
                  CSCL = CMPLX(AX,0.0E0)
                  CSCR = CMPLX(AY,0.0E0)
                  S1 = S1*CSCL
                  S2 = S2*CSCL
               endif
            endif
   20          continue
         Y(N) = S2*CSCR
         if (N /= 1) then
            NL = N - 1
            FNUI = NL
            K = NL
            DO 40 I = 1, NL
               ST = S2
               S2 = CMPLX(FNU+FNUI,0.0E0)*RZ*S2 + S1
               S1 = ST
               ST = S2*CSCR
               Y(K) = ST
               FNUI = FNUI - 1.0E0
               K = K - 1
               if (IFLAG < 3) then
                  STR = REAL(ST)
                  STI = AIMAG(ST)
                  STR = ABS(STR)
                  STI = ABS(STI)
                  STM = MAX(STR,STI)
                  if (STM > ASCLE) then
                     IFLAG = IFLAG + 1
                     ASCLE = BRY(IFLAG)
                     S1 = S1*CSCR
                     S2 = ST
                     AX = AX*TOL
                     AY = 1.0E0/AX
                     CSCL = CMPLX(AX,0.0E0)
                     CSCR = CMPLX(AY,0.0E0)
                     S1 = S1*CSCL
                     S2 = S2*CSCL
                  endif
               endif
   40             continue
         endif
      endif
      return
   endif
  endif
  NZ = -1
  if (NW == (-2)) NZ = -2
  return
  END
  subroutine DEZS17(Z,FNU,KODE,N,CY,NZ,RL,FNUL,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-769 (DEC 1989).
!
!     Original name: CBINU
!
!     DEZS17 COMPUTES THE I function IN THE RIGHT HALF Z PLANE
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, FNUL, RL, TOL
  INTEGER           KODE, N, NZ
!     .. Array Arguments ..
  COMPLEX           CY(N)
!     .. Local Scalars ..
  COMPLEX           CZERO
  REAL              AZ, DFNU
  INTEGER           I, INW, NLAST, NN, NUI, NW
!     .. Local Arrays ..
  COMPLEX           CW(2)
!     .. External subroutines ..
  EXTERNAL          DESS17, DEVS17, DEYS17, DGRS17, DGTS17, DGYS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, INT, MAX
!     .. Data statements ..
  DATA              CZERO/(0.0E0,0.0E0)/
!     .. Executable Statements ..
!
  NZ = 0
  AZ = ABS(Z)
  NN = N
  DFNU = FNU + N - 1
  if (AZ > 2.0E0) then
   if (AZ*AZ*0.25E0 > DFNU+1.0E0) goto 20
  endif
!     ------------------------------------------------------------------
!     POWER SERIES
!     ------------------------------------------------------------------
  CALL DGRS17(Z,FNU,KODE,NN,CY,NW,TOL,ELIM,ALIM)
  INW = ABS(NW)
  NZ = NZ + INW
  NN = NN - INW
  if (NN == 0) then
   return
  else if (NW >= 0) then
   return
  ELSE
   DFNU = FNU + NN - 1
  endif
   20 if (AZ >= RL) then
   if (DFNU > 1.0E0) then
      if (AZ+AZ < DFNU*DFNU) goto 40
   endif
!        ---------------------------------------------------------------
!        ASYMPTOTIC EXPANSION FOR LARGE Z
!        ---------------------------------------------------------------
   CALL DGYS17(Z,FNU,KODE,NN,CY,NW,RL,TOL,ELIM,ALIM)
   if (NW < 0) then
      goto 120
   ELSE
      return
   endif
  else if (DFNU <= 1.0E0) then
   goto 100
  endif
!     ------------------------------------------------------------------
!     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
!     ------------------------------------------------------------------
   40 CALL DEVS17(Z,FNU,KODE,1,NN,CY,NW,TOL,ELIM,ALIM)
  if (NW < 0) then
   goto 120
  ELSE
   NZ = NZ + NW
   NN = NN - NW
   if (NN == 0) then
      return
   ELSE
      DFNU = FNU + NN - 1
      if (DFNU <= FNUL) then
         if (AZ <= FNUL) goto 60
      endif
!           ------------------------------------------------------------
!           INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
!           ------------------------------------------------------------
      NUI = INT(FNUL-DFNU) + 1
      NUI = MAX(NUI,0)
      CALL DEYS17(Z,FNU,KODE,NN,CY,NW,NUI,NLAST,FNUL,TOL,ELIM, &
                    ALIM)
      if (NW < 0) then
         goto 120
      ELSE
         NZ = NZ + NW
         if (NLAST == 0) then
            return
         ELSE
            NN = NLAST
         endif
      endif
   60       if (AZ > RL) then
!              ---------------------------------------------------------
!              MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
!              ---------------------------------------------------------
!              ---------------------------------------------------------
!              OVERFLOW TEST ON K functionS USED IN WRONSKIAN
!              ---------------------------------------------------------
         CALL DEVS17(Z,FNU,KODE,2,2,CW,NW,TOL,ELIM,ALIM)
         if (NW < 0) then
            NZ = NN
            DO 80 I = 1, NN
               CY(I) = CZERO
   80             continue
            return
         else if (NW > 0) then
            goto 120
         ELSE
            CALL DESS17(Z,FNU,KODE,NN,CY,NW,CW,TOL,ELIM,ALIM)
            if (NW < 0) then
               goto 120
            ELSE
               return
            endif
         endif
      endif
   endif
  endif
!     ------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES
!     ------------------------------------------------------------------
  100 CALL DGTS17(Z,FNU,KODE,NN,CY,NW,TOL)
  if (NW >= 0) return
  120 NZ = -1
  if (NW == (-2)) NZ = -2
  if (NW == (-3)) NZ = -3
  return
  END
  subroutine DGRS17(Z,FNU,KODE,N,Y,NZ,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-771 (DEC 1989).
!
!     Original name: CSERI
!
!     DGRS17 COMPUTES THE I BESSEL function FOR REAL(Z) >= 0.0 BY
!     MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE
!     REGION CABS(Z) <= 2*SQRT(FNU+1). NZ=0 IS A NORMAL return.
!     NZ>0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
!     DUE TO UNDERFLOW. NZ < 0 MEANS UNDERFLOW OCCURRED, BUT THE
!     CONDITION CABS(Z) <= 2*SQRT(FNU+1) WAS VIOLATED AND THE
!     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, TOL
  INTEGER           KODE, N, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           AK1, CK, COEF, CONE, CRSC, CZ, CZERO, HZ, RZ, &
                    S1, S2
  REAL              AA, ACZ, AK, ARM, ASCLE, ATOL, AZ, DFNU, FNUP, &
                    RAK1, RS, RTR1, S, SS, X
  INTEGER           I, IB, IDUM, IFLAG, IL, K, L, M, NN, NW
!     .. Local Arrays ..
  COMPLEX           W(2)
!     .. External functions ..
  REAL              S14ABE, X02AME
  EXTERNAL          S14ABE, X02AME
!     .. External subroutines ..
  EXTERNAL          DGVS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, COS, EXP, LOG, MIN, REAL, &
                    SIN, SQRT
!     .. Data statements ..
  DATA              CZERO, CONE/(0.0E0,0.0E0), (1.0E0,0.0E0)/
!     .. Executable Statements ..
!
  NZ = 0
  AZ = ABS(Z)
  if (AZ /= 0.0E0) then
   X = REAL(Z)
   ARM = 1.0E+3*X02AME()
   RTR1 = SQRT(ARM)
   CRSC = CMPLX(1.0E0,0.0E0)
   IFLAG = 0
   if (AZ < ARM) then
      NZ = N
      if (FNU == 0.0E0) NZ = NZ - 1
   ELSE
      HZ = Z*CMPLX(0.5E0,0.0E0)
      CZ = CZERO
      if (AZ > RTR1) CZ = HZ*HZ
      ACZ = ABS(CZ)
      NN = N
      CK = LOG(HZ)
   20       continue
      DFNU = FNU + NN - 1
      FNUP = DFNU + 1.0E0
!           ------------------------------------------------------------
!           UNDERFLOW TEST
!           ------------------------------------------------------------
      AK1 = CK*CMPLX(DFNU,0.0E0)
      IDUM = 0
!           S14ABE assumed not to fail, therefore IDUM set to zero.
      AK = S14ABE(FNUP,IDUM)
      AK1 = AK1 - CMPLX(AK,0.0E0)
      if (KODE == 2) AK1 = AK1 - CMPLX(X,0.0E0)
      RAK1 = REAL(AK1)
      if (RAK1 > (-ELIM)) then
         if (RAK1 <= (-ALIM)) then
            IFLAG = 1
            SS = 1.0E0/TOL
            CRSC = CMPLX(TOL,0.0E0)
            ASCLE = ARM*SS
         endif
         AK = AIMAG(AK1)
         AA = EXP(RAK1)
         if (IFLAG == 1) AA = AA*SS
         COEF = CMPLX(AA,0.0E0)*CMPLX(COS(AK),SIN(AK))
         ATOL = TOL*ACZ/FNUP
         IL = MIN(2,NN)
         DO 60 I = 1, IL
            DFNU = FNU + NN - I
            FNUP = DFNU + 1.0E0
            S1 = CONE
            if (ACZ >= TOL*FNUP) then
               AK1 = CONE
               AK = FNUP + 2.0E0
               S = FNUP
               AA = 2.0E0
   40                continue
               RS = 1.0E0/S
               AK1 = AK1*CZ*CMPLX(RS,0.0E0)
               S1 = S1 + AK1
               S = S + AK
               AK = AK + 2.0E0
               AA = AA*ACZ*RS
               if (AA > ATOL) goto 40
            endif
            M = NN - I + 1
            S2 = S1*COEF
            W(I) = S2
            if (IFLAG /= 0) then
               CALL DGVS17(S2,NW,ASCLE,TOL)
               if (NW /= 0) goto 80
            endif
            Y(M) = S2*CRSC
            if (I /= IL) COEF = COEF*CMPLX(DFNU,0.0E0)/HZ
   60          continue
         goto 100
      endif
   80       NZ = NZ + 1
      Y(NN) = CZERO
      if (ACZ > DFNU) then
         goto 180
      ELSE
         NN = NN - 1
         if (NN == 0) then
            return
         ELSE
            goto 20
         endif
      endif
  100       if (NN > 2) then
         K = NN - 2
         AK = K
         RZ = (CONE+CONE)/Z
         if (IFLAG == 1) then
!                 ------------------------------------------------------
!                 RECUR BACKWARD WITH SCALED VALUES
!                 ------------------------------------------------------
!                 ------------------------------------------------------
!                 EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE
!                 THE UNDERFLOW LIMIT = ASCLE = X02AME()*CSCL*1.0E+3
!                 ------------------------------------------------------
            S1 = W(1)
            S2 = W(2)
            DO 120 L = 3, NN
               CK = S2
               S2 = S1 + CMPLX(AK+FNU,0.0E0)*RZ*S2
               S1 = CK
               CK = S2*CRSC
               Y(K) = CK
               AK = AK - 1.0E0
               K = K - 1
               if (ABS(CK) > ASCLE) goto 140
  120             continue
            return
  140             IB = L + 1
            if (IB > NN) return
         ELSE
            IB = 3
         endif
         DO 160 I = IB, NN
            Y(K) = CMPLX(AK+FNU,0.0E0)*RZ*Y(K+1) + Y(K+2)
            AK = AK - 1.0E0
            K = K - 1
  160          continue
      endif
      return
!           ------------------------------------------------------------
!           return WITH NZ < 0 IF CABS(Z*Z/4)>FNU+N-NZ-1 COMPLETE
!           THE CALCULATION IN DEZS17 WITH N=N-IABS(NZ)
!           ------------------------------------------------------------
  180       continue
      NZ = -NZ
      return
   endif
  endif
  Y(1) = CZERO
  if (FNU == 0.0E0) Y(1) = CONE
  if (N /= 1) then
   DO 200 I = 2, N
      Y(I) = CZERO
  200    continue
  endif
  return
  END
  subroutine DGSS17(ZR,S1,S2,NZ,ASCLE,ALIM,IUF)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-772 (DEC 1989).
!
!     Original name: CS1S2
!
!     DGSS17 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE
!     ADDITION OF THE I AND K functionS IN THE ANALYTIC CON-
!     TINUATION FORMULA WHERE S1=K function AND S2=I function.
!     ON KODE=1 THE I AND K functionS ARE DIFFERENT ORDERS OF
!     MAGNITUDE, BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER
!     OF MAGNITUDE AND THE MAXIMUM MUST BE AT LEAST ONE
!     PRECISION ABOVE THE UNDERFLOW LIMIT.
!
!     .. Scalar Arguments ..
  COMPLEX           S1, S2, ZR
  REAL              ALIM, ASCLE
  INTEGER           IUF, NZ
!     .. Local Scalars ..
  COMPLEX           C1, CZERO, S1D
  REAL              AA, ALN, AS1, AS2, XX
  INTEGER           IF1
!     .. External functions ..
  COMPLEX           S01EAE
  EXTERNAL          S01EAE
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, LOG, MAX, REAL
!     .. Data statements ..
  DATA              CZERO/(0.0E0,0.0E0)/
!     .. Executable Statements ..
!
  NZ = 0
  AS1 = ABS(S1)
  AS2 = ABS(S2)
  AA = REAL(S1)
  ALN = AIMAG(S1)
  if (AA /= 0.0E0 .or. ALN /= 0.0E0) then
   if (AS1 /= 0.0E0) then
      XX = REAL(ZR)
      ALN = -XX - XX + LOG(AS1)
      S1D = S1
      S1 = CZERO
      AS1 = 0.0E0
      if (ALN >= (-ALIM)) then
         C1 = LOG(S1D) - ZR - ZR
!               S1 = EXP(C1)
         IF1 = 1
         S1 = S01EAE(C1,IF1)
         AS1 = ABS(S1)
         IUF = IUF + 1
      endif
   endif
  endif
  AA = MAX(AS1,AS2)
  if (AA <= ASCLE) then
   S1 = CZERO
   S2 = CZERO
   NZ = 1
   IUF = 0
  endif
  return
  END
  subroutine DGTS17(Z,FNU,KODE,N,Y,NZ,TOL)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-773 (DEC 1989).
!     Mark 17 REVISED. IER-1703 (JUN 1995).
!
!     Original name: CMLRI
!
!     DGTS17 COMPUTES THE I BESSEL function FOR RE(Z) >= 0.0 BY THE
!     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              FNU, TOL
  INTEGER           KODE, N, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           CK, CNORM, CONE, CTWO, CZERO, P1, P2, PT, RZ, &
                    SUM
  REAL              ACK, AK, AP, AT, AZ, BK, FKAP, FKK, FLAM, FNF, &
                    RHO, RHO2, SCLE, TFNF, TST, X
  INTEGER           I, IAZ, IDUM, IFL, IFNU, INU, ITIME, K, KK, KM, &
                    M
!     .. External functions ..
  COMPLEX           S01EAE
  REAL              S14ABE, X02ANE
  EXTERNAL          S14ABE, S01EAE, X02ANE
!     .. Intrinsic functions ..
  INTRINSIC         ABS, CMPLX, CONJG, EXP, INT, LOG, MAX, MIN, &
                    REAL, SQRT
!     .. Data statements ..
  DATA              CZERO, CONE, CTWO/(0.0E0,0.0E0), (1.0E0,0.0E0), &
                    (2.0E0,0.0E0)/
!     .. Executable Statements ..
!
  SCLE = (1.0E+3*X02ANE())/TOL
  NZ = 0
  AZ = ABS(Z)
  X = REAL(Z)
  IAZ = INT(AZ)
  IFNU = INT(FNU)
  INU = IFNU + N - 1
  AT = IAZ + 1.0E0
  CK = CMPLX(AT,0.0E0)/Z
  RZ = CTWO/Z
  P1 = CZERO
  P2 = CONE
  ACK = (AT+1.0E0)/AZ
  RHO = ACK + SQRT(ACK*ACK-1.0E0)
  RHO2 = RHO*RHO
  TST = (RHO2+RHO2)/((RHO2-1.0E0)*(RHO-1.0E0))
  TST = TST/TOL
!     ------------------------------------------------------------------
!     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
!     ------------------------------------------------------------------
  AK = AT
  DO 20 I = 1, 80
   PT = P2
   P2 = P1 - CK*P2
   P1 = PT
   CK = CK + RZ
   AP = ABS(P2)
   if (AP > TST*AK*AK) then
      goto 40
   ELSE
      AK = AK + 1.0E0
   endif
   20 continue
  goto 180
   40 I = I + 1
  K = 0
  if (INU >= IAZ) then
!        ---------------------------------------------------------------
!        COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
!        ---------------------------------------------------------------
   P1 = CZERO
   P2 = CONE
   AT = INU + 1.0E0
   CK = CMPLX(AT,0.0E0)/Z
   ACK = AT/AZ
   TST = SQRT(ACK/TOL)
   ITIME = 1
   DO 60 K = 1, 80
      PT = P2
      P2 = P1 - CK*P2
      P1 = PT
      CK = CK + RZ
      AP = ABS(P2)
      if (AP >= TST) then
         if (ITIME == 2) then
            goto 80
         ELSE
            ACK = ABS(CK)
            FLAM = ACK + SQRT(ACK*ACK-1.0E0)
            FKAP = AP/ABS(P1)
            RHO = MIN(FLAM,FKAP)
            TST = TST*SQRT(RHO/(RHO*RHO-1.0E0))
            ITIME = 2
         endif
      endif
   60    continue
   goto 180
  endif
!     ------------------------------------------------------------------
!     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
!     ------------------------------------------------------------------
   80 K = K + 1
  KK = MAX(I+IAZ,K+INU)
  FKK = KK
  P1 = CZERO
!     ------------------------------------------------------------------
!     SCALE P2 AND SUM BY SCLE
!     ------------------------------------------------------------------
  P2 = CMPLX(SCLE,0.0E0)
  FNF = FNU - IFNU
  TFNF = FNF + FNF
  IDUM = 0
!     S14ABE assumed not to fail, therefore IDUM set to zero.
  BK = S14ABE(FKK+TFNF+1.0E0,IDUM) - S14ABE(FKK+1.0E0,IDUM) - &
       S14ABE(TFNF+1.0E0,IDUM)
  BK = EXP(BK)
  SUM = CZERO
  KM = KK - INU
  DO 100 I = 1, KM
   PT = P2
   P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
   P1 = PT
   AK = 1.0E0 - TFNF/(FKK+TFNF)
   ACK = BK*AK
   SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
   BK = ACK
   FKK = FKK - 1.0E0
  100 continue
  Y(N) = P2
  if (N /= 1) then
   DO 120 I = 2, N
      PT = P2
      P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
      P1 = PT
      AK = 1.0E0 - TFNF/(FKK+TFNF)
      ACK = BK*AK
      SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
      BK = ACK
      FKK = FKK - 1.0E0
      M = N - I + 1
      Y(M) = P2
  120    continue
  endif
  if (IFNU > 0) then
   DO 140 I = 1, IFNU
      PT = P2
      P2 = P1 + CMPLX(FKK+FNF,0.0E0)*RZ*P2
      P1 = PT
      AK = 1.0E0 - TFNF/(FKK+TFNF)
      ACK = BK*AK
      SUM = SUM + CMPLX(ACK+BK,0.0E0)*P1
      BK = ACK
      FKK = FKK - 1.0E0
  140    continue
  endif
  PT = Z
  if (KODE == 2) PT = PT - CMPLX(X,0.0E0)
  P1 = -CMPLX(FNF,0.0E0)*LOG(RZ) + PT
  IDUM = 0
!     S14ABE assumed not to fail, therefore IDUM set to zero.
  AP = S14ABE(1.0E0+FNF,IDUM)
  PT = P1 - CMPLX(AP,0.0E0)
!     ------------------------------------------------------------------
!     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
!     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
!     ------------------------------------------------------------------
  P2 = P2 + SUM
  AP = ABS(P2)
  P1 = CMPLX(1.0E0/AP,0.0E0)
!      CK = EXP(PT)*P1
  IFL = 1
  CK = S01EAE(PT,IFL)*P1
  if ((IFL >= 1 .and. IFL <= 3) .or. IFL == 5) goto 200
  PT = CONJG(P2)*P1
  CNORM = CK*PT
  DO 160 I = 1, N
   Y(I) = Y(I)*CNORM
  160 continue
  return
  180 NZ = -2
  return
  200 NZ = -3
  return
  END
  subroutine DGUS17(Z,CSH,CCH)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-774 (DEC 1989).
!
!     Original name: CSHCH
!
!     DGUS17 COMPUTES THE COMPLEX HYPERBOLIC functionS CSH=SINH(X+I*Y)
!     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
!
!     .. Scalar Arguments ..
  COMPLEX           CCH, CSH, Z
!     .. Local Scalars ..
  REAL              CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, X, Y
!     .. Intrinsic functions ..
  INTRINSIC         AIMAG, CMPLX, COS, COSH, REAL, SIN, SINH
!     .. Executable Statements ..
!
  X = REAL(Z)
  Y = AIMAG(Z)
  SH = SINH(X)
  CH = COSH(X)
  SN = SIN(Y)
  CN = COS(Y)
  CSHR = SH*CN
  CSHI = CH*SN
  CSH = CMPLX(CSHR,CSHI)
  CCHR = CH*CN
  CCHI = SH*SN
  CCH = CMPLX(CCHR,CCHI)
  return
  END
  subroutine DGVS17(Y,NZ,ASCLE,TOL)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-775 (DEC 1989).
!
!     Original name: CUCHK
!
!      Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
!      EXP(-ALIM)=ASCLE=1.0E+3*X02AME()/TOL. THE TEST IS MADE TO SEE
!      IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW
!      WHEN Y IS SCALED (BY TOL) TO ITS PROPER VALUE. Y IS ACCEPTED
!      IF THE UNDERFLOW IS AT LEAST ONE PRECISION BELOW THE MAGNITUDE
!      OF THE LARGEST COMPONENT; OTHERWISE THE PHASE ANGLE DOES NOT HAVE
!      ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.
!
!     .. Scalar Arguments ..
  COMPLEX           Y
  REAL              ASCLE, TOL
  INTEGER           NZ
!     .. Local Scalars ..
  REAL              SS, ST, YI, YR
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, MAX, MIN, REAL
!     .. Executable Statements ..
!
  NZ = 0
  YR = REAL(Y)
  YI = AIMAG(Y)
  YR = ABS(YR)
  YI = ABS(YI)
  ST = MIN(YR,YI)
  if (ST <= ASCLE) then
   SS = MAX(YR,YI)
   ST = ST/TOL
   if (SS < ST) NZ = 1
  endif
  return
  END
  subroutine DGWS17(ZR,FNU,N,Y,NZ,RZ,ASCLE,TOL,ELIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-776 (DEC 1989).
!
!     Original name: CKSCL
!
!     SET K functionS TO ZERO ON UNDERFLOW, continue RECURRENCE
!     ON SCALED functionS UNTIL TWO MEMBERS COME ON SCALE, THEN
!     return WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
!
!     .. Scalar Arguments ..
  COMPLEX           RZ, ZR
  REAL              ASCLE, ELIM, FNU, TOL
  INTEGER           N, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           CELM, CK, CS, CZERO, S1, S2, ZD
  REAL              AA, ACS, ALAS, AS, CSI, CSR, ELM, FN, HELIM, XX, &
                    ZRI
  INTEGER           I, IC, K, KK, NN, NW
!     .. Local Arrays ..
  COMPLEX           CY(2)
!     .. External subroutines ..
  EXTERNAL          DGVS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, COS, EXP, LOG, MIN, REAL, SIN
!     .. Data statements ..
  DATA              CZERO/(0.0E0,0.0E0)/
!     .. Executable Statements ..
!
  NZ = 0
  IC = 0
  XX = REAL(ZR)
  NN = MIN(2,N)
  DO 20 I = 1, NN
   S1 = Y(I)
   CY(I) = S1
   AS = ABS(S1)
   ACS = -XX + LOG(AS)
   NZ = NZ + 1
   Y(I) = CZERO
   if (ACS >= (-ELIM)) then
      CS = -ZR + LOG(S1)
      CSR = REAL(CS)
      CSI = AIMAG(CS)
      AA = EXP(CSR)/TOL
      CS = CMPLX(AA,0.0E0)*CMPLX(COS(CSI),SIN(CSI))
      CALL DGVS17(CS,NW,ASCLE,TOL)
      if (NW == 0) then
         Y(I) = CS
         NZ = NZ - 1
         IC = I
      endif
   endif
   20 continue
  if (N /= 1) then
   if (IC <= 1) then
      Y(1) = CZERO
      NZ = 2
   endif
   if (N /= 2) then
      if (NZ /= 0) then
         FN = FNU + 1.0E0
         CK = CMPLX(FN,0.0E0)*RZ
         S1 = CY(1)
         S2 = CY(2)
         HELIM = 0.5E0*ELIM
         ELM = EXP(-ELIM)
         CELM = CMPLX(ELM,0.0E0)
         ZRI = AIMAG(ZR)
         ZD = ZR
!
!              FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE
!              RECURRENCE IF S2 GETS LARGER THAN EXP(ELIM/2)
!
         DO 40 I = 3, N
            KK = I
            CS = S2
            S2 = CK*S2 + S1
            S1 = CS
            CK = CK + RZ
            AS = ABS(S2)
            ALAS = LOG(AS)
            ACS = -XX + ALAS
            NZ = NZ + 1
            Y(I) = CZERO
            if (ACS >= (-ELIM)) then
               CS = -ZD + LOG(S2)
               CSR = REAL(CS)
               CSI = AIMAG(CS)
               AA = EXP(CSR)/TOL
               CS = CMPLX(AA,0.0E0)*CMPLX(COS(CSI),SIN(CSI))
               CALL DGVS17(CS,NW,ASCLE,TOL)
               if (NW == 0) then
                  Y(I) = CS
                  NZ = NZ - 1
                  if (IC == (KK-1)) then
                     goto 60
                  ELSE
                     IC = KK
                     goto 40
                  endif
               endif
            endif
            if (ALAS >= HELIM) then
               XX = XX - ELIM
               S1 = S1*CELM
               S2 = S2*CELM
               ZD = CMPLX(XX,ZRI)
            endif
   40          continue
         NZ = N
         if (IC == N) NZ = N - 1
         goto 80
   60          NZ = KK - 2
   80          DO 100 K = 1, NZ
            Y(K) = CZERO
  100          continue
      endif
   endif
  endif
  return
  END
  subroutine DGXS17(Z,FNU,KODE,N,Y,NZ,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-777 (DEC 1989).
!
!     Original name: CBKNU
!
!     DGXS17 COMPUTES THE K BESSEL function IN THE RIGHT HALF Z PLANE
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, TOL
  INTEGER           KODE, N, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           CCH, CELM, CK, COEF, CONE, CRSC, CS, CSCL, CSH, &
                    CTWO, CZ, CZERO, F, FMU, P, P1, P2, PT, Q, RZ, &
                    S1, S2, SMU, ST, ZD
  REAL              A1, A2, AA, AK, ALAS, AS, ASCLE, BB, BK, CAZ, &
                    DNU, DNU2, ELM, ETEST, FC, FHS, FK, FKS, FPI, &
                    G1, G2, HELIM, HPI, P2I, P2M, P2R, PI, R1, RK, &
                    RTHPI, S, SPI, T1, T2, TM, TTH, XD, XX, YD, YY
  INTEGER           I, IC, IDUM, IFL, IFLAG, INU, INUB, J, K, KFLAG, &
                    KK, KMAX, KODED, NW
!     .. Local Arrays ..
  COMPLEX           CSR(3), CSS(3), CY(2)
  REAL              BRY(3), CC(8)
!     .. External functions ..
  COMPLEX           S01EAE
  REAL              S14ABE, X02AME, X02ALE
  INTEGER           X02BHE, X02BJE
  EXTERNAL          S14ABE, S01EAE, X02AME, X02ALE, X02BHE, X02BJE
!     .. External subroutines ..
  EXTERNAL          DGUS17, DGVS17, DGWS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, ATAN, CMPLX, CONJG, COS, EXP, INT, &
                    LOG, LOG10, MAX, MIN, REAL, SIN, SQRT
!     .. Data statements ..
!
!
!
  DATA              KMAX/30/
  DATA              R1/2.0E0/
  DATA              CZERO, CONE, CTWO/(0.0E0,0.0E0), (1.0E0,0.0E0), &
                    (2.0E0,0.0E0)/
  DATA              PI, RTHPI, SPI, HPI, FPI, &
                    TTH/3.14159265358979324E0, &
                    1.25331413731550025E0, 1.90985931710274403E0, &
                    1.57079632679489662E0, 1.89769999331517738E0, &
                    6.66666666666666666E-01/
  DATA              CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), &
                    CC(8)/5.77215664901532861E-01, &
                    -4.20026350340952355E-02, &
                    -4.21977345555443367E-02, &
                    7.21894324666309954E-03, &
                    -2.15241674114950973E-04, &
                    -2.01348547807882387E-05, &
                    1.13302723198169588E-06, &
                    6.11609510448141582E-09/
!     .. Executable Statements ..
!
  XX = REAL(Z)
  YY = AIMAG(Z)
  CAZ = ABS(Z)
  CSCL = CMPLX(1.0E0/TOL,0.0E0)
  CRSC = CMPLX(TOL,0.0E0)
  CSS(1) = CSCL
  CSS(2) = CONE
  CSS(3) = CRSC
  CSR(1) = CRSC
  CSR(2) = CONE
  CSR(3) = CSCL
  BRY(1) = (1.0E+3*X02AME())/TOL
  BRY(2) = 1.0E0/BRY(1)
  BRY(3) = X02ALE()
  NZ = 0
  IFLAG = 0
  KODED = KODE
  RZ = CTWO/Z
  INU = INT(FNU+0.5E0)
  DNU = FNU - INU
  if (ABS(DNU) /= 0.5E0) then
   DNU2 = 0.0E0
   if (ABS(DNU) > TOL) DNU2 = DNU*DNU
   if (CAZ <= R1) then
!           ------------------------------------------------------------
!           SERIES FOR CABS(Z) <= R1
!           ------------------------------------------------------------
      FC = 1.0E0
      SMU = LOG(RZ)
      FMU = SMU*CMPLX(DNU,0.0E0)
      CALL DGUS17(FMU,CSH,CCH)
      if (DNU /= 0.0E0) then
         FC = DNU*PI
         FC = FC/SIN(FC)
         SMU = CSH*CMPLX(1.0E0/DNU,0.0E0)
      endif
      A2 = 1.0E0 + DNU
!           ------------------------------------------------------------
!           GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU),
!           T2=1/GAM(1+DNU)
!           ------------------------------------------------------------
      IDUM = 0
!           S14ABE assumed not to fail, therefore IDUM set to zero.
      T2 = EXP(-S14ABE(A2,IDUM))
      T1 = 1.0E0/(T2*FC)
      if (ABS(DNU) > 0.1E0) then
         G1 = (T1-T2)/(DNU+DNU)
      ELSE
!              ---------------------------------------------------------
!              SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
!              ---------------------------------------------------------
         AK = 1.0E0
         S = CC(1)
         DO 20 K = 2, 8
            AK = AK*DNU2
            TM = CC(K)*AK
            S = S + TM
            if (ABS(TM) < TOL) goto 40
   20          continue
   40          G1 = -S
      endif
      G2 = 0.5E0*(T1+T2)*FC
      G1 = G1*FC
      F = CMPLX(G1,0.0E0)*CCH + SMU*CMPLX(G2,0.0E0)
      IFL = 1
      PT = S01EAE(FMU,IFL)
      if ((IFL >= 1 .and. IFL <= 3) .or. IFL == 5) goto 320
      P = CMPLX(0.5E0/T2,0.0E0)*PT
      Q = CMPLX(0.5E0/T1,0.0E0)/PT
      S1 = F
      S2 = P
      AK = 1.0E0
      A1 = 1.0E0
      CK = CONE
      BK = 1.0E0 - DNU2
      if (INU > 0 .or. N > 1) then
!              ---------------------------------------------------------
!              GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
!              ---------------------------------------------------------
         if (CAZ >= TOL) then
            CZ = Z*Z*CMPLX(0.25E0,0.0E0)
            T1 = 0.25E0*CAZ*CAZ
   60             continue
            F = (F*CMPLX(AK,0.0E0)+P+Q)*CMPLX(1.0E0/BK,0.0E0)
            P = P*CMPLX(1.0E0/(AK-DNU),0.0E0)
            Q = Q*CMPLX(1.0E0/(AK+DNU),0.0E0)
            RK = 1.0E0/AK
            CK = CK*CZ*CMPLX(RK,0.0E0)
            S1 = S1 + CK*F
            S2 = S2 + CK*(P-F*CMPLX(AK,0.0E0))
            A1 = A1*T1*RK
            BK = BK + AK + AK + 1.0E0
            AK = AK + 1.0E0
            if (A1 > TOL) goto 60
         endif
         KFLAG = 2
         BK = REAL(SMU)
         A1 = FNU + 1.0E0
         AK = A1*ABS(BK)
         if (AK > ALIM) KFLAG = 3
         P2 = S2*CSS(KFLAG)
         S2 = P2*RZ
         S1 = S1*CSS(KFLAG)
         if (KODED /= 1) then
!                  F = EXP(Z)
            IFL = 1
            F = S01EAE(Z,IFL)
            if ((IFL >= 1 .and. IFL <= 3) .or. IFL == 5) goto 320
            S1 = S1*F
            S2 = S2*F
         endif
         goto 160
      ELSE
!              ---------------------------------------------------------
!              GENERATE K(FNU,Z), 0.0D0 <= FNU < 0.5D0 AND N=1
!              ---------------------------------------------------------
         if (CAZ >= TOL) then
            CZ = Z*Z*CMPLX(0.25E0,0.0E0)
            T1 = 0.25E0*CAZ*CAZ
   80             continue
            F = (F*CMPLX(AK,0.0E0)+P+Q)*CMPLX(1.0E0/BK,0.0E0)
            P = P*CMPLX(1.0E0/(AK-DNU),0.0E0)
            Q = Q*CMPLX(1.0E0/(AK+DNU),0.0E0)
            RK = 1.0E0/AK
            CK = CK*CZ*CMPLX(RK,0.0E0)
            S1 = S1 + CK*F
            A1 = A1*T1*RK
            BK = BK + AK + AK + 1.0E0
            AK = AK + 1.0E0
            if (A1 > TOL) goto 80
         endif
         Y(1) = S1
!               if (KODED /= 1) Y(1) = S1*EXP(Z)
         if (KODED /= 1) then
            IFL = 1
            Y(1) = S01EAE(Z,IFL)
            if ((IFL >= 1 .and. IFL <= 3) .or. IFL == 5) goto 320
            Y(1) = S1*Y(1)
         endif
         return
      endif
   endif
  endif
!     ------------------------------------------------------------------
!     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
!     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
!     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
!     RECURSION
!     ------------------------------------------------------------------
  COEF = CMPLX(RTHPI,0.0E0)/SQRT(Z)
  KFLAG = 2
  if (KODED /= 2) then
   if (XX > ALIM) then
!           ------------------------------------------------------------
!           SCALE BY EXP(Z), IFLAG = 1 CASES
!           ------------------------------------------------------------
      KODED = 2
      IFLAG = 1
      KFLAG = 2
   ELSE
!           BLANK LINE
!            A1 = EXP(-XX)*REAL(CSS(KFLAG))
!            PT = CMPLX(A1,0.0E0)*CMPLX(COS(YY),-SIN(YY))
      IFL = 1
      PT = S01EAE(CMPLX(-XX,-YY),IFL)
      if ((IFL >= 1 .and. IFL <= 3) .or. IFL == 5) goto 320
      PT = PT*REAL(CSS(KFLAG))
      COEF = COEF*PT
   endif
  endif
  if (ABS(DNU) /= 0.5E0) then
!        ---------------------------------------------------------------
!        MILLER ALGORITHM FOR CABS(Z)>R1
!        ---------------------------------------------------------------
   AK = COS(PI*DNU)
   AK = ABS(AK)
   if (AK /= 0.0E0) then
      FHS = ABS(0.25E0-DNU2)
      if (FHS /= 0.0E0) then
!              ---------------------------------------------------------
!              COMPUTE R2=F(E). IF CABS(Z) >= R2, USE FORWARD RECURRENCE
!              TO DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT
!              LINE ON 12 <= E <= 60. E IS COMPUTED FROM
!              2**(-E)=B**(1-X02BJE())=TOL WHERE B IS THE BASE OF THE
!              ARITHMETIC.
!              ---------------------------------------------------------
         T1 = (X02BJE()-1)*LOG10(REAL(X02BHE()))*3.321928094E0
         T1 = MAX(T1,12.0E0)
         T1 = MIN(T1,60.0E0)
         T2 = TTH*T1 - 6.0E0
         if (XX /= 0.0E0) then
            T1 = ATAN(YY/XX)
            T1 = ABS(T1)
         ELSE
            T1 = HPI
         endif
         if (T2 > CAZ) then
!                 ------------------------------------------------------
!                 COMPUTE BACKWARD INDEX K FOR CABS(Z) < R2
!                 ------------------------------------------------------
            A2 = SQRT(CAZ)
            AK = FPI*AK/(TOL*SQRT(A2))
            AA = 3.0E0*T1/(1.0E0+CAZ)
            BB = 14.7E0*T1/(28.0E0+CAZ)
            AK = (LOG(AK)+CAZ*COS(AA)/(1.0E0+0.008E0*CAZ))/COS(BB)
            FK = 0.12125E0*AK*AK/CAZ + 1.5E0
         ELSE
!                 ------------------------------------------------------
!                 FORWARD RECURRENCE LOOP WHEN CABS(Z) >= R2
!                 ------------------------------------------------------
            ETEST = AK/(PI*CAZ*TOL)
            FK = 1.0E0
            if (ETEST >= 1.0E0) then
               FKS = 2.0E0
               RK = CAZ + CAZ + 2.0E0
               A1 = 0.0E0
               A2 = 1.0E0
               DO 100 I = 1, KMAX
                  AK = FHS/FKS
                  BK = RK/(FK+1.0E0)
                  TM = A2
                  A2 = BK*A2 - AK*A1
                  A1 = TM
                  RK = RK + 2.0E0
                  FKS = FKS + FK + FK + 2.0E0
                  FHS = FHS + FK + FK
                  FK = FK + 1.0E0
                  TM = ABS(A2)*FK
                  if (ETEST < TM) goto 120
  100                continue
               NZ = -2
               return
  120                FK = FK + SPI*T1*SQRT(T2/CAZ)
               FHS = ABS(0.25E0-DNU2)
            endif
         endif
         K = INT(FK)
!              ---------------------------------------------------------
!              BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
!              ---------------------------------------------------------
         FK = K
         FKS = FK*FK
         P1 = CZERO
         P2 = CMPLX(TOL,0.0E0)
         CS = P2
         DO 140 I = 1, K
            A1 = FKS - FK
            A2 = (FKS+FK)/(A1+FHS)
            RK = 2.0E0/(FK+1.0E0)
            T1 = (FK+XX)*RK
            T2 = YY*RK
            PT = P2
            P2 = (P2*CMPLX(T1,T2)-P1)*CMPLX(A2,0.0E0)
            P1 = PT
            CS = CS + P2
            FKS = A1 - FK + 1.0E0
            FK = FK - 1.0E0
  140          continue
!              ---------------------------------------------------------
!              COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR
!              BETTER SCALING
!              ---------------------------------------------------------
         TM = ABS(CS)
         PT = CMPLX(1.0E0/TM,0.0E0)
         S1 = PT*P2
         CS = CONJG(CS)*PT
         S1 = COEF*S1*CS
         if (INU > 0 .or. N > 1) then
!                 ------------------------------------------------------
!                 COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR
!                 SCALING
!                 ------------------------------------------------------
            TM = ABS(P2)
            PT = CMPLX(1.0E0/TM,0.0E0)
            P1 = PT*P1
            P2 = CONJG(P2)*PT
            PT = P1*P2
            S2 = S1*(CONE+(CMPLX(DNU+0.5E0,0.0E0)-PT)/Z)
            goto 160
         ELSE
            ZD = Z
            if (IFLAG == 1) then
               goto 240
            ELSE
               goto 260
            endif
         endif
      endif
   endif
  endif
!     ------------------------------------------------------------------
!     FNU=HALF ODD INTEGER CASE, DNU=-0.5
!     ------------------------------------------------------------------
  S1 = COEF
  S2 = COEF
!     ------------------------------------------------------------------
!     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION WITH
!     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
!     ------------------------------------------------------------------
  160 continue
  CK = CMPLX(DNU+1.0E0,0.0E0)*RZ
  if (N == 1) INU = INU - 1
  if (INU > 0) then
   INUB = 1
   if (IFLAG == 1) then
!           ------------------------------------------------------------
!           IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON
!           UNDERFLOW
!           ------------------------------------------------------------
      HELIM = 0.5E0*ELIM
      ELM = EXP(-ELIM)
      CELM = CMPLX(ELM,0.0E0)
      ASCLE = BRY(1)
      ZD = Z
      XD = XX
      YD = YY
      IC = -1
      J = 2
      DO 180 I = 1, INU
         ST = S2
         S2 = CK*S2 + S1
         S1 = ST
         CK = CK + RZ
         AS = ABS(S2)
         ALAS = LOG(AS)
         P2R = -XD + ALAS
         if (P2R >= (-ELIM)) then
            P2 = -ZD + LOG(S2)
            P2R = REAL(P2)
            P2I = AIMAG(P2)
            P2M = EXP(P2R)/TOL
            P1 = CMPLX(P2M,0.0E0)*CMPLX(COS(P2I),SIN(P2I))
            CALL DGVS17(P1,NW,ASCLE,TOL)
            if (NW == 0) then
               J = 3 - J
               CY(J) = P1
               if (IC == (I-1)) then
                  goto 200
               ELSE
                  IC = I
                  goto 180
               endif
            endif
         endif
         if (ALAS >= HELIM) then
            XD = XD - ELIM
            S1 = S1*CELM
            S2 = S2*CELM
            ZD = CMPLX(XD,YD)
         endif
  180       continue
      if (N == 1) S1 = S2
      goto 240
  200       KFLAG = 1
      INUB = I + 1
      S2 = CY(J)
      J = 3 - J
      S1 = CY(J)
      if (INUB > INU) then
         if (N == 1) S1 = S2
         goto 260
      endif
   endif
   P1 = CSR(KFLAG)
   ASCLE = BRY(KFLAG)
   DO 220 I = INUB, INU
      ST = S2
      S2 = CK*S2 + S1
      S1 = ST
      CK = CK + RZ
      if (KFLAG < 3) then
         P2 = S2*P1
         P2R = REAL(P2)
         P2I = AIMAG(P2)
         P2R = ABS(P2R)
         P2I = ABS(P2I)
         P2M = MAX(P2R,P2I)
         if (P2M > ASCLE) then
            KFLAG = KFLAG + 1
            ASCLE = BRY(KFLAG)
            S1 = S1*P1
            S2 = P2
            S1 = S1*CSS(KFLAG)
            S2 = S2*CSS(KFLAG)
            P1 = CSR(KFLAG)
         endif
      endif
  220    continue
   if (N == 1) S1 = S2
   goto 260
  ELSE
   if (N == 1) S1 = S2
   ZD = Z
   if (IFLAG /= 1) goto 260
  endif
  240 Y(1) = S1
  if (N /= 1) Y(2) = S2
  ASCLE = BRY(1)
  CALL DGWS17(ZD,FNU,N,Y,NZ,RZ,ASCLE,TOL,ELIM)
  INU = N - NZ
  if (INU <= 0) then
   return
  ELSE
   KK = NZ + 1
   S1 = Y(KK)
   Y(KK) = S1*CSR(1)
   if (INU == 1) then
      return
   ELSE
      KK = NZ + 2
      S2 = Y(KK)
      Y(KK) = S2*CSR(1)
      if (INU == 2) then
         return
      ELSE
         T2 = FNU + KK - 1
         CK = CMPLX(T2,0.0E0)*RZ
         KFLAG = 1
         goto 280
      endif
   endif
  endif
  260 Y(1) = S1*CSR(KFLAG)
  if (N == 1) then
   return
  ELSE
   Y(2) = S2*CSR(KFLAG)
   if (N == 2) then
      return
   ELSE
      KK = 2
   endif
  endif
  280 KK = KK + 1
  if (KK <= N) then
   P1 = CSR(KFLAG)
   ASCLE = BRY(KFLAG)
   DO 300 I = KK, N
      P2 = S2
      S2 = CK*S2 + S1
      S1 = P2
      CK = CK + RZ
      P2 = S2*P1
      Y(I) = P2
      if (KFLAG < 3) then
         P2R = REAL(P2)
         P2I = AIMAG(P2)
         P2R = ABS(P2R)
         P2I = ABS(P2I)
         P2M = MAX(P2R,P2I)
         if (P2M > ASCLE) then
            KFLAG = KFLAG + 1
            ASCLE = BRY(KFLAG)
            S1 = S1*P1
            S2 = P2
            S1 = S1*CSS(KFLAG)
            S2 = S2*CSS(KFLAG)
            P1 = CSR(KFLAG)
         endif
      endif
  300    continue
  endif
  return
  320 NZ = -3
  return
  END
  subroutine DGYS17(Z,FNU,KODE,N,Y,NZ,RL,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-778 (DEC 1989).
!
!     Original name: CASYI
!
!     DGYS17 COMPUTES THE I BESSEL function FOR REAL(Z) >= 0.0 BY
!     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z) IN THE
!     REGION CABS(Z)>MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL return.
!     NZ < 0 INDICATES AN OVERFLOW ON KODE=1.
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, RL, TOL
  INTEGER           KODE, N, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           AK1, CK, CONE, CS1, CS2, CZ, CZERO, DK, EZ, P1, &
                    RZ, S2
  REAL              AA, ACZ, AEZ, AK, ARG, ARM, ATOL, AZ, BB, BK, &
                    DFNU, DNU2, FDN, PI, RTPI, RTR1, S, SGN, SQK, X, &
                    YY
  INTEGER           I, IB, IERR1, IL, INU, J, JL, K, KODED, M, NN
!     .. External functions ..
  COMPLEX           S01EAE
  REAL              X02AME
  EXTERNAL          S01EAE, X02AME
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, COS, EXP, INT, MIN, MOD, &
                    REAL, SIN, SQRT
!     .. Data statements ..
  DATA              PI, RTPI/3.14159265358979324E0, &
                    0.159154943091895336E0/
  DATA              CZERO, CONE/(0.0E0,0.0E0), (1.0E0,0.0E0)/
!     .. Executable Statements ..
!
  NZ = 0
  AZ = ABS(Z)
  X = REAL(Z)
  ARM = 1.0E+3*X02AME()
  RTR1 = SQRT(ARM)
  IL = MIN(2,N)
  DFNU = FNU + N - IL
!     ------------------------------------------------------------------
!     OVERFLOW TEST
!     ------------------------------------------------------------------
  AK1 = CMPLX(RTPI,0.0E0)/Z
  AK1 = SQRT(AK1)
  CZ = Z
  if (KODE == 2) CZ = Z - CMPLX(X,0.0E0)
  ACZ = REAL(CZ)
  if (ABS(ACZ) > ELIM) then
   NZ = -1
  ELSE
   DNU2 = DFNU + DFNU
   KODED = 1
   if ((ABS(ACZ) <= ALIM) .or. (N <= 2)) then
      KODED = 0
      IERR1 = 1
      AK1 = AK1*S01EAE(CZ,IERR1)
!        Allow reduced precision from S01EAE, but disallow other errors.
      if ((IERR1 >= 1 .and. IERR1 <= 3) .or. IERR1 == 5) goto 140
   endif
   FDN = 0.0E0
   if (DNU2 > RTR1) FDN = DNU2*DNU2
   EZ = Z*CMPLX(8.0E0,0.0E0)
!        ---------------------------------------------------------------
!        WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO
!        THE FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF
!        THE EXPANSION FOR THE IMAGINARY PART.
!        ---------------------------------------------------------------
   AEZ = 8.0E0*AZ
   S = TOL/AEZ
   JL = INT(RL+RL) + 2
   YY = AIMAG(Z)
   P1 = CZERO
   if (YY /= 0.0E0) then
!           ------------------------------------------------------------
!           CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
!           SIGNIFICANCE WHEN FNU OR N IS LARGE
!           ------------------------------------------------------------
      INU = INT(FNU)
      ARG = (FNU-INU)*PI
      INU = INU + N - IL
      AK = -SIN(ARG)
      BK = COS(ARG)
      if (YY < 0.0E0) BK = -BK
      P1 = CMPLX(AK,BK)
      if (MOD(INU,2) == 1) P1 = -P1
   endif
   DO 60 K = 1, IL
      SQK = FDN - 1.0E0
      ATOL = S*ABS(SQK)
      SGN = 1.0E0
      CS1 = CONE
      CS2 = CONE
      CK = CONE
      AK = 0.0E0
      AA = 1.0E0
      BB = AEZ
      DK = EZ
      DO 20 J = 1, JL
         CK = CK*CMPLX(SQK,0.0E0)/DK
         CS2 = CS2 + CK
         SGN = -SGN
         CS1 = CS1 + CK*CMPLX(SGN,0.0E0)
         DK = DK + EZ
         AA = AA*ABS(SQK)/BB
         BB = BB + AEZ
         AK = AK + 8.0E0
         SQK = SQK - AK
         if (AA <= ATOL) goto 40
   20       continue
      goto 120
   40       S2 = CS1
      if (X+X < ELIM) then
         IERR1 = 1
         S2 = S2 + P1*CS2*S01EAE(-Z-Z,IERR1)
         if ((IERR1 >= 1 .and. IERR1 <= 3) .or. IERR1 == 5) &
               goto 140
      endif
      FDN = FDN + 8.0E0*DFNU + 4.0E0
      P1 = -P1
      M = N - IL + K
      Y(M) = S2*AK1
   60    continue
   if (N > 2) then
      NN = N
      K = NN - 2
      AK = K
      RZ = (CONE+CONE)/Z
      IB = 3
      DO 80 I = IB, NN
         Y(K) = CMPLX(AK+FNU,0.0E0)*RZ*Y(K+1) + Y(K+2)
         AK = AK - 1.0E0
         K = K - 1
   80       continue
      if (KODED /= 0) then
         IERR1 = 1
         CK = S01EAE(CZ,IERR1)
         if ((IERR1 >= 1 .and. IERR1 <= 3) .or. IERR1 == 5) &
               goto 140
         DO 100 I = 1, NN
            Y(I) = Y(I)*CK
  100          continue
      endif
   endif
   return
  120    NZ = -2
   return
  140    NZ = -3
  endif
  return
  END
  subroutine DGZS17(Z,FNU,KODE,MR,N,Y,NZ,RL,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-779 (DEC 1989).
!
!     Original name: CACAI
!
!     DGZS17 APPLIES THE ANALYTIC CONTINUATION FORMULA
!
!         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!                 MP=PI*MR*CMPLX(0.0,1.0)
!
!     TO continue THE K function FROM THE RIGHT HALF TO THE LEFT
!     HALF Z PLANE FOR USE WITH S17DGE WHERE FNU=1/3 OR 2/3 AND N=1.
!     DGZS17 IS THE SAME AS DLZS17 WITH THE PARTS FOR LARGER ORDERS AND
!     RECURRENCE REMOVED. A RECURSIVE CALL TO DLZS17 CAN RESULT IF S17DL
!     IS CALLED FROM S17DGE.
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, RL, TOL
  INTEGER           KODE, MR, N, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           C1, C2, CSGN, CSPN, ZN
  REAL              ARG, ASCLE, AZ, CPN, DFNU, FMR, PI, SGN, SPN, YY
  INTEGER           INU, IUF, NN, NW
!     .. Local Arrays ..
  COMPLEX           CY(2)
!     .. External functions ..
  REAL              X02AME
  EXTERNAL          X02AME
!     .. External subroutines ..
  EXTERNAL          DGRS17, DGSS17, DGTS17, DGXS17, DGYS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, COS, INT, MOD, SIGN, SIN
!     .. Data statements ..
  DATA              PI/3.14159265358979324E0/
!     .. Executable Statements ..
!
  NZ = 0
  ZN = -Z
  AZ = ABS(Z)
  NN = N
  DFNU = FNU + N - 1
  if (AZ > 2.0E0) then
   if (AZ*AZ*0.25E0 > DFNU+1.0E0) then
      if (AZ < RL) then
!              ---------------------------------------------------------
!              MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I
!              function
!              ---------------------------------------------------------
         CALL DGTS17(ZN,FNU,KODE,NN,Y,NW,TOL)
         if (NW < 0) then
            goto 40
         ELSE
            goto 20
         endif
      ELSE
!              ---------------------------------------------------------
!              ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I function
!              ---------------------------------------------------------
         CALL DGYS17(ZN,FNU,KODE,NN,Y,NW,RL,TOL,ELIM,ALIM)
         if (NW < 0) then
            goto 40
         ELSE
            goto 20
         endif
      endif
   endif
  endif
!     ------------------------------------------------------------------
!     POWER SERIES FOR THE I function
!     ------------------------------------------------------------------
  CALL DGRS17(ZN,FNU,KODE,NN,Y,NW,TOL,ELIM,ALIM)
!     ------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K function
!     ------------------------------------------------------------------
   20 CALL DGXS17(ZN,FNU,KODE,1,CY,NW,TOL,ELIM,ALIM)
  if (NW == 0) then
   FMR = MR
   SGN = -SIGN(PI,FMR)
   CSGN = CMPLX(0.0E0,SGN)
   if (KODE /= 1) then
      YY = -AIMAG(ZN)
      CPN = COS(YY)
      SPN = SIN(YY)
      CSGN = CSGN*CMPLX(CPN,SPN)
   endif
!        ---------------------------------------------------------------
!        CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!        WHEN FNU IS LARGE
!        ---------------------------------------------------------------
   INU = INT(FNU)
   ARG = (FNU-INU)*SGN
   CPN = COS(ARG)
   SPN = SIN(ARG)
   CSPN = CMPLX(CPN,SPN)
   if (MOD(INU,2) == 1) CSPN = -CSPN
   C1 = CY(1)
   C2 = Y(1)
   if (KODE /= 1) then
      IUF = 0
      ASCLE = (1.0E+3*X02AME())/TOL
      CALL DGSS17(ZN,C1,C2,NW,ASCLE,ALIM,IUF)
      NZ = NZ + NW
   endif
   Y(1) = CSPN*C1 + CSGN*C2
   return
  endif
   40 NZ = -1
  if (NW == (-2)) NZ = -2
  if (NW == (-3)) NZ = -3
  return
  END
  subroutine DLYS17(Z,FNU,KODE,MR,N,Y,NZ,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-782 (DEC 1989).
!
!     Original name: CBUNK
!
!     DLYS17 COMPUTES THE K BESSEL function FOR FNU>FNUL.
!     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
!     IN DCZS18 AND THE EXPANSION FOR H(2,FNU,Z) IN DCYS18
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, TOL
  INTEGER           KODE, MR, N, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  REAL              AX, AY, XX, YY
!     .. External subroutines ..
  EXTERNAL          DCYS18, DCZS18
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, REAL
!     .. Executable Statements ..
!
  NZ = 0
  XX = REAL(Z)
  YY = AIMAG(Z)
  AX = ABS(XX)*1.7321E0
  AY = ABS(YY)
  if (AY > AX) then
!        ---------------------------------------------------------------
!        ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!        APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I
!        AND HPI=PI/2
!        ---------------------------------------------------------------
   CALL DCYS18(Z,FNU,KODE,MR,N,Y,NZ,TOL,ELIM,ALIM)
  ELSE
!        ---------------------------------------------------------------
!        ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
!        -PI/3 <= ARG(Z) <= PI/3
!        ---------------------------------------------------------------
   CALL DCZS18(Z,FNU,KODE,MR,N,Y,NZ,TOL,ELIM,ALIM)
  endif
  return
  END
  subroutine DLZS17(Z,FNU,KODE,MR,N,Y,NZ,RL,FNUL,TOL,ELIM,ALIM)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-783 (DEC 1989).
!
!     Original name: CACON
!
!     DLZS17 APPLIES THE ANALYTIC CONTINUATION FORMULA
!
!         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!                 MP=PI*MR*CMPLX(0.0,1.0)
!
!     TO continue THE K function FROM THE RIGHT HALF TO THE LEFT
!     HALF Z PLANE
!
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              ALIM, ELIM, FNU, FNUL, RL, TOL
  INTEGER           KODE, MR, N, NZ
!     .. Array Arguments ..
  COMPLEX           Y(N)
!     .. Local Scalars ..
  COMPLEX           C1, C2, CK, CONE, CS, CSCL, CSCR, CSGN, CSPN, &
                    RZ, S1, S2, SC1, SC2, ST, ZN
  REAL              ARG, AS2, ASCLE, BSCLE, C1I, C1M, C1R, CPN, FMR, &
                    PI, SGN, SPN, YY
  INTEGER           I, INU, IUF, KFLAG, NN, NW
!     .. Local Arrays ..
  COMPLEX           CSR(3), CSS(3), CY(2)
  REAL              BRY(3)
!     .. External functions ..
  REAL              X02AME, X02ALE
  EXTERNAL          X02AME, X02ALE
!     .. External subroutines ..
  EXTERNAL          DEZS17, DGSS17, DGXS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, COS, INT, MAX, MIN, MOD, &
                    REAL, SIGN, SIN
!     .. Data statements ..
  DATA              PI/3.14159265358979324E0/
  DATA              CONE/(1.0E0,0.0E0)/
!     .. Executable Statements ..
!
  NZ = 0
  ZN = -Z
  NN = N
  CALL DEZS17(ZN,FNU,KODE,NN,Y,NW,RL,FNUL,TOL,ELIM,ALIM)
  if (NW >= 0) then
!        ---------------------------------------------------------------
!        ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K function
!        ---------------------------------------------------------------
   NN = MIN(2,N)
   CALL DGXS17(ZN,FNU,KODE,NN,CY,NW,TOL,ELIM,ALIM)
   if (NW == 0) then
      S1 = CY(1)
      FMR = MR
      SGN = -SIGN(PI,FMR)
      CSGN = CMPLX(0.0E0,SGN)
      if (KODE /= 1) then
         YY = -AIMAG(ZN)
         CPN = COS(YY)
         SPN = SIN(YY)
         CSGN = CSGN*CMPLX(CPN,SPN)
      endif
!           ------------------------------------------------------------
!           CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF
!           SIGNIFICANCE WHEN FNU IS LARGE
!           ------------------------------------------------------------
      INU = INT(FNU)
      ARG = (FNU-INU)*SGN
      CPN = COS(ARG)
      SPN = SIN(ARG)
      CSPN = CMPLX(CPN,SPN)
      if (MOD(INU,2) == 1) CSPN = -CSPN
      IUF = 0
      C1 = S1
      C2 = Y(1)
      ASCLE = (1.0E+3*X02AME())/TOL
      if (KODE /= 1) then
         CALL DGSS17(ZN,C1,C2,NW,ASCLE,ALIM,IUF)
         NZ = NZ + NW
         SC1 = C1
      endif
      Y(1) = CSPN*C1 + CSGN*C2
      if (N /= 1) then
         CSPN = -CSPN
         S2 = CY(2)
         C1 = S2
         C2 = Y(2)
         if (KODE /= 1) then
            CALL DGSS17(ZN,C1,C2,NW,ASCLE,ALIM,IUF)
            NZ = NZ + NW
            SC2 = C1
         endif
         Y(2) = CSPN*C1 + CSGN*C2
         if (N /= 2) then
            CSPN = -CSPN
            RZ = CMPLX(2.0E0,0.0E0)/ZN
            CK = CMPLX(FNU+1.0E0,0.0E0)*RZ
!                 ------------------------------------------------------
!                 SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON
!                 K functionS
!                 ------------------------------------------------------
            CSCL = CMPLX(1.0E0/TOL,0.0E0)
            CSCR = CMPLX(TOL,0.0E0)
            CSS(1) = CSCL
            CSS(2) = CONE
            CSS(3) = CSCR
            CSR(1) = CSCR
            CSR(2) = CONE
            CSR(3) = CSCL
            BRY(1) = ASCLE
            BRY(2) = 1.0E0/ASCLE
            BRY(3) = X02ALE()
            AS2 = ABS(S2)
            KFLAG = 2
            if (AS2 <= BRY(1)) then
               KFLAG = 1
            else if (AS2 >= BRY(2)) then
               KFLAG = 3
            endif
            BSCLE = BRY(KFLAG)
            S1 = S1*CSS(KFLAG)
            S2 = S2*CSS(KFLAG)
            CS = CSR(KFLAG)
            DO 20 I = 3, N
               ST = S2
               S2 = CK*S2 + S1
               S1 = ST
               C1 = S2*CS
               ST = C1
               C2 = Y(I)
               if (KODE /= 1) then
                  if (IUF >= 0) then
                     CALL DGSS17(ZN,C1,C2,NW,ASCLE,ALIM,IUF)
                     NZ = NZ + NW
                     SC1 = SC2
                     SC2 = C1
                     if (IUF == 3) then
                        IUF = -4
                        S1 = SC1*CSS(KFLAG)
                        S2 = SC2*CSS(KFLAG)
                        ST = SC2
                     endif
                  endif
               endif
               Y(I) = CSPN*C1 + CSGN*C2
               CK = CK + RZ
               CSPN = -CSPN
               if (KFLAG < 3) then
                  C1R = REAL(C1)
                  C1I = AIMAG(C1)
                  C1R = ABS(C1R)
                  C1I = ABS(C1I)
                  C1M = MAX(C1R,C1I)
                  if (C1M > BSCLE) then
                     KFLAG = KFLAG + 1
                     BSCLE = BRY(KFLAG)
                     S1 = S1*CS
                     S2 = ST
                     S1 = S1*CSS(KFLAG)
                     S2 = S2*CSS(KFLAG)
                     CS = CSR(KFLAG)
                  endif
               endif
   20             continue
         endif
      endif
      return
   endif
  endif
  NZ = -1
  if (NW == (-2)) NZ = -2
  if (NW == (-3)) NZ = -3
  return
  END
  INTEGER function P01ABE(IFAIL,IERROR,SRNAME,NREC,REC)
!     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
!     MARK 13 REVISED. IER-621 (APR 1988).
!     MARK 13B REVISED. IER-668 (AUG 1988).
!
!     P01ABE is the error-handling routine for the NAG Library.
!
!     P01ABE either returns the value of IERROR through the routine
!     name (soft failure), or terminates execution of the program
!     (hard failure). Diagnostic messages may be output.
!
!     If IERROR = 0 (successful exit from the calling routine),
!     the value 0 is returned through the routine name, and no
!     message is output
!
!     If IERROR is non-zero (abnormal exit from the calling routine),
!     the action taken depends on the value of IFAIL.
!
!     IFAIL =  1: soft failure, silent exit (i.e. no messages are
!                 output)
!     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
!     IFAIL =-13: soft failure, noisy exit but standard messages from
!                 P01ABE are suppressed
!     IFAIL =  0: hard failure, noisy exit
!
!     For compatibility with certain routines included before Mark 12
!     P01ABE also allows an alternative specification of IFAIL in which
!     it is regarded as a decimal integer with least significant digits
!     cba. Then
!
!     a = 0: hard failure  a = 1: soft failure
!     b = 0: silent exit   b = 1: noisy exit
!
!     except that hard failure now always implies a noisy exit.
!
!     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
!
!     .. Scalar Arguments ..
  INTEGER                 IERROR, IFAIL, NREC
  CHARACTER*(*)           SRNAME
!     .. Array Arguments ..
  CHARACTER*(*)           REC(*)
!     .. Local Scalars ..
  INTEGER                 I, NERR
  CHARACTER*72            MESS
!     .. External subroutines ..
  EXTERNAL                ABZP01, X04AAE, X04BAE
!     .. Intrinsic functions ..
  INTRINSIC               ABS, MOD
!     .. Executable Statements ..
  if (IERROR /= 0) then
!        Abnormal exit from calling routine
   if (IFAIL == -1 .or. IFAIL == 0 .or. IFAIL == -13 .or. &
         (IFAIL > 0 .and. MOD(IFAIL/10,10) /= 0)) then
!           Noisy exit
      CALL X04AAE(0,NERR)
      DO 20 I = 1, NREC
         CALL X04BAE(NERR,REC(I))
   20       continue
      if (IFAIL /= -13) then
         WRITE (MESS,FMT=99999) SRNAME, IERROR
         CALL X04BAE(NERR,MESS)
         if (ABS(MOD(IFAIL,10)) /= 1) then
!                 Hard failure
            CALL X04BAE(NERR, &
                       ' ** NAG hard failure - execution terminated' &
                          )
            CALL ABZP01
         ELSE
!                 Soft failure
            CALL X04BAE(NERR, &
                          ' ** NAG soft failure - control returned')
         endif
      endif
   endif
  endif
  P01ABE = IERROR
  return
!
  99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL', &
         ' =',I6)
  END
  COMPLEX function S01EAE(Z,IFAIL)
!     MARK 14 RELEASE. NAG COPYRIGHT 1989.
!     returns exp(Z) for complex Z.
!     .. Parameters ..
  REAL                    ONE, ZERO
  PARAMETER               (ONE=1.0E0,ZERO=0.0E0)
  CHARACTER*6             SRNAME
  PARAMETER               (SRNAME='S01EAE')
!     .. Scalar Arguments ..
  COMPLEX                 Z
  INTEGER                 IFAIL
!     .. Local Scalars ..
  REAL                    COSY, EXPX, LNSAFE, RECEPS, RESI, RESR, &
                          RTSAFS, SAFE, SAFSIN, SINY, X, XPLNCY, &
                          XPLNSY, Y
  INTEGER                 IER, NREC
  LOGICAL                 FIRST
!     .. Local Arrays ..
  CHARACTER*80            REC(2)
!     .. External functions ..
  REAL                    X02AHE, X02AJE, X02AME
  INTEGER                 P01ABE
  EXTERNAL                X02AHE, X02AJE, X02AME, P01ABE
!     .. Intrinsic functions ..
  INTRINSIC               ABS, AIMAG, CMPLX, COS, EXP, LOG, MIN, &
                          REAL, SIGN, SIN, SQRT
!     .. Save statement ..
  SAVE                    SAFE, LNSAFE, SAFSIN, RTSAFS, FIRST
!     .. Data statements ..
  DATA                    FIRST/.true./
!     .. Executable Statements ..
  if (FIRST) then
   FIRST = .false.
   SAFE = ONE/X02AME()
   LNSAFE = LOG(SAFE)
   RECEPS = ONE/X02AJE()
   SAFSIN = MIN(X02AHE(ONE),RECEPS)
   if (SAFSIN < RECEPS**0.75E0) then
!         Assume that SAFSIN is approximately sqrt(RECEPS), in which
!         case IFAIL=4 cannot occur.
      RTSAFS = SAFSIN
   ELSE
!         Set RTSAFS to the argument above which SINE and COSINE will
!         return results of less than half precision, assuming that
!         SAFSIN is approximately equal to RECEPS.
      RTSAFS = SQRT(SAFSIN)
   endif
  endif
  NREC = 0
  IER = 0
  X = REAL(Z)
  Y = AIMAG(Z)
  if (ABS(Y) > SAFSIN) then
   IER = 5
   NREC = 2
   WRITE (REC,FMT=99995) Z
   S01EAE = ZERO
  ELSE
   COSY = COS(Y)
   SINY = SIN(Y)
   if (X > LNSAFE) then
      if (COSY == ZERO) then
         RESR = ZERO
      ELSE
         XPLNCY = X + LOG(ABS(COSY))
         if (XPLNCY > LNSAFE) then
            IER = 1
            RESR = SIGN(SAFE,COSY)
         ELSE
            RESR = SIGN(EXP(XPLNCY),COSY)
         endif
      endif
      if (SINY == ZERO) then
         RESI = ZERO
      ELSE
         XPLNSY = X + LOG(ABS(SINY))
         if (XPLNSY > LNSAFE) then
            IER = IER + 2
            RESI = SIGN(SAFE,SINY)
         ELSE
            RESI = SIGN(EXP(XPLNSY),SINY)
         endif
      endif
   ELSE
      EXPX = EXP(X)
      RESR = EXPX*COSY
      RESI = EXPX*SINY
   endif
   S01EAE = CMPLX(RESR,RESI)
   if (IER == 3) then
      NREC = 2
      WRITE (REC,FMT=99997) Z
   else if (ABS(Y) > RTSAFS) then
      IER = 4
      NREC = 2
      WRITE (REC,FMT=99996) Z
   else if (IER == 1) then
      NREC = 2
      WRITE (REC,FMT=99999) Z
   else if (IER == 2) then
      NREC = 2
      WRITE (REC,FMT=99998) Z
   endif
  endif
  IFAIL = P01ABE(IFAIL,IER,SRNAME,NREC,REC)
  return
!
  99999 FORMAT (1X,'** Argument Z causes overflow in real part of result:' &
         ,/4X,'Z = (',1P,E13.5,',',E13.5,')')
  99998 FORMAT (1X,'** Argument Z causes overflow in imaginary part of r', &
         'esult:',/4X,'Z = (',1P,E13.5,',',E13.5,')')
  99997 FORMAT (1X,'** Argument Z causes overflow in both real and imagi', &
         'nary parts of result:',/4X,'Z = (',1P,E13.5,',',E13.5,')')
  99996 FORMAT (1X,'** The imaginary part of argument Z is so large that', &
         ' the result is',/4X,'accurate to less than half precisio', &
         'n: Z = (',1P,E13.5,',',E13.5,')')
  99995 FORMAT (1X,'** The imaginary part of argument Z is so large that', &
         ' the result has no',/4X,'precision: Z = (',1P,E13.5,',', &
         E13.5,')')
  END
  REAL function S14ABE(X,IFAIL)
!     MARK 8 RELEASE. NAG COPYRIGHT 1979.
!     MARK 11.5(F77) REVISED. (SEPT 1985.)
!        LNGAMMA(X) function
!        ABRAMOWITZ AND STEGUN  CH.6
!
!     **************************************************************
!
!     TO EXTRACT THE CORRECT CODE FOR A PARTICULAR MACHINE-RANGE,
!     ACTIVATE THE STATEMENTS CONTAINED IN COMMENTS BEGINNING  CDD ,
!     WHERE  DD  IS THE APPROXIMATE NUMBER OF SIGNIFICANT DECIMAL
!     DIGITS REPRESENTED BY THE MACHINE
!     DELETE THE ILLEGAL DUMMY STATEMENTS OF THE FORM
!     * EXPANSION (NNNN) *
!
!     ALSO INSERT APPROPRIATE DATA STATEMENTS TO DEFINE CONSTANTS
!     WHICH DEPEND ON THE RANGE OF NUMBERS REPRESENTED BY THE
!     MACHINE, RATHER THAN THE PRECISION (SUITABLE STATEMENTS FOR
!     SOME MACHINES ARE CONTAINED IN COMMENTS BEGINNING CRD WHERE
!     D IS A DIGIT WHICH SIMPLY DISTINGUISHES A GROUP OF MACHINES).
!     DELETE THE ILLEGAL DUMMY DATA STATEMENTS WITH VALUES WRITTEN
!     *VALUE*
!
!     **************************************************************
!
!        IMPLEMENTATION DEPENDENT CONSTANTS
!
!        if (X < XSMALL)GAMMA(X)=1/X
!             I.E.   XSMALL*EULGAM <= XRELPR
!        LNGAM(XVBIG)=GBIG <= XOVFLO
!        LNR2PI=LN(SQRT(2*PI))
!        if (X>XBIG)LNGAM(X)=(X-0.5)LN(X)-X+LNR2PI
!
!     .. Parameters ..
  CHARACTER*6          SRNAME
  PARAMETER            (SRNAME='S14ABE')
!     .. Scalar Arguments ..
  REAL                 X
  INTEGER              IFAIL
!     .. Local Scalars ..
  REAL                 G, GBIG, LNR2PI, T, XBIG, XSMALL, XVBIG, Y
  INTEGER              I, M
!     .. Local Arrays ..
  CHARACTER*1          P01REC(1)
!     .. External functions ..
  INTEGER              P01ABE
  EXTERNAL             P01ABE
!     .. Intrinsic functions ..
  INTRINSIC            LOG, REAL
!     .. Data statements ..
!08   DATA XSMALL,XBIG,LNR2PI/
!08  *1.0E-8,1.2E+3,9.18938533E-1/
!09   DATA XSMALL,XBIG,LNR2PI/
!09  *1.0E-9,4.8E+3,9.189385332E-1/
!12   DATA XSMALL,XBIG,LNR2PI/
!12  *1.0E-12,3.7E+5,9.189385332047E-1/
  DATA XSMALL,XBIG,LNR2PI/ &
  1.0E-15,6.8E+6,9.189385332046727E-1/
!17   DATA XSMALL,XBIG,LNR2PI/
!17  *1.0E-17,7.7E+7,9.18938533204672742E-1/
!19   DATA XSMALL,XBIG,LNR2PI/
!19  *1.0E-19,3.1E+8,9.189385332046727418E-1/
!
!     RANGE DEPENDENT CONSTANTS
! DK DK      DATA XVBIG,GBIG/4.81E+2461,2.72E+2465/
  DATA XVBIG,GBIG/4.08E+36,3.40E+38/
!     FOR IEEE SINGLE PRECISION
!R0   DATA XVBIG,GBIG/4.08E+36,3.40E+38/
!     FOR IBM 360/370 AND SIMILAR MACHINES
!R1   DATA XVBIG,GBIG/4.29E+73,7.231E+75/
!     FOR DEC10, HONEYWELL, UNIVAC 1100 (S.P.)
!R2   DATA XVBIG,GBIG/2.05E36,1.69E38/
!     FOR ICL 1900
!R3   DATA XVBIG,GBIG/3.39E+74,5.784E+76/
!     FOR CDC 7600/CYBER
!R4   DATA XVBIG,GBIG/1.72E+319,1.26E+322/
!     FOR UNIVAC 1100 (D.P.)
!R5   DATA XVBIG,GBIG/1.28E305,8.98E+307/
!     FOR IEEE DOUBLE PRECISION
!R7   DATA XVBIG,GBIG/2.54D+305,1.79D+308/
!     .. Executable Statements ..
  if (X > XSMALL) goto 20
!        VERY SMALL RANGE
  if (X <= 0.0) goto 160
  IFAIL = 0
  S14ABE = -LOG(X)
  goto 200
!
   20 if (X > 15.0) goto 120
!        MAIN SMALL X RANGE
  M = X
  T = X - FLOAT(M)
  M = M - 1
  G = 1.0
  if (M) 40, 100, 60
   40 G = G/X
  goto 100
   60 DO 80 I = 1, M
   G = (X-FLOAT(I))*G
   80 continue
  100 T = 2.0*T - 1.0
!
!      * EXPANSION (0026) *
!
!     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 08E.09
!08   Y = (((((((((((+1.88278283E-6*T-5.48272091E-6)*T+1.03144033E-5)
!08  *    *T-3.13088821E-5)*T+1.01593694E-4)*T-2.98340924E-4)
!08  *    *T+9.15547391E-4)*T-2.42216251E-3)*T+9.04037536E-3)
!08  *    *T-1.34119055E-2)*T+1.03703361E-1)*T+1.61692007E-2)*T +
!08  *    8.86226925E-1
!
!     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 09E.10
!09   Y = ((((((((((((-6.463247484E-7*T+1.882782826E-6)
!09  *    *T-3.382165478E-6)*T+1.031440334E-5)*T-3.393457634E-5)
!09  *    *T+1.015936944E-4)*T-2.967655076E-4)*T+9.155473906E-4)
!09  *    *T-2.422622002E-3)*T+9.040375355E-3)*T-1.341184808E-2)
!09  *    *T+1.037033609E-1)*T+1.616919866E-2)*T + 8.862269255E-1
!
!     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 12E.13
!12   Y = ((((((((((((((((-8.965837291520E-9*T+2.612707393536E-8)
!12  *    *T-3.802866827264E-8)*T+1.173294768947E-7)
!12  *    *T-4.275076254106E-7)*T+1.276176602829E-6)
!12  *    *T-3.748495971011E-6)*T+1.123829871408E-5)
!12  *    *T-3.364018663166E-5)*T+1.009331480887E-4)
!12  *    *T-2.968895120407E-4)*T+9.157850115110E-4)
!12  *    *T-2.422595461409E-3)*T+9.040335037321E-3)
!12  *    *T-1.341185056618E-2)*T+1.037033634184E-1)
!12  *    *T+1.616919872437E-2)*T + 8.862269254528E-1
!
!     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 15E.16
  Y = (((((((((((((((-1.243191705600000E-10*T+ &
      3.622882508800000E-10)*T-4.030909644800000E-10) &
      *T+1.265236705280000E-9)*T-5.419466096640000E-9) &
      *T+1.613133578240000E-8)*T-4.620920340480000E-8) &
      *T+1.387603440435200E-7)*T-4.179652784537600E-7) &
      *T+1.253148247777280E-6)*T-3.754930502328320E-6) &
      *T+1.125234962812416E-5)*T-3.363759801664768E-5) &
      *T+1.009281733953869E-4)*T-2.968901194293069E-4) &
      *T+9.157859942174304E-4)*T-2.422595384546340E-3
  Y = ((((Y*T+9.040334940477911E-3)*T-1.341185057058971E-2) &
      *T+1.037033634220705E-1)*T+1.616919872444243E-2)*T + &
      8.862269254527580E-1
!
!     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 17E.18
!17   Y = (((((((((((((((-1.46381209600000000E-11*T+
!17  *    4.26560716800000000E-11)*T-4.01499750400000000E-11)
!17  *    *T+1.27679856640000000E-10)*T-6.13513953280000000E-10)
!17  *    *T+1.82243164160000000E-9)*T-5.11961333760000000E-9)
!17  *    *T+1.53835215257600000E-8)*T-4.64774927155200000E-8)
!17  *    *T+1.39383522590720000E-7)*T-4.17808776355840000E-7)
!17  *    *T+1.25281466396672000E-6)*T-3.75499034136576000E-6)
!17  *    *T+1.12524642975590400E-5)*T-3.36375833240268800E-5)
!17  *    *T+1.00928148823365120E-4)*T-2.96890121633200000E-4
!17   Y = ((((((Y*T+9.15785997288933120E-4)*T-2.42259538436268176E-3)
!17  *    *T+9.04033494028101968E-3)*T-1.34118505705967765E-2)
!17  *    *T+1.03703363422075456E-1)*T+1.61691987244425092E-2)*T +
!17  *    8.86226925452758013E-1
!
!     EXPANSION (0026) EVALUATED AS Y(T)  --PRECISION 19E.19
!19   Y = (((((((((((((((+6.710886400000000000E-13*T-
!19  *    1.677721600000000000E-12)*T+6.710886400000000000E-13)
!19  *    *T-4.152360960000000000E-12)*T+2.499805184000000000E-11)
!19  *    *T-6.898581504000000000E-11)*T+1.859597107200000000E-10)
!19  *    *T-5.676387532800000000E-10)*T+1.725556326400000000E-9)
!19  *    *T-5.166307737600000000E-9)*T+1.548131827712000000E-8)
!19  *    *T-4.644574052352000000E-8)*T+1.393195837030400000E-7)
!19  *    *T-4.178233990758400000E-7)*T+1.252842254950400000E-6)
!19  *    *T-3.754985815285760000E-6)*T+1.125245651030528000E-5
!19   Y = (((((((((Y*T-3.363758423922688000E-5)
!19  *    *T+1.009281502108083200E-4)
!19  *    *T-2.968901215188000000E-4)*T+9.157859971435078400E-4)
!19  *    *T-2.422595384370689760E-3)*T+9.040334940288877920E-3)
!19  *    *T-1.341185057059651648E-2)*T+1.037033634220752902E-1)
!19  *    *T+1.616919872444250674E-2)*T + 8.862269254527580137E-1
!
  S14ABE = LOG(Y*G)
  IFAIL = 0
  goto 200
!
  120 if (X > XBIG) goto 140
!        MAIN LARGE X RANGE
  T = 450.0/(X*X) - 1.0
!
!      * EXPANSION (0059) *
!
!     EXPANSION (0059) EVALUATED AS Y(T)  --PRECISION 08E.09
!08   Y = (+3.89980902E-9*T-6.16502533E-6)*T + 8.33271644E-2
!
!     EXPANSION (0059) EVALUATED AS Y(T)  --PRECISION 09E.10
!09   Y = (+3.899809019E-9*T-6.165025333E-6)*T + 8.332716441E-2
!
!     EXPANSION (0059) EVALUATED AS Y(T)  --PRECISION 12E.13
!12   Y = ((-6.451144077930E-12*T+3.899809018958E-9)
!12  *    *T-6.165020494506E-6)*T + 8.332716440658E-2
!
!     EXPANSION (0059) EVALUATED AS Y(T)  --PRECISION 15E.16
  Y = (((+2.002019273379824E-14*T-6.451144077929628E-12) &
      *T+3.899788998764847E-9)*T-6.165020494506090E-6)*T + &
      8.332716440657866E-2
!
!     EXPANSION (0059) EVALUATED AS Y(T)  --PRECISION 17E.18
!17   Y = ((((-9.94561064728159347E-17*T+2.00201927337982364E-14)
!17  *    *T-6.45101975779653651E-12)*T+3.89978899876484712E-9)
!17  *    *T-6.16502049453716986E-6)*T + 8.33271644065786580E-2
!
!     EXPANSION (0059) EVALUATED AS Y(T)  --PRECISION 19E.19
!19   Y = (((((+7.196406678180202240E-19*T-9.945610647281593472E-17)
!19  *    *T+2.001911327279650935E-14)*T-6.451019757796536510E-12)
!19  *    *T+3.899788999169644998E-9)*T-6.165020494537169862E-6)*T +
!19  *    8.332716440657865795E-2
!
  S14ABE = (X-0.5)*LOG(X) - X + LNR2PI + Y/X
  IFAIL = 0
  goto 200
!
  140 if (X > XVBIG) goto 180
!        ASYMPTOTIC LARGE X RANGE
  S14ABE = (X-0.5)*LOG(X) - X + LNR2PI
  IFAIL = 0
  goto 200
!
!        FAILURE EXITS
  160 IFAIL = P01ABE(IFAIL,1,SRNAME,0,P01REC)
  S14ABE = 0.0
  goto 200
  180 IFAIL = P01ABE(IFAIL,2,SRNAME,0,P01REC)
  S14ABE = GBIG
!
  200 return
!
  END
  subroutine S17DGE(DERIV,Z,SCALE,AI,NZ,IFAIL)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-770 (DEC 1989).
!
!     Original name: CAIRY
!
!     PURPOSE  TO COMPUTE AIRY functionS AI(Z) AND DAI(Z) FOR COMPLEX Z
!
!     DESCRIPTION
!     ===========
!
!         ON SCALE='U', S17DGE COMPUTES THE COMPLEX AIRY function AI(Z)
!         OR ITS DERIVATIVE DAI(Z)/DZ ON DERIV='F' OR DERIV='D'
!         RESPECTIVELY. ON SCALE='S', A SCALING OPTION
!         CEXP(ZTA)*AI(Z) OR CEXP(ZTA)*DAI(Z)/DZ IS PROVIDED TO REMOVE
!         THE EXPONENTIAL DECAY IN -PI/3 < ARG(Z) < PI/3 AND THE
!         EXPONENTIAL GROWTH IN PI/3 < ABS(ARG(Z)) < PI WHERE
!         ZTA=(2/3)*Z*CSQRT(Z)
!
!         WHILE THE AIRY functionS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
!         THE WHOLE Z PLANE, THE CORRESPONDING SCALED functionS DEFINED
!         FOR SCALE='S' HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
!         DEFINITIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
!         MATHEMATICAL functionS (REF. 1).
!
!         INPUT
!           Z      - Z=CMPLX(X,Y)
!           DERIV  - return function (DERIV='F') OR DERIVATIVE
!                    (DERIV='D')
!           SCALE  - A PARAMETER TO INDICATE THE SCALING OPTION
!                    SCALE = 'U' OR 'u' returnS
!                             AI=AI(Z)                ON DERIV='F' OR
!                             AI=DAI(Z)/DZ            ON DERIV='D'
!                    SCALE = 'S' OR 's' returnS
!                             AI=CEXP(ZTA)*AI(Z)      ON DERIV='F' OR
!                             AI=CEXP(ZTA)*DAI(Z)/DZ  ON DERIV='D' WHERE
!                             ZTA=(2/3)*Z*CSQRT(Z)
!
!         OUTPUT
!           AI     - COMPLEX ANSWER DEPENDING ON THE CHOICES FOR DERIV
!                    AND SCALE
!           NZ     - UNDERFLOW INDICATOR
!                    NZ= 0   , NORMAL return
!                    NZ= 1   , AI=CMPLX(0.0,0.0) DUE TO UNDERFLOW IN
!                              -PI/3 < ARG(Z) < PI/3 ON SCALE='U'
!           IFAIL  - ERROR FLAG
!                   IFAIL=0, NORMAL return - COMPUTATION COMPLETED
!                   IFAIL=1, INPUT ERROR   - NO COMPUTATION
!                   IFAIL=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA)
!                            TOO LARGE WITH SCALE = 'U'
!                   IFAIL=3, CABS(Z) LARGE      - COMPUTATION COMPLETED
!                            LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION
!                            PRODUCE LESS THAN HALF OF MACHINE ACCURACY
!                   IFAIL=4, CABS(Z) TOO LARGE  - NO COMPUTATION
!                            COMPLETE LOSS OF ACCURACY BY ARGUMENT
!                            REDUCTION
!                   IFAIL=5, ERROR              - NO COMPUTATION,
!                            ALGORITHM TERMINATION CONDITION NOT MET
!
!     LONG DESCRIPTION
!     ================
!
!         AI AND DAI ARE COMPUTED FOR CABS(Z)>1.0 FROM THE K BESSEL
!         functionS BY
!
!            AI(Z)=C*SQRT(Z)*K(1/3,ZTA) , DAI(Z)=-C*Z*K(2/3,ZTA)
!                           C=1.0/(PI*SQRT(3.0))
!                           ZTA=(2/3)*Z**(3/2)
!
!         WITH THE POWER SERIES FOR CABS(Z) <= 1.0.
!
!         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
!         MENTARY functionS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
!         OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
!         THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
!         THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
!         FLAG IFAIL=3 IS TRIGGERED WHERE UR=X02AJE()=UNIT ROUNDOFF.
!         ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
!         ALL SIGNIFICANCE IS LOST AND IFAIL=4. IN ORDER TO USE THE INT
!         function, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
!         LARGEST INTEGER, U3=X02BBE(). THUS, THE MAGNITUDE OF ZETA
!         MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
!         AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
!         PRECISION ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE
!         PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
!         ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
!         NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
!         DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
!         EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
!         NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
!         PRECISION ARITHMETIC. SIMILAR CONSIDERATIONS HOLD FOR OTHER
!         MACHINES.
!
!         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
!         BESSEL function CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
!         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
!         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
!         ELEMENTARY functionS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
!         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
!         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
!         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
!         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
!         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
!         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
!         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
!         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
!         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
!         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
!         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
!         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
!         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
!         OR -PI/2+P.
!
!     REFERENCES
!     ==========
!               HANDBOOK OF MATHEMATICAL functionS BY M. ABRAMOWITZ
!                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
!                 COMMERCE, 1955.
!
!               COMPUTATION OF BESSEL functionS OF COMPLEX ARGUMENT
!                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
!
!               A subroutine PACKAGE FOR BESSEL functionS OF A COMPLEX
!                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
!                 1018, MAY, 1985
!
!               A PORTABLE PACKAGE FOR BESSEL functionS OF A COMPLEX
!                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
!                 MATH. SOFTWARE, 1986
!
!     DATE WRITTEN   830501   (YYMMDD)
!     REVISION DATE  830501   (YYMMDD)
!     AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!
!     .. Parameters ..
  CHARACTER*6       SRNAME
  PARAMETER         (SRNAME='S17DGE')
!     .. Scalar Arguments ..
  COMPLEX           AI, Z
  INTEGER           IFAIL, NZ
  CHARACTER         DERIV, SCALE
!     .. Local Scalars ..
  COMPLEX           CONE, CSQ, S1, S2, TRM1, TRM2, Z3, ZTA
  REAL              AA, AD, AK, ALAZ, ALIM, ATRM, AZ, AZ3, BB, BK, &
                    C1, C2, CK, COEF, D1, D2, DIG, DK, ELIM, FID, &
                    FNU, R1M5, RL, SAVAA, SFAC, TOL, TTH, Z3I, Z3R, &
                    ZI, ZR
  INTEGER           ID, IERR, IFL, IFLAG, K, K1, K2, KODE, MR, NN, &
                    NREC
!     .. Local Arrays ..
  COMPLEX           CY(1)
  CHARACTER*80      REC(1)
!     .. External functions ..
  COMPLEX           S01EAE
  REAL              X02AHE, X02AJE, X02AME
  INTEGER           P01ABE, X02BBE, X02BHE, X02BJE, X02BKE, X02BLE
  EXTERNAL          S01EAE, X02AHE, X02AJE, X02AME, P01ABE, X02BBE, &
                    X02BHE, X02BJE, X02BKE, X02BLE
!     .. External subroutines ..
  EXTERNAL          DGXS17, DGZS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, LOG, LOG10, MAX, MIN, REAL, &
                    SQRT
!     .. Data statements ..
  DATA              TTH, C1, C2, COEF/6.66666666666666667E-01, &
                    3.55028053887817240E-01, &
                    2.58819403792806799E-01, &
                    1.83776298473930683E-01/
  DATA              CONE/(1.0E0,0.0E0)/
!     .. Executable Statements ..
  IERR = 0
  NREC = 0
  NZ = 0
  if (DERIV == 'F' .or. DERIV == 'f') then
   ID = 0
  else if (DERIV == 'D' .or. DERIV == 'd') then
   ID = 1
  ELSE
   ID = -1
  endif
  if (SCALE == 'U' .or. SCALE == 'u') then
   KODE = 1
  else if (SCALE == 'S' .or. SCALE == 's') then
   KODE = 2
  ELSE
   KODE = -1
  endif
  if (ID == -1) then
   IERR = 1
   NREC = 1
   WRITE (REC,FMT=99999) DERIV
  else if (KODE == -1) then
   IERR = 1
   NREC = 1
   WRITE (REC,FMT=99998) SCALE
  endif
  if (IERR == 0) then
   AZ = ABS(Z)
   TOL = MAX(X02AJE(),1.0E-18)
   FID = ID
   if (AZ > 1.0E0) then
!           ------------------------------------------------------------
!           CASE FOR CABS(Z)>1.0
!           ------------------------------------------------------------
      FNU = (1.0E0+FID)/3.0E0
!           ------------------------------------------------------------
!           SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!           TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!           ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW
!           LIMIT.
!           EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!           EXP(ELIM)>EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS
!           NEAR UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC
!           IS DONE.
!           RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR
!           LARGE Z.
!           DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!           ------------------------------------------------------------
      K1 = X02BKE()
      K2 = X02BLE()
      R1M5 = LOG10(REAL(X02BHE()))
      K = MIN(ABS(K1),ABS(K2))
      ELIM = 2.303E0*(K*R1M5-3.0E0)
      K1 = X02BJE() - 1
      AA = R1M5*K1
      DIG = MIN(AA,18.0E0)
      AA = AA*2.303E0
      ALIM = ELIM + MAX(-AA,-41.45E0)
      RL = 1.2E0*DIG + 3.0E0
      ALAZ = LOG(AZ)
!           ------------------------------------------------------------
!           TEST FOR RANGE
!           ------------------------------------------------------------
      AA = 0.5E0/TOL
      BB = X02BBE(1.0E0)*0.5E0
      AA = MIN(AA,BB,X02AHE(1.0E0))
      AA = AA**TTH
      if (AZ > AA) then
         NZ = 0
         IERR = 4
         NREC = 1
         WRITE (REC,FMT=99997) AZ, AA
      ELSE
         AA = SQRT(AA)
         SAVAA = AA
         if (AZ > AA) then
            IERR = 3
            NREC = 1
            WRITE (REC,FMT=99996) AZ, AA
         endif
         CSQ = SQRT(Z)
         ZTA = Z*CSQ*CMPLX(TTH,0.0E0)
!              ---------------------------------------------------------
!              RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS
!              SMALL
!              ---------------------------------------------------------
         IFLAG = 0
         SFAC = 1.0E0
         ZI = AIMAG(Z)
         ZR = REAL(Z)
         AK = AIMAG(ZTA)
         if (ZR < 0.0E0) then
            BK = REAL(ZTA)
            CK = -ABS(BK)
            ZTA = CMPLX(CK,AK)
         endif
         if (ZI == 0.0E0) then
            if (ZR <= 0.0E0) ZTA = CMPLX(0.0E0,AK)
         endif
         AA = REAL(ZTA)
         if (AA >= 0.0E0 .and. ZR > 0.0E0) then
            if (KODE /= 2) then
!                    ---------------------------------------------------
!                    UNDERFLOW TEST
!                    ---------------------------------------------------
               if (AA >= ALIM) then
                  AA = -AA - 0.25E0*ALAZ
                  IFLAG = 2
                  SFAC = 1.0E0/TOL
                  if (AA < (-ELIM)) then
                     NZ = 1
                     AI = CMPLX(0.0E0,0.0E0)
                     IFAIL = P01ABE(IFAIL,IERR,SRNAME,NREC,REC)
                     return
                  endif
               endif
            endif
            CALL DGXS17(ZTA,FNU,KODE,1,CY,NZ,TOL,ELIM,ALIM)
         ELSE
            if (KODE /= 2) then
!                    ---------------------------------------------------
!                    OVERFLOW TEST
!                    ---------------------------------------------------
               if (AA <= (-ALIM)) then
                  AA = -AA + 0.25E0*ALAZ
                  IFLAG = 1
                  SFAC = TOL
                  if (AA > ELIM) goto 20
               endif
            endif
!                 ------------------------------------------------------
!                 DGXS17 AND DGZS17 return EXP(ZTA)*K(FNU,ZTA) ON KODE=2
!                 ------------------------------------------------------
            MR = 1
            if (ZI < 0.0E0) MR = -1
            CALL DGZS17(ZTA,FNU,KODE,MR,1,CY,NN,RL,TOL,ELIM,ALIM)
            if (NN >= 0) then
               NZ = NZ + NN
               goto 40
            else if (NN == (-3)) then
               NZ = 0
               IERR = 4
               NREC = 1
               WRITE (REC,FMT=99997) AZ, SAVAA
               IFAIL = P01ABE(IFAIL,IERR,SRNAME,NREC,REC)
               return
            else if (NN /= (-1)) then
               NZ = 0
               IERR = 5
               NREC = 1
               WRITE (REC,FMT=99995)
               IFAIL = P01ABE(IFAIL,IERR,SRNAME,NREC,REC)
               return
            endif
   20             NZ = 0
            IERR = 2
            NREC = 1
            WRITE (REC,FMT=99994)
            IFAIL = P01ABE(IFAIL,IERR,SRNAME,NREC,REC)
            return
         endif
   40          S1 = CY(1)*CMPLX(COEF,0.0E0)
         if (IFLAG /= 0) then
            S1 = S1*CMPLX(SFAC,0.0E0)
            if (ID == 1) then
               S1 = -S1*Z
               AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
            ELSE
               S1 = S1*CSQ
               AI = S1*CMPLX(1.0E0/SFAC,0.0E0)
            endif
         else if (ID == 1) then
            AI = -Z*S1
         ELSE
            AI = CSQ*S1
         endif
      endif
   ELSE
!           ------------------------------------------------------------
!           POWER SERIES FOR CABS(Z) <= 1.
!           ------------------------------------------------------------
      S1 = CONE
      S2 = CONE
      if (AZ < TOL) then
         AA = 1.0E+3*X02AME()
         S1 = CMPLX(0.0E0,0.0E0)
         if (ID == 1) then
            AI = -CMPLX(C2,0.0E0)
            AA = SQRT(AA)
            if (AZ > AA) S1 = Z*Z*CMPLX(0.5E0,0.0E0)
            AI = AI + S1*CMPLX(C1,0.0E0)
         ELSE
            if (AZ > AA) S1 = CMPLX(C2,0.0E0)*Z
            AI = CMPLX(C1,0.0E0) - S1
         endif
      ELSE
         AA = AZ*AZ
         if (AA >= TOL/AZ) then
            TRM1 = CONE
            TRM2 = CONE
            ATRM = 1.0E0
            Z3 = Z*Z*Z
            AZ3 = AZ*AA
            AK = 2.0E0 + FID
            BK = 3.0E0 - FID - FID
            CK = 4.0E0 - FID
            DK = 3.0E0 + FID + FID
            D1 = AK*DK
            D2 = BK*CK
            AD = MIN(D1,D2)
            AK = 24.0E0 + 9.0E0*FID
            BK = 30.0E0 - 9.0E0*FID
            Z3R = REAL(Z3)
            Z3I = AIMAG(Z3)
            DO 60 K = 1, 25
               TRM1 = TRM1*CMPLX(Z3R/D1,Z3I/D1)
               S1 = S1 + TRM1
               TRM2 = TRM2*CMPLX(Z3R/D2,Z3I/D2)
               S2 = S2 + TRM2
               ATRM = ATRM*AZ3/AD
               D1 = D1 + AK
               D2 = D2 + BK
               AD = MIN(D1,D2)
               if (ATRM < TOL*AD) then
                  goto 80
               ELSE
                  AK = AK + 18.0E0
                  BK = BK + 18.0E0
               endif
   60             continue
         endif
   80          if (ID == 1) then
            AI = -S2*CMPLX(C2,0.0E0)
            if (AZ > TOL) AI = AI + Z*Z*S1*CMPLX(C1/(1.0E0+FID), &
                                  0.0E0)
            if (KODE /= 1) then
               ZTA = Z*SQRT(Z)*CMPLX(TTH,0.0E0)
!                     AI = AI*EXP(ZTA)
               IFL = 1
               AI = AI*S01EAE(ZTA,IFL)
            endif
         ELSE
            AI = S1*CMPLX(C1,0.0E0) - Z*S2*CMPLX(C2,0.0E0)
            if (KODE /= 1) then
               ZTA = Z*SQRT(Z)*CMPLX(TTH,0.0E0)
!                     AI = AI*EXP(ZTA)
               IFL = 1
               AI = AI*S01EAE(ZTA,IFL)
            endif
         endif
      endif
   endif
  endif
  IFAIL = P01ABE(IFAIL,IERR,SRNAME,NREC,REC)
  return
!
  99999 FORMAT (1X,'** On entry, DERIV has illegal value: DERIV = ''',A, &
         '''')
  99998 FORMAT (1X,'** On entry, SCALE has illegal value: SCALE = ''',A, &
         '''')
  99997 FORMAT (1X,'** No computation because abs(Z) =',1P,E13.5,'>', &
         E13.5)
  99996 FORMAT (1X,'** Results lack precision because abs(Z) =',1P,E13.5, &
         '>',E13.5)
  99995 FORMAT (1X,'** No computation - algorithm termination condition ', &
         'not met.')
  99994 FORMAT (1X,'** No computation because real(ZTA) too large, where', &
         ' ZTA = (2/3)*Z**(3/2).')
  END
  subroutine S17DLE(M,FNU,Z,N,SCALE,CY,NZ,IFAIL)
!     MARK 13 RELEASE. NAG COPYRIGHT 1988.
!     MARK 14 REVISED. IER-781 (DEC 1989).
!
!     Original name: CBESH
!
!     PURPOSE  TO COMPUTE THE H-BESSEL functionS OF A COMPLEX ARGUMENT
!
!     DESCRIPTION
!     ===========
!
!         ON SCALE='U', S17DLE COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
!         HANKEL (BESSEL) functionS CY(J)=H(M,FNU+J-1,Z) FOR KINDS M=1
!         OR 2, REAL, NONNEGATIVE ORDERS FNU+J-1, J=1,...,N, AND COMPLEX
!         Z /= CMPLX(0.0E0,0.0E0) IN THE CUT PLANE -PI < ARG(Z) <= PI.
!         ON SCALE='S', S17DLE COMPUTES THE SCALED HANKEL functionS
!
!         CY(I)=H(M,FNU+J-1,Z)*EXP(-MM*Z*I)       MM=3-2M,      I**2=-1.
!
!         WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE UPPER
!         AND LOWER HALF PLANES. DEFINITIONS AND NOTATION ARE FOUND IN
!         THE NBS HANDBOOK OF MATHEMATICAL functionS (REF. 1).
!
!         INPUT
!           Z      - Z=CMPLX(X,Y), Z /= CMPLX(0.,0.),-PI < ARG(Z) <= PI
!           FNU    - ORDER OF INITIAL H function, FNU >= 0.0E0
!           SCALE  - A PARAMETER TO INDICATE THE SCALING OPTION
!                    SCALE = 'U' OR SCALE = 'u' returnS
!                             CY(J)=H(M,FNU+J-1,Z),      J=1,...,N
!                          = 'S' OR SCALE = 's' returnS
!                             CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))
!                                  J=1,...,N  ,  I**2=-1
!           M      - KIND OF HANKEL function, M=1 OR 2
!           N      - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1
!
!         OUTPUT
!           CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
!                    VALUES FOR THE SEQUENCE
!                    CY(J)=H(M,FNU+J-1,Z)  OR
!                    CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))  J=1,...,N
!                    DEPENDING ON SCALE, I**2=-1.
!           NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
!                    NZ= 0   , NORMAL return
!                    NZ>0 , FIRST NZ COMPONENTS OF CY SET TO ZERO
!                              DUE TO UNDERFLOW, CY(J)=CMPLX(0.0,0.0)
!                              J=1,...,NZ WHEN Y>0.0 AND M=1 OR
!                              Y < 0.0 AND M=2. FOR THE COMPLMENTARY
!                              HALF PLANES, NZ STATES ONLY THE NUMBER
!                              OF UNDERFLOWS.
!           IERR    -ERROR FLAG
!                    IERR=0, NORMAL return - COMPUTATION COMPLETED
!                    IERR=1, INPUT ERROR   - NO COMPUTATION
!                    IERR=2, OVERFLOW      - NO COMPUTATION,
!                            CABS(Z) TOO SMALL
!                    IERR=3  OVERFLOW      - NO COMPUTATION,
!                            FNU+N-1 TOO LARGE
!                    IERR=4, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
!                            BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
!                            REDUCTION PRODUCE LESS THAN HALF OF MACHINE
!                            ACCURACY
!                    IERR=5, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
!                            TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
!                            CANCE BY ARGUMENT REDUCTION
!                    IERR=6, ERROR              - NO COMPUTATION,
!                            ALGORITHM TERMINATION CONDITION NOT MET
!
!     LONG DESCRIPTION
!     ================
!
!         THE COMPUTATION IS CARRIED OUT BY THE RELATION
!
!         H(M,FNU,Z)=(1/MP)*EXP(-MP*FNU)*K(FNU,Z*EXP(-MP))
!             MP=MM*HPI*I,  MM=3-2*M,  HPI=PI/2,  I**2=-1
!
!         FOR M=1 OR 2 WHERE THE K BESSEL function IS COMPUTED FOR THE
!         RIGHT HALF PLANE RE(Z) >= 0.0. THE K function IS continueD
!         TO THE LEFT HALF PLANE BY THE RELATION
!
!         K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
!         MP=MR*PI*I, MR=+1 OR -1, RE(Z)>0, I**2=-1
!
!         WHERE I(FNU,Z) IS THE I BESSEL function.
!
!         EXPONENTIAL DECAY OF H(M,FNU,Z) OCCURS IN THE UPPER HALF Z
!         PLANE FOR M=1 AND THE LOWER HALF Z PLANE FOR M=2.  EXPONENTIAL
!         GROWTH OCCURS IN THE COMPLEMENTARY HALF PLANES.  SCALING
!         BY EXP(-MM*Z*I) REMOVES THE EXPONENTIAL BEHAVIOR IN THE
!         WHOLE Z PLANE FOR Z TO INFINITY.
!
!         FOR NEGATIVE ORDERS,THE FORMULAE
!
!               H(1,-FNU,Z) = H(1,FNU,Z)*CEXP( PI*FNU*I)
!               H(2,-FNU,Z) = H(2,FNU,Z)*CEXP(-PI*FNU*I)
!                         I**2=-1
!
!         CAN BE USED.
!
!         IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
!         MENTARY functionS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
!         LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
!         CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
!         LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
!         IERR=4 IS TRIGGERED WHERE UR=X02AJE()=UNIT ROUNDOFF. ALSO
!         IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
!         LOST AND IERR=5. IN ORDER TO USE THE INT function, ARGUMENTS
!         MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
!         INTEGER, U3=X02BBE(). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
!         RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
!         ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
!         ARITHMETIC AND 1.3E+8, 1.8E+16, 2.1E+9 IN DOUBLE PRECISION
!         ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
!         THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
!         TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
!         IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
!         SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.
!
!         THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
!         BESSEL function CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
!         ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
!         SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
!         ELEMENTARY functionS. HERE, S=MAX(1,ABS(LOG10(CABS(Z))),
!         ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
!         CABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
!         HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
!         ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
!         SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
!         THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
!         0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
!         THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
!         COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
!         BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
!         COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
!         MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
!         THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
!         OR -PI/2+P.
!
!     REFERENCES
!     ==========
!               HANDBOOK OF MATHEMATICAL functionS BY M. ABRAMOWITZ
!                 AND I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF
!                 COMMERCE, 1955.
!
!               COMPUTATION OF BESSEL functionS OF COMPLEX ARGUMENT
!                 BY D. E. AMOS, SAND83-0083, MAY, 1983.
!
!               COMPUTATION OF BESSEL functionS OF COMPLEX ARGUMENT
!                 AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY, 1983
!
!               A subroutine PACKAGE FOR BESSEL functionS OF A COMPLEX
!                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-
!                 1018, MAY, 1985
!
!               A PORTABLE PACKAGE FOR BESSEL functionS OF A COMPLEX
!                 ARGUMENT AND NONNEGATIVE ORDER BY D. E. AMOS, TRANS.
!                 MATH. SOFTWARE, 1986
!
!     DATE WRITTEN   830501   (YYMMDD)
!     REVISION DATE  830501   (YYMMDD)
!     AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!
!     .. Parameters ..
  CHARACTER*6       SRNAME
  PARAMETER         (SRNAME='S17DLE')
!     .. Scalar Arguments ..
  COMPLEX           Z
  REAL              FNU
  INTEGER           IFAIL, M, N, NZ
  CHARACTER*1       SCALE
!     .. Array Arguments ..
  COMPLEX           CY(N)
!     .. Local Scalars ..
  COMPLEX           CSGN, ZN, ZT
  REAL              AA, ALIM, ALN, ARG, ASCLE, ATOL, AZ, BB, CPN, &
                    DIG, ELIM, FMM, FN, FNUL, HPI, R1M5, RHPI, RL, &
                    RTOL, SGN, SPN, TOL, UFL, XN, XX, YN, YY
  INTEGER           I, IERR, INU, INUH, IR, K, K1, K2, KODE, MM, MR, &
                    NN, NREC, NUF, NW
!     .. Local Arrays ..
  CHARACTER*80      REC(1)
!     .. External functions ..
  REAL              X02AHE, X02AJE
  INTEGER           P01ABE, X02BBE, X02BHE, X02BJE, X02BKE, X02BLE
  EXTERNAL          X02AHE, X02AJE, P01ABE, X02BBE, X02BHE, X02BJE, &
                    X02BKE, X02BLE
!     .. External subroutines ..
  EXTERNAL          DEVS17, DGXS17, DLYS17, DLZS17
!     .. Intrinsic functions ..
  INTRINSIC         ABS, AIMAG, CMPLX, COS, EXP, INT, LOG, LOG10, &
                    MAX, MIN, MOD, REAL, SIGN, SIN, SQRT
!     .. Data statements ..
!
  DATA              HPI/1.57079632679489662E0/
!     .. Executable Statements ..
  NZ = 0
  NREC = 0
  XX = REAL(Z)
  YY = AIMAG(Z)
  IERR = 0
  if (SCALE == 'U' .or. SCALE == 'u') then
   KODE = 1
  else if (SCALE == 'S' .or. SCALE == 's') then
   KODE = 2
  ELSE
   KODE = -1
  endif
  if (XX == 0.0E0 .and. YY == 0.0E0) then
   IERR = 1
   NREC = 1
   WRITE (REC,FMT=99999)
  else if (FNU < 0.0E0) then
   IERR = 1
   NREC = 1
   WRITE (REC,FMT=99998) FNU
  else if (KODE == -1) then
   IERR = 1
   NREC = 1
   WRITE (REC,FMT=99997) SCALE
  else if (N < 1) then
   IERR = 1
   NREC = 1
   WRITE (REC,FMT=99996) N
  else if (M < 1 .or. M > 2) then
   IERR = 1
   NREC = 1
   WRITE (REC,FMT=99995) M
  endif
  if (IERR == 0) then
   NN = N
!        ---------------------------------------------------------------
!        SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!        TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!        ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!        EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!        EXP(ELIM)>EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!        UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!        RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR
!        LARGE Z.
!        DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!        FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE
!        FNU
!        ---------------------------------------------------------------
   TOL = MAX(X02AJE(),1.0E-18)
   K1 = X02BKE()
   K2 = X02BLE()
   R1M5 = LOG10(REAL(X02BHE()))
   K = MIN(ABS(K1),ABS(K2))
   ELIM = 2.303E0*(K*R1M5-3.0E0)
   K1 = X02BJE() - 1
   AA = R1M5*K1
   DIG = MIN(AA,18.0E0)
   AA = AA*2.303E0
   ALIM = ELIM + MAX(-AA,-41.45E0)
   FNUL = 10.0E0 + 6.0E0*(DIG-3.0E0)
   RL = 1.2E0*DIG + 3.0E0
   FN = FNU + NN - 1
   MM = 3 - M - M
   FMM = MM
   ZN = Z*CMPLX(0.0E0,-FMM)
   XN = REAL(ZN)
   YN = AIMAG(ZN)
   AZ = ABS(Z)
!        ---------------------------------------------------------------
!        TEST FOR RANGE
!        ---------------------------------------------------------------
   AA = 0.5E0/TOL
   BB = X02BBE(1.0E0)*0.5E0
   AA = MIN(AA,BB,X02AHE(1.0E0))
   if (AZ <= AA) then
      if (FN <= AA) then
         AA = SQRT(AA)
         if (AZ > AA) then
            IERR = 4
            NREC = 1
            WRITE (REC,FMT=99994) AZ, AA
         else if (FN > AA) then
            IERR = 4
            NREC = 1
            WRITE (REC,FMT=99993) FN, AA
         endif
!              ---------------------------------------------------------
!              OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
!              ---------------------------------------------------------
         UFL = EXP(-ELIM)
         if (AZ >= UFL) then
            if (FNU > FNUL) then
!                    ---------------------------------------------------
!                    UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU>FNUL
!                    ---------------------------------------------------
               MR = 0
               if ((XN < 0.0E0) .or. (XN == 0.0E0 .and. YN <&
                     0.0E0 .and. M == 2)) then
                  MR = -MM
                  if (XN == 0.0E0 .and. YN < 0.0E0) ZN = -ZN
               endif
               CALL DLYS17(ZN,FNU,KODE,MR,NN,CY,NW,TOL,ELIM,ALIM)
               if (NW < 0) then
                  goto 40
               ELSE
                  NZ = NZ + NW
               endif
            ELSE
               if (FN > 1.0E0) then
                  if (FN > 2.0E0) then
                     CALL DEVS17(ZN,FNU,KODE,2,NN,CY,NUF,TOL,ELIM, &
                                   ALIM)
                     if (NUF < 0) then
                        goto 60
                     ELSE
                        NZ = NZ + NUF
                        NN = NN - NUF
!                             ------------------------------------------
!                             HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1
!                             ON return FROM DEVS17
!                             IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
!                             ------------------------------------------
                        if (NN == 0) then
                           if (XN < 0.0E0) then
                              goto 60
                           ELSE
                              IFAIL = P01ABE(IFAIL,IERR,SRNAME, &
                                        NREC,REC)
                              return
                           endif
                        endif
                     endif
                  else if (AZ <= TOL) then
                     ARG = 0.5E0*AZ
                     ALN = -FN*LOG(ARG)
                     if (ALN > ELIM) goto 60
                  endif
               endif
               if ((XN < 0.0E0) .or. (XN == 0.0E0 .and. YN <&
                     0.0E0 .and. M == 2)) then
!                       ------------------------------------------------
!                       LEFT HALF PLANE COMPUTATION
!                       ------------------------------------------------
                  MR = -MM
                  CALL DLZS17(ZN,FNU,KODE,MR,NN,CY,NW,RL,FNUL,TOL, &
                                ELIM,ALIM)
                  if (NW < 0) then
                     goto 40
                  ELSE
                     NZ = NW
                  endif
               ELSE
!                       ------------------------------------------------
!                       RIGHT HALF PLANE COMPUTATION, XN >= 0. .and.
!                       (XN /= 0. .or. YN >= 0. .or. M=1)
!                       ------------------------------------------------
                  CALL DGXS17(ZN,FNU,KODE,NN,CY,NZ,TOL,ELIM,ALIM)
               endif
            endif
!                 ------------------------------------------------------
!                 H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)
!
!                 ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
!                 ------------------------------------------------------
            SGN = SIGN(HPI,-FMM)
!                 ------------------------------------------------------
!                 CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF
!                 SIGNIFICANCE WHEN FNU IS LARGE
!                 ------------------------------------------------------
            INU = INT(FNU)
            INUH = INU/2
            IR = INU - 2*INUH
            ARG = (FNU-INU+IR)*SGN
            RHPI = 1.0E0/SGN
            CPN = RHPI*COS(ARG)
            SPN = RHPI*SIN(ARG)
!                 ZN = CMPLX(-SPN,CPN)
            CSGN = CMPLX(-SPN,CPN)
!                 if (MOD(INUH,2)==1) ZN = -ZN
            if (MOD(INUH,2) == 1) CSGN = -CSGN
            ZT = CMPLX(0.0E0,-FMM)
            RTOL = 1.0E0/TOL
            ASCLE = UFL*RTOL
            DO 20 I = 1, NN
!                    CY(I) = CY(I)*ZN
!                    ZN = ZN*ZT
               ZN = CY(I)
               AA = REAL(ZN)
               BB = AIMAG(ZN)
               ATOL = 1.0E0
               if (MAX(ABS(AA),ABS(BB)) <= ASCLE) then
                  ZN = ZN*RTOL
                  ATOL = TOL
               endif
               ZN = ZN*CSGN
               CY(I) = ZN*ATOL
               CSGN = CSGN*ZT
   20             continue
            IFAIL = P01ABE(IFAIL,IERR,SRNAME,NREC,REC)
            return
   40             if (NW == (-3)) then
               NZ = 0
               IERR = 5
               NREC = 1
               WRITE (REC,FMT=99988) AZ, AA
               IFAIL = P01ABE(IFAIL,IERR,SRNAME,NREC,REC)
               return
            else if (NW /= (-1)) then
               NZ = 0
               IERR = 6
               NREC = 1
               WRITE (REC,FMT=99992)
               IFAIL = P01ABE(IFAIL,IERR,SRNAME,NREC,REC)
               return
            endif
   60             IERR = 3
            NZ = 0
            NREC = 1
            WRITE (REC,FMT=99991) FN
            IFAIL = P01ABE(IFAIL,IERR,SRNAME,NREC,REC)
            return
         ELSE
            IERR = 2
            NZ = 0
            NREC = 1
            WRITE (REC,FMT=99990) AZ, UFL
            IFAIL = P01ABE(IFAIL,IERR,SRNAME,NREC,REC)
            return
         endif
      ELSE
         NZ = 0
         IERR = 5
         NREC = 1
         WRITE (REC,FMT=99989) FN, AA
      endif
   ELSE
      NZ = 0
      IERR = 5
      NREC = 1
      WRITE (REC,FMT=99988) AZ, AA
   endif
  endif
  IFAIL = P01ABE(IFAIL,IERR,SRNAME,NREC,REC)
  return
!
  99999 FORMAT (1X,'** On entry, Z = (0.0,0.0)')
  99998 FORMAT (1X,'** On entry, FNU < 0: FNU = ',E13.5)
  99997 FORMAT (1X,'** On entry, SCALE has an illegal value: SCALE = ''', &
         A,'''')
  99996 FORMAT (1X,'** On entry, N <= 0: N = ',I16)
  99995 FORMAT (1X,'** On entry, M has illegal value: M = ',I16)
  99994 FORMAT (1X,'** Results lack precision because abs(Z) =',1P,E13.5, &
         '>',E13.5)
  99993 FORMAT (1X,'** Results lack precision, FNU+N-1 =',1P,E13.5, &
         '>',E13.5)
  99992 FORMAT (1X,'** No computation - algorithm termination condition ', &
         'not met.')
  99991 FORMAT (1X,'** No computation because FNU+N-1 =',1P,E13.5,' is t', &
         'oo large.')
  99990 FORMAT (1X,'** No computation because abs(Z) =',1P,E13.5,'<', &
         E13.5)
  99989 FORMAT (1X,'** No computation because FNU+N-1 =',1P,E13.5,'>', &
         E13.5)
  99988 FORMAT (1X,'** No computation because abs(Z) =',1P,E13.5,'>', &
         E13.5)
  END
  REAL function X02AHE(X)
!     MARK 9 RELEASE. NAG COPYRIGHT 1981.
!     MARK 11.5(F77) REVISED. (SEPT 1985.)
!
!     * MAXIMUM ARGUMENT FOR SIN AND COS *
!     returnS THE LARGEST POSITIVE REAL NUMBER MAXSC SUCH THAT
!     SIN(MAXSC) AND COS(MAXSC) CAN BE SUCCESSFULLY COMPUTED
!     BY THE COMPILER SUPPLIED SIN AND COS ROUTINES.
!
!     .. Scalar Arguments ..
  REAL                 X
  REAL CONX02
  DATA CONX02 /1.677721600000E+7 /
!     .. Executable Statements ..
  X02AHE = CONX02
  return
  END
  REAL function X02AJE()
!     MARK 12 RELEASE. NAG COPYRIGHT 1986.
!
!     returnS  (1/2)*B**(1-P)  IF ROUNDS IS .true.
!     returnS  B**(1-P)  OTHERWISE
!
  REAL CONX02
  DATA CONX02 /1.4210854715202E-14 /
!bc      DATA CONX02 /1.421090000020E-14 /
!     .. Executable Statements ..
  X02AJE = CONX02
  return
  END
  REAL function X02ALE()
!     MARK 12 RELEASE. NAG COPYRIGHT 1986.
!
!     returnS  (1 - B**(-P)) * B**EMAX  (THE LARGEST POSITIVE MODEL
!     NUMBER)
!
  REAL CONX02
! DK DK DK      DATA CONX02 /0577757777777777777777B /
  DATA CONX02 /1.e30/
!     .. Executable Statements ..
  X02ALE = CONX02
  return
  END
  REAL function X02AME()
!     MARK 12 RELEASE. NAG COPYRIGHT 1986.
!
!     returnS THE 'SAFE RANGE' PARAMETER
!     I.E. THE SMALLEST POSITIVE MODEL NUMBER Z SUCH THAT
!     FOR ANY X WHICH SATISFIES X >= Z AND X <= 1/Z
!     THE FOLLOWING CAN BE COMPUTED WITHOUT OVERFLOW, UNDERFLOW OR OTHER
!     ERROR
!
!        -X
!        1.0/X
!        SQRT(X)
!        LOG(X)
!        EXP(LOG(X))
!        Y**(LOG(X)/LOG(Y)) FOR ANY Y
!
  REAL CONX02
! DK DK DK     DATA CONX02 /0200044000000000000004B /
  DATA CONX02 /1.e-27/
!     .. Executable Statements ..
  X02AME = CONX02
  return
  END
  REAL function X02ANE()
!     MARK 15 RELEASE. NAG COPYRIGHT 1991.
!
!     returns the 'safe range' parameter for complex numbers,
!     i.e. the smallest positive model number Z such that
!     for any X which satisfies X >= Z and X <= 1/Z
!     the following can be computed without overflow, underflow or other
!     error
!
!        -W
!        1.0/W
!        SQRT(W)
!        LOG(W)
!        EXP(LOG(W))
!        Y**(LOG(W)/LOG(Y)) for any Y
!        ABS(W)
!
!     where W is any of cmplx(X,0), cmplx(0,X), cmplx(X,X),
!                   cmplx(1/X,0), cmplx(0,1/X), cmplx(1/X,1/X).
!
  REAL CONX02
!bc      DATA CONX02 /0000006220426276611547B /
!! DK DK      DATA CONX02 / 2.708212596942E-1233 /
  DATA CONX02 / 2.708212596942E-30 /
!     .. Executable Statements ..
  X02ANE = CONX02
  return
  END
  INTEGER function X02BBE(X)
!     NAG COPYRIGHT 1975
!     MARK 4.5 RELEASE
!     MARK 11.5(F77) REVISED. (SEPT 1985.)
!     * MAXINT *
!     returnS THE LARGEST INTEGER REPRESENTABLE ON THE COMPUTER
!     THE X PARAMETER IS NOT USED
!     .. Scalar Arguments ..
  REAL                    X
!     .. Executable Statements ..
!     FOR ICL 1900
!     X02BBE = 8388607
! DK DK DK      X02BBE =       70368744177663
  X02BBE =       744177663
  return
  END
  INTEGER function X02BHE()
!     MARK 12 RELEASE. NAG COPYRIGHT 1986.
!
!     returnS THE MODEL PARAMETER, B.
!
!     .. Executable Statements ..
  X02BHE =     2
  return
  END
  INTEGER function X02BJE()
!     MARK 12 RELEASE. NAG COPYRIGHT 1986.
!
!     returnS THE MODEL PARAMETER, p.
!
!     .. Executable Statements ..
  X02BJE =    47
  return
  END
  INTEGER function X02BKE()
!     MARK 12 RELEASE. NAG COPYRIGHT 1986.
!
!     returnS THE MODEL PARAMETER, EMIN.
!
!     .. Executable Statements ..
  X02BKE =  -8192
  return
  END
  INTEGER function X02BLE()
!     MARK 12 RELEASE. NAG COPYRIGHT 1986.
!
!     returnS THE MODEL PARAMETER, EMAX.
!
!     .. Executable Statements ..
  X02BLE =  8189
  return
  END
  subroutine X04AAE(I,NERR)
!     MARK 7 RELEASE. NAG COPYRIGHT 1978
!     MARK 7C REVISED IER-190 (MAY 1979)
!     MARK 11.5(F77) REVISED. (SEPT 1985.)
!     MARK 14 REVISED. IER-829 (DEC 1989).
!     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
!     (STORED IN NERR1).
!     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
!     VALUE SPECIFIED BY NERR.
!
!     .. Scalar Arguments ..
  INTEGER           I, NERR
!     .. Local Scalars ..
  INTEGER           NERR1
!     .. Save statement ..
  SAVE              NERR1
!     .. Data statements ..
  DATA              NERR1/0/
!     .. Executable Statements ..
  if (I == 0) NERR = NERR1
  if (I == 1) NERR1 = NERR
  return
  END
  subroutine X04BAE(NOUT,REC)
!     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
!
!     X04BAE writes the contents of REC to the unit defined by NOUT.
!
!     Trailing blanks are not output, except that if REC is entirely
!     blank, a single blank character is output.
!     If NOUT < 0, i.e. if NOUT is not a valid Fortran unit identifier,
!     then no output occurs.
!
!     .. Scalar Arguments ..
  INTEGER           NOUT
  CHARACTER*(*)     REC
!     .. Local Scalars ..
  INTEGER           I
!     .. Intrinsic functions ..
  INTRINSIC         LEN
!     .. Executable Statements ..
  if (NOUT >= 0) then
!        Remove trailing blanks
   DO 20 I = LEN(REC), 2, -1
      if (REC(I:I) /= ' ') goto 40
   20    continue
!        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
  endif
  return
!
  99999 FORMAT (A)
  END

