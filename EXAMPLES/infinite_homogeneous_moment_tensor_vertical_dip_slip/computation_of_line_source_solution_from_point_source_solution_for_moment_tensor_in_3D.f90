! The algorithm for computation of convolution
! y_k = sumi_p1^p2 (h_i*x_{k-i})
! The meaning of the variable
! Function:  stf[t] source time function
!
! Variable:  A:   the amplitude of the source time function
!            ts_stf: start time of the source time function
!            te_stf: end time of the source time function
!            ts_integrate: start time of the integration
!            te_integrate: end time of the integration
!            deltat: time step used in output
!
! we compute the seismic response for a line source from the integration of the seismic responce of the point source along the
! line.
! thus V_2D(x,z,t) = Int_{y_min}^{y_max} V_3D(x,y,z,t) dy

Implicit none

double precision :: amplitude, ts_stf, te_stf, ts_integrate, te_integrate
double precision, external :: stf, stf_integration, stf_convolution, dirac_function, stf_differential
integer :: i,j,k,p,q,ns_integration,ne_integration
double precision :: time, deltat,numerical_convolution_temp
double precision :: c_p, c_s, rho
double precision, dimension(3) :: gamma

double precision :: dominant_frequency
double precision, dimension(3) :: position_source
double precision, dimension(3) :: position_receiver
double precision :: distance_source_receiver
double precision :: ts_receiver,te_receiver
double precision, dimension(3,3) :: moment_tensor

double precision :: y_direction_step
double precision :: ys_direction_integration, ye_direction_integration
integer :: nys_direction_integration, nye_direction_integration

double precision :: distance_divided_by_c_p, distance_divided_by_c_s
integer :: ndistance_divided_by_c_p, ndistance_divided_by_c_s

double precision :: displacement, time_p_wave, time_s_wave
double precision :: pi = 3.141592653589793238462643d0
integer :: displacement_direction_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!media property
c_p = 3200.0d0
c_s = 1847.5d0
rho = 2200.0d0  ! in order to match the result of EX2DVAEL_Berg&If
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! moment_tensor
moment_tensor = 0.0d0
moment_tensor(1,3) = 1.0d0
moment_tensor(3,1) = 1.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Ricker wavelet source
Amplitude = -1.d8
dominant_frequency = 10.0d0
ts_stf = - 1.d0/dominant_frequency * 2.d0
te_stf = 1.d0/dominant_frequency * 5.0 + ts_stf
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!position_receiver
displacement_direction_index = 1
position_receiver(1) = 500.0
position_receiver(2) = 0.0
position_receiver(3) = 500.0
ts_receiver = ts_stf
te_receiver = 2.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Ricker wavelet source
deltat = 0.002d0
y_direction_step = 1.0d0/(4.d0*dominant_frequency) * c_s / 5.d0  ! we ensure that the smallest wavelength include four element
ys_direction_integration = - c_p * te_receiver * 1.5d0
ye_direction_integration = c_p * te_receiver * 1.5d0
nys_direction_integration = int(ys_direction_integration/y_direction_step) - 1
nye_direction_integration = int(ye_direction_integration/y_direction_step) + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!

ts_integrate = ts_receiver
te_integrate = te_receiver
ns_integration = int(ts_integrate/deltat) - 1
ne_integration = int(te_integrate/deltat) + 1

open(10, file='source_time_function.dat', status='unknown')
do i = ns_integration, ne_integration
   time = (dble(i)-1.d0)*deltat
   write(10,*) time, stf(ts_stf,te_stf,amplitude,time,dominant_frequency), &
                     stf_differential(ts_stf,te_stf,amplitude,time,dominant_frequency)
enddo
close(10)

open(20, file='x_displacement.dat', status='unknown')

! y_k = sumi_1^N (h_i*x_{k-i})

do k = ns_integration, ne_integration
  displacement = 0.d0
  do j = nys_direction_integration, nye_direction_integration
    position_source(1) = 0.d0
    position_source(2) = (dble(j)-1.d0) * y_direction_step
    position_source(3) = 0.d0

    distance_source_receiver = 0.d0
    do i = 1, 3
      distance_source_receiver = distance_source_receiver + (position_receiver(i) - position_source(i))**2
    enddo
    distance_source_receiver = dsqrt(distance_source_receiver)

    do i = 1,3
      gamma(i) = (position_receiver(i) - position_source(i)) / distance_source_receiver
    enddo

    distance_divided_by_c_p = distance_source_receiver / c_p
    distance_divided_by_c_s = distance_source_receiver / c_s

    ndistance_divided_by_c_p = int(distance_divided_by_c_p / deltat) + 1
    ndistance_divided_by_c_s = int(distance_divided_by_c_s / deltat) + 1

    numerical_convolution_temp = 0.d0
    do i = ndistance_divided_by_c_p, ndistance_divided_by_c_s
      time = (dble(k-i)-1.d0)*deltat
      numerical_convolution_temp = numerical_convolution_temp + &
                                 stf(ts_stf,te_stf,amplitude,time,dominant_frequency)*(dble(i)-1.d0)*deltat*deltat
    enddo
    time_p_wave = (dble(k)-1.d0)*deltat - distance_source_receiver / c_p
  time_s_wave = (dble(k)-1.d0)*deltat - distance_source_receiver / c_s

    do p = 1,3
    do q = 1,3
        displacement = displacement + &
                   !!!!!!!!!!
( 1.d0/(4.d0*rho*pi) * (15.0d0*gamma(displacement_direction_index)*gamma(p)*gamma(q) -  &
      3.0d0*gamma(displacement_direction_index)*dirac_function(p,q) - &
      3.0d0*gamma(p)*dirac_function(displacement_direction_index,q) - &
      3.0d0*gamma(q)*dirac_function(displacement_direction_index,p)) &
                         * 1.0d0/(distance_source_receiver**4) * numerical_convolution_temp * moment_tensor(p,q) + &
1.d0/(4.d0*rho*pi*c_p**2) * (6.0d0*gamma(displacement_direction_index)*gamma(p)*gamma(q) -  &
           1.0d0*gamma(displacement_direction_index)*dirac_function(p,q) - &
           1.0d0*gamma(p)*dirac_function(displacement_direction_index,q) - &
           1.0d0*gamma(q)*dirac_function(displacement_direction_index,p)) &
                          * 1.0d0/(distance_source_receiver**2) * moment_tensor(p,q) &
        * stf(ts_stf,te_stf,amplitude,time_p_wave,dominant_frequency) - &
1.d0/(4.d0*rho*pi*c_s**2) * (6.0d0*gamma(displacement_direction_index)*gamma(p)*gamma(q) -  &
           1.0d0*gamma(displacement_direction_index)*dirac_function(p,q) - &
           1.0d0*gamma(p)*dirac_function(displacement_direction_index,q) - &
           2.0d0*gamma(q)*dirac_function(displacement_direction_index,p)) &
                          * 1.0d0/(distance_source_receiver**2) * moment_tensor(p,q) &
        * stf(ts_stf,te_stf,amplitude,time_s_wave,dominant_frequency) + &
1.d0/(4.d0*rho*pi*c_p**3) * 1.0d0/(distance_source_receiver) * gamma(displacement_direction_index)*gamma(p)*gamma(q)  &
                    * moment_tensor(p,q) * stf_differential(ts_stf,te_stf,amplitude,time_p_wave,dominant_frequency) - &
1.d0/(4.d0*rho*pi*c_s**3) * (gamma(displacement_direction_index)*gamma(p) - dirac_function(displacement_direction_index,p)) &
                    * gamma(q) * 1.0d0/(distance_source_receiver) * moment_tensor(p,q) &
        * stf_differential(ts_stf,te_stf,amplitude,time_s_wave,dominant_frequency) ) * y_direction_step
    enddo
  enddo

  enddo

  write(20,*)(dble(k)-1.d0)*deltat + ts_stf, displacement
enddo

close(20)

end


double precision function dirac_function(i,j)
 integer :: i,j
 if(i /= j) then
   dirac_function = 0.0d0
 else
   dirac_function = 1.0d0
 endif
end function dirac_function

double precision function stf(ts_stf,te_stf,amplitude,t,dp)

 double precision :: ts_stf,te_stf,amplitude,t,dp
 double precision :: pi = 3.141592653589793238462643d0

! Ricker wavelet source
 if( t>= ts_stf .and. t <= te_stf ) then
   stf = amplitude * (1.d0 - 2.0*(pi*dp*(t + ts_stf))**2)*dexp(-(pi*dp*(t + ts_stf))**2)
 else
   stf = 0.0d0
 endif

end function stf

double precision function stf_differential(ts_stf,te_stf,amplitude,t,dp)

 double precision :: ts_stf,te_stf,amplitude,t,dp
 double precision :: pi = 3.141592653589793238462643d0

! Ricker wavelet source
 if( t>= ts_stf .and. t <= te_stf ) then
   stf_differential = amplitude * (-2.d0) * (pi*dp)**2 * (t + ts_stf) * &
                      (3.d0 - 2.0*(pi*dp*(t + ts_stf))**2)*dexp(-(pi*dp*(t + ts_stf))**2)
 else
   stf_differential = 0.0d0
 endif

end function stf_differential





