
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!                         Dimitri Komatitsch
!                     University of Pau, France
!
!                          (c) April 2007
!
!========================================================================

  subroutine attenuation_model(Qp_attenuation,Qs_attenuation,f0_attenuation, &
       inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,Mu_nu1,Mu_nu2)

! define the attenuation constants

  implicit none

  include "constants.h"

  double precision :: Qp_attenuation,Qs_attenuation,f0_attenuation
  double precision, dimension(N_SLS) :: inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2
  double precision :: Mu_nu1,Mu_nu2

  integer :: i_sls

  double precision, dimension(N_SLS) :: tau_epsilon_nu1,tau_sigma_nu1,tau_epsilon_nu2,tau_sigma_nu2

  double precision :: f1_attenuation, f2_attenuation
  
  
! f1 and f2 are computed as : f2/f1=12 and (log(f1)+log(f2))/2 = log(f0)
  f1_attenuation = exp(log(f0_attenuation)-log(12.d0)/2.d0)
  f2_attenuation = 12.d0 * f1_attenuation

! Call of C function that computes attenuation parameters (function in file "attenuation_compute_param.c"; 
! a main can be found in UTILS/attenuation directory).
! Beware of underscores in this function name; depending on your compiler and compilation options, you will have to add or 
! delete underscores. Also look in file "attenuation_compute_param.c" for this issue.
  call attenuation_compute_param(N_SLS, Qp_attenuation, Qs_attenuation, &
       f1_attenuation,f2_attenuation, &
       tau_sigma_nu1, tau_sigma_nu2, tau_epsilon_nu1, tau_epsilon_nu2)

! attenuation constants for standard linear solids

! nu1 is the dilatation mode
! nu2 is the shear mode

! array index (1) is the first standard linear solid, (2) is the second etc.

! from J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).
! Beware: these values implement specific values of the quality factors:
! Qp approximately equal to 13 and Qs approximately equal to 10,
! which means very high attenuation, see that paper for details.
! tau_epsilon_nu1(1) = 0.0334d0
! tau_sigma_nu1(1)   = 0.0303d0
! tau_epsilon_nu2(1) = 0.0352d0
! tau_sigma_nu2(1)   = 0.0287d0

! tau_epsilon_nu1(2) = 0.0028d0
! tau_sigma_nu1(2)   = 0.0025d0
! tau_epsilon_nu2(2) = 0.0029d0
! tau_sigma_nu2(2)   = 0.0024d0

! from J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
! in a linear viscoelastic medium, Geophysical Journal International,
! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).
! Beware: these values implement specific values of the quality factors:
! Qp approximately equal to 27 and Qs approximately equal to 20,
! which means very high attenuation, see that paper for details.
!  tau_epsilon_nu1(1) = 0.0325305d0
!  tau_sigma_nu1(1)   = 0.0311465d0
!  tau_epsilon_nu2(1) = 0.0332577d0
!  tau_sigma_nu2(1)   = 0.0304655d0

!  tau_epsilon_nu1(2) = 0.0032530d0
!  tau_sigma_nu1(2)   = 0.0031146d0
!  tau_epsilon_nu2(2) = 0.0033257d0
!  tau_sigma_nu2(2)   = 0.0030465d0

! values for Paul Cristini for fluid-solid ocean acoustics simulations

! for N_SLS = 2
! frequency range: 1.500000 Hz - 18.000000 Hz
! central frequency in log scale in Hz = 5.196152422706633
! target constant attenuation factor Q = 136.4376068115
! tau sigma evenly spaced in log frequency, do not depend on value of Q

! tau_sigma_nu1(1) = 0.10610329539459699422d0
! tau_sigma_nu1(2) = 0.00884194128288308401d0

! tau_epsilon_nu1(1) = 0.10754721280605997191d0
! tau_epsilon_nu1(2) = 0.00895488050110176612d0

! tau_epsilon_nu2(1) = tau_epsilon_nu1(1)
! tau_epsilon_nu2(2) = tau_epsilon_nu1(2)
! tau_sigma_nu2(1)   = tau_sigma_nu1(1)
! tau_sigma_nu2(2)   = tau_sigma_nu1(2)

!
!--- other constants computed from the parameters above, do not modify
!
  inv_tau_sigma_nu1(:) = ONE / tau_sigma_nu1(:)
  inv_tau_sigma_nu2(:) = ONE / tau_sigma_nu2(:)

  phi_nu1(:) = (ONE - tau_epsilon_nu1(:)/tau_sigma_nu1(:)) / tau_sigma_nu1(:)
  phi_nu2(:) = (ONE - tau_epsilon_nu2(:)/tau_sigma_nu2(:)) / tau_sigma_nu2(:)

  Mu_nu1 = ONE
  Mu_nu2 = ONE

  do i_sls = 1,N_SLS
    Mu_nu1 = Mu_nu1 - (ONE - tau_epsilon_nu1(i_sls)/tau_sigma_nu1(i_sls))
    Mu_nu2 = Mu_nu2 - (ONE - tau_epsilon_nu2(i_sls)/tau_sigma_nu2(i_sls))
  enddo

  end subroutine attenuation_model

