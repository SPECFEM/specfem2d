
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
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

  subroutine attenuation_model(N_SLS,QKappa_attenuation,Qmu_attenuation,f0_attenuation, &
       tau_epsilon_nu1_custom_real,tau_epsilon_nu2_custom_real,inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,Mu_nu1,Mu_nu2)

! define the attenuation constants

  implicit none

  include "constants.h"

  integer :: N_SLS
  double precision :: QKappa_attenuation,Qmu_attenuation,f0_attenuation
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: tau_epsilon_nu1_custom_real,tau_epsilon_nu2_custom_real, &
              inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2
  real(kind=CUSTOM_REAL) :: Mu_nu1,Mu_nu2

  integer :: i_sls

  double precision, dimension(N_SLS) :: tau_epsilon_nu1,tau_sigma_nu1,tau_epsilon_nu2,tau_sigma_nu2

  double precision :: f1_attenuation, f2_attenuation

! f1 and f2 are computed as : f2/f1=12 and (log(f1)+log(f2))/2 = log(f0)
  f1_attenuation = exp(log(f0_attenuation)-log(12.d0)/2.d0)
  f2_attenuation = 12.d0 * f1_attenuation

! Call of C function that computes attenuation parameters (function in file "attenuation_compute_param.c";
! a main can be found in UTILS/attenuation directory).
  call attenuation_compute_param(N_SLS, QKappa_attenuation, Qmu_attenuation, &
       f1_attenuation,f2_attenuation, &
       tau_sigma_nu1, tau_sigma_nu2, tau_epsilon_nu1, tau_epsilon_nu2)

! attenuation constants for standard linear solids

! nu1 is the dilatation/incompressibility mode (QKappa)
! nu2 is the shear mode (Qmu)

! array index (1) is the first standard linear solid, (2) is the second etc.

! When implementing viscoelasticity according to the Carcione 1993 paper, attenuation is
! non-causal rather than causal. We fixed the problem by using equations in Carcione's
! 2004 paper and his 2007 book. See also file doc/problem_attenuation_reference_Specfem2D_fixed_by_Xie_Zhinan.pdf

! from J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).
! Beware: these values implement specific values of the quality factors:
! Qp approximately equal to 13, Qkappa approximately to 20 and Qmu / Qs approximately to 10,
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
! Qkappa approximately to 27 and Qmu / Qs approximately to 20,
! which means very high attenuation, see that paper for details.
!  tau_epsilon_nu1(1) = 0.0325305d0
!  tau_sigma_nu1(1)   = 0.0311465d0
!  tau_epsilon_nu2(1) = 0.0332577d0
!  tau_sigma_nu2(1)   = 0.0304655d0

!  tau_epsilon_nu1(2) = 0.0032530d0
!  tau_sigma_nu1(2)   = 0.0031146d0
!  tau_epsilon_nu2(2) = 0.0033257d0
!  tau_sigma_nu2(2)   = 0.0030465d0

!
!--- other constants computed from the parameters above, do not modify
!
  if(CUSTOM_REAL == SIZE_REAL) then

   tau_epsilon_nu1_custom_real(:) = sngl(tau_epsilon_nu1(:))
   tau_epsilon_nu2_custom_real(:) = sngl(tau_epsilon_nu2(:))

   inv_tau_sigma_nu1(:) = sngl(dble(ONE) / tau_sigma_nu1(:))
   inv_tau_sigma_nu2(:) = sngl(dble(ONE) / tau_sigma_nu2(:))

   phi_nu1(:) = sngl((dble(ONE) - tau_epsilon_nu1(:)/tau_sigma_nu1(:)) / tau_sigma_nu1(:))
   phi_nu2(:) = sngl((dble(ONE) - tau_epsilon_nu2(:)/tau_sigma_nu2(:)) / tau_sigma_nu2(:))

   Mu_nu1 = dble(ONE)
   Mu_nu2 = dble(ONE)

   do i_sls = 1,N_SLS
     Mu_nu1 = sngl(dble(Mu_nu1) - (dble(ONE) - tau_epsilon_nu1(i_sls)/tau_sigma_nu1(i_sls)))
     Mu_nu2 = sngl(dble(Mu_nu2) - (dble(ONE) - tau_epsilon_nu2(i_sls)/tau_sigma_nu2(i_sls)))
   enddo

  else

   tau_epsilon_nu1_custom_real(:) = tau_epsilon_nu1(:)
   tau_epsilon_nu2_custom_real(:) = tau_epsilon_nu2(:)

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

  endif

  end subroutine attenuation_model

!--------------------------------------------------------------------------------

  subroutine shift_velocities_from_f0(vp,vs,rho,mu,lambda,f0_attenuation, &
                              tau_epsilon_nu1,tau_epsilon_nu2,inv_tau_sigma_nu1,inv_tau_sigma_nu2,N_SLS)

! From Emmanuel Chaljub, CNRS Grenoble, France:

! shift (i.e. change) velocities read from the input file to take average physical dispersion into account,
! i.e. if needed change the reference frequency at which these velocities are defined internally in the code:

!  by default, the velocity values that are read in the Par_file of the code are supposed to be the unrelaxed values,
!  i.e. the velocities at infinite frequency.
!  We may want to change this and impose that the values read are those for a given frequency (here f0_attenuation).
!
!  The unrelaxed values are then defined from the values read at frequency f0_attenuation as follows:
!
!     mu_unrelaxed = mu (w_ref) / [ 1 - ( Sum_k ak/(1+(w_ref*tau_k)**2) )/( 1 + Sum_k ak ) ]
!     where ak = tau_epsilon_k/tau_sigma_k - 1
!     and the ak are the solutions of the linear system:
!
!     Sum_k   { [  w*tau_k*Q(w) - (w*tau_k)^2 ] / [ 1 + (w*tau_k)^2 ] }  ak  = 1
!                  where tau_k = tau_epsilon_k

  implicit none

  include "constants.h"

! arguments
  integer, intent(in) :: N_SLS
  double precision, intent(in) :: rho,f0_attenuation
  double precision, intent(inout) :: vp,vs
  double precision, intent(out) :: mu,lambda
  real(kind=CUSTOM_REAL), dimension(N_SLS), intent(in) :: tau_epsilon_nu1,tau_epsilon_nu2,inv_tau_sigma_nu1,inv_tau_sigma_nu2

! local variables
  integer :: i_sls
  double precision :: xtmp1_nu1,xtmp1_nu2,xtmp2_nu1,xtmp2_nu2,xtmp_ak_nu1,xtmp_ak_nu2,delta_mu,kappa,delta_kappa

  mu = rho * vs*vs
  lambda = rho * vp*vp - TWO * mu
  kappa  = lambda + TWO_THIRDS*mu

  xtmp1_nu1 = ZERO
  xtmp2_nu1 = ONE
  xtmp1_nu2 = ZERO
  xtmp2_nu2 = ONE

  do i_sls = 1,N_SLS
!! DK DK changed this to the pre-computed inverse     xtmp_ak_nu2 = tau_epsilon_nu2(i_sls)/tau_sigma_nu2(i_sls) - ONE
     xtmp_ak_nu2 = tau_epsilon_nu2(i_sls)*inv_tau_sigma_nu2(i_sls) - ONE
     xtmp1_nu2 = xtmp1_nu2 + xtmp_ak_nu2/(ONE + (TWO * PI * f0_attenuation * tau_epsilon_nu2(i_sls))**2)
     xtmp2_nu2 = xtmp2_nu2 + xtmp_ak_nu2

!! DK DK changed this to the pre-computed inverse     xtmp_ak_nu1 = tau_epsilon_nu1(i_sls)/tau_sigma_nu1(i_sls) - ONE
     xtmp_ak_nu1 = tau_epsilon_nu1(i_sls)*inv_tau_sigma_nu1(i_sls) - ONE
     xtmp1_nu1 = xtmp1_nu1 + xtmp_ak_nu1/(ONE + (TWO * PI * f0_attenuation * tau_epsilon_nu1(i_sls))**2)
     xtmp2_nu1 = xtmp2_nu1 + xtmp_ak_nu1
  enddo

  delta_mu = ONE - xtmp1_nu2/xtmp2_nu2
  mu    = mu    / delta_mu

  delta_kappa = ONE - xtmp1_nu1/xtmp2_nu1
  kappa = kappa / delta_kappa

  lambda = kappa - TWO_THIRDS*mu
  vp = dsqrt((lambda + TWO * mu) / rho)
  vs = dsqrt(mu / rho)

  end subroutine shift_velocities_from_f0

