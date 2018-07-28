!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

!----
!---- subroutine to compute poroelastic velocities cpI, cpII, & cs as a function of the dominant frequency
!----

!
! This routine calculates the poroelastic wave speeds, cpI, cpII and cs.
!
! We have two frequencies:
!
! - fi is the dominant frequency of the source time function
!
! - f0 is a frequency used when we consider the frequency dependence of the viscous term.
!   Namely without viscous attenuation b = eta_f k^-1
!   While when we do introduce attenuation mechanism using memory variable b=b(t) and the expression
!   is that of our paper equation (161) - Zener model. It is the same model used for viscoelasticity in the code.
!
!   Again, fi is set by the source time function.
!   f0 is set in the Par_file as freq0_poroelastic, as well as Q0_poroelastic, to define the strain and stress relaxation times.
!
!   If you are running simulations without the fluid viscous dissipation,
!   then in Par_file you have the flag ATTENUATION_PORO_FLUID_PART     = .false.
!   and fi and Q0 are obsolete, and b=eta_f k^-1 simply.
!
!   If the flag is set to true, then you have to provide the quality factor and frequency range you want to work with.
!
!   One can refer to Section 8.3 of Morency, C. and Tromp, J., Spectral-element simulations of wave propagation
!   in poroelastic media, Geophys. J. Int., vol. 175, p. 301-345 (2008).
!

  subroutine get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare, &
                                        H_biot,C_biot,M_biot,mu_fr,phi, &
                                        tort,rho_s,rho_f,eta_f,perm_xx, &
                                        fi,f0,Q0_poroelastic,w_c,ATTENUATION_PORO_FLUID_PART)

  use constants, only: PI

  implicit none

  double precision,intent(out) :: cpIsquare,cpIIsquare,cssquare

  double precision,intent(in) :: H_biot,C_biot,M_biot
  double precision,intent(in) :: eta_f,rho_f,rho_s,perm_xx
  double precision,intent(in) :: mu_fr,phi,tort

  double precision,intent(in) :: fi,f0,Q0_poroelastic
  double precision,intent(out) :: w_c

  logical,intent(in) :: ATTENUATION_PORO_FLUID_PART

  ! local parameters
  double precision :: rho_bar
  double precision :: a_r,a_i,b_r,b_i,cc,alpha,aa1,aa2
  double precision :: xx,yy, gXI, gYI,gXII,gYII,f_c
  double precision :: taus,taue,bbr,bbi
  double precision :: wi,w0il
  double precision :: gA,gB,sa,sb,xxs,yys
  double precision :: att_I,att_II

  rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

  w_c = eta_f*phi/(tort*rho_f*perm_xx)

  f_c = w_c/(2*PI)

  wi = 2.d0 * PI * fi
  w0il =  2.d0 * PI * f0

  alpha = 10.d0**dlog10(wi)

  taue = (sqrt(Q0_poroelastic*Q0_poroelastic+1) + 1)/(w0il*Q0_poroelastic)
  taus = (sqrt(Q0_poroelastic*Q0_poroelastic+1) - 1)/(w0il*Q0_poroelastic)

  if (ATTENUATION_PORO_FLUID_PART) then
    ! high frequency, with memory variables
    bbr = eta_f/perm_xx * (1.d0+alpha*alpha*taus*taue)/(1.d0 + alpha*alpha*taus*taus)
    bbi = eta_f/perm_xx * alpha*(taue-taus)/(1.d0 + alpha*alpha*taus*taus)
  else
    ! low frequency
    bbr = eta_f/perm_xx
    bbi = 0.d0
  endif

  ! cs
  gA = (rho_f*tort*rho_bar-phi*rho_f**2)**2/(phi*rho_bar)**2 - (bbr**2-bbi**2)/alpha**2*&
      (phi*rho_f/(rho_bar*tort) -1.d0) - bbi/alpha*phi*rho_f/(rho_bar*tort)*&
      (rho_f*tort*rho_bar-phi*rho_f**2)/(phi*rho_bar)
  gB = -2.d0*bbr*bbi/alpha**2*(phi*rho_f/(rho_bar*tort) -1.d0) + bbr/alpha*phi*rho_f/&
      (rho_bar*tort)*(rho_f*tort*rho_bar-phi*rho_f**2)/(phi*rho_bar)
  !
  sa = (rho_f*tort*rho_bar-phi*rho_f**2)**2/(phi*rho_bar)**2 + (bbr**2-bbi**2)/alpha**2
  sb = 2.d0*bbr*bbi/alpha**2
  !
  xxs = sa*gA + sb*gB
  yys = gA*sb - sa*gB

  cssquare = mu_fr/(rho_bar-phi*rho_f/tort) * 2.d0*(gA**2+gB**2)/(sqrt(xxs**2+yys**2)+xxs)


  ! cpI & cpII
  a_r = rho_bar - phi*rho_f/tort - phi*rho_bar/(tort*rho_f)*bbi/alpha
  a_i = phi*rho_bar/(tort*rho_f)*bbr
  b_r = H_biot + M_biot*phi*rho_bar/(tort*rho_f) - 2.d0*phi*C_biot/tort - &
      phi*H_biot/(tort*rho_f)*bbi/alpha
  b_i = phi*H_biot/(tort*rho_f)*bbr
  cc = phi/(tort*rho_f)*(H_biot*M_biot - C_biot*C_biot)
  !
  xx = b_r*b_r - b_i*b_i/(alpha*alpha) - 4.d0*a_r*cc
  yy = 2.d0*b_r*b_i/alpha - 4.d0*a_i/alpha*cc
  !
  gXI = a_r*(b_r + sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx)) + &
        a_i/alpha*(b_i/alpha + sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))

  gYI = a_i/alpha*(b_r + sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx)) - &
        a_r*(b_i/alpha + sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))
  gYI = -gYI

  gXII = a_r*(b_r - sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx)) + &
        a_i/alpha*(b_i/alpha - sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))

  gYII = a_i/alpha*(b_r - sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx)) - &
        a_r*(b_i/alpha - sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))
  gYII = -gYII
  !
  !
  !
  cpIsquare = ((b_r + sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx))**2 + &
              (b_i/alpha + sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))**2)/&
              (sqrt(gXI**2+gYI**2) + gXI)

  cpIIsquare = ((b_r - sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx))**2 + &
              (b_i/alpha - sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))**2)/&
              (sqrt(gXII**2+gYII**2) + gXII)

  ! attenuation factors
  att_I = -alpha*sign(1.d0,yy)*sqrt(sqrt(gXI**2+gYI**2)-gXI) / &
           sqrt((b_r + sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx))**2+&
               (b_i/alpha + sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))**2)
  att_II = -alpha*sign(1.d0,yy)*sqrt(sqrt(gXII**2+gYII**2)-gXII) / &
           sqrt((b_r - sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx))**2+&
               (b_i/alpha - sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))**2)

  ! inverse quality factors
  aa1 = -gYI/gXI
  aa2 = -gYII/gXII

  end subroutine get_poroelastic_velocities


!
!-------------------------------------------------------------------------------
!

  subroutine get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)

  use constants, only: FOUR_THIRDS

  use specfem_par, only: porosity,tortuosity,poroelastcoef,density,kmato ! AXISYM

  implicit none

  integer,intent(in) :: ispec

  double precision,intent(out) :: phi,tort
  double precision,intent(out) :: kappa_s,kappa_f,kappa_fr,mu_s,mu_fr,rho_s,rho_f,rho_bar,eta_f

  ! local parameters
  integer :: material

  ! gets associated material
  material = kmato(ispec)

  ! gets material properties
  phi = porosity(material)
  tort = tortuosity(material)

  ! solid properties
  mu_s = poroelastcoef(2,1,material)
  !if (AXISYM) then ! ABAB !! Warning !! This is false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
    kappa_s = poroelastcoef(3,1,material) - FOUR_THIRDS * mu_s
  !else
  !kappa_s = poroelastcoef(3,1,material) - mu_s
  !endif
  rho_s = density(1,material)

  ! fluid properties
  kappa_f = poroelastcoef(1,2,material)
  rho_f = density(2,material)
  eta_f = poroelastcoef(2,2,material)

  ! frame properties
  mu_fr = poroelastcoef(2,3,material)
  !if (AXISYM) then ! ABAB !! Warning !! This is false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
    kappa_fr = poroelastcoef(3,3,material) - FOUR_THIRDS * mu_fr
  !else
  !  kappa_fr = poroelastcoef(3,3,material) - mu_fr
  !endif
  ! rho bar
  rho_bar =  (1.d0 - phi) * rho_s + phi * rho_f

  end subroutine get_poroelastic_material

!
!-------------------------------------------------------------------------------
!

  subroutine get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

  ! use specfem_par, only: AXISYM
  use constants, only: FOUR_THIRDS

  implicit none

  double precision,intent(in) :: phi,kappa_s,kappa_f,kappa_fr,mu_fr
  double precision,intent(out) :: D_biot,H_biot,C_biot,M_biot

  ! Biot coefficients for the input phi
  D_biot = kappa_s*(1.d0 + phi*(kappa_s/kappa_f - 1.d0))
  !if (AXISYM) then ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
    H_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr + FOUR_THIRDS*mu_fr
  !else
  !  H_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr + mu_fr
  !endif
  C_biot = kappa_s*(kappa_s - kappa_fr)/(D_biot - kappa_fr)
  M_biot = kappa_s*kappa_s/(D_biot - kappa_fr)

  end subroutine get_poroelastic_Biot_coeff
