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

  subroutine attenuation_model(QKappa_att,QMu_att,ATTENUATION_f0_REFERENCE,N_SLS, &
                               tau_epsilon_nu1_sent,inv_tau_sigma_nu1_sent,phi_nu1_sent,Mu_nu1_sent, &
                               tau_epsilon_nu2_sent,inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu2_sent)

! define the attenuation constants

  use constants, only: CUSTOM_REAL,ONE,USE_OLD_C_ATTENUATION_ROUTINE_INSTEAD
  use shared_input_parameters, only: USE_SOLVOPT

  implicit none

  double precision,intent(in) :: QKappa_att,QMu_att
  double precision,intent(in) :: ATTENUATION_f0_REFERENCE

  integer,intent(in) :: N_SLS
  real(kind=CUSTOM_REAL), dimension(N_SLS),intent(out) :: inv_tau_sigma_nu1_sent,phi_nu1_sent
  real(kind=CUSTOM_REAL), dimension(N_SLS),intent(out) :: inv_tau_sigma_nu2_sent,phi_nu2_sent
  real(kind=CUSTOM_REAL), dimension(N_SLS),intent(out) :: tau_epsilon_nu1_sent,tau_epsilon_nu2_sent
  real(kind=CUSTOM_REAL),intent(out) :: Mu_nu1_sent,Mu_nu2_sent

  ! local parameters
  double precision, dimension(N_SLS) :: tau_epsilon_nu1_d,tau_epsilon_nu2_d
  double precision, dimension(N_SLS) :: tau_sigma_nu1,tau_sigma_nu2
  double precision :: f_min_attenuation, f_max_attenuation

  ! safety check
  if (N_SLS < 1) call stop_the_code('Invalid N_SLS value, must be at least 1')

! attenuation constants for standard linear solids
! nu1 is the dilatation/incompressibility mode (QKappa)
! nu2 is the shear mode (Qmu)
! array index (1) is the first standard linear solid, (2) is the second etc.

! can instead call the old C routine when USE_SOLVOPT is false, for compatibility with older benchmarks
  if (USE_OLD_C_ATTENUATION_ROUTINE_INSTEAD .and. .not. USE_SOLVOPT) then

  ! f_min and f_max are computed as : f_max/f_min=12 and (log(f_min)+log(f_max))/2 = log(f0)
    f_min_attenuation = exp(log(ATTENUATION_f0_REFERENCE)-log(12.d0)/2.d0)
    f_max_attenuation = 12.d0 * f_min_attenuation

  ! call of C function that computes attenuation parameters (function in file "attenuation_compute_param.c";
  ! a main can be found in UTILS/attenuation directory).
    call attenuation_compute_param(N_SLS,QKappa_att,QMu_att,f_min_attenuation,f_max_attenuation, &
                                   tau_sigma_nu1,tau_sigma_nu2,tau_epsilon_nu1_d,tau_epsilon_nu2_d)

  else

  ! use a wide bandwidth (always OK when using three or more Standard Linear Solids, can be a bit inaccurate if using only two)
    f_min_attenuation = ATTENUATION_f0_REFERENCE / 10.d0
    f_max_attenuation = ATTENUATION_f0_REFERENCE * 10.d0

    call compute_attenuation_coeffs(N_SLS,QKappa_att,ATTENUATION_f0_REFERENCE,f_min_attenuation,f_max_attenuation, &
                                          tau_epsilon_nu1_d,tau_sigma_nu1)

    call compute_attenuation_coeffs(N_SLS,QMu_att,ATTENUATION_f0_REFERENCE,f_min_attenuation,f_max_attenuation, &
                                          tau_epsilon_nu2_d,tau_sigma_nu2)

  endif

! print *
! print *,'N_SLS, QKappa_att, QMu_att = ',N_SLS, QKappa_att, QMu_att
! print *,'ATTENUATION_f0_REFERENCE,f_min_attenuation,f_max_attenuation = ',ATTENUATION_f0_REFERENCE, &
!                              f_min_attenuation,f_max_attenuation
! print *,'tau_epsilon_nu1 = ',tau_epsilon_nu1_d
! print *,'tau_sigma_nu1 = ',tau_sigma_nu1
! print *,'tau_epsilon_nu2 = ',tau_epsilon_nu2_d
! print *,'tau_sigma_nu2 = ',tau_sigma_nu2
! print *

  if (any(tau_sigma_nu1 < 0.d0) .or. any(tau_sigma_nu2 < 0.d0) .or. &
      any(tau_epsilon_nu1_d < 0.d0) .or. any(tau_epsilon_nu2_d < 0.d0)) &
       call stop_the_code('error: negative relaxation time found for a viscoelastic material')

! in the old formulation of Carcione 1993, which is based on Liu et al. 1976, the 1/N factor is missing
! and thus this does not apply; it only applies to the right formula with 1/N included
  if (minval(tau_epsilon_nu1_d / tau_sigma_nu1) < 0.999d0 .or. &
                                    minval(tau_epsilon_nu2_d / tau_sigma_nu2) < 0.999d0) then
       print *
       print *,'*******************************************************************************'
       print *,'*******************************************************************************'
       print *,'*******************************************************************************'
       print *,'*******************************************************************************'
       print *,'*******************************************************************************'
       print *,'minval(tau_epsilon_nu1 / tau_sigma_nu1) = ',minval(tau_epsilon_nu1_d / tau_sigma_nu1)
       print *,'minval(tau_epsilon_nu2 / tau_sigma_nu2) = ',minval(tau_epsilon_nu2_d / tau_sigma_nu2)
       print *,'WARNING: tau_epsilon should never be smaller than tau_sigma for viscoelasticity'
       print *,'*******************************************************************************'
       print *,'*******************************************************************************'
       print *,'*******************************************************************************'
       print *,'*******************************************************************************'
       print *,'*******************************************************************************'
  endif

! When implementing viscoelasticity according to the Carcione 1993 paper, attenuation is
! non-causal rather than causal i.e. wave speed up instead of slowing down
! when attenuation is turned on. We fixed that issue (which is not incorrect but non traditional)
! by taking the unrelaxed state (infinite frequency) as a reference instead of the relaxed state (zero frequency)
! and also using equations in Carcione's 2007 book.
! See file doc/old_problem_attenuation_reference_Specfem2D_fixed_by_Xie_Zhinan.pdf
! and doc/how_we_modified_Carcione_1993_to_make_it_causal_and_include_the_missing_1_over_L_factor.pdf

! See also J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).

! and J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
! in a linear viscoelastic medium, Geophysical Journal International,
! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).

!
!--- other constants computed from the parameters above, do not modify
!

  ! SLS parameters
  ! for Qkappa (set by tau**1,phi**1,Mu_nu1**) and Qmu (set by tau**2,phi**2,Mu_nu2**)
  tau_epsilon_nu1_sent(:) = real(tau_epsilon_nu1_d(:),kind=CUSTOM_REAL)
  tau_epsilon_nu2_sent(:) = real(tau_epsilon_nu2_d(:),kind=CUSTOM_REAL)

  inv_tau_sigma_nu1_sent(:) = real(1.d0 / tau_sigma_nu1(:),kind=CUSTOM_REAL)
  inv_tau_sigma_nu2_sent(:) = real(1.d0 / tau_sigma_nu2(:),kind=CUSTOM_REAL)

  ! use the right formula with 1/N included
  phi_nu1_sent(:) = real((1.d0 - tau_epsilon_nu1_d(:)/tau_sigma_nu1(:)) / tau_sigma_nu1(:) &
                                                                           / sum(tau_epsilon_nu1_d/tau_sigma_nu1),kind=CUSTOM_REAL)
  phi_nu2_sent(:) = real((1.d0 - tau_epsilon_nu2_d(:)/tau_sigma_nu2(:)) / tau_sigma_nu2(:) &
                                                                           / sum(tau_epsilon_nu2_d/tau_sigma_nu2),kind=CUSTOM_REAL)

  Mu_nu1_sent = real(sum(tau_epsilon_nu1_d/tau_sigma_nu1) / dble(N_SLS),kind=CUSTOM_REAL)
  Mu_nu2_sent = real(sum(tau_epsilon_nu2_d/tau_sigma_nu2) / dble(N_SLS),kind=CUSTOM_REAL)

  if (Mu_nu1_sent < 1. .or. Mu_nu2_sent < 1.) &
    call stop_the_code('error in Zener viscoelasticity: must have Mu_nu1 and Mu_nu2 both greater than one')

  end subroutine attenuation_model

!
!--------------------------------------------------------------------------------
!

  subroutine shift_velocities_from_f0(vp,vs,rho, &
                                      ATTENUATION_f0_REFERENCE,N_SLS, &
                                      tau_epsilon_nu1,tau_epsilon_nu2,inv_tau_sigma_nu1,inv_tau_sigma_nu2)

! From Emmanuel Chaljub, ISTerre, OSU Grenoble, France:

! shift (i.e. change) velocities read from the input file to take average physical dispersion into account,
! i.e. if needed change the reference frequency at which these velocities are defined internally in the code:

!  by default, the velocity values that are read in the Par_file of the code are supposed to be the unrelaxed values,
!  i.e. the velocities at infinite frequency.
!  We may want to change this and impose that the values read are those for a given frequency (here ATTENUATION_f0_REFERENCE).
!
!  The unrelaxed values are then defined from the reference values read at frequency ATTENUATION_f0_REFERENCE as follows:
!
!     mu_unrelaxed = mu (w_ref) * [ ( 1 + (1/N) Sum_k ak ) / (1 + (1/N) Sum_k ak/(1+1/(w_ref*tau_sigma_k)**2) ) ]
!
!     where w_ref = 2*PI*ATTENUATION_f0_REFERENCE
!           tau_k = tau_epsilon_k is the strain relaxation time of the k-th SLS mechanism
!              ak = tau_k/tau_sigma_k - 1
!     The ak are the solutions of the linear system:
!     (1/N) Sum_k   { [  w*tau_k*Q(w) - (w*tau_k)^2 ] / [ 1 + (w*tau_k)^2 ] }  ak  = 1
!
!  To see how to obtain these formulas, see for instance equations (8), (9) and (10) of
!  P. Moczo, E. Bystricky, J. Kristek, J. M. Carcione and M. Bouchon,
!  Hybrid modeling of P-SV seismic motion at inhomogeneous viscoelastic topographic structures,
!  Bulletin of the Seismological Society of Americal, vol. 87, p. 1305-1323 (1997).
!
!  See also some of the PDF files in the "doc/" directory of the code.
!
!  The above formulas are for a Zener-body approximation of a Q(f) model, which is what SPECFEM uses.
!  If one wants to use exact formulas for a truly constant Q model instead (i.e., not approximated by a set of Zener bodies),
!  the expression of the scaling can be found for instance in equation (49) of Komatitsch and Tromp,
!  Spectral-Element Simulations of Global Seismic Wave Propagation-I. Validation, Geophys. J. Int. vol. 149 p. 390-412 (2002),
!  which is also in file "doc/exact_formula_to_scale_mu_from_a_frequency_f0_to_a_frequency_f1_from_Komatitsch_Tromp_GJI_2002.pdf":
!
!       mu(omega_c) = mu(omega_0) * (1 + 2 * ln(omega_c/omega_0) / (PI * Q_mu))
!
!  (the formula to scale Vs is the same except for the factor of 2, which needs to be removed because mu is related to Vs squared.
!
!  A similar expression can then be established for Q_Kappa, and conversion from Q_Kappa and Q_mu to Q_P and Q_S (if needed)
!  can be found for instance in file "Qkappa_Qmu_versus_Qp_Qs_relationship_in_2D_plane_strain.pdf"
!  in the "doc" directory of the code.

  use constants, only: ONE,TWO,PI,TWO_THIRDS,CUSTOM_REAL

  use specfem_par, only: AXISYM

  implicit none

! arguments
  double precision, intent(inout) :: vp,vs

  double precision, intent(in) :: rho
  double precision, intent(in) :: ATTENUATION_f0_REFERENCE

  integer,intent(in) :: N_SLS
  real(kind=CUSTOM_REAL), dimension(N_SLS),intent(in) :: tau_epsilon_nu1,tau_epsilon_nu2
  real(kind=CUSTOM_REAL), dimension(N_SLS),intent(in) :: inv_tau_sigma_nu1,inv_tau_sigma_nu2

! local variables
  integer :: i_sls
  double precision :: xtmp1_nu1,xtmp1_nu2,xtmp2_nu1,xtmp2_nu2,xtmp_ak_nu1,xtmp_ak_nu2
  double precision :: factor_mu,factor_kappa
  double precision :: kappa,mu,lambda

  mu = rho * vs*vs
  lambda = rho * vp*vp - TWO * mu

  if (AXISYM) then ! CHECK kappa
    kappa  = lambda + TWO_THIRDS*mu
  else
    kappa  = lambda + mu
  endif

  xtmp1_nu1 = ONE
  xtmp2_nu1 = ONE
  xtmp1_nu2 = ONE
  xtmp2_nu2 = ONE

  do i_sls = 1,N_SLS
!! DK DK changed this to the pre-computed inverse     xtmp_ak_nu2 = tau_epsilon_nu2(i_sls)/tau_sigma_nu2(i_sls) - ONE
     xtmp_ak_nu2 = tau_epsilon_nu2(i_sls)*inv_tau_sigma_nu2(i_sls) - ONE
     xtmp1_nu2 = xtmp1_nu2 + xtmp_ak_nu2/N_SLS
     xtmp2_nu2 = xtmp2_nu2 + xtmp_ak_nu2/(ONE + ONE/(TWO * PI * ATTENUATION_f0_REFERENCE / inv_tau_sigma_nu2(i_sls))**2)/N_SLS

!! DK DK changed this to the pre-computed inverse     xtmp_ak_nu1 = tau_epsilon_nu1(i_sls)/tau_sigma_nu1(i_sls) - ONE
     xtmp_ak_nu1 = tau_epsilon_nu1(i_sls)*inv_tau_sigma_nu1(i_sls) - ONE
     xtmp1_nu1 = xtmp1_nu1 + xtmp_ak_nu1/N_SLS
     xtmp2_nu1 = xtmp2_nu1 + xtmp_ak_nu1/(ONE + ONE/(TWO * PI * ATTENUATION_f0_REFERENCE / inv_tau_sigma_nu1(i_sls))**2)/N_SLS
  enddo

  factor_mu = xtmp1_nu2/xtmp2_nu2
  mu    = mu    * factor_mu

  factor_kappa = xtmp1_nu1/xtmp2_nu1
  kappa = kappa * factor_kappa

  if (AXISYM) then ! CHECK kappa
    lambda = kappa - TWO_THIRDS*mu
  else
    lambda = kappa - mu
  endif

  ! returns shifted vp,vs
  vp = dsqrt((lambda + TWO * mu) / rho)
  vs = dsqrt(mu / rho)

  end subroutine shift_velocities_from_f0

!
!--------------------------------------------------------------------------------
!

! use of SolvOpt to compute attenuation relaxation mechanisms,
! from Emilie Blanc, Bruno Lombard and Dimitri Komatitsch, CNRS Marseille, France, for a Generalized Zener body model.

! If you use this code for your own research, please cite some (or all) of these articles:
!
! @Article{BlKoChLoXi16,
! Title   = {Highly accurate stability-preserving optimization of the {Z}ener viscoelastic model,
!            with application to wave propagation in the presence of strong attenuation},
! Author  = {\'Emilie Blanc and Dimitri Komatitsch and Emmanuel Chaljub and Bruno Lombard and Zhinan Xie},
! Journal = {Geophysical Journal International},
! Year    = {2016},
! Number  = {1},
! Pages   = {427-439},
! Volume  = {205},
! Doi     = {10.1093/gji/ggw024}}

! The SolvOpt algorithm was developed by Franz Kappel and Alexei V. Kuntsevich
! and is available open source at http://www.uni-graz.at/imawww/kuntsevich/solvopt
!
! It is described in Kappel and Kuntsevich, An Implementation of Shor's r-Algorithm,
! Computational Optimization and Applications, vol. 15, p. 193-205 (2000).

!-------------------------------------------------------------------------

! From Bruno Lombard, May 2014:

! En interne dans le code ci-dessous on travaille en (Theta, Kappa).
! Les Theta sont les points et les Kappa sont les poids.
! Pour repasser en (Tau_Sigma, Tau_Epsilon), on doit appliquer les formules:
!
! Tau_Sigma = 1 / Theta
! Tau_Epsilon = (1 / Theta) * (1 + Nrelax * Kappa) = Tau_Sigma * (1 + Nrelax * Kappa)

! The system to solve can be found in equation (7) of:
! Lombard and Piraux, Numerical modeling of transient two-dimensional viscoelastic waves,
! Journal of Computational Physics, Volume 230, Issue 15, Pages 6099-6114 (2011)

! Suivant les compilateurs et les options de compilation utilisees,
! il peut y avoir des differences au 4eme chiffre significatif. C'est sans consequences sur la precision du calcul :
! l'erreur est de 0.015 % avec optimization non-lineaire, a comparer a 1.47 % avec Emmerich and Korn (1987).
! Si je relance le calcul en initialisant avec le resultat precedent, ce chiffre varie a nouveau tres legerement.

!-------------------------------------------------------------------------

! From Bruno Lombard, June 2014:

! j'ai relu en detail :

! [1] Carcione, Kosslof, Kosslof, "Viscoacoustic wave propagation simulation in the Earth",
!            Geophysics 53-6 (1988), 769-777
!
! [2] Carcione, Kosslof, Kosslof, "Wave propagation simulation in a linear viscoelastic medium",
!            Geophysical Journal International 95 (1988), 597-611
!
! [3] Moczo, Kristek, "On the rheological models used for time-domain methods of seismic wave propagation",
!            Geophysical Research Letters 32 (2005).

! Le probleme provient probablement d'une erreur recurrente dans [1,2] et datant de Liu et al 1976 :
! l'oubli du facteur 1/N dans la fonction de relaxation d'un modele de Zener a N elements.
! Il est effectivement facile de faire l'erreur. Voir l'equation (12) de [3], et le paragraphe qui suit.

! Du coup le module de viscoelasticite est faux dans [1,2], et donc le facteur de qualite,
! et donc les temps de relaxation tau_sigma...

! Apres, [2] calcule une solution analytique juste, mais avec des coefficients sans sens physique.
! Et quand SPECFEM2D obtient un bon accord avec cette solution analytique, ca valide SPECFEM, mais pas le calcul des coefficients.

! Il y a donc une erreur dans [1,2], et [3] a raison.

! Sa solution analytique decoule d'un travail sur ses fonctions de relaxation (A4),
! qu'il injecte ensuite dans la relation de comportement (A1) et travaille en Fourier.

! Le probleme est que sa fonction de relaxation (A4) est fausse : il manque 1/N.
! De ce fait, sa solution analytique est coherente avec sa solution numerique.
! Dans les deux cas, ce sont les memes temps de relaxation qui sont utilises. Mais ces temps sont calcules de facon erronee.

!-------------------------------------------------------------------------

! From Dimitri Komatitsch, June 2014:

! In [2] Carcione, Kosslof, Kosslof, "Wave propagation simulation in a linear viscoelastic medium",
!            Geophysical Journal International 95 (1988), 597-611
! there is another mistake: in Appendix B page 611 Carcione writes omega/(r*v),
! but that is not correct, it should be omega*r/v instead.

!---------------------------------------------------

! From Emilie Blanc, April 2014:

! le programme SolvOpt d'optimization non-lineaire
! avec contrainte. Ce programme prend quatre fonctions en entree :

! - fun() est la fonction a minimiser

! - grad() est le gradient de la fonction a minimiser par rapport a chaque parametre

! - func() est le maximum des residus (= 0 si toutes les contraintes sont satisfaites)

! - gradc() est le gradient du maximum des residus (= 0 si toutes les
! contraintes sont satisfaites)

! Ce programme a ete developpe par Kappel et Kuntsevich. Leur article decrit l'algorithme.

! J'ai utilise ce code pour la poroelasticite haute-frequence, et aussi en
! viscoelasticite fractionnaire (modele d'Andrade, avec Bruno Lombard et
! Cedric Bellis). Nous pouvons interagir sur l'algorithme d'optimization
! pour votre modele visco, et etudier l'effet des coefficients ainsi obtenus.

!---------------------------------------------------

! From Emilie Blanc, March 2014:

! Les entrees du programme principal sont le nombre de variables
! diffusives, le facteur de qualite voulu Qref et la frequence centrale f0.

! Cependant, pour l'optimization non-lineaire, j'ai mis theta_max=100*f0
! et non pas theta_max=2*pi*100*f0. En effet, dans le programme, on
! travaille sur les frequences, et non pas sur les frequences angulaires.
! Cela dit, dans les deux cas j'obtiens les memes coefficients...

!---------------------------------------------------

  subroutine compute_attenuation_coeffs(N,Qref,f0,f_min,f_max,tau_epsilon,tau_sigma)

  implicit none

  integer, intent(in) :: N
  double precision, intent(in) :: Qref,f_min,f_max,f0
  double precision, dimension(1:N), intent(out) :: tau_epsilon,tau_sigma

  integer :: i
  double precision, dimension(1:N) :: point,weight

! nonlinear optimization with constraints
  call nonlinear_optimization(N,Qref,f0,point,weight,f_min,f_max)

  do i = 1,N
    tau_sigma(i) = 1.d0 / point(i)
    tau_epsilon(i) = tau_sigma(i) * (1.d0 + N * weight(i))
  enddo

! print *,'points = '
! do i = 1,N
!   print *,point(i)
! enddo
! print *

! print *,'weights = '
! do i = 1,N
!   print *,weight(i)
! enddo
! print *

! print *,'tau_epsilon = '
! do i = 1,N
!   print *,tau_epsilon(i)
! enddo
! print *

! print *,'tau_sigma = '
! do i = 1,N
!   print *,tau_sigma(i)
! enddo
! print *

  end subroutine compute_attenuation_coeffs

!---------------------------------------------------

! classical calculation of the coefficients based on linear least squares

  subroutine decomposition_LU(a,i_min,n,indx,d)

  implicit none

  integer, intent(in) :: i_min,n
  double precision, intent(out) :: d
  integer, dimension(i_min:n), intent(inout) :: indx
  double precision, dimension(i_min:n,i_min:n), intent(inout) :: a

  integer i,imax,j,k
  double precision big,dum,somme,eps
  double precision, dimension(i_min:n) :: vv

  imax = 0
  d = 1.
  eps = 1.e-20

  do i = i_min,n
    big = 0.
    do j = i_min,n
      if (abs(a(i,j)) > big) then
        big = abs(a(i,j))
      endif
    enddo
    if (big == 0.) then
      print *,'Singular matrix in routine decomposition_LU'
    endif
    vv(i) = 1./big
  enddo

  do j = i_min,n
    do i = i_min,j-1
      somme = a(i,j)
      do k = i_min,i-1
        somme = somme - a(i,k)*a(k,j)
      enddo
      a(i,j) = somme
    enddo

    big = 0.

    do i = j,n
      somme = a(i,j)
      do k = i_min,j-1
        somme = somme - a(i,k)*a(k,j)
      enddo
      a(i,j) = somme
      dum = vv(i)*abs(somme)
      if (dum >= big) then
        big = dum
        imax = i
      endif
    enddo

    if (j /= imax) then
      do k = i_min,n
        dum = a(imax,k)
        a(imax,k) = a(j,k)
        a(j,k) = dum
      enddo
      d = -d
      vv(imax) = vv(j)
    endif

    indx(j) = imax
    if (a(j,j) == 0.) then
      a(j,j) = eps
    endif
    if (j /= n) then
      dum = 1./a(j,j)
      do i = j+1,n
        a(i,j) = a(i,j)*dum
      enddo
    endif
  enddo

  end subroutine decomposition_LU

  subroutine LUbksb(a,i_min,n,indx,b,m)

  implicit none

  integer, intent(in) :: i_min,n,m
  integer, dimension(i_min:n), intent(in) :: indx
  double precision, dimension(i_min:n,i_min:n), intent(in) :: a
  double precision, dimension(i_min:n,i_min:m), intent(inout) :: b

  integer i,ip,j,ii,k
  double precision somme

  do k = i_min,m

    ii = -1

    do i = i_min,n
      ip = indx(i)
      somme = b(ip,k)
      b(ip,k) = b(i,k)
      if (ii /= -1) then
        do j = ii,i-1
          somme = somme - a(i,j)*b(j,k)
        enddo
      else if (somme /= 0.) then
        ii = i
      endif
      b(i,k) = somme
    enddo

    do i = n,i_min,-1
      somme = b(i,k)
      do j = i+1,n
        somme = somme - a(i,j)*b(j,k)
      enddo
      b(i,k) = somme/a(i,i)
    enddo
  enddo

  end subroutine LUbksb

  subroutine syst_LU(a,i_min,n,b,m)

  implicit none

  integer, intent(in) :: i_min,n,m
  double precision, dimension(i_min:n,i_min:n), intent(in) :: a
  double precision, dimension(i_min:n,i_min:m), intent(inout) :: b

  integer i,j
  integer, dimension(i_min:n) :: indx
  double precision d
  double precision, dimension(i_min:n,i_min:n) :: aux

  do j = i_min,n
    indx(j) = 0
    do i = i_min,n
      aux(i,j) = a(i,j)
    enddo
  enddo

  call decomposition_LU(aux,i_min,n,indx,d)
  call LUbksb(aux,i_min,n,indx,b,m)

  end subroutine syst_LU

  subroutine lfit_zener(x,y,sig,ndat,poids,ia,covar,chisq,ma,Qref,point)
! ma = nombre de variable diffusive
! ndat = m = K nombre d'abcisse freq_k

  implicit none

  integer, intent(in) :: ndat,ma
  logical, dimension(1:ma), intent(in) :: ia
  double precision, intent(in) :: Qref
  double precision, intent(out) :: chisq
  double precision, dimension(1:ndat), intent(in) :: x,y,sig
  double precision, dimension(1:ma), intent(in) :: point
  double precision, dimension(1:ma), intent(out) :: poids
  double precision, dimension(1:ma,1:ma), intent(out) :: covar

  integer i,j,k,l,mfit
  double precision ym,wt,sig2i
  double precision, dimension(1:ma) :: afunc
  double precision, dimension(1:ma,1:1) :: beta

  mfit = 0

  do j = 1,ma
    if (ia(j)) then
      mfit = mfit + 1
    endif
  enddo
  if (mfit == 0) then
    print *,'lfit: no parameters to be fitted'
  endif

  do j = 1,mfit
    beta(j,1) = 0.
    do k = 1,mfit
      covar(j,k) = 0.
    enddo
  enddo

  do i = 1,ndat
    call func_zener(x(i),afunc,ma,Qref,point)
    ym = y(i)
    if (mfit < ma) then
      do j = 1,ma
        if (.not. ia(j)) then
          ym = ym - poids(j) * afunc(j)
        endif
      enddo
    endif
    sig2i = 1. / (sig(i) * sig(i))
    j = 0
    do l= 1,ma
      if (ia(l)) then
        j = j+1
        wt = afunc(l) * sig2i
        k = count(ia(1:l))
        covar(j,1:k) = covar(j,1:k) + wt * pack(afunc(1:l),ia(1:l))
        beta(j,1) = beta(j,1) + ym * wt
      endif
    enddo
  enddo

  do j = 2,mfit,1
  do k = 1,j-1,1
    covar(k,j) = covar(j,k)
  enddo
  enddo

  if (ma == 1) then
    poids(1) = beta(1,1)/covar(1,1)
  else if (ma > 1) then
    call syst_LU(covar,1,mfit,beta,1)
    poids(1:ma) = unpack(beta(1:ma,1),ia,poids(1:ma))
  endif

  chisq = 0.
  do i = 1,ndat
    call func_zener(x(i),afunc,ma,Qref,point)
    chisq=chisq+((y(i)-dot_product(poids(1:ma),afunc(1:ma)))/sig(i))**2
  enddo

  end subroutine lfit_zener

  subroutine func_zener(x,afunc,N,Qref,point)

  implicit none

  integer, intent(in) :: N
  double precision, intent(in) :: x,Qref
  double precision, dimension(1:N), intent(in) :: point
  double precision, dimension(1:N), intent(out) :: afunc

  integer k
  double precision num,deno

  do k = 1,N
    num  = x * (point(k) - x / Qref)
    deno = point(k) * point(k) + x * x
    afunc(k) = num / deno
  enddo

  end subroutine func_zener

  subroutine remplit_point(fmin,fmax,N,point)

  use constants, only: TWO_PI

  implicit none

  integer, intent(in) :: N
  double precision, intent(in) :: fmin,fmax
  double precision, dimension(1:N), intent(out) :: point

  integer l

  if (N == 1) then
    point(1) = sqrt(fmin * fmax)
  ELSE
    do l = 1, N, 1
      point(l) = (fmax/fmin) ** ((l-1.)/(N-1.))
! we work in angular frequency, not frequency
      point(l) = TWO_PI * point(l) * fmin
    enddo
  endif

  end subroutine remplit_point

  subroutine classical_linear_least_squares(Qref,poids,point,N,fmin,fmax)

  use constants, only: TWO_PI

  implicit none

  integer, intent(in) :: N
  double precision, intent(in) :: Qref,fmin,fmax
  double precision, dimension(1:N), intent(out) :: point,poids

  integer k,m
  logical, dimension(1:N) :: ia
  double precision ref,freq,chi2
  double precision, dimension(1:N,1:N) :: covar
  double precision, dimension(1:2*N-1) :: x,y_ref,sig

  m = 2*N-1

  call remplit_point(fmin,fmax,N,point)

  ref = 1.0 / Qref

  do k = 1,m
    freq = (fmax/fmin) ** ((k - 1.)/(m - 1.))
! we work in angular frequency, not frequency
    freq = TWO_PI * fmin * freq
    x(k) = freq
    y_ref(k) = ref
    sig(k) = 1.
  enddo

  do k = 1,N
    ia(k) = .true.
  enddo

  call lfit_zener(x,y_ref,sig,m,poids,ia,covar,chi2,N,Qref,point)

  end subroutine classical_linear_least_squares

! Calcul des coefficients par optimization non-lineaire avec contraintes

  subroutine solvopt(n,x,f,fun,flg,grad,options,flfc,func,flgc,gradc,Qref,Kopt,theta_min,theta_max,f_min,f_max)
!-----------------------------------------------------------------------------
! The subroutine SOLVOPT performs a modified version of Shor's r-algorithm in
! order to find a local minimum resp. maximum of a nonlinear function
! defined on the n-dimensional Euclidean space
! or
! a local minimum for a nonlinear constrained problem:
! min { f(x): g(x) ( < )= 0, g(x) in R(m), x in R(n) }.
! Arguments:
! n       is the space dimension (integer*4),
! x       is the n-vector, the coordinates of the starting point
!         at a call to the subroutine and the optimizer at regular return
!         (double precision),
! f       returns the optimum function value
!         (double precision),
! fun     is the entry name of a subroutine which computes the value
!         of the function fun at a point x, should be declared as external
!         in a calling routine,
!        synopsis: fun(x,f)
! grad    is the entry name of a subroutine which computes the gradient
!         vector of the function fun at a point x, should be declared as
!         external in a calling routine,
!         synopsis: grad(x,g)
! func    is the entry name of a subroutine which computes the MAXIMAL
!         RESIDIAL!!! (a scalar) for a set of constraints at a point x,
!         should be declared as external in a calling routine,
!         synopsis: func(x,fc)
! gradc   is the entry name of a subroutine which computes the gradient
!         vector for a constraint with the MAXIMAL RESIDUAL at a point x,
!         should be declared as external in a calling routine,
!        synopsis: gradc(x,gc)
! flg,    (logical) is a flag for the use of a subroutine grad:
!         .true. means gradients are calculated by the user-supplied routine.
! flfc,   (logical) is a flag for a constrained problem:
!         .true. means the maximal residual for a set of constraints
!         is calculated by func.
! flgc,   (logical) is a flag for the use of a subroutine gradc:
!         .true. means gradients of the constraints are calculated
!         by the user-supplied routine.
! options is a vector of optional parameters (double precision):
!     options(1)= H, where sign(H)=-1 resp. sign(H)=+1 means minimize resp.
!         maximize fun (valid only for an unconstrained problem) and
!         H itself is a factor for the initial trial step size
!         (options(1)=-1.d0 by default),
!     options(2)= relative error for the argument in terms of the infinity-norm
!         (1.d-4 by default),
!     options(3)= relative error for the function value (1.d-6 by default),
!     options(4)= limit for the number of iterations (1.5d4 by default),
!     options(5)= control of the display of intermediate results and error
!         resp. warning messages (default value is 0.d0, i.e., no intermediate
!         output but error and warning messages, see the manual for more),
!     options(6)= maximal admissible residual for a set of constraints
!         (options(6)=1.d-8 by default, see the manual for more),
!    *options(7)= the coefficient of space dilation (2.5d0 by default),
!    *options(8)= lower bound for the stepsize used for the difference
!        approximation of gradients (1.d-11 by default,see the manual for more).
!   (* ... changes should be done with care)
! returned optional values:
!     options(9),  the number of iterations, if positive,
!         or an abnormal stop code, if negative (see manual for more),
!                -1: allocation error,
!                -2: improper space dimension,
!                -3: fun returns an improper value,
!                -4: grad returns a zero vector or improper value at the
!                    starting point,
!                -5: func returns an improper value,
!                -6: gradc returns an improper value,
!                -7: function is unbounded,
!                -8: gradient is zero at the point,
!                    but stopping criteria are not fulfilled,
!                -9: iterations limit exceeded,
!               -11: Premature stop is possible,
!               -12: Result may not provide the true optimum,
!               -13: function is flat: result may be inaccurate
!                   in view of a point.
!               -14: function is steep: result may be inaccurate
!                    in view of a function value,
!       options(10), the number of objective function evaluations, and
!       options(11), the number of gradient evaluations.
!       options(12), the number of constraint function evaluations, and
!       options(13), the number of constraint gradient evaluations.
! ____________________________________________________________________________
!
      implicit none
      !include 'messages.inc'

      integer, intent(in) :: Kopt
      double precision, intent(in) :: Qref,theta_min,theta_max,f_min,f_max

      logical flg,flgc,flfc, constr, app, appconstr
      logical FsbPnt, FsbPnt1, termflag, stopf
      logical stopping, dispwarn, Reset, ksm,knan,obj
      integer n, kstore, ajp,ajpp,knorms, k, kcheck, numelem
      integer dispdata, ld, mxtc, termx, limxterm, nzero, krerun
      integer warnno, kflat, stepvanish, i,j,ni,ii, kd,kj,kc,ip
      integer iterlimit, kg,k1,k2, kless,   allocerr
      double precision options(13),doptions(13)
      double precision x(n),f
      double precision nsteps(3), gnorms(10), kk, nx
      double precision ajb,ajs, des, dq,du20,du10,du03
      double precision n_float, cnteps
      double precision low_bound, ZeroGrad, ddx, y
      double precision lowxbound, lowfbound, detfr, detxr, grbnd
      double precision fp,fp1,fc,f1,f2,fm,fopt,frec,fst, fp_rate
      double precision PenCoef, PenCoefNew
      double precision gamma,w,wdef,h1,h,hp
      double precision dx,ng,ngc,nng,ngt,nrmz,ng1,d,dd, laststep
      double precision zero,one,two,three,four,five,six,seven
      double precision eight,nine,ten,hundr
      double precision infty, epsnorm,epsnorm2,powerm12
      double precision, dimension(:,:), allocatable :: B
      double precision, dimension(:), allocatable :: g
      double precision, dimension(:), allocatable :: g0
      double precision, dimension(:), allocatable :: g1
      double precision, dimension(:), allocatable :: gt
      double precision, dimension(:), allocatable :: gc
      double precision, dimension(:), allocatable :: z
      double precision, dimension(:), allocatable :: x1
      double precision, dimension(:), allocatable :: xopt
      double precision, dimension(:), allocatable :: xrec
      double precision, dimension(:), allocatable :: grec
      double precision, dimension(:), allocatable :: xx
      double precision, dimension(:), allocatable :: deltax
      integer, dimension(:), allocatable :: idx
      character(len=100) :: endwarn
      character(len=19) :: allocerrstr
      external fun,grad,func,gradc

      data zero/0.d0/, one/1.d0/, two/2.d0/, three/3.d0/, four/4.d0/, &
         five/5.d0/, six/6.d0/, seven/7.d0/, eight/8.d0/, nine/9.d0/, &
         ten/1.d1/,  hundr/1.d2/, powerm12/1.d-12/, &
         infty /1.d100/, epsnorm /1.d-15/,  epsnorm2 /1.d-30/, &
         allocerrstr/'Allocation Error = '/
! Check the dimension:
      if (n < 2) then
          print *, 'SolvOpt error:'
          print *, 'Improper space dimension.'
         call stop_the_code('error in allocate statement in SolvOpt')
        options(9)=-one
        goto 999
      endif
      n_float=dble(n)
! allocate working arrays:
      allocate (B(n,n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (g(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (g0(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (g1(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (gt(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (gc(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (z(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (x1(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (xopt(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (xrec(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (grec(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (xx(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (deltax(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif
      allocate (idx(n),stat=allocerr)
      if (allocerr /= 0) then
         options(9)=-one
         print *,allocerrstr,allocerr
         call stop_the_code('error in allocate statement in SolvOpt')
      endif

! store flags:
      app= .not. flg
      constr=flfc
      appconstr= .not. flgc
! Default values for options:
      call soptions(doptions)
      do i =1,8
            if (options(i) == zero) then
               options(i)=doptions(i)
            else if (i == 2 .or. i == 3 .or. i == 6) then
               options(i)=dmax1(options(i),powerm12)
               options(i)=dmin1(options(i),one)
               if (i == 2)options(i)=dmax1(options(i),options(8)*hundr)
            else if (i == 7) then
               options(7)=dmax1(options(i),1.5d0)
            endif
      enddo

! WORKING CONSTANTS AND COUNTERS ----{

      options(10)=zero    !! counter for function calculations
      options(11)=zero    !! counter for gradient calculations
      options(12)=zero    !! counter for constraint function calculations
      options(13)=zero    !! counter for constraint gradient calculations
      iterlimit=idint(options(4))
      if (constr) then
        h1=-one           !! NLP: restricted to minimization
        cnteps=options(6)
      else
        h1=dsign(one,options(1))  !! Minimize resp. maximize a function
      endif
      k=0                         !! Iteration counter
      wdef=one/options(7)-one     !! Default space transf. coeff.

! Gamma control ---{
      ajb=one+1.d-1/n_float**2    !! Base I
      ajp=20
      ajpp=ajp                    !! Start value for the power
      ajs=1.15d0                  !! Base II
      knorms=0
      do i =1,10
       gnorms(i)=zero
      enddo
!---}
! Display control ---{
      if (options(5) <= zero) then
         dispdata=0
         if (options(5) == -one) then
            dispwarn=.false.
         else
            dispwarn=.true.
         endif
      else
         dispdata=idnint(options(5))
         dispwarn=.true.
      endif
      ld=dispdata
!---}

! Stepsize control ---{
      dq=5.1d0           !! Step divider (at f_{i+1} > gamma*f_{i})
      du20=two
      du10=1.5d0
      du03=1.05d0        !! Step multipliers (at certain steps made)
      kstore=3
      do i = 1,kstore
       nsteps(i)=zero    !! Steps made at the last 'kstore' iterations
      enddo
      if (app) then
        des=6.3d0        !! Desired number of steps per 1-D search
      else
        des=3.3d0
      endif
      mxtc=3             !! Number of trial cycles (steep wall detect)
!---}
      termx=0
      limxterm=50        !! Counter and limit for x-criterion
! stepsize for gradient approximation
      ddx=dmax1(1.d-11,options(8))

      low_bound=-one+1.d-4     !! Lower bound cosine used to detect a ravine
      ZeroGrad=n_float*1.d-16  !! Lower bound for a gradient norm
      nzero=0                  !! Zero-gradient events counter
! Low bound for the values of variables to take into account
      lowxbound=dmax1(options(2),1.d-3)
! Lower bound for function values to be considered as making difference
      lowfbound=options(3)**2
      krerun=0                 !! Re-run events counter
      detfr=options(3)*hundr   !! Relative error for f/f_{record}
      detxr=options(2)*ten     !! Relative error for norm(x)/norm(x_{record})
      warnno=0                 !! the number of warn.mess. to end with
      kflat=0                  !! counter for points of flatness
      stepvanish=0             !! counter for vanished steps
      stopf=.false.
! ----}  End of setting constants
! ----}  End of the preamble
!--------------------------------------------------------------------
! Compute the function  ( first time ) ----{
      call fun(x,f,Qref,n/2,n,Kopt,f_min,f_max)
      options(10)=options(10)+one
      if (dabs(f) >= infty) then
         if (dispwarn) then
            print *,'SolvOpt error: function equals infinity at the point.'
            print *,'Choose another starting point.'
         endif
         options(9)=-three
         goto 999
      endif
      do i = 1,n
        xrec(i)=x(i)
      enddo
      frec=f     !! record point and function value
! Constrained problem
      if (constr) then
          kless=0
          fp=f
          call func(x,fc,n/2,n,theta_min,theta_max)
          options(12)=options(12)+one
          if (dabs(fc) >= infty) then
             if (dispwarn) then
                print *,'SolvOpt error: FUNC returns infinite value at the point.'
                print *,'Choose another starting point.'
             endif
             options(9)=-five
             goto 999
          endif
        PenCoef=one          !! first rough approximation
        if (fc <= cnteps) then
         FsbPnt=.true.       !! feasible point
         fc=zero
        else
         FsbPnt=.false.
        endif
        f=f+PenCoef*fc
      endif
! ----}
! COMPUTE THE GRADIENT ( FIRST TIME ) ----{
      if (app) then
        do i = 1,n
         deltax(i)=h1*ddx
        enddo
        obj=.true.
        !if (constr) then
           !call apprgrdn()
        !else
           !call apprgrdn()
        !endif
        options(10)=options(10)+n_float
      else
        call grad(x,g,Qref,n/2,n,Kopt,f_min,f_max)
        options(11)=options(11)+one
      endif
      ng=zero
      do i = 1,n
         ng=ng+g(i)*g(i)
      enddo
      ng=dsqrt(ng)
      if (ng >= infty) then
         if (dispwarn) then
            print *,'SolvOpt error:'
            print *,'Gradient equals infinity at the starting point.'
            print *,'Choose another starting point.'
         endif
         options(9)=-four
         goto 999
      else if (ng < ZeroGrad) then
         if (dispwarn) then
            print *,'SolvOpt error:'
            print *,'Gradient equals zero at the starting point.'
            print *,'Choose another starting point.'
         endif
         options(9)=-four
         goto 999
      endif
      if (constr) then
       if (.not. FsbPnt) then
         !if (appconstr) then
            !do j = 1,n
              !if (x(j) >= zero) then
                 !deltax(j)=ddx
              !else
                 !deltax(j)=-ddx
              !endif
            !enddo
            !obj=.false.
            !call apprgrdn()
         if (.not. appconstr) then
            call gradc(x,gc,n/2,n,theta_min,theta_max)
         endif
         ngc=zero
         do i = 1,n
           ngc=ngc+gc(i)*gc(i)
         enddo
         ngc=dsqrt(ngc)
         if (ng >= infty) then
            if (dispwarn) then
               print *,'SolvOpt error: GRADC returns infinite vector at the point.'
               print *,'Choose another starting point.'
            endif
            options(9)=-six
            goto 999
         else if (ng < ZeroGrad) then
            if (dispwarn) then
               print *,'SolvOpt error: GRADC returns zero vector at an infeasible point.'
            endif
            options(9)=-six
            goto 999
         endif
         do i = 1,n
           g(i)=g(i)+PenCoef*gc(i)
         enddo
         ng=zero
         do i = 1,n
           ng=ng+g(i)*g(i)
           grec(i)=g(i)
         enddo
         ng=dsqrt(ng)
       endif
      endif
      do i = 1,n
        grec(i)=g(i)
      enddo
      nng=ng
! ----}
! INITIAL STEP SIZE
      d=zero
      do i = 1,n
        if (d < dabs(x(i))) d=dabs(x(i))
      enddo
      h=h1*dsqrt(options(2))*d                  !! smallest possible stepsize
      if (dabs(options(1)) /= one) then
        h=h1*dmax1(dabs(options(1)),dabs(h))    !! user-supplied stepsize
      else
          h=h1*dmax1(one/dlog(ng+1.1d0),dabs(h)) !! calculated stepsize
      endif

! RESETTING LOOP ----{
      do while (.true.)
        kcheck=0                       !! Set checkpoint counter.
        kg=0                           !! stepsizes stored
        kj=0                           !! ravine jump counter
        do i = 1,n
          do j = 1,n
            B(i,j)=zero
          enddo
          B(i,i)=one                   !! re-set transf. matrix to identity
          g1(i)=g(i)
        enddo
        fst=f
        dx=0
! ----}

! MAIN ITERATIONS ----{

        do while (.true.)
          k=k+1
          kcheck=kcheck+1
          laststep=dx
! ADJUST GAMMA --{
           gamma=one+dmax1(ajb**((ajp-kcheck)*n),two*options(3))
           gamma=dmin1 ( gamma,ajs**dmax1(one,dlog10(nng+one)) )
! --}
       ngt=zero
       ng1=zero
       dd=zero
       do i = 1,n
         d=zero
         do j = 1,n
            d=d+B(j,i)*g(j)
         enddo
         gt(i)=d
         dd=dd+d*g1(i)
         ngt=ngt+d*d
         ng1=ng1+g1(i)*g1(i)
       enddo
       ngt=dsqrt(ngt)
       ng1=dsqrt(ng1)
       dd=dd/ngt/ng1

       w=wdef
! JUMPING OVER A RAVINE ----{
       if (dd < low_bound) then
        if (kj == 2) then
          do i = 1,n
           xx(i)=x(i)
          enddo
        endif
        if (kj == 0) kd=4
        kj=kj+1
        w=-.9d0              !! use large coef. of space dilation
        h=h*two
        if (kj > 2*kd) then
          kd=kd+1
          warnno=1
          endwarn='Premature stopping is possible. Try to re-run the routine from the obtained point.'
          do i = 1,n
            if (dabs(x(i)-xx(i)) < epsnorm*dabs(x(i))) then
             if (dispwarn) then
                print *,'SolvOpt warning:'
                print *,'Ravine with a flat bottom is detected.'
             endif
            endif
          enddo
        endif
       else
        kj=0
       endif
! ----}
! DILATION ----{
       nrmz=zero
       do i = 1,n
         z(i)=gt(i)-g1(i)
         nrmz=nrmz+z(i)*z(i)
       enddo
       nrmz=dsqrt(nrmz)
       if (nrmz > epsnorm*ngt) then
        do i = 1,n
         z(i)=z(i)/nrmz
        enddo
! New direction in the transformed space: g1=gt+w*(z*gt')*z and
! new inverse matrix: B = B ( I + (1/alpha -1)zz' )
        d = zero
        do i = 1,n
          d=d+z(i)*gt(i)
        enddo
        ng1=zero
        d = d*w
        do i = 1,n
          dd=zero
          g1(i)=gt(i)+d*z(i)
          ng1=ng1+g1(i)*g1(i)
          do j = 1,n
             dd=dd+B(i,j)*z(j)
          enddo
          dd=w*dd
          do j = 1,n
            B(i,j)=B(i,j)+dd*z(j)
          enddo
        enddo
        ng1=dsqrt(ng1)
       else
        do i = 1,n
         z(i)=zero
         g1(i)=gt(i)
        enddo
        nrmz=zero
       endif
       do i = 1,n
           gt(i)=g1(i)/ng1
       enddo
        do i = 1,n
          d=zero
            do j = 1,n
               d=d+B(i,j)*gt(j)
            enddo
          g0(i)=d
        enddo
! ----}
! RESETTING ----{
        if (kcheck > 1) then
           numelem=0
           do i = 1,n
              if (dabs(g(i)) > ZeroGrad) then
                 numelem=numelem+1
                 idx(numelem)=i
              endif
           enddo
           if (numelem > 0) then
              grbnd=epsnorm*dble(numelem**2)
              ii=0
              do i = 1,numelem
                 j=idx(i)
                 if (dabs(g1(j)) <= dabs(g(j))*grbnd) ii=ii+1
              enddo
              if (ii == n .or. nrmz == zero) then
                if (dispwarn) then
                  print *,'SolvOpt warning: Normal re-setting of a transformation matrix'
                endif
                if (dabs(fst-f) < dabs(f)*1.d-2) then
                   ajp=ajp-10*n
                else
                   ajp=ajpp
                endif
                h=h1*dx/three
                k=k-1
                exit
              endif
           endif
        endif
! ----}
! STORE THE CURRENT VALUES AND SET THE COUNTERS FOR 1-D SEARCH
        do i = 1,n
         xopt(i)=x(i)
        enddo
        fopt=f
        k1=0
        k2=0
        ksm=.false.
        kc=0
        knan=.false.
        hp=h
        if (constr) Reset=.false.
! 1-D SEARCH ----{
        do while (.true.)
         do i = 1,n
          x1(i)=x(i)
         enddo
         f1=f
         if (constr) then
           FsbPnt1=FsbPnt
           fp1=fp
         endif
! NEW POINT
         do i = 1,n
            x(i)=x(i)+hp*g0(i)
         enddo
           ii=0
           do i = 1,n
            if (dabs(x(i)-x1(i)) < dabs(x(i))*epsnorm) ii=ii+1
           enddo
! function value
         call fun(x,f,Qref,n/2,n,Kopt,f_min,f_max)
         options(10)=options(10)+one
         if (h1*f >= infty) then
            if (dispwarn) then
              print *,'SolvOpt error: function is unbounded.'
            endif
            options(9)=-seven
            goto 999
         endif
         if (constr) then
           fp=f
           call func(x,fc,n/2,n,theta_min,theta_max)
           options(12)=options(12)+one
           if (dabs(fc) >= infty) then
               if (dispwarn) then
                  print *,'SolvOpt error: FUNC returns infinite value at the point.'
                  print *,'Choose another starting point.'
               endif
               options(9)=-five
               goto 999
           endif
           if (fc <= cnteps) then
              FsbPnt=.true.
              fc=zero
           else
              FsbPnt=.false.
              fp_rate=fp-fp1
              if (fp_rate < -epsnorm) then
               if (.not. FsbPnt1) then
                d=zero
                do i = 1,n
                  d=d+(x(i)-x1(i))**2
                enddo
                d=dsqrt(d)
                PenCoefNew=-1.5d1*fp_rate/d
                if (PenCoefNew > 1.2d0*PenCoef) then
                  PenCoef=PenCoefNew
                  Reset=.true.
                  kless=0
                  f=f+PenCoef*fc
                  exit
                endif
               endif
              endif
           endif
           f=f+PenCoef*fc
         endif
         if (dabs(f) >= infty) then
             if (dispwarn) then
               print *,'SolvOpt warning: function equals infinity at the point.'
             endif
             if (ksm .or. kc >= mxtc) then
                options(9)=-three
                goto 999
             else
                k2=k2+1
                k1=0
                hp=hp/dq
                do i = 1,n
                 x(i)=x1(i)
                enddo
                f=f1
                knan=.true.
                if (constr) then
                  FsbPnt=FsbPnt1
                  fp=fp1
                endif
             endif
! STEP SIZE IS ZERO TO THE EXTENT OF EPSNORM
         else if (ii == n) then
                stepvanish=stepvanish+1
                if (stepvanish >= 5) then
                    options(9)=-ten-four
                    if (dispwarn) then
                       print *,'SolvOpt: Termination warning:'
                       print *,'Stopping criteria are not fulfilled. The function is very steep at the solution.'
                    endif
                    goto 999
                else
                    do i = 1,n
                     x(i)=x1(i)
                    enddo
                    f=f1
                    hp=hp*ten
                    ksm=.true.
                    if (constr) then
                       FsbPnt=FsbPnt1
                       fp=fp1
                    endif
                endif
! USE SMALLER STEP
         else if (h1*f < h1*gamma**idint(dsign(one,f1))*f1) then
             if (ksm) exit
             k2=k2+1
             k1=0
             hp=hp/dq
             do i = 1,n
              x(i)=x1(i)
             enddo
             f=f1
             if (constr) then
                FsbPnt=FsbPnt1
                fp=fp1
             endif
             if (kc >= mxtc) exit
! 1-D OPTIMIZER IS LEFT BEHIND
         else
             if (h1*f <= h1*f1) exit
! USE LARGER STEP
             k1=k1+1
             if (k2 > 0) kc=kc+1
             k2=0
             if (k1 >= 20) then
                 hp=du20*hp
             else if (k1 >= 10) then
                 hp=du10*hp
             else if (k1 >= 3) then
                 hp=du03*hp
             endif
         endif
        enddo
! ----}  End of 1-D search
! ADJUST THE TRIAL STEP SIZE ----{
        dx=zero
        do i = 1,n
           dx=dx+(xopt(i)-x(i))**2
        enddo
        dx=dsqrt(dx)
        if (kg < kstore)  kg=kg+1
        if (kg >= 2) then
           do i =kg,2,-1
             nsteps(i)=nsteps(i-1)
           enddo
        endif
        d=zero
        do i = 1,n
           d=d+g0(i)*g0(i)
        enddo
        d=dsqrt(d)
        nsteps(1)=dx/(dabs(h)*d)
        kk=zero
        d=zero
        do i = 1,kg
           dd=dble(kg-i+1)
           d=d+dd
           kk=kk+nsteps(i)*dd
        enddo
        kk=kk/d
        if (kk > des) then
             if (kg == 1) then
                h=h*(kk-des+one)
             else
                h=h*dsqrt(kk-des+one)
             endif
        else if (kk < des) then
             h=h*dsqrt(kk/des)
        endif

        if (ksm) stepvanish=stepvanish+1
! ----}
! COMPUTE THE GRADIENT ----{
        if (app) then
          do j = 1,n
            if (g0(j) >= zero) then
               deltax(j)=h1*ddx
            else
               deltax(j)=-h1*ddx
            endif
          enddo
          obj=.true.
          !if (constr) then
             !call apprgrdn()
          !else
             !call apprgrdn()
          !endif
          !options(10)=options(10)+n_float
        else
          call grad(x,g,Qref,n/2,n,Kopt,f_min,f_max)
          options(11)=options(11)+one
        endif
        ng=zero
        do i = 1,n
          ng=ng+g(i)*g(i)
        enddo
        ng=dsqrt(ng)
        if (ng >= infty) then
         if (dispwarn) then
           print *,'SolvOpt error:'
           print *,'Gradient equals infinity at the starting point.'
         endif
         options(9)=-four
         goto 999
        else if (ng < ZeroGrad) then
         if (dispwarn) then
           print *,'SolvOpt warning:'
           print *,'Gradient is zero, but stopping criteria are not fulfilled.'
         endif
         ng=ZeroGrad
        endif
! Constraints:
        if (constr) then
         if (.not. FsbPnt) then
           if (ng < 1.d-2*PenCoef) then
              kless=kless+1
              if (kless >= 20) then
                 PenCoef=PenCoef/ten
                 Reset=.true.
                 kless=0
              endif
           else
              kless=0
           endif
           !if (appconstr) then
                 !do j = 1,n
                   !if (x(j) >= zero) then
                      !deltax(j)=ddx
                   !else
                      !deltax(j)=-ddx
                   !endif
                 !enddo
                 !obj=.false.
                 !call apprgrdn()
                 !options(12)=options(12)+n_float
           if (.not. appconstr) then
                 call gradc(x,gc,n/2,n,theta_min,theta_max)
                 options(13)=options(13)+one
           endif
           ngc=zero
           do i = 1,n
              ngc=ngc+gc(i)*gc(i)
           enddo
           ngc=dsqrt(ngc)
           if (ngc >= infty) then
                  if (dispwarn) then
                     print *,'SolvOpt error: GRADC returns infinite vector at the point.'
                  endif
                  options(9)=-six
                  goto 999
           else if (ngc < ZeroGrad .and. .not. appconstr) then
                  if (dispwarn) then
                     print *,'SolvOpt error: GRADC returns zero vector at an infeasible point.'
                  endif
                  options(9)=-six
                  goto 999
           endif
           do i = 1,n
             g(i)=g(i)+PenCoef*gc(i)
           enddo
           ng=zero
           do i = 1,n
              ng=ng+g(i)*g(i)
           enddo
           ng=dsqrt(ng)
           if (Reset) then
              if (dispwarn) then
                 print *,'SolvOpt warning:'
                 print *,'Re-setting due to the use of a new penalty coefficient.'
              endif
              h=h1*dx/three
              k=k-1
              nng=ng
              exit
           endif
         endif
        endif
        if (h1*f > h1*frec) then
          frec=f
          do i = 1,n
            xrec(i)=x(i)
            grec(i)=g(i)
          enddo
        endif
! ----}
       if (ng > ZeroGrad) then
        if (knorms < 10)  knorms=knorms+1
        if (knorms >= 2) then
          do i =knorms,2,-1
           gnorms(i)=gnorms(i-1)
          enddo
        endif
        gnorms(1)=ng
        nng=one
          do i = 1,knorms
            nng=nng*gnorms(i)
          enddo
        nng=nng**(one/dble(knorms))
       endif
! Norm X:
       nx=zero
       do i = 1,n
        nx=nx+x(i)*x(i)
       enddo
       nx=dsqrt(nx)

! DISPLAY THE CURRENT VALUES ----{
       if (k == ld) then
         print *, &
             'Iteration # ..... function value ..... ', &
             'Step Value ..... Gradient Norm'
         print '(5x,i5,7x,g13.5,6x,g13.5,7x,g13.5)', k,f,dx,ng
         ld=k+dispdata
       endif
!----}
! CHECK THE STOPPING CRITERIA ----{
      termflag=.true.
      if (constr) then
        if (.not. FsbPnt) termflag=.false.
      endif
      if (kcheck <= 5 .or. kcheck <= 12 .and. ng > one)termflag=.false.
      if (kc >= mxtc .or. knan)termflag=.false.
! ARGUMENT
       if (termflag) then
           ii=0
           stopping=.true.
           do i = 1,n
             if (dabs(x(i)) >= lowxbound) then
                ii=ii+1
                idx(ii)=i
                if (dabs(xopt(i)-x(i)) > options(2)*dabs(x(i))) then
                  stopping=.false.
                endif
             endif
           enddo
           if (ii == 0 .or. stopping) then
                stopping=.true.
                termx=termx+1
                d=zero
                do i = 1,n
                  d=d+(x(i)-xrec(i))**2
                enddo
                d=dsqrt(d)
! function
                if (dabs(f-frec) > detfr*dabs(f) .and. &
                  dabs(f-fopt) <= options(3)*dabs(f) .and. &
                  krerun <= 3 .and. .not. constr) then
                   stopping=.false.
                   if (ii > 0) then
                    do i = 1,ii
                     j=idx(i)
                     if (dabs(xrec(j)-x(j)) > detxr*dabs(x(j))) then
                       stopping=.true.
                       exit
                     endif
                    enddo
                   endif
                   if (stopping) then
                      if (dispwarn) then
                        print *,'SolvOpt warning:'
                        print *,'Re-run from recorded point.'
                      endif
                      ng=zero
                      do i = 1,n
                       x(i)=xrec(i)
                       g(i)=grec(i)
                       ng=ng+g(i)*g(i)
                      enddo
                      ng=dsqrt(ng)
                      f=frec
                      krerun=krerun+1
                      h=h1*dmax1(dx,detxr*nx)/dble(krerun)
                      warnno=2
                      endwarn='Result may not provide the optimum. The function apparently has many extremum points.'
                      exit
                   else
                      h=h*ten
                   endif
                else if (dabs(f-frec) > options(3)*dabs(f) .and. &
                  d < options(2)*nx .and. constr) then
                   continue
                else if (dabs(f-fopt) <= options(3)*dabs(f) .or. &
                   dabs(f) <= lowfbound .or. &
                   (dabs(f-fopt) <= options(3) .and. &
                    termx >= limxterm )) then
                  if (stopf) then
                   if (dx <= laststep) then
                    if (warnno == 1 .and. ng < dsqrt(options(3))) then
                       warnno=0
                    endif
                    if (.not. app) then
                      do i = 1,n
                       if (dabs(g(i)) <= epsnorm2) then
                         warnno=3
                         endwarn='Result may be inaccurate in the coordinates. The function is flat at the solution.'
                         exit
                       endif
                      enddo
                    endif
                    if (warnno /= 0) then
                       options(9)=-dble(warnno)-ten
                       if (dispwarn) then
                         print *,'SolvOpt: Termination warning:'
                         print *,endwarn
                         if (app) print *,'The above warning may be reasoned by inaccurate gradient approximation'
                       endif
                    else
                       options(9)=dble(k)
!! DK DK               if (dispwarn) print *,'SolvOpt: Normal termination.'
                    endif
                    goto 999
                   endif
                  else
                   stopf=.true.
                  endif
                else if (dx < powerm12*dmax1(nx,one) .and. &
                       termx >= limxterm) then
                     options(9)=-four-ten
                     if (dispwarn) then
                       print *,'SolvOpt: Termination warning:'
                       print *,'Stopping criteria are not fulfilled. The function is very steep at the solution.'
                       if (app) print *,'The above warning may be reasoned by inaccurate gradient approximation'
                       f=frec
                       do i = 1,n
                        x(i)=xrec(i)
                       enddo
                     endif
                     goto 999
                endif
           endif
       endif
! ITERATIONS LIMIT
            if (k == iterlimit) then
                options(9)=-nine
                if (dispwarn) then
                  print *,'SolvOpt warning:'
                  print *,'Iterations limit exceeded.'
                endif
                goto 999
            endif
! ----}
! ZERO GRADIENT ----{
          if (constr) then
            if (ng <= ZeroGrad) then
                if (dispwarn) then
                  print *,'SolvOpt: Termination warning:'
                  print *,'Gradient is zero, but stopping criteria are not fulfilled.'
                endif
                options(9)=-eight
                goto 999
            endif
          else
            if (ng <= ZeroGrad) then
             nzero=nzero+1
             if (dispwarn) then
               print *,'SolvOpt warning:'
               print *,'Gradient is zero, but stopping criteria are not fulfilled.'
             endif
             if (nzero >= 3) then
               options(9)=-eight
               goto 999
             endif
             do i = 1,n
               g0(i)=-h*g0(i)/two
             enddo
             do i =1,10
               do j = 1,n
                x(j)=x(j)+g0(j)
               enddo
               call fun(x,f,Qref,n/2,n,Kopt,f_min,f_max)
               options(10)=options(10)+one
               if (dabs(f) >= infty) then
                 if (dispwarn) then
                   print *,'SolvOpt error: function equals infinity at the point.'
                 endif
                 options(9)=-three
                 goto 999
               endif
               !if (app) then
                   !do j = 1,n
                     !if (g0(j) >= zero) then
                        !deltax(j)=h1*ddx
                     !else
                        !deltax(j)=-h1*ddx
                     !endif
                   !enddo
                   !obj=.true.
                   !call apprgrdn()
                   !options(10)=options(10)+n_float
               if (.not. app) then
                   call grad(x,g,Qref,n/2,n,Kopt,f_min,f_max)
                   options(11)=options(11)+one
               endif
               ng=zero
               do j = 1,n
                  ng=ng+g(j)*g(j)
               enddo
               ng=dsqrt(ng)
               if (ng >= infty) then
                    if (dispwarn) then
                      print *,'SolvOpt error:'
                      print *,'Gradient equals infinity at the starting point.'
                    endif
                    options(9)=-four
                    goto 999
               endif
               if (ng > ZeroGrad) exit
             enddo
             if (ng <= ZeroGrad) then
                if (dispwarn) then
                  print *,'SolvOpt: Termination warning:'
                  print *,'Gradient is zero, but stopping criteria are not fulfilled.'
                endif
                options(9)=-eight
                goto 999
             endif
             h=h1*dx
             exit
            endif
          endif
! ----}
! function is flat at the point ----{
          if (.not. constr .and. &
             dabs(f-fopt) < dabs(fopt)*options(3) .and. kcheck > 5 .and. ng < one) then

           ni=0
           do i = 1,n
             if (dabs(g(i)) <= epsnorm2) then
               ni=ni+1
               idx(ni)=i
             endif
           enddo
           if (ni >= 1 .and. ni <= n/2 .and. kflat <= 3) then
             kflat=kflat+1
             if (dispwarn) then
                print *,'SolvOpt warning:'
                print *,'The function is flat in certain directions.'
             endif
             warnno=1
             endwarn='Premature stopping is possible. Try to re-run the routine from the obtained point.'
             do i = 1,n
               x1(i)=x(i)
             enddo
             fm=f
             do i = 1,ni
              j=idx(i)
              f2=fm
              y=x(j)
              if (y == zero) then
                x1(j)=one
              else if (dabs(y) < one) then
                x1(j)=dsign(one,y)
              else
                x1(j)=y
              endif
              do ip=1,20
               x1(j)=x1(j)/1.15d0
               call fun(x1,f1,Qref,n/2,n,Kopt,f_min,f_max)
               options(10)=options(10)+one
               if (dabs(f1) < infty) then
                 if (h1*f1 > h1*fm) then
                   y=x1(j)
                   fm=f1
                 else if (h1*f2 > h1*f1) then
                   exit
                 else if (f2 == f1) then
                   x1(j)=x1(j)/1.5d0
                 endif
                 f2=f1
               endif
              enddo
              x1(j)=y
             enddo
             if (h1*fm > h1*f) then
              !if (app) then
                !do j = 1,n
                  !deltax(j)=h1*ddx
                !enddo
                !obj=.true.
                !call apprgrdn()
                !options(10)=options(10)+n_float
              if (.not. app) then
                call grad(x1,gt,Qref,n/2,n,Kopt,f_min,f_max)
                options(11)=options(11)+one
              endif
              ngt=zero
              do i = 1,n
                ngt=ngt+gt(i)*gt(i)
              enddo
              if (ngt > epsnorm2 .and. ngt < infty) then
                if (dispwarn) print *,'Trying to recover by shifting insensitive variables.'
                do i = 1,n
                 x(i)=x1(i)
                 g(i)=gt(i)
                enddo
                ng=ngt
                f=fm
                h=h1*dx/three
                options(3)=options(3)/five
                exit
              endif   !! regular gradient
             endif   !! a better value has been found
           endif   !! function is flat
          endif   !! pre-conditions are fulfilled
! ----}
       enddo   !! iterations
      enddo   !! restart

999   continue

! deallocate working arrays:
      deallocate (idx,deltax,xx,grec,xrec,xopt,x1,z,gc,gt,g1,g0,g,B)

  end subroutine solvopt

  subroutine soptions(default)
! SOPTIONS returns the default values for the optional parameters
! used by SolvOpt.

  implicit none

  double precision default(13)

  default(1)  = -1.d0
  default(2)  = 1.d-4
  default(3)  = 1.d-6
  default(4)  = 15.d3
  default(5)  = 0.d0
  default(6)  = 1.d-8
  default(7)  = 2.5d0
  default(8)  = 1.d-12
  default(9)  = 0.d0
  default(10) = 0.d0
  default(11) = 0.d0
  default(12) = 0.d0
  default(13) = 0.d0

  end subroutine soptions

  subroutine func_objective(x,res,freq,Qref,N,Nopt)

  implicit none

  integer, intent(in) :: N,Nopt
  double precision, intent(in) :: freq,Qref
  double precision, intent(out) :: res
  double precision, dimension(1:Nopt), intent(in) :: x

  integer i
  double precision num,deno

  res = 0.d0
  do i = 1,N
    num = x(N+i)*x(N+i)*freq*Qref*(x(i)*x(i) - freq/qref)
    deno = (x(i) ** 4.) + freq*freq
    res = res + num/deno
  enddo

  end subroutine func_objective

  subroutine func_mini(x,res,Qref,N,Nopt,K,f_min,f_max)

! Nopt = 2*N : nombre de coefficients a optimiser

  use constants, only: TWO_PI

  implicit none

  integer, intent(in) :: N,Nopt,K
  double precision, intent(in) :: Qref,f_min,f_max
  double precision, intent(out) :: res
  double precision, dimension(1:Nopt), intent(in) :: x

  integer i
  double precision d,freq,aux

  res = 0.d0
  do i = 1,K
! we work in angular frequency, not frequency
    freq = TWO_PI * f_min*((f_max/f_min)**((i-1.d0)/(K-1.d0)))
    call func_objective(x,aux,freq,Qref,N,Nopt)
    d = aux - 1.d0
    res = res + d*d
  enddo

  end subroutine func_mini

  subroutine grad_func_mini(x,grad,Qref,N,Nopt,K,f_min,f_max)

  use constants, only: TWO_PI

  implicit none

  integer, intent(in) :: N,Nopt,K
  double precision, intent(in) :: Qref,f_min,f_max
  double precision, dimension(1:Nopt), intent(in) :: x
  double precision, dimension(1:Nopt), intent(out) :: grad

  integer i,l
  double precision R,temp0,temp1,temp2,temp3,tamp,aux1,aux2,aux3,aux4
  double precision, dimension(1:N) :: point,poids
  double precision, dimension(1:K) :: freq

  do i = 1,K
! we work in angular frequency, not frequency
    freq(i) = TWO_PI * f_min*((f_max/f_min)**((i-1.d0)/(K-1.d0)))
  enddo

  do l= 1,N
    point(l) = x(l)
    poids(l) = x(N+l)
  enddo

  do l= 1,N
    grad(l) = 0.d0
    grad(N+l) = 0.d0

    do i = 1,K
      call func_objective(x,R,freq(i),Qref,N,Nopt)
      temp3 = R - 1.d0
      temp0 = freq(i)*Qref

      ! derivee par rapport aux poids
      temp1 = temp0*(point(l)*point(l) - freq(i)/qref)
      temp1 = temp1*2.d0*poids(l)
      temp2 = (point(l)**4.d0) + freq(i)*freq(i)
      temp1 = temp1/temp2
      tamp = 2.d0*temp3*temp1
      grad(N+l) = grad(N+l) + tamp

      ! derivee par rapport aux points
      aux1 = -2.d0*(point(l)**5.d0) + 2.d0*point(l)*freq(i)*freq(i) + 4.d0*(point(l)**3.d0)*freq(i)/Qref
      aux3 = temp2*temp2
      aux4 = aux1/aux3
      aux4 = aux4*temp0
      aux2 = aux4*poids(l)*poids(l)
      tamp = 2.d0*temp3*aux2
      grad(l) = grad(l) + tamp
    enddo
  enddo

  end subroutine grad_func_mini

  subroutine max_residu(x,res,N,Nopt,theta_min,theta_max)

  implicit none

  integer, intent(in) :: N,Nopt
  double precision, intent(in) :: theta_min,theta_max
  double precision, intent(out) :: res
  double precision, dimension(1:Nopt), intent(in) :: x

  integer l
  double precision temp,aux

  temp = 0.d0
  res = 0.d0

  do l= 1,N
    aux = res
    temp = max(0.d0,x(l)*x(l)-(theta_max-theta_min))
    res = max(temp,aux)
  enddo

  end subroutine max_residu

  subroutine grad_max_residu(x,grad,N,Nopt,theta_min,theta_max)

  implicit none

  integer, intent(in) :: N,Nopt
  double precision, intent(in) :: theta_min,theta_max
  double precision, dimension(1:Nopt), intent(in) :: x
  double precision, dimension(1:Nopt), intent(out) :: grad

  integer l,l0
  double precision temp,res,aux,temp2
  double precision, dimension(1:N) :: point

  temp = 0.d0
  res = 0.d0

  do l= 1,N
    point(l) = x(l)
  enddo

  l0 = 1
  do l= 1,N
    aux = res
    temp = max(0.d0,point(l)*point(l) - (theta_max-theta_min))
    res = max(temp,aux)
    if (temp > aux) then
      l0 = l
    endif
  enddo

  do l= 1,N
    grad(N+l) = 0.d0
    if (l /= l0) then
      grad(l) = 0.d0
    else
      call max_residu(x,temp2,N,Nopt,theta_min,theta_max)
      if (temp2 == 0.d0) then
        grad(l0) = 0.d0
      else
        grad(l0) = 2.d0*point(l0)
      endif
    endif
  enddo

  end subroutine grad_max_residu

  subroutine nonlinear_optimization(N,Qref,f0,point,poids,f_min,f_max)

  use shared_input_parameters, only: USE_SOLVOPT

  implicit none

  integer, intent(in) :: N
  double precision, intent(in) :: Qref,f0,f_min,f_max
  double precision, dimension(1:N), intent(out) :: point,poids

  external func_mini,grad_func_mini,max_residu,grad_max_residu

  integer K,i
  logical flg,flfc,flgc
  double precision theta_min,theta_max,res
  double precision, dimension(1:2*N) :: x
  double precision, dimension(1:13) :: options

  flg = .true.
  flgc = .true.
  flfc = .true.

  K = 4*N
  theta_min = 0.d0        ! arbitrary lower limit from Bruno Lombard to make sure points never become negative
  theta_max = 1000.d0*f0  ! arbitrary upper limit from Bruno Lombard to make sure points never tend to infinity

  ! this is used as a first guess
  call classical_linear_least_squares(Qref,poids,point,N,f_min,f_max)
  if (.not. USE_SOLVOPT) return

  ! what follows is the nonlinear optimization part

  do i = 1,N
    x(i)   = sqrt(abs(point(i)) - theta_min)
    x(N+i) = sqrt(abs(poids(i)))
  enddo

  call soptions(options)
  call solvopt(2*N,x,res,func_mini,flg,grad_func_mini,options,flfc, &
      max_residu,flgc,grad_max_residu,Qref,K,theta_min,theta_max,f_min,f_max)

  do i = 1,N
    point(i) = theta_min + x(i)*x(i)
    poids(i) = x(N+i)*x(N+i)
  enddo

  end subroutine nonlinear_optimization

