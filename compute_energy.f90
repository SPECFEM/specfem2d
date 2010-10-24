
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
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

  subroutine compute_energy(displ_elastic,veloc_elastic,displs_poroelastic,velocs_poroelastic, &
         displw_poroelastic,velocw_poroelastic, &
         xix,xiz,gammax,gammaz,jacobian,ibool,elastic,poroelastic,hprime_xx,hprime_zz, &
         nspec,npoin,assign_external_model,it,deltat,t0,kmato,elastcoef,density, &
         porosity,tortuosity, &
         vpext,vsext,rhoext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,anisotropic,anisotropy,wxgll,wzgll,numat, &
         pressure_element,vector_field_element,e1,e11, &
         potential_dot_acoustic,potential_dot_dot_acoustic,TURN_ATTENUATION_ON,Mu_nu1,Mu_nu2,N_SLS)

! compute kinetic and potential energy in the solid (acoustic elements are excluded)

  implicit none

  include "constants.h"

  integer :: nspec,npoin,numat

! vector field in an element
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLX) :: vector_field_element

! pressure in an element
  integer :: N_SLS
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: pressure_element

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec,N_SLS) :: e1,e11
  double precision, dimension(NGLLX,NGLLZ,nspec) :: Mu_nu1,Mu_nu2

  real(kind=CUSTOM_REAL), dimension(npoin) :: potential_dot_acoustic,potential_dot_dot_acoustic

  logical :: TURN_ATTENUATION_ON

  integer :: it
  double precision :: t0,deltat

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

  logical, dimension(nspec) :: elastic,poroelastic,anisotropic

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz,jacobian

  integer, dimension(nspec) :: kmato

  logical :: assign_external_model

  double precision, dimension(2,numat) :: density
  double precision, dimension(numat) :: porosity,tortuosity
  double precision, dimension(6,numat) :: anisotropy
  double precision, dimension(4,3,numat) :: elastcoef
  double precision, dimension(NGLLX,NGLLZ,nspec) :: vpext,vsext,rhoext
  double precision, dimension(NGLLX,NGLLZ,nspec) ::  c11ext,c15ext,c13ext,c33ext,c35ext,c55ext

  real(kind=CUSTOM_REAL), dimension(NDIM,npoin) :: displ_elastic,veloc_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin) :: displs_poroelastic,velocs_poroelastic
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin) :: displw_poroelastic,velocw_poroelastic

! Gauss-Lobatto-Legendre points and weights
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  real(kind=CUSTOM_REAL), dimension(NGLLZ) :: wzgll

! array with derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

! local variables
  integer :: i,j,k,ispec

! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  real(kind=CUSTOM_REAL) :: dwx_dxl,dwz_dxl,dwx_dzl,dwz_dzl

! jacobian
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: kinetic_energy,potential_energy
  real(kind=CUSTOM_REAL) :: cpl,csl,rhol,mul_relaxed,lambdal_relaxed,lambdalplus2mul_relaxed,kappal
  real(kind=CUSTOM_REAL) :: mul_s,kappal_s,rhol_s
  real(kind=CUSTOM_REAL) :: kappal_f,rhol_f
  real(kind=CUSTOM_REAL) :: mul_fr,kappal_fr,phil,tortl
  real(kind=CUSTOM_REAL) :: D_biot,H_biot,C_biot,M_biot,rhol_bar
  real(kind=CUSTOM_REAL) :: mul_G,lambdal_G,lambdalplus2mul_G

  kinetic_energy = ZERO
  potential_energy = ZERO

! loop over spectral elements
  do ispec = 1,nspec

!---
!--- elastic spectral element
!---
    if(elastic(ispec)) then

! get relaxed elastic parameters of current spectral element
      lambdal_relaxed = elastcoef(1,1,kmato(ispec))
      mul_relaxed = elastcoef(2,1,kmato(ispec))
      lambdalplus2mul_relaxed = elastcoef(3,1,kmato(ispec))
      rhol  = density(1,kmato(ispec))

! double loop over GLL points
      do j = 1,NGLLZ
        do i = 1,NGLLX

!--- if external medium, get elastic parameters of current grid point
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
            mul_relaxed = rhol*csl*csl
            lambdal_relaxed = rhol*cpl*cpl - TWO*mul_relaxed
            lambdalplus2mul_relaxed = lambdal_relaxed + TWO*mul_relaxed
          endif

! derivative along x and along z
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
          jacobianl = jacobian(i,j,ispec)

! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

! compute kinetic energy
          kinetic_energy = kinetic_energy + &
              rhol*(veloc_elastic(1,ibool(i,j,ispec))**2 + &
                    veloc_elastic(2,ibool(i,j,ispec))**2) *wxgll(i)*wzgll(j)*jacobianl / TWO

! compute potential energy
          potential_energy = potential_energy + (lambdalplus2mul_relaxed*dux_dxl**2 &
              + lambdalplus2mul_relaxed*duz_dzl**2 &
              + two*lambdal_relaxed*dux_dxl*duz_dzl + mul_relaxed*(dux_dzl + duz_dxl)**2)*wxgll(i)*wzgll(j)*jacobianl / TWO

        enddo
      enddo

!---
!--- poroelastic spectral element
!---
    elseif(poroelastic(ispec)) then

! get relaxed elastic parameters of current spectral element
!for now replaced by solid, fluid, and frame parameters of current spectral element
    phil = porosity(kmato(ispec))
    tortl = tortuosity(kmato(ispec))
!solid properties
    mul_s = elastcoef(2,1,kmato(ispec))
    kappal_s = elastcoef(3,1,kmato(ispec)) - FOUR_THIRDS*mul_s
    rhol_s = density(1,kmato(ispec))
!fluid properties
    kappal_f = elastcoef(1,2,kmato(ispec))
    rhol_f = density(2,kmato(ispec))
!frame properties
    mul_fr = elastcoef(2,3,kmato(ispec))
    kappal_fr = elastcoef(3,3,kmato(ispec)) - FOUR_THIRDS*mul_fr
    rhol_bar =  (1.d0 - phil)*rhol_s + phil*rhol_f
!Biot coefficients for the input phi
      D_biot = kappal_s*(1.d0 + phil*(kappal_s/kappal_f - 1.d0))
      H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) + kappal_fr + FOUR_THIRDS*mul_fr
      C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
      M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)
!The RHS has the form : div T -phi/c div T_f + phi/ceta_fk^-1.partial t w
!where T = G:grad u_s + C div w I
!and T_f = C div u_s I + M div w I
!we are expressing lambdaplus2mu, lambda, and mu for G, C, and M
      mul_G = mul_fr
      lambdal_G = H_biot - TWO*mul_fr
      lambdalplus2mul_G = lambdal_G + TWO*mul_G

! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX

! derivative along x and along z
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          dwx_dxi = ZERO
          dwz_dxi = ZERO

          dwx_dgamma = ZERO
          dwz_dgamma = ZERO

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)


            dwx_dxi = dwx_dxi + displw_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            dwz_dxi = dwz_dxi + displw_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dwx_dgamma = dwx_dgamma + displw_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            dwz_dgamma = dwz_dgamma + displw_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
          jacobianl = jacobian(i,j,ispec)

! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          dwx_dxl = dwx_dxi*xixl + dwx_dgamma*gammaxl
          dwx_dzl = dwx_dxi*xizl + dwx_dgamma*gammazl

          dwz_dxl = dwz_dxi*xixl + dwz_dgamma*gammaxl
          dwz_dzl = dwz_dxi*xizl + dwz_dgamma*gammazl

! compute potential energy
          potential_energy = potential_energy + ( lambdalplus2mul_G*dux_dxl**2 &
              + lambdalplus2mul_G*duz_dzl**2 &
              + two*lambdal_G*dux_dxl*duz_dzl + mul_G*(dux_dzl + duz_dxl)**2 &
              + two*C_biot*dwx_dxl*dux_dxl + two*C_biot*dwz_dzl*duz_dzl &
              + two*C_biot*(dwx_dxl*duz_dzl + dwz_dzl*dux_dxl) &
              + M_biot*dwx_dxl**2 + M_biot*dwz_dzl**2 &
              + two*M_biot*dwx_dxl*dwz_dzl )*wxgll(i)*wzgll(j)*jacobianl / TWO

! compute kinetic energy
         if(phil > 0.0d0) then
          kinetic_energy = kinetic_energy + ( &
              rhol_bar*(velocs_poroelastic(1,ibool(i,j,ispec))**2 + velocs_poroelastic(2,ibool(i,j,ispec))**2) &
              + rhol_f*tortl/phil*(velocw_poroelastic(1,ibool(i,j,ispec))**2 + velocw_poroelastic(2,ibool(i,j,ispec))**2) &
              + rhol_f*(velocs_poroelastic(1,ibool(i,j,ispec))*velocw_poroelastic(1,ibool(i,j,ispec)) &
              + velocs_poroelastic(2,ibool(i,j,ispec))*velocw_poroelastic(2,ibool(i,j,ispec))) &
                 )*wxgll(i)*wzgll(j)*jacobianl / TWO
         else
                   kinetic_energy = kinetic_energy +  &
              rhol_s*(velocs_poroelastic(1,ibool(i,j,ispec))**2 + velocs_poroelastic(2,ibool(i,j,ispec))**2) &
                 *wxgll(i)*wzgll(j)*jacobianl / TWO
         endif
        enddo
      enddo

!---
!--- acoustic spectral element
!---
    else

! for the definition of potential energy in an acoustic fluid, see for instance
! equation (23) of M. Maess et al., Journal of Sound and Vibration 296 (2006) 264-276

! in case of an acoustic medium, a potential Chi of (density * displacement) is used as in Chaljub and Valette,
! Geophysical Journal International, vol. 158, p. 131-141 (2004) and *NOT* a velocity potential
! as in Komatitsch and Tromp, Geophysical Journal International, vol. 150, p. 303-318 (2002).
! This permits acoustic-elastic coupling based on a non-iterative time scheme.
! Displacement is then: u = grad(Chi) / rho
! Velocity is then: v = grad(Chi_dot) / rho (Chi_dot being the time derivative of Chi)
! and pressure is: p = - Chi_dot_dot  (Chi_dot_dot being the time second derivative of Chi).

! compute pressure in this element
    call compute_pressure_one_element(pressure_element,potential_dot_dot_acoustic,displ_elastic,elastic, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,assign_external_model, &
         numat,kmato,elastcoef,vpext,vsext,rhoext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,anisotropic,anisotropy, &
         ispec,e1,e11, &
         TURN_ATTENUATION_ON,Mu_nu1,Mu_nu2,N_SLS)

! compute velocity vector field in this element
    call compute_vector_one_element(vector_field_element,potential_dot_acoustic,veloc_elastic,elastic, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,ispec,numat,kmato,density,rhoext,assign_external_model)

! get density of current spectral element
      lambdal_relaxed = elastcoef(1,1,kmato(ispec))
      mul_relaxed = elastcoef(2,1,kmato(ispec))
      rhol  = density(1,kmato(ispec))
      kappal  = lambdal_relaxed + TWO*mul_relaxed/3._CUSTOM_REAL
      cpl = sqrt((kappal + 4._CUSTOM_REAL*mul_relaxed/3._CUSTOM_REAL)/rhol)

! double loop over GLL points
      do j = 1,NGLLZ
        do i = 1,NGLLX

!--- if external medium, get density of current grid point
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

! compute kinetic energy
          kinetic_energy = kinetic_energy + &
              rhol*(vector_field_element(1,i,j)**2 + &
                    vector_field_element(2,i,j)**2) *wxgll(i)*wzgll(j)*jacobianl / TWO

! compute potential energy
          potential_energy = potential_energy + (pressure_element(i,j)**2)*wxgll(i)*wzgll(j)*jacobianl / (TWO * rhol * cpl**2)

        enddo
      enddo

    endif

  enddo

! save kinetic, potential and total energy for this time step in external file
  write(IOUT_ENERGY,*) real(dble(it-1)*deltat - t0,4),real(kinetic_energy,4), &
                     real(potential_energy,4),real(kinetic_energy + potential_energy,4)

  end subroutine compute_energy

