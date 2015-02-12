
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

  subroutine compute_energy()

! compute kinetic and potential energy in the solid (acoustic elements are excluded)

  use specfem_par

  implicit none


! local variables
  integer :: i,j,k,ispec

  real(kind=CUSTOM_REAL) :: cpl,csl,kappal


  kinetic_energy = ZERO
  potential_energy = ZERO

! loop over spectral elements
  do ispec = 1,nspec

    !---
    !--- elastic spectral element
    !---
    if(elastic(ispec)) then

      ! checks wave type
      if( .not. p_sv ) then
        call exit_MPI('output energy for SH waves not implemented yet')
      endif

      ! get relaxed elastic parameters of current spectral element
      lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
      mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
      lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))
      rhol  = density(1,kmato(ispec))

      ! double loop over GLL points
      do j = 1,NGLLZ
        do i = 1,NGLLX

          !--- if external medium, get elastic parameters of current grid point
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
            mul_unrelaxed_elastic = rhol*csl*csl
            lambdal_unrelaxed_elastic = rhol*cpl*cpl - TWO*mul_unrelaxed_elastic
            lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic
          endif

          ! derivative along x and along z
          dux_dxi = 0._CUSTOM_REAL
          duz_dxi = 0._CUSTOM_REAL

          dux_dgamma = 0._CUSTOM_REAL
          duz_dgamma = 0._CUSTOM_REAL

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)
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
          kinetic_energy = kinetic_energy  &
              + rhol*(veloc_elastic(1,ibool(i,j,ispec))**2  &
              + veloc_elastic(3,ibool(i,j,ispec))**2) *wxgll(i)*wzgll(j)*jacobianl / TWO

          ! compute potential energy
          potential_energy = potential_energy &
              + (lambdaplus2mu_unrelaxed_elastic*dux_dxl**2 &
              + lambdaplus2mu_unrelaxed_elastic*duz_dzl**2 &
              + two*lambdal_unrelaxed_elastic*dux_dxl*duz_dzl &
              + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2)*wxgll(i)*wzgll(j)*jacobianl / TWO

        enddo
      enddo

    !---
    !--- poroelastic spectral element
    !---
    else if(poroelastic(ispec)) then

      ! get unrelaxed elastic parameters of current spectral element
      !for now replaced by solid, fluid, and frame parameters of current spectral element
      phil = porosity(kmato(ispec))
      tortl = tortuosity(kmato(ispec))
      !solid properties
      mul_s = poroelastcoef(2,1,kmato(ispec))
      kappal_s = poroelastcoef(3,1,kmato(ispec)) - FOUR_THIRDS*mul_s
      rhol_s = density(1,kmato(ispec))
      !fluid properties
      kappal_f = poroelastcoef(1,2,kmato(ispec))
      rhol_f = density(2,kmato(ispec))
      !frame properties
      mul_fr = poroelastcoef(2,3,kmato(ispec))
      kappal_fr = poroelastcoef(3,3,kmato(ispec)) - FOUR_THIRDS*mul_fr
      rhol_bar =  (1.d0 - phil)*rhol_s + phil*rhol_f
      !Biot coefficients for the input phi
      D_biot = kappal_s*(1.d0 + phil*(kappal_s/kappal_f - 1.d0))
      H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) &
              + kappal_fr + FOUR_THIRDS*mul_fr
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
          potential_energy = potential_energy &
              + ( lambdalplus2mul_G*dux_dxl**2 &
              + lambdalplus2mul_G*duz_dzl**2 &
              + two*lambdal_G*dux_dxl*duz_dzl + mul_G*(dux_dzl + duz_dxl)**2 &
              + two*C_biot*dwx_dxl*dux_dxl + two*C_biot*dwz_dzl*duz_dzl &
              + two*C_biot*(dwx_dxl*duz_dzl + dwz_dzl*dux_dxl) &
              + M_biot*dwx_dxl**2 + M_biot*dwz_dzl**2 &
              + two*M_biot*dwx_dxl*dwz_dzl )*wxgll(i)*wzgll(j)*jacobianl / TWO

          ! compute kinetic energy
          if(phil > 0.0d0) then
            kinetic_energy = kinetic_energy &
              + ( rhol_bar*(velocs_poroelastic(1,ibool(i,j,ispec))**2 &
              + velocs_poroelastic(2,ibool(i,j,ispec))**2) &
              + rhol_f*tortl/phil*(velocw_poroelastic(1,ibool(i,j,ispec))**2 &
              + velocw_poroelastic(2,ibool(i,j,ispec))**2) &
              + rhol_f*(velocs_poroelastic(1,ibool(i,j,ispec))*velocw_poroelastic(1,ibool(i,j,ispec)) &
              + velocs_poroelastic(2,ibool(i,j,ispec))*velocw_poroelastic(2,ibool(i,j,ispec))) &
                 )*wxgll(i)*wzgll(j)*jacobianl / TWO
          else
            kinetic_energy = kinetic_energy  &
              + rhol_s*(velocs_poroelastic(1,ibool(i,j,ispec))**2 &
              + velocs_poroelastic(2,ibool(i,j,ispec))**2)*wxgll(i)*wzgll(j)*jacobianl / TWO
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
      call compute_pressure_one_element(ispec)

      ! compute velocity vector field in this element
      call compute_vector_one_element(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                              potential_dot_gravito,veloc_elastic,velocs_poroelastic,ispec)

      ! get density of current spectral element
      lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
      rhol  = density(1,kmato(ispec))
      kappal  = lambdal_unrelaxed_elastic
      cpl = sqrt(kappal/rhol)

      ! double loop over GLL points
      do j = 1,NGLLZ
        do i = 1,NGLLX

          !--- if external medium, get density of current grid point
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          jacobianl = jacobian(i,j,ispec)

          ! compute kinetic energy
          kinetic_energy = kinetic_energy &
              + rhol*(vector_field_element(1,i,j)**2 &
              + vector_field_element(2,i,j)**2) *wxgll(i)*wzgll(j)*jacobianl / TWO

          ! compute potential energy
          potential_energy = potential_energy &
              + (pressure_element(i,j)**2)*wxgll(i)*wzgll(j)*jacobianl / (TWO * rhol * cpl**2)

        enddo
      enddo

    endif

  enddo

  end subroutine compute_energy

