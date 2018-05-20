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

  subroutine compute_and_output_energy()

  use constants, only: IOUT_ENERGY,CUSTOM_REAL

  use specfem_par, only: GPU_MODE,myrank,it,deltat,kinetic_energy,potential_energy,t0

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: kinetic_energy_total,potential_energy_total

  ! safety check
  if (GPU_MODE) call stop_the_code('Error computing energy for output is not implemented on GPUs yet')

  ! computes energy
  call compute_energy()

  ! computes total for all processes
  call sum_all_cr(kinetic_energy,kinetic_energy_total)
  call sum_all_cr(potential_energy,potential_energy_total)

  ! saves kinetic, potential and total energy for this time step in external file
  if (myrank == 0) then
    write(IOUT_ENERGY,*) real(dble(it-1)*deltat - t0,4),real(kinetic_energy_total,4), &
                         real(potential_energy_total,4),real(kinetic_energy_total + potential_energy_total,4)
  endif

  end subroutine compute_and_output_energy

!
!----------------------------------------------------------------------------------------
!

  subroutine compute_energy()

! compute kinetic and potential energy in the solid (acoustic elements are excluded)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM,ZERO,TWO


  use specfem_par, only: AXISYM,is_on_the_axis,myrank,nspec,kinetic_energy,potential_energy, &
    ibool,hprime_xx,hprime_zz,hprimeBar_xx,xix,xiz,gammax,gammaz,jacobian,wxgll,wzgll, &
    displ_elastic,veloc_elastic, &
    displs_poroelastic,displw_poroelastic,velocs_poroelastic,velocw_poroelastic,potential_acoustic, &
    potential_dot_acoustic,potential_dot_dot_acoustic,vsext,vpext,rhoext,poroelastcoef,density,kmato,assign_external_model, &
    ispec_is_poroelastic,ispec_is_elastic,P_SV,ispec_is_PML
  implicit none

! local variables
  integer :: i,j,k,ispec
  real(kind=CUSTOM_REAL) :: cpl,csl,kappal
  real(kind=CUSTOM_REAL) :: mu_G,lambdal_G,lambdalplus2mul_G

  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl
  double precision :: rhol
  ! to evaluate cpI, cpII, and cs, and rI (poroelastic medium)
  double precision :: phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot
  double precision :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic

  double precision :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  double precision :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  double precision :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  double precision :: dwx_dxl,dwz_dxl,dwx_dzl,dwz_dzl
  ! vector field in an element
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: vector_field_element
  ! pressure in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: pressure_element

!! DK DK March 2018: restored this initialization that someone had removed for some reason, thus making the calculation wrong
!! DK DK March 2018: please do *NOT* remove it
  kinetic_energy = 0._CUSTOM_REAL
  potential_energy = 0._CUSTOM_REAL

  ! loop over spectral elements
  do ispec = 1,nspec

!! DK DK March 2018: only compute energy in the main domain, not in the PMLs, in which it is not physical
    if (ispec_is_PML(ispec)) cycle

    !---
    !--- elastic spectral element
    !---
    if (ispec_is_elastic(ispec)) then

      ! checks wave type
      if (.not. P_SV) then
        call exit_MPI(myrank,'output energy for SH waves not implemented yet')
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
          if (assign_external_model) then
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

          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              do k = 1,NGLJ
                dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprimeBar_xx(i,k)
              enddo
            else
              do k = 1,NGLJ
                dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
                duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
              enddo
            endif
          else
            do k = 1,NGLLX
              dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
              duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            enddo
          endif

          do k = 1,NGLLX
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

          ! compute kinetic energy ! TODO ABAB This integral over space should be adapted for axisym geometries if needed
          kinetic_energy = kinetic_energy  &
              + rhol * (veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) &
              *wxgll(i)*wzgll(j)*jacobianl / TWO

          ! compute potential energy ! TODO ABAB This integral over space should be adapted for axisym geometries if needed
          potential_energy = potential_energy &
              + (lambdaplus2mu_unrelaxed_elastic*dux_dxl**2 &
              + lambdaplus2mu_unrelaxed_elastic*duz_dzl**2 &
              + TWO*lambdal_unrelaxed_elastic*dux_dxl*duz_dzl &
              + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2)*wxgll(i)*wzgll(j)*jacobianl / TWO

        enddo
      enddo

    !---
    !--- poroelastic spectral element
    !---
    else if (ispec_is_poroelastic(ispec)) then

      ! get unrelaxed elastic parameters of current spectral element
      !for now replaced by solid, fluid, and frame parameters of current spectral element

      ! gets poroelastic material
      call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)

      ! Biot coefficients for the input phi
      call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

      !The RHS has the form : div T -phi/c div T_f + phi/ceta_fk^-1.partial t w
      !where T = G:grad u_s + C div w I
      !and T_f = C div u_s I + M div w I
      !we are expressing lambdaplus2mu, lambda, and mu for G, C, and M
      mu_G = mu_fr
      lambdal_G = H_biot - TWO*mu_fr
      lambdalplus2mul_G = lambdal_G + TWO*mu_G

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
              + TWO*lambdal_G*dux_dxl*duz_dzl + mu_G*(dux_dzl + duz_dxl)**2 &
              + TWO*C_biot*dwx_dxl*dux_dxl + TWO*C_biot*dwz_dzl*duz_dzl &
              + TWO*C_biot*(dwx_dxl*duz_dzl + dwz_dzl*dux_dxl) &
              + M_biot*dwx_dxl**2 + M_biot*dwz_dzl**2 &
              + TWO*M_biot*dwx_dxl*dwz_dzl )*wxgll(i)*wzgll(j)*jacobianl / TWO

          ! compute kinetic energy
          if (phi > 0.0d0) then
            kinetic_energy = kinetic_energy &
              + ( rho_bar*(velocs_poroelastic(1,ibool(i,j,ispec))**2 &
              + velocs_poroelastic(2,ibool(i,j,ispec))**2) &
              + rho_f*tort/phi*(velocw_poroelastic(1,ibool(i,j,ispec))**2 &
              + velocw_poroelastic(2,ibool(i,j,ispec))**2) &
              + rho_f*(velocs_poroelastic(1,ibool(i,j,ispec))*velocw_poroelastic(1,ibool(i,j,ispec)) &
              + velocs_poroelastic(2,ibool(i,j,ispec))*velocw_poroelastic(2,ibool(i,j,ispec))) &
                 )*wxgll(i)*wzgll(j)*jacobianl / TWO
          else
            kinetic_energy = kinetic_energy  &
              + rho_s*(velocs_poroelastic(1,ibool(i,j,ispec))**2 &
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
      call compute_pressure_one_element(ispec,pressure_element,displ_elastic,displs_poroelastic,displw_poroelastic, &
                                        potential_dot_dot_acoustic,potential_acoustic)

      ! compute velocity vector field in this element
      call compute_vector_one_element(potential_dot_acoustic,veloc_elastic,velocs_poroelastic,ispec,vector_field_element)

      ! get velocity and density in current spectral element
      if (.not. assign_external_model) then
        lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
        rhol  = density(1,kmato(ispec))
        kappal  = lambdal_unrelaxed_elastic
        cpl = sqrt(kappal/rhol)
      endif

      ! double loop over GLL points
      do j = 1,NGLLZ
        do i = 1,NGLLX

          !--- if external medium, get density of current grid point
          if (assign_external_model) then
            cpl = vpext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          jacobianl = jacobian(i,j,ispec)

          ! compute kinetic energy
          kinetic_energy = kinetic_energy &
              + rhol*(vector_field_element(1,i,j)**2 + vector_field_element(2,i,j)**2) *wxgll(i)*wzgll(j)*jacobianl / TWO

          ! compute potential energy
          potential_energy = potential_energy &
              + (pressure_element(i,j)**2)*wxgll(i)*wzgll(j)*jacobianl / (TWO * rhol * cpl**2)

        enddo
      enddo

    endif
  enddo

  end subroutine compute_energy

!
!----------------------------------------------------------------------------------------
!

  subroutine compute_energy_fields()

  ! computes maximum, integrated energy and duration fields

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,TWO,ZERO

  use specfem_par, only: AXISYM,is_on_the_axis,nspec,ibool,deltat,veloc_elastic,potential_dot_acoustic,ispec_is_elastic, &
                        ispec_is_poroelastic,integrated_kinetic_energy_field,max_kinetic_energy_field, &
                        integrated_potential_energy_field,max_potential_energy_field,kinetic_effective_duration_field, &
                        potential_effective_duration_field,total_integrated_energy_field,max_total_energy_field, &
                        total_effective_duration_field,velocs_poroelastic, &
                        poroelastcoef,vsext,vpext,rhoext,density,kmato,assign_external_model,jacobian,displ_elastic, &
                        hprime_xx,hprime_zz,hprimeBar_xx,xix,xiz,gammax,gammaz, &
                        displs_poroelastic,displw_poroelastic, &
                        potential_dot_dot_acoustic,potential_acoustic

  implicit none

  ! local variables
  integer :: ispec,i,j,k
  real(kind=CUSTOM_REAL) :: cpl,csl

  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl
  double precision :: rhol
  double precision :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic
  double precision :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  double precision :: dux_dxl,duz_dxl,dux_dzl,duz_dzl

  ! vector field in an element
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: vector_field_element
  ! pressure in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: pressure_element

  ! We save the value at the GLL point:
  i=2
  j=2

  ! loop over spectral elements
  do ispec = 1,nspec

    !---
    !--- elastic spectral element
    !---
    if (ispec_is_elastic(ispec)) then

      ! get relaxed elastic parameters of current spectral element
      lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
      mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
      lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))

      rhol  = density(1,kmato(ispec))

      !--- if external medium, get elastic parameters of current grid point
      if (assign_external_model) then
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
      if (AXISYM) then
        if (is_on_the_axis(ispec)) then
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprimeBar_xx(i,k)
          enddo
        else
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          enddo
        endif
      else
        do k = 1,NGLLX
          dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
        enddo
      endif
      do k = 1,NGLLX
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

      ! compute total integrated energy
      ! We record just one point per element (i=2, j=2)
      integrated_kinetic_energy_field(ispec) = integrated_kinetic_energy_field(ispec)  &
          +  rhol*(veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) * deltat / TWO

      if (max_kinetic_energy_field(ispec) < rhol*(veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) / &
          TWO) then
        max_kinetic_energy_field(ispec) = rhol*(veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) / TWO
      endif

      if (max_kinetic_energy_field(ispec) > ZERO) then
        kinetic_effective_duration_field(ispec) = TWO*integrated_kinetic_energy_field(ispec)/max_kinetic_energy_field(ispec)
      endif

      integrated_potential_energy_field(ispec) = integrated_potential_energy_field(ispec) &
              + (lambdaplus2mu_unrelaxed_elastic*dux_dxl**2 &
              + lambdaplus2mu_unrelaxed_elastic*duz_dzl**2 &
              + TWO*lambdal_unrelaxed_elastic*dux_dxl*duz_dzl &
              + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2) * deltat / TWO
      if (max_potential_energy_field(ispec) < (lambdaplus2mu_unrelaxed_elastic*dux_dxl**2 &
            + lambdaplus2mu_unrelaxed_elastic*duz_dzl**2 &
            + TWO*lambdal_unrelaxed_elastic*dux_dxl*duz_dzl &
            + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2) / TWO) then
        max_potential_energy_field(ispec) = (lambdaplus2mu_unrelaxed_elastic*dux_dxl**2 &
          + lambdaplus2mu_unrelaxed_elastic*duz_dzl**2 &
          + TWO*lambdal_unrelaxed_elastic*dux_dxl*duz_dzl &
          + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2) / TWO
      endif

      if (max_potential_energy_field(ispec) > ZERO) then
        potential_effective_duration_field(ispec) = TWO*integrated_potential_energy_field(ispec) / &
                                                    max_potential_energy_field(ispec)
      endif

      if (max_total_energy_field(ispec) < (lambdaplus2mu_unrelaxed_elastic*dux_dxl**2 &
            + lambdaplus2mu_unrelaxed_elastic*duz_dzl**2 &
            + TWO*lambdal_unrelaxed_elastic*dux_dxl*duz_dzl &
            + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2) / TWO + &
            rhol*(veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) / TWO) then
        max_total_energy_field(ispec) = (lambdaplus2mu_unrelaxed_elastic*dux_dxl**2 &
            + lambdaplus2mu_unrelaxed_elastic*duz_dzl**2 &
            + TWO*lambdal_unrelaxed_elastic*dux_dxl*duz_dzl &
            + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2) / TWO + &
            rhol*(veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) / TWO
      endif

    !---
    !--- poroelastic spectral element
    !---
    else if (ispec_is_poroelastic(ispec)) then
       ! safety check
       call stop_the_code( &
'COMPUTE_INTEGRATED_ENERGY_FIELD is not available for poroelastic media yet (but it would be very easy to implement)')

    !---
    !--- acoustic spectral element
    !---
    else

      ! compute velocity vector field in this element
      call compute_vector_one_element(potential_dot_acoustic,veloc_elastic,velocs_poroelastic,ispec,vector_field_element)

      ! compute pressure in this element
      call compute_pressure_one_element(ispec,pressure_element,displ_elastic,displs_poroelastic,displw_poroelastic, &
                                        potential_dot_dot_acoustic,potential_acoustic)

      !--- if external medium, get density of current grid point
      if (assign_external_model) then
        rhol = rhoext(i,j,ispec)
        cpl = vpext(i,j,ispec)
      else
        lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
        rhol  = density(1,kmato(ispec))
        cpl = sqrt(lambdal_unrelaxed_elastic/rhol) !lambdal_unrelaxed_elastic = kappal
      endif

      jacobianl = jacobian(i,j,ispec)

      ! compute total integrated energy ! = int_0^t v^2 dt
      integrated_kinetic_energy_field(ispec) = integrated_kinetic_energy_field(ispec)  &
           +  rhol * (vector_field_element(1,i,j)**2 + vector_field_element(2,i,j)**2) * deltat / TWO

      if (max_kinetic_energy_field(ispec) < rhol * (vector_field_element(1,i,j)**2 + vector_field_element(2,i,j)**2) / TWO) then
        max_kinetic_energy_field(ispec) = rhol * (vector_field_element(1,i,j)**2 + vector_field_element(2,i,j)**2) / TWO
      endif

      if (max_kinetic_energy_field(ispec) > ZERO) then
        kinetic_effective_duration_field(ispec) = TWO*integrated_kinetic_energy_field(ispec)/max_kinetic_energy_field(ispec)
      endif

      ! compute potential energy
      integrated_potential_energy_field(ispec) = integrated_potential_energy_field(ispec) + pressure_element(i,j)**2 * deltat &
                                                 / (TWO*rhol*cpl**2)
      if (max_potential_energy_field(ispec) < pressure_element(i,j)**2 / (TWO*rhol*cpl**2)) then
        max_potential_energy_field(ispec) = pressure_element(i,j)**2 / (TWO*rhol*cpl**2)
      endif

      if (max_potential_energy_field(ispec) > ZERO) then
        potential_effective_duration_field(ispec) = TWO*integrated_potential_energy_field(ispec) / &
                                                    max_potential_energy_field(ispec)
      endif

      if (max_total_energy_field(ispec) < pressure_element(i,j)**2 / (TWO*rhol*cpl**2) + &
            rhol * (vector_field_element(1,i,j)**2 + vector_field_element(2,i,j)**2) / TWO) then
        max_total_energy_field(ispec) = pressure_element(i,j)**2 / (TWO*rhol*cpl**2) + &
            rhol * (vector_field_element(1,i,j)**2 + vector_field_element(2,i,j)**2) / TWO
      endif

    endif

    total_integrated_energy_field(ispec) = integrated_kinetic_energy_field(ispec) + integrated_potential_energy_field(ispec)

    if (max_total_energy_field(ispec) > ZERO) then
      total_effective_duration_field(ispec) = TWO*total_integrated_energy_field(ispec) / max_total_energy_field(ispec)
    endif

  enddo

  end subroutine compute_energy_fields


