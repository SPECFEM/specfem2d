!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour and CNRS, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
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

  subroutine compute_forces_poro_viscous_damping()

! viscous damping

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,ZERO,USE_PORO_VISCOUS_DAMPING

  use specfem_par, only: nspec,it,SIMULATION_TYPE,SAVE_FORWARD, &
                         ATTENUATION_PORO_FLUID_PART, &
                         ibool,kmato,ispec_is_poroelastic, &
                         accels_poroelastic, &
                         accelw_poroelastic, &
                         velocw_poroelastic, &
                         porosity,tortuosity,permeability,poroelastcoef, &
                         jacobian,wxgll,wzgll, &
                         rx_viscous,rz_viscous,theta_e,theta_s, &
                         b_viscodampx,b_viscodampz

  implicit none

  ! local variables
  integer :: ispec,i,j,iglob

  double precision, dimension(3):: bl_relaxed_viscoelastic,bl_unrelaxed_elastic
  double precision :: eta_f,phi,tort
  double precision :: permlxx,permlxz,permlzz,invpermlxx,invpermlxz,invpermlzz,detk
  double precision :: viscodampx,viscodampz
  double precision :: visc_x,visc_z

  ! checks if anything to do
  if (.not. USE_PORO_VISCOUS_DAMPING) return

  ! loop over spectral elements
  do ispec = 1,nspec

    eta_f = poroelastcoef(2,2,kmato(ispec))

    if (ispec_is_poroelastic(ispec) .and. eta_f > 0.d0) then
      phi = porosity(kmato(ispec))
      tort = tortuosity(kmato(ispec))

      permlxx = permeability(1,kmato(ispec))
      permlxz = permeability(2,kmato(ispec))
      permlzz = permeability(3,kmato(ispec))

      ! calcul of the inverse of k
      detk = permlxx*permlzz - permlxz*permlxz

      if (detk /= ZERO) then
        invpermlxx = permlzz/detk
        invpermlxz = -permlxz/detk
        invpermlzz = permlxx/detk
      else
        call stop_the_code('Permeability matrix is not invertible')
      endif

      if (ATTENUATION_PORO_FLUID_PART) then
        ! viscous attenuation is implemented following the memory variable formulation of
        ! J. M. Carcione Wave fields in real media: wave propagation in anisotropic,
        ! anelastic and porous media, Elsevier, p. 304-305, 2007
        bl_relaxed_viscoelastic(1) = eta_f*invpermlxx*theta_e/theta_s
        bl_relaxed_viscoelastic(2) = eta_f*invpermlxz*theta_e/theta_s
        bl_relaxed_viscoelastic(3) = eta_f*invpermlzz*theta_e/theta_s
      else
        ! relaxed viscous coef
        bl_unrelaxed_elastic(1) = eta_f*invpermlxx
        bl_unrelaxed_elastic(2) = eta_f*invpermlxz
        bl_unrelaxed_elastic(3) = eta_f*invpermlzz
      endif

      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)

          ! computes the viscous damping term
          if (ATTENUATION_PORO_FLUID_PART) then
            ! with the unrelaxed viscous coef and add memory variable
            viscodampx = velocw_poroelastic(1,iglob)*bl_relaxed_viscoelastic(1) &
                       + velocw_poroelastic(2,iglob)*bl_relaxed_viscoelastic(2) - rx_viscous(i,j,ispec)
            viscodampz = velocw_poroelastic(1,iglob)*bl_relaxed_viscoelastic(2) &
                       + velocw_poroelastic(2,iglob)*bl_relaxed_viscoelastic(3) - rz_viscous(i,j,ispec)
          else
            ! no viscous attenuation
            viscodampx = velocw_poroelastic(1,iglob)*bl_unrelaxed_elastic(1) &
                       + velocw_poroelastic(2,iglob)*bl_unrelaxed_elastic(2)
            viscodampz = velocw_poroelastic(1,iglob)*bl_unrelaxed_elastic(2) &
                       + velocw_poroelastic(2,iglob)*bl_unrelaxed_elastic(3)
          endif

          ! viscous term eta_f k^-1 dot(w)
          visc_x = wxgll(i) * wzgll(j) * jacobian(i,j,ispec) * viscodampx
          visc_z = wxgll(i) * wzgll(j) * jacobian(i,j,ispec) * viscodampz

          ! solid contribution
          ! adds term + phi/tort eta_f k^-1 dot(w)
          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) + phi/tort * visc_x
          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + phi/tort * visc_z

          ! fluid contribution
          ! add - eta_f k^-1 dot(w)
          accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - visc_x
          accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - visc_z

          ! reconstructed/backward wavefield
          ! viscous damping contribution
          if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
            ! stores contribution
            b_viscodampx(i,j,ispec) = visc_x
            b_viscodampz(i,j,ispec) = visc_z
          endif

        enddo
      enddo

    endif ! end of test if poroelastic element

  enddo ! end of loop over all spectral elements

  ! saves viscous contribution to disk
  if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
    ! writes damping contributions to file
    write(23,rec=it) b_viscodampx(:,:,:)
    write(24,rec=it) b_viscodampz(:,:,:)
  endif

  end subroutine compute_forces_poro_viscous_damping


!
!-------------------------------------------------------------------------------------
!

  subroutine compute_forces_poro_viscous_damping_backward()

! viscous damping

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,ZERO,USE_PORO_VISCOUS_DAMPING

  use specfem_par, only: nspec,NSTEP,it,SIMULATION_TYPE, &
                         ibool,kmato,ispec_is_poroelastic, &
                         b_accels_poroelastic, &
                         b_accelw_poroelastic, &
                         porosity,tortuosity,poroelastcoef, &
                         b_viscodampx,b_viscodampz

  implicit none

  ! local variables
  integer :: ispec,i,j,iglob

  double precision :: eta_f,phi,tort
  double precision :: visc_x,visc_z

  ! checks if anything to do
  if (.not. USE_PORO_VISCOUS_DAMPING) return
  if (SIMULATION_TYPE /= 3) return

  ! reads in viscous contributions for reconstructed/backward wavefield
  read(23,rec=NSTEP-it+1) b_viscodampx(:,:,:)
  read(24,rec=NSTEP-it+1) b_viscodampz(:,:,:)

  ! loop over spectral elements
  do ispec = 1,nspec

    eta_f = poroelastcoef(2,2,kmato(ispec))

    if (ispec_is_poroelastic(ispec) .and. eta_f > 0.d0) then

      phi = porosity(kmato(ispec))
      tort = tortuosity(kmato(ispec))

      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)

          ! reconstructed/backward wavefield
          ! viscous damping contribution
          ! kernels simulation uses previously stored contributions
          visc_x = b_viscodampx(i,j,ispec)
          visc_z = b_viscodampz(i,j,ispec)

          ! solid
          b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) + phi/tort * visc_x
          b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + phi/tort * visc_z
          ! fluid
          b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) - visc_x
          b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) - visc_z

        enddo
      enddo

    endif ! end of test if poroelastic element

  enddo ! end of loop over all spectral elements

  end subroutine compute_forces_poro_viscous_damping_backward

