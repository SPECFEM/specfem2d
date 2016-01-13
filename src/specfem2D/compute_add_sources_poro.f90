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
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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
!=====================================================================

! for poro solver

  subroutine compute_add_sources_poro(accels_poroelastic,accelw_poroelastic,it,i_stage)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM

  use specfem_par, only: ispec_is_poroelastic,nglob_poroelastic, &
                         NSOURCES,source_type,anglesource,source_time_function,sourcearray, &
                         is_proc_source,ispec_selected_source, &
                         hxis_store,hgammas_store,ibool, &
                         porosity,tortuosity,density,kmato
  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic) :: accels_poroelastic,accelw_poroelastic

  integer :: it,i_stage

  !local variables
  integer :: i_source,i,j,iglob
  double precision :: hlagrange
  double precision :: phi,tort,rho_s,rho_f,rho_bar
  real(kind=CUSTOM_REAL) :: fac_s,fac_w
  integer :: material

  do i_source = 1,NSOURCES

    ! if this processor core carries the source and the source element is poroelastic
    if (is_proc_source(i_source) == 1 .and. ispec_is_poroelastic(ispec_selected_source(i_source))) then

      material = kmato(ispec_selected_source(i_source))
      phi = porosity(material)
      tort = tortuosity(material)
      rho_s = density(1,material)
      rho_f = density(2,material)

      rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

      fac_s = real((1.d0 - phi/tort),kind=CUSTOM_REAL)
      fac_w = real((1.d0 - rho_f/rho_bar),kind=CUSTOM_REAL)

      select case (source_type(i_source))

      case (1)
        ! collocated force
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec_selected_source(i_source))
            hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
            ! s
            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - hlagrange * &
                   fac_s * sin(anglesource(i_source))*source_time_function(i_source,it,i_stage)
            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + hlagrange * &
                   fac_s * cos(anglesource(i_source))*source_time_function(i_source,it,i_stage)
            ! w
            accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - hlagrange * &
                   fac_w * sin(anglesource(i_source))*source_time_function(i_source,it,i_stage)
            accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) + hlagrange * &
                   fac_w * cos(anglesource(i_source))*source_time_function(i_source,it,i_stage)
          enddo
        enddo

      case(2)
        ! moment tensor
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec_selected_source(i_source))

            ! solid contribution
            accels_poroelastic(:,iglob) = accels_poroelastic(:,iglob) + &
                        (1._CUSTOM_REAL - phi/tort)*sourcearray(i_source,:,i,j) &
                        * source_time_function(i_source,it,i_stage)

            ! fluid contribution
            accelw_poroelastic(:,iglob) = accelw_poroelastic(:,iglob) + &
                    (1._CUSTOM_REAL - rho_f/rho_bar)*sourcearray(i_source,:,i,j) * &
                    source_time_function(i_source,it,i_stage)
          enddo
        enddo
      end select

    endif ! if this processor core carries the source and the source element is elastic
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_poro

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_add_sources_poro_adjoint()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: myrank,nrec,NSTEP,it, &
                         ispec_is_poroelastic,which_proc_receiver,ispec_selected_rec, &
                         ibool,adj_sourcearrays,initialfield,SIMULATION_TYPE, &
                         kmato,porosity,density, &
                         accels_poroelastic,accelw_poroelastic

  implicit none

  ! local parameters
  integer :: irec_local,irec,i,j,iglob
  double precision :: phi,rho_s,rho_f,rho_bar

  ! checks if anything to do
  if (initialfield) return

  ! only for adjoint/kernel simulations
  if (.not. SIMULATION_TYPE == 3) return

  ! adjoint wavefield
  irec_local = 0
  do irec = 1,nrec
    ! add the source (only if this proc carries the source)
    if (myrank == which_proc_receiver(irec)) then
      irec_local = irec_local + 1

      if (ispec_is_poroelastic(ispec_selected_rec(irec))) then

        phi = porosity(kmato(ispec_selected_rec(irec)))
        rho_s = density(1,kmato(ispec_selected_rec(irec)))
        rho_f = density(2,kmato(ispec_selected_rec(irec)))

        rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

        ! add source array
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec_selected_rec(irec))

            ! solid contribution
            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) + adj_sourcearrays(irec_local,NSTEP-it+1,1,i,j)
            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + adj_sourcearrays(irec_local,NSTEP-it+1,2,i,j)

            ! fluid
            accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - &
                  rho_f/rho_bar*adj_sourcearrays(irec_local,NSTEP-it+1,1,i,j)
            accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - &
                  rho_f/rho_bar*adj_sourcearrays(irec_local,NSTEP-it+1,2,i,j)
          enddo
        enddo

      endif ! if element is poroelastic
    endif ! if this processor core carries the adjoint source and the source element is poroelastic
  enddo ! irec = 1,nrec

  end subroutine compute_add_sources_poro_adjoint
