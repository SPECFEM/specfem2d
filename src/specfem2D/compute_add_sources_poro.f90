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
                         NSOURCES,source_type,anglesource,source_time_function, &
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

  do i_source= 1,NSOURCES
    ! if this processor core carries the source and the source element is elastic
    if (is_proc_source(i_source) == 1 .and. ispec_is_poroelastic(ispec_selected_source(i_source))) then
      material = kmato(ispec_selected_source(i_source))
      phi = porosity(material)
      tort = tortuosity(material)
      rho_s = density(1,material)
      rho_f = density(2,material)

      rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

      fac_s = real((1.d0 - phi/tort),kind=CUSTOM_REAL)
      fac_w = real((1.d0 - rho_f/rho_bar),kind=CUSTOM_REAL)

      ! collocated force
      if (source_type(i_source) == 1) then
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
      endif
    endif ! if this processor core carries the source and the source element is elastic
  enddo ! do i_source= 1,NSOURCES

 end subroutine compute_add_sources_poro

