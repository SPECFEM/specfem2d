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
!=====================================================================

  subroutine add_acoustic_forcing_at_rigid_boundary(potential_dot_dot_acoustic,dot_e1)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: nglob_acoustic,nelem_acforcing,codeacforcing,numacforcing,ispec_is_acoustic, &
                         ibool,xix,xiz,jacobian,gammax,gammaz,wxgll,wzgll, &
                         ATTENUATION_VISCOACOUSTIC,N_SLS,nglob_att

  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(nglob_att,N_SLS) :: dot_e1

  !local variables
  integer :: inum,ispec,i,j,iglob
  real(kind=CUSTOM_REAL) :: xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight, &
                            displ_x,displ_z,displ_n

  ! loop on all the forced edges

  do inum = 1,nelem_acforcing

    ispec = numacforcing(inum)
    if (.not. ispec_is_acoustic(ispec)) cycle ! acoustic spectral element

    !--- left acoustic forcing boundary
    if (codeacforcing(IEDGE4,inum)) then
      i = 1
      do j = 1,NGLLZ
        iglob = ibool(i,j,ispec)
        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)

        ! define displacement components which will force the boundary
        if (PML_BOUNDARY_CONDITIONS) then
          if (ispec_is_PML(ispec)) then
            displ_x = 0.0d0
            displ_z = 0.0d0
          else
            call acoustic_forcing_boundary(iglob,displ_x,displ_z)
          endif
        else
          call acoustic_forcing_boundary(iglob,displ_x,displ_z)
        endif
        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n
        if (ATTENUATION_VISCOACOUSTIC) dot_e1(iglob,:) = dot_e1(iglob,:) + weight*displ_n
      enddo
    endif  !  end of left acoustic forcing boundary

    !--- right acoustic forcing boundary
    if (codeacforcing(IEDGE2,inum)) then
      i = NGLLX
      do j = 1,NGLLZ
        iglob = ibool(i,j,ispec)
        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)

        ! define displacement components which will force the boundary
        if (PML_BOUNDARY_CONDITIONS) then
          if (ispec_is_PML(ispec)) then
            displ_x = 0.0d0
            displ_z = 0.0d0
          else
            call acoustic_forcing_boundary(iglob,displ_x,displ_z)
          endif
        else
          call acoustic_forcing_boundary(iglob,displ_x,displ_z)
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n
        if (ATTENUATION_VISCOACOUSTIC) dot_e1(iglob,:) = dot_e1(iglob,:) + weight*displ_n

      enddo
    endif  !  end of right acoustic forcing boundary

    !--- bottom acoustic forcing boundary
    if (codeacforcing(IEDGE1,inum)) then
      j = 1
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        weight = jacobian1D * wxgll(i)

        ! define displacement components which will force the boundary
        if (PML_BOUNDARY_CONDITIONS) then
          if (ispec_is_PML(ispec)) then
            displ_x = 0.0d0
            displ_z = 0.0d0
          else
            call acoustic_forcing_boundary(iglob,displ_x,displ_z)
          endif
        else
          call acoustic_forcing_boundary(iglob,displ_x,displ_z)
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n
        if (ATTENUATION_VISCOACOUSTIC) dot_e1(iglob,:) = dot_e1(iglob,:) + weight*displ_n
      enddo
    endif  !  end of bottom acoustic forcing boundary

    !--- top acoustic forcing boundary
    if (codeacforcing(IEDGE3,inum)) then
      j = NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        weight = jacobian1D * wxgll(i)

        ! define displacement components which will force the boundary
        if (PML_BOUNDARY_CONDITIONS) then
          if (ispec_is_PML(ispec)) then
            displ_x = 0.0d0
            displ_z = 0.0d0
          else
            call acoustic_forcing_boundary(iglob,displ_x,displ_z)
          endif
        else
          call acoustic_forcing_boundary(iglob,displ_x,displ_z)
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n
        if (ATTENUATION_VISCOACOUSTIC) dot_e1(iglob,:) = dot_e1(iglob,:) + weight*displ_n
      enddo
    endif  !  end of top acoustic forcing boundary
  enddo

  end subroutine add_acoustic_forcing_at_rigid_boundary

