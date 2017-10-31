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

! for poro solver

 subroutine compute_coupling_poro_ac(potential_dot_dot_acoustic,accels_poroelastic,accelw_poroelastic,FORWARD_OR_ADJOINT)

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLZ,NGLJ,CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP,ONE

  use specfem_par, only: num_fluid_poro_edges,ibool,wxgll,wzgll,xix,xiz, &
                         gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse, &
                         fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge, &
                         fluid_poro_poroelastic_ispec,fluid_poro_poroelastic_iedge, &
                         porosity,tortuosity,density,kmato, &
                         potential_acoustic_adj_coupling, &
                         nglob_acoustic,nglob_poroelastic, &
                         SIMULATION_TYPE

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic), intent(in) :: potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic), intent(inout) :: accels_poroelastic,accelw_poroelastic

  ! FORWARD_OR_ADJOINT == 1 for forward wavefield, FORWARD_OR_ADJOINT == 3 for backward/reconstruction wavefield
  integer,intent(in) :: FORWARD_OR_ADJOINT

  !local variable
  integer :: inum,ispec_acoustic,iedge_acoustic,ispec_poroelastic,iedge_poroelastic, &
             ipoin1D,i,j,ii2,jj2,iglob
  real(kind=CUSTOM_REAL) :: pressure,xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight
  real(kind=CUSTOM_REAL) :: fac_s,fac_w
  double precision :: phi,tort,rho_f,rho_s,rho_bar
  integer :: material

  ! loop on all the coupling edges
  do inum = 1,num_fluid_poro_edges
    ! get the edge of the acoustic element
    ispec_acoustic = fluid_poro_acoustic_ispec(inum)
    iedge_acoustic = fluid_poro_acoustic_iedge(inum)
    ! get the corresponding edge of the poroelastic element
    ispec_poroelastic = fluid_poro_poroelastic_ispec(inum)
    iedge_poroelastic = fluid_poro_poroelastic_iedge(inum)

    ! implement 1D coupling along the edge
    do ipoin1D = 1,NGLLX
      ! get point values for the acoustic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_acoustic)
      j = jvalue_inverse(ipoin1D,iedge_acoustic)
      iglob = ibool(i,j,ispec_acoustic)

      ! get poroelastic parameters
      material = kmato(ispec_poroelastic)
      phi = porosity(material)
      tort = tortuosity(material)
      rho_f = density(2,material)
      rho_s = density(1,material)
      rho_bar = (1.d0-phi)*rho_s + phi*rho_f

      fac_s = real((1.d0 - phi/tort),kind=CUSTOM_REAL)
      fac_w = real((1.d0 - rho_f/rho_bar),kind=CUSTOM_REAL)

      ! compute pressure on the fluid/porous medium edge
      pressure = - potential_dot_dot_acoustic(iglob)

      if (SIMULATION_TYPE == 3 .and. FORWARD_OR_ADJOINT == 1) then
        ! new definition of adjoint displacement and adjoint potential
        ! adjoint definition: pressure^\dagger = potential^\dagger
        pressure = potential_acoustic_adj_coupling(iglob)
      endif

      ! get point values for the poroelastic side
      ii2 = ivalue(ipoin1D,iedge_poroelastic)
      jj2 = jvalue(ipoin1D,iedge_poroelastic)
      iglob = ibool(ii2,jj2,ispec_poroelastic)

      ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
      ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
      ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
      ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
      ! Blackwell Science, page 110, equation (4.60).
      if (iedge_acoustic == ITOP) then
        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        weight = jacobian1D * wxgll(i)
      else if (iedge_acoustic == IBOTTOM) then
        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        weight = jacobian1D * wxgll(i)
      else if (iedge_acoustic == ILEFT) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)
      else if (iedge_acoustic == IRIGHT) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)
      endif

      ! contribution to the solid phase
      accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) + weight*nx*pressure * fac_s
      accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + weight*nz*pressure * fac_s
      ! contribution to the fluid phase
      accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) + weight*nx*pressure * fac_w
      accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) + weight*nz*pressure * fac_w
    enddo
  enddo

 end subroutine compute_coupling_poro_ac
