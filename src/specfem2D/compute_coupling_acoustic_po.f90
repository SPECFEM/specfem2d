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

! for acoustic solver

 subroutine compute_coupling_acoustic_po(dot_e1)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP,ONE, &
                       USE_A_STRONG_FORMULATION_FOR_E1

  use specfem_par, only: num_fluid_poro_edges,ibool,wxgll,wzgll,xix,xiz, &
                         gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse, &
                         fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge, &
                         fluid_poro_poroelastic_ispec,fluid_poro_poroelastic_iedge, &
                         displs_poroelastic,displw_poroelastic, &
                         accels_poroelastic_adj_coupling,accelw_poroelastic_adj_coupling, &
                         potential_dot_dot_acoustic,SIMULATION_TYPE, &
                         ATTENUATION_VISCOACOUSTIC,N_SLS,nglob_att

  implicit none

  real(kind=CUSTOM_REAL),dimension(nglob_att,N_SLS) :: dot_e1

  ! local variables
  integer :: inum,ispec_acoustic,iedge_acoustic,ispec_poroelastic,iedge_poroelastic, &
             ipoin1D,i,j,iglob
  real(kind=CUSTOM_REAL) :: displ_x,displ_z,displw_x,displw_z,displ_n, &
                            xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight


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
      ! get point values for the poroelastic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_poroelastic)
      j = jvalue_inverse(ipoin1D,iedge_poroelastic)
      iglob = ibool(i,j,ispec_poroelastic)

      displ_x = displs_poroelastic(1,iglob)
      displ_z = displs_poroelastic(2,iglob)

      displw_x = displw_poroelastic(1,iglob)
      displw_z = displw_poroelastic(2,iglob)

      if (SIMULATION_TYPE == 3) then
        ! new definition of adjoint displacement and adjoint potential
        displ_x = accels_poroelastic_adj_coupling(1,iglob)
        displ_z = accels_poroelastic_adj_coupling(2,iglob)

        displw_x = accelw_poroelastic_adj_coupling(1,iglob)
        displw_z = accelw_poroelastic_adj_coupling(2,iglob)
      endif

      ! get point values for the acoustic side
      i = ivalue(ipoin1D,iedge_acoustic)
      j = jvalue(ipoin1D,iedge_acoustic)
      iglob = ibool(i,j,ispec_acoustic)

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

      ! compute dot product [u_s + u_w]*n
      displ_n = (displ_x + displw_x)*nx + (displ_z + displw_z)*nz
      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

      if (ATTENUATION_VISCOACOUSTIC .and. .not. USE_A_STRONG_FORMULATION_FOR_E1) &
          dot_e1(iglob,:) = dot_e1(iglob,:) + weight*displ_n

    enddo
  enddo

 end subroutine compute_coupling_acoustic_po

!========================================================================

! for acoustic solver

 subroutine compute_coupling_acoustic_po_backward()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP,ONE

  use specfem_par, only: num_fluid_poro_edges,ibool,wxgll,wzgll,xix,xiz, &
                         gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse, &
                         fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge, &
                         fluid_poro_poroelastic_ispec,fluid_poro_poroelastic_iedge, &
                         b_displs_poroelastic,b_displw_poroelastic, &
                         b_potential_dot_dot_acoustic

  implicit none

  !local variable
  integer :: inum,ispec_acoustic,iedge_acoustic,ispec_poroelastic,iedge_poroelastic, &
             ipoin1D,i,j,iglob
  real(kind=CUSTOM_REAL) :: b_displ_x,b_displ_z,b_displw_x,b_displw_z,b_displ_n, &
                            xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight


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
      ! get point values for the poroelastic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_poroelastic)
      j = jvalue_inverse(ipoin1D,iedge_poroelastic)
      iglob = ibool(i,j,ispec_poroelastic)

      b_displ_x = b_displs_poroelastic(1,iglob)
      b_displ_z = b_displs_poroelastic(2,iglob)

      b_displw_x = b_displw_poroelastic(1,iglob)
      b_displw_z = b_displw_poroelastic(2,iglob)

      ! get point values for the acoustic side
      i = ivalue(ipoin1D,iedge_acoustic)
      j = jvalue(ipoin1D,iedge_acoustic)
      iglob = ibool(i,j,ispec_acoustic)

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

      ! compute dot product [u_s + u_w]*n
      b_displ_n = (b_displ_x + b_displw_x)*nx + (b_displ_z + b_displw_z)*nz
      b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) + weight*b_displ_n
    enddo
  enddo

 end subroutine compute_coupling_acoustic_po_backward
