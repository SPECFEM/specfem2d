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
!========================================================================

! for viscoelastic solver

  subroutine compute_coupling_viscoelastic_ac()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,ZERO,ONE,TWO, &
    IRIGHT,ILEFT,IBOTTOM,ITOP,ALPHA_LDDRK,BETA_LDDRK

  use specfem_par, only: SIMULATION_TYPE,num_fluid_solid_edges,&
                         ibool,wxgll,wzgll,xix,xiz,gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse,&
                         minus_int_int_pressure_acoustic,minus_int_pressure_acoustic,minus_pressure_acoustic,&
                         accel_elastic,fluid_solid_acoustic_ispec, &
                         fluid_solid_acoustic_iedge,fluid_solid_elastic_ispec,fluid_solid_elastic_iedge,&
                         minus_int_int_pressure_acoustic_adj_coupling, &
                         AXISYM,coord,is_on_the_axis,xiglj,wxglj, &
                         rmemory_sfb_minus_pressure_acoustic,timeval,deltat,&
                         rmemory_sfb_minus_pressure_acoustic_LDDRK,i_stage,stage_time_scheme
  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,nspec_PML,ispec_is_PML,spec_to_PML,region_CPML, &
                K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store,minus_int_int_pressure_acoustic_old

  implicit none

  !local variable
  integer :: inum,ispec_acoustic,ispec_elastic,iedge_acoustic,iedge_elastic,ipoin1D,i,j,iglob,ii2,jj2,&
             ispec_PML,CPML_region_local,singularity_type
  real(kind=CUSTOM_REAL) :: pressure,xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1
  double precision :: kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z,&
                      A0,A1,A2,A3,A4,bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2

  ! loop on all the coupling edges
  do inum = 1,num_fluid_solid_edges

    ! get the edge of the acoustic element
    ispec_acoustic = fluid_solid_acoustic_ispec(inum)
    iedge_acoustic = fluid_solid_acoustic_iedge(inum)

    ! get the corresponding edge of the elastic element
    ispec_elastic = fluid_solid_elastic_ispec(inum)
    iedge_elastic = fluid_solid_elastic_iedge(inum)

    ! implement 1D coupling along the edge
    do ipoin1D = 1,NGLLX

      ! get point values for the acoustic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_acoustic)
      j = jvalue_inverse(ipoin1D,iedge_acoustic)
      iglob = ibool(i,j,ispec_acoustic)

      ! compute pressure on the fluid/solid edge
      pressure = - minus_pressure_acoustic(iglob)

      if (SIMULATION_TYPE == 3) then
        pressure = minus_int_int_pressure_acoustic_adj_coupling(iglob)
      endif

      ! PML
      ! (overwrites pressure value if needed)
      if (PML_BOUNDARY_CONDITIONS ) then
        if (ispec_is_PML(ispec_acoustic) .and. nspec_PML > 0) then
          ispec_PML = spec_to_PML(ispec_acoustic)
          CPML_region_local = region_CPML(ispec_acoustic)
          kappa_x = K_x_store(i,j,ispec_PML)
          kappa_z = K_z_store(i,j,ispec_PML)
          d_x = d_x_store(i,j,ispec_PML)
          d_z = d_z_store(i,j,ispec_PML)
          alpha_x = alpha_x_store(i,j,ispec_PML)
          alpha_z = alpha_z_store(i,j,ispec_PML)
          beta_x = alpha_x + d_x / kappa_x
          beta_z = alpha_z + d_z / kappa_z
          call l_parameter_computation(timeval,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                       CPML_region_local,A0,A1,A2,A3,A4,singularity_type,&
                                       bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2)

          if (stage_time_scheme == 1) then
            rmemory_sfb_minus_pressure_acoustic(1,i,j,inum) = &
                        coef0_1 * rmemory_sfb_minus_pressure_acoustic(1,i,j,inum) + &
                        coef1_1 * minus_int_int_pressure_acoustic(iglob) + coef2_1 * minus_int_int_pressure_acoustic_old(iglob)
          endif

          if (stage_time_scheme == 6) then
            rmemory_sfb_minus_pressure_acoustic_LDDRK(1,i,j,inum) = &
                    ALPHA_LDDRK(i_stage) * rmemory_sfb_minus_pressure_acoustic_LDDRK(1,i,j,inum) + &
                    deltat * (-bb_1 * rmemory_sfb_minus_pressure_acoustic(1,i,j,inum) + minus_int_int_pressure_acoustic(iglob))
            rmemory_sfb_minus_pressure_acoustic(1,i,j,inum) = rmemory_sfb_minus_pressure_acoustic(1,i,j,inum) + &
                    BETA_LDDRK(i_stage) * rmemory_sfb_minus_pressure_acoustic_LDDRK(1,i,j,inum)
          endif

          pressure = - (A0 * minus_pressure_acoustic(iglob) + A1 * minus_int_pressure_acoustic(iglob) + &
                        A2 * minus_int_int_pressure_acoustic(iglob) + A3 * rmemory_sfb_minus_pressure_acoustic(1,i,j,inum))
        endif
      endif

      ! get point values for the elastic side
      ii2 = ivalue(ipoin1D,iedge_elastic)
      jj2 = jvalue(ipoin1D,iedge_elastic)
      iglob = ibool(ii2,jj2,ispec_elastic)

      ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
      ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
      ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
      ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
      ! Blackwell Science, page 110, equation (4.60).

      if (AXISYM) then
! axial elements are always rotated by the mesher to make sure it is their i == 1 edge that is on the axis
! and thus the case of an edge located along j == constant does not need to be considered here
        if (is_on_the_axis(ispec_acoustic) .and. i == 1) then
          xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
          r_xiplus1(i,j) = xxi
        else if (is_on_the_axis(ispec_acoustic)) then
           r_xiplus1(i,j) = coord(1,ibool(i,j,ispec_acoustic))/(xiglj(i)+ONE)
        endif
      endif

      if (iedge_acoustic == ITOP ) then

        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
          else
            weight = jacobian1D * wxgll(i) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wxgll(i)
        endif

      else if (iedge_acoustic == IBOTTOM ) then

        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
          else
            weight = jacobian1D * wxgll(i) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wxgll(i)
        endif

      else if (iedge_acoustic == ILEFT ) then

        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            stop 'error: rotated element detected on the symmetry axis, this should not happen'
          else
            weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wzgll(j)
        endif

      else if (iedge_acoustic ==IRIGHT ) then

        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            stop 'error: rotated element detected on the symmetry axis, this should not happen'
          else
            weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wzgll(j)
        endif

      endif

      accel_elastic(1,iglob) = accel_elastic(1,iglob) + weight*nx*pressure
      accel_elastic(2,iglob) = accel_elastic(2,iglob) + weight*nz*pressure

    enddo
  enddo

  end subroutine compute_coupling_viscoelastic_ac

!========================================================================
! for viscoelastic solver

  subroutine compute_coupling_viscoelastic_ac_backward()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,ZERO,ONE,TWO,IRIGHT,ILEFT,IBOTTOM,ITOP

  use specfem_par, only: num_fluid_solid_edges,ibool,wxgll,wzgll,xix,xiz,gammax,gammaz, &
                         jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse, &
                         b_minus_pressure_acoustic,b_accel_elastic,fluid_solid_acoustic_ispec, &
                         fluid_solid_acoustic_iedge,fluid_solid_elastic_ispec,fluid_solid_elastic_iedge,&
                         AXISYM,coord,is_on_the_axis,xiglj,wxglj
  implicit none

  !local variable
  integer :: inum,ispec_acoustic,ispec_elastic,iedge_acoustic,iedge_elastic,ipoin1D,i,j,iglob,ii2,jj2
  real(kind=CUSTOM_REAL) :: b_pressure,xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  ! loop on all the coupling edges
  do inum = 1,num_fluid_solid_edges

    ! get the edge of the acoustic element
    ispec_acoustic = fluid_solid_acoustic_ispec(inum)
    iedge_acoustic = fluid_solid_acoustic_iedge(inum)

    ! get the corresponding edge of the elastic element
    ispec_elastic = fluid_solid_elastic_ispec(inum)
    iedge_elastic = fluid_solid_elastic_iedge(inum)

    ! implement 1D coupling along the edge
    do ipoin1D = 1,NGLLX

      ! get point values for the acoustic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_acoustic)
      j = jvalue_inverse(ipoin1D,iedge_acoustic)
      iglob = ibool(i,j,ispec_acoustic)

      b_pressure = - b_minus_pressure_acoustic(iglob)

      ! get point values for the elastic side
      ii2 = ivalue(ipoin1D,iedge_elastic)
      jj2 = jvalue(ipoin1D,iedge_elastic)
      iglob = ibool(ii2,jj2,ispec_elastic)

      ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
      ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
      ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
      ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
      ! Blackwell Science, page 110, equation (4.60).

      if (AXISYM) then
        if (is_on_the_axis(ispec_acoustic) .and. i == 1) then
          xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
          r_xiplus1(i,j) = xxi
        else if (is_on_the_axis(ispec_acoustic)) then
           r_xiplus1(i,j) = coord(1,ibool(i,j,ispec_acoustic))/(xiglj(i)+ONE)
        endif
      endif

      if (iedge_acoustic == ITOP ) then
        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
          else
            weight = jacobian1D * wxgll(i) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wxgll(i)
        endif
      else if (iedge_acoustic == IBOTTOM ) then
        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
          else
            weight = jacobian1D * wxgll(i) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wxgll(i)
        endif
      else if (iedge_acoustic == ILEFT ) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            stop 'error: rotated element detected on the symmetry axis, this should not happen'
          else
            weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wzgll(j)
        endif
      else if (iedge_acoustic ==IRIGHT ) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            stop 'error: rotated element detected on the symmetry axis, this should not happen'
          else
            weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wzgll(j)
        endif
      endif

      b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) + weight*nx*b_pressure
      b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) + weight*nz*b_pressure

    enddo
  enddo

  end subroutine compute_coupling_viscoelastic_ac_backward


