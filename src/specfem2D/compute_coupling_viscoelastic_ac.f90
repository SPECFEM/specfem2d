!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, INRIA and University of Pau, France,
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
! for viscoelastic solver

  subroutine compute_coupling_viscoelastic_ac()

  use specfem_par, only: SIMULATION_TYPE,num_fluid_solid_edges,&
                         ibool,wxgll,wzgll,xix,xiz,gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse,&
                         potential_acoustic,potential_acoustic_old,potential_dot_acoustic,potential_dot_dot_acoustic,&
                         accel_elastic,fluid_solid_acoustic_ispec, &
                         fluid_solid_acoustic_iedge,fluid_solid_elastic_ispec,fluid_solid_elastic_iedge,&
                         potential_acoustic_adj_coupling,AXISYM,coord,is_on_the_axis,xiglj,wxglj, &
                         PML_BOUNDARY_CONDITIONS,nspec_PML,K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,&
                         alpha_z_store,is_PML,spec_to_PML,region_CPML,rmemory_sfb_potential_ddot_acoustic,timeval,deltat,&
                         rmemory_sfb_potential_ddot_acoustic_LDDRK,i_stage,stage_time_scheme,alpha_LDDRK,beta_LDDRK
  implicit none
  include 'constants.h'

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
      pressure = - potential_dot_dot_acoustic(iglob)

      if( PML_BOUNDARY_CONDITIONS  ) then
        if( is_PML(ispec_acoustic) .and. nspec_PML > 0 ) then
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

          if( stage_time_scheme == 1 ) then
            rmemory_sfb_potential_ddot_acoustic(1,i,j,inum) = &
                        coef0_1 * rmemory_sfb_potential_ddot_acoustic(1,i,j,inum) + &
                        coef1_1 * potential_acoustic(iglob) + coef2_1 * potential_acoustic_old(iglob)
          endif

          if( stage_time_scheme == 6 ) then
            rmemory_sfb_potential_ddot_acoustic_LDDRK(1,i,j,inum) = &
                    alpha_LDDRK(i_stage) * rmemory_sfb_potential_ddot_acoustic_LDDRK(1,i,j,inum) + &
                    deltat * (-bb_1 * rmemory_sfb_potential_ddot_acoustic(1,i,j,inum) + potential_acoustic(iglob))
            rmemory_sfb_potential_ddot_acoustic(1,i,j,inum) = rmemory_sfb_potential_ddot_acoustic(1,i,j,inum) + &
                    beta_LDDRK(i_stage) * rmemory_sfb_potential_ddot_acoustic_LDDRK(1,i,j,inum)
          endif

          pressure = - (A0 * potential_dot_dot_acoustic(iglob) + A1 * potential_dot_acoustic(iglob) + &
                        A2 * potential_acoustic(iglob) + A3 * rmemory_sfb_potential_ddot_acoustic(1,i,j,inum))
        else
          pressure = - potential_dot_dot_acoustic(iglob)
        endif
      else
        pressure = - potential_dot_dot_acoustic(iglob)
      endif

      if( SIMULATION_TYPE == 3 ) then
        !<YANGL
        ! new definition of adjoint displacement and adjoint potential
        pressure = potential_acoustic_adj_coupling(iglob)
        !>YANGL
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

      if( AXISYM ) then
! axial elements are always rotated by the mesher to make sure it is their i == 1 edge that is on the axis
! and thus the case of an edge located along j == constant does not need to be considered here
        if( is_on_the_axis(ispec_acoustic) .and. i == 1 ) then
          xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
          r_xiplus1(i,j) = xxi
        else if(is_on_the_axis(ispec_acoustic) ) then
           r_xiplus1(i,j) = coord(1,ibool(i,j,ispec_acoustic))/(xiglj(i)+ONE)
        endif
      endif

      if( iedge_acoustic == ITOP  ) then

        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        if( AXISYM ) then
          if( is_on_the_axis(ispec_acoustic) ) then
            weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
          else
            weight = jacobian1D * wxgll(i) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wxgll(i)
        endif

      else if( iedge_acoustic == IBOTTOM  ) then

        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        if( AXISYM ) then
          if( is_on_the_axis(ispec_acoustic) ) then
            weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
          else
            weight = jacobian1D * wxgll(i) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wxgll(i)
        endif

      else if( iedge_acoustic == ILEFT  ) then

        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        if( AXISYM ) then
          if( is_on_the_axis(ispec_acoustic) ) then
            stop 'error: rotated element detected on the symmetry axis, this should not happen'
          else
            weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wzgll(j)
        endif

      else if( iedge_acoustic ==IRIGHT  ) then

        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        if( AXISYM ) then
          if( is_on_the_axis(ispec_acoustic) ) then
            stop 'error: rotated element detected on the symmetry axis, this should not happen'
          else
            weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wzgll(j)
        endif

      endif

      accel_elastic(1,iglob) = accel_elastic(1,iglob) + weight*nx*pressure
      accel_elastic(3,iglob) = accel_elastic(3,iglob) + weight*nz*pressure

    enddo
  enddo

  end subroutine compute_coupling_viscoelastic_ac

!========================================================================
! for viscoelastic solver

  subroutine compute_coupling_viscoelastic_ac_backward()

  use specfem_par, only: num_fluid_solid_edges,ibool,wxgll,wzgll,xix,xiz,gammax,gammaz, &
                         jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse, &
                         b_potential_dot_dot_acoustic,b_accel_elastic,fluid_solid_acoustic_ispec, &
                         fluid_solid_acoustic_iedge,fluid_solid_elastic_ispec,fluid_solid_elastic_iedge,&
                         AXISYM,coord,is_on_the_axis,xiglj,wxglj
  implicit none
  include 'constants.h'

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

      b_pressure = - b_potential_dot_dot_acoustic(iglob)

      ! get point values for the elastic side
      ii2 = ivalue(ipoin1D,iedge_elastic)
      jj2 = jvalue(ipoin1D,iedge_elastic)
      iglob = ibool(ii2,jj2,ispec_elastic)

      ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
      ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
      ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
      ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
      ! Blackwell Science, page 110, equation (4.60).

      if( AXISYM ) then
        if( is_on_the_axis(ispec_acoustic) .and. i == 1 ) then
          xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
          r_xiplus1(i,j) = xxi
        else if(is_on_the_axis(ispec_acoustic) ) then
           r_xiplus1(i,j) = coord(1,ibool(i,j,ispec_acoustic))/(xiglj(i)+ONE)
        endif
      endif

      if( iedge_acoustic == ITOP  ) then
        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        if( AXISYM ) then
          if( is_on_the_axis(ispec_acoustic) ) then
            weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
          else
            weight = jacobian1D * wxgll(i) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wxgll(i)
        endif
      else if( iedge_acoustic == IBOTTOM  ) then
        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        if( AXISYM ) then
          if( is_on_the_axis(ispec_acoustic) ) then
            weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
          else
            weight = jacobian1D * wxgll(i) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wxgll(i)
        endif
      else if( iedge_acoustic == ILEFT  ) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        if( AXISYM ) then
          if( is_on_the_axis(ispec_acoustic) ) then
            stop 'error: rotated element detected on the symmetry axis, this should not happen'
          else
            weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wzgll(j)
        endif
      else if( iedge_acoustic ==IRIGHT  ) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        if( AXISYM ) then
          if( is_on_the_axis(ispec_acoustic) ) then
            stop 'error: rotated element detected on the symmetry axis, this should not happen'
          else
            weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wzgll(j)
        endif
      endif

      b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) + weight*nx*b_pressure
      b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) + weight*nz*b_pressure

    enddo
  enddo

  end subroutine compute_coupling_viscoelastic_ac_backward


