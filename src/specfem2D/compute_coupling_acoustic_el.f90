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

  subroutine compute_coupling_acoustic_el(displ_elastic,displ_elastic_old,potential_dot_dot_acoustic,dot_e1)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM, &
    CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP,ONE,ALPHA_LDDRK,BETA_LDDRK,USE_A_STRONG_FORMULATION_FOR_E1

  use specfem_par, only: num_fluid_solid_edges,ibool,wxgll,wzgll,xix,xiz, &
                         gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse, &
                         fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
                         fluid_solid_elastic_ispec,fluid_solid_elastic_iedge, &
                         AXISYM,coord,is_on_the_axis,xiglj,wxglj, &
                         rmemory_fsb_displ_elastic,deltat, &
                         rmemory_fsb_displ_elastic_LDDRK,i_stage,time_stepping_scheme, &
                         nglob_acoustic,nglob_elastic,iglob_is_forced, &
                         ATTENUATION_VISCOACOUSTIC,N_SLS,nglob_att

  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML,nspec_PML,spec_to_PML,region_CPML, &
                         K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  implicit none

  real(kind=CUSTOM_REAL),dimension(NDIM,nglob_elastic) :: displ_elastic,displ_elastic_old
  real(kind=CUSTOM_REAL),dimension(nglob_acoustic) :: potential_dot_dot_acoustic

  real(kind=CUSTOM_REAL),dimension(nglob_att,N_SLS) :: dot_e1

  !local variable
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1
  integer :: inum,ispec_acoustic,ispec_elastic,iedge_acoustic,iedge_elastic,ipoin1D,i,j,iglob
  real(kind=CUSTOM_REAL) :: displ_x,displ_z,displ_n, &
                            xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight
  ! PML
  integer :: ispec_PML,CPML_region_local
  double precision :: kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z, &
                      A8,A9,A10,bb_xz_1,bb_xz_2,coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2

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

      ! get point values for the elastic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_elastic)
      j = jvalue_inverse(ipoin1D,iedge_elastic)
      iglob = ibool(i,j,ispec_elastic)

      displ_x = displ_elastic(1,iglob)
      displ_z = displ_elastic(2,iglob)

      ! PML elements
      ! overwrites displ_x and displ_z
      if (PML_BOUNDARY_CONDITIONS) then
        if (ispec_is_PML(ispec_elastic) .and. nspec_PML > 0) then
          ispec_PML = spec_to_PML(ispec_elastic)
          CPML_region_local = region_CPML(ispec_elastic)
          kappa_x = K_x_store(i,j,ispec_PML)
          kappa_z = K_z_store(i,j,ispec_PML)
          d_x = d_x_store(i,j,ispec_PML)
          d_z = d_z_store(i,j,ispec_PML)
          alpha_x = alpha_x_store(i,j,ispec_PML)
          alpha_z = alpha_z_store(i,j,ispec_PML)
          beta_x = alpha_x + d_x / kappa_x
          beta_z = alpha_z + d_z / kappa_z

          if (CPML_region_local == CPML_X_ONLY) then
            !ZN needed to be change in case of viscoelastic medium
            call lik_parameter_computation(deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                           CPML_region_local,13,A8,A9,A10,bb_xz_1,bb_xz_2, &
                                           coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2)
          else if (CPML_region_local == CPML_Z_ONLY) then
            call lik_parameter_computation(deltat,kappa_z,beta_z,alpha_z,kappa_x,beta_x,alpha_x, &
                                           CPML_region_local,31,A8,A9,A10,bb_xz_1,bb_xz_2, &
                                           coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2)

          else
            call stop_the_code('PML do not support a fluid-solid boundary in corner PML region')
          endif

          if (time_stepping_scheme == 1) then
            ! Newmark
            rmemory_fsb_displ_elastic(1,1,i,j,inum) = coef0_xz_1 * rmemory_fsb_displ_elastic(1,1,i,j,inum) + &
                   coef1_xz_1 * displ_elastic(1,iglob) + coef2_xz_1 * displ_elastic_old(1,iglob)
            rmemory_fsb_displ_elastic(1,2,i,j,inum) = coef0_xz_1 * rmemory_fsb_displ_elastic(1,2,i,j,inum) + &
                   coef1_xz_1 * displ_elastic(2,iglob) + coef2_xz_1 * displ_elastic_old(2,iglob)

          else if (time_stepping_scheme == 2) then
            ! LDDRK
            rmemory_fsb_displ_elastic_LDDRK(1,1,i,j,inum) = &
                   ALPHA_LDDRK(i_stage) * rmemory_fsb_displ_elastic_LDDRK(1,1,i,j,inum) + &
                   deltat * ( - bb_xz_1 * rmemory_fsb_displ_elastic(1,1,i,j,inum) + displ_elastic(1,iglob) )
            rmemory_fsb_displ_elastic(1,1,i,j,inum) = rmemory_fsb_displ_elastic(1,1,i,j,inum) + &
                   BETA_LDDRK(i_stage) * rmemory_fsb_displ_elastic_LDDRK(1,1,i,j,inum)

            rmemory_fsb_displ_elastic_LDDRK(1,2,i,j,inum) = &
                  ALPHA_LDDRK(i_stage) * rmemory_fsb_displ_elastic_LDDRK(1,2,i,j,inum) + &
                  deltat * ( - bb_xz_1 * rmemory_fsb_displ_elastic(1,2,i,j,inum) + displ_elastic(2,iglob) )
            rmemory_fsb_displ_elastic(1,2,i,j,inum) = rmemory_fsb_displ_elastic(1,2,i,j,inum) + &
                   BETA_LDDRK(i_stage) * rmemory_fsb_displ_elastic_LDDRK(1,2,i,j,inum)
          endif

          displ_x = A8 * displ_elastic(1,iglob) + A9 * rmemory_fsb_displ_elastic(1,1,i,j,inum)
          displ_z = A8 * displ_elastic(2,iglob) + A9 * rmemory_fsb_displ_elastic(1,2,i,j,inum)
        endif
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

      if (AXISYM) then
! axial elements are always rotated by the mesher to make sure it is their i == 1 edge that is on the axis
! and thus the case of an edge located along j == constant does not need to be considered here
        if (is_on_the_axis(ispec_acoustic) .and. i == 1) then
          xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
          r_xiplus1(i,j) = xxi
        else if (is_on_the_axis(ispec_acoustic)) then
          r_xiplus1(i,j) = coord(1,iglob)/(xiglj(i)+ONE)
        endif
      endif

      if (iedge_acoustic == ITOP) then

        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
             weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
          else
             weight = jacobian1D * wxgll(i) * coord(1,iglob)
          endif
        else
          weight = jacobian1D * wxgll(i)
        endif

      else if (iedge_acoustic == IBOTTOM) then

        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
          else
            weight = jacobian1D * wxgll(i) * coord(1,iglob)
          endif
        else
          weight = jacobian1D * wxgll(i)
        endif

      else if (iedge_acoustic == ILEFT) then

        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            call stop_the_code('error: rotated element detected on the symmetry axis, this should not happen')
          else
            weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wzgll(j)
        endif

      else if (iedge_acoustic == IRIGHT) then

        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            call stop_the_code('error: rotated element detected on the symmetry axis, this should not happen')
          else
            weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wzgll(j)
        endif

      endif

      ! compute dot product
      displ_n = displ_x*nx + displ_z*nz
!! DK DK QUENTIN visco begin ici comme tu disais il faudrait coupler la composante tangentielle
!! DK DK QUENTIN en plus de la composante normale je suppose
      if (.not. iglob_is_forced(iglob)) then
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n
        if (ATTENUATION_VISCOACOUSTIC .and. .not. USE_A_STRONG_FORMULATION_FOR_E1 ) &
          dot_e1(iglob,:) = dot_e1(iglob,:) + weight*displ_n
      endif
!! DK DK QUENTIN visco end

    enddo
  enddo

  end subroutine compute_coupling_acoustic_el

!========================================================================
! for acoustic solver: backward simulation in adjoint inversion

  subroutine compute_coupling_acoustic_el_backward(b_displ_elastic,b_potential_dot_dot_acoustic)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM, &
    CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP,ONE

  use specfem_par, only: num_fluid_solid_edges,ibool,wxgll,wzgll,xix,xiz, &
                         gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse, &
                         fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
                         fluid_solid_elastic_ispec,fluid_solid_elastic_iedge, &
                         AXISYM,coord,is_on_the_axis,xiglj,wxglj, &
                         nglob_acoustic,nglob_elastic

  implicit none

  real(kind=CUSTOM_REAL),dimension(NDIM,nglob_elastic) :: b_displ_elastic
  real(kind=CUSTOM_REAL),dimension(nglob_acoustic) :: b_potential_dot_dot_acoustic

  !local variable
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1
  integer :: inum,ispec_acoustic,ispec_elastic,iedge_acoustic,iedge_elastic,ipoin1D,i,j,iglob
  real(kind=CUSTOM_REAL) :: displ_x,displ_z,displ_n, &
                            xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight

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

      ! get point values for the elastic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_elastic)
      j = jvalue_inverse(ipoin1D,iedge_elastic)
      iglob = ibool(i,j,ispec_elastic)


      displ_x = b_displ_elastic(1,iglob)
      displ_z = b_displ_elastic(2,iglob)

      ! get point values for the acoustic side
      i = ivalue(ipoin1D,iedge_acoustic)
      j = jvalue(ipoin1D,iedge_acoustic)
      iglob = ibool(i,j,ispec_acoustic)

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
          r_xiplus1(i,j) = coord(1,iglob)/(xiglj(i)+ONE)
        endif
      endif

      if (iedge_acoustic == ITOP) then
        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
             weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
          else
             weight = jacobian1D * wxgll(i) * coord(1,iglob)
          endif
        else
          weight = jacobian1D * wxgll(i)
        endif
      else if (iedge_acoustic == IBOTTOM) then
        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            weight = jacobian1D * wxglj(i) * r_xiplus1(i,j)
          else
            weight = jacobian1D * wxgll(i) * coord(1,iglob)
          endif
        else
          weight = jacobian1D * wxgll(i)
        endif
      else if (iedge_acoustic == ILEFT) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            call stop_the_code('error: rotated element detected on the symmetry axis, this should not happen')
          else
            weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wzgll(j)
        endif
      else if (iedge_acoustic == IRIGHT) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        if (AXISYM) then
          if (is_on_the_axis(ispec_acoustic)) then
            call stop_the_code('error: rotated element detected on the symmetry axis, this should not happen')
          else
            weight = jacobian1D * wzgll(j) * coord(1,ibool(i,j,ispec_acoustic))
          endif
        else
          weight = jacobian1D * wzgll(j)
        endif
      endif

      ! compute dot product
      displ_n = displ_x*nx + displ_z*nz
      b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) + weight*displ_n

    enddo
  enddo

  end subroutine compute_coupling_acoustic_el_backward
