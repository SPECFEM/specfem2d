!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, Inria and University of Pau, France,
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
! license as circulated by CEA, CNRS and Inria at the following URL
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

! for acoustic solver

  subroutine compute_coupling_acoustic_el(nspec,nglob_elastic,nglob_acoustic,num_fluid_solid_edges,ibool,wxgll,wzgll,xix,xiz,&
                              gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse,displ_elastic,displ_elastic_old,&
                              potential_dot_dot_acoustic,fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
                              fluid_solid_elastic_ispec,fluid_solid_elastic_iedge,&
                              PML_BOUNDARY_CONDITIONS,nspec_PML,K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,&
                              alpha_z_store,is_PML,spec_to_PML,region_CPML,rmemory_fsb_displ_elastic,time,deltat)

   implicit none
   include 'constants.h'

   integer :: nspec,nglob_elastic,nglob_acoustic,num_fluid_solid_edges

   integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
   real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll,wzgll
   real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec)  :: xix,xiz,gammax,gammaz,jacobian
   integer, dimension(NGLLX,NEDGES) :: ivalue,jvalue,ivalue_inverse,jvalue_inverse

   real(kind=CUSTOM_REAL),dimension(3,nglob_elastic) :: displ_elastic,displ_elastic_old
   real(kind=CUSTOM_REAL),dimension(nglob_acoustic) :: potential_dot_dot_acoustic
   integer, dimension(num_fluid_solid_edges) :: fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
                                                fluid_solid_elastic_ispec,fluid_solid_elastic_iedge
   integer :: nspec_PML
   real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_PML) :: &
                  K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store
   logical:: PML_BOUNDARY_CONDITIONS
   logical, dimension(nspec) :: is_PML
   integer, dimension(nspec) :: spec_to_PML
   integer, dimension(nspec) :: region_CPML
   real(kind=CUSTOM_REAL),dimension(1,3,NGLLX,NGLLZ,num_fluid_solid_edges) :: rmemory_fsb_displ_elastic

!local variable

   integer :: inum,ispec_acoustic,ispec_elastic,iedge_acoustic,iedge_elastic,ipoin1D,i,j,iglob,&
              ispec_PML,CPML_region_local,singularity_type_xz
   real(kind=CUSTOM_REAL) :: displ_x,displ_z,displ_n,&
                             xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight
   real(kind=CUSTOM_REAL) :: time,deltat,kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z,&
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

          if(PML_BOUNDARY_CONDITIONS)then
             if(is_PML(ispec_elastic) .and. nspec_PML > 0) then
               ispec_PML = spec_to_PML(ispec_elastic)
               CPML_region_local = region_CPML(ispec_elastic)
               if(CPML_region_local == CPML_X_ONLY)then
                  kappa_x = K_x_store(i,j,ispec_PML)
                  kappa_z = K_z_store(i,j,ispec_PML)
                  d_x = d_x_store(i,j,ispec_PML)
                  d_z = d_z_store(i,j,ispec_PML)
                  alpha_x = alpha_x_store(i,j,ispec_PML)
                  alpha_z = alpha_z_store(i,j,ispec_PML)
                  beta_x = alpha_x + d_x / kappa_x
                  beta_z = alpha_z + d_z / kappa_z
                  call lik_parameter_computation(time,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z,&
                                           CPML_region_local,13,A8,A9,A10,singularity_type_xz,bb_xz_1,bb_xz_2,&
                                           coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2)
                  rmemory_fsb_displ_elastic(1,1,i,j,inum) = coef0_xz_1 * rmemory_fsb_displ_elastic(1,1,i,j,inum) + &
                                coef1_xz_1 * displ_elastic(1,iglob) + coef2_xz_1 * displ_elastic_old(1,iglob)
                  rmemory_fsb_displ_elastic(1,3,i,j,inum) = coef0_xz_1 * rmemory_fsb_displ_elastic(1,3,i,j,inum) + &
                                coef1_xz_1 * displ_elastic(3,iglob) + coef2_xz_1 * displ_elastic_old(3,iglob)
                  displ_x = A8 * displ_elastic(1,iglob) + A9 * rmemory_fsb_displ_elastic(1,1,i,j,inum)
                  displ_z = A8 * displ_elastic(3,iglob) + A9 * rmemory_fsb_displ_elastic(1,3,i,j,inum)
               else
                  stop 'PML do not support the fluid-solid boundary not inside CPML_X_ONLY'
               endif
             else
               displ_x = displ_elastic(1,iglob)
               displ_z = displ_elastic(3,iglob)
             endif
          else
            displ_x = displ_elastic(1,iglob)
            displ_z = displ_elastic(3,iglob)
          endif

!          if(SIMULATION_TYPE == 3) then
!            b_displ_x = b_displ_elastic(1,iglob)
!            b_displ_z = b_displ_elastic(3,iglob)
            !<YANGL
            ! new definition of adjoint displacement and adjoint potential
!            displ_x = -accel_elastic_adj_coupling(1,iglob)
!            displ_z = -accel_elastic_adj_coupling(3,iglob)
            !>YANGL
!          endif

          ! get point values for the acoustic side
          i = ivalue(ipoin1D,iedge_acoustic)
          j = jvalue(ipoin1D,iedge_acoustic)
          iglob = ibool(i,j,ispec_acoustic)

          ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
          ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
          ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
          ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
          ! Blackwell Science, page 110, equation (4.60).
          if(iedge_acoustic == ITOP)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = - zxi / jacobian1D
            nz = + xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_acoustic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
            weight = jacobian1D * wxgll(i)
          else if(iedge_acoustic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = - zgamma / jacobian1D
            nz = + xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          else if(iedge_acoustic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
            weight = jacobian1D * wzgll(j)
          endif

          ! compute dot product
          displ_n = displ_x*nx + displ_z*nz

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

!          if(SIMULATION_TYPE == 3) then
!          b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) + &
!                      weight*(b_displ_x*nx + b_displ_z*nz)
!          endif !if(SIMULATION_TYPE == 3) then

        enddo

      enddo

  end subroutine compute_coupling_acoustic_el


