
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
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
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

  subroutine invert_mass_matrix_init()

!  builds the global mass matrix

  use specfem_par, only: myrank,any_elastic,any_acoustic,any_gravitoacoustic,any_poroelastic, &
                                rmass_inverse_elastic_one, &
                                rmass_inverse_acoustic, &
                                rmass_inverse_gravitoacoustic, &
                                rmass_inverse_gravito, &
                                rmass_s_inverse_poroelastic, &
                                rmass_w_inverse_poroelastic, &
                                nspec,ibool,kmato,wxgll,wzgll,jacobian, &
                                elastic,acoustic,gravitoacoustic,poroelastic, &
                                assign_external_model, &
                                density,poroelastcoef,porosity,tortuosity, &
                                vpext,rhoext, &
                                anyabs,numabs,deltat,codeabs,codeabs_corner,&
                                ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
                                ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
                                rmass_inverse_elastic_three,&
                                nelemabs,vsext,xix,xiz,gammaz,gammax, &
                                K_x_store,K_z_store,is_PML,&
                                AXISYM,is_on_the_axis,coord,wxglj,xiglj, &
                                d_x_store,d_z_store,PML_BOUNDARY_CONDITIONS,region_CPML, &
                                spec_to_PML,time_stepping_scheme

  implicit none
  include 'constants.h'


  integer :: ibegin,iend,ispecabs,jbegin,jend,ispec,i,j,iglob


  ! local parameter
  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic,cpl,csl
  integer count_left,count_right,count_bottom
  real(kind=CUSTOM_REAL) :: nx,nz,vx,vy,vz,vn,rho_vp,rho_vs,tx,ty,tz,&
                            weight,xxi,zxi,xgamma,zgamma,jacobian1D
  double precision :: rhol,kappal,mul_relaxed,lambdal_relaxed
  double precision :: rhol_s,rhol_f,rhol_bar,phil,tortl
  integer :: ispec_PML
  logical :: this_element_has_PML

  if (myrank == 0) then
    write(IOUT,*) "  initializing mass matrices"
    call flush_IOUT()
  endif

  ! initialize mass matrix
  if(any_elastic) rmass_inverse_elastic_one(:) = 0._CUSTOM_REAL
  if(any_elastic) rmass_inverse_elastic_three(:) = 0._CUSTOM_REAL
  if(any_poroelastic) rmass_s_inverse_poroelastic(:) = 0._CUSTOM_REAL
  if(any_poroelastic) rmass_w_inverse_poroelastic(:) = 0._CUSTOM_REAL
  if(any_acoustic) rmass_inverse_acoustic(:) = 0._CUSTOM_REAL
  if (any_gravitoacoustic) then
      rmass_inverse_gravitoacoustic(:) = 0._CUSTOM_REAL
      rmass_inverse_gravito(:) = 0._CUSTOM_REAL
  endif

  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)

        ! if external density model (elastic or acoustic)
        if(assign_external_model) then
          rhol = rhoext(i,j,ispec)
          kappal = rhol * vpext(i,j,ispec)**2
        else
          rhol = density(1,kmato(ispec))
          lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
          mul_relaxed = poroelastcoef(2,1,kmato(ispec))
          kappal = lambdal_relaxed + 2.d0/3.d0*mul_relaxed
        endif

        if( poroelastic(ispec) ) then
          ! material is poroelastic

          rhol_s = density(1,kmato(ispec))
          rhol_f = density(2,kmato(ispec))
          phil = porosity(kmato(ispec))
          tortl = tortuosity(kmato(ispec))
          rhol_bar = (1.d0-phil)*rhol_s + phil*rhol_f

          ! for the solid mass matrix
          rmass_s_inverse_poroelastic(iglob) = rmass_s_inverse_poroelastic(iglob)  &
                  + wxgll(i)*wzgll(j)*jacobian(i,j,ispec)*(rhol_bar - phil*rhol_f/tortl)
          ! for the fluid mass matrix
          rmass_w_inverse_poroelastic(iglob) = rmass_w_inverse_poroelastic(iglob) &
                  + wxgll(i)*wzgll(j)*jacobian(i,j,ispec)*(rhol_bar*rhol_f*tortl  &
                  - phil*rhol_f*rhol_f)/(rhol_bar*phil)

          ! for elastic medium
        else if( elastic(ispec) ) then

          this_element_has_PML = .false.
          if(PML_BOUNDARY_CONDITIONS .and. size(is_PML) > 1) then
! do not merge this condition with the above line because array is_PML() sometimes has a dummy size of 1
            if (is_PML(ispec)) this_element_has_PML = .true.
          endif

          if(this_element_has_PML) then
            ispec_PML=spec_to_PML(ispec)
            if(time_stepping_scheme == 1)then

              if(region_CPML(ispec) == CPML_X_ONLY) then
                if (AXISYM) then  ! This PML can't be on the axis
                   rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                        + coord(1,iglob)*wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML) &
                        + d_x_store(i,j,ispec_PML) * deltat / 2.d0)
                 else ! not axisym
                   rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                        + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML)&
                        + d_x_store(i,j,ispec_PML) * deltat / 2.d0)
                 endif
                 rmass_inverse_elastic_three(iglob) = rmass_inverse_elastic_one(iglob)

              else if (region_CPML(ispec) == CPML_XZ_ONLY) then
                if (AXISYM) then  ! This corner can't be on the axis
                   rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                        + coord(1,iglob)*wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) &
                        * (K_x_store(i,j,ispec_PML) * K_z_store(i,j,ispec_PML) &
                        + (d_x_store(i,j,ispec_PML)*k_z_store(i,j,ispec_PML) + &
                        d_z_store(i,j,ispec_PML)*k_x_store(i,j,ispec_PML)) * deltat / 2.d0)
                 else ! not axisym
                   rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                        + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML) * K_z_store(i,j,ispec_PML)&
                        + (d_x_store(i,j,ispec_PML)*k_z_store(i,j,ispec_PML)+&
                          d_z_store(i,j,ispec_PML)*k_x_store(i,j,ispec_PML)) * deltat / 2.d0)
                 endif
                 rmass_inverse_elastic_three(iglob) = rmass_inverse_elastic_one(iglob)

              else if(region_CPML(ispec) == CPML_Z_ONLY) then
                if (AXISYM) then
                  if (is_on_the_axis(ispec)) then
                    if (is_on_the_axis(ispec) .and. i == 1) then ! First GLJ point
                      xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                      rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                         + xxi*wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec) &
                         * (K_z_store(i,j,ispec_PML) + d_z_store(i,j,ispec_PML)* deltat / 2.d0)
                    else
                      rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                         + coord(1,iglob)/(xiglj(i)+ONE)*wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec) &
                         * (K_z_store(i,j,ispec_PML) + d_z_store(i,j,ispec_PML)* deltat / 2.d0)
                    endif
                  else ! not on the axis
                    rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                         + coord(1,iglob)*wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) &
                         * (K_z_store(i,j,ispec_PML) + d_z_store(i,j,ispec_PML)* deltat / 2.d0)
                  endif
                else ! not axisym
                  rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                       + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_z_store(i,j,ispec_PML)&
                       + d_z_store(i,j,ispec_PML)* deltat / 2.d0)
                endif
                rmass_inverse_elastic_three(iglob) = rmass_inverse_elastic_one(iglob)
              endif

            else ! time_stepping_scheme /= 1
              if(region_CPML(ispec) == CPML_X_ONLY) then
                rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                     + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML))
                rmass_inverse_elastic_three(iglob) = rmass_inverse_elastic_one(iglob)
              else if (region_CPML(ispec) == CPML_XZ_ONLY) then
                rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                     + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML) * K_z_store(i,j,ispec_PML))
                rmass_inverse_elastic_three(iglob) = rmass_inverse_elastic_one(iglob)
              else if(region_CPML(ispec) == CPML_Z_ONLY) then
                rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                     + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * (K_z_store(i,j,ispec_PML))
                rmass_inverse_elastic_three(iglob) = rmass_inverse_elastic_one(iglob)
              endif
            endif

          else ! no PML

            if (AXISYM) then
              if (is_on_the_axis(ispec)) then
                if (is_on_the_axis(ispec) .and. i == 1) then ! First GLJ point
                  xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                  rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                      + xxi*wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec)
                else
                  rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                      + coord(1,iglob)/(xiglj(i)+ONE)*wxglj(i)*wzgll(j)*rhol*jacobian(i,j,ispec)
                endif
              else ! not on the axis
                rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                    + coord(1,iglob)*wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)
              endif
            else ! not axisym
              rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob)  &
                      + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)
            endif
            rmass_inverse_elastic_three(iglob) = rmass_inverse_elastic_one(iglob)
          endif

          ! for gravitoacoustic medium
          !!! PML NOT WORKING YET !!!
        else if( gravitoacoustic(ispec) ) then

        this_element_has_PML = .false.
        if(PML_BOUNDARY_CONDITIONS .and. size(is_PML) > 1) then
          if (is_PML(ispec)) stop 'PML not implemented yet for gravitoacoustic case'
        endif

          rmass_inverse_gravitoacoustic(iglob) = rmass_inverse_gravitoacoustic(iglob) &
                  + wxgll(i)*wzgll(j)*jacobian(i,j,ispec) / (kappal/rhol)
           rmass_inverse_gravito(iglob) = rmass_inverse_gravito(iglob) &
                  + wxgll(i)*wzgll(j)*jacobian(i,j,ispec)

         else
          ! for acoustic medium

          this_element_has_PML = .false.
          if(PML_BOUNDARY_CONDITIONS .and. size(is_PML) > 1) then
! do not merge this condition with the above line because array is_PML() sometimes has a dummy size of 1
            if (is_PML(ispec)) this_element_has_PML = .true.
          endif

          if(this_element_has_PML) then

            ispec_PML=spec_to_PML(ispec)
            if(time_stepping_scheme == 1)then
              if(region_CPML(ispec) == CPML_X_ONLY) then
                if (AXISYM) then   !! AB AB: This PML cannot be on the axis: it is a right PML
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob)  &
                       + coord(1,iglob)*wxgll(i)*wzgll(j)/ kappal*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML) &
                       + d_x_store(i,j,ispec_PML) * deltat / 2.d0)
                else ! not axisym
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob)  &
                       + wxgll(i)*wzgll(j)/ kappal*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML)&
                       + d_x_store(i,j,ispec_PML) * deltat / 2.d0)
                endif

              else if (region_CPML(ispec) == CPML_XZ_ONLY) then
                if (AXISYM) then   !! AB AB: This corner cannot be on the axis
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob)  &
                       + coord(1,iglob)*wxgll(i)*wzgll(j)/ kappal*jacobian(i,j,ispec) &
                       *  (K_x_store(i,j,ispec_PML) * K_z_store(i,j,ispec_PML) &
                       + (d_x_store(i,j,ispec_PML)*k_z_store(i,j,ispec_PML) &
                        + d_z_store(i,j,ispec_PML)*k_x_store(i,j,ispec_PML)) * deltat / 2.d0)
                else ! not axisym
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob)  &
                       + wxgll(i)*wzgll(j)/ kappal*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML) * K_z_store(i,j,ispec_PML) &
                       + (d_x_store(i,j,ispec_PML)*k_z_store(i,j,ispec_PML)&
                          + d_z_store(i,j,ispec_PML)*k_x_store(i,j,ispec_PML)) * deltat / 2.d0)
                endif

              else if(region_CPML(ispec) == CPML_Z_ONLY) then
                if (AXISYM) then
                  if (is_on_the_axis(ispec)) then
                    if (is_on_the_axis(ispec) .and. i == 1) then ! First GLJ point
                      xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                      rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob)  &
                         + xxi*wxglj(i)*wzgll(j)/kappal*jacobian(i,j,ispec) &
                         * (K_z_store(i,j,ispec_PML) + d_z_store(i,j,ispec_PML)* deltat / 2.d0)
                    else
                      rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob)  &
                         + coord(1,iglob)/(xiglj(i)+ONE)*wxglj(i)*wzgll(j)/kappal*jacobian(i,j,ispec) &
                         * (K_z_store(i,j,ispec_PML) + d_z_store(i,j,ispec_PML)* deltat / 2.d0)
                    endif
                  else ! not on the axis
                    rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob)  &
                         + coord(1,iglob)*wxgll(i)*wzgll(j)/kappal*jacobian(i,j,ispec) &
                         * (K_z_store(i,j,ispec_PML) + d_z_store(i,j,ispec_PML)* deltat / 2.d0)
                  endif
                else ! not axisym
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob)  &
                       + wxgll(i)*wzgll(j)/ kappal*jacobian(i,j,ispec) * (K_z_store(i,j,ispec_PML)&
                       + d_z_store(i,j,ispec_PML)* deltat / 2.d0)
                endif
              endif
            else  ! time_stepping_scheme /= 1
              if(region_CPML(ispec) == CPML_X_ONLY) then
                rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob)  &
                     + wxgll(i)*wzgll(j)/ kappal*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML))
              else if (region_CPML(ispec) == CPML_XZ_ONLY) then
                rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob)  &
                     + wxgll(i)*wzgll(j)/ kappal*jacobian(i,j,ispec) * (K_x_store(i,j,ispec_PML) * K_z_store(i,j,ispec_PML))
              else if(region_CPML(ispec) == CPML_Z_ONLY) then
                rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob)  &
                     + wxgll(i)*wzgll(j)/ kappal*jacobian(i,j,ispec) * (K_z_store(i,j,ispec_PML))
              endif
            endif
          else  ! no PML

            if (AXISYM) then
              if (is_on_the_axis(ispec)) then
                if (is_on_the_axis(ispec) .and. i == 1) then ! First GLJ point
                  xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                      + xxi*wxglj(i)*wzgll(j)*jacobian(i,j,ispec) / kappal
                else
                  rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                      + coord(1,iglob)/(xiglj(i)+ONE)*wxglj(i)*wzgll(j)*jacobian(i,j,ispec) / kappal
                endif
              else ! not on the axis
                rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                    + coord(1,iglob)*wxgll(i)*wzgll(j)*jacobian(i,j,ispec) / kappal
              endif
            else ! not axisym
              rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                   + wxgll(i)*wzgll(j)*jacobian(i,j,ispec) / kappal
            endif

          endif
        endif
      enddo
    enddo
  enddo ! of do ispec = 1,nspec

  !
  !--- DK and Zhinan Xie: add C Delta_t / 2 contribution to the mass matrix
  !--- DK and Zhinan Xie: in the case of Clayton-Engquist absorbing boundaries;
  !--- DK and Zhinan Xie: see for instance the book of Hughes (1987) chapter 9.
  !--- DK and Zhinan Xie: IMPORTANT: note that this implies that we must have two different mass matrices,
  !--- DK and Zhinan Xie: one per component of the wave field i.e. one per spatial dimension.
  !--- DK and Zhinan Xie: This was also suggested by Jean-Paul Ampuero in 2003.
  !

  if(.not. PML_BOUNDARY_CONDITIONS .and. anyabs .and. time_stepping_scheme == 1) then
     count_left=1
     count_right=1
     count_bottom=1
     do ispecabs = 1,nelemabs
        ispec = numabs(ispecabs)
        ! get elastic parameters of current spectral elemegammaznt
        lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
        mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
        rhol  = density(1,kmato(ispec))
        kappal  = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic/3._CUSTOM_REAL
        cpl = sqrt((kappal + 4._CUSTOM_REAL*mul_unrelaxed_elastic/3._CUSTOM_REAL)/rhol)
        csl = sqrt(mul_unrelaxed_elastic/rhol)

        !--- left absorbing boundary
        if(codeabs(IEDGE4,ispecabs)) then

           i = 1
           do j = 1,NGLLZ
              iglob = ibool(i,j,ispec)
              ! external velocity model
              if(assign_external_model) then
                 cpl = vpext(i,j,ispec)
                 csl = vsext(i,j,ispec)
                 rhol = rhoext(i,j,ispec)
              endif

              rho_vp = rhol*cpl
              rho_vs = rhol*csl

              xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
              zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
              jacobian1D = sqrt(xgamma**2 + zgamma**2)
              nx = - zgamma / jacobian1D
              nz = + xgamma / jacobian1D

              weight = jacobian1D * wzgll(j)

              ! Clayton-Engquist condition if elastic
              if(elastic(ispec)) then

                 vx = 1.0d0*deltat/2.0d0
                 vy = 1.0d0*deltat/2.0d0
                 vz = 1.0d0*deltat/2.0d0

                 vn = nx*vx+nz*vz

                 tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
                 ty = rho_vs*vy
                 tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

                rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob) + tx*weight
                rmass_inverse_elastic_three(iglob) = rmass_inverse_elastic_three(iglob) + tz*weight

              endif
           enddo

        endif  !  end of left absorbing boundary

        !--- right absorbing boundary
        if(codeabs(IEDGE2,ispecabs)) then

           i = NGLLX

           do j = 1,NGLLZ

              iglob = ibool(i,j,ispec)

              ! for analytical initial plane wave for Bielak's conditions
              ! left or right edge, horizontal normal vector

              ! external velocity model
              if(assign_external_model) then
                 cpl = vpext(i,j,ispec)
                 csl = vsext(i,j,ispec)
                 rhol = rhoext(i,j,ispec)
              endif

              rho_vp = rhol*cpl
              rho_vs = rhol*csl

              xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
              zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
              jacobian1D = sqrt(xgamma**2 + zgamma**2)
              nx = + zgamma / jacobian1D
              nz = - xgamma / jacobian1D

              weight = jacobian1D * wzgll(j)

              ! Clayton-Engquist condition if elastic
              if(elastic(ispec)) then

                 vx = 1.0d0*deltat/2.0d0
                 vy = 1.0d0*deltat/2.0d0
                 vz = 1.0d0*deltat/2.0d0

                 vn = nx*vx+nz*vz

                 tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
                 ty = rho_vs*vy
                 tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

                rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob) + tx*weight
                rmass_inverse_elastic_three(iglob) = rmass_inverse_elastic_three(iglob) + tz*weight

              endif

           enddo

        endif  !  end of right absorbing boundary

        !--- bottom absorbing boundary
        if(codeabs(IEDGE1,ispecabs)) then

           j = 1
           ibegin = 1
           iend = NGLLX

           do i = ibegin,iend

              iglob = ibool(i,j,ispec)
              ! external velocity model
              if(assign_external_model) then
                 cpl = vpext(i,j,ispec)
                 csl = vsext(i,j,ispec)
                 rhol = rhoext(i,j,ispec)
              endif

              rho_vp = rhol*cpl
              rho_vs = rhol*csl

              xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
              zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
              jacobian1D = sqrt(xxi**2 + zxi**2)
              nx = + zxi / jacobian1D
              nz = - xxi / jacobian1D

              weight = jacobian1D * wxgll(i)

              ! Clayton-Engquist condition if elastic
              if(elastic(ispec)) then

                 vx = 1.0d0*deltat/2.0d0
                 vy = 1.0d0*deltat/2.0d0
                 vz = 1.0d0*deltat/2.0d0

                 vn = nx*vx+nz*vz

                 tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
                 ty = rho_vs*vy
                 tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
                 if((codeabs_corner(1,ispecabs) .and. i == 1) .or. (codeabs_corner(2,ispecabs) .and. i == NGLLX)) then
                   tx = 0
                   ty = 0
                   tz = 0
                 else
                   rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob) + tx*weight
                   rmass_inverse_elastic_three(iglob) = rmass_inverse_elastic_three(iglob) + tz*weight
                 endif

             endif

           enddo

        endif  !  end of bottom absorbing boundary

        !--- top absorbing boundary
        if(codeabs(IEDGE3,ispecabs)) then

           j = NGLLZ

           ibegin = 1
           iend = NGLLX

           do i = ibegin,iend

              iglob = ibool(i,j,ispec)

              ! external velocity model
              if(assign_external_model) then
                 cpl = vpext(i,j,ispec)
                 csl = vsext(i,j,ispec)
                 rhol = rhoext(i,j,ispec)
              endif

              rho_vp = rhol*cpl
              rho_vs = rhol*csl

              xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
              zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
              jacobian1D = sqrt(xxi**2 + zxi**2)
              nx = - zxi / jacobian1D
              nz = + xxi / jacobian1D

              weight = jacobian1D * wxgll(i)

              ! Clayton-Engquist condition if elastic
              if(elastic(ispec)) then

                 vx = 1.0d0*deltat/2.0d0
                 vy = 1.0d0*deltat/2.0d0
                 vz = 1.0d0*deltat/2.0d0

                 vn = nx*vx+nz*vz

                 tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
                 ty = rho_vs*vy
                 tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
                 if((codeabs_corner(3,ispecabs) .and. i == 1) .or. (codeabs_corner(4,ispecabs) .and. i == NGLLX)) then
                   tx = 0
                   ty = 0
                   tz = 0
                 else
                   rmass_inverse_elastic_one(iglob) = rmass_inverse_elastic_one(iglob) + tx*weight
                   rmass_inverse_elastic_three(iglob) = rmass_inverse_elastic_three(iglob) + tz*weight
                 endif
            endif


           enddo
        endif  !  end of top absorbing boundary
     enddo
  endif  ! end of absorbing boundaries

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(.not. PML_BOUNDARY_CONDITIONS .and. anyabs .and. time_stepping_scheme /= 1) then

    do ispecabs=1,nelemabs

      ispec = numabs(ispecabs)

      ! Sommerfeld condition if acoustic
      if(acoustic(ispec)) then

        ! get elastic parameters of current spectral element
        lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
        mul_relaxed = poroelastcoef(2,1,kmato(ispec))
        kappal  = lambdal_relaxed + TWO*mul_relaxed/3._CUSTOM_REAL
        rhol = density(1,kmato(ispec))

        cpl = sqrt(kappal/rhol)

        !--- left absorbing boundary
        if(codeabs(IEDGE4,ispecabs)) then
          i = 1
          jbegin = ibegin_edge4(ispecabs)
          jend = iend_edge4(ispecabs)
          do j = jbegin,jend
            iglob = ibool(i,j,ispec)
            ! external velocity model
            if(assign_external_model) then
              cpl = vpext(i,j,ispec)
              rhol = rhoext(i,j,ispec)
            endif
            xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
            zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            weight = jacobian1D * wzgll(j)

            rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                  + 1.0d0*deltat/2.0d0*weight/cpl/rhol

          enddo

        endif  !  end of left absorbing boundary

        !--- right absorbing boundary
        if(codeabs(IEDGE2,ispecabs)) then
          i = NGLLX
          jbegin = ibegin_edge2(ispecabs)
          jend = iend_edge2(ispecabs)
          do j = jbegin,jend
            iglob = ibool(i,j,ispec)
            ! external velocity model
            if(assign_external_model) then
              cpl = vpext(i,j,ispec)
              rhol = rhoext(i,j,ispec)
            endif
            xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
            zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            weight = jacobian1D * wzgll(j)

            rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                  + 1.0d0*deltat/2.0d0*weight/cpl/rhol
          enddo
        endif  !  end of right absorbing boundary

        !--- bottom absorbing boundary
        if(codeabs(IEDGE1,ispecabs)) then
          j = 1
          ibegin = ibegin_edge1(ispecabs)
          iend = iend_edge1(ispecabs)
          ! exclude corners to make sure there is no contradiction on the normal
          if(codeabs_corner(1,ispecabs)) ibegin = 2
          if(codeabs_corner(2,ispecabs)) iend = NGLLX-1
          do i = ibegin,iend
            iglob = ibool(i,j,ispec)
            ! external velocity model
            if(assign_external_model) then
              cpl = vpext(i,j,ispec)
              rhol = rhoext(i,j,ispec)
            endif
            xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
            zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            weight = jacobian1D * wxgll(i)

            rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                 + 1.0d0*deltat/2.0d0*weight/cpl/rhol

          enddo
        endif  !  end of bottom absorbing boundary

        !--- top absorbing boundary
        if(codeabs(IEDGE3,ispecabs)) then
          j = NGLLZ
          ibegin = ibegin_edge3(ispecabs)
          iend = iend_edge3(ispecabs)
          ! exclude corners to make sure there is no contradiction on the normal
          if(codeabs_corner(3,ispecabs)) ibegin = 2
          if(codeabs_corner(4,ispecabs)) iend = NGLLX-1
          do i = ibegin,iend
            iglob = ibool(i,j,ispec)
            ! external velocity model
            if(assign_external_model) then
              cpl = vpext(i,j,ispec)
              rhol = rhoext(i,j,ispec)
            endif
            xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
            zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            weight = jacobian1D * wxgll(i)

            rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                 + 1.0d0*deltat/2.0d0*weight/cpl/rhol

          enddo
        endif  !  end of top absorbing boundary

      endif ! acoustic ispec
    enddo
  endif  ! end of absorbing boundaries

  end subroutine invert_mass_matrix_init
!
!-------------------------------------------------------------------------------------------------
!

  subroutine invert_mass_matrix()

! inverts the global mass matrix

  use specfem_par, only: myrank,any_elastic,any_acoustic,any_gravitoacoustic,any_poroelastic, &
                                rmass_inverse_elastic_one,rmass_inverse_elastic_three,&
                                rmass_inverse_acoustic, &
                                rmass_inverse_gravitoacoustic, &
                                rmass_inverse_gravito, &
                                rmass_s_inverse_poroelastic, &
                                rmass_w_inverse_poroelastic
  implicit none
  include 'constants.h'

  if (myrank == 0) then
    write(IOUT,*) "  inverting mass matrices"
    call flush_IOUT()
  endif

! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
! (this can happen when some degrees of freedom have been removed from some of the global arrays)
  if(any_elastic) &
    where(rmass_inverse_elastic_one <= 0._CUSTOM_REAL) rmass_inverse_elastic_one = 1._CUSTOM_REAL
  if(any_elastic) &
    where(rmass_inverse_elastic_three <= 0._CUSTOM_REAL) rmass_inverse_elastic_three = 1._CUSTOM_REAL
  if(any_poroelastic) &
    where(rmass_s_inverse_poroelastic <= 0._CUSTOM_REAL) rmass_s_inverse_poroelastic = 1._CUSTOM_REAL
  if(any_poroelastic) &
    where(rmass_w_inverse_poroelastic <= 0._CUSTOM_REAL) rmass_w_inverse_poroelastic = 1._CUSTOM_REAL
  if(any_acoustic) &
    where(rmass_inverse_acoustic <= 0._CUSTOM_REAL) rmass_inverse_acoustic = 1._CUSTOM_REAL
  if(any_gravitoacoustic) then
    where(rmass_inverse_gravitoacoustic <= 0._CUSTOM_REAL) rmass_inverse_gravitoacoustic = 1._CUSTOM_REAL
    where(rmass_inverse_gravito <= 0._CUSTOM_REAL) rmass_inverse_gravito = 1._CUSTOM_REAL
  endif

! compute the inverse of the mass matrix
  if(any_elastic) &
    rmass_inverse_elastic_one(:) = 1._CUSTOM_REAL / rmass_inverse_elastic_one(:)
  if(any_elastic) &
    rmass_inverse_elastic_three(:) = 1._CUSTOM_REAL / rmass_inverse_elastic_three(:)
  if(any_poroelastic) &
    rmass_s_inverse_poroelastic(:) = 1._CUSTOM_REAL / rmass_s_inverse_poroelastic(:)
  if(any_poroelastic) &
    rmass_w_inverse_poroelastic(:) = 1._CUSTOM_REAL / rmass_w_inverse_poroelastic(:)
  if(any_acoustic) &
    rmass_inverse_acoustic(:) = 1._CUSTOM_REAL / rmass_inverse_acoustic(:)
  if(any_gravitoacoustic) then
    rmass_inverse_gravitoacoustic(:) = 1._CUSTOM_REAL / rmass_inverse_gravitoacoustic(:)
    rmass_inverse_gravito(:) = 1._CUSTOM_REAL / rmass_inverse_gravito(:)
  endif

  end subroutine invert_mass_matrix

