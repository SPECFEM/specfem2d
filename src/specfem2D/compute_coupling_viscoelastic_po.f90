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

 subroutine compute_coupling_viscoelastic_po()

  use specfem_par, only: SIMULATION_TYPE,num_solid_poro_edges,&
                         ibool,wxgll,wzgll,xix,xiz,gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse, &
                         hprime_xx,hprime_zz, &
                         solid_poro_elastic_ispec,solid_poro_elastic_iedge, &
                         solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge,&
                         kmato,porosity,tortuosity,poroelastcoef,density, &
                         assign_external_model,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,anisotropy, &
                         displ_elastic,b_displ_elastic,displs_poroelastic,displw_poroelastic, &
                         b_displs_poroelastic,b_displw_poroelastic,&
                         accel_elastic,b_accel_elastic

  implicit none
  include 'constants.h'

  !local variables
  integer :: inum,ispec_elastic,iedge_elastic,ispec_poroelastic,iedge_poroelastic, &
             i,j,k,ipoin1D,iglob,ii2,jj2
  double precision :: phil,tortl,mul_s,kappal_s,rhol_s,kappal_f,rhol_f, &
                      mul_fr,kappal_fr,rhol_bar,D_biot,H_biot,C_biot,M_biot, &
                      mul_G,lambdal_G,lambdalplus2mul_G,lambdal_unrelaxed_elastic, &
                      mul_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic, &
                      c11,c13,c15,c33,c35,c55,c12,c23,c25

  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: dwx_dxl,dwz_dxl,dwx_dzl,dwz_dzl
  real(kind=CUSTOM_REAL) :: b_dux_dxi,b_dux_dgamma,b_duz_dxi,b_duz_dgamma
  real(kind=CUSTOM_REAL) :: b_dux_dxl,b_duz_dxl,b_dux_dzl,b_duz_dzl
  real(kind=CUSTOM_REAL) :: b_dwx_dxi,b_dwx_dgamma,b_dwz_dxi,b_dwz_dgamma
  real(kind=CUSTOM_REAL) :: b_dwx_dxl,b_dwz_dxl,b_dwx_dzl,b_dwz_dzl
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xz,sigma_zz
  real(kind=CUSTOM_REAL) :: b_sigma_xx,b_sigma_xz,b_sigma_zz
  real(kind=CUSTOM_REAL) :: xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl



  ! loop on all the coupling edges
  do inum = 1,num_solid_poro_edges

    ! get the edge of the elastic element
    ispec_elastic = solid_poro_elastic_ispec(inum)
    iedge_elastic = solid_poro_elastic_iedge(inum)

    ! get the corresponding edge of the poroelastic element
    ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
    iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

    ! implement 1D coupling along the edge
    do ipoin1D = 1,NGLLX

      ! get point values for the poroelastic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_poroelastic)
      j = jvalue_inverse(ipoin1D,iedge_poroelastic)
      iglob = ibool(i,j,ispec_poroelastic)

      ! get poroelastic domain paramters
      phil = porosity(kmato(ispec_poroelastic))
      tortl = tortuosity(kmato(ispec_poroelastic))
      !solid properties
      mul_s = poroelastcoef(2,1,kmato(ispec_poroelastic))
      kappal_s = poroelastcoef(3,1,kmato(ispec_poroelastic)) - 4._CUSTOM_REAL*mul_s/3._CUSTOM_REAL
      rhol_s = density(1,kmato(ispec_poroelastic))
      !fluid properties
      kappal_f = poroelastcoef(1,2,kmato(ispec_poroelastic))
      rhol_f = density(2,kmato(ispec_poroelastic))
      !frame properties
      mul_fr = poroelastcoef(2,3,kmato(ispec_poroelastic))
      kappal_fr = poroelastcoef(3,3,kmato(ispec_poroelastic)) - 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
      rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f
      !Biot coefficients for the input phi
      D_biot = kappal_s*(1._CUSTOM_REAL + phil*(kappal_s/kappal_f - 1._CUSTOM_REAL))
      H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) + &
               kappal_fr + 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
      C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
      M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)
      mul_G = mul_fr
      lambdal_G = H_biot - 2._CUSTOM_REAL*mul_fr
      lambdalplus2mul_G = lambdal_G + TWO*mul_G

      ! derivative along x and along z for u_s and w
      dux_dxi = ZERO
      duz_dxi = ZERO

      dux_dgamma = ZERO
      duz_dgamma = ZERO

      dwx_dxi = ZERO
      dwz_dxi = ZERO

      dwx_dgamma = ZERO
      dwz_dgamma = ZERO

      if(SIMULATION_TYPE == 3) then
        b_dux_dxi = ZERO
        b_duz_dxi = ZERO

        b_dux_dgamma = ZERO
        b_duz_dgamma = ZERO

        b_dwx_dxi = ZERO
        b_dwz_dxi = ZERO

        b_dwx_dgamma = ZERO
        b_dwz_dgamma = ZERO
      endif

      ! first double loop over GLL points to compute and store gradients
      ! we can merge the two loops because NGLLX == NGLLZ
      do k = 1,NGLLX
        dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
        duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
        dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
        duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)

        dwx_dxi = dwx_dxi + displw_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
        dwz_dxi = dwz_dxi + displw_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
        dwx_dgamma = dwx_dgamma + displw_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
        dwz_dgamma = dwz_dgamma + displw_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
        if(SIMULATION_TYPE == 3) then
          b_dux_dxi = b_dux_dxi + b_displs_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
          b_duz_dxi = b_duz_dxi + b_displs_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
          b_dux_dgamma = b_dux_dgamma + b_displs_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
          b_duz_dgamma = b_duz_dgamma + b_displs_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)

          b_dwx_dxi = b_dwx_dxi + b_displw_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
          b_dwz_dxi = b_dwz_dxi + b_displw_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
          b_dwx_dgamma = b_dwx_dgamma + b_displw_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
          b_dwz_dgamma = b_dwz_dgamma + b_displw_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
        endif
      enddo

      xixl = xix(i,j,ispec_poroelastic)
      xizl = xiz(i,j,ispec_poroelastic)
      gammaxl = gammax(i,j,ispec_poroelastic)
      gammazl = gammaz(i,j,ispec_poroelastic)

      ! derivatives of displacement
      dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
      dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

      duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
      duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

      dwx_dxl = dwx_dxi*xixl + dwx_dgamma*gammaxl
      dwx_dzl = dwx_dxi*xizl + dwx_dgamma*gammazl

      dwz_dxl = dwz_dxi*xixl + dwz_dgamma*gammaxl
      dwz_dzl = dwz_dxi*xizl + dwz_dgamma*gammazl

      if(SIMULATION_TYPE == 3) then
        b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
        b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

        b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
        b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl

        b_dwx_dxl = b_dwx_dxi*xixl + b_dwx_dgamma*gammaxl
        b_dwx_dzl = b_dwx_dxi*xizl + b_dwx_dgamma*gammazl

        b_dwz_dxl = b_dwz_dxi*xixl + b_dwz_dgamma*gammaxl
        b_dwz_dzl = b_dwz_dxi*xizl + b_dwz_dgamma*gammazl
      endif
      ! compute stress tensor (include attenuation or anisotropy if needed)

      ! no attenuation
      sigma_xx = lambdalplus2mul_G*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
      sigma_xz = mul_G*(duz_dxl + dux_dzl)
      sigma_zz = lambdalplus2mul_G*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

      if(SIMULATION_TYPE == 3) then
        b_sigma_xx = lambdalplus2mul_G*b_dux_dxl + lambdal_G*b_duz_dzl + C_biot*(b_dwx_dxl + b_dwz_dzl)
        b_sigma_xz = mul_G*(b_duz_dxl + b_dux_dzl)
        b_sigma_zz = lambdalplus2mul_G*b_duz_dzl + lambdal_G*b_dux_dxl + C_biot*(b_dwx_dxl + b_dwz_dzl)
      endif
      ! get point values for the elastic domain, which matches our side in the inverse direction
      ii2 = ivalue(ipoin1D,iedge_elastic)
      jj2 = jvalue(ipoin1D,iedge_elastic)
      iglob = ibool(ii2,jj2,ispec_elastic)

      ! get elastic properties
      lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec_elastic))
      mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec_elastic))
      lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec_elastic))

      ! derivative along x and along z for u_s and w
      dux_dxi = ZERO
      duz_dxi = ZERO

      dux_dgamma = ZERO
      duz_dgamma = ZERO

      if(SIMULATION_TYPE == 3) then
        b_dux_dxi = ZERO
        b_duz_dxi = ZERO

        b_dux_dgamma = ZERO
        b_duz_dgamma = ZERO
      endif

      ! first double loop over GLL points to compute and store gradients
      ! we can merge the two loops because NGLLX == NGLLZ
      do k = 1,NGLLX
        dux_dxi = dux_dxi + displ_elastic(1,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
        duz_dxi = duz_dxi + displ_elastic(3,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
        dux_dgamma = dux_dgamma + displ_elastic(1,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)
        duz_dgamma = duz_dgamma + displ_elastic(3,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)

        if(SIMULATION_TYPE == 3) then
          b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
          b_duz_dxi = b_duz_dxi + b_displ_elastic(3,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
          b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)
          b_duz_dgamma = b_duz_dgamma + b_displ_elastic(3,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)
        endif
      enddo

      xixl = xix(ii2,jj2,ispec_elastic)
      xizl = xiz(ii2,jj2,ispec_elastic)
      gammaxl = gammax(ii2,jj2,ispec_elastic)
      gammazl = gammaz(ii2,jj2,ispec_elastic)

      ! derivatives of displacement
      dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
      dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

      duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
      duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

      if(SIMULATION_TYPE == 3) then
        b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
        b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

        b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
        b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl
      endif
      ! compute stress tensor
      ! full anisotropy
      if(kmato(ispec_elastic) == 2) then
        ! implement anisotropy in 2D
        if(assign_external_model) then
          c11 = c11ext(ii2,jj2,ispec_elastic)
          c13 = c13ext(ii2,jj2,ispec_elastic)
          c15 = c15ext(ii2,jj2,ispec_elastic)
          c33 = c33ext(ii2,jj2,ispec_elastic)
          c35 = c35ext(ii2,jj2,ispec_elastic)
          c55 = c55ext(ii2,jj2,ispec_elastic)
          c12 = c12ext(ii2,jj2,ispec_elastic)
          c23 = c23ext(ii2,jj2,ispec_elastic)
          c25 = c25ext(ii2,jj2,ispec_elastic)
        else
          c11 = anisotropy(1,kmato(ispec_elastic))
          c13 = anisotropy(2,kmato(ispec_elastic))
          c15 = anisotropy(3,kmato(ispec_elastic))
          c33 = anisotropy(4,kmato(ispec_elastic))
          c35 = anisotropy(5,kmato(ispec_elastic))
          c55 = anisotropy(6,kmato(ispec_elastic))
          c12 = anisotropy(7,kmato(ispec_elastic))
          c23 = anisotropy(8,kmato(ispec_elastic))
          c25 = anisotropy(9,kmato(ispec_elastic))
        endif

        sigma_xx = sigma_xx + c11*dux_dxl + c15*(duz_dxl + dux_dzl) + c13*duz_dzl
        sigma_zz = sigma_zz + c13*dux_dxl + c35*(duz_dxl + dux_dzl) + c33*duz_dzl
        sigma_xz = sigma_xz + c15*dux_dxl + c55*(duz_dxl + dux_dzl) + c35*duz_dzl
      else
        ! no attenuation
        sigma_xx = sigma_xx + lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
        sigma_xz = sigma_xz + mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
        sigma_zz = sigma_zz + lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
      endif

      if(SIMULATION_TYPE == 3) then
        b_sigma_xx = b_sigma_xx + lambdaplus2mu_unrelaxed_elastic*b_dux_dxl + lambdal_unrelaxed_elastic*b_duz_dzl
        b_sigma_xz = b_sigma_xz + mul_unrelaxed_elastic*(b_duz_dxl + b_dux_dzl)
        b_sigma_zz = b_sigma_zz + lambdaplus2mu_unrelaxed_elastic*b_duz_dzl + lambdal_unrelaxed_elastic*b_dux_dxl
      endif

      ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
      ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
      ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
      ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
      ! Blackwell Science, page 110, equation (4.60).
      if(iedge_poroelastic == ITOP) then
        xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        weight = jacobian1D * wxgll(i)
      else if(iedge_poroelastic == IBOTTOM) then
        xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        weight = jacobian1D * wxgll(i)
      else if(iedge_poroelastic ==ILEFT) then
        xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)
      else if(iedge_poroelastic ==IRIGHT) then
        xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)
      endif

      accel_elastic(1,iglob) = accel_elastic(1,iglob) - weight * (sigma_xx*nx + sigma_xz*nz)/2.d0
      accel_elastic(3,iglob) = accel_elastic(3,iglob) - weight * (sigma_xz*nx + sigma_zz*nz)/2.d0

      if(SIMULATION_TYPE == 3) then
        b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - weight * (b_sigma_xx*nx + b_sigma_xz*nz)/2.d0
        b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) - weight * (b_sigma_xz*nx + b_sigma_zz*nz)/2.d0
      endif !if(SIMULATION_TYPE == 3) then

    enddo
  enddo
 end subroutine compute_coupling_viscoelastic_po

!========================================================================

! for viscoelastic solver

 subroutine compute_coupling_viscoelastic_po_backward()

  use specfem_par, only: SIMULATION_TYPE,num_solid_poro_edges,&
                         ibool,wxgll,wzgll,xix,xiz,gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse, &
                         hprime_xx,hprime_zz, &
                         solid_poro_elastic_ispec,solid_poro_elastic_iedge, &
                         solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge,&
                         kmato,porosity,tortuosity,poroelastcoef,density, &
                         assign_external_model,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,anisotropy, &
                         displ_elastic,b_displ_elastic,displs_poroelastic,displw_poroelastic, &
                         b_displs_poroelastic,b_displw_poroelastic,&
                         accel_elastic,b_accel_elastic

  implicit none
  include 'constants.h'

  !local variables
  integer :: inum,ispec_elastic,iedge_elastic,ispec_poroelastic,iedge_poroelastic, &
             i,j,k,ipoin1D,iglob,ii2,jj2
  double precision :: phil,tortl,mul_s,kappal_s,rhol_s,kappal_f,rhol_f, &
                      mul_fr,kappal_fr,rhol_bar,D_biot,H_biot,C_biot,M_biot, &
                      mul_G,lambdal_G,lambdalplus2mul_G,lambdal_unrelaxed_elastic, &
                      mul_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic, &
                      c11,c13,c15,c33,c35,c55,c12,c23,c25

  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: dwx_dxl,dwz_dxl,dwx_dzl,dwz_dzl
  real(kind=CUSTOM_REAL) :: b_dux_dxi,b_dux_dgamma,b_duz_dxi,b_duz_dgamma
  real(kind=CUSTOM_REAL) :: b_dux_dxl,b_duz_dxl,b_dux_dzl,b_duz_dzl
  real(kind=CUSTOM_REAL) :: b_dwx_dxi,b_dwx_dgamma,b_dwz_dxi,b_dwz_dgamma
  real(kind=CUSTOM_REAL) :: b_dwx_dxl,b_dwz_dxl,b_dwx_dzl,b_dwz_dzl
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xz,sigma_zz
  real(kind=CUSTOM_REAL) :: b_sigma_xx,b_sigma_xz,b_sigma_zz
  real(kind=CUSTOM_REAL) :: xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl



  ! loop on all the coupling edges
  do inum = 1,num_solid_poro_edges

    ! get the edge of the elastic element
    ispec_elastic = solid_poro_elastic_ispec(inum)
    iedge_elastic = solid_poro_elastic_iedge(inum)

    ! get the corresponding edge of the poroelastic element
    ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
    iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

    ! implement 1D coupling along the edge
    do ipoin1D = 1,NGLLX

      ! get point values for the poroelastic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_poroelastic)
      j = jvalue_inverse(ipoin1D,iedge_poroelastic)
      iglob = ibool(i,j,ispec_poroelastic)

      ! get poroelastic domain paramters
      phil = porosity(kmato(ispec_poroelastic))
      tortl = tortuosity(kmato(ispec_poroelastic))
      !solid properties
      mul_s = poroelastcoef(2,1,kmato(ispec_poroelastic))
      kappal_s = poroelastcoef(3,1,kmato(ispec_poroelastic)) - 4._CUSTOM_REAL*mul_s/3._CUSTOM_REAL
      rhol_s = density(1,kmato(ispec_poroelastic))
      !fluid properties
      kappal_f = poroelastcoef(1,2,kmato(ispec_poroelastic))
      rhol_f = density(2,kmato(ispec_poroelastic))
      !frame properties
      mul_fr = poroelastcoef(2,3,kmato(ispec_poroelastic))
      kappal_fr = poroelastcoef(3,3,kmato(ispec_poroelastic)) - 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
      rhol_bar =  (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f
      !Biot coefficients for the input phi
      D_biot = kappal_s*(1._CUSTOM_REAL + phil*(kappal_s/kappal_f - 1._CUSTOM_REAL))
      H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) + &
               kappal_fr + 4._CUSTOM_REAL*mul_fr/3._CUSTOM_REAL
      C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
      M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)
      mul_G = mul_fr
      lambdal_G = H_biot - 2._CUSTOM_REAL*mul_fr
      lambdalplus2mul_G = lambdal_G + TWO*mul_G

      ! derivative along x and along z for u_s and w
      b_dux_dxi = ZERO
      b_duz_dxi = ZERO

      b_dux_dgamma = ZERO
      b_duz_dgamma = ZERO

      b_dwx_dxi = ZERO
      b_dwz_dxi = ZERO

      b_dwx_dgamma = ZERO
      b_dwz_dgamma = ZERO

      ! first double loop over GLL points to compute and store gradients
      ! we can merge the two loops because NGLLX == NGLLZ
      do k = 1,NGLLX
        dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
        duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
        dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
        duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)

        dwx_dxi = dwx_dxi + displw_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
        dwz_dxi = dwz_dxi + displw_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
        dwx_dgamma = dwx_dgamma + displw_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
        dwz_dgamma = dwz_dgamma + displw_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
        if(SIMULATION_TYPE == 3) then
          b_dux_dxi = b_dux_dxi + b_displs_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
          b_duz_dxi = b_duz_dxi + b_displs_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
          b_dux_dgamma = b_dux_dgamma + b_displs_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
          b_duz_dgamma = b_duz_dgamma + b_displs_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)

          b_dwx_dxi = b_dwx_dxi + b_displw_poroelastic(1,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
          b_dwz_dxi = b_dwz_dxi + b_displw_poroelastic(2,ibool(k,j,ispec_poroelastic))*hprime_xx(i,k)
          b_dwx_dgamma = b_dwx_dgamma + b_displw_poroelastic(1,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
          b_dwz_dgamma = b_dwz_dgamma + b_displw_poroelastic(2,ibool(i,k,ispec_poroelastic))*hprime_zz(j,k)
        endif
      enddo

      xixl = xix(i,j,ispec_poroelastic)
      xizl = xiz(i,j,ispec_poroelastic)
      gammaxl = gammax(i,j,ispec_poroelastic)
      gammazl = gammaz(i,j,ispec_poroelastic)

      ! derivatives of displacement
      dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
      dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

      duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
      duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

      dwx_dxl = dwx_dxi*xixl + dwx_dgamma*gammaxl
      dwx_dzl = dwx_dxi*xizl + dwx_dgamma*gammazl

      dwz_dxl = dwz_dxi*xixl + dwz_dgamma*gammaxl
      dwz_dzl = dwz_dxi*xizl + dwz_dgamma*gammazl

      if(SIMULATION_TYPE == 3) then
        b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
        b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

        b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
        b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl

        b_dwx_dxl = b_dwx_dxi*xixl + b_dwx_dgamma*gammaxl
        b_dwx_dzl = b_dwx_dxi*xizl + b_dwx_dgamma*gammazl

        b_dwz_dxl = b_dwz_dxi*xixl + b_dwz_dgamma*gammaxl
        b_dwz_dzl = b_dwz_dxi*xizl + b_dwz_dgamma*gammazl
      endif
      ! compute stress tensor (include attenuation or anisotropy if needed)

      ! no attenuation
      sigma_xx = lambdalplus2mul_G*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
      sigma_xz = mul_G*(duz_dxl + dux_dzl)
      sigma_zz = lambdalplus2mul_G*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

      if(SIMULATION_TYPE == 3) then
        b_sigma_xx = lambdalplus2mul_G*b_dux_dxl + lambdal_G*b_duz_dzl + C_biot*(b_dwx_dxl + b_dwz_dzl)
        b_sigma_xz = mul_G*(b_duz_dxl + b_dux_dzl)
        b_sigma_zz = lambdalplus2mul_G*b_duz_dzl + lambdal_G*b_dux_dxl + C_biot*(b_dwx_dxl + b_dwz_dzl)
      endif
      ! get point values for the elastic domain, which matches our side in the inverse direction
      ii2 = ivalue(ipoin1D,iedge_elastic)
      jj2 = jvalue(ipoin1D,iedge_elastic)
      iglob = ibool(ii2,jj2,ispec_elastic)

      ! get elastic properties
      lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec_elastic))
      mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec_elastic))
      lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec_elastic))

      ! derivative along x and along z for u_s and w
      dux_dxi = ZERO
      duz_dxi = ZERO

      dux_dgamma = ZERO
      duz_dgamma = ZERO

      if(SIMULATION_TYPE == 3) then
        b_dux_dxi = ZERO
        b_duz_dxi = ZERO

        b_dux_dgamma = ZERO
        b_duz_dgamma = ZERO
      endif

      ! first double loop over GLL points to compute and store gradients
      ! we can merge the two loops because NGLLX == NGLLZ
      do k = 1,NGLLX
        dux_dxi = dux_dxi + displ_elastic(1,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
        duz_dxi = duz_dxi + displ_elastic(3,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
        dux_dgamma = dux_dgamma + displ_elastic(1,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)
        duz_dgamma = duz_dgamma + displ_elastic(3,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)

        if(SIMULATION_TYPE == 3) then
          b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
          b_duz_dxi = b_duz_dxi + b_displ_elastic(3,ibool(k,jj2,ispec_elastic))*hprime_xx(ii2,k)
          b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)
          b_duz_dgamma = b_duz_dgamma + b_displ_elastic(3,ibool(ii2,k,ispec_elastic))*hprime_zz(jj2,k)
        endif
      enddo

      xixl = xix(ii2,jj2,ispec_elastic)
      xizl = xiz(ii2,jj2,ispec_elastic)
      gammaxl = gammax(ii2,jj2,ispec_elastic)
      gammazl = gammaz(ii2,jj2,ispec_elastic)

      ! derivatives of displacement
      dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
      dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

      duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
      duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

      if(SIMULATION_TYPE == 3) then
        b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
        b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

        b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
        b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl
      endif
      ! compute stress tensor
      ! full anisotropy
      if(kmato(ispec_elastic) == 2) then
        ! implement anisotropy in 2D
        if(assign_external_model) then
          c11 = c11ext(ii2,jj2,ispec_elastic)
          c13 = c13ext(ii2,jj2,ispec_elastic)
          c15 = c15ext(ii2,jj2,ispec_elastic)
          c33 = c33ext(ii2,jj2,ispec_elastic)
          c35 = c35ext(ii2,jj2,ispec_elastic)
          c55 = c55ext(ii2,jj2,ispec_elastic)
          c12 = c12ext(ii2,jj2,ispec_elastic)
          c23 = c23ext(ii2,jj2,ispec_elastic)
          c25 = c25ext(ii2,jj2,ispec_elastic)
        else
          c11 = anisotropy(1,kmato(ispec_elastic))
          c13 = anisotropy(2,kmato(ispec_elastic))
          c15 = anisotropy(3,kmato(ispec_elastic))
          c33 = anisotropy(4,kmato(ispec_elastic))
          c35 = anisotropy(5,kmato(ispec_elastic))
          c55 = anisotropy(6,kmato(ispec_elastic))
          c12 = anisotropy(7,kmato(ispec_elastic))
          c23 = anisotropy(8,kmato(ispec_elastic))
          c25 = anisotropy(9,kmato(ispec_elastic))
        endif

        sigma_xx = sigma_xx + c11*dux_dxl + c15*(duz_dxl + dux_dzl) + c13*duz_dzl
        sigma_zz = sigma_zz + c13*dux_dxl + c35*(duz_dxl + dux_dzl) + c33*duz_dzl
        sigma_xz = sigma_xz + c15*dux_dxl + c55*(duz_dxl + dux_dzl) + c35*duz_dzl
      else
        ! no attenuation
        sigma_xx = sigma_xx + lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
        sigma_xz = sigma_xz + mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
        sigma_zz = sigma_zz + lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
      endif

      if(SIMULATION_TYPE == 3) then
        b_sigma_xx = b_sigma_xx + lambdaplus2mu_unrelaxed_elastic*b_dux_dxl + lambdal_unrelaxed_elastic*b_duz_dzl
        b_sigma_xz = b_sigma_xz + mul_unrelaxed_elastic*(b_duz_dxl + b_dux_dzl)
        b_sigma_zz = b_sigma_zz + lambdaplus2mu_unrelaxed_elastic*b_duz_dzl + lambdal_unrelaxed_elastic*b_dux_dxl
      endif

      ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
      ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
      ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
      ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
      ! Blackwell Science, page 110, equation (4.60).
      if(iedge_poroelastic == ITOP) then
        xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        weight = jacobian1D * wxgll(i)
      else if(iedge_poroelastic == IBOTTOM) then
        xxi = + gammaz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        zxi = - gammax(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        weight = jacobian1D * wxgll(i)
      else if(iedge_poroelastic ==ILEFT) then
        xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)
      else if(iedge_poroelastic ==IRIGHT) then
        xgamma = - xiz(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        zgamma = + xix(i,j,ispec_poroelastic) * jacobian(i,j,ispec_poroelastic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)
      endif

      accel_elastic(1,iglob) = accel_elastic(1,iglob) - weight * (sigma_xx*nx + sigma_xz*nz)/2.d0
      accel_elastic(3,iglob) = accel_elastic(3,iglob) - weight * (sigma_xz*nx + sigma_zz*nz)/2.d0

      if(SIMULATION_TYPE == 3) then
        b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - weight * (b_sigma_xx*nx + b_sigma_xz*nz)/2.d0
        b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) - weight * (b_sigma_xz*nx + b_sigma_zz*nz)/2.d0
      endif !if(SIMULATION_TYPE == 3) then

    enddo
  enddo
 end subroutine compute_coupling_viscoelastic_po_backward


