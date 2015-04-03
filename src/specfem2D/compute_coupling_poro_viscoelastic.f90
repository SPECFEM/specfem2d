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
! for poro solver

 subroutine compute_coupling_poro_viscoelastic()

  use specfem_par, only: SIMULATION_TYPE,num_solid_poro_edges,&
                         ibool,wxgll,wzgll,xix,xiz,gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse, &
                         hprime_xx,hprime_zz, &
                         solid_poro_elastic_ispec,solid_poro_elastic_iedge, &
                         solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge,&
                         kmato,porosity,tortuosity,poroelastcoef,density, &
                         assign_external_model,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,anisotropy, &
                         displ_elastic,b_displ_elastic,displs_poroelastic,displw_poroelastic, &
                         b_displs_poroelastic,b_displw_poroelastic, &
                         accels_poroelastic,b_accels_poroelastic

  implicit none
  include 'constants.h'

  !local variables
  integer :: inum,ispec_elastic,iedge_elastic,ispec_poroelastic,iedge_poroelastic, &
             i,j,k,ipoin1D,iglob
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
  real(kind=CUSTOM_REAL) :: sigmap,b_sigmap

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

      ! get point values for the elastic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_elastic)
      j = jvalue_inverse(ipoin1D,iedge_elastic)
      iglob = ibool(i,j,ispec_elastic)

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
        dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
        duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
        dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
        duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec_elastic))*hprime_zz(j,k)

        if(SIMULATION_TYPE == 3) then
          b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
          b_duz_dxi = b_duz_dxi + b_displ_elastic(3,ibool(k,j,ispec_elastic))*hprime_xx(i,k)
          b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
          b_duz_dgamma = b_duz_dgamma + b_displ_elastic(3,ibool(i,k,ispec_elastic))*hprime_zz(j,k)
        endif
      enddo

      xixl = xix(i,j,ispec_elastic)
      xizl = xiz(i,j,ispec_elastic)
      gammaxl = gammax(i,j,ispec_elastic)
      gammazl = gammaz(i,j,ispec_elastic)

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
          c11 = c11ext(i,j,ispec_elastic)
          c13 = c13ext(i,j,ispec_elastic)
          c15 = c15ext(i,j,ispec_elastic)
          c33 = c33ext(i,j,ispec_elastic)
          c35 = c35ext(i,j,ispec_elastic)
          c55 = c55ext(i,j,ispec_elastic)
          c12 = c12ext(i,j,ispec_elastic)
          c23 = c23ext(i,j,ispec_elastic)
          c25 = c25ext(i,j,ispec_elastic)
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
        sigma_xx = c11*dux_dxl + c15*(duz_dxl + dux_dzl) + c13*duz_dzl
        sigma_zz = c13*dux_dxl + c35*(duz_dxl + dux_dzl) + c33*duz_dzl
        sigma_xz = c15*dux_dxl + c55*(duz_dxl + dux_dzl) + c35*duz_dzl
      else
        ! no attenuation
        sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
        sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
        sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
      endif

      if(SIMULATION_TYPE == 3) then
        b_sigma_xx = lambdaplus2mu_unrelaxed_elastic*b_dux_dxl + lambdal_unrelaxed_elastic*b_duz_dzl
        b_sigma_xz = mul_unrelaxed_elastic*(b_duz_dxl + b_dux_dzl)
        b_sigma_zz = lambdaplus2mu_unrelaxed_elastic*b_duz_dzl + lambdal_unrelaxed_elastic*b_dux_dxl
      endif ! if(SIMULATION_TYPE == 3)

      ! get point values for the poroelastic side
      i = ivalue(ipoin1D,iedge_poroelastic)
      j = jvalue(ipoin1D,iedge_poroelastic)
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
      ! compute stress tensor

      ! no attenuation
      sigma_xx = sigma_xx + lambdalplus2mul_G*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
      sigma_xz = sigma_xz + mul_G*(duz_dxl + dux_dzl)
      sigma_zz = sigma_zz + lambdalplus2mul_G*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

      sigmap = C_biot*(dux_dxl + duz_dzl) + M_biot*(dwx_dxl + dwz_dzl)

      if(SIMULATION_TYPE == 3) then
        b_sigma_xx = b_sigma_xx + lambdalplus2mul_G*b_dux_dxl + lambdal_G*b_duz_dzl + C_biot*(b_dwx_dxl + b_dwz_dzl)
        b_sigma_xz = b_sigma_xz + mul_G*(b_duz_dxl + b_dux_dzl)
        b_sigma_zz = b_sigma_zz + lambdalplus2mul_G*b_duz_dzl + lambdal_G*b_dux_dxl + C_biot*(b_dwx_dxl + b_dwz_dzl)
        b_sigmap = C_biot*(b_dux_dxl + b_duz_dzl) + M_biot*(b_dwx_dxl + b_dwz_dzl)
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

      ! contribution to the solid phase
      accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) + &
        weight*((sigma_xx*nx + sigma_xz*nz)/2.d0 -phil/tortl*sigmap*nx)

      accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + &
        weight*((sigma_xz*nx + sigma_zz*nz)/2.d0 -phil/tortl*sigmap*nz)

      ! contribution to the fluid phase
      ! w = 0

      if(SIMULATION_TYPE == 3) then
        ! contribution to the solid phase
        b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) + &
        weight*((b_sigma_xx*nx + b_sigma_xz*nz)/2.d0 -phil/tortl*b_sigmap*nx)

        b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + &
        weight*((b_sigma_xz*nx + b_sigma_zz*nz)/2.d0 -phil/tortl*b_sigmap*nz)

        ! contribution to the fluid phase
        ! w = 0
      endif !if(SIMULATION_TYPE == 3) then

    enddo

  enddo

 end subroutine compute_coupling_poro_viscoelastic

!========================================================================

 subroutine compute_coupling_poro_viscoelastic_for_stabilization()

  use specfem_par, only: SIMULATION_TYPE,num_solid_poro_edges,ibool,ivalue,jvalue, &
                         solid_poro_elastic_ispec,solid_poro_elastic_iedge, &
                         solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge,&
                         veloc_elastic,b_veloc_elastic,accel_elastic,b_accel_elastic, &
                         accels_poroelastic,b_accels_poroelastic, &
                         velocs_poroelastic,b_velocs_poroelastic, &
                         accelw_poroelastic,b_accelw_poroelastic, &
                         velocw_poroelastic,b_velocw_poroelastic, &
                         icount,rmass_inverse_elastic_one,rmass_inverse_elastic_three, &
                         rmass_s_inverse_poroelastic,&
                         time_stepping_scheme,deltatover2,b_deltatover2

  implicit none
  include 'constants.h'

  !local variables
  integer :: inum,ispec_elastic,iedge_elastic,ispec_poroelastic,iedge_poroelastic, &
             i,j,ipoin1D,iglob

  icount(:)=ZERO

  ! loop on all the coupling edges
  do inum = 1,num_solid_poro_edges
     ! get the edge of the elastic element
     ispec_elastic = solid_poro_elastic_ispec(inum)
     iedge_elastic = solid_poro_elastic_iedge(inum)
     ! get the corresponding edge of the poroelastic element
     ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
     iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

     do ipoin1D = 1,NGLLX
       ! recovering original velocities and accelerations on boundaries (elastic side)
       i = ivalue(ipoin1D,iedge_poroelastic)
       j = jvalue(ipoin1D,iedge_poroelastic)
       iglob = ibool(i,j,ispec_poroelastic)
       icount(iglob) = icount(iglob) + 1

       if( icount(iglob) ==1 ) then
         if( time_stepping_scheme == 1 ) then
           veloc_elastic(1,iglob) = veloc_elastic(1,iglob) - deltatover2*accel_elastic(1,iglob)
           veloc_elastic(3,iglob) = veloc_elastic(3,iglob) - deltatover2*accel_elastic(3,iglob)
           accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic_one(iglob)
           accel_elastic(3,iglob) = accel_elastic(3,iglob) / rmass_inverse_elastic_three(iglob)
           ! recovering original velocities and accelerations on boundaries (poro side)
           velocs_poroelastic(1,iglob) = velocs_poroelastic(1,iglob) - deltatover2*accels_poroelastic(1,iglob)
           velocs_poroelastic(2,iglob) = velocs_poroelastic(2,iglob) - deltatover2*accels_poroelastic(2,iglob)
           accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
           accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)
           ! assembling accelerations
           accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
                                    ( 1.0/rmass_inverse_elastic_one(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
           accel_elastic(3,iglob) = ( accel_elastic(3,iglob) + accels_poroelastic(2,iglob) ) / &
                                    ( 1.0/rmass_inverse_elastic_three(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
           accels_poroelastic(1,iglob) = accel_elastic(1,iglob)
           accels_poroelastic(2,iglob) = accel_elastic(3,iglob)
           ! updating velocities
           velocs_poroelastic(1,iglob) = velocs_poroelastic(1,iglob) + deltatover2*accels_poroelastic(1,iglob)
           velocs_poroelastic(2,iglob) = velocs_poroelastic(2,iglob) + deltatover2*accels_poroelastic(2,iglob)
           veloc_elastic(1,iglob) = veloc_elastic(1,iglob) + deltatover2*accel_elastic(1,iglob)
           veloc_elastic(3,iglob) = veloc_elastic(3,iglob) + deltatover2*accel_elastic(3,iglob)
           ! zeros w
           accelw_poroelastic(1,iglob) = ZERO
           accelw_poroelastic(2,iglob) = ZERO
           velocw_poroelastic(1,iglob) = ZERO
           velocw_poroelastic(2,iglob) = ZERO
         endif

!         if(time_stepping_scheme == 2) then
!        recovering original velocities and accelerations on boundaries (elastic side)
!      veloc_elastic = veloc_elastic - beta_LDDRK(i_stage) * veloc_elastic_LDDRK
!      displ_elastic = displ_elastic - beta_LDDRK(i_stage) * displ_elastic_LDDRK
!      veloc_elastic_LDDRK = (veloc_elastic_LDDRK - deltat * accel_elastic) / alpha_LDDRK(i_stage)
!      displ_elastic_LDDRK = (displ_elastic_LDDRK - deltat * veloc_elastic) / alpha_LDDRK(i_stage)
!            accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic(iglob)
!            accel_elastic(3,iglob) = accel_elastic(3,iglob) / rmass_inverse_elastic(iglob)

            ! recovering original velocities and accelerations on boundaries (poro side)
!      velocs_poroelastic = velocs_poroelastic - beta_LDDRK(i_stage) * velocs_poroelastic_LDDRK
!      displs_poroelastic = displs_poroelastic - beta_LDDRK(i_stage) * displs_poroelastic_LDDRK
!      velocs_poroelastic_LDDRK = (velocs_poroelastic_LDDRK - deltat * accels_poroelastic) / alpha_LDDRK(i_stage)
!      displs_poroelastic_LDDRK = (velocs_poroelastic_LDDRK - deltat * velocs_poroelastic) / alpha_LDDRK(i_stage)
!            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
!            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)

            ! assembling accelerations
!            accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accel_elastic(3,iglob) = ( accel_elastic(3,iglob) + accels_poroelastic(2,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accels_poroelastic(1,iglob) = accel_elastic(1,iglob)
!            accels_poroelastic(2,iglob) = accel_elastic(3,iglob)

      ! updating velocities
            ! updating velocities(elastic side)
!      veloc_elastic_LDDRK = alpha_LDDRK(i_stage) * veloc_elastic_LDDRK + deltat * accel_elastic
!      displ_elastic_LDDRK = alpha_LDDRK(i_stage) * displ_elastic_LDDRK + deltat * veloc_elastic
!      veloc_elastic = veloc_elastic + beta_LDDRK(i_stage) * veloc_elastic_LDDRK
!      displ_elastic = displ_elastic + beta_LDDRK(i_stage) * displ_elastic_LDDRK
            ! updating velocities(poro side)
!      velocs_poroelastic_LDDRK = alpha_LDDRK(i_stage) * velocs_poroelastic_LDDRK + deltat * accels_poroelastic
!      displs_poroelastic_LDDRK = alpha_LDDRK(i_stage) * displs_poroelastic_LDDRK + deltat * velocs_poroelastic
!      velocs_poroelastic = velocs_poroelastic + beta_LDDRK(i_stage) * velocs_poroelastic_LDDRK
!      displs_poroelastic = displs_poroelastic + beta_LDDRK(i_stage) * displs_poroelastic_LDDRK

            ! zeros w
!            accelw_poroelastic(1,iglob) = ZERO
!            accelw_poroelastic(2,iglob) = ZERO
!            velocw_poroelastic(1,iglob) = ZERO
!            velocw_poroelastic(2,iglob) = ZERO
!            endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if(time_stepping_scheme == 3) then

        ! recovering original velocities and accelerations on boundaries (elastic side)
!        if(i_stage==1 .or. i_stage==2 .or. i_stage==3) then

!        if(i_stage == 1)weight_rk = 0.5d0
!        if(i_stage == 2)weight_rk = 0.5d0
!        if(i_stage == 3)weight_rk = 1.0d0

!        veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) - weight_rk * accel_elastic_rk(1,iglob,i_stage)
!  veloc_elastic(3,iglob) = veloc_elastic_initial_rk(3,iglob) - weight_rk * accel_elastic_rk(3,iglob,i_stage)
!        displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) - weight_rk * veloc_elastic_rk(1,iglob,i_stage)
!  displ_elastic(3,iglob) = displ_elastic_initial_rk(3,iglob) - weight_rk * veloc_elastic_rk(3,iglob,i_stage)


!        else if(i_stage==4) then

!        veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (accel_elastic_rk(1,iglob,1) + 2.0d0 * accel_elastic_rk(1,iglob,2) + &
!         2.0d0 * accel_elastic_rk(1,iglob,3) + accel_elastic_rk(1,iglob,4))

!        veloc_elastic(3,iglob) = veloc_elastic_initial_rk(3,iglob) - 1.0d0 / 6.0d0 * &
!        (accel_elastic_rk(3,iglob,1) + 2.0d0 * accel_elastic_rk(3,iglob,2) + &
!         2.0d0 * accel_elastic_rk(3,iglob,3) + accel_elastic_rk(3,iglob,4))

!        displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (veloc_elastic_rk(1,iglob,1) + 2.0d0 * veloc_elastic_rk(1,iglob,2) + &
!         2.0d0 * veloc_elastic_rk(1,iglob,3) + veloc_elastic_rk(1,iglob,4))

!        displ_elastic(3,iglob) = displ_elastic_initial_rk(3,iglob) - 1.0d0 / 6.0d0 * &
!        (veloc_elastic_rk(3,iglob,1) + 2.0d0 * veloc_elastic_rk(3,iglob,2) + &
!         2.0d0 * veloc_elastic_rk(3,iglob,3) + veloc_elastic_rk(3,iglob,4))

!        endif

!        accel_elastic(1,iglob) = accel_elastic(1,iglob) / rmass_inverse_elastic(iglob)
!        accel_elastic(3,iglob) = accel_elastic(3,iglob) / rmass_inverse_elastic(iglob)

!        accel_elastic_rk(1,iglob,i_stage) = accel_elastic(1,iglob) / deltat
!        accel_elastic_rk(3,iglob,i_stage) = accel_elastic(3,iglob) / deltat
!        veloc_elastic_rk(1,iglob,i_stage) = veloc_elastic(1,iglob) / deltat
!        veloc_elastic_rk(3,iglob,i_stage) = veloc_elastic(3,iglob) / deltat


        ! recovering original velocities and accelerations on boundaries (poro side)
!        if(i_stage==1 .or. i_stage==2 .or. i_stage==3) then

!        if(i_stage == 1)weight_rk = 0.5d0
!        if(i_stage == 2)weight_rk = 0.5d0
!        if(i_stage == 3)weight_rk = 1.0d0

!        velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) - weight_rk * accels_poroelastic_rk(1,iglob,i_stage)
!  velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) - weight_rk * accels_poroelastic_rk(2,iglob,i_stage)
!        displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) - weight_rk * velocs_poroelastic_rk(1,iglob,i_stage)
!  displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) - weight_rk * velocs_poroelastic_rk(2,iglob,i_stage)


!        else if(i_stage==4) then

!        velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (accels_poroelastic_rk(1,iglob,1) + 2.0d0 * accels_poroelastic_rk(1,iglob,2) + &
!         2.0d0 * accels_poroelastic_rk(1,iglob,3) + accels_poroelastic_rk(1,iglob,4))

!        velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) - 1.0d0 / 6.0d0 * &
!        (accels_poroelastic_rk(2,iglob,1) + 2.0d0 * accels_poroelastic_rk(2,iglob,2) + &
!         2.0d0 * accels_poroelastic_rk(2,iglob,3) + accels_poroelastic_rk(2,iglob,4))

!        displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) - 1.0d0 / 6.0d0 * &
!        (velocs_poroelastic_rk(1,iglob,1) + 2.0d0 * velocs_poroelastic_rk(1,iglob,2) + &
!         2.0d0 * velocs_poroelastic_rk(1,iglob,3) + velocs_poroelastic_rk(1,iglob,4))

!        displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) - 1.0d0 / 6.0d0 * &
!        (velocs_poroelastic_rk(2,iglob,1) + 2.0d0 * velocs_poroelastic_rk(2,iglob,2) + &
!         2.0d0 * velocs_poroelastic_rk(2,iglob,3) + velocs_poroelastic_rk(2,iglob,4))

!        endif

!        accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
!        accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)

!        accels_poroelastic_rk(1,iglob,i_stage) = accels_poroelastic(1,iglob) / deltat
!        accels_poroelastic_rk(2,iglob,i_stage) = accels_poroelastic(2,iglob) / deltat
!        velocs_poroelastic_rk(1,iglob,i_stage) = velocs_poroelastic(1,iglob) / deltat
!        velocs_poroelastic_rk(2,iglob,i_stage) = velocs_poroelastic(2,iglob) / deltat


        ! assembling accelerations
!            accel_elastic(1,iglob) = ( accel_elastic(1,iglob) + accels_poroelastic(1,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accel_elastic(3,iglob) = ( accel_elastic(3,iglob) + accels_poroelastic(2,iglob) ) / &
!                                   ( 1.0/rmass_inverse_elastic(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
!            accels_poroelastic(1,iglob) = accel_elastic(1,iglob)
!            accels_poroelastic(2,iglob) = accel_elastic(3,iglob)

   ! updating velocities
        ! updating velocities(elastic side)

 !       accel_elastic_rk(1,iglob,i_stage) = accel_elastic(1,iglob) * deltat
 !       accel_elastic_rk(3,iglob,i_stage) = accel_elastic(3,iglob) * deltat

 !       if(i_stage==1 .or. i_stage==2 .or. i_stage==3) then

 !       if(i_stage == 1)weight_rk = 0.5d0
 !       if(i_stage == 2)weight_rk = 0.5d0
 !       if(i_stage == 3)weight_rk = 1.0d0

 !       veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) + weight_rk * accel_elastic_rk(1,iglob,i_stage)
 ! veloc_elastic(3,iglob) = veloc_elastic_initial_rk(3,iglob) + weight_rk * accel_elastic_rk(3,iglob,i_stage)
 !       displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) + weight_rk * veloc_elastic_rk(1,iglob,i_stage)
 ! displ_elastic(3,iglob) = displ_elastic_initial_rk(3,iglob) + weight_rk * veloc_elastic_rk(3,iglob,i_stage)


 !       else if(i_stage==4) then

 !       veloc_elastic(1,iglob) = veloc_elastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (accel_elastic_rk(1,iglob,1) + 2.0d0 * accel_elastic_rk(1,iglob,2) + &
 !        2.0d0 * accel_elastic_rk(1,iglob,3) + accel_elastic_rk(1,iglob,4))
!
 !       veloc_elastic(3,iglob) = veloc_elastic_initial_rk(3,iglob) + 1.0d0 / 6.0d0 * &
 !       (accel_elastic_rk(3,iglob,1) + 2.0d0 * accel_elastic_rk(3,iglob,2) + &
 !        2.0d0 * accel_elastic_rk(3,iglob,3) + accel_elastic_rk(3,iglob,4))

 !       displ_elastic(1,iglob) = displ_elastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (veloc_elastic_rk(1,iglob,1) + 2.0d0 * veloc_elastic_rk(1,iglob,2) + &
 !        2.0d0 * veloc_elastic_rk(1,iglob,3) + veloc_elastic_rk(1,iglob,4))

 !       displ_elastic(3,iglob) = displ_elastic_initial_rk(3,iglob) + 1.0d0 / 6.0d0 * &
 !       (veloc_elastic_rk(3,iglob,1) + 2.0d0 * veloc_elastic_rk(3,iglob,2) + &
 !        2.0d0 * veloc_elastic_rk(3,iglob,3) + veloc_elastic_rk(3,iglob,4))

 !       endif
        ! updating velocities(poro side)

 !       accels_poroelastic_rk(1,iglob,i_stage) = deltat * accels_poroelastic(1,iglob)
 !       accels_poroelastic_rk(2,iglob,i_stage) = deltat * accels_poroelastic(2,iglob)
 !       velocs_poroelastic_rk(1,iglob,i_stage) = deltat * velocs_poroelastic(1,iglob)
 !       velocs_poroelastic_rk(2,iglob,i_stage) = deltat * velocs_poroelastic(2,iglob)


 !       if(i_stage==1 .or. i_stage==2 .or. i_stage==3) then

 !       if(i_stage == 1)weight_rk = 0.5d0
 !       if(i_stage == 2)weight_rk = 0.5d0
 !       if(i_stage == 3)weight_rk = 1.0d0

 !       velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) + weight_rk * accels_poroelastic_rk(1,iglob,i_stage)
 ! velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) + weight_rk * accels_poroelastic_rk(2,iglob,i_stage)
 !       displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) + weight_rk * velocs_poroelastic_rk(1,iglob,i_stage)
 ! displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) + weight_rk * velocs_poroelastic_rk(2,iglob,i_stage)


 !       else if(i_stage==4) then

 !       velocs_poroelastic(1,iglob) = velocs_poroelastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (accels_poroelastic_rk(1,iglob,1) + 2.0d0 * accels_poroelastic_rk(1,iglob,2) + &
 !        2.0d0 * accels_poroelastic_rk(1,iglob,3) + accels_poroelastic_rk(1,iglob,4))

 !       velocs_poroelastic(2,iglob) = velocs_poroelastic_initial_rk(2,iglob) + 1.0d0 / 6.0d0 * &
 !       (accels_poroelastic_rk(2,iglob,1) + 2.0d0 * accels_poroelastic_rk(2,iglob,2) + &
 !        2.0d0 * accels_poroelastic_rk(2,iglob,3) + accels_poroelastic_rk(2,iglob,4))
!
 !       displs_poroelastic(1,iglob) = displs_poroelastic_initial_rk(1,iglob) + 1.0d0 / 6.0d0 * &
 !       (velocs_poroelastic_rk(1,iglob,1) + 2.0d0 * velocs_poroelastic_rk(1,iglob,2) + &
 !        2.0d0 * velocs_poroelastic_rk(1,iglob,3) + velocs_poroelastic_rk(1,iglob,4))
!
 !       displs_poroelastic(2,iglob) = displs_poroelastic_initial_rk(2,iglob) + 1.0d0 / 6.0d0 * &
 !       (velocs_poroelastic_rk(2,iglob,1) + 2.0d0 * velocs_poroelastic_rk(2,iglob,2) + &
 !        2.0d0 * velocs_poroelastic_rk(2,iglob,3) + velocs_poroelastic_rk(2,iglob,4))

 !       endif

 !     endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if( SIMULATION_TYPE == 3 ) then
         b_veloc_elastic(1,iglob) = b_veloc_elastic(1,iglob) - b_deltatover2*b_accel_elastic(1,iglob)
         b_veloc_elastic(3,iglob) = b_veloc_elastic(3,iglob) - b_deltatover2*b_accel_elastic(3,iglob)
         b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) / rmass_inverse_elastic_one(iglob)
         b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) / rmass_inverse_elastic_three(iglob)
         ! recovering original velocities and accelerations on boundaries (poro side)
         b_velocs_poroelastic(1,iglob) = b_velocs_poroelastic(1,iglob) - b_deltatover2*b_accels_poroelastic(1,iglob)
         b_velocs_poroelastic(2,iglob) = b_velocs_poroelastic(2,iglob) - b_deltatover2*b_accels_poroelastic(2,iglob)
         b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) / rmass_s_inverse_poroelastic(iglob)
         b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) / rmass_s_inverse_poroelastic(iglob)
         ! assembling accelerations
         b_accel_elastic(1,iglob) = ( b_accel_elastic(1,iglob) + b_accels_poroelastic(1,iglob) ) / &
                        ( 1.0/rmass_inverse_elastic_one(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
         b_accel_elastic(3,iglob) = ( b_accel_elastic(3,iglob) + b_accels_poroelastic(2,iglob) ) / &
                        ( 1.0/rmass_inverse_elastic_three(iglob) +1.0/rmass_s_inverse_poroelastic(iglob) )
         b_accels_poroelastic(1,iglob) = b_accel_elastic(1,iglob)
         b_accels_poroelastic(2,iglob) = b_accel_elastic(3,iglob)
         ! updating velocities
         b_velocs_poroelastic(1,iglob) = b_velocs_poroelastic(1,iglob) + b_deltatover2*b_accels_poroelastic(1,iglob)
         b_velocs_poroelastic(2,iglob) = b_velocs_poroelastic(2,iglob) + b_deltatover2*b_accels_poroelastic(2,iglob)
         b_veloc_elastic(1,iglob) = b_veloc_elastic(1,iglob) + b_deltatover2*b_accel_elastic(1,iglob)
         b_veloc_elastic(3,iglob) = b_veloc_elastic(3,iglob) + b_deltatover2*b_accel_elastic(3,iglob)
         ! zeros w
         b_accelw_poroelastic(1,iglob) = ZERO
         b_accelw_poroelastic(2,iglob) = ZERO
         b_velocw_poroelastic(1,iglob) = ZERO
         b_velocw_poroelastic(2,iglob) = ZERO
       endif !if(SIMULATION_TYPE == 3)
     endif !if(icount(iglob) ==1)
   enddo
  enddo

 end subroutine compute_coupling_poro_viscoelastic_for_stabilization


