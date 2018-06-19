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

  subroutine compute_pressure_whole_medium(i_field)

! compute pressure in acoustic elements and in elastic elements

  use specfem_par, only: CUSTOM_REAL,NGLLX,NGLLZ,nspec,ibool,displ_elastic,displs_poroelastic,displw_poroelastic, &
                         potential_dot_dot_acoustic,potential_acoustic,b_displ_elastic,b_displs_poroelastic, &
                         b_displw_poroelastic,b_potential_dot_dot_acoustic,b_potential_acoustic
  use specfem_par_movie, only: vector_field_display

  implicit none

  ! forwrad or adjoint wavefield)
  integer :: i_field

  ! local parameters
  integer :: i,j,ispec,iglob
  ! pressure in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: pressure_element

  ! loop over spectral elements
  do ispec = 1,nspec

    ! compute pressure in this element
    if (i_field == 1) then
      call compute_pressure_one_element(ispec,pressure_element,displ_elastic,displs_poroelastic,displw_poroelastic, &
                                        potential_dot_dot_acoustic,potential_acoustic)
    else
      call compute_pressure_one_element(ispec,pressure_element,b_displ_elastic,b_displs_poroelastic,b_displw_poroelastic, &
                                        b_potential_dot_dot_acoustic,b_potential_acoustic)
    endif
    ! use vector_field_display as temporary storage, store pressure in its second component
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        vector_field_display(2,iglob) = pressure_element(i,j)
      enddo
    enddo

  enddo

  end subroutine compute_pressure_whole_medium

!
!=====================================================================
!

  subroutine compute_pressure_one_element(ispec,pressure_element,displ_elastic,displs_poroelastic,displw_poroelastic, &
                                          potential_dot_dot_acoustic,potential_acoustic)

! compute pressure in acoustic elements and in elastic elements

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,ZERO,TWO,NDIM

  use specfem_par, only: N_SLS,nglob,ispec_is_elastic,ispec_is_acoustic,ispec_is_poroelastic, &
    ispec_is_anisotropic,kmato,poroelastcoef,assign_external_model,vpext,vsext,rhoext, &
    ATTENUATION_VISCOELASTIC,AXISYM,is_on_the_axis, &
    anisotropy,c11ext,c12ext,c13ext,c15ext,c23ext,c25ext,c33ext,c35ext,c55ext, &
    hprimebar_xx,hprime_xx,hprime_zz,xix,xiz,gammax,gammaz,jacobian,ibool,coord,e1,e11,USE_TRICK_FOR_BETTER_PRESSURE

  implicit none

  integer,intent(in) :: ispec
  ! pressure in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: pressure_element
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: potential_dot_dot_acoustic,potential_acoustic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ_elastic,displs_poroelastic,displw_poroelastic

  ! local variables
  real(kind=CUSTOM_REAL) :: e1_sum,e11_sum
  integer :: i_sls
  integer :: i,j,k,iglob
  real(kind=CUSTOM_REAL) :: sigma_yy !! ,sigmap
  real(kind=CUSTOM_REAL) :: sigma_thetatheta
  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: denst
  real(kind=CUSTOM_REAL) :: cpl,csl
  real(kind=CUSTOM_REAL) :: mu_G,lambdal_G,lambdalplus2mul_G

  ! for anisotropy
  double precision ::  c11,c15,c13,c33,c35,c55,c12,c23,c25
  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl
  ! to evaluate cpI, cpII, and cs, and rI (poroelastic medium)
  double precision :: phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot
  double precision :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic
  double precision :: sigma_xx,sigma_zz

  double precision :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  double precision :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  double precision :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  double precision :: dwx_dxl,dwz_dzl

  double precision :: xxi

  ! initializes
  pressure_element(:,:) = 0._CUSTOM_REAL

! if elastic element
!
! from L. S. Bennethum, Compressibility Moduli for Porous Materials Incorporating Volume Fraction,
! J. Engrg. Mech., vol. 132(11), p. 1205-1214 (2006), below equation (5):
! for a 3D isotropic solid, pressure is defined in terms of the trace of the stress tensor as
! p = -1/3 (t11 + t22 + t33) where t is the Cauchy stress tensor.

! to compute pressure in 3D in an elastic solid, one uses pressure = - trace(sigma) / 3
! sigma_ij = lambda delta_ij trace(epsilon) + 2 mu epsilon_ij
!          = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_ij
! sigma_xx = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_xx
! sigma_yy = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_yy
! sigma_zz = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_zz
! pressure = - trace(sigma) / 3 = - (lambda + 2/3 mu) trace(epsilon) = - kappa * trace(epsilon)
!
! to compute pressure in 2D in an elastic solid in the plane strain convention i.e. in the P-SV case,
! one still uses pressure = - trace(sigma) / 3 but taking into account the fact
! that the off-plane strain epsilon_zz is zero by definition of the plane strain convention
! but thus the off-plane stress sigma_zz is not equal to zero,
! one has instead:  sigma_zz = lambda * (epsilon_xx + epsilon_yy), thus
! sigma_ij = lambda delta_ij trace(epsilon) + 2 mu epsilon_ij
!          = lambda (epsilon_xx + epsilon_yy) + 2 mu epsilon_ij
! sigma_xx = lambda (epsilon_xx + epsilon_yy) + 2 mu epsilon_xx
! sigma_yy = lambda (epsilon_xx + epsilon_yy) + 2 mu epsilon_yy
! sigma_zz = lambda * (epsilon_xx + epsilon_yy)
! pressure = - trace(sigma) / 3 = - (lambda + 2*mu/3) (epsilon_xx + epsilon_yy)

  if (ispec_is_elastic(ispec)) then
    ! elastic element

    ! get relaxed elastic parameters of current spectral element
    lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
    mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
    lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))

    do j = 1,NGLLZ
      do i = 1,NGLLX

        !--- if external medium, get elastic parameters of current grid point
        if (assign_external_model) then
          cpl = vpext(i,j,ispec)
          csl = vsext(i,j,ispec)
          denst = rhoext(i,j,ispec)
          mul_unrelaxed_elastic = denst*csl*csl
          lambdal_unrelaxed_elastic = denst*cpl*cpl - TWO*mul_unrelaxed_elastic
        endif

        ! derivative along x and along z
        dux_dxi = ZERO
        duz_dxi = ZERO

        dux_dgamma = ZERO
        duz_dgamma = ZERO

        ! first double loop over GLL points to compute and store gradients
        ! we can merge the two loops because NGLLX == NGLLZ

        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            do k = 1,NGLJ
              dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
              duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprimeBar_xx(i,k)
              dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
              duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
            enddo
          else
            do k = 1,NGLJ
              dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
              duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
              dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
              duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
            enddo
          endif
        else
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
          enddo
        endif

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! derivatives of displacement
        dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
        duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

        if (AXISYM .and. is_on_the_axis(ispec) .and. i == 1) then ! d_uz/dr=0 on the axis
          duz_dxl = 0.d0
        endif

! compute diagonal components of the stress tensor (include attenuation or anisotropy if needed)

        if (ATTENUATION_VISCOELASTIC) then

! attenuation is implemented following the memory variable formulation of
! J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993). More details can be found in
! J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation in a linear
! viscoelastic medium, Geophysical Journal International, vol. 95, p. 597-611 (1988).

! When implementing viscoelasticity according to the Carcione 1993 paper, attenuation is
! non-causal rather than causal i.e. wave speed up instead of slowing down
! when attenuation is turned on. We fixed that issue (which is not incorrect but non traditional)
! by taking the unrelaxed state (infinite frequency) as a reference instead of the relaxed state (zero frequency)
! and also using equations in Carcione's 2007 book.
! See file doc/old_problem_attenuation_reference_Specfem2D_fixed_by_Xie_Zhinan.pdf
! and doc/how_we_modified_Carcione_1993_to_make_it_causal_and_include_the_missing_1_over_L_factor.pdf

! See also J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).

! and J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
! in a linear viscoelastic medium, Geophysical Journal International,
! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).

          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              if (is_on_the_axis(ispec) .and. i == 1) then ! First GLJ point
                sigma_xx = 0._CUSTOM_REAL
                sigma_zz = 0._CUSTOM_REAL
                sigma_thetatheta = 0._CUSTOM_REAL
                xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                do k = 1,NGLJ
                  sigma_xx = sigma_xx + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  sigma_zz = sigma_zz + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  sigma_thetatheta = sigma_thetatheta + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                enddo
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                           + lambdal_unrelaxed_elastic*sigma_xx/xxi
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                           + lambdal_unrelaxed_elastic*sigma_zz/xxi
                sigma_thetatheta = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                   + lambdaplus2mu_unrelaxed_elastic*sigma_thetatheta/xxi
              else ! Not first GLJ point
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                           + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                           + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_thetatheta = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                        + lambdaplus2mu_unrelaxed_elastic &
                                        * displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              endif
            else ! Not on the axis
              sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                         + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                         + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              sigma_thetatheta = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                      + lambdaplus2mu_unrelaxed_elastic &
                                      * displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
            endif
            sigma_yy = sigma_thetatheta
          else ! Not axisym
            ! compute the stress using the unrelaxed Lame parameters (Carcione 2007 page 125)
            sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
            ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
            sigma_yy = lambdal_unrelaxed_elastic*(dux_dxl + duz_dzl)
            sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
          endif

          ! add the memory variables using the relaxed parameters (Carcione 2007 page 125)
          ! beware: there is a bug in Carcione's equation (2c) of his 1993 paper for sigma_zz, we fixed it in the code below
          e1_sum = 0._CUSTOM_REAL
          e11_sum = 0._CUSTOM_REAL

          do i_sls = 1,N_SLS
            e1_sum = e1_sum + e1(i_sls,i,j,ispec)
            e11_sum = e11_sum + e11(i_sls,i,j,ispec)
          enddo

          sigma_xx = sigma_xx + (lambdal_unrelaxed_elastic + mul_unrelaxed_elastic) * e1_sum &
                      + TWO * mul_unrelaxed_elastic * e11_sum
          ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
          sigma_yy = sigma_yy + (lambdal_unrelaxed_elastic + mul_unrelaxed_elastic) * e1_sum
          sigma_zz = sigma_zz + (lambdal_unrelaxed_elastic + mul_unrelaxed_elastic) * e1_sum &
                      - TWO * mul_unrelaxed_elastic * e11_sum

        else

          ! no attenuation

          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              if (is_on_the_axis(ispec) .and. i == 1) then ! First GLJ point
                sigma_xx = 0._CUSTOM_REAL
                sigma_zz = 0._CUSTOM_REAL
                sigma_thetatheta = 0._CUSTOM_REAL
                xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                do k = 1,NGLJ
                  sigma_xx = sigma_xx + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  sigma_zz = sigma_zz + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  sigma_thetatheta = sigma_thetatheta + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                enddo
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                           + lambdal_unrelaxed_elastic*sigma_xx/xxi
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                           + lambdal_unrelaxed_elastic*sigma_zz/xxi
                sigma_thetatheta = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                   + lambdaplus2mu_unrelaxed_elastic*sigma_thetatheta/xxi
              else ! Not first GLJ point
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                           + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                           + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_thetatheta = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                        + lambdaplus2mu_unrelaxed_elastic &
                                        * displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              endif
            else ! Not on the axis
              sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                         + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                         + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              sigma_thetatheta = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                      + lambdaplus2mu_unrelaxed_elastic &
                                      * displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
            endif
            sigma_yy = sigma_thetatheta
          else ! Not axisym
            sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
            ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
            sigma_yy = lambdal_unrelaxed_elastic*(dux_dxl + duz_dzl)
            sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
          endif

        endif

        ! full anisotropy
        if (ispec_is_anisotropic(ispec)) then
          if (assign_external_model) then
            c11 = c11ext(i,j,ispec)
            c15 = c15ext(i,j,ispec)
            c13 = c13ext(i,j,ispec)
            c33 = c33ext(i,j,ispec)
            c35 = c35ext(i,j,ispec)
            c55 = c55ext(i,j,ispec)
            c12 = c12ext(i,j,ispec)
            c23 = c23ext(i,j,ispec)
            c25 = c25ext(i,j,ispec)
          else
            c11 = anisotropy(1,kmato(ispec))
            c13 = anisotropy(2,kmato(ispec))
            c15 = anisotropy(3,kmato(ispec))
            c33 = anisotropy(4,kmato(ispec))
            c35 = anisotropy(5,kmato(ispec))
            c55 = anisotropy(6,kmato(ispec))
            c12 = anisotropy(7,kmato(ispec))
            c23 = anisotropy(8,kmato(ispec))
            c25 = anisotropy(9,kmato(ispec))
          endif

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          if (AXISYM .and. is_on_the_axis(ispec) .and. i == 1) then ! d_uz/dr=0 on the axis
            duz_dxl = 0.d0
          endif

          ! implement anisotropy in 2D
          sigma_xx = c11*dux_dxl + c13*duz_dzl + c15*(duz_dxl + dux_dzl)
          ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
          if (c12 < 1.e-7 .or. c23 < 1.e-7) call stop_the_code( &
              'cannot compute pressure for an anisotropic material if c12 or c23 are zero')
          sigma_yy = c12*dux_dxl + c23*duz_dzl + c25*(duz_dxl + dux_dzl)
          sigma_zz = c13*dux_dxl + c33*duz_dzl + c35*(duz_dxl + dux_dzl)

        endif

        ! store pressure
        ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
        pressure_element(i,j) = - (sigma_xx + sigma_yy + sigma_zz) / 3.d0

      enddo
    enddo

  else if (ispec_is_poroelastic(ispec)) then
    ! poro-elastic element

    lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
    mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))

    ! get poroelastic parameters of current spectral element
    call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)

    ! Biot coefficients for the input phi
    call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

    ! where T = G:grad u_s + C div w I
    ! and T_f = C div u_s I + M div w I
    ! we are expressing lambdaplus2mu, lambda, and mu for G, C, and M
    mu_G = mu_fr
    lambdal_G = H_biot - TWO*mu_fr
    lambdalplus2mul_G = lambdal_G + TWO*mu_G

    do j = 1,NGLLZ
      do i = 1,NGLLX

        ! derivative along x and along z
        dux_dxi = ZERO
        duz_dxi = ZERO

        dux_dgamma = ZERO
        duz_dgamma = ZERO

        dwx_dxi = ZERO
        dwz_dxi = ZERO

        dwx_dgamma = ZERO
        dwz_dgamma = ZERO

        ! first double loop over GLL points to compute and store gradients
        ! we can merge the two loops because NGLLX == NGLLZ
        do k = 1,NGLLX
          dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
          duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)

          dwx_dxi = dwx_dxi + displw_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          dwz_dxi = dwz_dxi + displw_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          dwx_dgamma = dwx_dgamma + displw_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
          dwz_dgamma = dwz_dgamma + displw_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! derivatives of displacement
        dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
        duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

        dwx_dxl = dwx_dxi*xixl + dwx_dgamma*gammaxl
        dwz_dzl = dwz_dxi*xizl + dwz_dgamma*gammazl

! compute diagonal components of the stress tensor (include attenuation if needed)

        if (ATTENUATION_VISCOELASTIC) then

! attenuation is implemented following the memory variable formulation of
! J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993). More details can be found in
! J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation in a linear
! viscoelastic medium, Geophysical Journal International, vol. 95, p. 597-611 (1988).

! When implementing viscoelasticity according to the Carcione 1993 paper, attenuation is
! non-causal rather than causal i.e. wave speed up instead of slowing down
! when attenuation is turned on. We fixed that issue (which is not incorrect but non traditional)
! by taking the unrelaxed state (infinite frequency) as a reference instead of the relaxed state (zero frequency)
! and also using equations in Carcione's 2007 book.
! See file doc/old_problem_attenuation_reference_Specfem2D_fixed_by_Xie_Zhinan.pdf
! and doc/how_we_modified_Carcione_1993_to_make_it_causal_and_include_the_missing_1_over_L_factor.pdf

! See also J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).

! and J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
! in a linear viscoelastic medium, Geophysical Journal International,
! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).

          ! compute the stress using the unrelaxed Lame parameters (Carcione 2007 page 125)
          sigma_xx = (lambdal_unrelaxed_elastic + 2.0 * mul_unrelaxed_elastic)*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
          ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
!         sigma_yy = ...  ! it is not zero because of the plane strain formulation, thus it should be computed here
          call stop_the_code('pressure calculation not implemented for poroelastic media yet, you should compute sigma_yy here')
          sigma_zz = (lambdal_unrelaxed_elastic + 2.0 * mul_unrelaxed_elastic)*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl

          ! add the memory variables using the relaxed parameters (Carcione 2007 page 125)
          ! beware: there is a bug in  Carcione's equation (2c) of his 1993 paper for sigma_zz, we fixed it in the code below
          e1_sum = 0._CUSTOM_REAL
          e11_sum = 0._CUSTOM_REAL

          do i_sls = 1,N_SLS
            e1_sum = e1_sum + e1(i_sls,i,j,ispec)
            e11_sum = e11_sum + e11(i_sls,i,j,ispec)
          enddo

          sigma_xx = sigma_xx + (lambdal_unrelaxed_elastic + mul_unrelaxed_elastic) * e1_sum &
                    + TWO * mul_unrelaxed_elastic * e11_sum
          ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
!         sigma_yy = ...  ! it is not zero because of the plane strain formulation, thus it should be computed here
          call stop_the_code('pressure calculation not implemented for poroelastic media yet, you should compute sigma_yy here')
          sigma_zz = sigma_zz + (lambdal_unrelaxed_elastic + mul_unrelaxed_elastic) * e1_sum &
                    - TWO * mul_unrelaxed_elastic * e11_sum

        else

          ! no attenuation
          sigma_xx = lambdalplus2mul_G*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
          ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
!         sigma_yy = ...  ! it is not zero because of the plane strain formulation, thus it should be computed here
!! DK DK we temporarily set it to zero here to avoid a warning by certain compilers about an unassigned variable
!! DK DK when computing pressure a few lines below; this should be removed and replaced with a true calculation
          sigma_yy = 0.d0
          call stop_the_code('pressure calculation not implemented for poroelastic media yet, you should compute sigma_yy here')
          sigma_zz = lambdalplus2mul_G*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

!         sigmap = C_biot*(dux_dxl + duz_dzl) + M_biot*(dwx_dxl + dwz_dzl)

        endif

        ! store pressure
        ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
        pressure_element(i,j) = - (sigma_xx + sigma_yy + sigma_zz) / 3.d0
!       pressure_element2(i,j) = - sigmap
      enddo
    enddo

! pressure = - Chi_dot_dot if acoustic element
  else if (ispec_is_acoustic(ispec)) then
    ! acoustic element

    if (USE_TRICK_FOR_BETTER_PRESSURE) then
      ! use a trick to increase accuracy of pressure seismograms in fluid (acoustic) elements:
      ! use the second derivative of the source for the source time function instead of the source itself,
      ! and then record -potential_acoustic() as pressure seismograms instead of -potential_dot_dot_acoustic();
      ! this is mathematically equivalent, but numerically significantly more accurate because in the explicit
      ! Newmark time scheme acceleration is accurate at zeroth order while displacement is accurate at second order,
      ! thus in fluid elements potential_dot_dot_acoustic() is accurate at zeroth order while potential_acoustic()
      ! is accurate at second order and thus contains significantly less numerical noise.
      ! compute pressure on the fluid/porous medium edge
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          ! store pressure
          pressure_element(i,j) = - potential_acoustic(iglob)
        enddo
      enddo
    else
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          ! store pressure
          pressure_element(i,j) = - potential_dot_dot_acoustic(iglob)
        enddo
      enddo
    endif
  endif ! end of test if acoustic or elastic element

  end subroutine compute_pressure_one_element

