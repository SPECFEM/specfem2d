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

  subroutine compute_forces_viscoelastic(accel_elastic,veloc_elastic,displ_elastic, &
                                         displ_elastic_old,dux_dxl_old,duz_dzl_old, &
                                         dux_dzl_plus_duz_dxl_old, &
                                         PML_BOUNDARY_CONDITIONS,e1,e11,e13,iphase)

  ! compute forces for the elastic elements
  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM, &
    ONE,TWO,PI,TINYVAL,FOUR_THIRDS

  use specfem_par, only: nglob,P_SV, &
                         ATTENUATION_VISCOELASTIC,nspec_ATT_el,N_SLS, &
                         ibool,ispec_is_elastic, &
                         xix,xiz,gammax,gammaz, &
                         jacobian,rho_vpstore,mustore,rhostore, &
                         qkappa_attenuation_store,qmu_attenuation_store, &
                         c11store,c12store,c13store,c15store,c22store,c23store,c25store,c33store,c35store,c55store, &
                         ispec_is_anisotropic, &
                         hprime_xx,hprimewgll_xx,hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         coord,iglob_is_forced

  ! overlapping communication
  use specfem_par, only: nspec_inner_elastic,nspec_outer_elastic,phase_ispec_inner_elastic

  ! AXISYM
  use specfem_par, only: AXISYM,is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj

  ! PML arrays
  use specfem_par, only: nspec_PML,ispec_is_PML,ROTATE_PML_ACTIVATE,ROTATE_PML_ANGLE

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ_elastic,veloc_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: accel_elastic

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ_elastic_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT_el),intent(inout) :: dux_dxl_old,duz_dzl_old,dux_dzl_plus_duz_dxl_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_ATT_el,N_SLS),intent(inout) :: e1,e11,e13

  ! CPML coefficients and memory variables
  logical,intent(in) :: PML_BOUNDARY_CONDITIONS

  integer,intent(in) :: iphase

  !---
  !--- local variables
  !---
  integer :: ispec,i,j,k,iglob

  ! spatial derivatives
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxl,dux_dzl,duz_dxl,duz_dzl
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ) :: dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: dummy_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempz1,tempz2
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: tempx3
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: sigma_thetatheta
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempz1l,tempz2l
  real(kind=CUSTOM_REAL) :: theta,ct,st
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xy,sigma_xz,sigma_zy,sigma_zz,sigma_zx
  real(kind=CUSTOM_REAL) :: sigma_xx_prime,sigma_xz_prime,sigma_zz_prime,sigma_zx_prime
  real(kind=CUSTOM_REAL) :: xxi


  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLZ) :: deriv
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: e1_sum,e11_sum,e13_sum

  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic, &
    lambdaplus2mu_unrelaxed_elastic,cpl,rhol,lambdalplusmul_unrelaxed_elastic

  ! for anisotropy
  double precision ::  c11,c15,c13,c33,c35,c55,c12,c23,c25,c22

  ! CPML coefficients and memory variables
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: accel_elastic_PML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old

  !additional variables for CPML implementation in viscoelastic simulation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: kappa_dux_dxl,kappa_duz_dzl,mu_dux_dxl,mu_duz_dzl, &
                                                    mu_dux_dzl,mu_duz_dxl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: kappa_pml_dux_dxl,kappa_pml_duz_dzl, &
                                                    mu_pml_dux_dxl,mu_pml_duz_dzl,mu_pml_dux_dzl,mu_pml_duz_dxl
  double precision :: qkappal,qmul

  integer :: num_elements,ispec_p

  ! this to avoid a warning at execution time about an undefined variable being used
  ! for the SH component in the case of a P-SV calculation, and vice versa
  ! P_SV-case
  sigma_xx = 0._CUSTOM_REAL
  sigma_xz = 0._CUSTOM_REAL
  sigma_zz = 0._CUSTOM_REAL
  sigma_zx = 0._CUSTOM_REAL
  ! SH-case
  sigma_xy = 0._CUSTOM_REAL
  sigma_zy = 0._CUSTOM_REAL

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_elastic
  else
    num_elements = nspec_inner_elastic
  endif

  ! loop over spectral elements
  do ispec_p = 1,num_elements

    ! returns element id from stored element list
    ispec = phase_ispec_inner_elastic(ispec_p,iphase)

    ! only for elastic spectral elements
    if (.not. ispec_is_elastic(ispec)) cycle

    if (ATTENUATION_VISCOELASTIC) then
      qkappal = maxval(qkappa_attenuation_store(:,:,ispec))
      qmul =  maxval(qmu_attenuation_store(:,:,ispec))
    endif

    ! gets local displacement for element
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        dummy_loc(1,i,j) = displ_elastic(1,iglob)
        dummy_loc(2,i,j) = displ_elastic(2,iglob)

        ! stores local array for element xi/gamma/jacobian (for better performance)
        deriv(1,i,j) = xix(i,j,ispec)
        deriv(2,i,j) = xiz(i,j,ispec)
        deriv(3,i,j) = gammax(i,j,ispec)
        deriv(4,i,j) = gammaz(i,j,ispec)
        deriv(5,i,j) = jacobian(i,j,ispec)
      enddo
    enddo

    ! first double loop over GLL points to compute and store gradients
    call mxm_4comp_singleA(dux_dxi,duz_dxi,dux_dgamma,duz_dgamma,dummy_loc,hprime_xx,hprime_zz)

    ! AXISYM case overwrites dux_dxi and duz_dxi
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            ! derivative along x and along z
            dux_dxi(i,j) = 0._CUSTOM_REAL
            duz_dxi(i,j) = 0._CUSTOM_REAL
            do k = 1,NGLJ
              dux_dxi(i,j) = dux_dxi(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
              duz_dxi(i,j) = duz_dxi(i,j) + dummy_loc(2,k,j)*hprimeBar_xx(i,k)
            enddo
          enddo
        enddo
      endif
    endif

    ! gets derivatives of ux and uz with respect to x and z
    do j = 1,NGLLZ
      do i = 1,NGLLX
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)

        ! derivatives of displacement
        dux_dxl(i,j) = dux_dxi(i,j)*xixl + dux_dgamma(i,j)*gammaxl
        dux_dzl(i,j) = dux_dxi(i,j)*xizl + dux_dgamma(i,j)*gammazl

        duz_dxl(i,j) = duz_dxi(i,j)*xixl + duz_dgamma(i,j)*gammaxl
        duz_dzl(i,j) = duz_dxi(i,j)*xizl + duz_dgamma(i,j)*gammazl
      enddo
    enddo

    ! AXISYM case overwrite duz_dxl
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        ! d_uz/dr=0 on the axis
        ! i == 1
        do j = 1,NGLLZ
          duz_dxl(1,j) = 0._CUSTOM_REAL
        enddo
      endif
    endif

    if (PML_BOUNDARY_CONDITIONS) then !viscoelastic PML
      if (ispec_is_PML(ispec)) then
        if (ATTENUATION_VISCOELASTIC) then
          if (qkappal < 9998.999d0 .and. qmul < 9998.999d0) then

            call pml_compute_memory_variables_viscoelastic(ispec,nglob,displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
                                                           kappa_pml_dux_dxl,kappa_pml_duz_dzl, &
                                                           mu_pml_dux_dxl,mu_pml_duz_dzl,mu_pml_dux_dzl, &
                                                           mu_pml_duz_dxl,kappa_dux_dxl,kappa_duz_dzl,mu_dux_dxl, &
                                                           mu_duz_dzl,mu_dux_dzl,mu_duz_dxl)
          else
            call pml_compute_memory_variables_elastic(ispec,nglob,displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
                                                     dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime, &
                                                     PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl, &
                                                     PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old)
          endif
        else
          call pml_compute_memory_variables_elastic(ispec,nglob,displ_elastic_old,dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
                                                   dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime, &
                                                   PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl, &
                                                   PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old)
        endif
      endif ! ispec_is_PML
    endif ! PML_BOUNDARY_CONDITIONS

    ! AXISYM case overwrite duz_dxl and duz_dxl_prime
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        ! d_uz/dr=0 on the axis
        !i == 1
        do j = 1,NGLLZ
          duz_dxl(1,j) = 0._CUSTOM_REAL
          duz_dxl_prime(1,j) = 0._CUSTOM_REAL
        enddo
      endif
    endif

    ! first double loop to compute gradient
    do j = 1,NGLLZ
      do i = 1,NGLLX
        ! get elastic parameters of current grid point
        mul_unrelaxed_elastic = mustore(i,j,ispec)
        rhol = rhostore(i,j,ispec)
        cpl = rho_vpstore(i,j,ispec)/rhol

        lambdal_unrelaxed_elastic = rhol*cpl*cpl - 2._CUSTOM_REAL * mul_unrelaxed_elastic
        lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + 2._CUSTOM_REAL * mul_unrelaxed_elastic
        lambdalplusmul_unrelaxed_elastic = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic

        ! compute stress tensor (include attenuation or anisotropy if needed)
        if (.not. ispec_is_anisotropic(ispec)) then
          ! isotropic elastic
          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              if (i == 1) then
                ! First GLJ point
                sigma_xx = 0._CUSTOM_REAL
                sigma_zz = 0._CUSTOM_REAL
                sigma_thetatheta(i,j) = 0._CUSTOM_REAL
                xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                r_xiplus1(i,j) = xxi
                do k = 1,NGLJ
                  sigma_xx = sigma_xx + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                  sigma_zz = sigma_zz + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                  sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                enddo
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
                          + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + sigma_xx/xxi)
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
                         + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + sigma_zz/xxi)
                sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
                sigma_zx = sigma_xz
                sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
                                   + lambdaplus2mu_unrelaxed_elastic*sigma_thetatheta(i,j)/xxi
              else
                ! Not first GLJ point
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
                           + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
                           + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
                sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
                sigma_zx = sigma_xz
                sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
                                        + lambdaplus2mu_unrelaxed_elastic * dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
                r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
              endif
            else
              ! Not on the axis
              sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) &
                         + lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
              sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) &
                         + lambdal_unrelaxed_elastic * (dux_dxl(i,j) + dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec)))
              sigma_xz = mul_unrelaxed_elastic*(duz_dxl(i,j) + dux_dzl(i,j))
              sigma_zx = sigma_xz
              sigma_thetatheta(i,j) = lambdal_unrelaxed_elastic * (duz_dzl(i,j) + dux_dxl(i,j)) &
                                      + lambdaplus2mu_unrelaxed_elastic * dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
            endif
          else
            ! Not AXISYM
            if (P_SV) then
              ! P_SV case
              ! 2D plane strain assumption:
              !   strain eps_yy == eps_xy == eps_zy == 0 (infinite medium in y-direction
              ! leads to:
              !   T_xx = (lambda + 2 mu) eps_xx + lambda eps_zz
              !   T_zz = (lambda + 2 mu) eps_zz + lambda eps_xx
              !   T_xz = T_zx = mu 2 eps_zx = mu (duz_dx + dux_dz)
              sigma_xx = lambdaplus2mu_unrelaxed_elastic * dux_dxl(i,j) + lambdal_unrelaxed_elastic * duz_dzl(i,j)
              sigma_zz = lambdaplus2mu_unrelaxed_elastic * duz_dzl(i,j) + lambdal_unrelaxed_elastic * dux_dxl(i,j)
              sigma_xz = mul_unrelaxed_elastic * (duz_dxl(i,j) + dux_dzl(i,j))
              sigma_zx = sigma_xz
            else
              ! SH-case
              ! displacement only in y-direction:
              !   ignores u_x and u_z, \partial_y == 0 (no derivatives in y-direction)
              ! leads to:
              !   T_xy = T_yx = mu 2 eps_yx = mu (duy_dx + dux_dy) = mu duy_dx     (sets dux_dy == 0)
              !   T_zy = T_yz = mu 2 eps_yz = mu (duy_dz + duz_dy) = mu duy_dz     (sets duz_dy == 0)
              ! instead of new variable names, we use the same variables from the P-SV case,
              ! but here dux_dxl == duy_dxl and dux_dzl == duy_dzl
              sigma_xy = mul_unrelaxed_elastic * dux_dxl(i,j)
              sigma_zy = mul_unrelaxed_elastic * dux_dzl(i,j)
            endif
          endif ! AXISYM

        else
          ! full anisotropy
          c11 = c11store(i,j,ispec)
          c13 = c13store(i,j,ispec)
          c15 = c15store(i,j,ispec)
          c33 = c33store(i,j,ispec)
          c35 = c35store(i,j,ispec)
          c55 = c55store(i,j,ispec)
          c12 = c12store(i,j,ispec)
          c23 = c23store(i,j,ispec)
          c25 = c25store(i,j,ispec)
          if (AXISYM) then
            c22 = c22store(i,j,ispec) ! This variable is used for axisym simulations only
          endif

          ! implement anisotropy in 2D
          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              if (i == 1) then
                ! first GLJ point, on the axis
                sigma_xx = 0._CUSTOM_REAL
                sigma_zz = 0._CUSTOM_REAL
                sigma_thetatheta(i,j) = 0._CUSTOM_REAL
                xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                r_xiplus1(i,j) = xxi
                do k = 1,NGLJ ! Compute the sum
                  sigma_xx = sigma_xx + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                  sigma_zz = sigma_zz + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                  sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + dummy_loc(1,k,j)*hprimeBar_xx(i,k)
                enddo
                sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*sigma_xx/xxi
                sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*sigma_zz/xxi
                sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
                sigma_zx = sigma_xz
                sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + c22*sigma_thetatheta(i,j)/xxi
              else
                ! not first GLJ point but not on the axis
                sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
                sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
                sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
                sigma_zx = sigma_xz
                sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + &
                                        c22*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
              endif
            else
              ! axisym but not on the axis
              sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c12*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
              sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c23*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
              sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
              sigma_zx = sigma_xz
              sigma_thetatheta(i,j) = c12*dux_dxl(i,j) + c23*duz_dzl(i,j) + c22*dummy_loc(1,i,j)/coord(1,ibool(i,j,ispec))
            endif
          else
            ! not AXISYM
            if (P_SV) then
              ! P_SV case
              ! P_SV case
              ! 2D plane strain assumption:
              !   strain eps_yy == eps_xy == eps_zy == 0 (infinite medium in y-direction
              ! leads to:
              !   T_xx = c11 eps_xx + c13 eps_zz + c15 2 eps_zx
              !   T_zz = c13 eps_xx + c33 eps_zz + c35 2 eps_zx
              !   T_xz = c15 eps_xx + c35 eps_zz + c55 2 eps_zx
              ! using the Voigt notation: T_i = sum_j=1^6 C_ij eps_j
              !                           T_1 = T_xx,  T_2 = T_yy,  T_3 = T_zz
              !                           T_3 = T_yz,  T_4 = T_xz,  T_6 = T_xy
              !                           and same for strain eps_1 = eps_xx, ..
              !
              ! the out-of-plane stress T_yy = c12 eps_xx + c23 eps_zz + c25 2 eps_zx is not computed
              ! as it is not used any further
              sigma_xx = c11*dux_dxl(i,j) + c13*duz_dzl(i,j) + c15*(duz_dxl(i,j) + dux_dzl(i,j))
              sigma_zz = c13*dux_dxl(i,j) + c33*duz_dzl(i,j) + c35*(duz_dxl(i,j) + dux_dzl(i,j))
              sigma_xz = c15*dux_dxl(i,j) + c35*duz_dzl(i,j) + c55*(duz_dxl(i,j) + dux_dzl(i,j))
              sigma_zx = sigma_xz
            else
              ! SH-case
              ! displacement only in y-direction:
              !   ignores u_x and u_z, \partial_y == 0 (no derivatives in y-direction)
              ! leads to:
              !   T_xy = T_yx = c66 2 eps_yx = c66 (duy_dx + dux_dy) = c66 duy_dx     (sets dux_dy == 0)
              !   T_zy = T_yz = c44 2 eps_yz = c44 (duy_dz + duz_dy) = c44 duy_dz     (sets duz_dy == 0)
              ! instead of new variable names, we use the same variables from the P-SV case,
              ! but here dux_dxl == duy_dxl and dux_dzl == duy_dzl
              !
              ! note: since c44 and c66 are not given in the Par_file input,
              !       uses c55 == mu in both directions, thus still isotropic for SH case
              sigma_xy = c55 * dux_dxl(i,j)
              sigma_zy = c55 * dux_dzl(i,j)
            endif
          endif ! AXISYM
        endif ! anisotropic

        ! add the memory variables using the relaxed parameters (Carcione 2007 page 125)
        ! beware: there is a bug in Carcione's equation (2c) of his 1993 paper for sigma_zz, we fixed it in the code below.

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
        if (ATTENUATION_VISCOELASTIC) then
          ! This routine updates the memory variables and with the current grad(displ)
          ! and gets the attenuation contribution (in e*_sum variables)
          call compute_attenuation_viscoelastic(e1,e11,e13,dux_dxl(i,j),dux_dzl(i,j),duz_dxl(i,j),duz_dzl(i,j), &
                                                dux_dxl_old,duz_dzl_old, &
                                                dux_dzl_plus_duz_dxl_old,PML_BOUNDARY_CONDITIONS,i,j,ispec, &
                                                e1_sum,e11_sum,e13_sum)

! use the right formula with 1/N included
! i.e. use the unrelaxed moduli here (see Carcione's book, third edition, equation (3.189))
          if (P_SV) then
            ! P_SV case
            sigma_xx = sigma_xx + lambdalplusmul_unrelaxed_elastic * e1_sum + TWO * mul_unrelaxed_elastic * e11_sum
            sigma_xz = sigma_xz + mul_unrelaxed_elastic * e13_sum
            sigma_zz = sigma_zz + lambdalplusmul_unrelaxed_elastic * e1_sum - TWO * mul_unrelaxed_elastic * e11_sum
            sigma_zx = sigma_xz
          else
            ! SH case
            ! attenuation not implemented yet for SH case ...
            sigma_xy = sigma_xy + 0._CUSTOM_REAL
            sigma_zy = sigma_zy + 0._CUSTOM_REAL
          endif
        endif ! ATTENUATION_VISCOELASTIC

        if (PML_BOUNDARY_CONDITIONS) then
          ! note: P-SV waves implement the PML equations:
          !         rho F^{-1}[-omega^2 s_x s_z] * u_x = d/dx{ (lambda+2mu) F^{-1}[s_z/s_x] * dux_dx + lambda duz_dz }
          !                                            + d/dz{ mu duz_dx + mu F^{-1}[s_x/s_z] * dux_dz }
          !         rho F^{-1}[-omega^2 s_x s_z] * u_z = d/dx{ mu F^{-1}[s_z/s_x] * duz_dx + mu dux_dz }
          !                                            + d/dz{ lambda dux_dx + (lambda+2mu) F^{-1}[s_x/s_z] * duz_dz }
          !       with F^{-1} being the inverse Fourier transform and '*' being convolution.
          !
          !       given the above expressions, we can write the stress components as
          !           sigma_xx = (lambda+2mu) F^{-1}[s_z/s_x] * dux_dx + lambda duz_dz
          !           sigma_zx = mu duz_dx + mu F^{-1}[s_x/s_z] * dux_dz
          !
          !           sigma_zz = (lambda+2mu) F^{-1}[s_x/s_z] * duz_dz + lambda dux_dx
          !           sigma_xz = mu F^{-1}[s_z/s_x] * duz_dx + mu dux_dz
          !
          !       note that the stress becomes non-symmetric.
          !       -> see Xie et al. (2014), appendix B, eq. (B6a) and (B6b), where component y == z here
          !
          !       SH waves implement the PML equation:
          !         rho F^{-1}[-omega**2 s_x s_z] * u_y = d/dx{ mu F^{-1}[s_z/s_x] * duy_dx }
          !                                             + d/dz{ mu F^{-1}[s_x/s_z] * duy_dz }
          !            with F^{-1} being the inverse Fourier transform.
          !
          !       for SH-waves, we only have F^{-1}[..] expressions for the stress components.
          !       these get computed by using the modified dux_dxl,dux_dzl,.. terms.

          !       given that PML_dux_dxl,PML_dux_dzl,.. arrays store the original, unmodified dux_dx,dux_dz,.. strain values,
          !       there are no such expressions for the stress components in the SH-wave case.
          !
          !       terms without an F^{-1}[..] expression still appear in the P-SV case, thus need to be added
          !       together with the PML-modified terms dux_dxl,dux_dzl,.. in those cases.
          !
          !       -> see in pml_compute_memory_variables_elastic() how the dux_dxl,dux_dzl,.. terms are modified
          !          using the recursive convolutional scheme with memory variables rmemory_dux_dz,.. arrays
          !
          !       thus, for the SH-case we already computed: sigma_xy = mul_unrelaxed_elastic * dux_dxl(i,j)
          !                                                  sigma_zy = mul_unrelaxed_elastic * dux_dzl(i,j)
          !       based on the convolutional dux_dxl,.. expressions, and therefore won't need to add anything here.
          !       plus, SH-waves ignore attenuation and the PML rotation cases.
          !
          if (ATTENUATION_VISCOELASTIC) then
            if (ispec_is_PML(ispec)) then
               if (qkappal < 9998.999d0 .and. qmul < 9998.999d0) then
                 sigma_xx = (lambdalplusmul_unrelaxed_elastic * kappa_pml_dux_dxl(i,j) + mul_unrelaxed_elastic * &
                            mu_pml_dux_dxl(i,j)) + (lambdalplusmul_unrelaxed_elastic * kappa_duz_dzl(i,j) - &
                            mul_unrelaxed_elastic * mu_duz_dzl(i,j))
                 sigma_xz = (mul_unrelaxed_elastic * mu_pml_duz_dxl(i,j) + mul_unrelaxed_elastic * mu_dux_dzl(i,j))
                 sigma_zx = (mul_unrelaxed_elastic * mu_duz_dxl(i,j) + mul_unrelaxed_elastic * mu_pml_dux_dzl(i,j))
                 sigma_zz = (lambdalplusmul_unrelaxed_elastic * kappa_dux_dxl(i,j) - mul_unrelaxed_elastic * &
                            mu_dux_dxl(i,j))+ (lambdalplusmul_unrelaxed_elastic * kappa_pml_duz_dzl(i,j) + &
                            mul_unrelaxed_elastic * mu_pml_duz_dzl(i,j))
               else
                 sigma_xx = lambdaplus2mu_unrelaxed_elastic * dux_dxl(i,j) + lambdal_unrelaxed_elastic * PML_duz_dzl(i,j)
                 sigma_zz = lambdaplus2mu_unrelaxed_elastic * duz_dzl(i,j) + lambdal_unrelaxed_elastic * PML_dux_dxl(i,j)
                 sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl(i,j))
                 sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl(i,j))
               endif
            endif
          else ! No ATTENUATION
            if (ispec_is_PML(ispec) .and. nspec_PML > 0) then
              if (ROTATE_PML_ACTIVATE) then
                theta = - ROTATE_PML_ANGLE/180._CUSTOM_REAL*PI
                !debug
                !if (it == 1) write(*,*) 'debug: rotate PML',theta,ROTATE_PML_ACTIVATE,cos(theta),sin(theta)
                ct = cos(theta)
                st = sin(theta)
                sigma_xx_prime = lambdaplus2mu_unrelaxed_elastic * (ct**2*dux_dxl(i,j) + ct*st*duz_dxl(i,j) + &
                                                                    ct*st*dux_dzl(i,j) + st**2*duz_dzl(i,j)) &
                                 + lambdal_unrelaxed_elastic*(st**2*PML_dux_dxl(i,j) - ct*st*PML_duz_dxl(i,j) - &
                                                              ct*st*PML_dux_dzl(i,j) + ct**2*PML_duz_dzl(i,j))

                sigma_xz_prime = mul_unrelaxed_elastic * (-ct*st*dux_dxl(i,j) + ct**2*duz_dxl(i,j) - &
                                                          st**2*dux_dzl(i,j) + ct*st*duz_dzl(i,j)) &
                                 + mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j) - st**2*PML_duz_dxl(i,j) + &
                                                            ct**2*PML_dux_dzl(i,j) + ct*st*PML_duz_dzl(i,j))

                sigma_zx_prime = mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j) + ct**2*PML_duz_dxl(i,j) - &
                                                          st**2*PML_dux_dzl(i,j) + ct*st*PML_duz_dzl(i,j)) &
                                 + mul_unrelaxed_elastic * (-ct*st*dux_dxl_prime(i,j) - st**2*duz_dxl_prime(i,j) + &
                                                          ct**2*dux_dzl_prime(i,j) + ct*st*duz_dzl_prime(i,j))

                sigma_zz_prime = lambdaplus2mu_unrelaxed_elastic*(st**2*dux_dxl_prime(i,j) - ct*st*duz_dxl_prime(i,j) - &
                                                                  ct*st*dux_dzl_prime(i,j) + ct**2*duz_dzl_prime(i,j)) &
                                 + lambdal_unrelaxed_elastic*(ct**2*PML_dux_dxl(i,j) + ct*st*PML_duz_dxl(i,j) + &
                                                              ct*st*PML_dux_dzl(i,j) + st**2*PML_duz_dzl(i,j))

                sigma_xx = ct**2*sigma_xx_prime - ct*st*sigma_xz_prime - ct*st*sigma_zx_prime + st**2*sigma_zz_prime
                sigma_xz = ct*st*sigma_xx_prime + ct**2*sigma_xz_prime - st**2*sigma_zx_prime - ct*st*sigma_zz_prime
                sigma_zx = ct*st*sigma_xx_prime - st**2*sigma_xz_prime + ct**2*sigma_zx_prime - ct*st*sigma_zz_prime
                sigma_zz = st**2*sigma_xx_prime + ct*st*sigma_xz_prime + ct*st*sigma_zx_prime + ct**2*sigma_zz_prime
              else
                ! stress components:
                !   sigma_xx = (lambda+2mu) F^{-1}[s_z/s_x] * dux_dx + lambda duz_dz
                !   sigma_zz = (lambda+2mu) F^{-1}[s_x/s_z] * duz_dz + lambda dux_dx
                !
                !   sigma_zx = mu duz_dx + mu F^{-1}[s_x/s_z] * dux_dz
                !            = mu ( duz_dx + F^{-1}[s_x/s_z] * dux_dz)
                !   sigma_xz = mu F^{-1}[s_z/s_x] * duz_dx + mu dux_dz
                !            = mu ( dux_dz + F^{-1}[s_z/s_x] * duz_dx )
                !
                ! note that PML_dux_dxl,PML_dux_dzl,.. arrays contain the original, unmodified dux_dx,dux_dz,.. strain values.
                !
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl(i,j) + lambdal_unrelaxed_elastic*PML_duz_dzl(i,j)
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl(i,j) + lambdal_unrelaxed_elastic*PML_dux_dxl(i,j)
                sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j) + dux_dzl(i,j))
                sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j) + duz_dxl(i,j))
              endif
            endif
          endif ! ATTENUATION
        endif ! PML_BOUNDARY_CONDITION

        ! weak formulation term based on stress tensor (non-symmetric form)
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)
        jacobianl = deriv(5,i,j)

        !! ABAB with the notations of Komatitsch & Tromp 1999 (with 3 -> 2) :
        ! tempx1(i,j) = w.J.F_{11}^{ij}
        ! tempz1(i,j) = w.J.F_{21}^{ij}
        ! tempx2(i,j) = w.J.F_{12}^{ij}
        ! tempz2(i,j) = w.J.F_{22}^{ij}

        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            ! This is normal, we always add a contribution depending on the value on the axis
            ! i.e. we purposely sum something at point (i,j) with something at point (1,j)
            tempx3(i,j) = jacobian(1,j,ispec) * sigma_thetatheta(1,j)*hprimeBarwglj_xx(1,i)

            ! not first GLJ point
            if (abs(coord(1,ibool(i,j,ispec))) > TINYVAL) then
              if (i == 1) then
                write(*,*) "Element number:",ispec
                call stop_the_code("Error: an axial element is rotated. The code should have been stopped before. Check that your &
                 &coordinates are greater than TINYVAL. Maybe you should also have a look to &
                 &doc/problematic_case_that_we_exclude_for_axisymmetric.pdf")
              endif
              tempx3(i,j) = tempx3(i,j) + wxglj(i) * jacobianl &
                            * sigma_thetatheta(i,j)/(xiglj(i)+ONE) ! this goes to accel_x
            endif

            tempx2(i,j) = r_xiplus1(i,j) * jacobianl &
                          * (sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
            tempz2(i,j) = r_xiplus1(i,j) * jacobianl &
                          * (sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
            tempx1(i,j) = r_xiplus1(i,j) * jacobianl &
                          * (sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
            tempz1(i,j) = r_xiplus1(i,j) * jacobianl &
                          * (sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z
          else
            ! axisym but not on the axis
            tempx2(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
                          *(sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
            tempz2(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
                          *(sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
            tempx1(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
                          *(sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
            tempz1(i,j) = coord(1,ibool(i,j,ispec)) * jacobianl &
                          *(sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z
            tempx3(i,j) = wxgll(i) * jacobianl * sigma_thetatheta(i,j) ! this goes to accel_x
          endif
        else
          ! default (not axisym case)
          if (P_SV) then
            ! P_SV case
            tempx1(i,j) = jacobianl * (sigma_xx*xixl + sigma_zx*xizl) ! this goes to accel_x
            tempz1(i,j) = jacobianl * (sigma_xz*xixl + sigma_zz*xizl) ! this goes to accel_z

            tempx2(i,j) = jacobianl * (sigma_xx*gammaxl + sigma_zx*gammazl) ! this goes to accel_x
            tempz2(i,j) = jacobianl * (sigma_xz*gammaxl + sigma_zz*gammazl) ! this goes to accel_z
          else
            ! SH-case
            tempx1(i,j) = jacobianl * (sigma_xy*xixl + sigma_zy*xizl) ! this goes to accel_x
            tempx2(i,j) = jacobianl * (sigma_xy*gammaxl + sigma_zy*gammazl) ! this goes to accel_x
            tempz1(i,j) = 0._CUSTOM_REAL
            tempz2(i,j) = 0._CUSTOM_REAL
          endif
        endif ! AXISYM
      enddo
    enddo  ! end of the loops on the collocation points i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update the displacement memory variable
    if (PML_BOUNDARY_CONDITIONS) then
      ! calculates contribution from each C-PML element to update acceleration
      call pml_compute_accel_contribution_elastic(ispec,nglob, &
                                                   dummy_loc,displ_elastic_old,veloc_elastic, &
                                                   accel_elastic_PML,r_xiplus1)
    endif

    !
    ! second double-loop over GLL to compute all the terms
    !
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. iglob_is_forced(iglob)) then
              ! along x direction and z direction
              ! and assemble the contributions
              ! we can merge the two loops because NGLLX == NGLLZ

              ! assembles the contributions
              tempx1l = 0._CUSTOM_REAL
              tempx2l = 0._CUSTOM_REAL
              tempz1l = 0._CUSTOM_REAL
              tempz2l = 0._CUSTOM_REAL
              do k = 1,NGLJ
                tempx1l = tempx1l + tempx1(k,j) * hprimeBarwglj_xx(k,i)
                tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
                tempz1l = tempz1l + tempz1(k,j) * hprimeBarwglj_xx(k,i)
                tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
              enddo
              ! sums contributions from each element to the global values
              accel_elastic(1,iglob) = accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxglj(i) * tempx2l)
              accel_elastic(2,iglob) = accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxglj(i) * tempz2l)

              accel_elastic(1,iglob) = accel_elastic(1,iglob) - wzgll(j) * tempx3(i,j)
            endif
          enddo
        enddo
      else
        ! Axisym but not on the axis
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. iglob_is_forced(iglob)) then
              ! assembles the contributions
              tempx1l = 0._CUSTOM_REAL
              tempx2l = 0._CUSTOM_REAL
              tempz1l = 0._CUSTOM_REAL
              tempz2l = 0._CUSTOM_REAL
              do k = 1,NGLLX
                tempx1l = tempx1l + tempx1(k,j) * hprimewgll_xx(k,i)
                tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
                tempz1l = tempz1l + tempz1(k,j) * hprimewgll_xx(k,i)
                tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
              enddo
              ! sums contributions from each element to the global values
              accel_elastic(1,iglob) = accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxgll(i) * tempx2l)
              accel_elastic(2,iglob) = accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxgll(i) * tempz2l)

              accel_elastic(1,iglob) = accel_elastic(1,iglob) - wzgll(j) * tempx3(i,j)
            endif
          enddo
        enddo
      endif
    else
      ! default (not AXISYM case)
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          if (.not. iglob_is_forced(iglob)) then
            ! assembles the contributions
            tempx1l = 0._CUSTOM_REAL
            tempx2l = 0._CUSTOM_REAL
            tempz1l = 0._CUSTOM_REAL
            tempz2l = 0._CUSTOM_REAL
            do k = 1,NGLLX
              tempx1l = tempx1l + tempx1(k,j) * hprimewgll_xx(k,i)
              tempx2l = tempx2l + tempx2(i,k) * hprimewgll_zz(k,j)
              tempz1l = tempz1l + tempz1(k,j) * hprimewgll_xx(k,i)
              tempz2l = tempz2l + tempz2(i,k) * hprimewgll_zz(k,j)
            enddo
            ! sums contributions from each element to the global values
            accel_elastic(1,iglob) = accel_elastic(1,iglob) - (wzgll(j) * tempx1l + wxgll(i) * tempx2l)
            accel_elastic(2,iglob) = accel_elastic(2,iglob) - (wzgll(j) * tempz1l + wxgll(i) * tempz2l)
          endif
        enddo
      enddo
    endif ! AXISYM

    ! adds PML_BOUNDARY_CONDITIONS contribution
    if (PML_BOUNDARY_CONDITIONS) then
      if (ispec_is_PML(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. iglob_is_forced(iglob)) then
              accel_elastic(1,iglob) = accel_elastic(1,iglob) - accel_elastic_PML(1,i,j)
              accel_elastic(2,iglob) = accel_elastic(2,iglob) - accel_elastic_PML(2,i,j)
            endif
          enddo
        enddo
      endif
    endif

  enddo ! end of loop over all spectral elements

  contains

!---------------------------------------------------------------------------------------

  subroutine mxm_4comp_singleA(x1,x2,z1,z2,A,B,C)

! matrix x matrix multiplication, merging 4 loops for x1,x2 = A^t B^t and z1,z2 = A C^t
!
! index notation:
! general matrix multiplication: uij = (A B)ij = Aik Bkj
!                          here: xij = (A^t B^t)ij = Akj Bik = (B A)ij
!                                zij = (A C^t)ij = Aik Cjk
!
! original loops:
!
!      do j = 1,NGLLZ
!        do i = 1,NGLLX
!          ! derivative along x and along z
!          dux_dxi(i,j) = 0._CUSTOM_REAL
!          duz_dxi(i,j) = 0._CUSTOM_REAL
!          dux_dgamma(i,j) = 0._CUSTOM_REAL
!          duz_dgamma(i,j) = 0._CUSTOM_REAL
!
!          ! first double loop over GLL points to compute and store gradients
!          ! we can merge the two loops because NGLLX == NGLLZ
!          do k = 1,NGLLX
!            dux_dxi(i,j) = dux_dxi(i,j) + dummy_loc(1,k,j)*hprime_xx(i,k)
!            duz_dxi(i,j) = duz_dxi(i,j) + dummy_loc(2,k,j)*hprime_xx(i,k)
!            dux_dgamma(i,j) = dux_dgamma(i,j) + dummy_loc(1,i,k)*hprime_zz(j,k)
!            duz_dgamma(i,j) = duz_dgamma(i,j) + dummy_loc(2,i,k)*hprime_zz(j,k)
!          enddo
!        enddo
!      enddo

  use constants, only: NDIM,NGLLX,NGLLZ,CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: x1,x2,z1,z2

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ),intent(in) :: A
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(in) :: B,C

  ! local parameters
  integer :: i,j,k

  select case(NGLLX)
  case (5)
    do j = 1,5
      do i = 1,5
        ! loop unrolling
        x1(i,j) = A(1,1,j) * B(i,1) + A(1,2,j) * B(i,2) + A(1,3,j) * B(i,3) + A(1,4,j) * B(i,4) + A(1,5,j) * B(i,5)
        x2(i,j) = A(2,1,j) * B(i,1) + A(2,2,j) * B(i,2) + A(2,3,j) * B(i,3) + A(2,4,j) * B(i,4) + A(2,5,j) * B(i,5)

        z1(i,j) = A(1,i,1) * C(j,1) + A(1,i,2) * C(j,2) + A(1,i,3) * C(j,3) + A(1,i,4) * C(j,4) + A(1,i,5) * C(j,5)
        z2(i,j) = A(2,i,1) * C(j,1) + A(2,i,2) * C(j,2) + A(2,i,3) * C(j,3) + A(2,i,4) * C(j,4) + A(2,i,5) * C(j,5)
      enddo
    enddo

  case default
    do j = 1,NGLLZ
      do i = 1,NGLLX
        x1(i,j) = 0._CUSTOM_REAL
        x2(i,j) = 0._CUSTOM_REAL
        z1(i,j) = 0._CUSTOM_REAL
        z2(i,j) = 0._CUSTOM_REAL
        do k = 1,NGLLX
          x1(i,j) = x1(i,j) + A(1,k,j) * B(i,k)
          x2(i,j) = x2(i,j) + A(2,k,j) * B(i,k)

          z1(i,j) = z1(i,j) + A(1,i,k) * C(j,k)
          z2(i,j) = z2(i,j) + A(2,i,k) * C(j,k)
        enddo
      enddo
    enddo
  end select

  end subroutine mxm_4comp_singleA

  end subroutine compute_forces_viscoelastic

