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
!=====================================================================

  subroutine compute_kernels()

! computes adjoint sensitivity kernel contributions
!
! see e.g. Tromp et al. (2005) for elastic calculation
! and Morency et al. (2009) for poroelastic calculation

  use specfem_par, only: any_acoustic,any_elastic,any_poroelastic,it,NTSTEP_BETWEEN_COMPUTE_KERNELS, &
    APPROXIMATE_HESS_KL

  implicit none

  ! checks if anything to do
  if (mod(it,NTSTEP_BETWEEN_COMPUTE_KERNELS) /= 0) return

  ! acoustic simulations
  if (any_acoustic) then
    call compute_kernels_ac()
  endif

  ! elastic simulations
  if (any_elastic) then
    call compute_kernels_el()
  endif

  ! poro-elastic simulations
  if (any_poroelastic) then
    call compute_kernels_po()
  endif

  ! computes an approximative Hessian for preconditioning kernels
  if (APPROXIMATE_HESS_KL) then
    call compute_kernels_Hessian()
  endif

  end subroutine compute_kernels

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_el()

! elastic kernel calculations
! see e.g. Tromp et al. (2005)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,HALF,TWO,FOUR_THIRDS

  use specfem_par, only: ispec_is_elastic, & ! AXISYM,
                         nspec,ibool,accel_elastic,b_displ_elastic, &
                         P_SV,displ_elastic, &
                         ibool,hprime_xx,hprime_zz,xix,xiz,gammax,gammaz, &
                         GPU_MODE

  use specfem_par_gpu, only: Mesh_pointer

  use specfem_par, only: deltat,ispec_is_anisotropic, &
                         rho_kl,mu_kl,kappa_kl, &
                         c11_kl,c13_kl,c15_kl,c33_kl,c35_kl,c55_kl

  ! AXISYM case
  use specfem_par, only: AXISYM,is_on_the_axis,hprimeBar_xx,NGLJ

  implicit none

  !local variables
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: b_dux_dxi,b_dux_dgamma,b_duz_dxi,b_duz_dgamma
  real(kind=CUSTOM_REAL) :: b_dux_dxl,b_duz_dxl,b_dux_dzl,b_duz_dzl
  real(kind=CUSTOM_REAL) :: dsxx,dsxz,dszz
  real(kind=CUSTOM_REAL) :: b_dsxx,b_dsxz,b_dszz
  real(kind=CUSTOM_REAL) :: rho_k_loc,kappa_k_loc,mu_k_loc
  real(kind=CUSTOM_REAL) :: c11_k_loc,c13_k_loc,c15_k_loc,c33_k_loc,c35_k_loc,c55_k_loc
  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl

  ! elastic kernels
  if (.not. GPU_MODE) then
    ! updates kernels on CPU
    do ispec = 1,nspec
      ! only elastic elements
      if (.not. ispec_is_elastic(ispec)) cycle

      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! derivative along x and along z
          dux_dxi = 0._CUSTOM_REAL
          duz_dxi = 0._CUSTOM_REAL
          dux_dgamma = 0._CUSTOM_REAL
          duz_dgamma = 0._CUSTOM_REAL

          b_dux_dxi = 0._CUSTOM_REAL
          b_duz_dxi = 0._CUSTOM_REAL
          b_dux_dgamma = 0._CUSTOM_REAL
          b_duz_dgamma = 0._CUSTOM_REAL

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)

            b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            b_duz_dxi = b_duz_dxi + b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            b_duz_dgamma = b_duz_dgamma + b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
          enddo

          ! AXISYM case overwrites dux_dxi and duz_dxi
          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              ! derivative along x and along z
              dux_dxi = 0._CUSTOM_REAL
              duz_dxi = 0._CUSTOM_REAL
              b_dux_dxi = 0._CUSTOM_REAL
              b_duz_dxi = 0._CUSTOM_REAL
              do k = 1,NGLJ
                dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                b_duz_dxi = b_duz_dxi + b_displ_elastic(2,ibool(k,j,ispec))*hprimeBar_xx(i,k)
              enddo
            endif
          endif

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

          ! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl
          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
          b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl
          b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
          b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl

          ! AXISYM case overwrite duz_dxl
          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              ! d_uz/dr=0 on the axis
              if (i == 1) then
                duz_dxl = 0._CUSTOM_REAL
                b_duz_dxl = 0._CUSTOM_REAL
              endif
            endif
          endif

          ! isotropic kernel contributions
          if (P_SV) then
            ! P-SV waves
            !   K_kappa = - kappa int_0^T [ div(s^adj) div(s) ] dt
            !   K_mu    = - 2 mu  int_0^T [ D^adj : D ] dt               (where D = 1/2[grad(s) + grad(s)^T] - 1/3 div(s) I )
            !
            ! the factors (- kappa dt) and (- 2 mu dt) will be added later in save_adjoint_kernels() routine
            dsxx =  dux_dxl
            dsxz = HALF * (duz_dxl + dux_dzl)
            dszz =  duz_dzl

            b_dsxx =  b_dux_dxl
            b_dsxz = HALF * (b_duz_dxl + b_dux_dzl)
            b_dszz =  b_duz_dzl

            kappa_k_loc = (dsxx + dszz) * (b_dsxx + b_dszz)
            mu_k_loc = dsxx * b_dsxx + dszz * b_dszz + &
                          2._CUSTOM_REAL * dsxz * b_dsxz - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k_loc
          else
            ! SH (membrane) waves
            !   K_kappa = 0
            !   K_mu    = - 2 mu  int_0^T [ D^adj : D ] dt
            ! where
            !   D = 1/2[grad(s) + grad(s)^T] - 1/3 div(s) I       deviatoric strain
            ! and the factor (- 2 mu dt) will be added later in save_adjoint_kernels() routine
            !
            ! note:
            !   D = [    dux_dx         1/2(dux_dy+duy_dx) 1/2(dux_dz+duz_dx)        [ 1/3 (dux_dx + duy_dy + duz_dz)
            !        1/2(duy_dx+dux_dy)     duy_dy         1/2(duy_dz+duz_dy)     -         1/3 (dux_dx + duy_dy + duz_dz)
            !        1/2(duz_dx+dux_dz) 1/2(duz_dy+duy_dz)     duz_dz   ]                       1/3 (dux_dx + duy_dy + duz_dx) ]
            !
            ! SH-waves: plane strain assumption ux==uz==0 and d/dy==0
            !   D = [   0             1/2 duy_dx       0                   [
            !          1/2 duy_dx       0             1/2 duy_dz       -             0
            !           0             1/2 duy_dz       0                                     ]
            !
            !  D^adj : D = sum_i sum_j D^adj_ij * D_ij
            !            = 1/2duy_dx^adj * 1/2duy_dx + 1/2duy_dx^adj * 1/2duy_dx
            !                + 1/2duy_dz^adj * 1/2duy_dz + 1/2duy_dz^adj * 1/2duy_dz
            !            = 1/2 ( duy_dx^adj * duy_dx) + 1/2 (duy_dz^adj * duy_dz)
            !            = 1/2 ( duy_dx^adj * duy_dx + duy_dz^adj * duy_dz )
            kappa_k_loc = 0.0_CUSTOM_REAL
            mu_k_loc = HALF * (dux_dxl * b_dux_dxl + dux_dzl * b_dux_dzl)
          endif


          ! adding contributions to sensitivity kernel values (on local basis (i,j,ispec))
          iglob = ibool(i,j,ispec)
          rho_k_loc =  accel_elastic(1,iglob)*b_displ_elastic(1,iglob) + accel_elastic(2,iglob)*b_displ_elastic(2,iglob)

          ! note: we will add the minus sign to the kernels when multiplying with the material properties
          !       at the end of the time looping in save_adjoint_kernels.f90

          ! isotropic kernels
          ! for parameterization (rho,mu,kappa): "primary" kernels
          ! density kernel
          rho_kl(i,j,ispec) = rho_kl(i,j,ispec) + rho_k_loc
          ! shear modulus kernel
          mu_kl(i,j,ispec) =  mu_kl(i,j,ispec) +  mu_k_loc
          ! bulk modulus kernel
          kappa_kl(i,j,ispec) = kappa_kl(i,j,ispec) +  kappa_k_loc

          ! Voigt kernels, e.g., see Sieminski, 2007a,b
          if (ispec_is_anisotropic(ispec)) then
            c11_k_loc = dux_dxl*b_dux_dxl
            c13_k_loc = dux_dxl*b_duz_dzl + duz_dzl*b_dux_dxl
            c15_k_loc = 2*(dux_dxl*HALF*(b_dux_dzl+b_duz_dxl)+&
                           HALF*(dux_dzl+duz_dxl)*b_dux_dxl)
            c33_k_loc = duz_dzl*b_duz_dzl
            c35_k_loc = 2*(duz_dzl*HALF*(b_dux_dzl+b_duz_dxl)+&
                           HALF*(dux_dzl+duz_dxl)*b_duz_dzl)
            c55_k_loc = 4*HALF*(dux_dzl+duz_dxl)*HALF*(b_dux_dzl+b_duz_dxl)

            ! "primary" kernels
            !rho_kl(i,j,ispec) = rho_kl(i,j,ispec) + rho_k_loc  ! taken from above
            c11_kl(i,j,ispec) = c11_kl(i,j,ispec) + c11_k_loc
            c13_kl(i,j,ispec) = c13_kl(i,j,ispec) + c13_k_loc
            c15_kl(i,j,ispec) = c15_kl(i,j,ispec) + c15_k_loc
            c33_kl(i,j,ispec) = c33_kl(i,j,ispec) + c33_k_loc
            c35_kl(i,j,ispec) = c35_kl(i,j,ispec) + c35_k_loc
            c55_kl(i,j,ispec) = c55_kl(i,j,ispec) + c55_k_loc
          endif
        enddo
      enddo
    enddo

  else
    ! updates kernels on GPU
    call compute_kernels_elastic_cuda(Mesh_pointer,deltat)
  endif

  end subroutine compute_kernels_el

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_ac()

! acoustic kernel calculations
! see e.g. Tromp et al. (2005)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,ZERO,HALF,TWO

  use specfem_par, only: nspec,ispec_is_acoustic,ibool, &
                         rhostore,kappastore,deltat, &
                         hprime_xx,hprime_zz,xix,xiz,gammax,gammaz, &
                         potential_acoustic,b_potential_acoustic,potential_dot_dot_acoustic, &
                         NTSTEP_BETWEEN_COMPUTE_KERNELS, &
                         rho_ac_kl,kappa_ac_kl,rhop_ac_kl,alpha_ac_kl,GPU_MODE

  use specfem_par_gpu, only: Mesh_pointer

  ! AXISYM case
  use specfem_par, only: AXISYM,is_on_the_axis,hprimeBar_xx,NGLJ

  implicit none

  !local variables
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,b_tempx1l,b_tempx2l
  real(kind=CUSTOM_REAL) :: rhol,kappal
  real(kind=CUSTOM_REAL) :: b_displ_loc(2),accel_loc(2)
  double precision :: xixl,xizl,gammaxl,gammazl

  if (.not. GPU_MODE) then
    ! kernels on CPU
    do ispec = 1, nspec
      ! acoustic kernels
      if (ispec_is_acoustic(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)

            rhol = rhostore(i,j,ispec)
            kappal = kappastore(i,j,ispec)

            ! calcul the displacement by computing the gradient of potential / rho
            ! and calcul the acceleration by computing the gradient of potential_dot_dot / rho
            !
            ! note: we use the displ potential_acoustic for the adjoint wavefield, and not the accel potential_dot_dot_acoustic,
            !       to match the rho-kernel expressions from Luo et al. (2013)
            tempx1l = ZERO
            tempx2l = ZERO
            b_tempx1l = ZERO
            b_tempx2l = ZERO

            !  bb_tempx1l = ZERO
            !  bb_tempx2l = ZERO
            do k = 1,NGLLX
              ! derivative along x
              !tempx1l = tempx1l + potential_dot_dot_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
              tempx1l = tempx1l + potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k) !!! YANGL
              b_tempx1l = b_tempx1l + b_potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
              ! derivative along z
              !tempx2l = tempx2l + potential_dot_dot_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
              tempx2l = tempx2l + potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k) !!! YANGL
              b_tempx2l = b_tempx2l + b_potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
            enddo

            ! AXISYM case overwrites dux_dxi
            if (AXISYM) then
              if (is_on_the_axis(ispec)) then
                ! derivative along x
                tempx1l = ZERO
                b_tempx1l = ZERO
                do k = 1,NGLJ
                  tempx1l = tempx1l + potential_acoustic(ibool(k,j,ispec)) * hprimeBar_xx(i,k)
                  b_tempx1l = b_tempx1l + b_potential_acoustic(ibool(k,j,ispec)) * hprimeBar_xx(i,k)
                enddo
              endif
            endif

            xixl = xix(i,j,ispec)
            xizl = xiz(i,j,ispec)
            gammaxl = gammax(i,j,ispec)
            gammazl = gammaz(i,j,ispec)

            ! derivatives of potential
            ! warning : the variable is named accel_loc but it is displ that is computed
            accel_loc(1) = (tempx1l*xixl + tempx2l*gammaxl)
            accel_loc(2) = (tempx1l*xizl + tempx2l*gammazl)

            b_displ_loc(1) = (b_tempx1l*xixl + b_tempx2l*gammaxl)
            b_displ_loc(2) = (b_tempx1l*xizl + b_tempx2l*gammazl)

            ! AXISYM case overwrite dux_dxl
            if (AXISYM) then
              if (is_on_the_axis(ispec)) then
                ! dchi/dr=rho * u_r=0 on the axis
                if (i == 1) then
                  accel_loc(1) = 0._CUSTOM_REAL
                  b_displ_loc(1) = 0._CUSTOM_REAL
                endif
              endif
            endif

            ! acoustic kernel integration
            ! YANGL
            !!!! old expression (from elastic kernels)
            !!!rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) - rhol  * &
            !!!      dot_product(accel_loc(:),b_displ_loc(:)) * deltat
            !!!kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) - kappal * &
            !!!      potential_dot_dot_acoustic(iglob) / kappal * &
            !!!      b_potential_dot_dot_acoustic(iglob) / kappal
            !!!      * deltat

            ! new expression (from PDE-constrained optimization, coupling terms changed as well)
            ! note: from Luo et al. (2013), kernels are given for absolute perturbations d(1/rho) and d(1/kappa)
            !       and only change from Peter et al. (2011) for acoustic kernels:
            !         K_kappa = - int_0^T [ phi^adj \partial_t^2 phi ] dt     see (A-27)
            !         K_rho   = - int_0^T [ grad(phi^adj) * grad(phi) ] dt        (A-28)
            !
            !       since we output relative perturbation kernels for elastic domains, we make here also use of the variation
            !         d(1/rho) = - 1/rho^2 d(rho) = - 1/rho d(rho)/rho = - 1/rho dln(rho)
            !
            !       to obtain relative kernels, we start from eq. (A-24)
            !         dChi = int_V [ d(1/kappa) K_kappa + d(1/rho) K_rho ] d^3x (+ elastic kernels)
            !              = int_V [ (-1/kappa K_kappa) dln(kappa) + (- 1/rho K_rho) dln(rho)
            !
            !       and see that relative perturbation kernels are given by
            !          \tilde{K}_kappa = - 1/kappa K_kappa
            !                          = + 1/kappa int_0^T [ phi^adj \partial_t^2 phi ] dt
            !                          = + 1/kappa int_0^T [ \partial_t^2 phi^adj phi ] dt              (equivalent)
            !
            !          \tilde{K}_rho   = - 1/rho   K_rho
            !                          = + 1/rho   int_0^T [ grad(phi^adj) * grad(phi) ] dt
            !                          = + rho     int_0^T [ 1/rho grad(phi^adj) * 1/rho grad(phi) ] dt   (equivalent)
            !
            ! density kernel
            rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) &
                                 + 1.0_CUSTOM_REAL/rhol * dot_product(accel_loc(:),b_displ_loc(:)) &
                                        * (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS)

            ! kappa (bulk modulus) kernel
            kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) &
                                   + 1.0_CUSTOM_REAL/kappal * potential_dot_dot_acoustic(iglob) * b_potential_acoustic(iglob) &
                                        * (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS)

            ! rho prime kernel
            rhop_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + kappa_ac_kl(i,j,ispec)

            ! alpha kernel
            alpha_ac_kl(i,j,ispec) = TWO *  kappa_ac_kl(i,j,ispec)
          enddo
        enddo
      endif
    enddo
  else
    ! on GPU
    call compute_kernels_acoustic_cuda(Mesh_pointer,deltat)
  endif

  end subroutine compute_kernels_ac

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_po()

! kernel calculations
! see e.g. Morency et al. (2009)

  use constants, only: CUSTOM_REAL,FOUR_THIRDS,NGLLX,NGLLZ,TWO,HALF

  use specfem_par, only: nspec,ispec_is_poroelastic,ibool,deltat, &
                         phistore,tortstore,kappaarraystore,mufr_store,rhoarraystore,permstore,etastore, &
                         accels_poroelastic,accelw_poroelastic,velocw_poroelastic, &
                         b_displs_poroelastic,b_displw_poroelastic, &
                         epsilondev_s,b_epsilondev_s, &
                         epsilondev_w,b_epsilondev_w, &
                         rhot_kl,rhof_kl,sm_kl,eta_kl,B_kl,C_kl,M_kl, &
                         mufr_kl,rhob_kl,rhofb_kl, &
                         mufrb_kl,phi_kl,rhobb_kl,rhofbb_kl,phib_kl,cpI_kl,cpII_kl,cs_kl,ratio_kl, &
                         GPU_MODE,NTSTEP_BETWEEN_COMPUTE_KERNELS
  implicit none

  !local variables
  integer :: i,j,ispec,iglob
  real(kind=CUSTOM_REAL) :: rholb,dd1
  real(kind=CUSTOM_REAL) :: ratio
  real(kind=CUSTOM_REAL) :: dsxx,dsxz,dszz,dszx_xz
  real(kind=CUSTOM_REAL) :: b_dsxx,b_dsxz,b_dszz,b_dszx_xz
  real(kind=CUSTOM_REAL) :: dwxx,dwzz,b_dwxx,b_dwzz
  real(kind=CUSTOM_REAL) :: rhot_k_loc,rhof_k_loc,sm_k_loc,eta_k_loc,B_k_loc,C_k_loc,M_k_loc,mufr_k_loc

  ! to evaluate cpI, cpII, and cs, and rI (poroelastic medium)
  double precision :: phi,tort,kappa_s,kappa_f,kappa_fr,mu_fr
  double precision :: rho_s,rho_f,rho_bar,eta_f
  double precision :: D_biot,H_biot,C_biot,M_biot
  double precision :: B_biot
  double precision :: perm_xx
  double precision :: afactor,bfactor,cfactor
  double precision :: gamma1,gamma2,gamma3,gamma4
  double precision :: cpIsquare,cpIIsquare,cssquare
  ! helper parameters
  double precision :: phi_over_tort,phi_over_tort_ratio,ratio_square

  ! safety check
  if (GPU_MODE) call stop_the_code('Error poroelastic kernels not implemented on GPUs yet')

  ! kernels on local nodes
  do ispec = 1, nspec
    ! only poroelastic elements
    if (.not. ispec_is_poroelastic(ispec)) cycle

    do j = 1, NGLLZ
      do i = 1, NGLLX
        iglob = ibool(i,j,ispec)

        ! kernel contributions
        rhot_k_loc = accels_poroelastic(1,iglob) * b_displs_poroelastic(1,iglob) + &
                        accels_poroelastic(2,iglob) * b_displs_poroelastic(2,iglob)

        rhof_k_loc = accelw_poroelastic(1,iglob) * b_displs_poroelastic(1,iglob) + &
                        accelw_poroelastic(2,iglob) * b_displs_poroelastic(2,iglob) + &
                        accels_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                        accels_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)

        sm_k_loc = accelw_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                       accelw_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)

        eta_k_loc = velocw_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                       velocw_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)

        ! gets poroelastic material
        phi = phistore(i,j,ispec)
        tort = tortstore(i,j,ispec)
        kappa_s = kappaarraystore(1,i,j,ispec)
        kappa_f = kappaarraystore(2,i,j,ispec)
        kappa_fr = kappaarraystore(3,i,j,ispec)
        mu_fr = mufr_store(i,j,ispec)

        ! Biot coefficients for the input phi
        call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

        B_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr

        ! permeability
        perm_xx = permstore(1,i,j,ispec)
        eta_f = etastore(i,j,ispec)

        rho_s = rhoarraystore(1,i,j,ispec)
        rho_f = rhoarraystore(2,i,j,ispec)
        rho_bar = (1.d0 - phi)*rho_s + phi * rho_f

        ! Approximated velocities (no viscous dissipation)
        phi_over_tort = phi / tort

        afactor = rho_bar - phi_over_tort * rho_f
        bfactor = H_biot + phi_over_tort * rho_bar/rho_f * M_biot - TWO * phi_over_tort * C_biot
        cfactor = phi_over_tort / rho_f * (H_biot * M_biot - C_biot * C_biot)

        cpIsquare = (bfactor + sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)
        cpIIsquare = (bfactor - sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)
        cssquare = mu_fr/afactor

        ! Approximated ratio r = amplitude "w" field/amplitude "s" field (no viscous dissipation)
        ! used later for wavespeed kernels calculation, which are presently implemented for inviscid case,
        ! contrary to primary and density-normalized kernels, which are consistent with viscous fluid case.
        gamma1 = H_biot - phi_over_tort * C_biot
        gamma2 = C_biot - phi_over_tort * M_biot
        gamma3 = phi_over_tort * ( M_biot * (afactor/rho_f + phi_over_tort) - C_biot)
        gamma4 = phi_over_tort * ( C_biot * (afactor/rho_f + phi_over_tort) - H_biot)

        ratio = HALF*(gamma1 - gamma3)/gamma4 + HALF*sqrt((gamma1-gamma3)**2/gamma4**2 + 4.d0 * gamma2/gamma4)

        ! helper
        phi_over_tort_ratio = phi_over_tort * ratio
        ratio_square = ratio**2

        rhot_kl(i,j,ispec) = rhot_kl(i,j,ispec) - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * rho_bar * rhot_k_loc
        rhof_kl(i,j,ispec) = rhof_kl(i,j,ispec) - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * rho_f * rhof_k_loc
        sm_kl(i,j,ispec) = sm_kl(i,j,ispec) - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * rho_f * tort/phi * sm_k_loc

        !at the moment works with constant permeability
        eta_kl(i,j,ispec) = eta_kl(i,j,ispec) - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * eta_f/perm_xx * eta_k_loc

        ! for B_k & mufr_k
        dsxx = epsilondev_s(1,i,j,ispec) ! dux_dxl
        dszz = epsilondev_s(2,i,j,ispec) ! duz_dzl
        dsxz = epsilondev_s(3,i,j,ispec) ! dux_dzl
        dszx_xz = epsilondev_s(4,i,j,ispec) ! 0.5_CUSTOM_REAL * (duz_dxl + dux_dzl)

        b_dsxx = b_epsilondev_s(1,i,j,ispec) ! b_dux_dxl
        b_dszz = b_epsilondev_s(2,i,j,ispec) ! b_duz_dzl
        b_dsxz = b_epsilondev_s(3,i,j,ispec) ! b_dux_dzl
        b_dszx_xz = b_epsilondev_s(4,i,j,ispec) ! 0.5_CUSTOM_REAL * (b_duz_dxl + b_dux_dzl)

        B_k_loc = (dsxx + dszz) *  (b_dsxx + b_dszz) * (H_biot - FOUR_THIRDS * mu_fr)

        mufr_k_loc = (dsxx * b_dsxx + dszz * b_dszz + &
                        2._CUSTOM_REAL * dszx_xz * b_dszx_xz - &
                        1._CUSTOM_REAL/3._CUSTOM_REAL * (dsxx + dszz) * (b_dsxx + b_dszz) ) * mu_fr

        ! from older compute_forces_poro_solid ...
        !  iglob = ibool(i,j,ispec)
        !  dsxx =  dux_dxl
        !  dsxz = HALF * (duz_dxl + dux_dzl)
        !  dszz =  duz_dzl
        !
        !  b_dsxx =  b_dux_dxl
        !  b_dsxz = HALF * (b_duz_dxl + b_dux_dzl)
        !  b_dszz =  b_duz_dzl
        !
        !  B_k(iglob) = (dux_dxl + duz_dzl) *  (b_dux_dxl + b_duz_dzl) * (H_biot - FOUR_THIRDS * mu_fr)
        !  mufr_k(iglob) = (dsxx * b_dsxx + dszz * b_dszz + &
        !                  2._CUSTOM_REAL * dsxz * b_dsxz - &
        !                  1._CUSTOM_REAL/3._CUSTOM_REAL * (dux_dxl + duz_dzl) * (b_dux_dxl + b_duz_dzl) ) * mu_fr

        ! for C_k & M_k
        dsxx = epsilondev_w(1,i,j,ispec) ! dux_dxl
        dszz = epsilondev_w(2,i,j,ispec) ! duz_dzl
        dwxx = epsilondev_w(3,i,j,ispec) ! dwx_dxl
        dwzz = epsilondev_w(4,i,j,ispec) ! dwz_dzl

        b_dsxx = b_epsilondev_w(1,i,j,ispec) ! b_dux_dxl
        b_dszz = b_epsilondev_w(2,i,j,ispec) ! b_duz_dzl
        b_dwxx = b_epsilondev_w(3,i,j,ispec) ! b_dwx_dxl
        b_dwzz = b_epsilondev_w(4,i,j,ispec) ! b_dwz_dzl

        C_k_loc = ((dsxx + dszz)*(b_dwxx + b_dwzz) + (dwxx + dwzz)*(b_dsxx + b_dszz)) * C_biot
        M_k_loc = (dwxx + dwzz)*(b_dwxx + b_dwzz) * M_biot

        ! from older compute_forces_poro_fluid ...
        !C_k(iglob) =  ((dux_dxl + duz_dzl) *  (b_dwx_dxl + b_dwz_dzl) + &
        !                (dwx_dxl + dwz_dzl) *  (b_dux_dxl + b_duz_dzl)) * C_biot
        !M_k(iglob) = (dwx_dxl + dwz_dzl) *  (b_dwx_dxl + b_dwz_dzl) * M_biot

        B_kl(i,j,ispec) = B_kl(i,j,ispec) - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * B_k_loc
        C_kl(i,j,ispec) = C_kl(i,j,ispec) - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * C_k_loc
        M_kl(i,j,ispec) = M_kl(i,j,ispec) - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * M_k_loc

        mufr_kl(i,j,ispec) = mufr_kl(i,j,ispec) - TWO * (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * mufr_k_loc

        ! density kernels
        rholb = rho_bar - phi_over_tort * rho_f
        rhob_kl(i,j,ispec) = rhot_kl(i,j,ispec) + B_kl(i,j,ispec) + mufr_kl(i,j,ispec)
        rhofb_kl(i,j,ispec) = rhof_kl(i,j,ispec) + C_kl(i,j,ispec) + M_kl(i,j,ispec) + sm_kl(i,j,ispec)

        mufrb_kl(i,j,ispec) = mufr_kl(i,j,ispec)
        phi_kl(i,j,ispec) = - sm_kl(i,j,ispec) - M_kl(i,j,ispec)

        ! wave speed kernels
        dd1 = (1._CUSTOM_REAL+rholb/rho_f) * ratio_square  + 2._CUSTOM_REAL*ratio + tort/phi

        rhobb_kl(i,j,ispec) = rhob_kl(i,j,ispec) &
              - phi_over_tort * rho_f/B_biot * &
                (cpIIsquare + (cpIsquare - cpIIsquare) * ( (phi_over_tort_ratio + 1._CUSTOM_REAL)/dd1 + &
                (rho_bar**2 * ratio_square / rho_f**2 * (phi_over_tort_ratio + 1._CUSTOM_REAL)*(phi_over_tort_ratio + &
                phi_over_tort * (1._CUSTOM_REAL + rho_f/rho_bar) - 1._CUSTOM_REAL) )/dd1**2 ) - FOUR_THIRDS * cssquare ) &
                * B_kl(i,j,ispec) &
              - rho_bar * ratio_square / M_biot * (cpIsquare - cpIIsquare) * &
                (phi_over_tort_ratio + 1._CUSTOM_REAL)**2 / dd1**2 &
                * M_kl(i,j,ispec) &
              + rho_bar * ratio / C_biot * (cpIsquare - cpIIsquare) * ( (phi_over_tort_ratio + 1._CUSTOM_REAL)/dd1 - &
                phi_over_tort_ratio * (phi_over_tort_ratio + 1._CUSTOM_REAL) * &
                (1._CUSTOM_REAL + rho_bar * ratio / rho_f)/dd1**2) &
                * C_kl(i,j,ispec) &
              + phi_over_tort * rho_f * cssquare / mu_fr &
                * mufrb_kl(i,j,ispec)

        rhofbb_kl(i,j,ispec) = rhofb_kl(i,j,ispec) &
              + phi_over_tort * rho_f/ B_biot * &
                (cpIIsquare + (cpIsquare - cpIIsquare) * ( (phi_over_tort_ratio  + 1._CUSTOM_REAL)/dd1 + &
                (rho_bar**2 * ratio_square / rho_f**2 * (phi_over_tort_ratio + 1._CUSTOM_REAL)*(phi_over_tort_ratio + &
                phi_over_tort * (1._CUSTOM_REAL + rho_f/rho_bar) - 1._CUSTOM_REAL) )/dd1**2 ) - FOUR_THIRDS * cssquare ) &
                * B_kl(i,j,ispec) &
              + rho_bar * ratio_square / M_biot * (cpIsquare - cpIIsquare) * &
                (phi_over_tort_ratio + 1._CUSTOM_REAL)**2 / dd1**2 &
                * M_kl(i,j,ispec) &
              - rho_bar * ratio / C_biot * (cpIsquare - cpIIsquare)* ( (phi_over_tort_ratio + 1._CUSTOM_REAL)/dd1 - &
                phi_over_tort_ratio * (phi_over_tort_ratio + 1._CUSTOM_REAL) * &
                (1._CUSTOM_REAL + rho_bar * ratio / rho_f)/dd1**2) &
                * C_kl(i,j,ispec) &
              - phi_over_tort * rho_f * cssquare/ mu_fr &
                * mufrb_kl(i,j,ispec)

        phib_kl(i,j,ispec) = phi_kl(i,j,ispec) &
              - phi_over_tort * rho_bar / B_biot * ( cpIsquare - rho_f/rho_bar * cpIIsquare - &
                (cpIsquare - cpIIsquare) * ( (TWO * ratio_square * phi_over_tort + (1._CUSTOM_REAL + rho_f/rho_bar) * &
                (TWO * ratio * phi_over_tort + 1._CUSTOM_REAL))/dd1 + (phi_over_tort_ratio + 1._CUSTOM_REAL) * &
                (phi_over_tort_ratio + phi_over_tort * (1._CUSTOM_REAL + rho_f/rho_bar) - 1._CUSTOM_REAL) * ((1._CUSTOM_REAL+ &
                rho_bar/rho_f - TWO * phi_over_tort) * ratio_square + TWO * ratio)/dd1**2 ) - &
                FOUR_THIRDS * rho_f * cssquare/rho_bar ) &
                * B_kl(i,j,ispec) &
              + rho_f/M_biot * (cpIsquare-cpIIsquare) * &
                (TWO * ratio * (phi_over_tort_ratio + 1._CUSTOM_REAL)/dd1 - (phi_over_tort_ratio + 1._CUSTOM_REAL)**2 * &
                ((1._CUSTOM_REAL + rho_bar/rho_f - TWO * phi_over_tort) * ratio_square + TWO * ratio) / dd1**2) &
                * M_kl(i,j,ispec) &
              + phi_over_tort * rho_f / C_biot * (cpIsquare - cpIIsquare) * ratio * (&
                (1._CUSTOM_REAL + rho_f/rho_bar * ratio)/dd1 - (phi_over_tort_ratio + 1._CUSTOM_REAL) * &
                (1._CUSTOM_REAL + rho_bar/rho_f * ratio) * &
                ((1._CUSTOM_REAL + rho_bar/rho_f - TWO * phi_over_tort) * ratio + TWO) / dd1**2) &
                * C_kl(i,j,ispec) &
              - phi_over_tort * rho_f * cssquare / mu_fr &
                * mufrb_kl(i,j,ispec)

        ! wavespeed kernels
        cpI_kl(i,j,ispec) = 2._CUSTOM_REAL * cpIsquare * rho_bar / B_biot * ( &
                1._CUSTOM_REAL - phi_over_tort + (phi_over_tort_ratio + 1._CUSTOM_REAL) * (phi_over_tort_ratio + &
                phi_over_tort * (1._CUSTOM_REAL + rho_f/rho_bar) - 1._CUSTOM_REAL) / dd1 ) &
                * B_kl(i,j,ispec) &
              + 2._CUSTOM_REAL * cpIsquare * rho_f * tort / (phi*M_biot) * &
                (phi_over_tort_ratio + 1._CUSTOM_REAL)**2 /dd1 &
                * M_kl(i,j,ispec) &
              + 2._CUSTOM_REAL * cpIsquare * rho_f / C_biot * &
                (phi_over_tort_ratio + 1._CUSTOM_REAL) * (1._CUSTOM_REAL + rho_bar/rho_f * ratio) / dd1 &
                * C_kl(i,j,ispec)
        cpII_kl(i,j,ispec) = 2._CUSTOM_REAL * cpIIsquare * rho_bar / B_biot * ( &
                phi_over_tort * rho_f / rho_bar - (phi_over_tort_ratio + 1._CUSTOM_REAL) * (phi_over_tort_ratio + &
                phi_over_tort * (1._CUSTOM_REAL + rho_f/rho_bar) - 1._CUSTOM_REAL) / dd1  ) &
                * B_kl(i,j,ispec) &
              + 2._CUSTOM_REAL * cpIIsquare * rho_f * tort / (phi*M_biot) * &
                (1._CUSTOM_REAL - (phi_over_tort_ratio + 1._CUSTOM_REAL)**2 / dd1 ) &
                * M_kl(i,j,ispec) &
              + 2._CUSTOM_REAL * cpIIsquare * rho_f / C_biot * &
                (1._CUSTOM_REAL - (phi_over_tort_ratio + 1._CUSTOM_REAL) * (1._CUSTOM_REAL + rho_bar / rho_f * ratio) / dd1 ) &
                * C_kl(i,j,ispec)

        cs_kl(i,j,ispec) = - 8._CUSTOM_REAL / 3._CUSTOM_REAL * cssquare * rho_bar / B_biot * &
                (1._CUSTOM_REAL - phi_over_tort * rho_f/rho_bar) &
                * B_kl(i,j,ispec) &
              + 2._CUSTOM_REAL * (rho_bar - rho_f * phi_over_tort) / mu_fr * cssquare &
                * mufrb_kl(i,j,ispec)

        ratio_kl(i,j,ispec) = ratio * rho_bar * phi_over_tort / B_biot * (cpIsquare - cpIIsquare) * &
                (phi_over_tort * (2._CUSTOM_REAL * ratio + 1._CUSTOM_REAL + rho_f / rho_bar) / dd1 - &
                (phi_over_tort_ratio + 1._CUSTOM_REAL) * &
                (phi_over_tort_ratio + phi_over_tort * ( 1._CUSTOM_REAL + rho_f / rho_bar) - 1._CUSTOM_REAL) * &
                (2._CUSTOM_REAL * ratio * (1._CUSTOM_REAL + rho_bar / rho_f - phi_over_tort) + 2._CUSTOM_REAL) / dd1**2  ) &
                * B_kl(i,j,ispec) &
              + ratio * rho_f * tort / (phi * M_biot) * (cpIsquare - cpIIsquare) * 2._CUSTOM_REAL * phi_over_tort * &
                ( (phi_over_tort_ratio + 1._CUSTOM_REAL) / dd1 - (phi_over_tort_ratio + 1._CUSTOM_REAL)**2 * &
                ((1._CUSTOM_REAL + rho_bar / rho_f -phi_over_tort) * ratio + 1._CUSTOM_REAL) / dd1**2 ) &
                * M_kl(i,j,ispec) &
              + ratio * rho_f / C_biot * (cpIsquare - cpIIsquare) * &
                ( (2._CUSTOM_REAL * phi_over_tort_ratio * rho_bar / rho_f + phi_over_tort + rho_bar / rho_f) / dd1 - &
                2._CUSTOM_REAL * phi_over_tort * (phi_over_tort_ratio + 1._CUSTOM_REAL) * &
                (1._CUSTOM_REAL + rho_bar / rho_f * ratio) * &
                ((1._CUSTOM_REAL + rho_bar / rho_f - phi_over_tort) * ratio + 1._CUSTOM_REAL) / dd1**2 ) &
                * C_kl(i,j,ispec)
      enddo
    enddo
  enddo

  end subroutine compute_kernels_po

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_Hessian()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,ZERO

  use specfem_par, only: nglob,nspec,ibool,ispec_is_acoustic,ispec_is_elastic, &
                         rhostore, &
                         hprime_xx,hprime_zz,xix,xiz,gammax,gammaz, &
                         any_elastic,any_acoustic, &
                         accel_elastic,b_accel_elastic, &
                         potential_acoustic,b_potential_acoustic, &
                         rhorho_el_Hessian_final1,rhorho_el_Hessian_final2, &
                         rhorho_ac_Hessian_final1,rhorho_ac_Hessian_final2, &
                         deltat,GPU_MODE,NTSTEP_BETWEEN_COMPUTE_KERNELS, &
                         APPROXIMATE_HESS_KL

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  !local variables
  real(kind=CUSTOM_REAL), dimension(nglob) :: rhorho_el_Hessian_temp1, rhorho_el_Hessian_temp2
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,b_tempx1l,b_tempx2l
  real(kind=CUSTOM_REAL) :: rhol
  real(kind=CUSTOM_REAL) :: b_accel_loc(2),accel_loc(2)
  double precision :: xixl,xizl,gammaxl,gammazl

  integer :: i,j,k,ispec,iglob

  ! checks if anything to do
  if (.not. APPROXIMATE_HESS_KL) return

  if (.not. GPU_MODE) then
    ! elastic domains
    if (any_elastic) then
      ! approximate Hessians
      ! pre-computes contributions on global points
      do iglob = 1,nglob
        rhorho_el_Hessian_temp1(iglob) = b_accel_elastic(1,iglob)*b_accel_elastic(1,iglob) + &
                                         b_accel_elastic(2,iglob)*b_accel_elastic(2,iglob)
        rhorho_el_Hessian_temp2(iglob) = accel_elastic(1,iglob)*b_accel_elastic(1,iglob) + &
                                         accel_elastic(2,iglob)*b_accel_elastic(2,iglob)
      enddo

      ! on local GLL basis
      do ispec = 1, nspec
        if (ispec_is_elastic(ispec)) then
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob = ibool(i,j,ispec)
              rhorho_el_Hessian_final1(i,j,ispec) = rhorho_el_Hessian_final1(i,j,ispec) + &
                                                    rhorho_el_Hessian_temp1(iglob) * (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS)
              rhorho_el_Hessian_final2(i,j,ispec) = rhorho_el_Hessian_final2(i,j,ispec) + &
                                                    rhorho_el_Hessian_temp2(iglob) * (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS)
            enddo
          enddo
        endif
      enddo
    endif

    ! acoustic domains
    if (any_acoustic) then
      ! on local GLL basis
      do ispec = 1,nspec
        if (ispec_is_acoustic(ispec)) then
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob = ibool(i,j,ispec)

              ! there are 2 different ways to compute an approximate Hessian,
              ! either using the potentials directly or accelerations - check which one is better...

              ! way 1:
              ! expressions using potentials
              !rhorho_ac_Hessian_final1(i,j,ispec) = rhorho_ac_Hessian_final1(i,j,ispec) &
              !                                    + potential_dot_dot_acoustic(iglob) * potential_dot_dot_acoustic(iglob) &
              !                                      * (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS)
              !
              !rhorho_ac_Hessian_final2(i,j,ispec) = rhorho_ac_Hessian_final2(i,j,ispec) &
              !                                    + potential_dot_dot_acoustic(iglob) * b_potential_dot_dot_acoustic(iglob) &
              !                                      * (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS)

              ! way 2:
              ! expressions using accelerations
              ! acceleration is defined as the gradient of potential_dot_dot / rho
              !
              ! note: we also use the displ potential for the adjoint wavefield, and not the accel potential_dot_dot,
              ! similar to the rho-kernel expressions from Luo et al. (2013)
              rhol = rhostore(i,j,ispec)

              tempx1l = ZERO
              tempx2l = ZERO
              b_tempx1l = ZERO
              b_tempx2l = ZERO
              do k = 1,NGLLX  ! merging loops NGLLX == NGLLZ
                ! derivative along x
                tempx1l = tempx1l + potential_acoustic(ibool(k,j,ispec)) * hprime_xx(i,k)
                b_tempx1l = b_tempx1l + b_potential_acoustic(ibool(k,j,ispec)) * hprime_xx(i,k)
                ! derivative along z
                tempx2l = tempx2l + potential_acoustic(ibool(i,k,ispec)) * hprime_zz(j,k)
                b_tempx2l = b_tempx2l + b_potential_acoustic(ibool(i,k,ispec)) * hprime_zz(j,k)
              enddo

              xixl = xix(i,j,ispec)
              xizl = xiz(i,j,ispec)
              gammaxl = gammax(i,j,ispec)
              gammazl = gammaz(i,j,ispec)

              ! derivatives of potential
              ! warning : the variable is named accel_loc but it is displ that is computed
              accel_loc(1) = (tempx1l*xixl + tempx2l*gammaxl) / rhol
              accel_loc(2) = (tempx1l*xizl + tempx2l*gammazl) / rhol

              b_accel_loc(1) = (b_tempx1l*xixl + b_tempx2l*gammaxl) / rhol
              b_accel_loc(2) = (b_tempx1l*xizl + b_tempx2l*gammazl) / rhol

              ! expressions using acceleration = 1/rho grad(potential_dot_dot_acoustic)
              rhorho_ac_Hessian_final1(i,j,ispec) = rhorho_ac_Hessian_final1(i,j,ispec) &
                                                  + dot_product(accel_loc(:),accel_loc(:)) &
                                                    * (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS)

              rhorho_ac_Hessian_final2(i,j,ispec) = rhorho_ac_Hessian_final2(i,j,ispec) &
                                                  + dot_product(accel_loc(:),b_accel_loc(:)) &
                                                    * (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS)
            enddo
          enddo
        endif
      enddo
    endif

  else
    ! on GPU
    ! computes contribution to density and bulk modulus kernel
    call compute_kernels_hess_cuda(Mesh_pointer,any_elastic,any_acoustic)
  endif

  end subroutine compute_kernels_Hessian

