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

  use constants, only: APPROXIMATE_HESS_KL

  use specfem_par, only: any_acoustic,any_elastic,any_poroelastic

  implicit none

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

  use specfem_par, only: ispec_is_elastic,rho_k, & ! AXISYM,
                         rho_kl,mu_kl,kappa_kl,rhop_kl,beta_kl,alpha_kl,bulk_c_kl,bulk_beta_kl, &
                         nglob,nspec,ibool,accel_elastic,b_displ_elastic, &
                         density,poroelastcoef,kmato,assign_external_model,rhoext,vsext,vpext, &
                         deltat,P_SV,displ_elastic, &
                         mu_k,kappa_k,ibool,hprime_xx,hprime_zz,xix,xiz,gammax,gammaz, &
                         GPU_MODE,it,NSTEP,NSTEP_BETWEEN_COMPUTE_KERNELS

  use specfem_par_gpu, only: Mesh_pointer,deltatf

  use specfem_par, only: c11_k,c13_k,c15_k,c33_k,c35_k,c55_k,ispec_is_anisotropic, &
                         rho_kl,c11_kl,c13_kl,c15_kl,c33_kl,c35_kl,c55_kl

  implicit none

  !local variables
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: b_dux_dxi,b_dux_dgamma,b_duz_dxi,b_duz_dgamma
  real(kind=CUSTOM_REAL) :: b_dux_dxl,b_duz_dxl,b_dux_dzl,b_duz_dzl
  real(kind=CUSTOM_REAL) :: dsxx,dsxz,dszz
  real(kind=CUSTOM_REAL) :: b_dsxx,b_dsxz,b_dszz
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal

  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl

  ! elastic kernels
  if (.not. GPU_MODE) then
    ! updates kernels on CPU
    do ispec = 1,nspec
      if (ispec_is_elastic(ispec)) then
        do j = 1,NGLLZ; do i = 1,NGLLX
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

          iglob = ibool(i,j,ispec)

          ! isotropic kernel contributions
          if (P_SV) then
            ! P-SV waves
            dsxx =  dux_dxl
            dsxz = HALF * (duz_dxl + dux_dzl)
            dszz =  duz_dzl

            b_dsxx =  b_dux_dxl
            b_dsxz = HALF * (b_duz_dxl + b_dux_dzl)
            b_dszz =  b_duz_dzl

            kappa_k(iglob) = (dsxx + dszz) *  (b_dsxx + b_dszz)
            mu_k(iglob) = dsxx * b_dsxx + dszz * b_dszz + &
                          2._CUSTOM_REAL * dsxz * b_dsxz - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k(iglob)
          else
            ! SH (membrane) waves
            mu_k(iglob) = dux_dxl * b_dux_dxl + dux_dzl * b_dux_dzl
          endif

          ! Voigt kernels, e.g., see Sieminski, 2007a,b
          if (ispec_is_anisotropic(ispec)) then
            c11_k(iglob) = dux_dxl*b_dux_dxl
            c13_k(iglob) = dux_dxl*b_duz_dzl + duz_dzl*b_dux_dxl
            c15_k(iglob) = 2*(dux_dxl*HALF*(b_dux_dzl+b_duz_dxl)+&
                           HALF*(dux_dzl+duz_dxl)*b_dux_dxl)
            c33_k(iglob) = duz_dzl*b_duz_dzl
            c35_k(iglob) = 2*(duz_dzl*HALF*(b_dux_dzl+b_duz_dxl)+&
                           HALF*(dux_dzl+duz_dxl)*b_duz_dzl)
            c55_k(iglob) = 4*HALF*(dux_dzl+duz_dxl)*HALF*(b_dux_dzl+b_duz_dxl)
          endif
        enddo; enddo
      endif
    enddo

    do iglob = 1,nglob
      rho_k(iglob) =  accel_elastic(1,iglob)*b_displ_elastic(1,iglob) + accel_elastic(2,iglob)*b_displ_elastic(2,iglob)
    enddo

  else
    ! updates kernels on GPU
    call compute_kernels_elastic_cuda(Mesh_pointer,deltatf)
  endif

 do ispec = 1, nspec
    if (ispec_is_elastic(ispec)) then
      ! isotropic kernels
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)

          ! for parameterization (rho,mu,kappa): "primary" kernels
          ! density kernel
          rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - rho_k(iglob)
          ! shear modulus kernel
          mu_kl(i,j,ispec) =  mu_kl(i,j,ispec) -  mu_k(iglob)
          ! bulk modulus kernel
          kappa_kl(i,j,ispec) = kappa_kl(i,j,ispec) -  kappa_k(iglob)

        enddo
      enddo
      ! Voigt kernels, e.g., see Sieminski, 2007a,b
      if (ispec_is_anisotropic(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - rho_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c11_kl(i,j,ispec) = c11_kl(i,j,ispec) - c11_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c13_kl(i,j,ispec) = c13_kl(i,j,ispec) - c13_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c15_kl(i,j,ispec) = c15_kl(i,j,ispec) - c15_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c33_kl(i,j,ispec) = c33_kl(i,j,ispec) - c33_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c35_kl(i,j,ispec) = c35_kl(i,j,ispec) - c35_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            c55_kl(i,j,ispec) = c55_kl(i,j,ispec) - c55_k(iglob) * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS)
            rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + c11_kl(i,j,ispec) + &
                                 c13_kl(i,j,ispec) + c15_kl(i,j,ispec) + c33_kl(i,j,ispec) + &
                                 c35_kl(i,j,ispec) + c55_kl(i,j,ispec)
          enddo
        enddo
      endif


    endif
  enddo

  ! only at the last time step we multiply by delta and parameter value, it is not necessary to do it at each iteration
  if (NSTEP - it == mod(NSTEP,NSTEP_BETWEEN_COMPUTE_KERNELS)) then

    do ispec = 1, nspec
      if (ispec_is_elastic(ispec)) then
        ! isotropic kernels
        do j = 1, NGLLZ
          do i = 1, NGLLX
            if (.not. assign_external_model) then
              rhol = density(1,kmato(ispec))
              mul = poroelastcoef(2,1,kmato(ispec))
              !if (AXISYM) then ! ABAB !!
              !Warning !! This is false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
                kappal = poroelastcoef(3,1,kmato(ispec)) - FOUR_THIRDS * mul
              !else
              !  kappal = poroelastcoef(3,1,kmato(ispec)) - mul
              !endif
            else
              rhol = rhoext(i,j,ispec)
              mul = rhoext(i,j,ispec)*vsext(i,j,ispec)*vsext(i,j,ispec)
              !if (AXISYM) then ! ABAB !!
              ! Warning !! This is false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
                kappal = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) - FOUR_THIRDS * mul
              !else
              !  kappal = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) - mul
              !endif
            endif

            ! for parameterization (rho,mu,kappa): "primary" kernels
            ! density kernel
            rho_kl(i,j,ispec) = rhol * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_kl(i,j,ispec)
            ! shear modulus kernel
            mu_kl(i,j,ispec) =  TWO * mul * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * mu_kl(i,j,ispec)
            ! bulk modulus kernel
            kappa_kl(i,j,ispec) = kappal * (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * kappa_kl(i,j,ispec)

            ! for parameterization (rho,beta,alpha):
            ! rho prime kernel
            rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + kappa_kl(i,j,ispec) + mu_kl(i,j,ispec)
            ! Vs kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
            beta_kl(i,j,ispec) = TWO * (mu_kl(i,j,ispec) - FOUR_THIRDS * mul/kappal * kappa_kl(i,j,ispec))
            ! Vp kernel
            ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) Check Kappa
            alpha_kl(i,j,ispec) = TWO * (1._CUSTOM_REAL + FOUR_THIRDS * mul/kappal) * kappa_kl(i,j,ispec)
            ! for bulk velocity c parameterization (rho,bulk_c,beta):
            bulk_c_kl(i,j,ispec) =  TWO * kappa_kl(i,j,ispec)
            bulk_beta_kl(i,j,ispec) =  TWO * mu_kl(i,j,ispec)
          enddo
        enddo

      endif ! elastic
    enddo !nspec loop

  endif ! it == NSTEP

  end subroutine compute_kernels_el

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_ac()

! acoustic kernel calculations
! see e.g. Tromp et al. (2005)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,ZERO,HALF,TWO

  use specfem_par, only: nspec,ispec_is_acoustic,ibool,kappal_ac_global,rhol_ac_global, &
                         poroelastcoef,density,kmato,assign_external_model,rhoext,vpext,deltat, &
                         hprime_xx,hprime_zz,xix,xiz,gammax,gammaz, &
                         potential_acoustic,b_potential_acoustic,potential_dot_dot_acoustic, &
                         accel_ac,b_displ_ac,NSTEP_BETWEEN_COMPUTE_KERNELS, &
                         rho_ac_kl,kappa_ac_kl,rhop_ac_kl,alpha_ac_kl,GPU_MODE

  use specfem_par_gpu, only: Mesh_pointer,deltatf

  implicit none

  !local variables
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,b_tempx1l,b_tempx2l
  double precision :: xixl,xizl,gammaxl,gammazl

  if (.not. GPU_MODE) then
    ! kernels on CPU
    do ispec = 1, nspec
      if (ispec_is_acoustic(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            if (.not. assign_external_model) then
              kappal_ac_global(iglob) = poroelastcoef(3,1,kmato(ispec))
              rhol_ac_global(iglob) = density(1,kmato(ispec))
            else
              kappal_ac_global(iglob) = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec)
              rhol_ac_global(iglob)   = rhoext(i,j,ispec)
            endif

            ! calcul the displacement by computing the gradient of potential / rho
            ! and calcul the acceleration by computing the gradient of potential_dot_dot / rho
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

            xixl = xix(i,j,ispec)
            xizl = xiz(i,j,ispec)
            gammaxl = gammax(i,j,ispec)
            gammazl = gammaz(i,j,ispec)

            ! derivatives of potential
            accel_ac(1,iglob) = (tempx1l*xixl + tempx2l*gammaxl) / rhol_ac_global(iglob)
            accel_ac(2,iglob) = (tempx1l*xizl + tempx2l*gammazl) / rhol_ac_global(iglob)
            b_displ_ac(1,iglob) = (b_tempx1l*xixl + b_tempx2l*gammaxl) / rhol_ac_global(iglob)
            b_displ_ac(2,iglob) = (b_tempx1l*xizl + b_tempx2l*gammazl) / rhol_ac_global(iglob)
          enddo !i = 1, NGLLX
        enddo !j = 1, NGLLZ
      endif
    enddo

    do ispec = 1,nspec
      if (ispec_is_acoustic(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            ! YANGL
            !!!! old expression (from elastic kernels)
            !!!rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) - rhol_ac_global(iglob)  * &
            !!!      dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * deltat
            !!!kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) - kappal_ac_global(iglob) * &
            !!!      potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * &
            !!!      b_potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob)&
            !!!      * deltat
            !!!! new expression (from PDE-constrained optimization, coupling terms changed as well)
            rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + rhol_ac_global(iglob) * &
                                   dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * &
                                   (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS)
                                   !warning : the variable is named accel_ac but it is displ_ac that is computed
            kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) + kappal_ac_global(iglob) * &
                                     potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * &
                                     b_potential_acoustic(iglob)/kappal_ac_global(iglob) * &
                                     (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS)
            ! YANGL
            rhop_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + kappa_ac_kl(i,j,ispec)
            alpha_ac_kl(i,j,ispec) = TWO *  kappa_ac_kl(i,j,ispec)
          enddo
        enddo
      endif
    enddo
  else
    ! on GPU
    call compute_kernels_acoustic_cuda(Mesh_pointer,deltatf)
  endif

  end subroutine compute_kernels_ac

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_po()

! kernel calculations
! see e.g. Morency et al. (2009)

  use constants, only: CUSTOM_REAL,FOUR_THIRDS,NGLLX,NGLLZ,TWO,HALF

  use specfem_par, only: nglob,nspec,ispec_is_poroelastic,ibool,deltat, &
                         kmato,permeability, &
                         accels_poroelastic,accelw_poroelastic,velocw_poroelastic, &
                         b_displs_poroelastic,b_displw_poroelastic, &
                         epsilondev_s,b_epsilondev_s, &
                         epsilondev_w,b_epsilondev_w, &
                         rhot_k,rhof_k,sm_k,eta_k,B_k,C_k, &
                         rhot_kl,rhof_kl,sm_kl,eta_kl,B_kl,C_kl,M_kl,M_k, &
                         mufr_kl,mufr_k,rhob_kl,rhofb_kl, &
                         mufrb_kl,phi_kl,rhobb_kl,rhofbb_kl,phib_kl,cpI_kl,cpII_kl,cs_kl,ratio_kl, &
                         GPU_MODE,NSTEP_BETWEEN_COMPUTE_KERNELS
  implicit none

  !local variables
  integer :: i,j,ispec,iglob
  real(kind=CUSTOM_REAL) :: rholb,dd1
  real(kind=CUSTOM_REAL) :: ratio
  real(kind=CUSTOM_REAL) :: dsxx,dsxz,dszz,dszx_xz
  real(kind=CUSTOM_REAL) :: b_dsxx,b_dsxz,b_dszz,b_dszx_xz
  real(kind=CUSTOM_REAL) :: dwxx,dwzz,b_dwxx,b_dwzz

  ! to evaluate cpI, cpII, and cs, and rI (poroelastic medium)
  double precision :: phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot
  double precision :: B_biot
  double precision :: perm_xx
  double precision :: afactor,bfactor,cfactor
  double precision :: gamma1,gamma2,gamma3,gamma4
  double precision :: cpIsquare,cpIIsquare,cssquare

  integer :: material

  ! safety check
  if (GPU_MODE) call stop_the_code('Error poroelastic kernels not implemented on GPUs yet')

  ! kernel contributions on global nodes
  do iglob = 1,nglob
    rhot_k(iglob) = accels_poroelastic(1,iglob) * b_displs_poroelastic(1,iglob) + &
                    accels_poroelastic(2,iglob) * b_displs_poroelastic(2,iglob)

    rhof_k(iglob) = accelw_poroelastic(1,iglob) * b_displs_poroelastic(1,iglob) + &
                    accelw_poroelastic(2,iglob) * b_displs_poroelastic(2,iglob) + &
                    accels_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                    accels_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)

    sm_k(iglob)  = accelw_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                   accelw_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)

    eta_k(iglob) = velocw_poroelastic(1,iglob) * b_displw_poroelastic(1,iglob) + &
                   velocw_poroelastic(2,iglob) * b_displw_poroelastic(2,iglob)
  enddo

  ! kernels on local nodes
  do ispec = 1, nspec
    if (ispec_is_poroelastic(ispec)) then

      ! gets poroelastic material
      call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)

      ! Biot coefficients for the input phi
      call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

      B_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr

      ! permeability
      material = kmato(ispec)
      perm_xx = permeability(1,material)

      ! Approximated velocities (no viscous dissipation)
      afactor = rho_bar - phi/tort*rho_f
      bfactor = H_biot + phi*rho_bar/(tort*rho_f)*M_biot - TWO*phi/tort*C_biot
      cfactor = phi/(tort*rho_f)*(H_biot*M_biot - C_biot*C_biot)

      cpIsquare = (bfactor + sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)
      cpIIsquare = (bfactor - sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)
      cssquare = mu_fr/afactor

      ! Approximated ratio r = amplitude "w" field/amplitude "s" field (no viscous dissipation)
      ! used later for wavespeed kernels calculation, which are presently implemented for inviscid case,
      ! contrary to primary and density-normalized kernels, which are consistent with viscous fluid case.
      gamma1 = H_biot - phi/tort*C_biot
      gamma2 = C_biot - phi/tort*M_biot
      gamma3 = phi/tort*( M_biot*(afactor/rho_f + phi/tort) - C_biot)
      gamma4 = phi/tort*( C_biot*(afactor/rho_f + phi/tort) - H_biot)

      ratio = HALF*(gamma1 - gamma3)/gamma4 + HALF*sqrt((gamma1-gamma3)**2/gamma4**2 + 4.d0 * gamma2/gamma4)

      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)

          rhot_kl(i,j,ispec) = rhot_kl(i,j,ispec) - (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_bar * rhot_k(iglob)
          rhof_kl(i,j,ispec) = rhof_kl(i,j,ispec) - (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_f * rhof_k(iglob)
          sm_kl(i,j,ispec) = sm_kl(i,j,ispec) - (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * rho_f*tort/phi * sm_k(iglob)

          !at the moment works with constant permeability
          eta_kl(i,j,ispec) = eta_kl(i,j,ispec) - (deltat * NSTEP_BETWEEN_COMPUTE_KERNELS) * eta_f/perm_xx * eta_k(iglob)

          ! for B_k & mufr_k
          dsxx = epsilondev_s(1,i,j,ispec) ! dux_dxl
          dszz = epsilondev_s(2,i,j,ispec) ! duz_dzl
          dsxz = epsilondev_s(3,i,j,ispec) ! dux_dzl
          dszx_xz = epsilondev_s(4,i,j,ispec) ! 0.5_CUSTOM_REAL * (duz_dxl + dux_dzl)

          b_dsxx = b_epsilondev_s(1,i,j,ispec) ! b_dux_dxl
          b_dszz = b_epsilondev_s(2,i,j,ispec) ! b_duz_dzl
          b_dsxz = b_epsilondev_s(3,i,j,ispec) ! b_dux_dzl
          b_dszx_xz = b_epsilondev_s(4,i,j,ispec) ! 0.5_CUSTOM_REAL * (b_duz_dxl + b_dux_dzl)

          B_k(iglob) = (dsxx + dszz) *  (b_dsxx + b_dszz) * (H_biot - FOUR_THIRDS * mu_fr)

          mufr_k(iglob) = (dsxx * b_dsxx + dszz * b_dszz + &
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

          C_k(iglob) =  ( (dsxx + dszz)*(b_dwxx + b_dwzz) + (dwxx + dwzz)*(b_dsxx + b_dszz) ) * C_biot
          M_k(iglob) = (dwxx + dwzz)*(b_dwxx + b_dwzz) * M_biot

          ! from older compute_forces_poro_fluid ...
          !C_k(iglob) =  ((dux_dxl + duz_dzl) *  (b_dwx_dxl + b_dwz_dzl) + &
          !                (dwx_dxl + dwz_dzl) *  (b_dux_dxl + b_duz_dzl)) * C_biot
          !M_k(iglob) = (dwx_dxl + dwz_dzl) *  (b_dwx_dxl + b_dwz_dzl) * M_biot

          B_kl(i,j,ispec) = B_kl(i,j,ispec) - (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS) * B_k(iglob)
          C_kl(i,j,ispec) = C_kl(i,j,ispec) - (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS) * C_k(iglob)
          M_kl(i,j,ispec) = M_kl(i,j,ispec) - (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS) * M_k(iglob)

          mufr_kl(i,j,ispec) = mufr_kl(i,j,ispec) - TWO * (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS) * mufr_k(iglob)

          ! density kernels
          rholb = rho_bar - phi*rho_f/tort
          rhob_kl(i,j,ispec) = rhot_kl(i,j,ispec) + B_kl(i,j,ispec) + mufr_kl(i,j,ispec)
          rhofb_kl(i,j,ispec) = rhof_kl(i,j,ispec) + C_kl(i,j,ispec) + M_kl(i,j,ispec) + sm_kl(i,j,ispec)

          mufrb_kl(i,j,ispec) = mufr_kl(i,j,ispec)
          phi_kl(i,j,ispec) = - sm_kl(i,j,ispec) - M_kl(i,j,ispec)

          ! wave speed kernels
          dd1 = (1._CUSTOM_REAL+rholb/rho_f)*ratio**2 + 2._CUSTOM_REAL*ratio + tort/phi

          rhobb_kl(i,j,ispec) = rhob_kl(i,j,ispec) &
                - phi*rho_f/(tort*B_biot) * &
                  (cpIIsquare + (cpIsquare - cpIIsquare)*( (phi / &
                  tort*ratio +1._CUSTOM_REAL)/dd1 + &
                  (rho_bar**2*ratio**2/rho_f**2*(phi / tort*ratio+1._CUSTOM_REAL)*(phi/tort*ratio + &
                  phi/tort * &
                  (1._CUSTOM_REAL+rho_f/rho_bar)-1._CUSTOM_REAL) )/dd1**2 ) - &
                  FOUR_THIRDS*cssquare ) &
                  * B_kl(i,j,ispec) &
                - rho_bar*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                  (phi/tort*ratio + &
                  1._CUSTOM_REAL)**2/dd1**2*M_kl(i,j,ispec) + &
                  rho_bar*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                  (phi/tort*ratio+1._CUSTOM_REAL)/dd1 - &
                  phi*ratio/tort*(phi / tort*ratio+1._CUSTOM_REAL)*&
                  (1._CUSTOM_REAL+rho_bar*ratio/rho_f)/dd1**2) &
                  * C_kl(i,j,ispec) &
                + phi*rho_f*cssquare / (tort*mu_fr) &
                  * mufrb_kl(i,j,ispec)

          rhofbb_kl(i,j,ispec) = rhofb_kl(i,j,ispec) &
                + phi*rho_f/(tort*B_biot) * (cpIIsquare + (cpIsquare - cpIIsquare)*( (phi/ &
                  tort*ratio +1._CUSTOM_REAL)/dd1+&
                  (rho_bar**2*ratio**2/rho_f**2*(phi/tort*ratio+1)*(phi/tort*ratio+ &
                  phi/tort*(1._CUSTOM_REAL+rho_f/rho_bar)-1._CUSTOM_REAL) )/dd1**2 )- &
                  FOUR_THIRDS*cssquare ) &
                  * B_kl(i,j,ispec) &
                + rho_bar*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                  (phi/tort*ratio + 1._CUSTOM_REAL)**2/dd1**2 &
                  * M_kl(i,j,ispec) &
                - rho_bar*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                  (phi/tort*ratio+1._CUSTOM_REAL)/dd1 - &
                  phi*ratio/tort*(phi/tort*ratio+1._CUSTOM_REAL)*&
                  (1._CUSTOM_REAL+rho_bar*ratio/rho_f)/dd1**2) &
                  * C_kl(i,j,ispec) &
                - phi*rho_f*cssquare/(tort*mu_fr) &
                  * mufrb_kl(i,j,ispec)

          phib_kl(i,j,ispec) = phi_kl(i,j,ispec) &
                - phi*rho_bar/(tort*B_biot) * ( cpIsquare - rho_f/rho_bar*cpIIsquare- &
                  (cpIsquare-cpIIsquare)*( (TWO*ratio**2*phi/tort + (1._CUSTOM_REAL+rho_f/rho_bar)* &
                  (TWO*ratio*phi/tort+1._CUSTOM_REAL))/dd1 + (phi/tort*ratio+1._CUSTOM_REAL)*(phi*&
                  ratio/tort+phi/tort*(1._CUSTOM_REAL+rho_f/rho_bar)-1._CUSTOM_REAL)*((1._CUSTOM_REAL+ &
                  rho_bar/rho_f-TWO*phi/tort)*ratio**2+TWO*ratio)/dd1**2 ) - &
                  FOUR_THIRDS*rho_f*cssquare/rho_bar ) &
                  * B_kl(i,j,ispec) &
                + rho_f/M_biot * (cpIsquare-cpIIsquare) &
                  *( TWO*ratio*(phi/tort*ratio+1._CUSTOM_REAL)/dd1 - &
                    (phi/tort*ratio+1._CUSTOM_REAL)**2 &
                    *((1._CUSTOM_REAL+rho_bar/rho_f-TWO*phi/tort)*ratio**2+TWO*ratio)/dd1**2) &
                  * M_kl(i,j,ispec) &
                + phi*rho_f/(tort*C_biot)* (cpIsquare-cpIIsquare)*ratio* (&
                  (1._CUSTOM_REAL+rho_f/rho_bar*ratio)/dd1 - (phi/tort*ratio+1._CUSTOM_REAL)* &
                  (1._CUSTOM_REAL+rho_bar/rho_f*ratio)*((1._CUSTOM_REAL+rho_bar/rho_f-TWO*phi/tort)*ratio+TWO)/dd1**2 ) &
                  * C_kl(i,j,ispec) &
                - phi*rho_f*cssquare /(tort*mu_fr) &
                  * mufrb_kl(i,j,ispec)

          ! wavespeed kernels
          cpI_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIsquare/B_biot*rho_bar*( &
                  1._CUSTOM_REAL-phi/tort + (phi/tort*ratio+ 1._CUSTOM_REAL)*(phi/tort*&
                  ratio+phi/tort* (1._CUSTOM_REAL+rho_f/rho_bar)-1._CUSTOM_REAL)/dd1 ) &
                  * B_kl(i,j,ispec) &
                + 2._CUSTOM_REAL*cpIsquare*rho_f*tort/(phi*M_biot) *&
                  (phi/tort*ratio+1._CUSTOM_REAL)**2/dd1 &
                  * M_kl(i,j,ispec) &
                + 2._CUSTOM_REAL*cpIsquare*rho_f/C_biot * &
                  (phi/tort*ratio+1._CUSTOM_REAL)* (1._CUSTOM_REAL+rho_bar/rho_f*ratio)/dd1 &
                  * C_kl(i,j,ispec)
          cpII_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIIsquare*rho_bar/B_biot * (&
                  phi*rho_f/(tort*rho_bar) - (phi/tort*ratio+ 1._CUSTOM_REAL)*(phi/tort*ratio+phi/tort* &
                  (1._CUSTOM_REAL+rho_f/rho_bar)-&
                  1._CUSTOM_REAL)/dd1  ) &
                  * B_kl(i,j,ispec) &
                + 2._CUSTOM_REAL*cpIIsquare*rho_f*tort/(phi*M_biot) * (&
                  1._CUSTOM_REAL - (phi/tort*ratio+ 1._CUSTOM_REAL)**2/dd1  ) &
                  * M_kl(i,j,ispec) &
                + 2._CUSTOM_REAL*cpIIsquare*rho_f/C_biot * (&
                  1._CUSTOM_REAL - (phi/tort*ratio+ 1._CUSTOM_REAL)*(1._CUSTOM_REAL+&
                  rho_bar/rho_f*ratio)/dd1  ) &
                  * C_kl(i,j,ispec)

          cs_kl(i,j,ispec) = - 8._CUSTOM_REAL/3._CUSTOM_REAL*cssquare* rho_bar/B_biot &
                  *(1._CUSTOM_REAL-phi*rho_f/(tort*rho_bar)) &
                  * B_kl(i,j,ispec) &
                + 2._CUSTOM_REAL*(rho_bar-rho_f*phi/tort)/mu_fr*cssquare &
                  * mufrb_kl(i,j,ispec)

          ratio_kl(i,j,ispec) = ratio*rho_bar*phi/(tort*B_biot) * (cpIsquare-cpIIsquare) &
                  * (phi/tort*(2._CUSTOM_REAL*ratio+1._CUSTOM_REAL+rho_f/rho_bar)/dd1 - (phi/tort*ratio+1._CUSTOM_REAL)*&
                    (phi/tort*ratio+phi/tort*( 1._CUSTOM_REAL+rho_f/rho_bar)-1._CUSTOM_REAL)*(2._CUSTOM_REAL*ratio*(&
                      1._CUSTOM_REAL+rho_bar/rho_f-phi/tort) + 2._CUSTOM_REAL)/dd1**2  ) &
                  * B_kl(i,j,ispec) &
                + ratio*rho_f*tort/(phi*M_biot)*(cpIsquare-cpIIsquare) * 2._CUSTOM_REAL*phi/tort &
                  * ( (phi/tort*ratio+1._CUSTOM_REAL)/dd1 - (phi/tort*ratio+1._CUSTOM_REAL)**2 &
                      * ((1._CUSTOM_REAL+rho_bar/rho_f-phi/tort)*ratio + 1._CUSTOM_REAL)/dd1**2 ) &
                  * M_kl(i,j,ispec) &
                + ratio*rho_f/C_biot*(cpIsquare-cpIIsquare) &
                  * ( (2._CUSTOM_REAL*phi*rho_bar*ratio/(tort*rho_f)+phi/tort+rho_bar/rho_f)/dd1 - &
                       2._CUSTOM_REAL*phi/tort*(phi/tort*ratio+1._CUSTOM_REAL)*(1._CUSTOM_REAL+rho_bar/rho_f*ratio) &
                      *((1._CUSTOM_REAL + rho_bar/rho_f - phi/tort)*ratio+1._CUSTOM_REAL)/dd1**2 ) &
                  * C_kl(i,j,ispec)
        enddo
      enddo
    endif
  enddo

  end subroutine compute_kernels_po

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_Hessian()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: nglob,nspec,ibool,ispec_is_acoustic,ispec_is_elastic, &
                         any_elastic,any_acoustic, &
                         accel_elastic,b_accel_elastic,accel_ac,b_accel_ac, &
                         rhorho_el_Hessian_final1,rhorho_el_Hessian_final2, &
                         rhorho_ac_Hessian_final1,rhorho_ac_Hessian_final2, &
                         deltat,GPU_MODE,NSTEP_BETWEEN_COMPUTE_KERNELS

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  !local variables
  real(kind=CUSTOM_REAL), dimension(nglob) :: rhorho_el_Hessian_temp1, rhorho_el_Hessian_temp2
  integer :: i,j,ispec,iglob


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
                                                    rhorho_el_Hessian_temp1(iglob) * (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS)
              rhorho_el_Hessian_final2(i,j,ispec) = rhorho_el_Hessian_final2(i,j,ispec) + &
                                                    rhorho_el_Hessian_temp2(iglob) * (deltat*NSTEP_BETWEEN_COMPUTE_KERNELS)
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
              rhorho_ac_Hessian_final1(i,j,ispec) = rhorho_ac_Hessian_final1(i,j,ispec) + &
                                                    dot_product(accel_ac(:,iglob),accel_ac(:,iglob)) * &
                                                    (deltat* NSTEP_BETWEEN_COMPUTE_KERNELS)
              rhorho_ac_Hessian_final2(i,j,ispec) = rhorho_ac_Hessian_final2(i,j,ispec) + &
                                                    dot_product(accel_ac(:,iglob),b_accel_ac(:,iglob)) * &
                                                    (deltat* NSTEP_BETWEEN_COMPUTE_KERNELS)
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

