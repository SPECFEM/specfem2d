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

! calculates C-PML acoustic/elastic memory variables and computes stress sigma

! second-order accurate convolution term calculation from equation (21) of
! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
! Anisotropic-medium PML for vector FETD with modified basis functions,
! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)


  subroutine pml_compute_memory_variables_acoustic(ispec,nglob,potential_acoustic_old,dux_dxl,dux_dzl)

! for acoustic elements

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,CPML_X_ONLY,CPML_Z_ONLY,ALPHA_LDDRK,BETA_LDDRK,C_LDDRK

  use specfem_par, only: time_stepping_scheme,i_stage,it,deltat, &
                         ibool,xix,xiz,gammax,gammaz, &
                         hprime_xx,hprime_zz, &
                         rmemory_acoustic_dux_dx,rmemory_acoustic_dux_dz, &
                         rmemory_acoustic_dux_dx_LDDRK,rmemory_acoustic_dux_dz_LDDRK, &
                         AXISYM,is_on_the_axis,hprimeBar_xx

  ! PML arrays
  use specfem_par, only: nspec_PML,ispec_is_PML,spec_to_PML,region_CPML, &
                K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  implicit none

  integer,intent(in) :: ispec

  integer,intent(in) :: nglob
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: potential_acoustic_old

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(inout) :: dux_dxl,dux_dzl

  ! local parameters
  integer :: i,j,k
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl

  ! PML
  integer :: ispec_PML
  integer :: CPML_region_local
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl,PML_dux_dzl,PML_dux_dxl_old,PML_dux_dzl_old

  double precision :: kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z,time_n,time_nsub1
  double precision :: A5,A6,A7
  double precision :: bb_zx_1,bb_zx_2
  double precision :: coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2
  double precision :: A8,A9,A10
  double precision :: bb_xz_1,bb_xz_2
  double precision :: coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2

  ! checks if anything to do in this slice
  if (nspec_PML == 0) return
  if (.not. ispec_is_PML(ispec)) return

  do j = 1,NGLLZ
    do i = 1,NGLLX
      ! stores initial derivatives
      PML_dux_dxl(i,j) = dux_dxl(i,j)
      PML_dux_dzl(i,j) = dux_dzl(i,j)

      select case (time_stepping_scheme)
      case (1)
        ! Newmark
        dux_dxi = 0._CUSTOM_REAL
        dux_dgamma = 0._CUSTOM_REAL

        do k = 1,NGLLX
          dux_dxi = dux_dxi + potential_acoustic_old(ibool(k,j,ispec)) * hprime_xx(i,k)
          dux_dgamma = dux_dgamma + potential_acoustic_old(ibool(i,k,ispec)) * hprime_zz(j,k)
        enddo

        ! AXISYM overwrites dux_dxi
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            dux_dxi = 0._CUSTOM_REAL
            do k = 1,NGLLX
              dux_dxi = dux_dxi + potential_acoustic_old(ibool(k,j,ispec)) * hprimeBar_xx(i,k)
            enddo
          endif
        endif

        ! derivatives of potential
        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        PML_dux_dxl_old(i,j) = dux_dxi * xixl + dux_dgamma * gammaxl
        PML_dux_dzl_old(i,j) = dux_dxi * xizl + dux_dgamma * gammazl

!        if (AXISYM) then ! TODO LDDRK
!          if (is_on_the_axis(ispec) .and. i == 1) then
!            ! dchi/dr=rho * u_r=0 on the axis
!            PML_dux_dxl_old(i,j) = 0._CUSTOM_REAL
!          endif
!        endif
      end select
    enddo
  enddo

  ! timing
  select case (time_stepping_scheme)
  case (1)
    ! Newmark
    time_n = (it-1) * deltat
    time_nsub1 = (it-2) * deltat
  case (2)
    ! LDDRK
    time_n = (it-1) * deltat + C_LDDRK(i_stage) * deltat
  case default
    call stop_the_code('Sorry, time stepping scheme not implemented yet in PML memory variable updates')
  end select

  ! local PML element
  ispec_PML = spec_to_PML(ispec)
  CPML_region_local = region_CPML(ispec)

  do j = 1,NGLLZ
    do i = 1,NGLLX
      kappa_x = K_x_store(i,j,ispec_PML)
      kappa_z = K_z_store(i,j,ispec_PML)

      d_x = d_x_store(i,j,ispec_PML)
      d_z = d_z_store(i,j,ispec_PML)

      alpha_x = alpha_x_store(i,j,ispec_PML)
      alpha_z = alpha_z_store(i,j,ispec_PML)

      beta_x = alpha_x + d_x / kappa_x
      beta_z = alpha_z + d_z / kappa_z

      ! gets PML coefficients
      ! the subroutine of lik_parameter_computation is located at the end of compute_forces_viscoelastic.F90
      call lik_parameter_computation(deltat,kappa_z,beta_z,alpha_z,kappa_x,beta_x,alpha_x, &
                                     CPML_region_local,31,A5,A6,A7,bb_zx_1,bb_zx_2, &
                                     coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2)

      call lik_parameter_computation(deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                     CPML_region_local,13,A8,A9,A10,bb_xz_1,bb_xz_2, &
                                     coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2)

      select case (time_stepping_scheme)
      case (1)
        ! Newmark
        rmemory_acoustic_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
                                                   coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
        rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                                                   coef1_zx_2 * PML_dux_dxl(i,j) + coef2_zx_2 * PML_dux_dxl_old(i,j)

        rmemory_acoustic_dux_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
                                                   coef1_xz_1 * PML_dux_dzl(i,j) + coef2_xz_1 * PML_dux_dzl_old(i,j)
        rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                                                   coef1_xz_2 * PML_dux_dzl(i,j) + coef2_xz_2 * PML_dux_dzl_old(i,j)
      case (2)
        ! LDDRK
        rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,1) = &
               ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,1) + &
               deltat * (-bb_zx_1 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + PML_dux_dxl(i,j))
        rmemory_acoustic_dux_dx(i,j,ispec_PML,1) = rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
               BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,1)

        rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) = &
               ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) + &
               deltat * (-bb_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + PML_dux_dxl(i,j))
        rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
               BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2)

        rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,1) = &
               ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,1) + &
               deltat * (-bb_xz_1 * rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + PML_dux_dzl(i,j))
        rmemory_acoustic_dux_dz(i,j,ispec_PML,1) = rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
               BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,1)

        rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) = &
               ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) + &
               deltat * (-bb_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + PML_dux_dzl(i,j))
        rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
               BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2)
      case default
        call stop_the_code('Time stepping scheme not implemented yet for PML memory variables')
      end select

      dux_dxl(i,j) = A5 * PML_dux_dxl(i,j) + A6 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
                                             A7 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2)
      dux_dzl(i,j) = A8 * PML_dux_dzl(i,j) + A9 *  rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
                                             A10 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2)
    enddo
  enddo

  end subroutine pml_compute_memory_variables_acoustic

!
!-------------------------------------------------------------------------------------
!

  subroutine pml_compute_memory_variables_elastic(ispec,nglob,displ_elastic_old, &
                                                  dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
                                                  dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime, &
                                                  PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl, &
                                                  PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old)

! for elastic elements

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLZ,NGLJ, &
    CPML_X_ONLY,CPML_Z_ONLY,ALPHA_LDDRK,BETA_LDDRK,C_LDDRK

  use specfem_par, only: time_stepping_scheme,i_stage,it,deltat, &
                         ibool,xix,xiz,gammax,gammaz, &
                         hprime_xx,hprime_zz, &
                         AXISYM,is_on_the_axis,hprimeBar_xx

  use specfem_par, only: rmemory_dux_dx,rmemory_dux_dz,rmemory_duz_dx,rmemory_duz_dz, &
                         rmemory_dux_dx_prime,rmemory_dux_dz_prime,rmemory_duz_dx_prime,rmemory_duz_dz_prime, &
                         rmemory_dux_dx_LDDRK,rmemory_dux_dz_LDDRK,rmemory_duz_dx_LDDRK,rmemory_duz_dz_LDDRK

  ! PML arrays
  use specfem_par, only: nspec_PML,ispec_is_PML,spec_to_PML,region_CPML, &
                K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store, &
                ROTATE_PML_ACTIVATE

  implicit none

  integer,intent(in) :: ispec

  integer,intent(in) :: nglob
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ_elastic_old

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(inout) :: dux_dxl,dux_dzl,duz_dxl,duz_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: dux_dxl_prime,dux_dzl_prime,duz_dxl_prime,duz_dzl_prime

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old

  ! local parameters
  integer :: i,j,k
  real(kind=CUSTOM_REAL) :: dux_dxi_old,dux_dgamma_old,duz_dxi_old,duz_dgamma_old
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl

  ! PML
  integer :: ispec_PML
  integer :: CPML_region_local
  double precision :: time_n,time_nsub1
  double precision :: kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z
  double precision :: A5,A6,A7
  double precision :: bb_zx_1,bb_zx_2,coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2
  double precision :: A8,A9,A10
  double precision :: bb_xz_1,bb_xz_2,coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2

  ! checks if anything to do in this slice
  if (nspec_PML == 0) return
  if (.not. ispec_is_PML(ispec)) return

  do j = 1,NGLLZ
    do i = 1,NGLLX
      ! stores initial derivatives
      PML_dux_dxl(i,j) = dux_dxl(i,j)
      PML_dux_dzl(i,j) = dux_dzl(i,j)
      PML_duz_dzl(i,j) = duz_dzl(i,j)
      PML_duz_dxl(i,j) = duz_dxl(i,j)

      if (AXISYM) then ! TODO LDDRK
        if (is_on_the_axis(ispec) .and. i == 1) then
          ! d_uz/dr=0 on the axis
          PML_duz_dxl(i,j) = 0._CUSTOM_REAL
        endif
      endif

      ! derivative along x and along z
      !
      ! first double loop over GLL points to compute and store gradients
      ! we can merge the two loops because NGLLX == NGLLZ
      dux_dxi_old = 0._CUSTOM_REAL
      duz_dxi_old = 0._CUSTOM_REAL
      dux_dgamma_old = 0._CUSTOM_REAL
      duz_dgamma_old = 0._CUSTOM_REAL
      do k = 1,NGLLX
        dux_dxi_old = dux_dxi_old + displ_elastic_old(1,ibool(k,j,ispec))*hprime_xx(i,k)
        duz_dxi_old = duz_dxi_old + displ_elastic_old(2,ibool(k,j,ispec))*hprime_xx(i,k)
        dux_dgamma_old = dux_dgamma_old + displ_elastic_old(1,ibool(i,k,ispec))*hprime_zz(j,k)
        duz_dgamma_old = duz_dgamma_old + displ_elastic_old(2,ibool(i,k,ispec))*hprime_zz(j,k)
      enddo

      ! AXISYM overwrite dux_dxi and duz_dxi
      if (AXISYM) then
        if (is_on_the_axis(ispec)) then
          dux_dxi_old = 0._CUSTOM_REAL
          duz_dxi_old = 0._CUSTOM_REAL
          do k = 1,NGLJ
            dux_dxi_old = dux_dxi_old + displ_elastic_old(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
            duz_dxi_old = duz_dxi_old + displ_elastic_old(2,ibool(k,j,ispec))*hprimeBar_xx(i,k)
          enddo
        endif
      endif

      ! derivatives of displacement
      xixl = xix(i,j,ispec)
      xizl = xiz(i,j,ispec)
      gammaxl = gammax(i,j,ispec)
      gammazl = gammaz(i,j,ispec)

      PML_dux_dxl_old(i,j) = dux_dxi_old*xixl + dux_dgamma_old*gammaxl !dux_dxl_old
      PML_dux_dzl_old(i,j) = dux_dxi_old*xizl + dux_dgamma_old*gammazl !dux_dzl_old
      PML_duz_dxl_old(i,j) = duz_dxi_old*xixl + duz_dgamma_old*gammaxl !duz_dxl_old
      PML_duz_dzl_old(i,j) = duz_dxi_old*xizl + duz_dgamma_old*gammazl !duz_dzl_old

      if (AXISYM) then ! TODO LDDRK
        if (is_on_the_axis(ispec) .and. i == 1) then
          ! d_uz/dr=0 on the axis
          PML_duz_dxl_old(i,j) = 0._CUSTOM_REAL
        endif
      endif
    enddo
  enddo

  !------------------------------------------------------------------------------
  !---------------------------- LEFT & RIGHT ------------------------------------
  !------------------------------------------------------------------------------
  ! timing
  select case (time_stepping_scheme)
  case (1)
    ! Newmark
    time_n = (it-1) * deltat
    time_nsub1 = (it-2) * deltat
  case (2)
    ! LDDRK
    time_n = (it-1) * deltat + C_LDDRK(i_stage) * deltat
  case default
    call stop_the_code('Sorry, time stepping scheme not implemented yet in PML memory variable updates')
  end select

  ! local PML element
  ispec_PML = spec_to_PML(ispec)
  CPML_region_local = region_CPML(ispec)

  do j = 1,NGLLZ
    do i = 1,NGLLX
      kappa_x = K_x_store(i,j,ispec_PML)
      kappa_z = K_z_store(i,j,ispec_PML)

      d_x = d_x_store(i,j,ispec_PML)
      d_z = d_z_store(i,j,ispec_PML)

      alpha_x = alpha_x_store(i,j,ispec_PML)
      alpha_z = alpha_z_store(i,j,ispec_PML)

      beta_x = alpha_x + d_x / kappa_x
      beta_z = alpha_z + d_z / kappa_z

      ! gets PML coefficients
      call lik_parameter_computation(deltat,kappa_z,beta_z,alpha_z,kappa_x,beta_x,alpha_x, &
                                     CPML_region_local,31,A5,A6,A7,bb_zx_1,bb_zx_2, &
                                     coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2)

      call lik_parameter_computation(deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                     CPML_region_local,13,A8,A9,A10,bb_xz_1,bb_xz_2, &
                                     coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2)

      select case (time_stepping_scheme)
      case (1)
        ! Newmark
        if (ROTATE_PML_ACTIVATE) then
          rmemory_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_dux_dx(i,j,ispec_PML,1) + &
                                            coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
          rmemory_dux_dz(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_dux_dz(i,j,ispec_PML,1) + &
                                            coef1_zx_1 * PML_dux_dzl(i,j) + coef2_zx_1 * PML_dux_dzl_old(i,j)
          rmemory_duz_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_duz_dx(i,j,ispec_PML,1) + &
                                            coef1_zx_1 * PML_duz_dxl(i,j) + coef2_zx_1 * PML_duz_dxl_old(i,j)
          rmemory_duz_dz(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_duz_dz(i,j,ispec_PML,1) + &
                                            coef1_zx_1 * PML_duz_dzl(i,j) + coef2_zx_1 * PML_duz_dzl_old(i,j)
          rmemory_dux_dx_prime(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_dux_dx_prime(i,j,ispec_PML,1) + &
                                                  coef1_xz_1 * PML_dux_dxl(i,j) + coef2_xz_1 * PML_dux_dxl_old(i,j)
          rmemory_dux_dz_prime(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_dux_dz_prime(i,j,ispec_PML,1) + &
                                                  coef1_xz_1 * PML_dux_dzl(i,j) + coef2_xz_1 * PML_dux_dzl_old(i,j)
          rmemory_duz_dx_prime(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_duz_dx_prime(i,j,ispec_PML,1) + &
                                                  coef1_xz_1 * PML_duz_dxl(i,j) + coef2_xz_1 * PML_duz_dxl_old(i,j)
          rmemory_duz_dz_prime(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_duz_dz_prime(i,j,ispec_PML,1) + &
                                                  coef1_xz_1 * PML_duz_dzl(i,j) + coef2_xz_1 * PML_duz_dzl_old(i,j)

          rmemory_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + &
                                            coef1_zx_2 * PML_dux_dxl(i,j) + coef2_zx_2 * PML_dux_dxl_old(i,j)
          rmemory_dux_dz(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dz(i,j,ispec_PML,2) + &
                                            coef1_zx_2 * PML_dux_dzl(i,j) + coef2_zx_2 * PML_dux_dzl_old(i,j)
          rmemory_duz_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + &
                                            coef1_zx_2 * PML_duz_dxl(i,j) + coef2_zx_2 * PML_duz_dxl_old(i,j)
          rmemory_duz_dz(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dz(i,j,ispec_PML,2) + &
                                            coef1_zx_2 * PML_duz_dzl(i,j) + coef2_zx_2 * PML_duz_dzl_old(i,j)

          rmemory_dux_dx_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dx_prime(i,j,ispec_PML,2) + &
                                                  coef1_xz_2 * PML_dux_dxl(i,j) + coef2_xz_2 * PML_dux_dxl_old(i,j)
          rmemory_dux_dz_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dz_prime(i,j,ispec_PML,2) + &
                                                  coef1_xz_2 * PML_dux_dzl(i,j) + coef2_xz_2 * PML_dux_dzl_old(i,j)
          rmemory_duz_dx_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dx_prime(i,j,ispec_PML,2) + &
                                                  coef1_xz_2 * PML_duz_dxl(i,j) + coef2_xz_2 * PML_duz_dxl_old(i,j)
          rmemory_duz_dz_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dz_prime(i,j,ispec_PML,2) + &
                                                  coef1_xz_2 * PML_duz_dzl(i,j) + coef2_xz_2 * PML_duz_dzl_old(i,j)

        else
          ! non-rotated, element aligns with x/y/z-coordinates
          rmemory_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_dux_dx(i,j,ispec_PML,1) + &
                                            coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
          rmemory_duz_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_duz_dx(i,j,ispec_PML,1) + &
                                            coef1_zx_1 * PML_duz_dxl(i,j) + coef2_zx_1 * PML_duz_dxl_old(i,j)
          rmemory_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + &
                                            coef1_zx_2 * PML_dux_dxl(i,j) + coef2_zx_2 * PML_dux_dxl_old(i,j)
          rmemory_duz_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + &
                                            coef1_zx_2 * PML_duz_dxl(i,j) + coef2_zx_2 * PML_duz_dxl_old(i,j)


          rmemory_dux_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_dux_dz(i,j,ispec_PML,1) + &
                                            coef1_xz_1 * PML_dux_dzl(i,j) + coef2_xz_1 * PML_dux_dzl_old(i,j)
          rmemory_duz_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_duz_dz(i,j,ispec_PML,1) + &
                                            coef1_xz_1 * PML_duz_dzl(i,j) + coef2_xz_1 * PML_duz_dzl_old(i,j)
          rmemory_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dz(i,j,ispec_PML,2) + &
                                            coef1_xz_2 * PML_dux_dzl(i,j) + coef2_xz_2 * PML_dux_dzl_old(i,j)
          rmemory_duz_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dz(i,j,ispec_PML,2) + &
                                            coef1_xz_2 * PML_duz_dzl(i,j) + coef2_xz_2 * PML_duz_dzl_old(i,j)
        endif ! ROTATE_PML_ACTIVATE

      case (2)
        ! LDDRK
        rmemory_dux_dx_LDDRK(i,j,ispec_PML,1) = ALPHA_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,1) + &
                                                deltat * (-bb_zx_1 * rmemory_dux_dx(i,j,ispec_PML,1) + PML_dux_dxl(i,j))
        rmemory_dux_dx(i,j,ispec_PML,1) = rmemory_dux_dx(i,j,ispec_PML,1) + &
                                          BETA_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,1)
        rmemory_duz_dx_LDDRK(i,j,ispec_PML,1) = ALPHA_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,1) + &
                                                deltat * (-bb_zx_1 * rmemory_duz_dx(i,j,ispec_PML,1) + PML_duz_dxl(i,j))
        rmemory_duz_dx(i,j,ispec_PML,1) = rmemory_duz_dx(i,j,ispec_PML,1) + &
                                          BETA_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,1)
        rmemory_dux_dx_LDDRK(i,j,ispec_PML,2) = ALPHA_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,2) + &
                                                deltat * (-bb_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + PML_dux_dxl(i,j))
        rmemory_dux_dx(i,j,ispec_PML,2) = rmemory_dux_dx(i,j,ispec_PML,2) + &
                                          BETA_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,2)
        rmemory_duz_dx_LDDRK(i,j,ispec_PML,2) = ALPHA_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,2) + &
                                                deltat * (-bb_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + PML_duz_dxl(i,j))
        rmemory_duz_dx(i,j,ispec_PML,2) = rmemory_duz_dx(i,j,ispec_PML,2) + &
                                          BETA_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,2)

        rmemory_dux_dz_LDDRK(i,j,ispec_PML,1) = ALPHA_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,1) + &
                                                deltat * (-bb_xz_1 * rmemory_dux_dz(i,j,ispec_PML,1) + PML_dux_dzl(i,j))
        rmemory_dux_dz(i,j,ispec_PML,1) = rmemory_dux_dz(i,j,ispec_PML,1) + &
                                          BETA_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,1)
        rmemory_duz_dz_LDDRK(i,j,ispec_PML,1) = ALPHA_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,1) + &
                                                deltat * (-bb_xz_1 * rmemory_duz_dz(i,j,ispec_PML,1) + PML_duz_dzl(i,j))
        rmemory_duz_dz(i,j,ispec_PML,1) = rmemory_duz_dz(i,j,ispec_PML,1) + &
                                          BETA_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,1)

        rmemory_dux_dz_LDDRK(i,j,ispec_PML,2) = ALPHA_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,2) + &
                                                deltat * (-bb_xz_2 * rmemory_dux_dz(i,j,ispec_PML,2) + PML_dux_dzl(i,j))
        rmemory_dux_dz(i,j,ispec_PML,2) = rmemory_dux_dz(i,j,ispec_PML,2) + &
                                          BETA_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,2)

        rmemory_duz_dz_LDDRK(i,j,ispec_PML,2) = ALPHA_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,2) + &
                                                deltat * (-bb_xz_2 * rmemory_duz_dz(i,j,ispec_PML,2) + PML_duz_dzl(i,j))
        rmemory_duz_dz(i,j,ispec_PML,2) = rmemory_duz_dz(i,j,ispec_PML,2) + &
                                          BETA_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,2)

      case default
        call stop_the_code('Time stepping scheme not implemented yet for elastic PML memory variable update')
      end select

      if (ROTATE_PML_ACTIVATE) then
        dux_dxl(i,j) = A5 * PML_dux_dxl(i,j)  + A6 * rmemory_dux_dx(i,j,ispec_PML,1) + A6 * rmemory_dux_dx(i,j,ispec_PML,2)
        dux_dzl(i,j) = A5 * PML_dux_dzl(i,j)  + A6 * rmemory_dux_dz(i,j,ispec_PML,1) + A6 * rmemory_dux_dz(i,j,ispec_PML,2)
        duz_dxl(i,j) = A5 * PML_duz_dxl(i,j)  + A6 * rmemory_duz_dx(i,j,ispec_PML,1) + A6 * rmemory_duz_dx(i,j,ispec_PML,2)
        duz_dzl(i,j) = A5 * PML_duz_dzl(i,j)  + A6 * rmemory_duz_dz(i,j,ispec_PML,1) + A6 * rmemory_duz_dz(i,j,ispec_PML,2)

        dux_dxl_prime(i,j) = A8 * PML_dux_dxl(i,j) + &
                             A9 * rmemory_dux_dx_prime(i,j,ispec_PML,1) + A10 * rmemory_dux_dx_prime(i,j,ispec_PML,2)
        dux_dzl_prime(i,j) = A8 * PML_dux_dzl(i,j) + &
                             A9 * rmemory_dux_dz_prime(i,j,ispec_PML,1) + A10 * rmemory_dux_dz_prime(i,j,ispec_PML,2)
        duz_dxl_prime(i,j) = A8 * PML_duz_dxl(i,j) + &
                             A9 * rmemory_duz_dx_prime(i,j,ispec_PML,1) + A10 * rmemory_duz_dx_prime(i,j,ispec_PML,2)
        duz_dzl_prime(i,j) = A8 * PML_duz_dzl(i,j) + &
                             A9 * rmemory_duz_dz_prime(i,j,ispec_PML,1) + A10 * rmemory_duz_dz_prime(i,j,ispec_PML,2)
      else
        dux_dxl(i,j) = A5 * PML_dux_dxl(i,j) + A6 * rmemory_dux_dx(i,j,ispec_PML,1) + A7 * rmemory_dux_dx(i,j,ispec_PML,2)
        duz_dxl(i,j) = A5 * PML_duz_dxl(i,j) + A6 * rmemory_duz_dx(i,j,ispec_PML,1) + A7 * rmemory_duz_dx(i,j,ispec_PML,2)
        dux_dzl(i,j) = A8 * PML_dux_dzl(i,j) + A9 * rmemory_dux_dz(i,j,ispec_PML,1) + A10 * rmemory_dux_dz(i,j,ispec_PML,2)
        duz_dzl(i,j) = A8 * PML_duz_dzl(i,j) + A9 * rmemory_duz_dz(i,j,ispec_PML,1) + A10 * rmemory_duz_dz(i,j,ispec_PML,2)
      endif
    enddo
  enddo

  ! AXISYM enforces zero derivatives on axis ! TODO LDDRK
  if (AXISYM) then
    if (is_on_the_axis(ispec)) then
      ! d_uz/dr=0 on the axis
      ! i == 1
      do j = 1,NGLLZ
        rmemory_duz_dx(1,j,ispec_PML,1) = 0._CUSTOM_REAL
        rmemory_duz_dx(1,j,ispec_PML,2) = 0._CUSTOM_REAL
        duz_dxl(1,j) = 0._CUSTOM_REAL
        if (time_stepping_scheme == 2) then
          rmemory_duz_dx_LDDRK(1,j,ispec_PML,1) = 0._CUSTOM_REAL
          rmemory_duz_dx_LDDRK(1,j,ispec_PML,2) = 0._CUSTOM_REAL
        endif
        if (ROTATE_PML_ACTIVATE) then
          duz_dxl_prime(1,j) = 0._CUSTOM_REAL
          rmemory_duz_dx_prime(1,j,ispec_PML,1) = 0._CUSTOM_REAL
          rmemory_duz_dx_prime(1,j,ispec_PML,2) = 0._CUSTOM_REAL
        endif
      enddo
    endif
  endif

  end subroutine pml_compute_memory_variables_elastic

!
!-------------------------------------------------------------------------------------
!

  subroutine pml_compute_memory_variables_viscoelastic(ispec,nglob,displ_elastic_old, &
                                                       dux_dxl,dux_dzl,duz_dxl,duz_dzl, &
                                                       kappa_pml_dux_dxl,kappa_pml_duz_dzl, &
                                                       mu_pml_dux_dxl,mu_pml_duz_dzl,mu_pml_dux_dzl, &
                                                       mu_pml_duz_dxl,kappa_dux_dxl,kappa_duz_dzl,mu_dux_dxl, &
                                                       mu_duz_dzl,mu_dux_dzl,mu_duz_dxl)

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLZ,CPML_X_ONLY,CPML_Z_ONLY

  use specfem_par, only: time_stepping_scheme,deltat, &
                         ibool,xix,xiz,gammax,gammaz,hprime_xx,hprime_zz

  use specfem_par, only: rmemory_dux_dx,rmemory_dux_dz,rmemory_duz_dx,rmemory_duz_dz, &
                         kaPML_rmemory_dux_dxl,kaPML_rmemory_duz_dzl,muPML_rmemory_dux_dxl, &
                         muPML_rmemory_duz_dzl,muPML_rmemory_duz_dxl,muPML_rmemory_dux_dzl, &
                         nspec_PML,ispec_is_PML

  use specfem_par, only: tau_epsilon_nu1,tau_epsilon_nu2,Mu_nu1,Mu_nu2, &
                         N_SLS,inv_tau_sigma_nu1,inv_tau_sigma_nu2,phi_nu1,phi_nu2,spec_to_PML, &
                         region_CPML,K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store, &
                         AXISYM,is_on_the_axis,hprimeBar_xx

  implicit none

  integer,intent(in) :: ispec
  integer,intent(in) :: nglob

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ_elastic_old

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(in) :: dux_dxl,dux_dzl,duz_dxl,duz_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: kappa_pml_dux_dxl,kappa_pml_duz_dzl,mu_pml_dux_dxl, &
                                                                mu_pml_duz_dzl,mu_pml_dux_dzl,mu_pml_duz_dxl

  ! local parameters
  integer :: i,j,k,i_sls
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl
  real(kind=CUSTOM_REAL) :: dux_dxi_old,dux_dgamma_old,duz_dxi_old,duz_dgamma_old
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: PML_dux_dxl_old,PML_dux_dzl_old,PML_duz_dxl_old,PML_duz_dzl_old

  ! PML
  integer :: ispec_PML
  integer :: CPML_region_local
  double precision :: kappa_x,kappa_z,d_x,d_z,alpha_x,alpha_z,beta_x,beta_z
  double precision :: tauinv_kappa,tauinv_mu,tao_epsilon_kappa,tao_epsilon_mu,Mu_kappa,Mu_mu
  double precision :: A_01,A_2_ka,A_3_ka,A_ik_ka,A_2_mu,A_3_mu,A_ik_mu, &
                      A_10,A_4_ka,A_5_ka,A_ki_ka,A_4_mu,A_5_mu,A_ki_mu, &
                      A_2_ka_sum,A_2_mu_sum,A_3_ka_sum,A_3_mu_sum,A_4_ka_sum,A_4_mu_sum,A_5_ka_sum,A_5_mu_sum
  double precision :: coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2, &
                      coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2, &
                      coef0_l_ka,coef1_l_ka,coef2_l_ka,coef0_l_mu,coef1_l_mu,coef2_l_mu

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: kappa_dux_dxl,kappa_duz_dzl,mu_dux_dxl,mu_duz_dzl, &
                                                                mu_dux_dzl,mu_duz_dxl
  real(kind=CUSTOM_REAL) :: sum_kappa_dux_dx,sum_kappa_duz_dz,sum_mu_dux_dx,sum_mu_duz_dz, &
                            sum_mu_dux_dz,sum_mu_duz_dx, kappa_sum_dux_dxl,kappa_sum_duz_dzl, &
                            mu_sum_dux_dxl,mu_sum_duz_dzl,mu_sum_dux_dzl,mu_sum_duz_dxl

  double precision :: A_kappa,A_mu

  ! checks if anything to do in this slice
  if (nspec_PML == 0) return
  if (.not. ispec_is_PML(ispec)) return

  ! compute the spatial derivatives

  ! AXISYM case
  if (AXISYM .and. is_on_the_axis(ispec)) then

  do j = 1,NGLLZ
    do i = 1,NGLLX

      ! stores initial derivatives
      PML_dux_dxl(i,j) = dux_dxl(i,j)
      PML_dux_dzl(i,j) = dux_dzl(i,j)
      PML_duz_dzl(i,j) = duz_dzl(i,j)
      PML_duz_dxl(i,j) = duz_dxl(i,j)

      dux_dxi_old = 0._CUSTOM_REAL
      duz_dxi_old = 0._CUSTOM_REAL
      dux_dgamma_old = 0._CUSTOM_REAL
      duz_dgamma_old = 0._CUSTOM_REAL
      do k = 1,NGLLX
        dux_dxi_old = dux_dxi_old + displ_elastic_old(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
        duz_dxi_old = duz_dxi_old + displ_elastic_old(2,ibool(k,j,ispec))*hprimeBar_xx(i,k)
        dux_dgamma_old = dux_dgamma_old + displ_elastic_old(1,ibool(i,k,ispec))*hprime_zz(j,k)
        duz_dgamma_old = duz_dgamma_old + displ_elastic_old(2,ibool(i,k,ispec))*hprime_zz(j,k)
      enddo

      ! derivatives of displacement
      xixl = xix(i,j,ispec)
      xizl = xiz(i,j,ispec)
      gammaxl = gammax(i,j,ispec)
      gammazl = gammaz(i,j,ispec)

      PML_dux_dxl_old(i,j) = dux_dxi_old*xixl + dux_dgamma_old*gammaxl ! this is dux_dxl_old
      PML_dux_dzl_old(i,j) = dux_dxi_old*xizl + dux_dgamma_old*gammazl ! this is dux_dzl_old
      PML_duz_dxl_old(i,j) = duz_dxi_old*xixl + duz_dgamma_old*gammaxl ! this is duz_dxl_old
      PML_duz_dzl_old(i,j) = duz_dxi_old*xizl + duz_dgamma_old*gammazl ! this is duz_dzl_old
    enddo
  enddo

  else ! non AXISYM case

  do j = 1,NGLLZ
    do i = 1,NGLLX

      ! stores initial derivatives
      PML_dux_dxl(i,j) = dux_dxl(i,j)
      PML_dux_dzl(i,j) = dux_dzl(i,j)
      PML_duz_dzl(i,j) = duz_dzl(i,j)
      PML_duz_dxl(i,j) = duz_dxl(i,j)

      dux_dxi_old = 0._CUSTOM_REAL
      duz_dxi_old = 0._CUSTOM_REAL
      dux_dgamma_old = 0._CUSTOM_REAL
      duz_dgamma_old = 0._CUSTOM_REAL
      do k = 1,NGLLX
        dux_dxi_old = dux_dxi_old + displ_elastic_old(1,ibool(k,j,ispec))*hprime_xx(i,k)
        duz_dxi_old = duz_dxi_old + displ_elastic_old(2,ibool(k,j,ispec))*hprime_xx(i,k)
        dux_dgamma_old = dux_dgamma_old + displ_elastic_old(1,ibool(i,k,ispec))*hprime_zz(j,k)
        duz_dgamma_old = duz_dgamma_old + displ_elastic_old(2,ibool(i,k,ispec))*hprime_zz(j,k)
      enddo

      ! derivatives of displacement
      xixl = xix(i,j,ispec)
      xizl = xiz(i,j,ispec)
      gammaxl = gammax(i,j,ispec)
      gammazl = gammaz(i,j,ispec)

      PML_dux_dxl_old(i,j) = dux_dxi_old*xixl + dux_dgamma_old*gammaxl ! this is dux_dxl_old
      PML_dux_dzl_old(i,j) = dux_dxi_old*xizl + dux_dgamma_old*gammazl ! this is dux_dzl_old
      PML_duz_dxl_old(i,j) = duz_dxi_old*xixl + duz_dgamma_old*gammaxl ! this is duz_dxl_old
      PML_duz_dzl_old(i,j) = duz_dxi_old*xizl + duz_dgamma_old*gammazl ! this is duz_dzl_old
    enddo
  enddo

  endif ! end of test on AXISYM

  ! local PML element
  ispec_PML = spec_to_PML(ispec)
  CPML_region_local = region_CPML(ispec)

  do j = 1,NGLLZ
    do i = 1,NGLLX
      kappa_x = K_x_store(i,j,ispec_PML)
      kappa_z = K_z_store(i,j,ispec_PML)

      d_x = d_x_store(i,j,ispec_PML)
      d_z = d_z_store(i,j,ispec_PML)

      alpha_x = alpha_x_store(i,j,ispec_PML)
      alpha_z = alpha_z_store(i,j,ispec_PML)

      beta_x = alpha_x + d_x / kappa_x
      beta_z = alpha_z + d_z / kappa_z

      kappa_sum_dux_dxl = 0._CUSTOM_REAL
      kappa_sum_duz_dzl = 0._CUSTOM_REAL
      mu_sum_dux_dxl = 0._CUSTOM_REAL
      mu_sum_duz_dzl = 0._CUSTOM_REAL
      mu_sum_dux_dzl = 0._CUSTOM_REAL
      mu_sum_duz_dxl = 0._CUSTOM_REAL

      A_2_ka_sum = 0._CUSTOM_REAL
      A_2_mu_sum = 0._CUSTOM_REAL
      A_3_ka_sum = 0._CUSTOM_REAL
      A_3_mu_sum = 0._CUSTOM_REAL
      A_4_ka_sum = 0._CUSTOM_REAL
      A_4_mu_sum = 0._CUSTOM_REAL
      A_5_ka_sum = 0._CUSTOM_REAL
      A_5_mu_sum = 0._CUSTOM_REAL

      if (inv_tau_sigma_nu1(i,j,ispec,1) < 0.) cycle

      sum_kappa_dux_dx = 0._CUSTOM_REAL
      sum_kappa_duz_dz = 0._CUSTOM_REAL
      sum_mu_dux_dx = 0._CUSTOM_REAL
      sum_mu_duz_dz = 0._CUSTOM_REAL
      sum_mu_dux_dz = 0._CUSTOM_REAL
      sum_mu_duz_dx = 0._CUSTOM_REAL

      ! loop on all the standard linear solids

      call compute_coef_convolution(alpha_z,deltat,coef0_zx_1,coef1_zx_1,coef2_zx_1)
      call compute_coef_convolution(beta_x,deltat,coef0_zx_2,coef1_zx_2,coef2_zx_2)
      call compute_coef_convolution(alpha_x,deltat,coef0_xz_1,coef1_xz_1,coef2_xz_1)
      call compute_coef_convolution(beta_z,deltat,coef0_xz_2,coef1_xz_2,coef2_xz_2)


      select case (time_stepping_scheme)
      case (1)
        ! Newmark
! alpha_z convolve dux_dx
        rmemory_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_dux_dx(i,j,ispec_PML,1) + &
                                          coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)

! alpha_z convolve duz_dx
        rmemory_duz_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_duz_dx(i,j,ispec_PML,1) + &
                                          coef1_zx_1 * PML_duz_dxl(i,j) + coef2_zx_1 * PML_duz_dxl_old(i,j)

! beta_x convolve dux_dx
        rmemory_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + &
                                          coef1_zx_2 * PML_dux_dxl(i,j) + coef2_zx_2 * PML_dux_dxl_old(i,j)

! beta_x convolve duz_dx
        rmemory_duz_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + &
                                          coef1_zx_2 * PML_duz_dxl(i,j) + coef2_zx_2 * PML_duz_dxl_old(i,j)

! alpha_x convolve dux_dz
        rmemory_dux_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_dux_dz(i,j,ispec_PML,1) + &
                                          coef1_xz_1 * PML_dux_dzl(i,j) + coef2_xz_1 * PML_dux_dzl_old(i,j)

! alpha_x convolve duz_dz
        rmemory_duz_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_duz_dz(i,j,ispec_PML,1) + &
                                          coef1_xz_1 * PML_duz_dzl(i,j) + coef2_xz_1 * PML_duz_dzl_old(i,j)

! beta_z convolve dux_dz
        rmemory_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dz(i,j,ispec_PML,2) + &
                                          coef1_xz_2 * PML_dux_dzl(i,j) + coef2_xz_2 * PML_dux_dzl_old(i,j)

! beta_z convolve duz_dz
        rmemory_duz_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dz(i,j,ispec_PML,2) + &
                                          coef1_xz_2 * PML_duz_dzl(i,j) + coef2_xz_2 * PML_duz_dzl_old(i,j)
      case (2)
      ! LDDRK
        call stop_the_code('Time stepping scheme LDDRK not implemented yet for viscoelastic PML memory variable update')

      case default
        call stop_the_code('Unknown time scheme for viscoelastic PML memory variable update')

      end select

      do i_sls = 1,N_SLS
        tauinv_kappa = inv_tau_sigma_nu1(i,j,ispec,i_sls)
        tauinv_mu = inv_tau_sigma_nu2(i,j,ispec,i_sls)
        tao_epsilon_kappa = tau_epsilon_nu1(i,j,ispec,i_sls)
        tao_epsilon_mu = tau_epsilon_nu2(i,j,ispec,i_sls)
        Mu_kappa = Mu_nu1(i,j,ispec) * dble(N_SLS)
        Mu_mu = Mu_nu2(i,j,ispec) * dble(N_SLS)

        A_kappa = phi_nu1(i,j,ispec,i_sls)
        A_mu = phi_nu2(i,j,ispec,i_sls)

        ! gets viscelastic PML coefficients
        call lik_parameter_computation_viscoelastic(kappa_z,beta_z,alpha_z,kappa_x,beta_x,alpha_x, &
                                              CPML_region_local,31,A_01,A_2_ka,A_3_ka,A_ki_ka,tauinv_kappa, &
                                              tao_epsilon_kappa)
        call lik_parameter_computation_viscoelastic(kappa_z,beta_z,alpha_z,kappa_x,beta_x,alpha_x, &
                                              CPML_region_local,31,A_01,A_2_mu,A_3_mu,A_ki_mu,tauinv_mu, &
                                              tao_epsilon_mu)

        call lik_parameter_computation_viscoelastic(kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                              CPML_region_local,13,A_10,A_4_ka,A_5_ka,A_ik_ka,tauinv_kappa, &
                                              tao_epsilon_kappa)
        call lik_parameter_computation_viscoelastic(kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                              CPML_region_local,13,A_10,A_4_mu,A_5_mu,A_ik_mu,tauinv_mu, &
                                              tao_epsilon_mu)

        call compute_coef_convolution(tauinv_kappa,deltat,coef0_l_ka,coef1_l_ka,coef2_l_ka)
        call compute_coef_convolution(tauinv_mu,deltat,coef0_l_mu,coef1_l_mu,coef2_l_mu)

        select case (time_stepping_scheme)

        case (1)
        ! Newmark

! inv_tau_sigma_nu1 convolve dux_dx
          kaPML_rmemory_dux_dxl(i,j,ispec_PML,i_sls) = coef0_l_ka * kaPML_rmemory_dux_dxl(i,j,ispec_PML,i_sls) + &
                                                       coef1_l_ka * PML_dux_dxl(i,j) + coef2_l_ka * PML_dux_dxl_old(i,j)

! inv_tau_sigma_nu1 convolve duz_dz
          kaPML_rmemory_duz_dzl(i,j,ispec_PML,i_sls) = coef0_l_ka * kaPML_rmemory_duz_dzl(i,j,ispec_PML,i_sls) + &
                                                       coef1_l_ka * PML_duz_dzl(i,j) + coef2_l_ka * PML_duz_dzl_old(i,j)

! inv_tau_sigma_nu2 convolve dux_dx
          muPML_rmemory_dux_dxl(i,j,ispec_PML,i_sls) = coef0_l_mu * muPML_rmemory_dux_dxl(i,j,ispec_PML,i_sls) + &
                                                       coef1_l_mu * PML_dux_dxl(i,j) + coef2_l_mu * PML_dux_dxl_old(i,j)

! inv_tau_sigma_nu2 convolve duz_dz
          muPML_rmemory_duz_dzl(i,j,ispec_PML,i_sls) = coef0_l_mu * muPML_rmemory_duz_dzl(i,j,ispec_PML,i_sls) + &
                                                       coef1_l_mu * PML_duz_dzl(i,j) + coef2_l_mu * PML_duz_dzl_old(i,j)

! inv_tau_sigma_nu2 convolve dux_dz
          muPML_rmemory_dux_dzl(i,j,ispec_PML,i_sls) = coef0_l_mu * muPML_rmemory_dux_dzl(i,j,ispec_PML,i_sls) + &
                                                       coef1_l_mu * PML_dux_dzl(i,j) + coef2_l_mu * PML_dux_dzl_old(i,j)

! inv_tau_sigma_nu2 convolve duz_dx
          muPML_rmemory_duz_dxl(i,j,ispec_PML,i_sls) = coef0_l_mu * muPML_rmemory_duz_dxl(i,j,ispec_PML,i_sls) + &
                                                       coef1_l_mu * PML_duz_dxl(i,j) + coef2_l_mu * PML_duz_dxl_old(i,j)
          A_2_ka_sum = A_2_ka_sum + A_2_ka
          A_2_mu_sum = A_2_mu_sum + A_2_mu
          A_3_ka_sum = A_3_ka_sum + A_3_ka
          A_3_mu_sum = A_3_mu_sum + A_3_mu
          A_4_ka_sum = A_4_ka_sum + A_4_ka
          A_4_mu_sum = A_4_mu_sum + A_4_mu
          A_5_ka_sum = A_5_ka_sum + A_5_ka
          A_5_mu_sum = A_5_mu_sum + A_5_mu

          kappa_sum_dux_dxl = kappa_sum_dux_dxl + A_ki_ka * kaPML_rmemory_dux_dxl(i,j,ispec_PML,i_sls) / Mu_kappa
          kappa_sum_duz_dzl = kappa_sum_duz_dzl + A_ik_ka * kaPML_rmemory_duz_dzl(i,j,ispec_PML,i_sls) / Mu_kappa

          mu_sum_dux_dxl = mu_sum_dux_dxl + A_ki_mu * muPML_rmemory_dux_dxl(i,j,ispec_PML,i_sls) / Mu_mu
          mu_sum_duz_dzl = mu_sum_duz_dzl + A_ik_mu * muPML_rmemory_duz_dzl(i,j,ispec_PML,i_sls) / Mu_mu
          mu_sum_dux_dzl = mu_sum_dux_dzl + A_ik_mu * muPML_rmemory_dux_dzl(i,j,ispec_PML,i_sls) / Mu_mu
          mu_sum_duz_dxl = mu_sum_duz_dxl + A_ki_mu * muPML_rmemory_duz_dxl(i,j,ispec_PML,i_sls) / Mu_mu

          sum_kappa_dux_dx = sum_kappa_dux_dx + A_kappa * kaPML_rmemory_dux_dxl(i,j,ispec_PML,i_sls)
          sum_kappa_duz_dz = sum_kappa_duz_dz + A_kappa * kaPML_rmemory_duz_dzl(i,j,ispec_PML,i_sls)

          sum_mu_dux_dx = sum_mu_dux_dx + A_mu * muPML_rmemory_dux_dxl(i,j,ispec_PML,i_sls)
          sum_mu_duz_dz = sum_mu_duz_dz + A_mu * muPML_rmemory_duz_dzl(i,j,ispec_PML,i_sls)
          sum_mu_dux_dz = sum_mu_dux_dz + A_mu * muPML_rmemory_dux_dzl(i,j,ispec_PML,i_sls)
          sum_mu_duz_dx = sum_mu_duz_dx + A_mu * muPML_rmemory_duz_dxl(i,j,ispec_PML,i_sls)

        case (2)
        ! LDDRK
          call stop_the_code('Time stepping scheme LDDRK not implemented yet for viscoelastic PML memory variable update')

        case default
          call stop_the_code('Unknown time scheme for viscoelastic PML memory variable update')

        end select

      enddo

      kappa_pml_dux_dxl(i,j) = A_01 * PML_dux_dxl(i,j) + kappa_sum_dux_dxl + &
                               A_2_ka_sum * rmemory_dux_dx(i,j,ispec_PML,1) / Mu_kappa + &
                               A_3_ka_sum * rmemory_dux_dx(i,j,ispec_PML,2) / Mu_kappa
      kappa_pml_duz_dzl(i,j) = A_10 * PML_duz_dzl(i,j) + kappa_sum_duz_dzl + &
                               A_4_ka_sum * rmemory_duz_dz(i,j,ispec_PML,1) / Mu_kappa + &
                               A_5_ka_sum * rmemory_duz_dz(i,j,ispec_PML,2) / Mu_kappa

      mu_pml_dux_dxl(i,j) = A_01 * PML_dux_dxl(i,j) + mu_sum_dux_dxl + &
                            A_2_mu_sum * rmemory_dux_dx(i,j,ispec_PML,1) / Mu_mu + &
                            A_3_mu_sum * rmemory_dux_dx(i,j,ispec_PML,2) / Mu_mu
      mu_pml_duz_dzl(i,j) = A_10 * PML_duz_dzl(i,j) + mu_sum_duz_dzl + &
                            A_4_mu_sum * rmemory_duz_dz(i,j,ispec_PML,1) / Mu_mu + &
                            A_5_mu_sum * rmemory_duz_dz(i,j,ispec_PML,2) / Mu_mu
      mu_pml_dux_dzl(i,j) = A_10 * PML_dux_dzl(i,j) + mu_sum_dux_dzl + &
                            A_4_mu_sum * rmemory_dux_dz(i,j,ispec_PML,1) / Mu_mu + &
                            A_5_mu_sum * rmemory_dux_dz(i,j,ispec_PML,2) / Mu_mu
      mu_pml_duz_dxl(i,j) = A_01 * PML_duz_dxl(i,j) + mu_sum_duz_dxl + &
                            A_2_mu_sum * rmemory_duz_dx(i,j,ispec_PML,1) / Mu_mu + &
                            A_3_mu_sum * rmemory_duz_dx(i,j,ispec_PML,2) / Mu_mu

      kappa_dux_dxl(i,j) = PML_dux_dxl(i,j) + sum_kappa_dux_dx
      kappa_duz_dzl(i,j) = PML_duz_dzl(i,j) + sum_kappa_duz_dz

      mu_dux_dxl(i,j) = PML_dux_dxl(i,j) + sum_mu_dux_dx
      mu_duz_dzl(i,j) = PML_duz_dzl(i,j) + sum_mu_duz_dz
      mu_dux_dzl(i,j) = PML_dux_dzl(i,j) + sum_mu_dux_dz
      mu_duz_dxl(i,j) = PML_duz_dxl(i,j) + sum_mu_duz_dx

    enddo
  enddo

  end subroutine pml_compute_memory_variables_viscoelastic

