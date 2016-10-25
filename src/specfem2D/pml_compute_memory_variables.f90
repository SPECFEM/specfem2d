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
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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
  integer :: singularity_type_zx,singularity_type_xz
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

        if (AXISYM) then
          if (is_on_the_axis(ispec) .and. i == 1) then
            ! dchi/dr=rho * u_r=0 on the axis
            PML_dux_dxl_old(i,j) = 0._CUSTOM_REAL
          endif
        endif
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
    stop 'Sorry, time stepping scheme not implemented yet in PML memory variable updates'
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
      call lik_parameter_computation(time_n,deltat,kappa_z,beta_z,alpha_z,kappa_x,beta_x,alpha_x, &
                                     CPML_region_local,31,A5,A6,A7,singularity_type_zx,bb_zx_1,bb_zx_2, &
                                     coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2)

      call lik_parameter_computation(time_n,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                     CPML_region_local,13,A8,A9,A10,singularity_type_xz,bb_xz_1,bb_xz_2, &
                                     coef0_xz_1,coef1_xz_1,coef2_xz_1,coef0_xz_2,coef1_xz_2,coef2_xz_2)

      select case (time_stepping_scheme)
      case (1)
        ! Newmark
        rmemory_acoustic_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
                                                   coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
        if (singularity_type_zx == 0) then
          rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                                                     coef1_zx_2 * PML_dux_dxl(i,j) + coef2_zx_2 * PML_dux_dxl_old(i,j)
        else
          rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                                                     coef1_zx_2 * time_n * PML_dux_dxl(i,j) + &
                                                     coef2_zx_2 * time_nsub1 * PML_dux_dxl_old(i,j)
        endif

        rmemory_acoustic_dux_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
                                                   coef1_xz_1 * PML_dux_dzl(i,j) + coef2_xz_1 * PML_dux_dzl_old(i,j)
        if (singularity_type_xz == 0) then
          rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                                                     coef1_xz_2 * PML_dux_dzl(i,j) + coef2_xz_2 * PML_dux_dzl_old(i,j)
        else
          rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                                                     coef1_xz_2 * time_n * PML_dux_dzl(i,j) + &
                                                     coef2_xz_2 * time_nsub1 * PML_dux_dzl_old(i,j)
        endif

      case (2)
        ! LDDRK
        rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,1) = &
               ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,1) + &
               deltat * (-bb_zx_1 * rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + PML_dux_dxl(i,j))
        rmemory_acoustic_dux_dx(i,j,ispec_PML,1) = rmemory_acoustic_dux_dx(i,j,ispec_PML,1) + &
               BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,1)

        if (singularity_type_zx == 0) then
          rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) = &
                 ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) + &
                 deltat * (-bb_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + PML_dux_dxl(i,j))
          rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                 BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2)
        else
          rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) = &
                 ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2) + &
                 deltat * (-bb_zx_2 * rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + PML_dux_dxl(i,j) * time_n)
          rmemory_acoustic_dux_dx(i,j,ispec_PML,2) = rmemory_acoustic_dux_dx(i,j,ispec_PML,2) + &
                 BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dx_LDDRK(i,j,ispec_PML,2)
        endif

        rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,1) = &
               ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,1) + &
               deltat * (-bb_xz_1 * rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + PML_dux_dzl(i,j))
        rmemory_acoustic_dux_dz(i,j,ispec_PML,1) = rmemory_acoustic_dux_dz(i,j,ispec_PML,1) + &
               BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,1)

        if (singularity_type_xz == 0) then
          rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) = &
                 ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) + &
                 deltat * (-bb_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + PML_dux_dzl(i,j))
          rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                 BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2)
        else
          rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) = &
                 ALPHA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2) + &
                 deltat * (-bb_xz_2 * rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + PML_dux_dzl(i,j) * time_n)
          rmemory_acoustic_dux_dz(i,j,ispec_PML,2) = rmemory_acoustic_dux_dz(i,j,ispec_PML,2) + &
                 BETA_LDDRK(i_stage) * rmemory_acoustic_dux_dz_LDDRK(i,j,ispec_PML,2)
        endif

      case default
        stop 'Time stepping scheme not implemented yet for PML memory variables'
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
  integer :: singularity_type_zx,singularity_type_xz
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

      if (AXISYM) then
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
    stop 'Sorry, time stepping scheme not implemented yet in PML memory variable updates'
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
      call lik_parameter_computation(time_n,deltat,kappa_z,beta_z,alpha_z,kappa_x,beta_x,alpha_x, &
                                     CPML_region_local,31,A5,A6,A7,singularity_type_zx,bb_zx_1,bb_zx_2, &
                                     coef0_zx_1,coef1_zx_1,coef2_zx_1,coef0_zx_2,coef1_zx_2,coef2_zx_2)

      call lik_parameter_computation(time_n,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                     CPML_region_local,13,A8,A9,A10,singularity_type_xz,bb_xz_1,bb_xz_2, &
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
          if (singularity_type_zx == 0) then
            rmemory_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + &
                                              coef1_zx_2 * PML_dux_dxl(i,j) + coef2_zx_2 * PML_dux_dxl_old(i,j)
            rmemory_dux_dz(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dz(i,j,ispec_PML,2) + &
                                              coef1_zx_2 * PML_dux_dzl(i,j) + coef2_zx_2 * PML_dux_dzl_old(i,j)
            rmemory_duz_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + &
                                              coef1_zx_2 * PML_duz_dxl(i,j) + coef2_zx_2 * PML_duz_dxl_old(i,j)
            rmemory_duz_dz(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dz(i,j,ispec_PML,2) + &
                                              coef1_zx_2 * PML_duz_dzl(i,j) + coef2_zx_2 * PML_duz_dzl_old(i,j)
          else
            rmemory_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + &
                                              coef1_zx_2 * time_n * PML_dux_dxl(i,j) + &
                                              coef2_zx_2 * time_nsub1 * PML_dux_dxl_old(i,j)
            rmemory_dux_dz(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dz(i,j,ispec_PML,2) + &
                                              coef1_zx_2 * time_n * PML_dux_dzl(i,j) + &
                                              coef2_zx_2 * time_nsub1 * PML_dux_dzl_old(i,j)
            rmemory_duz_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + &
                                              coef1_zx_2 * time_n * PML_duz_dxl(i,j) + &
                                              coef2_zx_2 * time_nsub1 * PML_duz_dxl_old(i,j)
            rmemory_duz_dz(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dz(i,j,ispec_PML,2) + &
                                              coef1_zx_2 * time_n * PML_duz_dzl(i,j) + &
                                              coef2_zx_2 * time_nsub1 * PML_duz_dzl_old(i,j)
          endif

          if (singularity_type_xz == 0) then
            rmemory_dux_dx_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dx_prime(i,j,ispec_PML,2) + &
                                                    coef1_xz_2 * PML_dux_dxl(i,j) + coef2_xz_2 * PML_dux_dxl_old(i,j)
            rmemory_dux_dz_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dz_prime(i,j,ispec_PML,2) + &
                                                    coef1_xz_2 * PML_dux_dzl(i,j) + coef2_xz_2 * PML_dux_dzl_old(i,j)
            rmemory_duz_dx_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dx_prime(i,j,ispec_PML,2) + &
                                                    coef1_xz_2 * PML_duz_dxl(i,j) + coef2_xz_2 * PML_duz_dxl_old(i,j)
            rmemory_duz_dz_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dz_prime(i,j,ispec_PML,2) + &
                                                    coef1_xz_2 * PML_duz_dzl(i,j) + coef2_xz_2 * PML_duz_dzl_old(i,j)
          else
            rmemory_dux_dx_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dx_prime(i,j,ispec_PML,2) + &
                                                    coef1_xz_2 * time_n * PML_dux_dxl(i,j) + &
                                                    coef2_xz_2 * time_nsub1 * PML_dux_dxl_old(i,j)
            rmemory_dux_dz_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dz_prime(i,j,ispec_PML,2) + &
                                                    coef1_xz_2 * time_n * PML_dux_dzl(i,j) + &
                                                    coef2_xz_2 * time_nsub1 * PML_dux_dzl_old(i,j)
            rmemory_duz_dx_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dx_prime(i,j,ispec_PML,2) + &
                                                    coef1_xz_2 * time_n * PML_duz_dxl(i,j) + &
                                                    coef2_xz_2 * time_nsub1 * PML_duz_dxl_old(i,j)
            rmemory_duz_dz_prime(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dz_prime(i,j,ispec_PML,2) + &
                                                    coef1_xz_2 * time_n * PML_duz_dzl(i,j) + &
                                                    coef2_xz_2 * time_nsub1 * PML_duz_dzl_old(i,j)
          endif

        else
          ! non-rotated, element aligns with x/y/z-coordinates
          rmemory_dux_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_dux_dx(i,j,ispec_PML,1) + &
                                            coef1_zx_1 * PML_dux_dxl(i,j) + coef2_zx_1 * PML_dux_dxl_old(i,j)
          rmemory_duz_dx(i,j,ispec_PML,1) = coef0_zx_1 * rmemory_duz_dx(i,j,ispec_PML,1) + &
                                            coef1_zx_1 * PML_duz_dxl(i,j) + coef2_zx_1 * PML_duz_dxl_old(i,j)
          if (singularity_type_zx == 0) then
            rmemory_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + &
                                              coef1_zx_2 * PML_dux_dxl(i,j) + coef2_zx_2 * PML_dux_dxl_old(i,j)
            rmemory_duz_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + &
                                              coef1_zx_2 * PML_duz_dxl(i,j) + coef2_zx_2 * PML_duz_dxl_old(i,j)
          else
            rmemory_dux_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + &
                                              coef1_zx_2 * time_n * PML_dux_dxl(i,j) + &
                                              coef2_zx_2 * time_nsub1 * PML_dux_dxl_old(i,j)
            rmemory_duz_dx(i,j,ispec_PML,2) = coef0_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + &
                                              coef1_zx_2 * time_n * PML_duz_dxl(i,j) + &
                                              coef2_zx_2 * time_nsub1 * PML_duz_dxl_old(i,j)
          endif

          rmemory_dux_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_dux_dz(i,j,ispec_PML,1) + &
                                            coef1_xz_1 * PML_dux_dzl(i,j) + coef2_xz_1 * PML_dux_dzl_old(i,j)
          rmemory_duz_dz(i,j,ispec_PML,1) = coef0_xz_1 * rmemory_duz_dz(i,j,ispec_PML,1) + &
                                            coef1_xz_1 * PML_duz_dzl(i,j) + coef2_xz_1 * PML_duz_dzl_old(i,j)
          if (singularity_type_xz == 0) then
            rmemory_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dz(i,j,ispec_PML,2) + &
                                              coef1_xz_2 * PML_dux_dzl(i,j) + coef2_xz_2 * PML_dux_dzl_old(i,j)
            rmemory_duz_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dz(i,j,ispec_PML,2) + &
                                              coef1_xz_2 * PML_duz_dzl(i,j) + coef2_xz_2 * PML_duz_dzl_old(i,j)
          else
            rmemory_dux_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_dux_dz(i,j,ispec_PML,2) + &
                                              coef1_xz_2 * time_n * PML_dux_dzl(i,j) + &
                                              coef2_xz_2 * time_nsub1 * PML_dux_dzl_old(i,j)
            rmemory_duz_dz(i,j,ispec_PML,2) = coef0_xz_2 * rmemory_duz_dz(i,j,ispec_PML,2) + &
                                              coef1_xz_2 * time_n * PML_duz_dzl(i,j) + &
                                              coef2_xz_2 * time_nsub1 * PML_duz_dzl_old(i,j)
          endif
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
        if (singularity_type_zx == 0) then
          rmemory_dux_dx_LDDRK(i,j,ispec_PML,2) = ALPHA_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,2) + &
                                                  deltat * (-bb_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + PML_dux_dxl(i,j))
          rmemory_dux_dx(i,j,ispec_PML,2) = rmemory_dux_dx(i,j,ispec_PML,2) + &
                                            BETA_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,2)
          rmemory_duz_dx_LDDRK(i,j,ispec_PML,2) = ALPHA_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,2) + &
                                                  deltat * (-bb_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + PML_duz_dxl(i,j))
          rmemory_duz_dx(i,j,ispec_PML,2) = rmemory_duz_dx(i,j,ispec_PML,2) + &
                                            BETA_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,2)
        else
          rmemory_dux_dx_LDDRK(i,j,ispec_PML,2) = ALPHA_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,2) + &
                deltat * (-bb_zx_2 * rmemory_dux_dx(i,j,ispec_PML,2) + PML_dux_dxl(i,j) * time_n)
          rmemory_dux_dx(i,j,ispec_PML,2) = rmemory_dux_dx(i,j,ispec_PML,2) + &
                                            BETA_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML,2)

          rmemory_duz_dx_LDDRK(i,j,ispec_PML,2) = ALPHA_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,2) + &
                deltat * (-bb_zx_2 * rmemory_duz_dx(i,j,ispec_PML,2) + PML_duz_dxl(i,j) * time_n)
          rmemory_duz_dx(i,j,ispec_PML,2) = rmemory_duz_dx(i,j,ispec_PML,2) + &
                                            BETA_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML,2)
        endif

        rmemory_dux_dz_LDDRK(i,j,ispec_PML,1) = ALPHA_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,1) + &
                                                deltat * (-bb_xz_1 * rmemory_dux_dz(i,j,ispec_PML,1) + PML_dux_dzl(i,j))
        rmemory_dux_dz(i,j,ispec_PML,1) = rmemory_dux_dz(i,j,ispec_PML,1) + &
                                          BETA_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,1)
        rmemory_duz_dz_LDDRK(i,j,ispec_PML,1) = ALPHA_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,1) + &
                                                deltat * (-bb_xz_1 * rmemory_duz_dz(i,j,ispec_PML,1) + PML_duz_dzl(i,j))
        rmemory_duz_dz(i,j,ispec_PML,1) = rmemory_duz_dz(i,j,ispec_PML,1) + &
                                          BETA_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,1)
        if (singularity_type_xz == 0) then
          rmemory_dux_dz_LDDRK(i,j,ispec_PML,2) = ALPHA_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,2) + &
                                                  deltat * (-bb_xz_2 * rmemory_dux_dz(i,j,ispec_PML,2) + PML_dux_dzl(i,j))
          rmemory_dux_dz(i,j,ispec_PML,2) = rmemory_dux_dz(i,j,ispec_PML,2) + &
                                            BETA_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,2)

          rmemory_duz_dz_LDDRK(i,j,ispec_PML,2) = ALPHA_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,2) + &
                                                  deltat * (-bb_xz_2 * rmemory_duz_dz(i,j,ispec_PML,2) + PML_duz_dzl(i,j))
          rmemory_duz_dz(i,j,ispec_PML,2) = rmemory_duz_dz(i,j,ispec_PML,2) + &
                                            BETA_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,2)
        else
          rmemory_dux_dz_LDDRK(i,j,ispec_PML,2) = ALPHA_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,2) + &
                deltat * (-bb_xz_2 * rmemory_dux_dz(i,j,ispec_PML,2) + PML_dux_dzl(i,j) * time_n)
          rmemory_dux_dz(i,j,ispec_PML,2) = rmemory_dux_dz(i,j,ispec_PML,2) + &
                                            BETA_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML,2)

          rmemory_duz_dz_LDDRK(i,j,ispec_PML,2) = ALPHA_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,2) + &
                deltat * (-bb_xz_2 * rmemory_duz_dz(i,j,ispec_PML,2) + PML_duz_dzl(i,j) * time_n)
          rmemory_duz_dz(i,j,ispec_PML,2) = rmemory_duz_dz(i,j,ispec_PML,2) + &
                                            BETA_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML,2)
        endif

      case default
        stop 'Time stepping scheme not implemented yet for elastic PML memory variable update'
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

  ! AXISYM enforces zero derivatives on axis
  if (AXISYM) then
    if (is_on_the_axis(ispec)) then
      ! d_uz/dr=0 on the axis
      ! i == 1
      do j = 1,NGLLZ
        rmemory_duz_dx(1,j,ispec_PML,1) = 0._CUSTOM_REAL
        rmemory_duz_dx(1,j,ispec_PML,2) = 0._CUSTOM_REAL
        if (time_stepping_scheme == 2) then
          rmemory_duz_dx_LDDRK(1,j,ispec_PML,1) = 0._CUSTOM_REAL
          rmemory_duz_dx_LDDRK(1,j,ispec_PML,2) = 0._CUSTOM_REAL
        endif
        if (ROTATE_PML_ACTIVATE) then
          rmemory_duz_dx_prime(1,j,ispec_PML,1) = 0._CUSTOM_REAL
          rmemory_duz_dx_prime(1,j,ispec_PML,2) = 0._CUSTOM_REAL
        endif
      enddo
    endif
  endif

  end subroutine pml_compute_memory_variables_elastic

