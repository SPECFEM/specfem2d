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

  subroutine compute_coef_convolution(bb,deltat,coef0,coef1,coef2)

  ! compute coefficient used in second order convolution scheme, from
  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-medium PML for vector FETD with modified basis functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

  use constants, only: CUSTOM_REAL

  implicit none

  double precision :: bb,coef0,coef1,coef2
  double precision :: deltat

  coef0 = exp(- bb * deltat)

  if (abs(bb) > 1e-5_CUSTOM_REAL) then
    coef1 = (1.0_CUSTOM_REAL - exp(-bb * deltat * 0.5_CUSTOM_REAL)) / bb
    coef2 = (1.0_CUSTOM_REAL - exp(-bb * deltat * 0.5_CUSTOM_REAL)) * exp(-bb * deltat * 0.5_CUSTOM_REAL) / bb
  else
    coef1 = deltat * 0.5_CUSTOM_REAL
    coef2 = deltat * 0.5_CUSTOM_REAL
  endif

  end subroutine compute_coef_convolution

!========================================================================

  subroutine lik_parameter_computation(timeval,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                      CPML_region_local,index_ik,A_0,A_1,A_2,singularity_type_2,bb_1,bb_2, &
                                      coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2)

  use constants, only: CUSTOM_REAL,CPML_X_ONLY,CPML_Z_ONLY,CPML_XZ_ONLY

  implicit none

  double precision, intent(in) :: timeval
  double precision, intent(in) :: deltat
  double precision, intent(in) :: kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z
  integer, intent(in) :: CPML_region_local,index_ik

  double precision, intent(out) :: A_0,A_1,A_2
  double precision, intent(out) :: coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2
  integer, intent(out) :: singularity_type_2

  ! local variables
  double precision :: bar_A_0,bar_A_1,bar_A_2,alpha_0,bb_1,bb_2
  integer :: CPML_X_ONLY_TEMP,CPML_Z_ONLY_TEMP,CPML_XZ_ONLY_TEMP

  logical,parameter :: FIRST_ORDER_CONVOLUTION = .false.

  if (index_ik == 13) then
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Z_ONLY_TEMP = CPML_Z_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
  else if (index_ik == 31) then
    CPML_X_ONLY_TEMP = CPML_Z_ONLY
    CPML_Z_ONLY_TEMP = CPML_X_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
  else
    stop 'In lik_parameter_computation index_ik must be equal to 13 or 31'
  endif

  if (CPML_region_local == CPML_XZ_ONLY_TEMP) then
  !----------------A0-------------------------
    bar_A_0 = kappa_x / kappa_z
    A_0 = bar_A_0

    if (abs(alpha_x-beta_z) >= 1.e-5_CUSTOM_REAL) then
      !----------------A1,2-------------------------
      bar_A_1 = - bar_A_0 * (alpha_x - alpha_z) * (alpha_x - beta_x) / (alpha_x-beta_z)
      bar_A_2 = - bar_A_0 * (beta_z - alpha_z) * (beta_z - beta_x) / (beta_z-alpha_x)

      A_1 = bar_A_1
      A_2 = bar_A_2

      singularity_type_2 = 0 ! 0 means no singularity, 1 means first order singularity

    else if (abs(alpha_x-beta_z) < 1.e-5_CUSTOM_REAL) then
      !----------------A1,2,3-------------------------
      alpha_0 = max(alpha_x,beta_z)

      bar_A_1 = bar_A_0 * (-2._CUSTOM_REAL * alpha_0 + (alpha_z + beta_x))
      bar_A_2 = bar_A_0 * (alpha_0 - alpha_z) * (alpha_0-beta_x)

      A_1 = bar_A_1 + timeval * bar_A_2
      A_2 = -bar_A_2

      singularity_type_2 = 1 ! 0 means no singularity, 1 means first order singularity

    else
      stop 'error in lik_parameter_computation'
    endif

  else if (CPML_region_local == CPML_X_ONLY_TEMP) then
  !----------------A0-------------------------
    bar_A_0 = kappa_x
    A_0 = bar_A_0
  !----------------A1,2,3-------------------------
    bar_A_1 = - bar_A_0 * (alpha_x - beta_x)
    bar_A_2 = 0._CUSTOM_REAL

    A_1 = bar_A_1
    A_2 = bar_A_2

    singularity_type_2 = 0 ! 0 means no singularity, 1 means first order singularity

  else if (CPML_region_local == CPML_Z_ONLY_TEMP) then
  !----------------A0-------------------------
    bar_A_0 = 1._CUSTOM_REAL / kappa_z
    A_0 = bar_A_0

    !----------------A1,2,3-------------------------
    bar_A_1 = 0._CUSTOM_REAL
    bar_A_2 = - bar_A_0 * (beta_z - alpha_z)

    A_1 = bar_A_1
    A_2 = bar_A_2

    singularity_type_2 = 0 ! 0 means no singularity, 1 means first order singularity
  endif

  bb_1 = alpha_x
  call compute_coef_convolution(bb_1,deltat,coef0_1,coef1_1,coef2_1)

  bb_2 = beta_z
  call compute_coef_convolution(bb_2,deltat,coef0_2,coef1_2,coef2_2)

 end subroutine lik_parameter_computation

!========================================================================

  subroutine l_parameter_computation(timeval,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                     CPML_region_local,A_0,A_1,A_2,A_3,A_4,singularity_type, &
                                     bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2)

  use constants, only: CUSTOM_REAL,CPML_X_ONLY,CPML_Z_ONLY,CPML_XZ_ONLY

  implicit none

  double precision, intent(in) :: timeval
  double precision, intent(in) :: deltat
  double precision, intent(in) :: kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z
  integer, intent(in) :: CPML_region_local

  double precision, intent(out) :: A_0, A_1, A_2, A_3, A_4
  double precision, intent(out) :: bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2
  integer, intent(out) :: singularity_type

  ! local variables
  double precision :: bar_A_0, bar_A_1, bar_A_2, bar_A_3, bar_A_4
  double precision :: alpha_0, beta_xyz_1, beta_xyz_2

  beta_xyz_1 = beta_x + beta_z
  beta_xyz_2 = beta_x * beta_z

  if (CPML_region_local == CPML_XZ_ONLY) then
    bar_A_0 = kappa_x * kappa_z
    bar_A_1 = bar_A_0 * (beta_x + beta_z - alpha_x - alpha_z)
    bar_A_2 = bar_A_0 * (beta_x - alpha_x) * (beta_z - alpha_z - alpha_x) &
            - bar_A_0 * (beta_z - alpha_z) * alpha_z

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2

    beta_xyz_1 = beta_x + beta_z
    beta_xyz_2 = beta_x * beta_z

    if (abs( alpha_x - alpha_z ) >= 1.e-5_CUSTOM_REAL) then
      bar_A_3 = bar_A_0 * alpha_x**2 * (beta_x - alpha_x) * (beta_z - alpha_x) / (alpha_z - alpha_x)
      bar_A_4 = bar_A_0 * alpha_z**2 * (beta_x - alpha_z) * (beta_z - alpha_z) / (alpha_x - alpha_z)

      A_3 = bar_A_3
      A_4 = bar_A_4

      singularity_type = 0
    else if (abs( alpha_x - alpha_z ) < 1.e-5_CUSTOM_REAL) then
      alpha_0 = alpha_x
      bar_A_3 = bar_A_0 * (- 4._CUSTOM_REAL * alpha_0**3  &
                           + 3._CUSTOM_REAL * alpha_0**2 * beta_xyz_1 - 2._CUSTOM_REAL * alpha_0 * beta_xyz_2)
      bar_A_4 = bar_A_0 * alpha_0**2 * (beta_x - alpha_0) * (beta_z - alpha_0)

      A_3 = bar_A_3 + timeval * bar_A_4
      A_4 = -bar_A_4

      singularity_type = 1
    endif
  else if (CPML_region_local == CPML_X_ONLY) then
    bar_A_0 = kappa_x
    bar_A_1 = bar_A_0 * (beta_x - alpha_x)
    bar_A_2 = - bar_A_0 * alpha_x * (beta_x - alpha_x)

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2

    bar_A_3 = bar_A_0 * alpha_x**2 * (beta_x - alpha_x)
    bar_A_4 = 0._CUSTOM_REAL

    A_3 = bar_A_3
    A_4 = bar_A_4

    singularity_type = 0
  else if (CPML_region_local == CPML_Z_ONLY) then
    bar_A_0 = kappa_z
    bar_A_1 = bar_A_0 * (beta_z - alpha_z)
    bar_A_2 = - bar_A_0 * alpha_z * (beta_z - alpha_z)

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2

    bar_A_3 = 0._CUSTOM_REAL
    bar_A_4 = bar_A_0 * alpha_z**2 * (beta_z - alpha_z)

    A_3 = bar_A_3
    A_4 = bar_A_4

    singularity_type = 0
  endif

  bb_1 = alpha_x
  call compute_coef_convolution(bb_1,deltat,coef0_1,coef1_1,coef2_1)

  bb_2 = alpha_z
  call compute_coef_convolution(bb_2,deltat,coef0_2,coef1_2,coef2_2)

  end subroutine l_parameter_computation

!=====================================================================

  subroutine rebuild_value_on_PML_interface_acoustic(it,b_potential_acoustic,b_potential_dot_acoustic)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: nspec,nglob,ispec_is_acoustic,any_acoustic,ibool,nglob_interface,point_interface

  ! PML arrays
  use specfem_par, only: ispec_is_PML,pml_interface_history_potential_dot,pml_interface_history_potential

  implicit none

  integer,intent(in) :: it

  real(kind=CUSTOM_REAL), dimension(nglob),intent(out) :: b_potential_acoustic,b_potential_dot_acoustic

  !local variables
  integer :: i,j,ispec

  if (any_acoustic) then
    do ispec = 1,nspec
      if (ispec_is_acoustic(ispec) .and. ispec_is_PML(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            b_potential_dot_acoustic(ibool(i,j,ispec)) = 0._CUSTOM_REAL
            b_potential_acoustic(ibool(i,j,ispec)) = 0._CUSTOM_REAL
          enddo
        enddo
      endif
    enddo

    if (nglob_interface > 0) then
      do i = 1, nglob_interface
        b_potential_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot(i,it)
        b_potential_acoustic(point_interface(i)) = pml_interface_history_potential(i,it)
      enddo
    endif
  endif

  end subroutine rebuild_value_on_PML_interface_acoustic

!=====================================================================

  subroutine rebuild_value_on_PML_interface_acoustic_accel(it,b_potential_dot_dot_acoustic)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: any_acoustic,nglob,nglob_interface,point_interface

  ! PML arrays
  use specfem_par, only: pml_interface_history_potential_dot_dot

  implicit none

  integer,intent(in) :: it

  real(kind=CUSTOM_REAL), dimension(nglob),intent(out) :: b_potential_dot_dot_acoustic

  !local variables
  integer :: i

  if (any_acoustic .and. nglob_interface > 0) then
    do i = 1, nglob_interface
      b_potential_dot_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot_dot(i,it)
    enddo
  endif

  end subroutine rebuild_value_on_PML_interface_acoustic_accel


!========================================================================

  subroutine rebuild_value_on_PML_interface_viscoelastic(it)

  use constants, only: NGLLX,NGLLZ,CUSTOM_REAL

  use specfem_par, only: nspec,ispec_is_elastic,any_elastic,ibool,nglob_interface,point_interface, &
                         b_veloc_elastic,b_displ_elastic
  ! PML arrays
  use specfem_par, only: ispec_is_PML,pml_interface_history_veloc,pml_interface_history_displ

  implicit none

  integer,intent(in) :: it

  !local variables
  integer :: i,j,ispec

  if (any_elastic) then
    do ispec = 1,nspec
      if (ispec_is_elastic(ispec) .and. ispec_is_PML(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            b_veloc_elastic(:,ibool(i,j,ispec)) = 0.0_CUSTOM_REAL
            b_displ_elastic(:,ibool(i,j,ispec)) = 0.0_CUSTOM_REAL
          enddo
        enddo
      endif
    enddo

    if (nglob_interface > 0) then
      do i = 1, nglob_interface
        b_veloc_elastic(1,point_interface(i)) = pml_interface_history_veloc(1,i,it)
        b_veloc_elastic(2,point_interface(i)) = pml_interface_history_veloc(2,i,it)

        b_displ_elastic(1,point_interface(i)) = pml_interface_history_displ(1,i,it)
        b_displ_elastic(2,point_interface(i)) = pml_interface_history_displ(2,i,it)
      enddo
    endif
  endif

  end subroutine rebuild_value_on_PML_interface_viscoelastic

!========================================================================

  subroutine rebuild_value_on_PML_interface_viscoelastic_accel(it)

  use constants, only: NGLLX,NGLLZ,CUSTOM_REAL

  use specfem_par, only: any_elastic,nglob_interface,point_interface, &
                         b_accel_elastic
  ! PML arrays
  use specfem_par, only: pml_interface_history_accel

  implicit none

  integer,intent(in) :: it

  !local variables
  integer :: i

  if (any_elastic .and. nglob_interface > 0) then
    do i = 1, nglob_interface
      b_accel_elastic(1,point_interface(i)) = pml_interface_history_accel(1,i,it)
      b_accel_elastic(2,point_interface(i)) = pml_interface_history_accel(2,i,it)
    enddo
  endif

  end subroutine rebuild_value_on_PML_interface_viscoelastic_accel

!========================================================================

  subroutine pml_boundary_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                   potential_acoustic,potential_acoustic_old)

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the potential
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the potential.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: nglob,ibool,nelemabs,codeabs,anyabs,numabs,ispec_is_PML

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob),intent(inout) :: potential_dot_dot_acoustic,potential_dot_acoustic, &
                                                            potential_acoustic,potential_acoustic_old

  ! local parameters
  integer :: i,j,ispecabs,ispec,iglob,ibegin,iend

  ! checks if anything to do
  if (.not. anyabs) return

  ! set Dirichelet boundary condition on outer boundary of CFS-PML
  do ispecabs = 1,nelemabs
    ispec = numabs(ispecabs)

    if (ispec_is_PML(ispec)) then
!--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then
        i = 1
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          potential_acoustic_old(iglob) = 0._CUSTOM_REAL
          potential_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
        enddo
      endif
!--- right absorbing boundary
      if (codeabs(IEDGE2,ispecabs)) then
        i = NGLLX
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          potential_acoustic_old(iglob) = 0._CUSTOM_REAL
          potential_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
        enddo
      endif
!--- bottom absorbing boundary
      if (codeabs(IEDGE1,ispecabs)) then
        j = 1
        ibegin = 1
        iend = NGLLX
        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          potential_acoustic_old(iglob) = 0._CUSTOM_REAL
          potential_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
        enddo
      endif
!--- top absorbing boundary
      if (codeabs(IEDGE3,ispecabs)) then
        j = NGLLZ
! exclude corners to make sure there is no contradiction on the normal
        ibegin = 1
        iend = NGLLX
        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          potential_acoustic_old(iglob) = 0._CUSTOM_REAL
          potential_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
          potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
        enddo
      endif  !  end of top absorbing boundary
    endif ! end of ispec_is_PML
  enddo ! end specabs loop

  end subroutine pml_boundary_acoustic

!========================================================================

  subroutine pml_boundary_elastic(accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: nglob,ibool,nelemabs,codeabs,anyabs,numabs,ispec_is_PML,nspec_PML

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: accel_elastic,veloc_elastic,displ_elastic,displ_elastic_old

  ! local parameters
  integer :: i,j,ispecabs,ispec,iglob,ibegin,iend

  ! checks if anything to do
  if (.not. anyabs) return
  if (nspec_PML == 0) return

  !set Dirichlet boundary condition on outer boundary of PML

  ! we have to put Dirichlet on the boundary of the PML
  do ispecabs = 1,nelemabs
    ispec = numabs(ispecabs)

    if (ispec_is_PML(ispec)) then

!--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then
        i = 1
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          displ_elastic_old(:,iglob) = 0._CUSTOM_REAL
          displ_elastic(:,iglob) = 0._CUSTOM_REAL
          veloc_elastic(:,iglob) = 0._CUSTOM_REAL
          accel_elastic(:,iglob) = 0._CUSTOM_REAL
        enddo
      endif

!--- right absorbing boundary
      if (codeabs(IEDGE2,ispecabs)) then
        i = NGLLX
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          displ_elastic_old(:,iglob) = 0._CUSTOM_REAL
          displ_elastic(:,iglob) = 0._CUSTOM_REAL
          veloc_elastic(:,iglob) = 0._CUSTOM_REAL
          accel_elastic(:,iglob) = 0._CUSTOM_REAL
        enddo
      endif

!--- bottom absorbing boundary
      if (codeabs(IEDGE1,ispecabs)) then
        j = 1
        ibegin = 1
        iend = NGLLX
        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          displ_elastic_old(:,iglob) = 0._CUSTOM_REAL
          displ_elastic(:,iglob) = 0._CUSTOM_REAL
          veloc_elastic(:,iglob) = 0._CUSTOM_REAL
          accel_elastic(:,iglob) = 0._CUSTOM_REAL
        enddo
      endif

!--- top absorbing boundary
      if (codeabs(IEDGE3,ispecabs)) then
        j = NGLLZ
        ibegin = 1
        iend = NGLLX
        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          displ_elastic_old(:,iglob) = 0._CUSTOM_REAL
          displ_elastic(:,iglob) = 0._CUSTOM_REAL
          veloc_elastic(:,iglob) = 0._CUSTOM_REAL
          accel_elastic(:,iglob) = 0._CUSTOM_REAL
        enddo
      endif

    endif ! end of ispec_is_PML
  enddo ! end specabs loop

  end subroutine pml_boundary_elastic
