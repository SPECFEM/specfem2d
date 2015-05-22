
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

 subroutine compute_coef_convolution(bb,deltat,coef0,coef1,coef2)
  ! compute coefficient used in second order convolution scheme, from
  ! second-order accurate convolution term calculation from equation (21) of
  ! Shumin Wang, Robert Lee, and Fernando L. Teixeira,
  ! Anisotropic-Medium PML for Vector FETD With Modified Basis Functions,
  ! IEEE Transactions on Antennas and Propagation, vol. 54, no. 1, (2006)

   implicit none
   include "constants.h"

   double precision :: bb,coef0,coef1,coef2
   double precision :: deltat

   coef0 = exp(- bb * deltat)

   if(  abs(bb) > 1e-5_CUSTOM_REAL ) then
     coef1 = (1._CUSTOM_REAL - exp(-bb * deltat / 2._CUSTOM_REAL)) / bb
     coef2 = (1._CUSTOM_REAL - exp(-bb * deltat / 2._CUSTOM_REAL)) * exp(-bb * deltat / 2._CUSTOM_REAL) / bb
   else
     coef1 = deltat / 2._CUSTOM_REAL
     coef2 = deltat / 2._CUSTOM_REAL
   endif

 end subroutine compute_coef_convolution

!========================================================================

 subroutine lik_parameter_computation(timeval,deltat,kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z, &
                                      CPML_region_local,index_ik,A_0,A_1,A_2,singularity_type_2,bb_1,bb_2, &
                                      coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2)
  implicit none
  include "constants.h"

  double precision, intent(in) :: timeval
  double precision :: deltat
  double precision, intent(in) :: kappa_x,beta_x,alpha_x,kappa_z,beta_z,alpha_z
  integer, intent(in) :: CPML_region_local,index_ik

  double precision, intent(out) :: A_0,A_1,A_2
  double precision, intent(out) :: coef0_1,coef1_1,coef2_1,coef0_2,coef1_2,coef2_2
  integer, intent(out) :: singularity_type_2

  ! local variables
  double precision :: bar_A_0,bar_A_1,bar_A_2,alpha_0,bb_1,bb_2
  integer :: CPML_X_ONLY_TEMP,CPML_Z_ONLY_TEMP,CPML_XZ_ONLY_TEMP

  logical,parameter :: FIRST_ORDER_CONVOLUTION = .false.

  if( index_ik == 13 ) then
    CPML_X_ONLY_TEMP = CPML_X_ONLY
    CPML_Z_ONLY_TEMP = CPML_Z_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
  else if( index_ik == 31 ) then
    CPML_X_ONLY_TEMP = CPML_Z_ONLY
    CPML_Z_ONLY_TEMP = CPML_X_ONLY
    CPML_XZ_ONLY_TEMP = CPML_XZ_ONLY
  else
    stop 'In lik_parameter_computation index_ik must be equal to 13 or 31'
  endif

  if( CPML_region_local == CPML_XZ_ONLY_TEMP ) then
  !----------------A0-------------------------
    bar_A_0 = kappa_x / kappa_z
    A_0 = bar_A_0

    if( abs(alpha_x-beta_z) >= 1.e-5_CUSTOM_REAL ) then
      !----------------A1,2-------------------------
      bar_A_1 = - bar_A_0 * (alpha_x - alpha_z) * (alpha_x - beta_x) / (alpha_x-beta_z)
      bar_A_2 = - bar_A_0 * (beta_z - alpha_z) * (beta_z - beta_x) / (beta_z-alpha_x)

      A_1 = bar_A_1
      A_2 = bar_A_2

      singularity_type_2 = 0 ! 0 means no singularity, 1 means first order singularity

    else if( abs(alpha_x-beta_z) < 1.e-5_CUSTOM_REAL ) then
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

  else if( CPML_region_local == CPML_X_ONLY_TEMP ) then
  !----------------A0-------------------------
    bar_A_0 = kappa_x
    A_0 = bar_A_0
  !----------------A1,2,3-------------------------
    bar_A_1 = - bar_A_0 * (alpha_x - beta_x)
    bar_A_2 = 0._CUSTOM_REAL

    A_1 = bar_A_1
    A_2 = bar_A_2

    singularity_type_2 = 0 ! 0 means no singularity, 1 means first order singularity

  else if( CPML_region_local == CPML_Z_ONLY_TEMP ) then
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
                                    CPML_region_local,A_0,A_1,A_2,A_3,A_4,singularity_type,&
                                    bb_1,coef0_1,coef1_1,coef2_1,bb_2,coef0_2,coef1_2,coef2_2)

  implicit none
  include "constants.h"

  double precision :: timeval
  double precision :: deltat
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

  if(  CPML_region_local == CPML_XZ_ONLY ) then
    bar_A_0 = kappa_x * kappa_z
    bar_A_1 = bar_A_0 * (beta_x + beta_z - alpha_x - alpha_z)
    bar_A_2 = bar_A_0 * (beta_x - alpha_x) * (beta_z - alpha_z - alpha_x) &
            - bar_A_0 * (beta_z - alpha_z) * alpha_z

    A_0 = bar_A_0
    A_1 = bar_A_1
    A_2 = bar_A_2

    beta_xyz_1 = beta_x + beta_z
    beta_xyz_2 = beta_x * beta_z

    if(  abs( alpha_x - alpha_z ) >= 1.e-5_CUSTOM_REAL ) then
      bar_A_3 = bar_A_0 * alpha_x**2 * (beta_x - alpha_x) * (beta_z - alpha_x) / (alpha_z - alpha_x)
      bar_A_4 = bar_A_0 * alpha_z**2 * (beta_x - alpha_z) * (beta_z - alpha_z) / (alpha_x - alpha_z)

      A_3 = bar_A_3
      A_4 = bar_A_4

      singularity_type = 0
    else if ( abs( alpha_x - alpha_z ) < 1.e-5_CUSTOM_REAL ) then
      alpha_0 = alpha_x
      bar_A_3 = bar_A_0 * (- 4._CUSTOM_REAL * alpha_0**3  &
                           + 3._CUSTOM_REAL * alpha_0**2 * beta_xyz_1 - 2._CUSTOM_REAL * alpha_0 * beta_xyz_2)
      bar_A_4 = bar_A_0 * alpha_0**2 * (beta_x - alpha_0) * (beta_z - alpha_0)

      A_3 = bar_A_3 + timeval * bar_A_4
      A_4 = -bar_A_4

      singularity_type = 1
    endif
  else if ( CPML_region_local == CPML_X_ONLY ) then
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
  else if ( CPML_region_local == CPML_Z_ONLY ) then
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

 subroutine rebuild_value_on_PML_interface_acoustic(it)

  use specfem_par, only: nspec,acoustic,any_acoustic,is_pml,ibool,nglob_interface,point_interface, &
                         b_potential_dot_acoustic,b_potential_acoustic,&
                         pml_interface_history_potential_dot,pml_interface_history_potential
  implicit none
  include "constants.h"
  integer :: it

  !local variables
  integer :: i,j,ispec

  do ispec = 1,nspec
    do j = 1, NGLLZ
      do i = 1, NGLLX
        if( acoustic(ispec) .and. is_pml(ispec) ) then
          b_potential_dot_acoustic(ibool(i,j,ispec)) = 0._CUSTOM_REAL
          b_potential_acoustic(ibool(i,j,ispec)) = 0._CUSTOM_REAL
        endif
      enddo
    enddo
  enddo

  if( any_acoustic .and. nglob_interface > 0 ) then
    do i = 1, nglob_interface
      b_potential_dot_acoustic(point_interface(i)) = pml_interface_history_potential_dot(i,it)
      b_potential_acoustic(point_interface(i)) = pml_interface_history_potential(i,it)
    enddo
  endif

 end subroutine rebuild_value_on_PML_interface_acoustic
!========================================================================

 subroutine rebuild_value_on_PML_interface_viscoelastic(it)

  use specfem_par, only: nspec,elastic,any_elastic,is_pml,ibool,nglob_interface,point_interface, &
                         b_veloc_elastic,b_displ_elastic,&
                         pml_interface_history_veloc,pml_interface_history_displ
  implicit none
  include "constants.h"
  integer :: it

  !local variables
  integer :: i,j,ispec

  do ispec = 1,nspec
    do j = 1, NGLLZ
      do i = 1, NGLLX
        if( elastic(ispec) .and. is_pml(ispec) ) then
           b_veloc_elastic(:,ibool(i,j,ispec)) = 0.
           b_displ_elastic(:,ibool(i,j,ispec)) = 0.
        endif
      enddo
    enddo
  enddo

  if(any_elastic .and. nglob_interface > 0) then
     do i = 1, nglob_interface
       b_veloc_elastic(1,point_interface(i)) = pml_interface_history_veloc(1,i,it)
       b_veloc_elastic(2,point_interface(i)) = pml_interface_history_veloc(2,i,it)
       b_veloc_elastic(3,point_interface(i)) = pml_interface_history_veloc(3,i,it)
       b_displ_elastic(1,point_interface(i)) = pml_interface_history_displ(1,i,it)
       b_displ_elastic(2,point_interface(i)) = pml_interface_history_displ(2,i,it)
       b_displ_elastic(3,point_interface(i)) = pml_interface_history_displ(3,i,it)
     enddo
  endif

 end subroutine rebuild_value_on_PML_interface_viscoelastic
!========================================================================
