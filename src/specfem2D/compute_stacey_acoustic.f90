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

  subroutine compute_stacey_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic)

! absorbing boundaries
! for Stacey paraxial absorbing conditions (more precisely: Sommerfeld in the case of a fluid) we implement them here

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,ZERO,ONE,TWO,TWO_THIRDS,IEDGE1,IEDGE2,IEDGE3,IEDGE4,USE_A_STRONG_FORMULATION_FOR_E1

  use specfem_par, only: nglob,num_abs_boundary_faces,anyabs,it,any_acoustic, &
                         ibool, &
                         abs_boundary_ispec,ispec_is_acoustic, &
                         codeabs,codeabs_corner, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         rho_vpstore, &
                         wxgll,wzgll, &
                         ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
                         ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
                         SAVE_FORWARD, &
                         ib_left,ib_right,ib_bottom,ib_top, &
                         b_absorb_acoustic_left,b_absorb_acoustic_right, &
                         b_absorb_acoustic_bottom,b_absorb_acoustic_top, &
                         STACEY_ABSORBING_CONDITIONS, &
                         ATTENUATION_VISCOACOUSTIC,dot_e1, &
                         NO_BACKWARD_RECONSTRUCTION

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob),intent(inout) :: potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: potential_dot_acoustic

  ! local parameters
  integer :: ispecabs,ispec,i,j,iglob
  integer :: ibegin,iend,jbegin,jend
  real(kind=CUSTOM_REAL) :: weight,xxi,zxi,xgamma,zgamma,jacobian1D
  ! material properties of the acoustic medium
  real(kind=CUSTOM_REAL) :: rho_vp

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. anyabs) return
  if (.not. any_acoustic) return

  do ispecabs = 1,num_abs_boundary_faces
    ispec = abs_boundary_ispec(ispecabs)

    ! Sommerfeld condition if acoustic
    if (ispec_is_acoustic(ispec)) then

      !--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then
        i = 1
        jbegin = ibegin_edge4(ispecabs)
        jend = iend_edge4(ispecabs)
        do j = jbegin,jend
          iglob = ibool(i,j,ispec)
          rho_vp = rho_vpstore(i,j,ispec)

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma ** 2 + zgamma ** 2)
          weight = jacobian1D * wzgll(j)

          ! adds absorbing boundary contribution
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob) * weight/rho_vp
          if (ATTENUATION_VISCOACOUSTIC .and. .not. USE_A_STRONG_FORMULATION_FOR_E1) &
                   dot_e1(iglob,:) = dot_e1(iglob,:) - potential_dot_acoustic(iglob) * weight/rho_vp

          if (SAVE_FORWARD) then
            ! saves contribution
            b_absorb_acoustic_left(j,ib_left(ispecabs),it) = potential_dot_acoustic(iglob) * weight/rho_vp
          endif
        enddo
      endif  !  end of left absorbing boundary

      !--- right absorbing boundary
      if (codeabs(IEDGE2,ispecabs)) then
        i = NGLLX
        jbegin = ibegin_edge2(ispecabs)
        jend = iend_edge2(ispecabs)
        do j = jbegin,jend
          iglob = ibool(i,j,ispec)
          rho_vp = rho_vpstore(i,j,ispec)

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma ** 2 + zgamma ** 2)
          weight = jacobian1D * wzgll(j)

          ! adds absorbing boundary contribution
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob) * weight/rho_vp
          if (ATTENUATION_VISCOACOUSTIC .and. .not. USE_A_STRONG_FORMULATION_FOR_E1) &
                 dot_e1(iglob,:) = dot_e1(iglob,:) - potential_dot_acoustic(iglob) * weight/rho_vp

          if (SAVE_FORWARD) then
            ! saves contribution
            b_absorb_acoustic_right(j,ib_right(ispecabs),it) = potential_dot_acoustic(iglob) * weight/rho_vp
          endif
        enddo
      endif  !  end of right absorbing boundary

      !--- bottom absorbing boundary
      if (codeabs(IEDGE1,ispecabs)) then
        j = 1
        ibegin = ibegin_edge1(ispecabs)
        iend = iend_edge1(ispecabs)
        ! exclude corners to make sure there is no contradiction on the normal
        if (codeabs_corner(1,ispecabs)) ibegin = 2
        if (codeabs_corner(2,ispecabs)) iend = NGLLX-1
        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          rho_vp = rho_vpstore(i,j,ispec)

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi ** 2 + zxi ** 2)
          weight = jacobian1D * wxgll(i)

          ! adds absorbing boundary contribution
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob) * weight/rho_vp
          if (ATTENUATION_VISCOACOUSTIC .and. .not. USE_A_STRONG_FORMULATION_FOR_E1) &
             dot_e1(iglob,:) = dot_e1(iglob,:) - potential_dot_acoustic(iglob) * weight/rho_vp

          if (SAVE_FORWARD) then
            ! saves contribution
            b_absorb_acoustic_bottom(i,ib_bottom(ispecabs),it) = potential_dot_acoustic(iglob) * weight/rho_vp
          endif
        enddo
      endif  !  end of bottom absorbing boundary

      !--- top absorbing boundary
      if (codeabs(IEDGE3,ispecabs)) then
        j = NGLLZ
        ibegin = ibegin_edge3(ispecabs)
        iend = iend_edge3(ispecabs)
        ! exclude corners to make sure there is no contradiction on the normal
        if (codeabs_corner(3,ispecabs)) ibegin = 2
        if (codeabs_corner(4,ispecabs)) iend = NGLLX-1
        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          rho_vp = rho_vpstore(i,j,ispec)

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi ** 2 + zxi ** 2)
          weight = jacobian1D * wxgll(i)

          ! adds absorbing boundary contribution
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob) * weight/rho_vp
          if (ATTENUATION_VISCOACOUSTIC .and. .not. USE_A_STRONG_FORMULATION_FOR_E1) &
               dot_e1(iglob,:) = dot_e1(iglob,:) - potential_dot_acoustic(iglob) * weight/rho_vp

          if (SAVE_FORWARD .and. (.not. NO_BACKWARD_RECONSTRUCTION)) then
            ! saves contribution
            b_absorb_acoustic_top(i,ib_top(ispecabs),it) = potential_dot_acoustic(iglob) * weight/rho_vp
          endif
        enddo
      endif  !  end of top absorbing boundary

    endif ! acoustic ispec
  enddo

  end subroutine compute_stacey_acoustic

!
!------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_acoustic_backward(b_potential_dot_dot_acoustic)

! absorbing boundaries
! uses contributions stored in forward simulation

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,ZERO,ONE,TWO,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: nglob,num_abs_boundary_faces,anyabs,it,NSTEP,any_acoustic, &
                         ibool,abs_boundary_ispec,ispec_is_acoustic, &
                         codeabs,codeabs_corner, &
                         ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
                         ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
                         ib_left,ib_right,ib_bottom,ib_top, &
                         b_absorb_acoustic_left,b_absorb_acoustic_right, &
                         b_absorb_acoustic_bottom,b_absorb_acoustic_top, &
                         STACEY_ABSORBING_CONDITIONS,NO_BACKWARD_RECONSTRUCTION

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob),intent(inout) :: b_potential_dot_dot_acoustic

  ! local parameters
  integer :: ispecabs,ispec,i,j,iglob
  integer :: ibegin,iend,jbegin,jend
  integer :: it_tmp

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. anyabs) return
  if (.not. any_acoustic) return
  if (NO_BACKWARD_RECONSTRUCTION) return

  ! time increment step
  it_tmp = NSTEP - it + 1

  do ispecabs = 1,num_abs_boundary_faces
    ispec = abs_boundary_ispec(ispecabs)

    ! Sommerfeld condition if acoustic
    if (ispec_is_acoustic(ispec)) then
      !--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then
        i = 1
        jbegin = ibegin_edge4(ispecabs)
        jend = iend_edge4(ispecabs)
        do j = jbegin,jend
          iglob = ibool(i,j,ispec)
          b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                                                b_absorb_acoustic_left(j,ib_left(ispecabs),it_tmp)

        enddo
      endif  !  end of left absorbing boundary

      !--- right absorbing boundary
      if (codeabs(IEDGE2,ispecabs)) then
        i = NGLLX
        jbegin = ibegin_edge2(ispecabs)
        jend = iend_edge2(ispecabs)
        do j = jbegin,jend
          iglob = ibool(i,j,ispec)
          ! adds (previously) stored contribution
          b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                                                b_absorb_acoustic_right(j,ib_right(ispecabs),it_tmp)
        enddo
      endif  !  end of right absorbing boundary

      !--- bottom absorbing boundary
      if (codeabs(IEDGE1,ispecabs)) then
        j = 1
        ibegin = ibegin_edge1(ispecabs)
        iend = iend_edge1(ispecabs)
        ! exclude corners to make sure there is no contradiction on the normal
        if (codeabs_corner(1,ispecabs)) ibegin = 2
        if (codeabs_corner(2,ispecabs)) iend = NGLLX-1
        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          ! adds (previously) stored contribution
          b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                                                b_absorb_acoustic_bottom(i,ib_bottom(ispecabs),it_tmp)
        enddo
      endif  !  end of bottom absorbing boundary

      !--- top absorbing boundary
      if (codeabs(IEDGE3,ispecabs)) then
        j = NGLLZ
        ibegin = ibegin_edge3(ispecabs)
        iend = iend_edge3(ispecabs)
        ! exclude corners to make sure there is no contradiction on the normal
        if (codeabs_corner(3,ispecabs)) ibegin = 2
        if (codeabs_corner(4,ispecabs)) iend = NGLLX-1
        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                                                b_absorb_acoustic_top(i,ib_top(ispecabs),it_tmp)
        enddo
      endif  !  end of top absorbing boundary
    endif ! acoustic ispec
  enddo

  end subroutine compute_stacey_acoustic_backward

