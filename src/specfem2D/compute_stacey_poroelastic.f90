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

!---------------------------------------------------------------------------------------------
!
! fluid part
!
!---------------------------------------------------------------------------------------------

  subroutine compute_stacey_poro_fluid(accelw_poroelastic,velocs_poroelastic,velocw_poroelastic)

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4,TWO,ZERO

  use specfem_par, only: it,STACEY_ABSORBING_CONDITIONS, &
                         anyabs,num_abs_boundary_faces,abs_boundary_ispec, &
                         ibool,ispec_is_poroelastic,any_poroelastic, &
                         codeabs,codeabs_corner, &
                         nglob_poroelastic, &
                         phistore,tortstore,rhoarraystore,vpIIstore,rho_vpstore,rho_vsstore, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         wxgll,wzgll, &
                         ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro,iend_edge3_poro, &
                         ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro, &
                         SIMULATION_TYPE,SAVE_FORWARD, &
                         b_absorb_poro_w_left,b_absorb_poro_w_right, &
                         b_absorb_poro_w_bottom,b_absorb_poro_w_top, &
                         ib_left,ib_right,ib_bottom,ib_top, &
                         NO_BACKWARD_RECONSTRUCTION

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic),intent(inout) :: accelw_poroelastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic),intent(in) :: velocs_poroelastic,velocw_poroelastic

  ! local parameters
  integer :: ispec,i,j,iglob
  integer :: ispecabs,ibegin,iend,jbegin,jend

  real(kind=CUSTOM_REAL) :: nx,nz,vx,vz,vn,vxf,vzf,vnf
  real(kind=CUSTOM_REAL) :: rho_vpI,rho_vpII,rho_vs
  real(kind=CUSTOM_REAL) :: tx,tz,weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: cpIl,cpIIl,csl
  real(kind=CUSTOM_REAL) :: phi,tort,rho_s,rho_f,rho_bar

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. anyabs) return
  if (.not. any_poroelastic) return

  ! absorbing boundaries
  do ispecabs = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(ispecabs)

    if (ispec_is_poroelastic(ispec)) then

      !--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then
        i = 1

        jbegin = ibegin_edge4_poro(ispecabs)
        jend = iend_edge4_poro(ispecabs)

        do j = jbegin,jend
          iglob = ibool(i,j,ispec)

          ! poroelastic material
          phi = phistore(i,j,ispec)
          tort = tortstore(i,j,ispec)

          rho_s = rhoarraystore(1,i,j,ispec)
          rho_f = rhoarraystore(2,i,j,ispec)

          cpIl = rho_vpstore(i,j,ispec) / rho_s
          cpIIl = vpIIstore(i,j,ispec)
          csl = rho_vsstore(i,j,ispec) / rho_s

          rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

          rho_vpI = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIl
          rho_vpII = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIIl
          rho_vs = rho_f/rho_bar*(rho_bar-rho_f*phi/tort)*csl

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = - zgamma / jacobian1D
          nz = + xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

          ! solid contribution
          vx = velocs_poroelastic(1,iglob)
          vz = velocs_poroelastic(2,iglob)

          ! fluid contribution
          vxf = velocw_poroelastic(1,iglob)
          vzf = velocw_poroelastic(2,iglob)

          vn = nx*vx+nz*vz
          vnf = nx*vxf+nz*vzf

          tx = rho_vpII*vnf*nx - rho_vs*(vx-vn*nx)
          tz = rho_vpII*vnf*nz - rho_vs*(vz-vn*nz)

          accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - tx*weight
          accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - tz*weight

          ! saves contribution
          if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. NO_BACKWARD_RECONSTRUCTION)) then
            b_absorb_poro_w_left(1,j,ib_left(ispecabs),it) = tx*weight
            b_absorb_poro_w_left(2,j,ib_left(ispecabs),it) = tz*weight
          endif

        enddo

      endif  !  end of left absorbing boundary

      !--- right absorbing boundary
      if (codeabs(IEDGE2,ispecabs)) then
        i = NGLLX

        jbegin = ibegin_edge2_poro(ispecabs)
        jend = iend_edge2_poro(ispecabs)

        do j = jbegin,jend
          iglob = ibool(i,j,ispec)

          ! poroelastic material
          phi = phistore(i,j,ispec)
          tort = tortstore(i,j,ispec)

          rho_s = rhoarraystore(1,i,j,ispec)
          rho_f = rhoarraystore(2,i,j,ispec)

          cpIl = rho_vpstore(i,j,ispec) / rho_s
          cpIIl = vpIIstore(i,j,ispec)
          csl = rho_vsstore(i,j,ispec) / rho_s

          rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

          rho_vpI = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIl
          rho_vpII = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIIl
          rho_vs = rho_f/rho_bar*(rho_bar-rho_f*phi/tort)*csl

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = + zgamma / jacobian1D
          nz = - xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

          ! solid contribution
          vx = velocs_poroelastic(1,iglob)
          vz = velocs_poroelastic(2,iglob)

          ! fluid contribution
          vxf = velocw_poroelastic(1,iglob)
          vzf = velocw_poroelastic(2,iglob)

          vn = nx*vx+nz*vz
          vnf = nx*vxf+nz*vzf

          tx = rho_vpII*vnf*nx - rho_vs*(vx-vn*nx)
          tz = rho_vpII*vnf*nz - rho_vs*(vz-vn*nz)

          accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - tx*weight
          accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - tz*weight

          ! saves contribution
          if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. NO_BACKWARD_RECONSTRUCTION)) then
            b_absorb_poro_w_right(1,j,ib_right(ispecabs),it) = tx*weight
            b_absorb_poro_w_right(2,j,ib_right(ispecabs),it) = tz*weight
          endif

        enddo

      endif  !  end of right absorbing boundary

      !--- bottom absorbing boundary
      if (codeabs(IEDGE1,ispecabs)) then
        j = 1

        ibegin = ibegin_edge1_poro(ispecabs)
        iend = iend_edge1_poro(ispecabs)

        ! exclude corners to make sure there is no contradiction on the normal
        if (codeabs_corner(1,ispecabs)) ibegin = 2
        if (codeabs_corner(2,ispecabs)) iend = NGLLX-1

        do i = ibegin,iend
          iglob = ibool(i,j,ispec)

          ! poroelastic material
          phi = phistore(i,j,ispec)
          tort = tortstore(i,j,ispec)

          rho_s = rhoarraystore(1,i,j,ispec)
          rho_f = rhoarraystore(2,i,j,ispec)

          cpIl = rho_vpstore(i,j,ispec) / rho_s
          cpIIl = vpIIstore(i,j,ispec)
          csl = rho_vsstore(i,j,ispec) / rho_s

          rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

          rho_vpI = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIl
          rho_vpII = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIIl
          rho_vs = rho_f/rho_bar*(rho_bar-rho_f*phi/tort)*csl

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = + zxi / jacobian1D
          nz = - xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

          ! solid contribution
          vx = velocs_poroelastic(1,iglob)
          vz = velocs_poroelastic(2,iglob)

          ! fluid contribution
          vxf = velocw_poroelastic(1,iglob)
          vzf = velocw_poroelastic(2,iglob)

          vn = nx*vx+nz*vz
          vnf = nx*vxf+nz*vzf

          tx = rho_vpII*vnf*nx - rho_vs*(vx-vn*nx)
          tz = rho_vpII*vnf*nz - rho_vs*(vz-vn*nz)

          accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - tx*weight
          accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - tz*weight

          ! saves contribution
          if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. NO_BACKWARD_RECONSTRUCTION)) then
            b_absorb_poro_w_bottom(1,i,ib_bottom(ispecabs),it) = tx*weight
            b_absorb_poro_w_bottom(2,i,ib_bottom(ispecabs),it) = tz*weight
          endif

        enddo

      endif  !  end of bottom absorbing boundary

      !--- top absorbing boundary
      if (codeabs(IEDGE3,ispecabs)) then
        j = NGLLZ

        ibegin = ibegin_edge3_poro(ispecabs)
        iend = iend_edge3_poro(ispecabs)

        ! exclude corners to make sure there is no contradiction on the normal
        if (codeabs_corner(3,ispecabs)) ibegin = 2
        if (codeabs_corner(4,ispecabs)) iend = NGLLX-1

        do i = ibegin,iend
          iglob = ibool(i,j,ispec)

          ! poroelastic material
          phi = phistore(i,j,ispec)
          tort = tortstore(i,j,ispec)

          rho_s = rhoarraystore(1,i,j,ispec)
          rho_f = rhoarraystore(2,i,j,ispec)

          cpIl = rho_vpstore(i,j,ispec) / rho_s
          cpIIl = vpIIstore(i,j,ispec)
          csl = rho_vsstore(i,j,ispec) / rho_s

          rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

          rho_vpI = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIl
          rho_vpII = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIIl
          rho_vs = rho_f/rho_bar*(rho_bar-rho_f*phi/tort)*csl

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = - zxi / jacobian1D
          nz = + xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

          ! solid contribution
          vx = velocs_poroelastic(1,iglob)
          vz = velocs_poroelastic(2,iglob)

          ! fluid contribution
          vxf = velocw_poroelastic(1,iglob)
          vzf = velocw_poroelastic(2,iglob)

          vn = nx*vx+nz*vz
          vnf = nx*vxf+nz*vzf

          tx = rho_vpII*vnf*nx - rho_vs*(vx-vn*nx)
          tz = rho_vpII*vnf*nz - rho_vs*(vz-vn*nz)

          accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - tx*weight
          accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - tz*weight

          ! saves contribution
          if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. NO_BACKWARD_RECONSTRUCTION)) then
            b_absorb_poro_w_top(1,i,ib_top(ispecabs),it) = tx*weight
            b_absorb_poro_w_top(2,i,ib_top(ispecabs),it) = tz*weight
          endif

        enddo

      endif  !  end of top absorbing boundary

    endif ! if ispec_is_poroelastic(ispec)

  enddo

  end subroutine compute_stacey_poro_fluid

!
!---------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_poro_fluid_backward(b_accelw_poroelastic)

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4,TWO,ZERO

  use specfem_par, only: it,NSTEP,STACEY_ABSORBING_CONDITIONS,NO_BACKWARD_RECONSTRUCTION, &
                         anyabs,num_abs_boundary_faces,abs_boundary_ispec, &
                         ibool,ispec_is_poroelastic,any_poroelastic, &
                         codeabs,codeabs_corner, &
                         nglob_poroelastic, &
                         ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro,iend_edge3_poro, &
                         ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro, &
                         SIMULATION_TYPE, &
                         b_absorb_poro_w_left,b_absorb_poro_w_right, &
                         b_absorb_poro_w_bottom,b_absorb_poro_w_top, &
                         ib_left,ib_right,ib_bottom,ib_top

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic),intent(inout) :: b_accelw_poroelastic

  ! local parameters
  integer :: ispec,i,j,iglob
  integer :: ispecabs,ibegin,iend,jbegin,jend
  integer :: it_tmp

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. anyabs) return
  if (.not. any_poroelastic) return

  if (NO_BACKWARD_RECONSTRUCTION) return
  if (SIMULATION_TYPE /= 3) return

  ! time increment index
  it_tmp = NSTEP - it + 1

  ! absorbing boundaries
  do ispecabs = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(ispecabs)

    if (ispec_is_poroelastic(ispec)) then

      !--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then
        i = 1

        jbegin = ibegin_edge4_poro(ispecabs)
        jend = iend_edge4_poro(ispecabs)

        do j = jbegin,jend
          iglob = ibool(i,j,ispec)
          b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) - &
                                          b_absorb_poro_w_left(1,j,ib_left(ispecabs),it_tmp)
          b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) - &
                                          b_absorb_poro_w_left(2,j,ib_left(ispecabs),it_tmp)
        enddo
      endif  !  end of left absorbing boundary

      !--- right absorbing boundary
      if (codeabs(IEDGE2,ispecabs)) then
        i = NGLLX

        jbegin = ibegin_edge2_poro(ispecabs)
        jend = iend_edge2_poro(ispecabs)

        do j = jbegin,jend
          iglob = ibool(i,j,ispec)
          b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) - &
                                          b_absorb_poro_w_right(1,j,ib_right(ispecabs),it_tmp)
          b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) - &
                                          b_absorb_poro_w_right(2,j,ib_right(ispecabs),it_tmp)
        enddo
      endif  !  end of right absorbing boundary

      !--- bottom absorbing boundary
      if (codeabs(IEDGE1,ispecabs)) then
        j = 1

        ibegin = ibegin_edge1_poro(ispecabs)
        iend = iend_edge1_poro(ispecabs)

        ! exclude corners to make sure there is no contradiction on the normal
        if (codeabs_corner(1,ispecabs)) ibegin = 2
        if (codeabs_corner(2,ispecabs)) iend = NGLLX-1

        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) - &
                                          b_absorb_poro_w_bottom(1,i,ib_bottom(ispecabs),it_tmp)
          b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) - &
                                          b_absorb_poro_w_bottom(2,i,ib_bottom(ispecabs),it_tmp)
        enddo
      endif  !  end of bottom absorbing boundary

      !--- top absorbing boundary
      if (codeabs(IEDGE3,ispecabs)) then
        j = NGLLZ

        ibegin = ibegin_edge3_poro(ispecabs)
        iend = iend_edge3_poro(ispecabs)

        ! exclude corners to make sure there is no contradiction on the normal
        if (codeabs_corner(3,ispecabs)) ibegin = 2
        if (codeabs_corner(4,ispecabs)) iend = NGLLX-1

        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) - &
                                          b_absorb_poro_w_top(1,i,ib_top(ispecabs),it_tmp)
          b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) - &
                                          b_absorb_poro_w_top(2,i,ib_top(ispecabs),it_tmp)
        enddo
      endif  !  end of top absorbing boundary

    endif ! if ispec_is_poroelastic(ispec)

  enddo

  end subroutine compute_stacey_poro_fluid_backward

!---------------------------------------------------------------------------------------------
!
! solid part
!
!---------------------------------------------------------------------------------------------

  subroutine compute_stacey_poro_solid(accels_poroelastic,velocs_poroelastic,velocw_poroelastic)

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4,TWO,ZERO

  use specfem_par, only: it,STACEY_ABSORBING_CONDITIONS, &
                         anyabs,num_abs_boundary_faces,abs_boundary_ispec, &
                         ibool,ispec_is_poroelastic,any_poroelastic, &
                         codeabs,codeabs_corner, &
                         nglob_poroelastic, &
                         phistore,tortstore,rhoarraystore,vpIIstore,rho_vpstore,rho_vsstore, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         wxgll,wzgll, &
                         ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro,iend_edge3_poro, &
                         ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro, &
                         SIMULATION_TYPE,SAVE_FORWARD, &
                         b_absorb_poro_s_left,b_absorb_poro_s_right, &
                         b_absorb_poro_s_bottom,b_absorb_poro_s_top, &
                         ib_left,ib_right,ib_bottom,ib_top, &
                         NO_BACKWARD_RECONSTRUCTION

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic),intent(inout) :: accels_poroelastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic),intent(in) :: velocs_poroelastic,velocw_poroelastic

  ! local parameters
  integer :: ispec,i,j,iglob
  integer :: ispecabs,ibegin,iend,jbegin,jend

  real(kind=CUSTOM_REAL) :: nx,nz,vx,vz,vn,vxf,vzf,vnf
  real(kind=CUSTOM_REAL) :: rho_vpI,rho_vpII,rho_vs
  real(kind=CUSTOM_REAL) :: tx,tz,weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: cpIl,cpIIl,csl
  real(kind=CUSTOM_REAL) :: phi,tort,rho_s,rho_f,rho_bar

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. anyabs) return
  if (.not. any_poroelastic) return

  ! absorbing boundaries
  do ispecabs = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(ispecabs)

    if (ispec_is_poroelastic(ispec)) then

      !--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then

        i = 1

        jbegin = ibegin_edge4_poro(ispecabs)
        jend = iend_edge4_poro(ispecabs)

        do j = jbegin,jend
          iglob = ibool(i,j,ispec)

          ! poroelastic material
          phi = phistore(i,j,ispec)
          tort = tortstore(i,j,ispec)

          rho_s = rhoarraystore(1,i,j,ispec)
          rho_f = rhoarraystore(2,i,j,ispec)

          cpIl = rho_vpstore(i,j,ispec) / rho_s
          cpIIl = vpIIstore(i,j,ispec)
          csl = rho_vsstore(i,j,ispec) / rho_s

          rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

          rho_vpI = (rho_bar - phi/tort*rho_f)*cpIl
          rho_vpII = (rho_bar - phi/tort*rho_f)*cpIIl
          rho_vs = (rho_bar - phi/tort*rho_f)*csl

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = - zgamma / jacobian1D
          nz = + xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

          ! solid contribution
          vx = velocs_poroelastic(1,iglob)
          vz = velocs_poroelastic(2,iglob)

          ! fluid contribution
          vxf = velocw_poroelastic(1,iglob)
          vzf = velocw_poroelastic(2,iglob)

          vn = nx*vx+nz*vz
          vnf = nx*vxf+nz*vzf

          tx = rho_vpI*vn*nx + rho_vs*(vx-vn*nx)
          tz = rho_vpI*vn*nz + rho_vs*(vz-vn*nz)

          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - tx*weight
          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) - tz*weight

          ! saves contribution for reconstructed/backward wavefield
          if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. NO_BACKWARD_RECONSTRUCTION)) then
            b_absorb_poro_s_left(1,j,ib_left(ispecabs),it) = tx*weight
            b_absorb_poro_s_left(2,j,ib_left(ispecabs),it) = tz*weight
          endif

        enddo

      endif  !  end of left absorbing boundary

      !--- right absorbing boundary
      if (codeabs(IEDGE2,ispecabs)) then

        i = NGLLX

        jbegin = ibegin_edge2_poro(ispecabs)
        jend = iend_edge2_poro(ispecabs)

        do j = jbegin,jend
          iglob = ibool(i,j,ispec)

          ! poroelastic material
          phi = phistore(i,j,ispec)
          tort = tortstore(i,j,ispec)

          rho_s = rhoarraystore(1,i,j,ispec)
          rho_f = rhoarraystore(2,i,j,ispec)

          cpIl = rho_vpstore(i,j,ispec) / rho_s
          cpIIl = vpIIstore(i,j,ispec)
          csl = rho_vsstore(i,j,ispec) / rho_s

          rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

          rho_vpI = (rho_bar - phi/tort*rho_f)*cpIl
          rho_vpII = (rho_bar - phi/tort*rho_f)*cpIIl
          rho_vs = (rho_bar - phi/tort*rho_f)*csl

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = + zgamma / jacobian1D
          nz = - xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

          ! solid contribution
          vx = velocs_poroelastic(1,iglob)
          vz = velocs_poroelastic(2,iglob)

          ! fluid contribution
          vxf = velocw_poroelastic(1,iglob)
          vzf = velocw_poroelastic(2,iglob)

          vn = nx*vx+nz*vz
          vnf = nx*vxf+nz*vzf

          tx = rho_vpI*vn*nx + rho_vs*(vx-vn*nx)
          tz = rho_vpI*vn*nz + rho_vs*(vz-vn*nz)

          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - tx*weight
          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) - tz*weight

          ! saves contribution for reconstructed/backward wavefield
          if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. NO_BACKWARD_RECONSTRUCTION)) then
            b_absorb_poro_s_right(1,j,ib_right(ispecabs),it) = tx*weight
            b_absorb_poro_s_right(2,j,ib_right(ispecabs),it) = tz*weight
          endif

        enddo

      endif  !  end of right absorbing boundary

      !--- bottom absorbing boundary
      if (codeabs(IEDGE1,ispecabs)) then

        j = 1

        ibegin = ibegin_edge1_poro(ispecabs)
        iend = iend_edge1_poro(ispecabs)

        ! exclude corners to make sure there is no contradiction on the normal
        if (codeabs_corner(1,ispecabs)) ibegin = 2
        if (codeabs_corner(2,ispecabs)) iend = NGLLX-1

        do i = ibegin,iend
          iglob = ibool(i,j,ispec)

          ! poroelastic material
          phi = phistore(i,j,ispec)
          tort = tortstore(i,j,ispec)

          rho_s = rhoarraystore(1,i,j,ispec)
          rho_f = rhoarraystore(2,i,j,ispec)

          cpIl = rho_vpstore(i,j,ispec) / rho_s
          cpIIl = vpIIstore(i,j,ispec)
          csl = rho_vsstore(i,j,ispec) / rho_s

          rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

          rho_vpI = (rho_bar - phi/tort*rho_f)*cpIl
          rho_vpII = (rho_bar - phi/tort*rho_f)*cpIIl
          rho_vs = (rho_bar - phi/tort*rho_f)*csl

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = + zxi / jacobian1D
          nz = - xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

          ! solid contribution
          vx = velocs_poroelastic(1,iglob)
          vz = velocs_poroelastic(2,iglob)

          ! fluid contribution
          vxf = velocw_poroelastic(1,iglob)
          vzf = velocw_poroelastic(2,iglob)

          vn = nx*vx+nz*vz
          vnf = nx*vxf+nz*vzf

          tx = rho_vpI*vn*nx + rho_vs*(vx-vn*nx)
          tz = rho_vpI*vn*nz + rho_vs*(vz-vn*nz)

          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - tx*weight
          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) - tz*weight

          ! saves contribution for reconstructed/backward wavefield
          if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. NO_BACKWARD_RECONSTRUCTION)) then
            b_absorb_poro_s_bottom(1,i,ib_bottom(ispecabs),it) = tx*weight
            b_absorb_poro_s_bottom(2,i,ib_bottom(ispecabs),it) = tz*weight
          endif

        enddo

      endif  !  end of bottom absorbing boundary

      !--- top absorbing boundary
      if (codeabs(IEDGE3,ispecabs)) then

        j = NGLLZ

        ibegin = ibegin_edge3_poro(ispecabs)
        iend = iend_edge3_poro(ispecabs)

        ! exclude corners to make sure there is no contradiction on the normal
        if (codeabs_corner(3,ispecabs)) ibegin = 2
        if (codeabs_corner(4,ispecabs)) iend = NGLLX-1

        do i = ibegin,iend
          iglob = ibool(i,j,ispec)

          ! poroelastic material
          phi = phistore(i,j,ispec)
          tort = tortstore(i,j,ispec)

          rho_s = rhoarraystore(1,i,j,ispec)
          rho_f = rhoarraystore(2,i,j,ispec)

          cpIl = rho_vpstore(i,j,ispec) / rho_s
          cpIIl = vpIIstore(i,j,ispec)
          csl = rho_vsstore(i,j,ispec) / rho_s

          rho_bar = (1.d0 - phi)*rho_s + phi*rho_f

          rho_vpI = (rho_bar - phi/tort*rho_f)*cpIl
          rho_vpII = (rho_bar - phi/tort*rho_f)*cpIIl
          rho_vs = (rho_bar - phi/tort*rho_f)*csl

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = - zxi / jacobian1D
          nz = + xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

          ! solid contribution
          vx = velocs_poroelastic(1,iglob)
          vz = velocs_poroelastic(2,iglob)

          ! fluid contribution
          vxf = velocw_poroelastic(1,iglob)
          vzf = velocw_poroelastic(2,iglob)

          vn = nx*vx+nz*vz
          vnf = nx*vxf+nz*vzf

          tx = rho_vpI*vn*nx + rho_vs*(vx-vn*nx)
          tz = rho_vpI*vn*nz + rho_vs*(vz-vn*nz)

          accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - tx*weight
          accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) - tz*weight

          ! saves contribution for reconstructed/backward wavefield
          if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. NO_BACKWARD_RECONSTRUCTION)) then
            b_absorb_poro_s_top(1,i,ib_top(ispecabs),it) = tx*weight
            b_absorb_poro_s_top(2,i,ib_top(ispecabs),it) = tz*weight
          endif

        enddo

      endif  !  end of top absorbing boundary

    endif ! if ispec_is_poroelastic(ispec)

  enddo

  end subroutine compute_stacey_poro_solid

!
!---------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_poro_solid_backward(b_accels_poroelastic)

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4,TWO,ZERO

  use specfem_par, only: it,NSTEP,STACEY_ABSORBING_CONDITIONS,NO_BACKWARD_RECONSTRUCTION, &
                         anyabs,num_abs_boundary_faces,abs_boundary_ispec, &
                         ibool,ispec_is_poroelastic,any_poroelastic, &
                         codeabs,codeabs_corner, &
                         nglob_poroelastic, &
                         ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro,iend_edge3_poro, &
                         ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro, &
                         SIMULATION_TYPE, &
                         b_absorb_poro_s_left,b_absorb_poro_s_right, &
                         b_absorb_poro_s_bottom,b_absorb_poro_s_top, &
                         ib_left,ib_right,ib_bottom,ib_top

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic),intent(inout) :: b_accels_poroelastic

  ! local parameters
  integer :: ispec,i,j,iglob
  integer :: ispecabs,ibegin,iend,jbegin,jend
  integer :: it_tmp

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. anyabs) return
  if (.not. any_poroelastic) return

  if (NO_BACKWARD_RECONSTRUCTION) return
  if (SIMULATION_TYPE /= 3) return

  ! time increment index
  it_tmp = NSTEP - it + 1

  ! absorbing boundaries
  do ispecabs = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(ispecabs)

    if (ispec_is_poroelastic(ispec)) then

      !--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then
        i = 1

        jbegin = ibegin_edge4_poro(ispecabs)
        jend = iend_edge4_poro(ispecabs)

        do j = jbegin,jend
          iglob = ibool(i,j,ispec)
          ! reconstructed/backward wavefield
          b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - &
                                          b_absorb_poro_s_left(1,j,ib_left(ispecabs),it_tmp)
          b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - &
                                          b_absorb_poro_s_left(2,j,ib_left(ispecabs),it_tmp)
        enddo
      endif  !  end of left absorbing boundary

      !--- right absorbing boundary
      if (codeabs(IEDGE2,ispecabs)) then
        i = NGLLX

        jbegin = ibegin_edge2_poro(ispecabs)
        jend = iend_edge2_poro(ispecabs)

        do j = jbegin,jend
          iglob = ibool(i,j,ispec)
          ! reconstructed/backward wavefield
          b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - &
                                          b_absorb_poro_s_right(1,j,ib_right(ispecabs),it_tmp)
          b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - &
                                          b_absorb_poro_s_right(2,j,ib_right(ispecabs),it_tmp)
        enddo
      endif  !  end of right absorbing boundary

      !--- bottom absorbing boundary
      if (codeabs(IEDGE1,ispecabs)) then
        j = 1

        ibegin = ibegin_edge1_poro(ispecabs)
        iend = iend_edge1_poro(ispecabs)

        ! exclude corners to make sure there is no contradiction on the normal
        if (codeabs_corner(1,ispecabs)) ibegin = 2
        if (codeabs_corner(2,ispecabs)) iend = NGLLX-1

        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          ! reconstructed/backward wavefield
          b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - &
                                          b_absorb_poro_s_bottom(1,i,ib_bottom(ispecabs),it_tmp)
          b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - &
                                          b_absorb_poro_s_bottom(2,i,ib_bottom(ispecabs),it_tmp)
        enddo
      endif  !  end of bottom absorbing boundary

      !--- top absorbing boundary
      if (codeabs(IEDGE3,ispecabs)) then
        j = NGLLZ

        ibegin = ibegin_edge3_poro(ispecabs)
        iend = iend_edge3_poro(ispecabs)

        ! exclude corners to make sure there is no contradiction on the normal
        if (codeabs_corner(3,ispecabs)) ibegin = 2
        if (codeabs_corner(4,ispecabs)) iend = NGLLX-1

        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          ! reconstructed/backward wavefield
          b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - &
                                          b_absorb_poro_s_top(1,i,ib_top(ispecabs),it_tmp)
          b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - &
                                          b_absorb_poro_s_top(2,i,ib_top(ispecabs),it_tmp)
        enddo
      endif  !  end of top absorbing boundary

    endif ! if ispec_is_poroelastic(ispec)

  enddo

  end subroutine compute_stacey_poro_solid_backward


