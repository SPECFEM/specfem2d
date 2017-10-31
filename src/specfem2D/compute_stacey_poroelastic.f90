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

  subroutine compute_stacey_poro_fluid(f0)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4,TWO,ZERO

  use specfem_par, only: it,STACEY_ABSORBING_CONDITIONS, &
                         anyabs,nelemabs,numabs, &
                         ATTENUATION_PORO_FLUID_PART, &
                         ibool,kmato,ispec_is_poroelastic, &
                         codeabs,codeabs_corner, &
                         accelw_poroelastic, &
                         velocw_poroelastic,velocs_poroelastic, &
                         permeability,xix,xiz,gammax,gammaz, &
                         jacobian, &
                         wxgll,wzgll, &
                         ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro,iend_edge3_poro, &
                         ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro, &
                         SIMULATION_TYPE,SAVE_FORWARD, &
                         b_absorb_poro_w_left,b_absorb_poro_w_right, &
                         b_absorb_poro_w_bottom,b_absorb_poro_w_top, &
                         ib_left,ib_right,ib_bottom,ib_top, &
                         freq0_poroelastic,Q0_poroelastic

  implicit none

  double precision,intent(in) :: f0

  ! local parameters
  double precision :: w_c

  integer :: ispec,i,j,iglob
  integer :: ispecabs,ibegin,iend,jbegin,jend

  real(kind=CUSTOM_REAL) :: nx,nz,vx,vz,vn,vxf,vzf,vnf
  real(kind=CUSTOM_REAL) :: rho_vpI,rho_vpII,rho_vs
  real(kind=CUSTOM_REAL) :: tx,tz,weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: cpIl,cpIIl,csl

  double precision :: cpIsquare,cpIIsquare,cssquare
  double precision :: phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot
  double precision :: permlxx

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. anyabs) return

  ! absorbing boundaries
  do ispecabs= 1,nelemabs

    ispec = numabs(ispecabs)

    if (ispec_is_poroelastic(ispec)) then

      ! get poroelastic parameters of current spectral element
      call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)

      ! Biot coefficients for the input phi
      call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

      permlxx = permeability(1,kmato(ispec))

      call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mu_fr,phi, &
                   tort,rho_s,rho_f,eta_f,permlxx,f0,freq0_poroelastic,Q0_poroelastic,w_c,ATTENUATION_PORO_FLUID_PART)

      cpIl = sqrt(cpIsquare)
      cpIIl = sqrt(cpIIsquare)
      csl = sqrt(cssquare)

!--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then

        i = 1

        jbegin = ibegin_edge4_poro(ispecabs)
        jend = iend_edge4_poro(ispecabs)

        do j = jbegin,jend
          iglob = ibool(i,j,ispec)

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = - zgamma / jacobian1D
          nz = + xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

          rho_vpI = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIl
          rho_vpII = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIIl
          rho_vs = rho_f/rho_bar*(rho_bar-rho_f*phi/tort)*csl

          if (ispec_is_poroelastic(ispec)) then
            vx = velocs_poroelastic(1,iglob)
            vz = velocs_poroelastic(2,iglob)
            vxf = velocw_poroelastic(1,iglob)
            vzf = velocw_poroelastic(2,iglob)

            vn = nx*vx+nz*vz
            vnf = nx*vxf+nz*vzf

            tx = rho_vpII*vnf*nx - rho_vs*(vx-vn*nx)
            tz = rho_vpII*vnf*nz - rho_vs*(vz-vn*nz)

            accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - tx*weight
            accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - tz*weight

            if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
              b_absorb_poro_w_left(1,j,ib_left(ispecabs),it) = tx*weight
              b_absorb_poro_w_left(2,j,ib_left(ispecabs),it) = tz*weight
            endif

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

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = + zgamma / jacobian1D
          nz = - xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

          rho_vpI = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIl
          rho_vpII = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIIl
          rho_vs = rho_f/rho_bar*(rho_bar-rho_f*phi/tort)*csl

          if (ispec_is_poroelastic(ispec)) then
            vx = velocs_poroelastic(1,iglob)
            vz = velocs_poroelastic(2,iglob)
            vxf = velocw_poroelastic(1,iglob)
            vzf = velocw_poroelastic(2,iglob)

            vn = nx*vx+nz*vz
            vnf = nx*vxf+nz*vzf

            tx = rho_vpII*vnf*nx - rho_vs*(vx-vn*nx)
            tz = rho_vpII*vnf*nz - rho_vs*(vz-vn*nz)

            accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - tx*weight
            accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - tz*weight

            if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
              b_absorb_poro_w_right(1,j,ib_right(ispecabs),it) = tx*weight
              b_absorb_poro_w_right(2,j,ib_right(ispecabs),it) = tz*weight
            endif

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

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = + zxi / jacobian1D
          nz = - xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

          rho_vpI = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIl
          rho_vpII = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIIl
          rho_vs = rho_f/rho_bar*(rho_bar-rho_f*phi/tort)*csl

          if (ispec_is_poroelastic(ispec)) then
            vx = velocs_poroelastic(1,iglob)
            vz = velocs_poroelastic(2,iglob)
            vxf = velocw_poroelastic(1,iglob)
            vzf = velocw_poroelastic(2,iglob)

            vn = nx*vx+nz*vz
            vnf = nx*vxf+nz*vzf

            tx = rho_vpII*vnf*nx - rho_vs*(vx-vn*nx)
            tz = rho_vpII*vnf*nz - rho_vs*(vz-vn*nz)

            accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - tx*weight
            accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - tz*weight

            if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
              b_absorb_poro_w_bottom(1,i,ib_bottom(ispecabs),it) = tx*weight
              b_absorb_poro_w_bottom(2,i,ib_bottom(ispecabs),it) = tz*weight
            endif

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

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = - zxi / jacobian1D
          nz = + xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

          rho_vpI = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIl
          rho_vpII = (rho_f*tort*rho_bar - phi*rho_f*rho_f)/(phi*rho_bar)*cpIIl
          rho_vs = rho_f/rho_bar*(rho_bar-rho_f*phi/tort)*csl

          if (ispec_is_poroelastic(ispec)) then
            vx = velocs_poroelastic(1,iglob)
            vz = velocs_poroelastic(2,iglob)
            vxf = velocw_poroelastic(1,iglob)
            vzf = velocw_poroelastic(2,iglob)

            vn = nx*vx+nz*vz
            vnf = nx*vxf+nz*vzf

            tx = rho_vpII*vnf*nx - rho_vs*(vx-vn*nx)
            tz = rho_vpII*vnf*nz - rho_vs*(vz-vn*nz)

            accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - tx*weight
            accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) - tz*weight

            if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
              b_absorb_poro_w_top(1,i,ib_top(ispecabs),it) = tx*weight
              b_absorb_poro_w_top(2,i,ib_top(ispecabs),it) = tz*weight
            endif

          endif

        enddo

      endif  !  end of top absorbing boundary

    endif ! if ispec_is_poroelastic(ispec)

  enddo

  end subroutine compute_stacey_poro_fluid

!
!---------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_poro_fluid_backward()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4,TWO,ZERO

  use specfem_par, only: it,NSTEP,STACEY_ABSORBING_CONDITIONS, &
                         anyabs,nelemabs,numabs, &
                         ibool,ispec_is_poroelastic, &
                         codeabs,codeabs_corner, &
                         b_accelw_poroelastic, &
                         ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro,iend_edge3_poro, &
                         ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro, &
                         SIMULATION_TYPE, &
                         b_absorb_poro_w_left,b_absorb_poro_w_right, &
                         b_absorb_poro_w_bottom,b_absorb_poro_w_top, &
                         ib_left,ib_right,ib_bottom,ib_top

  implicit none

  ! local parameters
  integer :: ispec,i,j,iglob
  integer :: ispecabs,ibegin,iend,jbegin,jend

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. anyabs) return
  if (SIMULATION_TYPE /= 3) return

  ! absorbing boundaries
  do ispecabs= 1,nelemabs

    ispec = numabs(ispecabs)

    if (ispec_is_poroelastic(ispec)) then

!--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then
        i = 1

        jbegin = ibegin_edge4_poro(ispecabs)
        jend = iend_edge4_poro(ispecabs)

        do j = jbegin,jend
          iglob = ibool(i,j,ispec)

          if (ispec_is_poroelastic(ispec)) then
            b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) - &
                                            b_absorb_poro_w_left(1,j,ib_left(ispecabs),NSTEP-it+1)
            b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) - &
                                            b_absorb_poro_w_left(2,j,ib_left(ispecabs),NSTEP-it+1)
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

          if (ispec_is_poroelastic(ispec)) then
            b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) - &
                                            b_absorb_poro_w_right(1,j,ib_right(ispecabs),NSTEP-it+1)
            b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) - &
                                            b_absorb_poro_w_right(2,j,ib_right(ispecabs),NSTEP-it+1)
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

          if (ispec_is_poroelastic(ispec)) then
            b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) - &
                                            b_absorb_poro_w_bottom(1,i,ib_bottom(ispecabs),NSTEP-it+1)
            b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) - &
                                            b_absorb_poro_w_bottom(2,i,ib_bottom(ispecabs),NSTEP-it+1)
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

          if (ispec_is_poroelastic(ispec)) then
            b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) - &
                                            b_absorb_poro_w_top(1,i,ib_top(ispecabs),NSTEP-it+1)
            b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) - &
                                            b_absorb_poro_w_top(2,i,ib_top(ispecabs),NSTEP-it+1)
          endif
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

  subroutine compute_stacey_poro_solid(f0)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4,TWO,ZERO

  use specfem_par, only: it,NSTEP,STACEY_ABSORBING_CONDITIONS, &
                         anyabs,nelemabs,numabs, &
                         ATTENUATION_PORO_FLUID_PART, &
                         ibool,kmato,ispec_is_poroelastic, &
                         codeabs,codeabs_corner, &
                         accels_poroelastic,b_accels_poroelastic, &
                         velocw_poroelastic,velocs_poroelastic, &
                         permeability,xix,xiz,gammax,gammaz, &
                         jacobian, &
                         wxgll,wzgll, &
                         ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro,iend_edge3_poro, &
                         ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro, &
                         SIMULATION_TYPE,SAVE_FORWARD, &
                         b_absorb_poro_s_left,b_absorb_poro_s_right, &
                         b_absorb_poro_s_bottom,b_absorb_poro_s_top, &
                         ib_left,ib_right,ib_bottom,ib_top, &
                         freq0_poroelastic,Q0_poroelastic

  implicit none

  double precision,intent(in) :: f0

  ! local parameters
  double precision :: w_c

  integer :: ispec,i,j,iglob
  integer :: ispecabs,ibegin,iend,jbegin,jend

  real(kind=CUSTOM_REAL) :: nx,nz,vx,vz,vn,vxf,vzf,vnf
  real(kind=CUSTOM_REAL) :: rho_vpI,rho_vpII,rho_vs
  real(kind=CUSTOM_REAL) :: tx,tz,weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: cpIl,cpIIl,csl

  double precision :: cpIsquare,cpIIsquare,cssquare
  double precision :: phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot
  double precision :: permlxx

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. anyabs) return

  ! absorbing boundaries
  do ispecabs = 1,nelemabs

    ispec = numabs(ispecabs)

    if (ispec_is_poroelastic(ispec)) then

      ! get poroelastic parameters of current spectral element
      call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)

      ! Biot coefficients for the input phi
      call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

      permlxx = permeability(1,kmato(ispec))

      call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mu_fr,phi, &
                    tort,rho_s,rho_f,eta_f,permlxx,f0,freq0_poroelastic,Q0_poroelastic,w_c,ATTENUATION_PORO_FLUID_PART)

      cpIl = sqrt(cpIsquare)
      cpIIl = sqrt(cpIIsquare)
      csl = sqrt(cssquare)

!--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then

        i = 1

        jbegin = ibegin_edge4_poro(ispecabs)
        jend = iend_edge4_poro(ispecabs)

        do j = jbegin,jend
          iglob = ibool(i,j,ispec)

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = - zgamma / jacobian1D
          nz = + xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

          rho_vpI = (rho_bar - phi/tort*rho_f)*cpIl
          rho_vpII = (rho_bar - phi/tort*rho_f)*cpIIl
          rho_vs = (rho_bar - phi/tort*rho_f)*csl

          if (ispec_is_poroelastic(ispec)) then
            vx = velocs_poroelastic(1,iglob)
            vz = velocs_poroelastic(2,iglob)

            vxf = velocw_poroelastic(1,iglob)
            vzf = velocw_poroelastic(2,iglob)

            vn = nx*vx+nz*vz
            vnf = nx*vxf+nz*vzf

            tx = rho_vpI*vn*nx + rho_vs*(vx-vn*nx)
            tz = rho_vpI*vn*nz + rho_vs*(vz-vn*nz)

            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - tx*weight
            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) - tz*weight

            ! reconstructed/backward wavefield
            if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
              b_absorb_poro_s_left(1,j,ib_left(ispecabs),it) = tx*weight
              b_absorb_poro_s_left(2,j,ib_left(ispecabs),it) = tz*weight
            else if (SIMULATION_TYPE == 3) then
              b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - &
                                              b_absorb_poro_s_left(1,j,ib_left(ispecabs),NSTEP-it+1)
              b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - &
                                              b_absorb_poro_s_left(2,j,ib_left(ispecabs),NSTEP-it+1)
            endif

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

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = + zgamma / jacobian1D
          nz = - xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

          rho_vpI = (rho_bar - phi/tort*rho_f)*cpIl
          rho_vpII = (rho_bar - phi/tort*rho_f)*cpIIl
          rho_vs = (rho_bar - phi/tort*rho_f)*csl

          if (ispec_is_poroelastic(ispec)) then
            vx = velocs_poroelastic(1,iglob)
            vz = velocs_poroelastic(2,iglob)
            vxf = velocw_poroelastic(1,iglob)
            vzf = velocw_poroelastic(2,iglob)

            vn = nx*vx+nz*vz
            vnf = nx*vxf+nz*vzf

            tx = rho_vpI*vn*nx + rho_vs*(vx-vn*nx)
            tz = rho_vpI*vn*nz + rho_vs*(vz-vn*nz)

            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - tx*weight
            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) - tz*weight

            ! reconstructed/backward wavefield
            if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
              b_absorb_poro_s_right(1,j,ib_right(ispecabs),it) = tx*weight
              b_absorb_poro_s_right(2,j,ib_right(ispecabs),it) = tz*weight
            else if (SIMULATION_TYPE == 3) then
              b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - &
                                              b_absorb_poro_s_right(1,j,ib_right(ispecabs),NSTEP-it+1)
              b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - &
                                              b_absorb_poro_s_right(2,j,ib_right(ispecabs),NSTEP-it+1)
            endif

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

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = + zxi / jacobian1D
          nz = - xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

          rho_vpI = (rho_bar - phi/tort*rho_f)*cpIl
          rho_vpII = (rho_bar - phi/tort*rho_f)*cpIIl
          rho_vs = (rho_bar - phi/tort*rho_f)*csl

          if (ispec_is_poroelastic(ispec)) then
            vx = velocs_poroelastic(1,iglob)
            vz = velocs_poroelastic(2,iglob)
            vxf = velocw_poroelastic(1,iglob)
            vzf = velocw_poroelastic(2,iglob)

            vn = nx*vx+nz*vz
            vnf = nx*vxf+nz*vzf

            tx = rho_vpI*vn*nx + rho_vs*(vx-vn*nx)
            tz = rho_vpI*vn*nz + rho_vs*(vz-vn*nz)

            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - tx*weight
            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) - tz*weight

            ! reconstructed/backward wavefield
            if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
              b_absorb_poro_s_bottom(1,i,ib_bottom(ispecabs),it) = tx*weight
              b_absorb_poro_s_bottom(2,i,ib_bottom(ispecabs),it) = tz*weight
            else if (SIMULATION_TYPE == 3) then
              b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - &
                                              b_absorb_poro_s_bottom(1,i,ib_bottom(ispecabs),NSTEP-it+1)
              b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - &
                                              b_absorb_poro_s_bottom(2,i,ib_bottom(ispecabs),NSTEP-it+1)
            endif

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

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = - zxi / jacobian1D
          nz = + xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

          rho_vpI = (rho_bar - phi/tort*rho_f)*cpIl
          rho_vpII = (rho_bar - phi/tort*rho_f)*cpIIl
          rho_vs = (rho_bar - phi/tort*rho_f)*csl

          if (ispec_is_poroelastic(ispec)) then
            vx = velocs_poroelastic(1,iglob)
            vz = velocs_poroelastic(2,iglob)
            vxf = velocw_poroelastic(1,iglob)
            vzf = velocw_poroelastic(2,iglob)

            vn = nx*vx+nz*vz
            vnf = nx*vxf+nz*vzf

            tx = rho_vpI*vn*nx + rho_vs*(vx-vn*nx)
            tz = rho_vpI*vn*nz + rho_vs*(vz-vn*nz)

            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - tx*weight
            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) - tz*weight

            ! reconstructed/backward wavefield
            if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
              b_absorb_poro_s_top(1,i,ib_top(ispecabs),it) = tx*weight
              b_absorb_poro_s_top(2,i,ib_top(ispecabs),it) = tz*weight
            else if (SIMULATION_TYPE == 3) then
              b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - &
                                              b_absorb_poro_s_top(1,i,ib_top(ispecabs),NSTEP-it+1)
              b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - &
                                              b_absorb_poro_s_top(2,i,ib_top(ispecabs),NSTEP-it+1)
            endif

          endif

        enddo

      endif  !  end of top absorbing boundary

    endif ! if ispec_is_poroelastic(ispec)

  enddo

  end subroutine compute_stacey_poro_solid

!
!---------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_poro_solid_backward()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4,TWO,ZERO

  use specfem_par, only: it,NSTEP,STACEY_ABSORBING_CONDITIONS, &
                         anyabs,nelemabs,numabs, &
                         ibool,ispec_is_poroelastic, &
                         codeabs,codeabs_corner, &
                         b_accels_poroelastic, &
                         ibegin_edge1_poro,iend_edge1_poro,ibegin_edge3_poro,iend_edge3_poro, &
                         ibegin_edge4_poro,iend_edge4_poro,ibegin_edge2_poro,iend_edge2_poro, &
                         SIMULATION_TYPE, &
                         b_absorb_poro_s_left,b_absorb_poro_s_right, &
                         b_absorb_poro_s_bottom,b_absorb_poro_s_top, &
                         ib_left,ib_right,ib_bottom,ib_top

  implicit none

  ! local parameters
  integer :: ispec,i,j,iglob
  integer :: ispecabs,ibegin,iend,jbegin,jend

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. anyabs) return
  if (SIMULATION_TYPE /= 3) return

  ! absorbing boundaries
  do ispecabs = 1,nelemabs

    ispec = numabs(ispecabs)

    if (ispec_is_poroelastic(ispec)) then

!--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then
        i = 1

        jbegin = ibegin_edge4_poro(ispecabs)
        jend = iend_edge4_poro(ispecabs)

        do j = jbegin,jend
          iglob = ibool(i,j,ispec)

          if (ispec_is_poroelastic(ispec)) then
            ! reconstructed/backward wavefield
            b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - &
                                            b_absorb_poro_s_left(1,j,ib_left(ispecabs),NSTEP-it+1)
            b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - &
                                            b_absorb_poro_s_left(2,j,ib_left(ispecabs),NSTEP-it+1)
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

          if (ispec_is_poroelastic(ispec)) then
            ! reconstructed/backward wavefield
            b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - &
                                            b_absorb_poro_s_right(1,j,ib_right(ispecabs),NSTEP-it+1)
            b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - &
                                            b_absorb_poro_s_right(2,j,ib_right(ispecabs),NSTEP-it+1)
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

          if (ispec_is_poroelastic(ispec)) then
            ! reconstructed/backward wavefield
            b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - &
                                            b_absorb_poro_s_bottom(1,i,ib_bottom(ispecabs),NSTEP-it+1)
            b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - &
                                            b_absorb_poro_s_bottom(2,i,ib_bottom(ispecabs),NSTEP-it+1)
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

          if (ispec_is_poroelastic(ispec)) then
            ! reconstructed/backward wavefield
            b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) - &
                                            b_absorb_poro_s_top(1,i,ib_top(ispecabs),NSTEP-it+1)
            b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) - &
                                            b_absorb_poro_s_top(2,i,ib_top(ispecabs),NSTEP-it+1)
          endif
        enddo
      endif  !  end of top absorbing boundary

    endif ! if ispec_is_poroelastic(ispec)

  enddo

  end subroutine compute_stacey_poro_solid_backward


