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

  subroutine compute_stacey_elastic(accel_elastic,veloc_elastic)

! absorbing boundaries
!
! Clayton-Engquist condition if elastic

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM, &
    ZERO,ONE,TWO,TWO_THIRDS,FOUR_THIRDS,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: nglob,num_abs_boundary_faces,anyabs,it,any_elastic, &
                         ibool, &
                         abs_boundary_ispec,ispec_is_elastic, &
                         codeabs,codeabs_corner, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         rho_vpstore,rho_vsstore,rhostore,mustore, &
                         wxgll,wzgll,P_SV, &
                         SIMULATION_TYPE,SAVE_FORWARD, &
                         b_absorb_elastic_left,b_absorb_elastic_right, &
                         b_absorb_elastic_bottom,b_absorb_elastic_top, &
                         ib_left,ib_right,ib_bottom,ib_top, &
                         STACEY_ABSORBING_CONDITIONS, &
                         NO_BACKWARD_RECONSTRUCTION

  ! initialfield
  use specfem_par, only: v0x_left,v0z_left,v0x_right,v0z_right,v0x_bot,v0z_bot, &
                        t0x_left,t0z_left,t0x_right,t0z_right,t0x_bot,t0z_bot, &
                        add_Bielak_conditions_bottom,add_Bielak_conditions_right, &
                        add_Bielak_conditions_top,add_Bielak_conditions_left, &
                        initialfield,over_critical_angle, &
                        anglesource,anglesource_refl,A_plane,B_plane,C_plane,c_inc,c_refl,time_offset

  ! for Bielak
  use specfem_par, only: x_source,z_source,f0_source

  use specfem_par, only: displ_elastic,hprime_xx,hprime_zz

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: accel_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: veloc_elastic

  ! local parameters
  integer :: ispecabs,ispec,i,j,k,iglob,iiglob
  real(kind=CUSTOM_REAL) :: weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: nx,nz,vx,vy,vz,vn,tx,ty,tz
  real(kind=CUSTOM_REAL) :: dn,dx_n,dz_n,dx_t,dz_t,tx_n,tz_n,tx_t,tz_t
  real(kind=CUSTOM_REAL) :: unit_tx,unit_tz,v_norm
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl
  real(kind=CUSTOM_REAL) :: tempx1l_n,tempz1l_n,tempx2l_n,tempz2l_n
  real(kind=CUSTOM_REAL) :: tempx1l_t,tempz1l_t,tempx2l_t,tempz2l_t
  real(kind=CUSTOM_REAL) :: dux_dxl,dux_dzl,duz_dxl,duz_dzl,deriv_x,deriv_z

  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: rho_vp,rho_vs ! cpl,csl
  real(kind=CUSTOM_REAL) :: rhol,cpl,csl
  real(kind=CUSTOM_REAL) :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic

  ! for analytical initial plane wave for Bielak's conditions
  double precision :: veloc_horiz,veloc_vert,dxUx,dzUx,dxUz,dzUz,traction_x_t0,traction_z_t0
  integer :: count_left,count_right,count_bottom

  ! Stacey second order absorbing boundary scheme called P3
  ! implemented for spectral elements, see section 5.3 in:
  !    Casadei et al. (2002),
  !    A mortar spectral/finite element method for complex 2D and 3D elastodynamic problems,
  !    Comput. Methods Appl. Mech. Engrg. 191, 5119-5148.
  logical, parameter :: STACEY_SECOND_ORDER_P3 = .false.

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. anyabs) return
  if (.not. any_elastic) return

  ! Clayton-Engquist condition if elastic
  count_left = 1
  count_right = 1
  count_bottom = 1

  do ispecabs = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(ispecabs)

    ! only for elastic elements, skip others
    if (.not. ispec_is_elastic(ispec) ) cycle

    !--- left absorbing boundary
    if (codeabs(IEDGE4,ispecabs)) then
      i = 1
      do j = 1,NGLLZ
        ! Clayton-Engquist condition if elastic
        iglob = ibool(i,j,ispec)

        veloc_horiz = 0._CUSTOM_REAL
        veloc_vert = 0._CUSTOM_REAL
        traction_x_t0 = 0._CUSTOM_REAL
        traction_z_t0 = 0._CUSTOM_REAL

        ! for analytical initial plane wave for Bielak's conditions
        ! left or right edge, horizontal normal vector
        if (add_Bielak_conditions_left .and. initialfield) then
          if (.not. over_critical_angle) then
            call compute_Bielak_conditions(iglob,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                                           x_source(1), z_source(1), A_plane, B_plane, C_plane, &
                                           anglesource(1), anglesource_refl, &
                                           c_inc, c_refl, time_offset,f0_source(1))

            ! get elastic parameters of current grid point
            mul_unrelaxed_elastic = mustore(i,j,ispec)
            rhol = rhostore(i,j,ispec)
            cpl = rho_vpstore(i,j,ispec)/rhol

            lambdal_unrelaxed_elastic = rhol*cpl*cpl - 2._CUSTOM_REAL * mul_unrelaxed_elastic
            lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + 2._CUSTOM_REAL * mul_unrelaxed_elastic

            traction_x_t0 = lambdaplus2mu_unrelaxed_elastic * dxUx + lambdal_unrelaxed_elastic * dzUz
            traction_z_t0 = mul_unrelaxed_elastic * (dxUz + dzUx)
          else
            veloc_horiz = v0x_left(count_left,it)
            veloc_vert = v0z_left(count_left,it)
            traction_x_t0 = t0x_left(count_left,it)
            traction_z_t0 = t0z_left(count_left,it)
            count_left = count_left+1
          endif
        endif

        rho_vp = rho_vpstore(i,j,ispec)
        rho_vs = rho_vsstore(i,j,ispec)

        ! normal pointing left
        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D

        if (P_SV) then
          ! P_SV case
          vx = veloc_elastic(1,iglob) - veloc_horiz
          vz = veloc_elastic(2,iglob) - veloc_vert
          vn = nx*vx + nz*vz

          ! first-order Clayton-Engquist (1977) uses tractions
          !   T_normal = - rho * vp * veloc_normal
          ! with velocity's normal component: veloc_normal = vn * normal
          !   T_tangential = - rho * vs * veloc_tangential
          ! with velocity's tangential component: veloc_tangential = v - vn * normal
          ! total traction
          !    T = T_normal + T_tangential
          tx = rho_vp*vn*nx + rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz + rho_vs*(vz-vn*nz)

          ! second-order Stacey (1988), named P3
          !   T_tangential  = - rho * vs * veloc_tangential
          !   T_tangential2 = + rho * vs * (2 * vs - vp) * d/dx_tangential( displ_normal )
          ! and
          !   T_normal  = - rho * vp * veloc_normal
          !   T_normal2 = - rho * vs * (2 * vs - vp) * d/dx_tangential( displ_tangential )
          ! and total traction is
          !   T_normal + Tnormal2 + T_tangential + T_tangential2
          !
          ! derivatives: grad_tangential(f) is the dot-product grad(f) * vector_tangential
          if (STACEY_SECOND_ORDER_P3) then
            ! normal displacement
            dn = nx*displ_elastic(1,iglob) + nz*displ_elastic(2,iglob)
            dx_n = dn * nx
            dz_n = dn * nz

            ! tangential displacement
            dx_t = displ_elastic(1,iglob) - dx_n
            dz_t = displ_elastic(2,iglob) - dz_n

            ! tangential vector norm
            v_norm = sqrt(dx_t*dx_t + dz_t*dz_t)

            if (v_norm > 0.0_CUSTOM_REAL) then
              ! unit vector in tangential direction
              unit_tx = dx_t / v_norm
              unit_tz = dz_t / v_norm
            else
              ! unit vector in tangential direction
              ! t_tangential = (tx tz)^-1 = (nz -nx)^-1
              !   component tangential tx -> normal nz
              !   component tangential tz -> normal -nx
              unit_tx = nz
              unit_tz = - nx
            endif

            ! additional terms
            rhol = rhostore(i,j,ispec)
            cpl = rho_vpstore(i,j,ispec) / rhol
            csl = rho_vsstore(i,j,ispec) / rhol

            xixl = xix(i,j,ispec)
            xizl = xiz(i,j,ispec)
            gammaxl = gammax(i,j,ispec)
            gammazl = gammaz(i,j,ispec)

            ! displ_normal
            ! derivative along x
            tempx1l_n = 0._CUSTOM_REAL ! u_normal
            tempx2l_n = 0._CUSTOM_REAL
            tempx1l_t = 0._CUSTOM_REAL ! u_tangential
            tempx2l_t = 0._CUSTOM_REAL
            do k = 1,NGLLX
              iiglob = ibool(k,j,ispec)
              ! normal displacement
              dn = nx*displ_elastic(1,iiglob) + nz*displ_elastic(2,iiglob)
              dx_n = dn * nx
              dz_n = dn * nz
              ! tangential displacement
              dx_t = displ_elastic(1,iiglob) - dx_n
              dz_t = displ_elastic(2,iiglob) - dz_n

              tempx1l_n = tempx1l_n + dx_n * hprime_xx(i,k)
              tempx2l_n = tempx2l_n + dz_n * hprime_xx(i,k)

              tempx1l_t = tempx1l_t + dx_t * hprime_xx(i,k)
              tempx2l_t = tempx2l_t + dz_t * hprime_xx(i,k)
            enddo
            ! derivative along z
            tempz1l_n = 0._CUSTOM_REAL ! u_normal
            tempz2l_n = 0._CUSTOM_REAL
            tempz1l_t = 0._CUSTOM_REAL ! u_tangential
            tempz2l_t = 0._CUSTOM_REAL
            do k = 1,NGLLX
              iiglob = ibool(i,k,ispec)
              ! normal displacement
              dn = nx*displ_elastic(1,iiglob) + nz*displ_elastic(2,iiglob)
              dx_n = dn * nx
              dz_n = dn * nz
              ! tangential displacement
              dx_t = displ_elastic(1,iiglob) - dx_n
              dz_t = displ_elastic(2,iiglob) - dz_n

              tempz1l_n = tempz1l_n + dx_n * hprime_zz(j,k)
              tempz2l_n = tempz2l_n + dz_n * hprime_zz(j,k)

              tempz1l_t = tempz1l_t + dx_t * hprime_zz(j,k)
              tempz2l_t = tempz2l_t + dz_t * hprime_zz(j,k)
            enddo

            ! derivatives of displacement du_normal / du
            dux_dxl = (tempx1l_n*xixl + tempz1l_n*gammaxl)
            dux_dzl = (tempx1l_n*xizl + tempz1l_n*gammazl)

            duz_dxl = (tempx2l_n*xixl + tempz2l_n*gammaxl)
            duz_dzl = (tempx2l_n*xizl + tempz2l_n*gammazl)

            ! derivative du_normal / dx_tangential (along tangential vector)
            deriv_x = (dux_dxl * unit_tx + dux_dzl * unit_tz) * unit_tx
            deriv_z = (duz_dxl * unit_tx + duz_dzl * unit_tz) * unit_tz

            ! tangential traction: additional P3 term (minus sign added as tx will be subtracted)
            tx_t = - rho_vs * (2.0 * csl - cpl) * deriv_x
            tz_t = - rho_vs * (2.0 * csl - cpl) * deriv_z

            ! derivatives of displacement du_tangential / du
            dux_dxl = (tempx1l_t*xixl + tempz1l_t*gammaxl)
            dux_dzl = (tempx1l_t*xizl + tempz1l_t*gammazl)

            duz_dxl = (tempx2l_t*xixl + tempz2l_t*gammaxl)
            duz_dzl = (tempx2l_t*xizl + tempz2l_t*gammazl)

            ! derivative du_tangential / dx_tangential (along tangential vector)
            deriv_x = (dux_dxl * unit_tx + dux_dzl * unit_tz) * unit_tx
            deriv_z = (duz_dxl * unit_tx + duz_dzl * unit_tz) * unit_tz

            ! normal traction: additional P3 term
            tx_n = - rho_vs * (2.0 * csl - cpl) * deriv_x
            tz_n = - rho_vs * (2.0 * csl - cpl) * deriv_z

            !if (j==1 .and. ispec==100) print *,'debug: ',ispecabs,ispec,tx,tz,'n',tx_n,tz_n,'t',tx_t,tz_t!,'v',v_norm

            ! P3: adds second-order terms
            tx = tx + tx_t + tx_n
            tz = tz + tz_t + tz_n
          endif  ! STACEY_SECOND_ORDER_P3

        else
          ! SH case
          vy = veloc_elastic(1,iglob)
          ty = rho_vs*vy
        endif

        weight = jacobian1D * wzgll(j)

        if (P_SV) then
          ! P_SV case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + traction_x_t0)*weight
          accel_elastic(2,iglob) = accel_elastic(2,iglob) - (tz + traction_z_t0)*weight
        else
          ! SH case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - ty*weight
        endif

        if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
          if (P_SV) then
            ! P-SV waves
            b_absorb_elastic_left(1,j,ib_left(ispecabs),it) = (tx + traction_x_t0)*weight
            b_absorb_elastic_left(2,j,ib_left(ispecabs),it) = (tz + traction_z_t0)*weight
          else
            ! SH (membrane) waves
            b_absorb_elastic_left(1,j,ib_left(ispecabs),it) = ty*weight
          endif
        endif
      enddo
    endif  !  end of left absorbing boundary

    !--- right absorbing boundary
    if (codeabs(IEDGE2,ispecabs)) then
      i = NGLLX
      do j = 1,NGLLZ
        ! Clayton-Engquist condition if elastic
        iglob = ibool(i,j,ispec)

        veloc_horiz = 0._CUSTOM_REAL
        veloc_vert = 0._CUSTOM_REAL
        traction_x_t0 = 0._CUSTOM_REAL
        traction_z_t0 = 0._CUSTOM_REAL

        ! for analytical initial plane wave for Bielak's conditions
        ! left or right edge, horizontal normal vector
        if (add_Bielak_conditions_right .and. initialfield) then
          if (.not. over_critical_angle) then
            call compute_Bielak_conditions(iglob,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                                           x_source(1), z_source(1), A_plane, B_plane, C_plane, &
                                           anglesource(1), anglesource_refl, &
                                           c_inc, c_refl, time_offset,f0_source(1))

            ! get elastic parameters of current grid point
            mul_unrelaxed_elastic = mustore(i,j,ispec)
            rhol = rhostore(i,j,ispec)
            cpl = rho_vpstore(i,j,ispec)/rhol

            lambdal_unrelaxed_elastic = rhol*cpl*cpl - 2._CUSTOM_REAL * mul_unrelaxed_elastic
            lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + 2._CUSTOM_REAL * mul_unrelaxed_elastic

            traction_x_t0 = lambdaplus2mu_unrelaxed_elastic * dxUx + lambdal_unrelaxed_elastic * dzUz
            traction_z_t0 = mul_unrelaxed_elastic * (dxUz + dzUx)
          else
            veloc_horiz = v0x_right(count_right,it)
            veloc_vert = v0z_right(count_right,it)
            traction_x_t0 = t0x_right(count_right,it)
            traction_z_t0 = t0z_right(count_right,it)
            count_right = count_right+1
          endif
        endif

        rho_vp = rho_vpstore(i,j,ispec)
        rho_vs = rho_vsstore(i,j,ispec)

        ! normal pointing right
        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D

        if (P_SV) then
          ! P_SV case
          vx = veloc_elastic(1,iglob) - veloc_horiz
          vz = veloc_elastic(2,iglob) - veloc_vert
          vn = nx*vx + nz*vz
          tx = rho_vp*vn*nx + rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz + rho_vs*(vz-vn*nz)

          ! second-order Stacey (1988), named P3
          !   T_tangential  = - rho * vs * veloc_tangential
          !   T_tangential2 = + rho * vs * (2 * vs - vp) * d/dx_tangential( displ_normal )
          ! and
          !   T_normal  = - rho * vp * veloc_normal
          !   T_normal2 = - rho * vs * (2 * vs - vp) * d/dx_tangential( displ_tangential )
          ! and total traction is
          !   T_normal + Tnormal2 + T_tangential + T_tangential2
          !
          ! derivatives: grad_tangential(f) is the dot-product grad(f) * vector_tangential
          if (STACEY_SECOND_ORDER_P3) then
            ! normal displacement
            dn = nx*displ_elastic(1,iglob) + nz*displ_elastic(2,iglob)
            dx_n = dn * nx
            dz_n = dn * nz

            ! tangential displacement
            dx_t = displ_elastic(1,iglob) - dx_n
            dz_t = displ_elastic(2,iglob) - dz_n

            ! tangential vector norm
            v_norm = sqrt(dx_t*dx_t + dz_t*dz_t)

            if (v_norm > 0.0_CUSTOM_REAL) then
              ! unit vector in tangential direction
              unit_tx = dx_t / v_norm
              unit_tz = dz_t / v_norm
            else
              ! unit vector in tangential direction
              ! t_tangential = (tx tz)^-1 = (nz -nx)^-1
              !   component tangential tx -> normal nz
              !   component tangential tz -> normal -nx
              unit_tx = nz
              unit_tz = - nx
            endif

            ! additional terms
            rhol = rhostore(i,j,ispec)
            cpl = rho_vpstore(i,j,ispec) / rhol
            csl = rho_vsstore(i,j,ispec) / rhol

            xixl = xix(i,j,ispec)
            xizl = xiz(i,j,ispec)
            gammaxl = gammax(i,j,ispec)
            gammazl = gammaz(i,j,ispec)

            ! displ_normal
            ! derivative along x
            tempx1l_n = 0._CUSTOM_REAL ! u_normal
            tempx2l_n = 0._CUSTOM_REAL
            tempx1l_t = 0._CUSTOM_REAL ! u_tangential
            tempx2l_t = 0._CUSTOM_REAL
            do k = 1,NGLLX
              iiglob = ibool(k,j,ispec)
              ! normal displacement
              dn = nx*displ_elastic(1,iiglob) + nz*displ_elastic(2,iiglob)
              dx_n = dn * nx
              dz_n = dn * nz
              ! tangential displacement
              dx_t = displ_elastic(1,iiglob) - dx_n
              dz_t = displ_elastic(2,iiglob) - dz_n

              tempx1l_n = tempx1l_n + dx_n * hprime_xx(i,k)
              tempx2l_n = tempx2l_n + dz_n * hprime_xx(i,k)

              tempx1l_t = tempx1l_t + dx_t * hprime_xx(i,k)
              tempx2l_t = tempx2l_t + dz_t * hprime_xx(i,k)
            enddo
            ! derivative along z
            tempz1l_n = 0._CUSTOM_REAL ! u_normal
            tempz2l_n = 0._CUSTOM_REAL
            tempz1l_t = 0._CUSTOM_REAL ! u_tangential
            tempz2l_t = 0._CUSTOM_REAL
            do k = 1,NGLLX
              iiglob = ibool(i,k,ispec)
              ! normal displacement
              dn = nx*displ_elastic(1,iiglob) + nz*displ_elastic(2,iiglob)
              dx_n = dn * nx
              dz_n = dn * nz
              ! tangential displacement
              dx_t = displ_elastic(1,iiglob) - dx_n
              dz_t = displ_elastic(2,iiglob) - dz_n

              tempz1l_n = tempz1l_n + dx_n * hprime_zz(j,k)
              tempz2l_n = tempz2l_n + dz_n * hprime_zz(j,k)

              tempz1l_t = tempz1l_t + dx_t * hprime_zz(j,k)
              tempz2l_t = tempz2l_t + dz_t * hprime_zz(j,k)
            enddo

            ! derivatives of displacement du_normal / du
            dux_dxl = (tempx1l_n*xixl + tempz1l_n*gammaxl)
            dux_dzl = (tempx1l_n*xizl + tempz1l_n*gammazl)

            duz_dxl = (tempx2l_n*xixl + tempz2l_n*gammaxl)
            duz_dzl = (tempx2l_n*xizl + tempz2l_n*gammazl)

            ! derivative du_normal / dx_tangential (along tangential vector)
            deriv_x = (dux_dxl * unit_tx + dux_dzl * unit_tz) * unit_tx
            deriv_z = (duz_dxl * unit_tx + duz_dzl * unit_tz) * unit_tz

            ! tangential traction: additional P3 term (minus sign added as tx will be subtracted)
            tx_t = - rho_vs * (2.0 * csl - cpl) * deriv_x
            tz_t = - rho_vs * (2.0 * csl - cpl) * deriv_z

            ! derivatives of displacement du_tangential / du
            dux_dxl = (tempx1l_t*xixl + tempz1l_t*gammaxl)
            dux_dzl = (tempx1l_t*xizl + tempz1l_t*gammazl)

            duz_dxl = (tempx2l_t*xixl + tempz2l_t*gammaxl)
            duz_dzl = (tempx2l_t*xizl + tempz2l_t*gammazl)

            ! derivative du_tangential / dx_tangential (along tangential vector)
            deriv_x = (dux_dxl * unit_tx + dux_dzl * unit_tz) * unit_tx
            deriv_z = (duz_dxl * unit_tx + duz_dzl * unit_tz) * unit_tz

            ! normal traction: additional P3 term
            tx_n = - rho_vs * (2.0 * csl - cpl) * deriv_x
            tz_n = - rho_vs * (2.0 * csl - cpl) * deriv_z

            !if (j==1 .and. ispec==100) print *,'debug: ',ispecabs,ispec,tx,tz,'n',tx_n,tz_n,'t',tx_t,tz_t!,'v',v_norm

            ! P3: adds second-order terms
            tx = tx + tx_t + tx_n
            tz = tz + tz_t + tz_n
          endif  ! STACEY_SECOND_ORDER_P3

        else
          ! SH case
          vy = veloc_elastic(1,iglob)
          ty = rho_vs*vy
        endif

        weight = jacobian1D * wzgll(j)

        if (P_SV) then
          ! P_SV case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx - traction_x_t0)*weight
          accel_elastic(2,iglob) = accel_elastic(2,iglob) - (tz - traction_z_t0)*weight
        else
          ! SH case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - ty*weight
        endif

        if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. NO_BACKWARD_RECONSTRUCTION)) then
          if (P_SV) then
            ! P-SV waves
            b_absorb_elastic_right(1,j,ib_right(ispecabs),it) = (tx - traction_x_t0)*weight
            b_absorb_elastic_right(2,j,ib_right(ispecabs),it) = (tz - traction_z_t0)*weight
          else
            ! SH (membrane) waves
            b_absorb_elastic_right(1,j,ib_right(ispecabs),it) = ty*weight
          endif
        endif
      enddo
    endif  !  end of right absorbing boundary

    !--- bottom absorbing boundary
    if (codeabs(IEDGE1,ispecabs)) then
      j = 1
      do i = 1,NGLLX
        ! Clayton-Engquist condition if elastic
        iglob = ibool(i,j,ispec)

        veloc_horiz = 0._CUSTOM_REAL
        veloc_vert = 0._CUSTOM_REAL
        traction_x_t0 = 0._CUSTOM_REAL
        traction_z_t0 = 0._CUSTOM_REAL

        ! for analytical initial plane wave for Bielak's conditions
        ! top or bottom edge, vertical normal vector
        if (add_Bielak_conditions_bottom .and. initialfield) then
          if (.not. over_critical_angle) then
            call compute_Bielak_conditions(iglob,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                                           x_source(1), z_source(1), A_plane, B_plane, C_plane, &
                                           anglesource(1), anglesource_refl, &
                                           c_inc, c_refl, time_offset,f0_source(1))

            ! get elastic parameters of current grid point
            mul_unrelaxed_elastic = mustore(i,j,ispec)
            rhol = rhostore(i,j,ispec)
            cpl = rho_vpstore(i,j,ispec)/rhol

            lambdal_unrelaxed_elastic = rhol*cpl*cpl - 2._CUSTOM_REAL * mul_unrelaxed_elastic
            lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + 2._CUSTOM_REAL * mul_unrelaxed_elastic

            traction_x_t0 = mul_unrelaxed_elastic * (dxUz + dzUx)
            traction_z_t0 = lambdal_unrelaxed_elastic * dxUx + lambdaplus2mu_unrelaxed_elastic * dzUz
          else
            veloc_horiz = v0x_bot(count_bottom,it)
            veloc_vert = v0z_bot(count_bottom,it)
            traction_x_t0 = t0x_bot(count_bottom,it)
            traction_z_t0 = t0z_bot(count_bottom,it)
            count_bottom = count_bottom+1
          endif
        endif

        rho_vp = rho_vpstore(i,j,ispec)
        rho_vs = rho_vsstore(i,j,ispec)

        ! normal pointing down
        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D

        if (P_SV) then
          ! P_SV case
          vx = veloc_elastic(1,iglob) - veloc_horiz
          vz = veloc_elastic(2,iglob) - veloc_vert
          vn = nx*vx+nz*vz
          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

          ! second-order Stacey (1988), named P3
          !   T_tangential  = - rho * vs * veloc_tangential
          !   T_tangential2 = + rho * vs * (2 * vs - vp) * d/dx_tangential( displ_normal )
          ! and
          !   T_normal  = - rho * vp * veloc_normal
          !   T_normal2 = - rho * vs * (2 * vs - vp) * d/dx_tangential( displ_tangential )
          ! and total traction is
          !   T_normal + Tnormal2 + T_tangential + T_tangential2
          !
          ! derivatives: grad_tangential(f) is the dot-product grad(f) * vector_tangential
          if (STACEY_SECOND_ORDER_P3) then
            ! normal displacement
            dn = nx*displ_elastic(1,iglob) + nz*displ_elastic(2,iglob)
            dx_n = dn * nx
            dz_n = dn * nz

            ! tangential displacement
            dx_t = displ_elastic(1,iglob) - dx_n
            dz_t = displ_elastic(2,iglob) - dz_n

            ! tangential vector norm
            v_norm = sqrt(dx_t*dx_t + dz_t*dz_t)

            if (v_norm > 0.0_CUSTOM_REAL) then
              ! unit vector in tangential direction
              unit_tx = dx_t / v_norm
              unit_tz = dz_t / v_norm
            else
              ! unit vector in tangential direction
              ! t_tangential = (tx tz)^-1 = (nz -nx)^-1
              !   component tangential tx -> normal nz
              !   component tangential tz -> normal -nx
              unit_tx = nz
              unit_tz = - nx
            endif

            ! additional terms
            rhol = rhostore(i,j,ispec)
            cpl = rho_vpstore(i,j,ispec) / rhol
            csl = rho_vsstore(i,j,ispec) / rhol

            xixl = xix(i,j,ispec)
            xizl = xiz(i,j,ispec)
            gammaxl = gammax(i,j,ispec)
            gammazl = gammaz(i,j,ispec)

            ! displ_normal
            ! derivative along x
            tempx1l_n = 0._CUSTOM_REAL ! u_normal
            tempx2l_n = 0._CUSTOM_REAL
            tempx1l_t = 0._CUSTOM_REAL ! u_tangential
            tempx2l_t = 0._CUSTOM_REAL
            do k = 1,NGLLX
              iiglob = ibool(k,j,ispec)
              ! normal displacement
              dn = nx*displ_elastic(1,iiglob) + nz*displ_elastic(2,iiglob)
              dx_n = dn * nx
              dz_n = dn * nz
              ! tangential displacement
              dx_t = displ_elastic(1,iiglob) - dx_n
              dz_t = displ_elastic(2,iiglob) - dz_n

              tempx1l_n = tempx1l_n + dx_n * hprime_xx(i,k)
              tempx2l_n = tempx2l_n + dz_n * hprime_xx(i,k)

              tempx1l_t = tempx1l_t + dx_t * hprime_xx(i,k)
              tempx2l_t = tempx2l_t + dz_t * hprime_xx(i,k)
            enddo
            ! derivative along z
            tempz1l_n = 0._CUSTOM_REAL ! u_normal
            tempz2l_n = 0._CUSTOM_REAL
            tempz1l_t = 0._CUSTOM_REAL ! u_tangential
            tempz2l_t = 0._CUSTOM_REAL
            do k = 1,NGLLX
              iiglob = ibool(i,k,ispec)
              ! normal displacement
              dn = nx*displ_elastic(1,iiglob) + nz*displ_elastic(2,iiglob)
              dx_n = dn * nx
              dz_n = dn * nz
              ! tangential displacement
              dx_t = displ_elastic(1,iiglob) - dx_n
              dz_t = displ_elastic(2,iiglob) - dz_n

              tempz1l_n = tempz1l_n + dx_n * hprime_zz(j,k)
              tempz2l_n = tempz2l_n + dz_n * hprime_zz(j,k)

              tempz1l_t = tempz1l_t + dx_t * hprime_zz(j,k)
              tempz2l_t = tempz2l_t + dz_t * hprime_zz(j,k)
            enddo

            ! derivatives of displacement du_normal / du
            dux_dxl = (tempx1l_n*xixl + tempz1l_n*gammaxl)
            dux_dzl = (tempx1l_n*xizl + tempz1l_n*gammazl)

            duz_dxl = (tempx2l_n*xixl + tempz2l_n*gammaxl)
            duz_dzl = (tempx2l_n*xizl + tempz2l_n*gammazl)

            ! derivative du_normal / dx_tangential (along tangential vector)
            deriv_x = (dux_dxl * unit_tx + dux_dzl * unit_tz) * unit_tx
            deriv_z = (duz_dxl * unit_tx + duz_dzl * unit_tz) * unit_tz

            ! tangential traction: additional P3 term (minus sign added as tx will be subtracted)
            tx_t = - rho_vs * (2.0 * csl - cpl) * deriv_x
            tz_t = - rho_vs * (2.0 * csl - cpl) * deriv_z

            ! derivatives of displacement du_tangential / du
            dux_dxl = (tempx1l_t*xixl + tempz1l_t*gammaxl)
            dux_dzl = (tempx1l_t*xizl + tempz1l_t*gammazl)

            duz_dxl = (tempx2l_t*xixl + tempz2l_t*gammaxl)
            duz_dzl = (tempx2l_t*xizl + tempz2l_t*gammazl)

            ! derivative du_tangential / dx_tangential (along tangential vector)
            deriv_x = (dux_dxl * unit_tx + dux_dzl * unit_tz) * unit_tx
            deriv_z = (duz_dxl * unit_tx + duz_dzl * unit_tz) * unit_tz

            ! normal traction: additional P3 term
            tx_n = - rho_vs * (2.0 * csl - cpl) * deriv_x
            tz_n = - rho_vs * (2.0 * csl - cpl) * deriv_z

            !if (j==1 .and. ispec==100) print *,'debug: ',ispecabs,ispec,tx,tz,'n',tx_n,tz_n,'t',tx_t,tz_t!,'v',v_norm

            ! P3: adds second-order terms
            tx = tx + tx_t + tx_n
            tz = tz + tz_t + tz_n
          endif  ! STACEY_SECOND_ORDER_P3

        else
          ! SH case
          vy = veloc_elastic(1,iglob)
          ty = rho_vs*vy
        endif

! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
        if ((codeabs_corner(1,ispecabs) .and. i == 1) .or. (codeabs_corner(2,ispecabs) .and. i == NGLLX)) then
          tx = 0._CUSTOM_REAL
          ty = 0._CUSTOM_REAL
          tz = 0._CUSTOM_REAL
        endif

        weight = jacobian1D * wxgll(i)

        if (P_SV) then
          ! P_SV case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + traction_x_t0)*weight
          accel_elastic(2,iglob) = accel_elastic(2,iglob) - (tz + traction_z_t0)*weight
        else
          ! SH case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - ty*weight
        endif

        if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
          if (P_SV) then
            ! P-SV waves
            b_absorb_elastic_bottom(1,i,ib_bottom(ispecabs),it) = (tx + traction_x_t0)*weight
            b_absorb_elastic_bottom(2,i,ib_bottom(ispecabs),it) = (tz + traction_z_t0)*weight
          else
            ! SH (membrane) waves
            b_absorb_elastic_bottom(1,i,ib_bottom(ispecabs),it) = ty*weight
          endif
        endif
      enddo
    endif  !  end of bottom absorbing boundary

    !--- top absorbing boundary
    if (codeabs(IEDGE3,ispecabs)) then
      j = NGLLZ
      do i = 1,NGLLX
        ! Clayton-Engquist condition if elastic
        iglob = ibool(i,j,ispec)

        veloc_horiz = 0._CUSTOM_REAL
        veloc_vert = 0._CUSTOM_REAL
        traction_x_t0 = 0._CUSTOM_REAL
        traction_z_t0 = 0._CUSTOM_REAL

        ! for analytical initial plane wave for Bielak's conditions
        ! top or bottom edge, vertical normal vector
        if (add_Bielak_conditions_top .and. initialfield) then
        ! at the top there is no test for whether we are above the critical angle
        ! because a critical angle can only exist when the top edge is a free surface, not in an infinite medium
          call compute_Bielak_conditions(iglob,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                                         x_source(1), z_source(1), A_plane, B_plane, C_plane, &
                                         anglesource(1), anglesource_refl, &
                                         c_inc, c_refl, time_offset,f0_source(1))

          ! get elastic parameters of current grid point
          mul_unrelaxed_elastic = mustore(i,j,ispec)
          rhol = rhostore(i,j,ispec)
          cpl = rho_vpstore(i,j,ispec)/rhol

          lambdal_unrelaxed_elastic = rhol*cpl*cpl - 2._CUSTOM_REAL * mul_unrelaxed_elastic
          lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + 2._CUSTOM_REAL * mul_unrelaxed_elastic

          traction_x_t0 = mul_unrelaxed_elastic * (dxUz + dzUx)
          traction_z_t0 = lambdal_unrelaxed_elastic * dxUx + lambdaplus2mu_unrelaxed_elastic * dzUz
        endif

        rho_vp = rho_vpstore(i,j,ispec)
        rho_vs = rho_vsstore(i,j,ispec)

        ! normal pointing up
        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D

        if (P_SV) then
          ! P_SV case
          vx = veloc_elastic(1,iglob) - veloc_horiz
          vz = veloc_elastic(2,iglob) - veloc_vert
          vn = nx*vx+nz*vz
          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

          ! second-order Stacey (1988), named P3
          !   T_tangential  = - rho * vs * veloc_tangential
          !   T_tangential2 = + rho * vs * (2 * vs - vp) * d/dx_tangential( displ_normal )
          ! and
          !   T_normal  = - rho * vp * veloc_normal
          !   T_normal2 = - rho * vs * (2 * vs - vp) * d/dx_tangential( displ_tangential )
          ! and total traction is
          !   T_normal + Tnormal2 + T_tangential + T_tangential2
          !
          ! derivatives: grad_tangential(f) is the dot-product grad(f) * vector_tangential
          if (STACEY_SECOND_ORDER_P3) then
            ! normal displacement
            dn = nx*displ_elastic(1,iglob) + nz*displ_elastic(2,iglob)
            dx_n = dn * nx
            dz_n = dn * nz

            ! tangential displacement
            dx_t = displ_elastic(1,iglob) - dx_n
            dz_t = displ_elastic(2,iglob) - dz_n

            ! tangential vector norm
            v_norm = sqrt(dx_t*dx_t + dz_t*dz_t)

            if (v_norm > 0.0_CUSTOM_REAL) then
              ! unit vector in tangential direction
              unit_tx = dx_t / v_norm
              unit_tz = dz_t / v_norm
            else
              ! unit vector in tangential direction
              ! t_tangential = (tx tz)^-1 = (nz -nx)^-1
              !   component tangential tx -> normal nz
              !   component tangential tz -> normal -nx
              unit_tx = nz
              unit_tz = - nx
            endif

            ! additional terms
            rhol = rhostore(i,j,ispec)
            cpl = rho_vpstore(i,j,ispec) / rhol
            csl = rho_vsstore(i,j,ispec) / rhol

            xixl = xix(i,j,ispec)
            xizl = xiz(i,j,ispec)
            gammaxl = gammax(i,j,ispec)
            gammazl = gammaz(i,j,ispec)

            ! displ_normal
            ! derivative along x
            tempx1l_n = 0._CUSTOM_REAL ! u_normal
            tempx2l_n = 0._CUSTOM_REAL
            tempx1l_t = 0._CUSTOM_REAL ! u_tangential
            tempx2l_t = 0._CUSTOM_REAL
            do k = 1,NGLLX
              iiglob = ibool(k,j,ispec)
              ! normal displacement
              dn = nx*displ_elastic(1,iiglob) + nz*displ_elastic(2,iiglob)
              dx_n = dn * nx
              dz_n = dn * nz
              ! tangential displacement
              dx_t = displ_elastic(1,iiglob) - dx_n
              dz_t = displ_elastic(2,iiglob) - dz_n

              tempx1l_n = tempx1l_n + dx_n * hprime_xx(i,k)
              tempx2l_n = tempx2l_n + dz_n * hprime_xx(i,k)

              tempx1l_t = tempx1l_t + dx_t * hprime_xx(i,k)
              tempx2l_t = tempx2l_t + dz_t * hprime_xx(i,k)
            enddo
            ! derivative along z
            tempz1l_n = 0._CUSTOM_REAL ! u_normal
            tempz2l_n = 0._CUSTOM_REAL
            tempz1l_t = 0._CUSTOM_REAL ! u_tangential
            tempz2l_t = 0._CUSTOM_REAL
            do k = 1,NGLLX
              iiglob = ibool(i,k,ispec)
              ! normal displacement
              dn = nx*displ_elastic(1,iiglob) + nz*displ_elastic(2,iiglob)
              dx_n = dn * nx
              dz_n = dn * nz
              ! tangential displacement
              dx_t = displ_elastic(1,iiglob) - dx_n
              dz_t = displ_elastic(2,iiglob) - dz_n

              tempz1l_n = tempz1l_n + dx_n * hprime_zz(j,k)
              tempz2l_n = tempz2l_n + dz_n * hprime_zz(j,k)

              tempz1l_t = tempz1l_t + dx_t * hprime_zz(j,k)
              tempz2l_t = tempz2l_t + dz_t * hprime_zz(j,k)
            enddo

            ! derivatives of displacement du_normal / du
            dux_dxl = (tempx1l_n*xixl + tempz1l_n*gammaxl)
            dux_dzl = (tempx1l_n*xizl + tempz1l_n*gammazl)

            duz_dxl = (tempx2l_n*xixl + tempz2l_n*gammaxl)
            duz_dzl = (tempx2l_n*xizl + tempz2l_n*gammazl)

            ! derivative du_normal / dx_tangential (along tangential vector)
            deriv_x = (dux_dxl * unit_tx + dux_dzl * unit_tz) * unit_tx
            deriv_z = (duz_dxl * unit_tx + duz_dzl * unit_tz) * unit_tz

            ! tangential traction: additional P3 term (minus sign added as tx will be subtracted)
            tx_t = - rho_vs * (2.0 * csl - cpl) * deriv_x
            tz_t = - rho_vs * (2.0 * csl - cpl) * deriv_z

            ! derivatives of displacement du_tangential / du
            dux_dxl = (tempx1l_t*xixl + tempz1l_t*gammaxl)
            dux_dzl = (tempx1l_t*xizl + tempz1l_t*gammazl)

            duz_dxl = (tempx2l_t*xixl + tempz2l_t*gammaxl)
            duz_dzl = (tempx2l_t*xizl + tempz2l_t*gammazl)

            ! derivative du_tangential / dx_tangential (along tangential vector)
            deriv_x = (dux_dxl * unit_tx + dux_dzl * unit_tz) * unit_tx
            deriv_z = (duz_dxl * unit_tx + duz_dzl * unit_tz) * unit_tz

            ! normal traction: additional P3 term
            tx_n = - rho_vs * (2.0 * csl - cpl) * deriv_x
            tz_n = - rho_vs * (2.0 * csl - cpl) * deriv_z

            !if (j==1 .and. ispec==100) print *,'debug: ',ispecabs,ispec,tx,tz,'n',tx_n,tz_n,'t',tx_t,tz_t!,'v',v_norm

            ! P3: adds second-order terms
            tx = tx + tx_t + tx_n
            tz = tz + tz_t + tz_n
          endif  ! STACEY_SECOND_ORDER_P3

        else
          ! SH case
          vy = veloc_elastic(1,iglob)
          ty = rho_vs*vy
        endif

! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
        if ((codeabs_corner(3,ispecabs) .and. i == 1) .or. (codeabs_corner(4,ispecabs) .and. i == NGLLX)) then
          tx = 0._CUSTOM_REAL
          ty = 0._CUSTOM_REAL
          tz = 0._CUSTOM_REAL
        endif

        weight = jacobian1D * wxgll(i)

        if (P_SV) then
          ! P_SV case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx - traction_x_t0)*weight
          accel_elastic(2,iglob) = accel_elastic(2,iglob) - (tz - traction_z_t0)*weight
        else
          ! SH case
          accel_elastic(1,iglob) = accel_elastic(1,iglob) - ty*weight
        endif

        if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
          if (P_SV) then
            ! P-SV waves
            b_absorb_elastic_top(1,i,ib_top(ispecabs),it) = (tx - traction_x_t0)*weight
            b_absorb_elastic_top(2,i,ib_top(ispecabs),it) = (tz - traction_z_t0)*weight
          else
            ! SH (membrane) waves
            b_absorb_elastic_top(1,i,ib_top(ispecabs),it) = ty*weight
          endif
        endif
      enddo
    endif  !  end of top absorbing boundary

  enddo

  end subroutine compute_stacey_elastic

!
!------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_elastic_backward(b_accel_elastic)

! absorbing boundaries

! Clayton-Engquist condition if elastic

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: nglob,any_elastic,ibool,ispec_is_elastic, &
                         NSTEP,it,num_abs_boundary_faces,anyabs, &
                         abs_boundary_ispec,codeabs, &
                         b_absorb_elastic_left,b_absorb_elastic_right, &
                         b_absorb_elastic_bottom,b_absorb_elastic_top, &
                         ib_left,ib_right,ib_bottom,ib_top, &
                         STACEY_ABSORBING_CONDITIONS,P_SV,NO_BACKWARD_RECONSTRUCTION

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(inout) :: b_accel_elastic

  ! local parameters
  integer :: ispecabs,ispec,i,j,iglob,it_tmp

  ! checks if anything to do
  if (.not. STACEY_ABSORBING_CONDITIONS) return
  if (.not. anyabs) return
  if (.not. any_elastic) return
  if (NO_BACKWARD_RECONSTRUCTION) return

  ! time increment index
  it_tmp = NSTEP - it + 1

  ! Clayton-Engquist condition if elastic
  do ispecabs = 1,num_abs_boundary_faces

    ispec = abs_boundary_ispec(ispecabs)
    if (.not. ispec_is_elastic(ispec) ) cycle

    !--- left absorbing boundary
    if (codeabs(IEDGE4,ispecabs)) then
      i = 1
      do j = 1,NGLLZ
        ! Clayton-Engquist condition if elastic
        iglob = ibool(i,j,ispec)
        if (P_SV) then
          ! P-SV waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_left(1,j,ib_left(ispecabs),it_tmp)
          b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - b_absorb_elastic_left(2,j,ib_left(ispecabs),it_tmp)
        else
          ! SH (membrane) waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_left(1,j,ib_left(ispecabs),it_tmp)
        endif
      enddo
    endif  !  end of left absorbing boundary

    !--- right absorbing boundary
    if (codeabs(IEDGE2,ispecabs)) then
      i = NGLLX
      do j = 1,NGLLZ
        iglob = ibool(i,j,ispec)
        if (P_SV) then
          ! P-SV waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_right(1,j,ib_right(ispecabs),it_tmp)
          b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - b_absorb_elastic_right(2,j,ib_right(ispecabs),it_tmp)
        else
          ! SH (membrane) waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_right(1,j,ib_right(ispecabs),it_tmp)
        endif
      enddo
    endif  !  end of right absorbing boundary

    !--- bottom absorbing boundary
    if (codeabs(IEDGE1,ispecabs)) then
      j = 1
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        if (P_SV) then
          ! P-SV waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_bottom(1,i,ib_bottom(ispecabs),it_tmp)
          b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - b_absorb_elastic_bottom(2,i,ib_bottom(ispecabs),it_tmp)
        else
          ! SH (membrane) waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_bottom(1,i,ib_bottom(ispecabs),it_tmp)
        endif
      enddo
    endif  !  end of bottom absorbing boundary

    !--- top absorbing boundary
    if (codeabs(IEDGE3,ispecabs)) then
      j = NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        if (P_SV) then
          ! P-SV waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_top(1,i,ib_top(ispecabs),it_tmp)
          b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - b_absorb_elastic_top(2,i,ib_top(ispecabs),it_tmp)
        else
          ! SH (membrane) waves
          b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - b_absorb_elastic_top(1,i,ib_top(ispecabs),it_tmp)
        endif
      enddo
    endif  !  end of top absorbing boundary

  enddo

  end subroutine compute_stacey_elastic_backward

!
!------------------------------------------------------------------------------------------
!
!
! this routine is in principle unused... left here for reference
!
!  subroutine UNUSED_compute_gradient_field_element(ispec,field,dux_dxl,dux_dzl,duz_dxl,duz_dzl)
!
!  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM
!
!  use specfem_par, only: nglob,ibool,xix,xiz,gammax,gammaz
!  use specfem_par, only: hprime_xx,hprime_zz
!  use specfem_par, only: AXISYM,is_on_the_axis,hprimeBar_xx,NGLJ
!
!  implicit none
!
!  integer,intent(in) :: ispec
!  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: field
!
!  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ),intent(out) :: dux_dxl,dux_dzl,duz_dxl,duz_dzl
!
!  ! local parameters
!  integer :: i,j,k,iglob
!  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl
!  real(kind=CUSTOM_REAL) :: dux_dxi,duz_dxi,dux_dgamma,duz_dgamma
!
!
!  do j = 1,NGLLZ
!    do i = 1,NGLLX
!      ! derivative along x and along z
!      dux_dxi = 0._CUSTOM_REAL
!      duz_dxi = 0._CUSTOM_REAL
!
!      dux_dgamma = 0._CUSTOM_REAL
!      duz_dgamma = 0._CUSTOM_REAL
!
!      ! first double loop over GLL points to compute and store gradients
!      ! we can merge the two loops because NGLLX == NGLLZ
!      if (AXISYM) then
!        if (is_on_the_axis(ispec)) then
!          do k = 1,NGLJ
!            ! derivative along x
!            iglob = ibool(k,j,ispec)
!            dux_dxi = dux_dxi + field(1,iglob)*hprimeBar_xx(i,k)
!            duz_dxi = duz_dxi + field(2,iglob)*hprimeBar_xx(i,k)
!            ! derivative along z
!            iglob = ibool(i,k,ispec)
!            dux_dgamma = dux_dgamma + field(1,iglob)*hprime_zz(j,k)
!            duz_dgamma = duz_dgamma + field(2,iglob)*hprime_zz(j,k)
!          enddo
!        else
!          do k = 1,NGLJ
!            ! derivative along x
!            iglob = ibool(k,j,ispec)
!            dux_dxi = dux_dxi + field(1,iglob)*hprime_xx(i,k)
!            duz_dxi = duz_dxi + field(2,iglob)*hprime_xx(i,k)
!            ! derivative along z
!            iglob = ibool(i,k,ispec)
!            dux_dgamma = dux_dgamma + field(1,iglob)*hprime_zz(j,k)
!            duz_dgamma = duz_dgamma + field(2,iglob)*hprime_zz(j,k)
!          enddo
!        endif
!      else
!        do k = 1,NGLLX
!          ! derivative along x
!          iglob = ibool(k,j,ispec)
!          dux_dxi = dux_dxi + field(1,iglob)*hprime_xx(i,k)
!          duz_dxi = duz_dxi + field(2,iglob)*hprime_xx(i,k)
!          ! derivative along z
!          iglob = ibool(i,k,ispec)
!          dux_dgamma = dux_dgamma + field(1,iglob)*hprime_zz(j,k)
!          duz_dgamma = duz_dgamma + field(2,iglob)*hprime_zz(j,k)
!        enddo
!      endif
!
!      xixl = xix(i,j,ispec)
!      xizl = xiz(i,j,ispec)
!      gammaxl = gammax(i,j,ispec)
!      gammazl = gammaz(i,j,ispec)
!
!      ! derivatives in x and z directions
!      dux_dxl(i,j) = dux_dxi*xixl + dux_dgamma*gammaxl
!      dux_dzl(i,j) = dux_dxi*xizl + dux_dgamma*gammazl
!
!      duz_dxl(i,j) = duz_dxi*xixl + duz_dgamma*gammaxl
!      duz_dzl(i,j) = duz_dxi*xizl + duz_dgamma*gammazl
!
!    enddo
!  enddo
!
!  end subroutine UNUSED_compute_gradient_field_element

