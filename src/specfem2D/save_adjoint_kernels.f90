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

  subroutine save_adjoint_kernels()

! saves adjoint sensitivity kernels to file

  use constants, only: NGLLX,NGLLZ,IMAIN,APPROXIMATE_HESS_KL

  use specfem_par, only: myrank, nspec, ibool, coord, save_ASCII_kernels, &
                         any_acoustic, any_elastic, any_poroelastic, &
                         rho_ac_kl, kappa_ac_kl, alpha_ac_kl, rhop_ac_kl, &
                         rho_kl, kappa_kl, mu_kl, rhop_kl, alpha_kl, beta_kl, &
                         bulk_c_kl, bulk_beta_kl, &
                         rhorho_ac_Hessian_final1, rhorho_ac_Hessian_final2, &
                         rhorho_el_Hessian_final1, rhorho_el_Hessian_final2, &
                         rhot_kl, rhof_kl, sm_kl, eta_kl, mufr_kl, B_kl, &
                         C_kl, M_kl, rhob_kl, rhofb_kl, phi_kl, mufrb_kl, &
                         rhobb_kl, rhofbb_kl, phib_kl, cpI_kl, cpII_kl, cs_kl, ratio_kl, GPU_MODE

  use specfem_par, only: ispec_is_anisotropic, c11_kl, c13_kl, c15_kl, c33_kl, c35_kl, c55_kl

  implicit none

  integer :: i, j, ispec, iglob
  double precision :: xx, zz

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Writing Kernels file'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (any_acoustic) then
    ! acoustic domain
    if (save_ASCII_kernels) then
      ! ascii format
      do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            xx = coord(1,iglob)
            zz = coord(2,iglob)
            write(95,'(4e15.5e4)') xx,zz,rho_ac_kl(i,j,ispec),kappa_ac_kl(i,j,ispec)
            write(96,'(4e15.5e4)') xx,zz,rhop_ac_kl(i,j,ispec),alpha_ac_kl(i,j,ispec)
          enddo
        enddo
      enddo
      close(95)
      close(96)
    else
      ! binary format
      write(200) rho_ac_kl
      write(201) kappa_ac_kl
      write(202) rhop_ac_kl
      write(203) alpha_ac_kl
      close(200)
      close(201)
      close(202)
      close(203)
      ! Hessian kernels
      if (APPROXIMATE_HESS_KL) then
        write(212) rhorho_ac_Hessian_final1
        write(213) rhorho_ac_Hessian_final2
        close(212)
        close(213)
      endif
    endif
  endif

  if (any_elastic) then
    ! elastic domains
    if (save_ASCII_kernels) then
      ! ascii format
      do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            xx = coord(1,iglob)
            zz = coord(2,iglob)
            if (count(ispec_is_anisotropic(:) .eqv. .true.) >= 1) then
              ! anisotropic
              write(97,'(9e15.5e4)') xx, zz, rho_kl(i,j,ispec), c11_kl(i,j,ispec), &
                                     c13_kl(i,j,ispec), c15_kl(i,j,ispec), c33_kl(i,j,ispec), c35_kl(i,j,ispec), &
                                     c55_kl(i,j,ispec)
            else
              ! isotropic
              ! parameterization (rho,kappa,mu) "primary" kernels
              write(97,'(5e15.5e4)') xx,zz,rho_kl(i,j,ispec),kappa_kl(i,j,ispec),mu_kl(i,j,ispec)
              ! parameterization (rho,Vp,Vs)
              write(98,'(5e15.5e4)') xx,zz,rhop_kl(i,j,ispec),alpha_kl(i,j,ispec),beta_kl(i,j,ispec)
            endif
          enddo
        enddo
      enddo
      close(97)
      close(98)
    else
      ! binary format
      write(204) rho_kl
      if (count(ispec_is_anisotropic(:) .eqv. .true.) >= 1) then ! anisotropic
        write(205) c11_kl
        write(206) c13_kl
        write(207) c15_kl
        write(208) c33_kl
        write(209) c35_kl
        write(210) c55_kl
      else
        write(205) kappa_kl
        write(206) mu_kl
        write(207) rhop_kl
        write(208) alpha_kl
        write(209) beta_kl
        write(210) bulk_c_kl
        write(211) bulk_beta_kl
        close(211)
      endif
      close(204)
      close(205)
      close(202)
      close(207)
      close(208)
      close(209)
      close(210)
      ! Hessian kernels
      if (APPROXIMATE_HESS_KL) then
        write(214) rhorho_el_Hessian_final1
        write(215) rhorho_el_Hessian_final2
        close(214)
        close(215)
      endif
    endif
  endif

  if (any_poroelastic) then
    ! poro-elastic domains
    ! checks
    if (GPU_MODE) call stop_the_code('poroelastic kernel output not implemented on GPUs yet')
    if (.not. SAVE_ASCII_KERNELS) call stop_the_code('poroelastic simulations must use SAVE_ASCII_KERNELS')
    ! ascii format
    do ispec = 1, nspec
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          xx = coord(1,iglob)
          zz = coord(2,iglob)
          write(144,'(5e11.3)') xx,zz,mufr_kl(i,j,ispec),B_kl(i,j,ispec),C_kl(i,j,ispec)
          write(155,'(5e11.3)') xx,zz,M_kl(i,j,ispec),rhot_kl(i,j,ispec),rhof_kl(i,j,ispec)

          write(16,'(5e11.3)') xx,zz,sm_kl(i,j,ispec),eta_kl(i,j,ispec)

          write(17,'(5e11.3)') xx,zz,mufrb_kl(i,j,ispec),B_kl(i,j,ispec),C_kl(i,j,ispec)
          write(18,'(5e11.3)') xx,zz,M_kl(i,j,ispec),rhob_kl(i,j,ispec),rhofb_kl(i,j,ispec)

          write(19,'(5e11.3)') xx,zz,phi_kl(i,j,ispec),eta_kl(i,j,ispec)
          write(20,'(5e11.3)') xx,zz,cpI_kl(i,j,ispec),cpII_kl(i,j,ispec),cs_kl(i,j,ispec)
          write(21,'(5e11.3)') xx,zz,rhobb_kl(i,j,ispec),rhofbb_kl(i,j,ispec),ratio_kl(i,j,ispec)
          write(22,'(5e11.3)') xx,zz,phib_kl(i,j,ispec),eta_kl(i,j,ispec)
        enddo
      enddo
    enddo
    close(144)
    close(155)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)
  endif

  end subroutine save_adjoint_kernels

