
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
!=====================================================================
!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_el()

! elastic kernel calculations
! see e.g. Tromp et al. (2005)

  use specfem_par, only: elastic,rho_k,rhorho_el_hessian_temp1,rhorho_el_hessian_temp2, &
                         rho_kl,mu_kl,kappa_kl,rhop_kl,beta_kl,alpha_kl,bulk_c_kl,bulk_beta_kl, &
                         rhorho_el_hessian_final1,rhorho_el_hessian_final2, &
                         nglob,nspec,ibool,accel_elastic,b_displ_elastic,b_accel_elastic, &
                         rhol_global,mul_global,kappal_global, &
                         density,poroelastcoef,kmato,assign_external_model,rhoext,vsext,vpext,&
                         deltat,p_sv,displ_elastic,&
                         mu_k,kappa_k,elastic,ibool,hprime_xx,hprime_zz,xix,xiz,gammax,gammaz
  implicit none
  include "constants.h"

  !local variables
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duy_dxi,duy_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duy_dxl,duz_dxl,dux_dzl,duy_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: b_dux_dxi,b_dux_dgamma,b_duy_dxi,b_duy_dgamma,b_duz_dxi,b_duz_dgamma
  real(kind=CUSTOM_REAL) :: b_dux_dxl,b_duy_dxl,b_duz_dxl,b_dux_dzl,b_duy_dzl,b_duz_dzl
  real(kind=CUSTOM_REAL) :: dsxx,dsxz,dszz
  real(kind=CUSTOM_REAL) :: b_dsxx,b_dsxz,b_dszz

  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl

  do ispec = 1,nspec
    if( elastic(ispec) ) then
      do j=1,NGLLZ; do i=1,NGLLX
        ! derivative along x and along z
        dux_dxi = 0._CUSTOM_REAL; duy_dxi = 0._CUSTOM_REAL; duz_dxi = 0._CUSTOM_REAL
        dux_dgamma = 0._CUSTOM_REAL; duy_dgamma = 0._CUSTOM_REAL; duz_dgamma = 0._CUSTOM_REAL
        b_dux_dxi = 0._CUSTOM_REAL; b_duy_dxi = 0._CUSTOM_REAL; b_duz_dxi = 0._CUSTOM_REAL
        b_dux_dgamma = 0._CUSTOM_REAL; b_duy_dgamma = 0._CUSTOM_REAL; b_duz_dgamma = 0._CUSTOM_REAL

        ! first double loop over GLL points to compute and store gradients
        ! we can merge the two loops because NGLLX == NGLLZ
        do k = 1,NGLLX
          dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          duy_dxi = duy_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
          dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
          duy_dgamma = duy_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
          duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)


          b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          b_duy_dxi = b_duy_dxi + b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          b_duz_dxi = b_duz_dxi + b_displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
          b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
          b_duy_dgamma = b_duy_dgamma + b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
          b_duz_dgamma = b_duz_dgamma + b_displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)

        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! derivatives of displacement
        dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
        dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl
        duy_dxl = duy_dxi*xixl + duy_dgamma*gammaxl
        duy_dzl = duy_dxi*xizl + duy_dgamma*gammazl
        duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
        duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl


        b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
        b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

        b_duy_dxl = b_duy_dxi*xixl + b_duy_dgamma*gammaxl
        b_duy_dzl = b_duy_dxi*xizl + b_duy_dgamma*gammazl

        b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
        b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl

        iglob = ibool(i,j,ispec)
        if( p_sv ) then !P-SV waves
          dsxx =  dux_dxl
          dsxz = HALF * (duz_dxl + dux_dzl)
          dszz =  duz_dzl

          b_dsxx =  b_dux_dxl
          b_dsxz = HALF * (b_duz_dxl + b_dux_dzl)
          b_dszz =  b_duz_dzl

          kappa_k(iglob) = (dux_dxl + duz_dzl) *  (b_dux_dxl + b_duz_dzl)
          mu_k(iglob) = dsxx * b_dsxx + dszz * b_dszz + &
                        2._CUSTOM_REAL * dsxz * b_dsxz - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k(iglob)
        else !SH (membrane) waves
          mu_k(iglob) = duy_dxl * b_duy_dxl + duy_dzl * b_duy_dzl
        endif
      enddo; enddo
    endif
  enddo

  do iglob = 1,nglob
    rho_k(iglob) =  accel_elastic(1,iglob)*b_displ_elastic(1,iglob) + &
                    accel_elastic(2,iglob)*b_displ_elastic(2,iglob) + &
                    accel_elastic(3,iglob)*b_displ_elastic(3,iglob)
    rhorho_el_hessian_temp1(iglob) = b_accel_elastic(1,iglob)*b_accel_elastic(1,iglob) + &
                                     b_accel_elastic(2,iglob)*b_accel_elastic(2,iglob) + &
                                     b_accel_elastic(3,iglob)*b_accel_elastic(3,iglob)
    rhorho_el_hessian_temp2(iglob) = accel_elastic(1,iglob)*b_accel_elastic(1,iglob) + &
                                     accel_elastic(2,iglob)*b_accel_elastic(2,iglob) + &
                                     accel_elastic(3,iglob)*b_accel_elastic(3,iglob)
  enddo

  do ispec = 1, nspec
    if( elastic(ispec) ) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          if( .not. assign_external_model ) then
            rhol_global(iglob) = density(1,kmato(ispec))
            mul_global(iglob) = poroelastcoef(2,1,kmato(ispec))
            kappal_global(iglob) = poroelastcoef(3,1,kmato(ispec)) - &
                                   4._CUSTOM_REAL*mul_global(iglob) / 3._CUSTOM_REAL
          else
            rhol_global(iglob)   = rhoext(i,j,ispec)
            mul_global(iglob)    = rhoext(i,j,ispec)*vsext(i,j,ispec)*vsext(i,j,ispec)
            kappal_global(iglob) = rhoext(i,j,ispec)*vpext(i,j,ispec)*vpext(i,j,ispec) - &
                                   4._CUSTOM_REAL*mul_global(iglob) / 3._CUSTOM_REAL
          endif

          rho_kl(i,j,ispec) = rho_kl(i,j,ispec) - rhol_global(iglob)  * rho_k(iglob) * deltat
          mu_kl(i,j,ispec) =  mu_kl(i,j,ispec) - TWO * mul_global(iglob) * mu_k(iglob) * deltat
          kappa_kl(i,j,ispec) = kappa_kl(i,j,ispec) - kappal_global(iglob) * kappa_k(iglob) * deltat
          !
          rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + kappa_kl(i,j,ispec) + mu_kl(i,j,ispec)
          beta_kl(i,j,ispec) = TWO * (mu_kl(i,j,ispec) - 4._CUSTOM_REAL * mul_global(iglob)/&
                        (3._CUSTOM_REAL * kappal_global(iglob)) * kappa_kl(i,j,ispec))
          alpha_kl(i,j,ispec) = TWO * (1._CUSTOM_REAL + 4._CUSTOM_REAL * mul_global(iglob)/&
                         (3._CUSTOM_REAL * kappal_global(iglob))) * kappa_kl(i,j,ispec)
          !
          bulk_c_kl(i,j,ispec) =  TWO * kappa_kl(i,j,ispec)
          bulk_beta_kl(i,j,ispec) =  TWO * mu_kl(i,j,ispec)

          rhorho_el_hessian_final1(i,j,ispec) = rhorho_el_hessian_final1(i,j,ispec) + &
                                  rhorho_el_hessian_temp1(iglob) * deltat
          rhorho_el_hessian_final2(i,j,ispec) = rhorho_el_hessian_final2(i,j,ispec) + &
                                  rhorho_el_hessian_temp2(iglob) * deltat

        enddo
      enddo
    endif
  enddo

  end subroutine compute_kernels_el

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_ac()

! acoustic kernel calculations
! see e.g. Tromp et al. (2005)

  use specfem_par, only: nspec,acoustic,ibool,kappal_ac_global,rhol_ac_global,&
                         poroelastcoef,density,kmato,assign_external_model,rhoext,vpext,deltat,&
                         hprime_xx,hprime_zz,xix,xiz,gammax,gammaz,&
                         potential_acoustic,b_potential_acoustic,b_potential_dot_dot_acoustic,&
                         accel_ac,b_accel_ac,b_displ_ac,&
                         rho_ac_kl,kappa_ac_kl,rhop_ac_kl,alpha_ac_kl,rhorho_ac_hessian_final1,&
                         rhorho_ac_hessian_final2
  implicit none
  include "constants.h"

  !local variables
  integer :: i,j,k,ispec,iglob
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,b_tempx1l,b_tempx2l,bb_tempx1l,bb_tempx2l
  double precision :: xixl,xizl,gammaxl,gammazl

  do ispec = 1, nspec
    if( acoustic(ispec) ) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          if( .not. assign_external_model ) then
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
          bb_tempx1l = ZERO
          bb_tempx2l = ZERO
          do k = 1,NGLLX
            ! derivative along x
            !tempx1l = tempx1l + potential_dot_dot_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
            tempx1l = tempx1l + potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k) !!! YANGL
            b_tempx1l = b_tempx1l + b_potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
            bb_tempx1l = bb_tempx1l + b_potential_dot_dot_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
            ! derivative along z
            !tempx2l = tempx2l + potential_dot_dot_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
            tempx2l = tempx2l + potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k) !!! YANGL
            b_tempx2l = b_tempx2l + b_potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
            bb_tempx2l = bb_tempx2l + b_potential_dot_dot_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
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
          b_accel_ac(1,iglob) = (bb_tempx1l*xixl + bb_tempx2l*gammaxl) / rhol_ac_global(iglob)
          b_accel_ac(2,iglob) = (bb_tempx1l*xizl + bb_tempx2l*gammazl) / rhol_ac_global(iglob)
        enddo !i = 1, NGLLX
      enddo !j = 1, NGLLZ
    endif
  enddo

  do ispec = 1,nspec
    if( acoustic(ispec) ) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          !<YANGL
          !!!! old expression (from elastic kernels)
          !!!rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) - rhol_ac_global(iglob)  * &
          !!!      dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * deltat
          !!!kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) - kappal_ac_global(iglob) * &
          !!!      potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * &
          !!!      b_potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob)&
          !!!      * deltat
          !!!! new expression (from PDE-constrained optimization, coupling terms changed as well)
          rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + rhol_ac_global(iglob) * &
                                 dot_product(accel_ac(:,iglob),b_displ_ac(:,iglob)) * deltat
          kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) + kappal_ac_global(iglob) * &
                                   potential_acoustic(iglob)/kappal_ac_global(iglob) * &
                                   b_potential_dot_dot_acoustic(iglob)/kappal_ac_global(iglob) * deltat
          !>YANGL
          rhop_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) + kappa_ac_kl(i,j,ispec)
          alpha_ac_kl(i,j,ispec) = TWO *  kappa_ac_kl(i,j,ispec)
          rhorho_ac_hessian_final1(i,j,ispec) =  rhorho_ac_hessian_final1(i,j,ispec) + &
                                                 dot_product(accel_ac(:,iglob),accel_ac(:,iglob)) * deltat
          rhorho_ac_hessian_final2(i,j,ispec) =  rhorho_ac_hessian_final2(i,j,ispec) + &
                                                 dot_product(accel_ac(:,iglob),b_accel_ac(:,iglob)) * deltat
        enddo
      enddo
    endif
  enddo

  end subroutine compute_kernels_ac

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_po()

! kernel calculations
! see e.g. Morency et al. (2009)

  use specfem_par, only: nglob,nspec,poroelastic,ibool,deltat, &
                         kmato,porosity,tortuosity,density,permeability,poroelastcoef, &
                         ratio,B_biot,M_biot,C_biot,cpIsquare,cpIIsquare,cssquare, &
                         accels_poroelastic,accelw_poroelastic,velocw_poroelastic, &
                         b_displs_poroelastic,b_displw_poroelastic, &
                         rhot_k,rhof_k,sm_k,eta_k,B_k,C_k, &
                         phil_global,tortl_global,rhol_s_global,rhol_f_global,rhol_bar_global, &
                         etal_f_global,permlxx_global,permlxz_global,permlzz_global,mulfr_global, &
                         rhot_kl,rhof_kl,sm_kl,eta_kl,B_kl,C_kl,M_kl,M_k, &
                         mufr_kl,mufr_k,rhol_bar_global,rhob_kl,rhofb_kl,Bb_kl,Cb_kl,Mb_kl, &
                         mufrb_kl,phi_kl,rhobb_kl,rhofbb_kl,phib_kl,cpI_kl,cpII_kl,cs_kl,ratio_kl
  implicit none
  include "constants.h"

  !local variables
  integer :: i,j,ispec,iglob
  real(kind=CUSTOM_REAL) :: rholb,dd1

  do iglob =1,nglob
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

  do ispec = 1, nspec
    if( poroelastic(ispec) ) then
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          phil_global(iglob) = porosity(kmato(ispec))
          tortl_global(iglob) = tortuosity(kmato(ispec))
          rhol_s_global(iglob) = density(1,kmato(ispec))
          rhol_f_global(iglob) = density(2,kmato(ispec))
          rhol_bar_global(iglob) =  (1._CUSTOM_REAL - phil_global(iglob))*rhol_s_global(iglob) + &
                                    phil_global(iglob)*rhol_f_global(iglob)
          etal_f_global(iglob) = poroelastcoef(2,2,kmato(ispec))
          permlxx_global(iglob) = permeability(1,kmato(ispec))
          permlxz_global(iglob) = permeability(2,kmato(ispec))
          permlzz_global(iglob) = permeability(3,kmato(ispec))
          mulfr_global(iglob) = poroelastcoef(2,3,kmato(ispec))

          rhot_kl(i,j,ispec) = rhot_kl(i,j,ispec) - deltat * rhol_bar_global(iglob) * rhot_k(iglob)
          rhof_kl(i,j,ispec) = rhof_kl(i,j,ispec) - deltat * rhol_f_global(iglob) * rhof_k(iglob)
          sm_kl(i,j,ispec) = sm_kl(i,j,ispec) - &
                             deltat * rhol_f_global(iglob)*tortl_global(iglob)/phil_global(iglob) * sm_k(iglob)
          !at the moment works with constant permeability
          eta_kl(i,j,ispec) = eta_kl(i,j,ispec) - deltat * etal_f_global(iglob)/permlxx_global(iglob) * eta_k(iglob)
          B_kl(i,j,ispec) = B_kl(i,j,ispec) - deltat * B_k(iglob)
          C_kl(i,j,ispec) = C_kl(i,j,ispec) - deltat * C_k(iglob)
          M_kl(i,j,ispec) = M_kl(i,j,ispec) - deltat * M_k(iglob)
          mufr_kl(i,j,ispec) = mufr_kl(i,j,ispec) - TWO * deltat * mufr_k(iglob)
          ! density kernels
          rholb = rhol_bar_global(iglob) - phil_global(iglob)*rhol_f_global(iglob)/tortl_global(iglob)
          rhob_kl(i,j,ispec) = rhot_kl(i,j,ispec) + B_kl(i,j,ispec) + mufr_kl(i,j,ispec)
          rhofb_kl(i,j,ispec) = rhof_kl(i,j,ispec) + C_kl(i,j,ispec) + M_kl(i,j,ispec) + sm_kl(i,j,ispec)
          Bb_kl(i,j,ispec) = B_kl(i,j,ispec)
          Cb_kl(i,j,ispec) = C_kl(i,j,ispec)
          Mb_kl(i,j,ispec) = M_kl(i,j,ispec)
          mufrb_kl(i,j,ispec) = mufr_kl(i,j,ispec)
          phi_kl(i,j,ispec) = - sm_kl(i,j,ispec) - M_kl(i,j,ispec)
          ! wave speed kernels
          dd1 = (1._CUSTOM_REAL+rholb/rhol_f_global(iglob))*ratio**2 + &
                2._CUSTOM_REAL*ratio + tortl_global(iglob)/phil_global(iglob)

          rhobb_kl(i,j,ispec) = rhob_kl(i,j,ispec) - &
                phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*B_biot) * &
                (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil_global(iglob) / &
                tortl_global(iglob)*ratio +1._CUSTOM_REAL)/dd1 + &
                (rhol_bar_global(iglob)**2*ratio**2/rhol_f_global(iglob)**2*(phil_global(iglob) / &
                tortl_global(iglob)*ratio+1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*ratio + &
                phil_global(iglob)/tortl_global(iglob) * &
                (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-1._CUSTOM_REAL) )/dd1**2 ) - &
                FOUR_THIRDS*cssquare )*Bb_kl(i,j,ispec) - &
                rhol_bar_global(iglob)*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                (phil_global(iglob)/tortl_global(iglob)*ratio + &
                1._CUSTOM_REAL)**2/dd1**2*Mb_kl(i,j,ispec) + &
                rhol_bar_global(iglob)*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                phil_global(iglob)*ratio/tortl_global(iglob)*(phil_global(iglob) / &
                tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                (1._CUSTOM_REAL+rhol_bar_global(iglob)*ratio/rhol_f_global(iglob))/dd1**2)*Cb_kl(i,j,ispec)+ &
                phil_global(iglob)*rhol_f_global(iglob)*cssquare / &
                (tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
          rhofbb_kl(i,j,ispec) = rhofb_kl(i,j,ispec) + &
                 phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*B_biot) * &
                 (cpIIsquare + (cpIsquare - cpIIsquare)*( (phil_global(iglob)/ &
                 tortl_global(iglob)*ratio +1._CUSTOM_REAL)/dd1+&
                 (rhol_bar_global(iglob)**2*ratio**2/rhol_f_global(iglob)**2*(phil_global(iglob)/ &
                 tortl_global(iglob)*ratio+1)*(phil_global(iglob)/tortl_global(iglob)*ratio+ &
                 phil_global(iglob)/tortl_global(iglob)*&
                 (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-1._CUSTOM_REAL) )/dd1**2 )- &
                 FOUR_THIRDS*cssquare )*Bb_kl(i,j,ispec) + &
                 rhol_bar_global(iglob)*ratio**2/M_biot * (cpIsquare - cpIIsquare)* &
                 (phil_global(iglob)/tortl_global(iglob)*ratio + &
                 1._CUSTOM_REAL)**2/dd1**2*Mb_kl(i,j,ispec) - &
                 rhol_bar_global(iglob)*ratio/C_biot * (cpIsquare - cpIIsquare)* (&
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                 phil_global(iglob)*ratio/tortl_global(iglob)*(phil_global(iglob)/ &
                 tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                 (1._CUSTOM_REAL+rhol_bar_global(iglob)*ratio/rhol_f_global(iglob))/dd1**2)*Cb_kl(i,j,ispec)- &
                 phil_global(iglob)*rhol_f_global(iglob)*cssquare/ &
                 (tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
          phib_kl(i,j,ispec) = phi_kl(i,j,ispec) - &
                 phil_global(iglob)*rhol_bar_global(iglob)/(tortl_global(iglob)*B_biot) &
                 * ( cpIsquare - rhol_f_global(iglob)/rhol_bar_global(iglob)*cpIIsquare- &
                 (cpIsquare-cpIIsquare)*( (TWO*ratio**2*phil_global(iglob)/ &
                 tortl_global(iglob) + (1._CUSTOM_REAL+&
                 rhol_f_global(iglob)/rhol_bar_global(iglob))* &
                 (TWO*ratio*phil_global(iglob)/tortl_global(iglob)+&
                 1._CUSTOM_REAL))/dd1 + (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                 1._CUSTOM_REAL)*(phil_global(iglob)*&
                 ratio/tortl_global(iglob)+phil_global(iglob)/tortl_global(iglob)* &
                 (1._CUSTOM_REAL+rhol_f_global(iglob)/&
                 rhol_bar_global(iglob))-1._CUSTOM_REAL)*((1._CUSTOM_REAL+ &
                 rhol_bar_global(iglob)/rhol_f_global(iglob)-&
                 TWO*phil_global(iglob)/tortl_global(iglob))*ratio**2+TWO*ratio)/dd1**2 ) - &
                 FOUR_THIRDS*rhol_f_global(iglob)*cssquare/rhol_bar_global(iglob) )*Bb_kl(i,j,ispec) + &
                 rhol_f_global(iglob)/M_biot * (cpIsquare-cpIIsquare)*(&
                 TWO*ratio*(phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2*( &
                 (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                 rhol_f_global(iglob)-TWO*phil_global(iglob)/tortl_global(iglob))*ratio**2+TWO*ratio)/dd1**2 &
                 )*Mb_kl(i,j,ispec) + &
                 phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*C_biot)* &
                 (cpIsquare-cpIIsquare)*ratio* (&
                 (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob)*ratio)/dd1 - &
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)* &
                 (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                 rhol_f_global(iglob)*ratio)*((1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-TWO*&
                 phil_global(iglob)/tortl_global(iglob))*ratio+TWO)/dd1**2&
                 )*Cb_kl(i,j,ispec) -&
                 phil_global(iglob)*rhol_f_global(iglob)*cssquare &
                 /(tortl_global(iglob)*mulfr_global(iglob))*mufrb_kl(i,j,ispec)
          cpI_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIsquare/B_biot*rhol_bar_global(iglob)*( &
                 1._CUSTOM_REAL-phil_global(iglob)/tortl_global(iglob) + &
                 (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                 1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*&
                 ratio+phil_global(iglob)/tortl_global(iglob)* &
                 (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-&
                 1._CUSTOM_REAL)/dd1 &
                 )* Bb_kl(i,j,ispec) +&
                 2._CUSTOM_REAL*cpIsquare*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot) *&
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2/dd1*Mb_kl(i,j,ispec)+&
                 2._CUSTOM_REAL*cpIsquare*rhol_f_global(iglob)/C_biot * &
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)* &
                 (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                 rhol_f_global(iglob)*ratio)/dd1*Cb_kl(i,j,ispec)
          cpII_kl(i,j,ispec) = 2._CUSTOM_REAL*cpIIsquare*rhol_bar_global(iglob)/B_biot * (&
                 phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)*rhol_bar_global(iglob)) - &
                 (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                 1._CUSTOM_REAL)*(phil_global(iglob)/tortl_global(iglob)*&
                 ratio+phil_global(iglob)/tortl_global(iglob)* &
                 (1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-&
                 1._CUSTOM_REAL)/dd1  ) * Bb_kl(i,j,ispec) +&
                 2._CUSTOM_REAL*cpIIsquare*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot) * (&
                 1._CUSTOM_REAL - (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                 1._CUSTOM_REAL)**2/dd1  )*Mb_kl(i,j,ispec) + &
                 2._CUSTOM_REAL*cpIIsquare*rhol_f_global(iglob)/C_biot * (&
                 1._CUSTOM_REAL - (phil_global(iglob)/tortl_global(iglob)*ratio+ &
                 1._CUSTOM_REAL)*(1._CUSTOM_REAL+&
                 rhol_bar_global(iglob)/rhol_f_global(iglob)*ratio)/dd1  )*Cb_kl(i,j,ispec)
          cs_kl(i,j,ispec) = - 8._CUSTOM_REAL/3._CUSTOM_REAL*cssquare* &
                 rhol_bar_global(iglob)/B_biot*(1._CUSTOM_REAL-&
                 phil_global(iglob)*rhol_f_global(iglob)/(tortl_global(iglob)* &
                 rhol_bar_global(iglob)))*Bb_kl(i,j,ispec) + &
                 2._CUSTOM_REAL*(rhol_bar_global(iglob)-rhol_f_global(iglob)*&
                 phil_global(iglob)/tortl_global(iglob))/&
                 mulfr_global(iglob)*cssquare*mufrb_kl(i,j,ispec)
          ratio_kl(i,j,ispec) = ratio*rhol_bar_global(iglob)*phil_global(iglob)/(tortl_global(iglob)*B_biot) * &
                 (cpIsquare-cpIIsquare) * ( &
                 phil_global(iglob)/tortl_global(iglob)*(2._CUSTOM_REAL*ratio+1._CUSTOM_REAL+rhol_f_global(iglob)/ &
                 rhol_bar_global(iglob))/dd1 - (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)*&
                 (phil_global(iglob)/tortl_global(iglob)*ratio+phil_global(iglob)/tortl_global(iglob)*(&
                 1._CUSTOM_REAL+rhol_f_global(iglob)/rhol_bar_global(iglob))-1._CUSTOM_REAL)*(2._CUSTOM_REAL*ratio*(&
                 1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob)) +&
                 2._CUSTOM_REAL)/dd1**2  )*Bb_kl(i,j,ispec) + &
                 ratio*rhol_f_global(iglob)*tortl_global(iglob)/(phil_global(iglob)*M_biot)*(cpIsquare-cpIIsquare) * &
                 2._CUSTOM_REAL*phil_global(iglob)/tortl_global(iglob) * (&
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)/dd1 - &
                 (phil_global(iglob)/tortl_global(iglob)*ratio+1._CUSTOM_REAL)**2*( &
                 (1._CUSTOM_REAL+rhol_bar_global(iglob)/&
                 rhol_f_global(iglob)-phil_global(iglob)/tortl_global(iglob))*ratio+ &
                 1._CUSTOM_REAL)/dd1**2 )*Mb_kl(i,j,ispec) +&
                 ratio*rhol_f_global(iglob)/C_biot*(cpIsquare-cpIIsquare) * (&
                 (2._CUSTOM_REAL*phil_global(iglob)*rhol_bar_global(iglob)* &
                 ratio/(tortl_global(iglob)*rhol_f_global(iglob))+&
                 phil_global(iglob)/tortl_global(iglob)+rhol_bar_global(iglob)/rhol_f_global(iglob))/dd1 - &
                 2._CUSTOM_REAL*phil_global(iglob)/tortl_global(iglob)*(phil_global(iglob)/tortl_global(iglob)*ratio+&
                 1._CUSTOM_REAL)*(1._CUSTOM_REAL+rhol_bar_global(iglob)/rhol_f_global(iglob)*ratio)*((1._CUSTOM_REAL+&
                 rhol_bar_global(iglob)/rhol_f_global(iglob)- &
                 phil_global(iglob)/tortl_global(iglob))*ratio+1._CUSTOM_REAL)/&
                 dd1**2 )*Cb_kl(i,j,ispec)
        enddo
      enddo
    endif
  enddo

  end subroutine compute_kernels_po

