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
! for poroelastic solver: update memory variables with fourth-order Runge-Kutta time scheme for attenuation

 subroutine compute_attenuation_poro_fluid_part()

  use specfem_par, only: nspec,poroelastic,poroelastcoef,kmato,permeability,ibool,viscox_loc,viscoz_loc, &
                         velocw_poroelastic,time_stepping_scheme,deltat,i_stage,stage_time_scheme, &
                         rx_viscous,rz_viscous,viscox,viscoz, &
                         rx_viscous_force_RK,rx_viscous_initial_rk,rz_viscous_force_RK,rz_viscous_initial_rk, &
                         rx_viscous_LDDRK,rz_viscous_LDDRK,alpha_LDDRK,beta_LDDRK, &
                         alphaval,betaval,gammaval,theta_e,theta_s,thetainv

  implicit none
  include "constants.h"

  ! local variables
  integer :: i,j,ispec,iglob
  double precision :: etal_f,permlxx,permlxz,permlzz,detk,invpermlxx,invpermlxz,invpermlzz, &
                      Sn,Snp1,weight_rk
  double precision, dimension(3) :: bl_unrelaxed_elastic

  ! loop over spectral elements
  do ispec = 1,nspec
    if( poroelastic(ispec) .and. poroelastcoef(2,2,kmato(ispec)) > 0.d0 ) then
      etal_f = poroelastcoef(2,2,kmato(ispec))
      permlxx = permeability(1,kmato(ispec))
      permlxz = permeability(2,kmato(ispec))
      permlzz = permeability(3,kmato(ispec))

      ! calcul of the inverse of k
      detk = permlxx * permlzz - permlxz * permlxz
      if( detk /= ZERO ) then
        invpermlxx = permlzz/detk
        invpermlxz = -permlxz/detk
        invpermlzz = permlxx/detk
      else
        stop 'Permeability matrix is not invertible'
      endif

      ! relaxed viscous coef
      bl_unrelaxed_elastic(1) = etal_f*invpermlxx
      bl_unrelaxed_elastic(2) = etal_f*invpermlxz
      bl_unrelaxed_elastic(3) = etal_f*invpermlzz

      do j=1,NGLLZ
        do i=1,NGLLX
          iglob = ibool(i,j,ispec)
          viscox_loc(i,j) = velocw_poroelastic(1,iglob) * bl_unrelaxed_elastic(1) + &
                            velocw_poroelastic(2,iglob) * bl_unrelaxed_elastic(2)
          viscoz_loc(i,j) = velocw_poroelastic(1,iglob) * bl_unrelaxed_elastic(2) + &
                            velocw_poroelastic(2,iglob) * bl_unrelaxed_elastic(3)

          if( time_stepping_scheme == 1 ) then
            ! evolution rx_viscous
            Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscox(i,j,ispec)
            Snp1 = - (1.d0 - theta_e/theta_s)/theta_s*viscox_loc(i,j)
            rx_viscous(i,j,ispec) = alphaval * rx_viscous(i,j,ispec) + betaval * Sn + gammaval * Snp1

            ! evolution rz_viscous
            Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscoz(i,j,ispec)
            Snp1 = - (1.d0 - theta_e/theta_s)/theta_s*viscoz_loc(i,j)
            rz_viscous(i,j,ispec) = alphaval * rz_viscous(i,j,ispec) + betaval * Sn + gammaval * Snp1
          endif

          if( time_stepping_scheme == 2 ) then
            Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscox(i,j,ispec)
            rx_viscous_LDDRK(i,j,ispec) = alpha_LDDRK(i_stage) * rx_viscous_LDDRK(i,j,ispec) + &
                                          deltat * (Sn + thetainv * rx_viscous(i,j,ispec))
            rx_viscous(i,j,ispec)= rx_viscous(i,j,ispec)+beta_LDDRK(i_stage) * rx_viscous_LDDRK(i,j,ispec)

            Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscoz(i,j,ispec)
            rz_viscous_LDDRK(i,j,ispec)= alpha_LDDRK(i_stage) * rz_viscous_LDDRK(i,j,ispec)+&
                                         deltat * (Sn + thetainv * rz_viscous(i,j,ispec))
            rz_viscous(i,j,ispec)= rz_viscous(i,j,ispec)+beta_LDDRK(i_stage) * rz_viscous_LDDRK(i,j,ispec)
          endif

          if(time_stepping_scheme == 3) then

            Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscox(i,j,ispec)
            rx_viscous_force_RK(i,j,ispec,i_stage) = deltat * (Sn + thetainv * rx_viscous(i,j,ispec))

            if( i_stage==1 .or. i_stage==2 .or. i_stage==3 ) then
              if( i_stage == 1 )weight_rk = 0.5d0
              if( i_stage == 2 )weight_rk = 0.5d0
              if( i_stage == 3 )weight_rk = 1.0d0

              if( i_stage==1 ) then
                rx_viscous_initial_rk(i,j,ispec) = rx_viscous(i,j,ispec)
              endif
                  rx_viscous(i,j,ispec) = rx_viscous_initial_rk(i,j,ispec) + &
                                          weight_rk * rx_viscous_force_RK(i,j,ispec,i_stage)
            else if( i_stage==4 ) then
                rx_viscous(i,j,ispec) = rx_viscous_initial_rk(i,j,ispec) + &
                                        1.0d0 / 6.0d0 * (rx_viscous_force_RK(i,j,ispec,i_stage) + &
                                        2.0d0 * rx_viscous_force_RK(i,j,ispec,i_stage) + &
                                        2.0d0 * rx_viscous_force_RK(i,j,ispec,i_stage) + &
                                        rx_viscous_force_RK(i,j,ispec,i_stage))
            endif

            Sn   = - (1.d0 - theta_e/theta_s)/theta_s*viscoz(i,j,ispec)
            rz_viscous_force_RK(i,j,ispec,i_stage) = deltat * (Sn + thetainv * rz_viscous(i,j,ispec))

            if( i_stage==1 .or. i_stage==2 .or. i_stage==3 ) then
              if( i_stage == 1 )weight_rk = 0.5d0
              if( i_stage == 2 )weight_rk = 0.5d0
              if( i_stage == 3 )weight_rk = 1.0d0

              if( i_stage==1 ) then
                rz_viscous_initial_rk(i,j,ispec) = rz_viscous(i,j,ispec)
              endif
              rz_viscous(i,j,ispec) = rz_viscous_initial_rk(i,j,ispec) + &
                                      weight_rk * rz_viscous_force_RK(i,j,ispec,i_stage)
            else if(i_stage==4) then
              rz_viscous(i,j,ispec) = rz_viscous_initial_rk(i,j,ispec) + &
                                      1.0d0 / 6.0d0 * (rz_viscous_force_RK(i,j,ispec,i_stage) + &
                                      2.0d0 * rz_viscous_force_RK(i,j,ispec,i_stage) + &
                                      2.0d0 * rz_viscous_force_RK(i,j,ispec,i_stage) + &
                                      rz_viscous_force_RK(i,j,ispec,i_stage))
            endif
          endif
        enddo
      enddo

      if( stage_time_scheme == 1 ) then
        ! save visco for Runge-Kutta scheme when used together with Newmark
        viscox(:,:,ispec) = viscox_loc(:,:)
        viscoz(:,:,ispec) = viscoz_loc(:,:)
      endif

    endif  ! end of poroelastic element loop
  enddo   ! end of spectral element loop

 end subroutine compute_attenuation_poro_fluid_part
