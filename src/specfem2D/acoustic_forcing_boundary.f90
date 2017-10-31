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

! define the forcing applied at the bottom boundary
! programmer Florian Cachoux and Raphael F. Garcia
! in collaboration with D. Komatitsch and R. Martin
! variable forcing_type should be passed as a parameter in future versions

  subroutine acoustic_forcing_boundary(iglob,displ_x,displ_z)

  use constants, only: TINYVAL,ZERO
  use specfem_par

  implicit none

  integer,intent(in) :: iglob

  real(kind=CUSTOM_REAL) :: displ_x,displ_z

! local variables
  real, parameter :: pigrec = 3.1415927

  real :: alpha,tho,A,c,delayed,delta_x
  integer :: forcing_type,k,ngoce_time_step,n_models,kk,ll

  double precision, dimension(:), allocatable :: goce_time,distance
  double precision, dimension(:,:), allocatable :: syn
  double precision :: t,t_used,signal_x1,signal_x2,fracx,fract
  double precision :: x,z

  double precision :: f0 = 6000.0

  forcing_type = 1

  ! length of a PML element along x-axis at the edge which will be forced
  delta_x = 2.4
  alpha = 2

  ! infrasounds / seismic
  tho = 30.0

  A = 1
  x = coord(1,iglob)
  z = coord(2,iglob)
  delayed = 0

  ! speed of light
  c = 300000000.

  if (forcing_type == 1) then !! First test function : same forcing for the whole boundary
    !if (ispec_is_PML(ispec_acoustic)) then
    !  displ_x = 0
    !  displ_z = 0
    !else
      ! infrasounds / seismic
      !displ_x = 0 !* Apo
      !displ_z = A * (exp(-(alpha*(deltat*it-40-t0)/tho)**2) &
      !             - exp(-(alpha*(deltat*it-70-t0)/tho)**2)) !* Apo
      t_used = deltat*(it-1) - 0.0007d0
      !sin(2.0d0*pigrec*f0*t_used)
      if ((z < -1.5d0) .and. (z > -3.5d0)) then
        displ_x =  real(2.0d0 * f0*f0 * (2.0d0 * f0*f0 * t_used**2 - 1.0d0) * &
                             exp(-f0*f0*t_used**2),kind=CUSTOM_REAL) !(z+2.5d0)**3
        displ_z = 0._CUSTOM_REAL
      else
        displ_x =  0._CUSTOM_REAL
        displ_z = 0._CUSTOM_REAL
      endif
  endif

  !! Second test function : moving forcing
  if (forcing_type == 2) then
    displ_x = 0._CUSTOM_REAL !* Apo
    displ_z = real(dble(A) * (exp(-(alpha*(deltat*it-40-t0-(x-delayed)/c)/tho)**2) &
                 - exp(-(alpha*(deltat*it-70-t0-(x-delayed)/c)/tho)**2)),kind=CUSTOM_REAL) !* Apo
  endif

  !! forcing external
  if (forcing_type == 3) then
    ngoce_time_step = 255
    n_models = 28
    t =it*deltat

    allocate(goce_time(ngoce_time_step))
    allocate(distance(n_models))
    allocate(syn(n_models,ngoce_time_step))

    open(1000,file='../../EXAMPLES/acoustic_forcing_bottom/distance.txt',form='formatted')
    open(1001,file='../../EXAMPLES/acoustic_forcing_bottom/forcing_signals.txt',form='formatted')

    read(1001,*) goce_time(:)

    do k = 1,n_models
      read(1001,*) syn(k,:)
      read(1000,*) distance(k)
    enddo

    close(1000)
    close(1001)

    kk = 1
    do while(x >= distance(kk) .and. kk /= n_models)
      kk = kk+1
    enddo

    ll = 1
    do while(t >= goce_time(ll) .and. ll /= ngoce_time_step)
      ll = ll+1
    enddo

    if (x == 0 .and. it == 1) then
      displ_z =  real(syn(1,1),kind=CUSTOM_REAL)
    else
      if (x == 0) then
        fract = (t-goce_time(ll-1))/(goce_time(ll)-goce_time(ll-1))
        displ_z =  real((syn(1,ll-1) + fract * (syn(1,ll)-syn(1,ll-1))),kind=CUSTOM_REAL)
      else
        if (it == 1) then
          fracx = (x-distance(kk-1))/(distance(kk)-distance(kk-1))
          displ_z =  real((syn(kk-1,1) + fracx * (syn(kk,1)-syn(kk-1,1))),kind=CUSTOM_REAL)
        else
          ! interpolation in time
          fract = (t-goce_time(ll-1))/(goce_time(ll)-goce_time(ll-1))
          ! in x1 = distance(kk-1)
          signal_x1 = syn(kk-1,ll-1) + fract * (syn(kk-1,ll)-syn(kk-1,ll-1))
          ! in x2 = distance(kk)
          signal_x2 = syn(kk,ll-1) + fract * (syn(kk,ll)-syn(kk,ll-1))
          ! spatial interpolation
          fracx = (x-distance(kk-1))/(distance(kk)-distance(kk-1))
          displ_z =  real((signal_x1 + fracx * (signal_x2 - signal_x1)),kind=CUSTOM_REAL)
        endif
      endif
    endif

    displ_x = 0._CUSTOM_REAL

  endif

  if (abs(displ_x) < TINYVAL) displ_x = 0._CUSTOM_REAL
  if (abs(displ_z) < TINYVAL) displ_z = 0._CUSTOM_REAL

  end subroutine acoustic_forcing_boundary
