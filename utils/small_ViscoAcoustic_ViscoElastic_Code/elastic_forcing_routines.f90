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

subroutine build_forced(xforced,coord,NGLOB,forced,nb_forced)
  ! This subroutine compute which global point is forced and which is not
  ! depending on what we want

  implicit none

  ! Input variables
  integer, intent(in) :: NGLOB
  double precision, intent(in) :: xforced
  double precision, dimension(2,NGLOB), intent(in) :: coord

  ! Output variables
  logical, dimension(NGLOB), intent(out) :: forced
  integer, intent(out) :: nb_forced

  ! Local variables
  integer :: iglob
  double precision :: TINYVAL = 1.0d-30

  nb_forced=0
  do iglob = 1,NGLOB
     if (abs(coord(1,iglob) - xforced) < TINYVAL) then
        forced(iglob) = .false.
        nb_forced=nb_forced+1
     endif
  enddo

end subroutine build_forced

! =======================================================================

subroutine enforce_elastic_forcing(NGLOB,iglob,it,deltat,deltatover2,deltatsqover2,coord,displ,veloc,accel)
  ! This subroutine set the acceleration at a GLL point iglob
  ! It then calculate numerically the velocity and displacement (which are also initialized here) from a Newmark scheme

  implicit none

  ! Input variables
  integer, intent(in) :: NGLOB,iglob,it
  real(kind=4), intent(in) :: deltat,deltatover2,deltatsqover2
  double precision, dimension(2,NGLOB), intent(in) :: coord

  ! Output variables
  real(kind=4), dimension(2,NGLOB), intent(inout) :: displ,veloc
  real(kind=4), dimension(2,NGLOB), intent(out) :: accel

  ! Local variables
  real(kind=4), dimension(2) :: accelOld,velocOld
  double precision :: factor

  factor = 1.0d0*(coord(2,iglob) - 1500.0d0)**2

  !if (abs(coord(2,iglob) - 1500.0d0) < 500.0d0) then
    if (it == 1) then ! We initialize the variables
      displ(1,iglob) = factor*(-1.0d0)
      displ(2,iglob) = factor*0.0d0
      veloc(1,iglob) = factor*0.0d0
      veloc(2,iglob) = factor*0.0d0
      accel(1,iglob) = factor*900.0d0
      accel(2,iglob) = factor*0.0d0
    else ! We set what we want
      accel(1,iglob) = factor*30.0**2*cos(30.0*(it-1)*deltat)
      accelOld(1) = factor*30.0**2*cos(30.0*(it-2)*deltat)
      accel(2,iglob) = 0.0d0
      accelOld(2) = 0.0d0
      ! Do not change anything below: we compute numerically the velocity and displacement
      velocOld(1) = veloc(1,iglob)
      !print *,iglob,veloc(1,iglob),accelOld(1),accel(1,iglob),deltatover2
      veloc(1,iglob) = veloc(1,iglob) + deltatover2*(accelOld(1) + accel(1,iglob))
      !print *,"OK"
      displ(1,iglob) = displ(1,iglob) + deltat*velocOld(1) + deltatsqover2*accelOld(1)
      velocOld(2) = veloc(2,iglob)
      veloc(2,iglob) = veloc(2,iglob) + deltatover2*(accelOld(2) + accel(2,iglob))
      displ(2,iglob) = displ(2,iglob) + deltat*velocOld(2) + deltatsqover2*accelOld(2)
    endif
  !else !nodes clamped
  !  displ(:,iglob) = 0.0d0
  !  veloc(:,iglob) = 0.0d0
  !  accel(:,iglob) = 0.0d0
  !endif

end subroutine enforce_elastic_forcing
!=====================================================================================================================

subroutine enforce_acoustic_forcing(NGLOB,iglob,it,deltat,deltatover2,deltatsqover2,coord, &
     potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic)
  ! This subroutine set the potential_dot_dot_acoustic at a GLL point iglob
  ! It then calculate numerically the potential_dott_acoustic and potential_acoustic (which are also initialized here) from a Newmark scheme

  implicit none

  ! Input variables
  integer, intent(in) :: NGLOB,iglob,it
  real(kind=4), intent(in) :: deltat,deltatover2,deltatsqover2
  double precision, dimension(2,NGLOB), intent(in) :: coord

  ! Output variables
  real(kind=4), dimension(NGLOB), intent(inout) :: potential_acoustic,potential_dot_acoustic
  real(kind=4), dimension(NGLOB), intent(out) :: potential_dot_dot_acoustic

  ! Local variables
  real(kind=4) :: potential_dot_dot_acousticOld,potential_dot_acousticOld
  double precision :: factor

  factor = 1.0d0 !*(coord(2,iglob) - 1500.0d0)**2

  !if (abs(coord(2,iglob) - 1500.0d0) < 500.0d0) then
  if (it == 1) then ! We initialize the variables
     potential_acoustic(iglob) = factor*(-1.0d0)
     potential_dot_acoustic(iglob) = factor*0.0d0
     potential_dot_dot_acoustic(iglob) = factor*900.0d0


  else ! We set what we want
     potential_dot_dot_acoustic(iglob) = factor*30.0**2*cos(30.0*(it-1)*deltat)
     potential_dot_dot_acousticOld = factor*30.0**2*cos(30.0*(it-2)*deltat)
     ! Do not change anything below: we compute numerically the velocity and displacement
     potential_dot_acousticOld = potential_dot_acoustic(iglob)
     potential_dot_acoustic(iglob) = potential_dot_acoustic(iglob) + &
          deltatover2*(potential_dot_dot_acousticOld + potential_dot_dot_acoustic(iglob))
     potential_acoustic(iglob) = potential_acoustic(iglob) + &
          deltat* potential_dot_acousticOld + deltatsqover2*potential_dot_dot_acousticOld
  endif
  !else !nodes clamped
  !  potential_acoustic(iglob) = 0.0d0
  !  potential_dot_acoustic(iglob) = 0.0d0
  !  potential_dot_dot_acoustic(iglob) = 0.0d0
  !endif

end subroutine enforce_acoustic_forcing
