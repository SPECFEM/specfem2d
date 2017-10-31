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

subroutine read_model_nspec()

! reads in nspec from database

  use tomography_par, only: MAX_STRING_LEN,IIN,NSPEC,myrank

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: prname
  character(len=MAX_STRING_LEN) :: dummy
  integer :: idummy
  logical :: ldummy

  ! opens database file
  write(prname,"('./OUTPUT_FILES/Database',i5.5,'.bin')") myrank
  open(IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(prname)
    call exit_MPI(myrank,'file not found')
  endif

  ! reads lines until nspec
  read(IIN) dummy  ! simulation_title
  read(IIN) idummy,idummy,ldummy,ldummy  ! SIMULATION_TYPE, NOISE_TOMOGRAPHY, SAVE_FORWARD, UNDO_ATTENUATION_AND_OR_PML
  read(IIN) nspec

  close(IIN)

  ! user output
  if (myrank == 0) then
    print *,'number of spectral elements: ',nspec
    print *
  endif

end subroutine read_model_nspec
