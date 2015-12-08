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


subroutine read_model_nspec()

! reads in nspec from database

  use tomography_par,only: MAX_STRING_LEN,IIN,NSPEC,myrank

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: m_file
  character(len=80) :: datlin

  ! opens database file
  write(m_file,'(a,i5.5)') './OUTPUT_FILES/'//'Database',myrank
  open(IIN,file=trim(m_file),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_MPI(myrank,'file not found')
  endif

  ! reads lines until nspec
  read(IIN,"(a80)") datlin  ! header
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin  ! simulation_title
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin ! AXISYM
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin  ! SIMULATION_TYPE, ...
  read(IIN,"(a80)") datlin
  read(IIN,*) nspec

  close(IIN)

  if (myrank == 0) then
    print *,'number of spectral-elements: ',nspec
    print *,''
  endif

end subroutine read_model_nspec
