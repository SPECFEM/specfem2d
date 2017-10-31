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

  subroutine datim(string_input)

! get date and time

  use constants, only: IMAIN,MAX_STRING_LEN

  implicit none

  character(len=MAX_STRING_LEN) :: string_input

  ! local parameters
  character(len=8) :: datein
  character(len=10) :: timein
  character(len=16) :: dateprint
  character(len=8) :: timeprint

  datein = ' '
  timein = ' '

  call date_and_time(datein,timein)

  dateprint = datein(7:8)//' - '//datein(5:6)//' - '//datein(1:4)
  timeprint = timein(1:2)//':'//timein(3:4)//':'//timein(5:6)

  write(IMAIN,"(//1x,79('-')/1x,79('-')/1x,'Program SPECFEM2D: ')")
  write(IMAIN,"(1x,79('-')/1x,79('-')/1x,a)") trim(string_input)
  write(IMAIN,"(1x,79('-')/,1x,79('-')/' D a t e : ',a16,30x,' T i m e  : ',a8/1x,79('-'),/1x,79('-'))") dateprint,timeprint
  call flush_IMAIN()

  end subroutine datim

