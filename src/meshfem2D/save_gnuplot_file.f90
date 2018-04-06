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

  subroutine save_gnuplot_file(ngnod,nx,nz,x,z)

! creates a Gnuplot file that displays the grid

  use constants, only: IMAIN,OUTPUT_FILES

  implicit none

  integer :: ngnod,nx,nz
  double precision, dimension(0:nx,0:nz) :: x,z

  ! local parameters
  integer :: ier,istepx,istepz,ili,icol

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'Saving the grid in Gnuplot format...'
  write(IMAIN,*)

  open(unit=20,file=trim(OUTPUT_FILES)//'gridfile.gnu',status='unknown',iostat=ier)
  if (ier /= 0 ) then
    print *,'Error opening gnuplot file for writing: ',trim(OUTPUT_FILES)//'gridfile.gnu'
    print *,'Please make sure directory ',trim(OUTPUT_FILES),' exists...'
    call stop_the_code('Error saving gnuplot file')
  endif

  ! draw horizontal lines of the grid
  write(IMAIN,*) 'drawing horizontal lines of the grid'
  istepx = 1
  if (ngnod == 4) then
    istepz = 1
  else
    istepz = 2
  endif
  do ili=0,nz,istepz
    do icol=0,nx-istepx,istepx
       write(20,*) sngl(x(icol,ili)),sngl(z(icol,ili))
       write(20,*) sngl(x(icol+istepx,ili)),sngl(z(icol+istepx,ili))
       write(20,10)
    enddo
  enddo

  ! draw vertical lines of the grid
  write(IMAIN,*) 'drawing vertical lines of the grid'
  if (ngnod == 4) then
    istepx = 1
  else
    istepx = 2
  endif
  istepz = 1
  do icol=0,nx,istepx
    do ili=0,nz-istepz,istepz
       write(20,*) sngl(x(icol,ili)),sngl(z(icol,ili))
       write(20,*) sngl(x(icol,ili+istepz)),sngl(z(icol,ili+istepz))
       write(20,10)
    enddo
  enddo

10   format('')

  close(20)

  ! create a Gnuplot script to display the grid
  open(unit=20,file=trim(OUTPUT_FILES)//'plot_gridfile.gnu',status='unknown',iostat=ier)
  if (ier /= 0 ) call stop_the_code('Error saving plotgnu file')

  write(20,*) '#set term wxt'
  write(20,*) 'set term postscript landscape monochrome solid "Helvetica" 22'
  write(20,*) 'set output "',trim(OUTPUT_FILES)//'gridfile.ps"'
  write(20,*) '#set xrange [',sngl(minval(x)),':',sngl(maxval(x)),']'
  write(20,*) '#set yrange [',sngl(minval(z)),':',sngl(maxval(z)),']'
  ! use same unit length on both X and Y axes
  write(20,*) 'set size ratio -1'
  write(20,*) 'set loadpath "'//trim(OUTPUT_FILES)//'"'
  write(20,*) 'plot "gridfile.gnu" title "Macrobloc mesh" w l'
  write(20,*) 'pause -1 "Hit any key..."'
  close(20)

  write(IMAIN,*) 'Grid saved in Gnuplot format...'
  write(IMAIN,*)

  end subroutine save_gnuplot_file
