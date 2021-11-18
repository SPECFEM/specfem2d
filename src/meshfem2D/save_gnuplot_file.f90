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

  subroutine save_gnuplot_file(NGNOD,nx,nz,x,z)

! creates a Gnuplot file that displays the grid

  use constants, only: IMAIN,IOUT_VIS,OUTPUT_FILES,myrank

  implicit none

  integer :: NGNOD,nx,nz
  double precision, dimension(0:nx,0:nz) :: x,z

  ! local parameters
  integer :: ier,istepx,istepz,ili,icol

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Saving the grid in Gnuplot format...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  open(unit=IOUT_VIS,file=trim(OUTPUT_FILES)//'gridfile.gnu',status='unknown',iostat=ier)
  if (ier /= 0 ) then
    print *,'Error opening gnuplot file for writing: ',trim(OUTPUT_FILES)//'gridfile.gnu'
    print *,'Please make sure directory ',trim(OUTPUT_FILES),' exists...'
    call stop_the_code('Error saving gnuplot file')
  endif

  ! draw horizontal lines of the grid
  if (myrank == 0) write(IMAIN,*) 'drawing horizontal lines of the grid'
  istepx = 1
  if (NGNOD == 4) then
    istepz = 1
  else
    istepz = 2
  endif
  do ili=0,nz,istepz
    do icol=0,nx-istepx,istepx
       write(IOUT_VIS,*) sngl(x(icol,ili)),sngl(z(icol,ili))
       write(IOUT_VIS,*) sngl(x(icol+istepx,ili)),sngl(z(icol+istepx,ili))
       write(IOUT_VIS,10)
    enddo
  enddo

  ! draw vertical lines of the grid
  if (myrank == 0) write(IMAIN,*) 'drawing vertical lines of the grid'
  if (NGNOD == 4) then
    istepx = 1
  else
    istepx = 2
  endif
  istepz = 1
  do icol=0,nx,istepx
    do ili=0,nz-istepz,istepz
       write(IOUT_VIS,*) sngl(x(icol,ili)),sngl(z(icol,ili))
       write(IOUT_VIS,*) sngl(x(icol,ili+istepz)),sngl(z(icol,ili+istepz))
       write(IOUT_VIS,10)
    enddo
  enddo

10   format('')

  close(IOUT_VIS)

  ! create a Gnuplot script to display the grid
  open(unit=IOUT_VIS,file=trim(OUTPUT_FILES)//'plot_gridfile.gnu',status='unknown',iostat=ier)
  if (ier /= 0 ) call stop_the_code('Error saving plotgnu file')

  write(IOUT_VIS,*) '#set term wxt'
  write(IOUT_VIS,*) 'set term postscript landscape monochrome solid "Helvetica" 22'
  write(IOUT_VIS,*) 'set output "',trim(OUTPUT_FILES)//'gridfile.ps"'
  write(IOUT_VIS,*) '#set xrange [',sngl(minval(x)),':',sngl(maxval(x)),']'
  write(IOUT_VIS,*) '#set yrange [',sngl(minval(z)),':',sngl(maxval(z)),']'
  ! use same unit length on both X and Y axes
  write(IOUT_VIS,*) 'set size ratio -1'
  write(IOUT_VIS,*) 'set loadpath "'//trim(OUTPUT_FILES)//'"'
  write(IOUT_VIS,*) 'plot "gridfile.gnu" title "Macrobloc mesh" w l'
  write(IOUT_VIS,*) 'pause -1 "Hit any key..."'
  close(IOUT_VIS)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Grid saved in Gnuplot format...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine save_gnuplot_file
