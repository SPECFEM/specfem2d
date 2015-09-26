
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


  subroutine save_gnuplot_file(ngnod,nx,nz,x,z)

! creates a Gnuplot file that displays the grid

  implicit none

  integer :: ngnod,nx,nz
  double precision, dimension(0:nx,0:nz) :: x,z

  ! local parameters
  integer :: ios,istepx,istepz,ili,icol

  print *
  print *,'Saving the grid in Gnuplot format...'

  open(unit=20,file='OUTPUT_FILES/gridfile.gnu',status='unknown',iostat=ios)
  if( ios /= 0 ) stop 'error saving gnuplot file'

  ! draw horizontal lines of the grid
  print *,'drawing horizontal lines of the grid'
  istepx = 1
  if(ngnod == 4) then
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
  print *,'drawing vertical lines of the grid'
  if(ngnod == 4) then
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
  open(unit=20,file='OUTPUT_FILES/plotgnu',status='unknown',iostat=ios)
  if( ios /= 0 ) stop 'error saving plotgnu file'

  write(20,*) '#set term wxt'
  write(20,*) 'set term postscript landscape monochrome solid "Helvetica" 22'
  write(20,*) 'set output "grid.ps"'
  write(20,*) '#set xrange [',sngl(minval(x)),':',sngl(maxval(x)),']'
  write(20,*) '#set yrange [',sngl(minval(z)),':',sngl(maxval(z)),']'
  ! use same unit length on both X and Y axes
  write(20,*) 'set size ratio -1'
  write(20,*) 'plot "gridfile.gnu" title "Macrobloc mesh" w l'
  write(20,*) 'pause -1 "Hit any key..."'
  close(20)

  print *,'Grid saved in Gnuplot format...'
  print *

  end subroutine save_gnuplot_file
