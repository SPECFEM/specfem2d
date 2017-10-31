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

  subroutine plot_gll()

! output the Gauss-Lobatto-Legendre mesh in a gnuplot file

  use constants, only: IMAIN,NGLLX,NGLLZ,OUTPUT_FILES

  use specfem_par, only: knods,ibool,coorg,coord,ngnod,nspec

  implicit none

  integer iy,ix,iglobnum,iglobnum2,ibloc,inode,ispec

! coordinates of the nodes for Gnuplot file
  integer, parameter :: MAXNGNOD = 9
  double precision xval(MAXNGNOD),zval(MAXNGNOD)

  character(len=70) name

!
!---- output the GLL mesh in a Gnuplot file
!

  write(IMAIN,*)
  write(IMAIN,*) 'Generating gnuplot meshes...'
  write(IMAIN,*)
  call flush_IMAIN()

! create non empty files for the case of 4-node elements

  name='macros1.gnu'
  open(unit=30,file=trim(OUTPUT_FILES)//trim(name),status='unknown')

  name='macros2.gnu'
  open(unit=31,file=trim(OUTPUT_FILES)//trim(name),status='unknown')
  write(31,"('')")

  name='gllmesh1.gnu'
  open(unit=20,file=trim(OUTPUT_FILES)//trim(name),status='unknown')

  name='gllmesh2.gnu'
  open(unit=21,file=trim(OUTPUT_FILES)//trim(name),status='unknown')
  write(21,"('')")

  do ispec = 1,nspec

!
!----    plot the lines in xi-direction
!
    do iy = 1,NGLLZ
      do ix = 1,NGLLX-1
!
!----   get the global point number
!
        iglobnum = ibool(ix,iy,ispec)
!
!----   do the same for next point on horizontal line
!
        iglobnum2 = ibool(ix+1,iy,ispec)

        write(20,*) coord(1,iglobnum),coord(2,iglobnum)
        write(20,*) coord(1,iglobnum2),coord(2,iglobnum2)
        write(20,"('')")

        if (iy == 1 .or. iy == NGLLZ) then
          write(21,*) coord(1,iglobnum),coord(2,iglobnum)
          write(21,*) coord(1,iglobnum2),coord(2,iglobnum2)
          write(21,"('')")
        endif

      enddo
    enddo

!
!----    plot the lines in eta-direction
!
    do ix = 1,NGLLX
      do iy = 1,NGLLZ-1
!
!----   get the global point number
!
        iglobnum = ibool(ix,iy,ispec)
!
!----   do the same for next point on vertical line
!
        iglobnum2 = ibool(ix,iy+1,ispec)

        write(20,*) coord(1,iglobnum),coord(2,iglobnum)
        write(20,*) coord(1,iglobnum2),coord(2,iglobnum2)
        write(20,"('')")

        if (ix == 1 .or. ix == NGLLX) then
          write(21,*) coord(1,iglobnum),coord(2,iglobnum)
          write(21,*) coord(1,iglobnum2),coord(2,iglobnum2)
          write(21,"('')")
        endif

      enddo
    enddo
  enddo

!
!----  plot the macrobloc mesh using Gnuplot
!
  do ibloc = 1,nspec
    do inode = 1,ngnod

      xval(inode) = coorg(1,knods(inode,ibloc))
      zval(inode) = coorg(2,knods(inode,ibloc))

    enddo

    if (ngnod == 4) then
!
!----  4-node rectangular element
!

! draw the edges of the element using one color
      write(30,*) xval(1),zval(1)
      write(30,*) xval(2),zval(2)
      write(30,"('')")
      write(30,*) xval(2),zval(2)
      write(30,*) xval(3),zval(3)
      write(30,"('')")
      write(30,*) xval(3),zval(3)
      write(30,*) xval(4),zval(4)
      write(30,"('')")
      write(30,*) xval(4),zval(4)
      write(30,*) xval(1),zval(1)
      write(30,"('')")

    else

!
!----  9-node rectangular element
!

! draw the edges of the element using one color
      write(30,*) xval(1),zval(1)
      write(30,*) xval(5),zval(5)
      write(30,"('')")
      write(30,*) xval(5),zval(5)
      write(30,*) xval(2),zval(2)
      write(30,"('')")
      write(30,*) xval(2),zval(2)
      write(30,*) xval(6),zval(6)
      write(30,"('')")
      write(30,*) xval(6),zval(6)
      write(30,*) xval(3),zval(3)
      write(30,"('')")
      write(30,*) xval(3),zval(3)
      write(30,*) xval(7),zval(7)
      write(30,"('')")
      write(30,*) xval(7),zval(7)
      write(30,*) xval(4),zval(4)
      write(30,"('')")
      write(30,*) xval(4),zval(4)
      write(30,*) xval(8),zval(8)
      write(30,"('')")
      write(30,*) xval(8),zval(8)
      write(30,*) xval(1),zval(1)
      write(30,"('')")

! draw middle lines using another color
      write(31,*) xval(5),zval(5)
      write(31,*) xval(9),zval(9)
      write(31,"('')")
      write(31,*) xval(9),zval(9)
      write(31,*) xval(7),zval(7)
      write(31,"('')")
      write(31,*) xval(8),zval(8)
      write(31,*) xval(9),zval(9)
      write(31,"('')")
      write(31,*) xval(9),zval(9)
      write(31,*) xval(6),zval(6)
      write(31,"('')")

    endif

  enddo

  close(20)
  close(21)

  close(30)
  close(31)

!
!----  generate the command file for Gnuplot
!
  open(unit=20,file=trim(OUTPUT_FILES)//'plotall_gll_mesh.gnu',status='unknown')
  write(20,*) 'set term wxt'
  write(20,*) '# set term postscript landscape color solid "Helvetica" 22'
  write(20,*) '# set output "gll_mesh.ps"'
  write(20,*) 'set xlabel "X"'
  write(20,*) 'set ylabel "Y"'
  write(20,*) 'set title "Gauss-Lobatto-Legendre Mesh"'
  write(20,*) 'set size ratio -1'
  write(20,*) 'set loadpath "'//trim(OUTPUT_FILES)//'"'
  write(20,*) 'plot "gllmesh1.gnu" title '''' w l lc 2, "gllmesh2.gnu" title '''' w l lc 3'
  write(20,*) 'pause -1 "Hit any key to exit..."'
  close(20)

  open(unit=20,file=trim(OUTPUT_FILES)//'plotall_macro_mesh.gnu',status='unknown')
  write(20,*) 'set term wxt'
  write(20,*) '# set term postscript landscape color solid "Helvetica" 22'
  write(20,*) '# set output "macro_mesh.ps"'
  write(20,*) 'set xlabel "X"'
  write(20,*) 'set ylabel "Y"'
  write(20,*) 'set title "Spectral Element (Macrobloc) Mesh"'
  write(20,*) 'set size ratio -1'
  write(20,*) 'set loadpath "'//trim(OUTPUT_FILES)//'"'
  write(20,*) 'plot "macros2.gnu" title '''' w l lc 2, "macros1.gnu" title '''' w l lc 3'
  write(20,*) 'pause -1 "Hit any key to exit..."'
  close(20)

  end subroutine plot_gll

