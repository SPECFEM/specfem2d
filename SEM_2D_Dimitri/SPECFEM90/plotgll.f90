!=====================================================================
!
!                 S p e c f e m  V e r s i o n  4 . 2
!                 -----------------------------------
!
!                         Dimitri Komatitsch
!    Department of Earth and Planetary Sciences - Harvard University
!                         Jean-Pierre Vilotte
!                 Departement de Sismologie - IPGP - Paris
!                           (c) June 1998
!
!=====================================================================

  subroutine plotgll(knods,ibool,coorg,coord)
!
!=======================================================================
!
!     "p l o t g l l" : Print the Gauss-Lobatto-Legendre mesh
!                       in a Gnuplot file
!
!=======================================================================
!
  use mesh01
  use spela202
  use iounit

  implicit none

  integer knods(ngnod,nspec),ibool(nxgll,nxgll,nspec)
  double precision coorg(ndime,npgeo),coord(ndime,npoin)

! coordinates of the nodes for Gnuplot file
  integer maxnnode
  parameter(maxnnode=9)
  real xval(maxnnode),zval(maxnnode)

  integer ispel,iy,ix,iglobnum,iglobnum2,ibloc,inode
  character(len=70) name

!
!---- print the GLL mesh in a Gnuplot file
!

  write(iout,*)
  write(iout,*) 'Generating gnuplot meshes...'
  write(iout,*)

! create non empty files for the case of 4-nodes elements

  name='macros1.gnu'
  open(unit=30,file=name,status='unknown')

  name='macros2.gnu'
  open(unit=31,file=name,status='unknown')
  write(31,10)

  name='gllmesh1.gnu'
  open(unit=20,file=name,status='unknown')

  name='gllmesh2.gnu'
  open(unit=21,file=name,status='unknown')
  write(21,10)

  do ispel = 1,nspec

!
!----    plot the lines in xi-direction
!
   do iy = 1,nygll
          do ix = 1,nxgll-1
!
!----   get the global point number
!
         iglobnum = ibool(ix,iy,ispel)
!
!----   do the same for next point on horizontal line
!
         iglobnum2 = ibool(ix+1,iy,ispel)

  write(20,15) sngl(coord(1,iglobnum)),sngl(coord(2,iglobnum))
  write(20,15) sngl(coord(1,iglobnum2)),sngl(coord(2,iglobnum2))
  write(20,10)

  if ((iy == 1).or.(iy == nygll)) then
  write(21,15) sngl(coord(1,iglobnum)),sngl(coord(2,iglobnum))
  write(21,15) sngl(coord(1,iglobnum2)),sngl(coord(2,iglobnum2))
  write(21,10)
  endif

          enddo
  enddo

!
!----    plot the lines in eta-direction
!
   do ix = 1,nxgll
          do iy = 1,nygll-1
!
!----   get the global point number
!
         iglobnum = ibool(ix,iy,ispel)
!
!----   do the same for next point on vertical line
!
         iglobnum2 = ibool(ix,iy+1,ispel)

  write(20,15) sngl(coord(1,iglobnum)),sngl(coord(2,iglobnum))
  write(20,15) sngl(coord(1,iglobnum2)),sngl(coord(2,iglobnum2))
  write(20,10)

  if ((ix == 1).or.(ix == nxgll)) then
  write(21,15) sngl(coord(1,iglobnum)),sngl(coord(2,iglobnum))
  write(21,15) sngl(coord(1,iglobnum2)),sngl(coord(2,iglobnum2))
  write(21,10)
  endif

          enddo
  enddo
  enddo

!
!----  Plot the macroblocs mesh using Gnuplot
!
  do ibloc = 1,nspec
  do inode = 1,ngnod

   xval(inode) =  sngl(coorg(1,knods(inode,ibloc)))
   zval(inode) =  sngl(coorg(2,knods(inode,ibloc)))

  enddo

  if(ngnod  ==  4) then
!
!----  4-noded rectangular element
!

! draw the edges of the element using one color
    write(30,15) xval(1),zval(1)
    write(30,15) xval(2),zval(2)
    write(30,10)
    write(30,15) xval(2),zval(2)
    write(30,15) xval(3),zval(3)
    write(30,10)
    write(30,15) xval(3),zval(3)
    write(30,15) xval(4),zval(4)
    write(30,10)
    write(30,15) xval(4),zval(4)
    write(30,15) xval(1),zval(1)
    write(30,10)

  else

!
!----  9-noded rectangular element
!

! draw the edges of the element using one color
    write(30,15) xval(1),zval(1)
    write(30,15) xval(5),zval(5)
    write(30,10)
    write(30,15) xval(5),zval(5)
    write(30,15) xval(2),zval(2)
    write(30,10)
    write(30,15) xval(2),zval(2)
    write(30,15) xval(6),zval(6)
    write(30,10)
    write(30,15) xval(6),zval(6)
    write(30,15) xval(3),zval(3)
    write(30,10)
    write(30,15) xval(3),zval(3)
    write(30,15) xval(7),zval(7)
    write(30,10)
    write(30,15) xval(7),zval(7)
    write(30,15) xval(4),zval(4)
    write(30,10)
    write(30,15) xval(4),zval(4)
    write(30,15) xval(8),zval(8)
    write(30,10)
    write(30,15) xval(8),zval(8)
    write(30,15) xval(1),zval(1)
    write(30,10)

! draw middle lines using another color
    write(31,15) xval(5),zval(5)
    write(31,15) xval(9),zval(9)
    write(31,10)
    write(31,15) xval(9),zval(9)
    write(31,15) xval(7),zval(7)
    write(31,10)
    write(31,15) xval(8),zval(8)
    write(31,15) xval(9),zval(9)
    write(31,10)
    write(31,15) xval(9),zval(9)
    write(31,15) xval(6),zval(6)
    write(31,10)

  endif

 enddo

  close(20)
  close(21)

  close(30)
  close(31)

!
!----  generate the command file for Gnuplot
!
  open(unit=20,file='plotmeshes',status='unknown')
  write(20,*) '#!/bin/sh'
  write(20,10)
  write(20,*) 'gnuplot macros_mesh.gnu'
  write(20,*) 'gnuplot gll_mesh.gnu'
  close(20)

  open(unit=20,file='gll_mesh.gnu',status='unknown')
  write(20,*) 'set term x11'
  write(20,*) 'set xlabel "X"'
  write(20,*) 'set ylabel "Y"'
  write(20,*) 'set title "Gauss-Lobatto-Legendre Mesh"'
  write(20,*) 'plot "gllmesh1.gnu" title '''' w l 2,', &
            ' "gllmesh2.gnu" title '''' w linesp 1 3'
  write(20,*) 'pause -1 "Hit any key to exit..."'
  close(20)

  open(unit=20,file='macros_mesh.gnu',status='unknown')
  write(20,*) 'set term x11'
  write(20,*) 'set xlabel "X"'
  write(20,*) 'set ylabel "Y"'
  write(20,*) 'set title "Spectral Elements (Macroblocs) Mesh"'
  write(20,*) 'plot "macros2.gnu" title '''' w l 2,', &
            ' "macros1.gnu" title '''' w linesp 1 3'
  write(20,*) 'pause -1 "Hit any key to exit..."'
  close(20)

!
!----
!

10 format('')
15 format(e10.5,1x,e10.5)

  return
  end subroutine plotgll
