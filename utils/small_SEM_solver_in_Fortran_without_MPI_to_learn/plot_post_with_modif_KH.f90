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

  subroutine plot_post(displ,coord,ibool,NGLOB,NSPEC,x_source,z_source,x_receiver,z_receiver,it,deltat,NGLLX,NGLLZ,NDIM, &
        is_on_left_edge_of_KH_contour,is_on_right_edge_of_KH_contour,is_on_bottom_edge_of_KH_contour,is_on_top_edge_of_KH_contour)

!
! PostScript display routine
!

  implicit none

  include "precision.h"

  double precision :: cutsnaps = 1.d0        ! minimum amplitude in % for snapshots
  double precision :: sizemax_arrows = 1.d0  ! maximum size of arrows on vector plots in cm
  double precision :: height = 0.25d0        ! height of domain numbers in centimeters

! US letter paper or European A4
  logical, parameter :: US_LETTER = .false.

! X and Z axis origin of PostScript plot in centimeters
  double precision, parameter :: ORIG_X = 2.4d0
  double precision, parameter :: ORIG_Z = 2.9d0

! dot to centimeter conversion for PostScript
  double precision, parameter :: CENTIM = 28.5d0

! parameters for arrows for PostScript snapshot
  double precision, parameter :: ARROW_ANGLE = 20.d0
  double precision, parameter :: ARROW_RATIO = 0.40d0

! size of frame used for Postscript display in percentage of the size of the page
  double precision, parameter :: RPERCENTX = 70.0d0,RPERCENTZ = 77.0d0

! pi
  double precision, parameter :: PI = 3.141592653589793d0

  integer :: it,NSPEC,NGLOB,NGLLX,NGLLZ,NDIM

  integer ibool(NGLLX,NGLLZ,NSPEC)

!! DK DK added this for KH
  logical, dimension(NSPEC) :: is_on_left_edge_of_KH_contour
  logical, dimension(NSPEC) :: is_on_right_edge_of_KH_contour
  logical, dimension(NSPEC) :: is_on_bottom_edge_of_KH_contour
  logical, dimension(NSPEC) :: is_on_top_edge_of_KH_contour

  real(kind=CUSTOM_REAL) :: deltat
  double precision :: timeval

  real(kind=CUSTOM_REAL) :: displ(NDIM,NGLOB)
  double precision :: coord(NDIM,NGLOB)

  double precision :: x_source,z_source
  double precision :: x_receiver,z_receiver

  double precision :: xmax,zmax,xw,zw,usoffset,sizex,sizez

! for the file name
  character(len=100) :: file_name

! to suppress useless white spaces in postscript lines
  character(len=100) :: postscript_line
  character(len=1), dimension(100) :: ch1,ch2
  equivalence (postscript_line,ch1)
  logical :: first

  double precision convert,x1,xa,za,xb,zb
  double precision z1,x2,z2,d,d1,d2,dummy,theta,thetaup,thetadown

  integer ispec,is,ir,line_length
  integer index_char,ii,ipoin

  double precision ratio_page,dispmax,xmin,zmin

! A4 or US letter paper
  if (US_LETTER) then
    usoffset = 1.75d0
    sizex = 27.94d0
    sizez = 21.59d0
  else
    usoffset = 0.d0
    sizex = 29.7d0
    sizez = 21.d0
  endif

! get minimum and maximum values of mesh coordinates
  xmin = minval(coord(1,:))
  zmin = minval(coord(2,:))
  xmax = maxval(coord(1,:))
  zmax = maxval(coord(2,:))
! write(*,*) 'X min, max in PostScript display = ',xmin,xmax
! write(*,*) 'Z min, max in PostScript display = ',zmin,zmax

! ratio of physical page size/size of the domain meshed
  ratio_page = min(rpercentz*sizez/(zmax-zmin),rpercentx*sizex/(xmax-xmin)) / 100.d0

! compute the maximum of the norm of the vector
  dispmax = maxval(sqrt(displ(1,:)**2 + displ(2,:)**2))
! write(*,*) 'Max norm of vector in PostScript display = ',dispmax
! write(*,*)

!
!---- open PostScript file
!
  write(file_name,"('vect',i6.6,'.ps')") it
  open(unit=24,file=file_name,status='unknown')

!
!---- write PostScript header
!
  write(24,10) 'SPECFEM2D demo code results'
  write(24,*) '/CM {28.5 mul} def'
  write(24,*) '/LR {rlineto} def'
  write(24,*) '/LT {lineto} def'
  write(24,*) '/L {lineto} def'
  write(24,*) '/MR {rmoveto} def'
  write(24,*) '/MV {moveto} def'
  write(24,*) '/M {moveto} def'
  write(24,*) '/ST {stroke} def'
  write(24,*) '/CP {closepath} def'
  write(24,*) '/RG {setrgbcolor} def'
  write(24,*) '/GF {gsave fill grestore} def'
  write(24,*) '% different useful symbols'
  write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
  write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/Cross {gsave 0.05 CM setlinewidth'
  write(24,*) 'gsave 3 3 MR -6. -6. LR ST grestore'
  write(24,*) 'gsave 3 -3 MR -6. 6. LR ST grestore'
  write(24,*) '0.01 CM setlinewidth} def'
  write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
  write(24,*) '/Diamond {gsave 0.05 CM setlinewidth 0 4.2 MR'
  write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
  write(24,*) 'grestore 0.01 CM setlinewidth} def'
  write(24,*) '%'
  write(24,*) '% gray levels for the velocity model'
  write(24,*) '/BK {setgray fill} def'
  write(24,*) '% black and white version'
  write(24,*) '%/BK {pop 1 setgray fill} def'
  write(24,*) '%'
  write(24,*) '% magenta for vectors'
  write(24,*) '/Colvects {0.01 CM setlinewidth 1. 0. 1. RG} def'
  write(24,*) '% black and white version'
  write(24,*) '%/Colvects {0.01 CM setlinewidth 0. setgray} def'
  write(24,*) '%'
  write(24,*) '% chartreuse for macrobloc mesh'
  write(24,*) '/Colmesh {0.02 CM setlinewidth 0.5 1. 0. RG} def'
  write(24,*) '% black and white version'
  write(24,*) '%/Colmesh {0.02 CM setlinewidth 0. setgray} def'
  write(24,*) '%'
  write(24,*) '% cyan for sources and receivers'
  write(24,*) '/Colreceiv {0. 1. 1. RG} def'
  write(24,*) '% black and white version'
  write(24,*) '%/Colreceiv {0. setgray} def'
  write(24,*) '%'
  write(24,*) '% macro to draw an arrow'
  write(24,*) '/F {MV LR gsave LR ST grestore LR ST} def'
  write(24,*) '% macro to draw the contour of the elements'
  write(24,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
  write(24,*) '%'
  write(24,*) '.01 CM setlinewidth'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.35 CM scalefont setfont'
  write(24,*) '%'
  write(24,*) '/vshift ',-height/2,' CM def'
  write(24,*) '/Rshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop neg vshift MR show } def'
  write(24,*) '/Cshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
  write(24,*) '/fN {/Helvetica-Bold findfont ',height,' CM scalefont setfont} def'
  write(24,*) '%'
  write(24,*) 'gsave newpath 90 rotate'
  write(24,*) '0 ',-sizez,' CM translate 1. 1. scale'
  write(24,*) '% uncomment this to zoom on parts of the mesh'
  write(24,*) '% -32 CM -21 CM translate 3. 3. scale'
  write(24,*) '% -52 CM -24 CM translate 4. 4. scale'
  write(24,*) '%'

!
!--- write captions of PostScript figure
!
  write(24,*) '0 setgray'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.5 CM scalefont setfont'

  write(24,*) '24. CM 1.2 CM MV'
  write(24,610) usoffset,it
  write(24,*) '%'

  write(24,*) '24. CM 1.95 CM MV'
  timeval = it*deltat
  if (timeval >= 1.d-3 .and. timeval < 1000.d0) then
    write(24,600) usoffset,timeval
  else
    write(24,601) usoffset,timeval
  endif
  write(24,*) '%'
  write(24,*) '24. CM 2.7 CM MV'
  write(24,640) usoffset,dispmax

  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.6 CM scalefont setfont'
  write(24,*) '11 CM 1.1 CM MV'
  write(24,*) '(X axis) show'
  write(24,*) '%'
  write(24,*) '1.4 CM 9.5 CM MV'
  write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
  write(24,*) '(Y axis) show'
  write(24,*) 'grestore'
  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.7 CM scalefont setfont'
  write(24,*) '24.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(Displacement vector field) show'
  write(24,*) 'grestore'
  write(24,*) '25.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(SPECFEM2D demo code results) show'
  write(24,*) 'grestore'
  write(24,*) '26.45 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'

  write(24,*) '(Elastic Wave 2D - Spectral Element Method) show'

  write(24,*) 'grestore'

  write(24,*) '%'
  write(24,*) '1 1 scale'
  write(24,*) '%'

!
!---- print the spectral elements mesh in PostScript
!

  convert = PI / 180.d0

!
!---- draw the spectral element mesh
!

  write(24,*) '%'
  write(24,*) '% spectral element mesh'
  write(24,*) '%'

  do ispec=1,NSPEC

  write(24,*) '% elem ',ispec

  is = 1
  ir = 1
  x1 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
  z1 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  write(24,*) 'mark'
  write(24,681) x1,z1

! draw straight lines for the element edges (if drawing curved elements with 9 nodes, they will thus appear with no curvature)

  ir=NGLLX
  x2 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
  z2 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  ir=NGLLX
  is=NGLLZ
  x2 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
  z2 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  is=NGLLZ
  ir=1
  x2 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
  z2 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  ir=1
  is=2
  x2 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
  z2 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  write(24,*) 'CO'
  write(24,*) '0 setgray ST'

  enddo

!
!---- draw the KH box on the mesh
!

!! DK DK added this for KH

  write(24,*) '%'
  write(24,*) '% KH box'
  write(24,*) '%'

  write(24,*) '0.05 CM setlinewidth'

  do ispec=1,NSPEC

  write(24,*) '% elem ',ispec

  if (is_on_left_edge_of_KH_contour(ispec)) then
    write(24,*) '1 0.6470 0 RG' ! Orange

    is = 1
    ir = 1
    x1 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
    z1 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
    x1 = x1 * centim
    z1 = z1 * centim

    is = NGLLZ
    ir = 1
    x2 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
    z2 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim

    write(24,602) x1,z1,x2,z2
  endif

  if (is_on_right_edge_of_KH_contour(ispec)) then
    write(24,*) '0 0 1 RG'  ! Blue

    is = 1
    ir = NGLLX
    x1 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
    z1 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
    x1 = x1 * centim
    z1 = z1 * centim

    is = NGLLZ
    ir = NGLLX
    x2 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
    z2 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim

    write(24,602) x1,z1,x2,z2
  endif

  if (is_on_bottom_edge_of_KH_contour(ispec)) then
    write(24,*) '0 1 0 RG'  ! Green

    is = 1
    ir = 1
    x1 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
    z1 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
    x1 = x1 * centim
    z1 = z1 * centim

    is = 1
    ir = NGLLX
    x2 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
    z2 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim

    write(24,602) x1,z1,x2,z2
  endif

  if (is_on_top_edge_of_KH_contour(ispec)) then
    write(24,*) '1 0.7529 0.7960 RG' ! Pink

    is = NGLLZ
    ir = 1
    x1 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
    z1 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
    x1 = x1 * centim
    z1 = z1 * centim

    is = NGLLZ
    ir = NGLLX
    x2 = (coord(1,ibool(ir,is,ispec))-xmin)*ratio_page + orig_x
    z2 = (coord(2,ibool(ir,is,ispec))-zmin)*ratio_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim

    write(24,602) x1,z1,x2,z2
  endif

  enddo

  write(24,*) '0.01 CM setlinewidth'
  write(24,*) '0 setgray'

!
!----  draw the normalized vector field
!

! return if the maximum vector equals zero (no source)
  if (dispmax == 0.d0) then
    write(*,*) 'null vector: returning!'
    return
  endif

  write(24,*) '%'
  write(24,*) '% vector field'
  write(24,*) '%'

  write(24,*) '0 setgray'

! draw the vectors at all the GLL nodes of the mesh
  do ipoin=1,NGLOB

  x1 =(coord(1,ipoin)-xmin)*ratio_page
  z1 =(coord(2,ipoin)-zmin)*ratio_page

  x2 = displ(1,ipoin)*sizemax_arrows/dispmax
  z2 = displ(2,ipoin)*sizemax_arrows/dispmax

  d = sqrt(x2**2 + z2**2)

! ignore if vector is too small
  if (d > cutsnaps*sizemax_arrows/100.d0) then

  d1 = d * ARROW_RATIO
  d2 = d1 * cos(ARROW_ANGLE*convert)

  dummy = x2/d
  if (dummy > 0.9999d0) dummy = 0.9999d0
  if (dummy < -0.9999d0) dummy = -0.9999d0
  theta = acos(dummy)

  if (z2 < 0.d0) theta = 360.d0*convert - theta
  thetaup = theta - ARROW_ANGLE*convert
  thetadown = theta + ARROW_ANGLE*convert

! draw the vector
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  x2 = x2 * centim
  z2 = z2 * centim
  xa = -d2*cos(thetaup)
  za = -d2*sin(thetaup)
  xa = xa * centim
  za = za * centim
  xb = -d2*cos(thetadown)
  zb = -d2*sin(thetadown)
  xb = xb * centim
  zb = zb * centim
  write(postscript_line,700) xb,zb,xa,za,x2,z2,x1,z1

! suppress useless white spaces to make PostScript file smaller

! suppress leading white spaces again, if any
  postscript_line = adjustl(postscript_line)

  line_length = len_trim(postscript_line)
  index_char = 1
  first = .false.
  do ii = 1,line_length-1
    if (ch1(ii) /= ' ' .or. first) then
      if (ch1(ii) /= ' ' .or. ch1(ii+1) /= ' ') then
        ch2(index_char) = ch1(ii)
        index_char = index_char + 1
        first = .true.
      endif
    endif
  enddo
  ch2(index_char) = ch1(line_length)
  write(24,"(100(a1))") (ch2(ii), ii=1,index_char)

  endif

  enddo

  write(24,*) '0 setgray'

!
!----  write position of the source
!
  xw = x_source
  zw = z_source
  xw = (xw-xmin)*ratio_page + orig_x
  zw = (zw-zmin)*ratio_page + orig_z
  xw = xw * centim
  zw = zw * centim
  write(24,500) xw,zw
  write(24,*) 'Cross'

!
!----  write position of the receiver
!
  write(24,*) '% position of the receiver'

  xw = x_receiver
  zw = z_receiver

  xw = (xw-xmin)*ratio_page + orig_x
  zw = (zw-zmin)*ratio_page + orig_z
  xw = xw * centim
  zw = zw * centim
  write(24,500) xw,zw
  write(24,*) 'Diamond'

  write(24,*) '%'
  write(24,*) 'grestore'
  write(24,*) 'showpage'

  close(24)

 10  format('%!PS-Adobe-2.0',/,'%%',/,'%% Title: ',a50,/,'%% Created by: Specfem2D',/,'%% Author: Dimitri Komatitsch',/,'%%')
 600 format(f6.3,' neg CM 0 MR (Time =',f8.3,' s) show')
 601 format(f6.3,' neg CM 0 MR (Time =',1pe12.3,' s) show')
 602 format(f6.2,1x,f6.2,' M ',f6.2,1x,f6.2,' L ST')
 610 format(f6.3,' neg CM 0 MR (Time step = ',i7,') show')
 640 format(f6.3,' neg CM 0 MR (Max norm =',1pe12.3,') show')

 500 format(f6.2,1x,f6.2,' M')
 681 format(f6.2,1x,f6.2)
 700 format(8(f6.2,1x),'F')

  end subroutine plot_post

