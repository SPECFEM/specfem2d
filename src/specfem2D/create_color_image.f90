
!========================================================================
!
!                   S P E C F E M 2 D  Version 6 . 2
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
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

  subroutine create_color_image(color_image_2D_data,iglob_image_color_2D, &
                                  NX,NY,it,cutsnaps,image_color_vp_display)

! display a given field as a red and blue color JPEG image

! to display the snapshots : display image*.jpg

  implicit none

  include "constants.h"

  integer :: NX,NY,it

  double precision :: cutsnaps

  integer, dimension(NX,NY) :: iglob_image_color_2D

  double precision, dimension(NX,NY) :: color_image_2D_data
  double precision, dimension(NX,NY) :: image_color_vp_display

! for JPEG
  character(len=1), dimension(3,NX,NY) :: JPEG_raw_image

  integer :: ix,iy,R,G,B

  double precision :: amplitude_max,normalized_value,vpmin,vpmax,x1

  character(len=100) :: filename

! open the image file
  write(filename,"('OUTPUT_FILES/image',i7.7,'.jpg')") it

! compute maximum amplitude
  amplitude_max = maxval(abs(color_image_2D_data))
  vpmin = HUGEVAL
  vpmax = TINYVAL
  do iy=1,NY
    do ix=1,NX
      if ( iglob_image_color_2D(ix,iy) > -1 ) then
        vpmin = min(vpmin,image_color_vp_display(ix,iy))
        vpmax = max(vpmax,image_color_vp_display(ix,iy))
      endif

    enddo
  enddo

! in the image format, the image starts in the upper-left corner
  do iy=NY,1,-1
    do ix=1,NX

! check if pixel is defined or not (can be above topography for instance)
      if(iglob_image_color_2D(ix,iy) == -1) then

!!! use light blue to display undefined region above topography
!!        R = 204
!!        G = 255
!!        B = 255
! now use white to display undefined region above topography to avoid visual confusion with a water layer
        R = 255
        G = 255
        B = 255

! suppress small amplitudes considered as noise
      else if (abs(color_image_2D_data(ix,iy)) < amplitude_max * cutsnaps) then

! use P velocity model as background where amplitude is negligible
        if((vpmax-vpmin)/vpmin > 0.02d0) then
          x1 = (image_color_vp_display(ix,iy)-vpmin)/(vpmax-vpmin)
        else
          x1 = 0.5d0
        endif

! rescale to avoid very dark gray levels
        x1 = x1*0.7 + 0.2
        if(x1 > 1.d0) x1=1.d0

! invert scale: white = vpmin, dark gray = vpmax
        x1 = 1.d0 - x1

! map to [0,255]
        x1 = x1 * 255.d0

        R = nint(x1)
        if(R < 0) R = 0
        if(R > 255) R = 255
        G = R
        B = R

      else

! define normalized image data in [-1:1] and convert to nearest integer
! keeping in mind that data values can be negative
        if( amplitude_max >= TINYVAL ) then
          normalized_value = color_image_2D_data(ix,iy) / amplitude_max
        else
          normalized_value = color_image_2D_data(ix,iy) / TINYVAL
        endif

! suppress values outside of [-1:+1]
        if(normalized_value < -1.d0) normalized_value = -1.d0
        if(normalized_value > 1.d0) normalized_value = 1.d0

! use red if positive value, blue if negative, no green
        if(normalized_value >= 0.d0) then
          R = nint(255.d0*normalized_value**POWER_DISPLAY_COLOR)
          G = 0
          B = 0
        else
          R = 0
          G = 0
          B = nint(255.d0*abs(normalized_value)**POWER_DISPLAY_COLOR)
        endif

      endif

! for JPEG
     JPEG_raw_image(1,ix,NY-iy+1) = char(R)
     JPEG_raw_image(2,ix,NY-iy+1) = char(G)
     JPEG_raw_image(3,ix,NY-iy+1) = char(B)

    enddo
  enddo

! for JPEG
  call write_jpeg_image(JPEG_raw_image,NX,NY,filename)

  end subroutine create_color_image

!========================================================================

! older (now unused) version that first creates a PNM image and then uses "convert" from ImageMagick to create a JPEG file
  subroutine create_color_image_from_PNM(color_image_2D_data,iglob_image_color_2D, &
                                  NX,NY,it,cutsnaps,image_color_vp_display)

! display a given field as a red and blue color JPEG image

! to display the snapshots : display image*.jpg

! when compiling with Intel ifort, use " -assume byterecl " option to create binary PNM images

  implicit none

  include "constants.h"

  integer :: NX,NY,it

  double precision :: cutsnaps

  integer, dimension(NX,NY) :: iglob_image_color_2D

  double precision, dimension(NX,NY) :: color_image_2D_data
  double precision, dimension(NX,NY) :: image_color_vp_display

  integer :: ix,iy,R,G,B,tenthousands,thousands,hundreds,tens,units,remainder,current_rec

  double precision :: amplitude_max,normalized_value,vpmin,vpmax,x1

  character(len=100) :: filename,system_command

! create temporary image files in binary PNM P6 format (smaller) or ASCII PNM P3 format (easier to edit)
  logical, parameter :: BINARY_FILE = .true.

! ASCII code of character '0' and of carriage return character
  integer, parameter :: ascii_code_of_zero = 48, ascii_code_of_carriage_return = 10

! open the image file
  write(filename,"('OUTPUT_FILES/image',i7.7,'.pnm')") it

  if(BINARY_FILE) then

    open(unit=27,file=filename,status='unknown',access='direct',recl=1)
    write(27,rec=1) 'P'
    write(27,rec=2) '6' ! write P6 = binary PNM image format
    write(27,rec=3) char(ascii_code_of_carriage_return)

! compute and write horizontal size
    remainder = NX

    tenthousands = remainder / 10000
    remainder = remainder - 10000 * tenthousands

    thousands = remainder / 1000
    remainder = remainder - 1000 * thousands

    hundreds = remainder / 100
    remainder = remainder - 100 * hundreds

    tens = remainder / 10
    remainder = remainder - 10 * tens

    units = remainder

    write(27,rec=4) char(tenthousands + ascii_code_of_zero)
    write(27,rec=5) char(thousands + ascii_code_of_zero)
    write(27,rec=6) char(hundreds + ascii_code_of_zero)
    write(27,rec=7) char(tens + ascii_code_of_zero)
    write(27,rec=8) char(units + ascii_code_of_zero)
    write(27,rec=9) ' '

! compute and write vertical size
    remainder = NY

    tenthousands = remainder / 10000
    remainder = remainder - 10000 * tenthousands

    thousands = remainder / 1000
    remainder = remainder - 1000 * thousands

    hundreds = remainder / 100
    remainder = remainder - 100 * hundreds

    tens = remainder / 10
    remainder = remainder - 10 * tens

    units = remainder

    write(27,rec=10) char(tenthousands + ascii_code_of_zero)
    write(27,rec=11) char(thousands + ascii_code_of_zero)
    write(27,rec=12) char(hundreds + ascii_code_of_zero)
    write(27,rec=13) char(tens + ascii_code_of_zero)
    write(27,rec=14) char(units + ascii_code_of_zero)
    write(27,rec=15) char(ascii_code_of_carriage_return)

! number of shades
    write(27,rec=16) '2'
    write(27,rec=17) '5'
    write(27,rec=18) '5'
    write(27,rec=19) char(ascii_code_of_carriage_return)

! block of image data starts at sixteenth character
    current_rec = 20

  else

    open(unit=27,file=filename,status='unknown')
    write(27,"('P3')") ! write P3 = ASCII PNM image format
    write(27,*) NX,NY  ! write image size
    write(27,*) '255'  ! number of shades

  endif

! compute maximum amplitude
  amplitude_max = maxval(abs(color_image_2D_data))
  vpmin = HUGEVAL
  vpmax = TINYVAL
  do iy=1,NY
    do ix=1,NX
      if ( iglob_image_color_2D(ix,iy) > -1 ) then
        vpmin = min(vpmin,image_color_vp_display(ix,iy))
        vpmax = max(vpmax,image_color_vp_display(ix,iy))
      endif

    enddo
  enddo

! in the PNM format, the image starts in the upper-left corner
  do iy=NY,1,-1
    do ix=1,NX

! check if pixel is defined or not (can be above topography for instance)
      if(iglob_image_color_2D(ix,iy) == -1) then

! use light blue to display undefined region above topography
        R = 204
        G = 255
        B = 255

! suppress small amplitudes considered as noise
      else if (abs(color_image_2D_data(ix,iy)) < amplitude_max * cutsnaps) then

! use P velocity model as background where amplitude is negligible
        if((vpmax-vpmin)/vpmin > 0.02d0) then
          x1 = (image_color_vp_display(ix,iy)-vpmin)/(vpmax-vpmin)
        else
          x1 = 0.5d0
        endif

! rescale to avoid very dark gray levels
        x1 = x1*0.7 + 0.2
        if(x1 > 1.d0) x1=1.d0

! invert scale: white = vpmin, dark gray = vpmax
        x1 = 1.d0 - x1

! map to [0,255]
        x1 = x1 * 255.d0

        R = nint(x1)
        if(R < 0) R = 0
        if(R > 255) R = 255
        G = R
        B = R

      else

! define normalized image data in [-1:1] and convert to nearest integer
! keeping in mind that data values can be negative
        if( amplitude_max >= TINYVAL ) then
          normalized_value = color_image_2D_data(ix,iy) / amplitude_max
        else
          normalized_value = color_image_2D_data(ix,iy) / TINYVAL
        endif

! suppress values outside of [-1:+1]
        if(normalized_value < -1.d0) normalized_value = -1.d0
        if(normalized_value > 1.d0) normalized_value = 1.d0

! use red if positive value, blue if negative, no green
        if(normalized_value >= 0.d0) then
          R = nint(255.d0*normalized_value**POWER_DISPLAY_COLOR)
          G = 0
          B = 0
        else
          R = 0
          G = 0
          B = nint(255.d0*abs(normalized_value)**POWER_DISPLAY_COLOR)
        endif

      endif

! write color image
      if(BINARY_FILE) then

! first write red
        write(27,rec=current_rec) char(R)
        current_rec = current_rec + 1

! then write green
        write(27,rec=current_rec) char(G)
        current_rec = current_rec + 1

! then write blue
        write(27,rec=current_rec) char(B)
        current_rec = current_rec + 1

      else

        write(27,"(i3,' ',i3,' ',i3)") R,G,B

      endif

    enddo
  enddo

! close the file
  close(27)

! open image file and create system command to convert image to more convenient format
! use the "convert" command from ImageMagick http://www.imagemagick.org
  write(system_command,"('cd OUTPUT_FILES ; convert image',i7.7,'.pnm image',i7.7,'.jpg ; rm -f image',i7.7,'.pnm')") it,it,it

! call the system to convert image to GIF
! this line can be safely commented out if your compiler does not implement "system()" for system calls;
! in such a case you will simply get images in PNM format in directory OUTPUT_FILES instead of GIF format
  call system(system_command)

  end subroutine create_color_image_from_PNM

