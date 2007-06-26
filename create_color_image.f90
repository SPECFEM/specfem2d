
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!                         Dimitri Komatitsch
!                     University of Pau, France
!
!                          (c) April 2007
!
!========================================================================

  subroutine create_color_image(color_image_2D_data,iglob_image_color_2D,NX,NY,it,cutsnaps)

! display a given field as a red and blue color image

! to display the snapshots : display image*.gif

! when compiling with Intel ifort, use " -assume byterecl " option to create binary PNM images

  implicit none

  include "constants.h"

  integer NX,NY,it

  double precision cutsnaps

  integer, dimension(NX,NY) :: iglob_image_color_2D

  double precision, dimension(NX,NY) :: color_image_2D_data

  integer ix,iy,R,G,B,tenthousands,thousands,hundreds,tens,units,remainder,current_rec

  double precision amplitude_max,normalized_value

  character(len=100) file_name,system_command

! create temporary image files in binary PNM P6 format (smaller) or ASCII PNM P3 format (easier to edit)
  logical, parameter :: BINARY_FILE = .true.

! ASCII code of character '0' and of carriage return character
  integer, parameter :: ascii_code_of_zero = 48, ascii_code_of_carriage_return = 10

! open the image file
  write(file_name,"('OUTPUT_FILES/image',i6.6,'.pnm')") it

  if(BINARY_FILE) then

    open(unit=27,file=file_name,status='unknown',access='direct',recl=1)
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

    open(unit=27,file=file_name,status='unknown')
    write(27,"('P3')") ! write P3 = ASCII PNM image format
    write(27,*) NX,NY  ! write image size
    write(27,*) '255'  ! number of shades

  endif

! compute maximum amplitude
  amplitude_max = maxval(abs(color_image_2D_data))

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

! use black background where amplitude is negligible
          R = 0
          G = 0
          B = 0

      else

! define normalized image data in [-1:1] and convert to nearest integer
! keeping in mind that data values can be negative
        normalized_value = color_image_2D_data(ix,iy) / amplitude_max

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
  write(system_command,"('cd OUTPUT_FILES ; convert image',i6.6,'.pnm image',i6.6,'.gif')") it,it
 
! call the system to convert image to GIF
  call system(system_command)

  end subroutine create_color_image

