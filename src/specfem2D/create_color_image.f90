
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently maNZ_IMAGE_color more authors!)
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

  subroutine create_color_image()


! display a given field as a red and blue color JPEG image

! to display the snapshots : display image*.jpg

  use specfem_par, only: image_color_data,iglob_image_color, &
                         NX_IMAGE_color,NZ_IMAGE_color,it,isnapshot_number,cutsnaps,image_color_vp_display, &
                         USE_SNAPSHOT_NUMBER_IN_FILENAME,POWER_DISPLAY_COLOR, &
                         DRAW_SOURCES_AND_RECEIVERS,NSOURCES,p_sv,nrec, &
                         ix_image_color_source,iy_image_color_source,ix_image_color_receiver,iy_image_color_receiver, &
                         USE_CONSTANT_MAX_AMPLITUDE,CONSTANT_MAX_AMPLITUDE_TO_USE

  implicit none

  include "constants.h"

  integer :: i

! for the JPEG library
  character(len=1), dimension(3,NX_IMAGE_color,NZ_IMAGE_color) :: JPEG_raw_image

  integer :: ix,iy,R,G,B

  double precision :: amplitude_max,normalized_value,vpmin,vpmax,x1

  character(len=100) :: filename


! size of cross and square in pixels drawn to represent the source and the receivers in JPEG pictures
  integer :: half_width_cross, thickness_cross, half_size_square

! make the size of the source and receiver symbols depend on the size of the picture
! using a rule of thumb
  thickness_cross = 1
  if(NX_IMAGE_color > 2000 .or. NZ_IMAGE_color > 2000) then
    half_width_cross = 6
    half_size_square = 4
  else if(NX_IMAGE_color <= 100 .or. NZ_IMAGE_color <= 100) then
    half_width_cross = 2
    half_size_square = 1
  else if(NX_IMAGE_color <= 250 .or. NZ_IMAGE_color <= 250) then
    half_width_cross = 3
    half_size_square = 2
  else
    half_width_cross = 5
    half_size_square = 3
  endif

! open the image file
! slightly change the beginning of the file name depending if we use the time step of the image number, to avoid confusion
  if(USE_SNAPSHOT_NUMBER_IN_FILENAME) then
    isnapshot_number = isnapshot_number + 1
    write(filename,"('OUTPUT_FILES/img',i7.7,'.jpg')") isnapshot_number
  else
    write(filename,"('OUTPUT_FILES/image',i7.7,'.jpg')") it
  endif

! compute maximum amplitude
  if(.not. USE_CONSTANT_MAX_AMPLITUDE) then
    amplitude_max = maxval(abs(image_color_data))
  else
    amplitude_max = CONSTANT_MAX_AMPLITUDE_TO_USE
!   in case of a pre-defined and constant maximum, truncate all values that are outside that constant range
    where(image_color_data > +CONSTANT_MAX_AMPLITUDE_TO_USE) image_color_data = +CONSTANT_MAX_AMPLITUDE_TO_USE
    where(image_color_data < -CONSTANT_MAX_AMPLITUDE_TO_USE) image_color_data = -CONSTANT_MAX_AMPLITUDE_TO_USE
  endif

  vpmin = HUGEVAL
  vpmax = TINYVAL
  do iy=1,NZ_IMAGE_color
    do ix=1,NX_IMAGE_color
! negative values in image_color_vp_display are a flag indicating a water layer to color in light blue later
      if ( iglob_image_color(ix,iy) > -1 .and. image_color_vp_display(ix,iy) >= 0) then
        vpmin = min(vpmin,image_color_vp_display(ix,iy))
        vpmax = max(vpmax,image_color_vp_display(ix,iy))
      endif

    enddo
  enddo

! in the image format, the image starts in the upper-left corner
  do iy=NZ_IMAGE_color,1,-1
    do ix=1,NX_IMAGE_color

! check if pixel is defined or not (can be above topography for instance)
      if(iglob_image_color(ix,iy) == -1) then

! use white to display undefined region above topography to avoid visual confusion with a water layer
        R = 255
        G = 255
        B = 255

! suppress small amplitudes considered as noise and display the background velocity model instead
      else if (abs(image_color_data(ix,iy)) < amplitude_max * cutsnaps) then

! use P velocity model as background where amplitude is negligible
        if((p_sv) .and. ((vpmax-vpmin)/max(vpmin, TINYVAL) > 0.02d0)) then
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

! negative values in image_color_vp_display are a flag indicating a water layer to color in light blue
        if (image_color_vp_display(ix,iy) < 0) then
! use light blue to display water
!!!!!!    R = 204
!!!!!!    G = 255
!!!!!!    B = 255
          R = 135 !!! LightSkyBlue
          G = 206
          B = 250
        endif

      else

! define normalized image data in [-1:1] and convert to nearest integer
! keeping in mind that data values can be negative
        if( amplitude_max >= TINYVAL ) then
          normalized_value = image_color_data(ix,iy) / amplitude_max
        else
          normalized_value = image_color_data(ix,iy) / TINYVAL
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
     JPEG_raw_image(1,ix,NZ_IMAGE_color-iy+1) = char(R)
     JPEG_raw_image(2,ix,NZ_IMAGE_color-iy+1) = char(G)
     JPEG_raw_image(3,ix,NZ_IMAGE_color-iy+1) = char(B)

    enddo
  enddo

!
!----  draw position of the sources and receivers
!
  if (DRAW_SOURCES_AND_RECEIVERS) then

! draw position of the sources with orange crosses
    do i=1,NSOURCES

! avoid edge effects for source or receiver symbols that can be partly outside of the image
      do iy = max(iy_image_color_source(i) - half_width_cross,1), min(iy_image_color_source(i) + half_width_cross,NZ_IMAGE_color)
        do ix = max(ix_image_color_source(i) - thickness_cross,1), min(ix_image_color_source(i) + thickness_cross,NX_IMAGE_color)
! use orange color
          R = 255
          G = 157
          B = 0
! for JPEG
          JPEG_raw_image(1,ix,NZ_IMAGE_color-iy+1) = char(R)
          JPEG_raw_image(2,ix,NZ_IMAGE_color-iy+1) = char(G)
          JPEG_raw_image(3,ix,NZ_IMAGE_color-iy+1) = char(B)
        enddo
      enddo

! avoid edge effects for source or receiver symbols that can be partly outside of the image
      do iy = max(iy_image_color_source(i) - thickness_cross,1), min(iy_image_color_source(i) + thickness_cross,NZ_IMAGE_color)
        do ix = max(ix_image_color_source(i) - half_width_cross,1), min(ix_image_color_source(i) + half_width_cross,NX_IMAGE_color)
! use orange color
          R = 255
          G = 157
          B = 0
! for JPEG
          JPEG_raw_image(1,ix,NZ_IMAGE_color-iy+1) = char(R)
          JPEG_raw_image(2,ix,NZ_IMAGE_color-iy+1) = char(G)
          JPEG_raw_image(3,ix,NZ_IMAGE_color-iy+1) = char(B)
        enddo
      enddo

    enddo

! draw position of the receivers with green squares
    do i=1,nrec
! avoid edge effects for source or receiver symbols that can be partly outside of the image
      do iy = max(iy_image_color_receiver(i) - half_size_square,1),&
                                          min(iy_image_color_receiver(i) + half_size_square,NZ_IMAGE_color)
        do ix = max(ix_image_color_receiver(i) - half_size_square,1), &
                                          min(ix_image_color_receiver(i) + half_size_square,NX_IMAGE_color)
! use dark green color
          R = 30
          G = 180
          B = 60
! for JPEG
          JPEG_raw_image(1,ix,NZ_IMAGE_color-iy+1) = char(R)
          JPEG_raw_image(2,ix,NZ_IMAGE_color-iy+1) = char(G)
          JPEG_raw_image(3,ix,NZ_IMAGE_color-iy+1) = char(B)
        enddo
      enddo
    enddo

  endif

! for JPEG
  call write_jpeg_image(JPEG_raw_image,NX_IMAGE_color,NZ_IMAGE_color,filename)

  end subroutine create_color_image

