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

  subroutine create_color_image(i_field,plot_b_wavefield_only)

! display a given field as a red and blue color JPEG image

! to display the snapshots : display image*.jpg

  use constants, only: TINYVAL,VERYTINYVAL,HUGEVAL,STABILITY_THRESHOLD,OUTPUT_FILES,IMAIN

  use specfem_par, only: myrank,it,NSOURCES,P_SV,nrec

  use specfem_par_movie, only: image_color_data,iglob_image_color,NX_IMAGE_color,NZ_IMAGE_color, &
    isnapshot_number,cutsnaps,image_color_vp_display, &
    USE_SNAPSHOT_NUMBER_IN_FILENAME,POWER_DISPLAY_COLOR, &
    DRAW_SOURCES_AND_RECEIVERS, &
    ix_image_color_source,iy_image_color_source,ix_image_color_receiver,iy_image_color_receiver, &
    USE_CONSTANT_MAX_AMPLITUDE,CONSTANT_MAX_AMPLITUDE_TO_USE,SIMULATION_TYPE

  implicit none

  integer :: i_field
  logical :: plot_b_wavefield_only

  ! local parameters
  integer :: i

  ! for the JPEG library
  character(len=1), dimension(3,NX_IMAGE_color,NZ_IMAGE_color) :: JPEG_raw_image
  integer :: ix,iy,R,G,B
  double precision :: amplitude_max,normalized_value,vpmin,vpmax,x1
  character(len=100) :: filename
  logical :: do_warning

  ! size of cross and square in pixels drawn to represent the source and the receivers in JPEG pictures
  integer :: half_width_cross, thickness_cross, half_size_square

! make the size of the source and receiver symbols depend on the size of the picture
! using a rule of thumb
  thickness_cross = 1
  if (NX_IMAGE_color > 2000 .or. NZ_IMAGE_color > 2000) then
    half_width_cross = 6
    half_size_square = 4
  else if (NX_IMAGE_color <= 100 .or. NZ_IMAGE_color <= 100) then
    half_width_cross = 2
    half_size_square = 1
  else if (NX_IMAGE_color <= 250 .or. NZ_IMAGE_color <= 250) then
    half_width_cross = 3
    half_size_square = 2
  else
    half_width_cross = 5
    half_size_square = 3
  endif

! open the image file
! slightly change the beginning of the file name depending if we use the time step of the image number, to avoid confusion
  if (USE_SNAPSHOT_NUMBER_IN_FILENAME) then
    isnapshot_number = isnapshot_number + 1
    if (i_field == 1 .and. SIMULATION_TYPE == 1) then
      write(filename,"(a,i7.7,a)") trim(OUTPUT_FILES)//'forward_img',isnapshot_number,'.jpg'
    else if (i_field == 1 .and. SIMULATION_TYPE == 3 .and. .not. plot_b_wavefield_only) then
      write(filename,"(a,i7.7,a)") trim(OUTPUT_FILES)//'adjoint_img',isnapshot_number,'.jpg'
    else
      write(filename,"(a,i7.7,a)") trim(OUTPUT_FILES)//'b_forward_img',isnapshot_number,'.jpg'
    endif
  else
    if (i_field == 1 .and. SIMULATION_TYPE == 1) then
      write(filename,"(a,i7.7,a)") trim(OUTPUT_FILES)//'forward_image',it,'.jpg'
    else if (i_field == 1 .and. SIMULATION_TYPE == 3 .and. .not. plot_b_wavefield_only) then
      write(filename,"(a,i7.7,a)") trim(OUTPUT_FILES)//'adjoint_image',it,'.jpg'
    else
      write(filename,"(a,i7.7,a)") trim(OUTPUT_FILES)//'b_forward_image',it,'.jpg'
    endif
  endif

! compute maximum amplitude
  if (.not. USE_CONSTANT_MAX_AMPLITUDE) then
    amplitude_max = maxval(abs(image_color_data))
  else
    amplitude_max = CONSTANT_MAX_AMPLITUDE_TO_USE
!   in case of a pre-defined and constant maximum, truncate all values that are outside that constant range
    where(image_color_data > +CONSTANT_MAX_AMPLITUDE_TO_USE) image_color_data = +CONSTANT_MAX_AMPLITUDE_TO_USE
    where(image_color_data < -CONSTANT_MAX_AMPLITUDE_TO_USE) image_color_data = -CONSTANT_MAX_AMPLITUDE_TO_USE
  endif
  ! user output
  write(IMAIN,*) 'Color image maximum amplitude = ',amplitude_max


! this trick checks for NaN (Not a Number), which is not even equal to itself
  if (amplitude_max > STABILITY_THRESHOLD .or. amplitude_max < 0 .or. amplitude_max /= amplitude_max) then
    print *,'Warning: failed creating color image, maximum value of amplitude in image color is invalid'
    print *,'amplitude max = ',amplitude_max,' with threshold at ', STABILITY_THRESHOLD
    print *,'Please check your simulation setup...'
    call exit_MPI(myrank,'code became unstable and blew up (image_color_data)')
  endif
  do_warning = .false.

  vpmin = HUGEVAL
  vpmax = TINYVAL
  do iy= 1,NZ_IMAGE_color
    do ix= 1,NX_IMAGE_color
! negative values in image_color_vp_display are a flag indicating a water layer to color in light blue later
      if (iglob_image_color(ix,iy) > -1 .and. image_color_vp_display(ix,iy) >= 0) then
        vpmin = min(vpmin,image_color_vp_display(ix,iy))
        vpmax = max(vpmax,image_color_vp_display(ix,iy))
      endif
    enddo
  enddo

! in the image format, the image starts in the upper-left corner
  do iy=NZ_IMAGE_color,1,-1
    do ix= 1,NX_IMAGE_color

! check if pixel is defined or not (can be above topography for instance)
      if (iglob_image_color(ix,iy) == -1) then

! use white to display undefined region above topography to avoid visual confusion with a water layer
        R = 255
        G = 255
        B = 255

! suppress small amplitudes considered as noise and display the background velocity model instead
      else if (abs(image_color_data(ix,iy)) < amplitude_max * cutsnaps) then

! use P velocity model as background where amplitude is negligible
        if ((P_SV) .and. ((vpmax-vpmin)/max(vpmin, TINYVAL) > 0.02d0)) then
          x1 = (image_color_vp_display(ix,iy)-vpmin)/(vpmax-vpmin)
        else
          x1 = 0.5d0
        endif

! rescale to avoid very dark gray levels
        x1 = x1*0.7 + 0.2
        if (x1 > 1.d0) x1=1.d0

! invert scale: white = vpmin, dark gray = vpmax
        x1 = 1.d0 - x1

! map to [0,255]
        x1 = x1 * 255.d0

        R = nint(x1)
        if (R < 0) R = 0
        if (R > 255) R = 255
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
        if (amplitude_max >= VERYTINYVAL) then
          normalized_value = image_color_data(ix,iy) / amplitude_max
        else
          normalized_value = image_color_data(ix,iy) / VERYTINYVAL
        endif
        ! check value (isNaN)
        if (normalized_value /= normalized_value) then
          ! will be set to zero
          normalized_value = 0.d0
          do_warning = .true.
        endif

! suppress values outside of [-1:+1]
        if (normalized_value < -1.d0) normalized_value = -1.d0
        if (normalized_value > 1.d0) normalized_value = 1.d0

! use red if positive value, blue if negative, no green
        if (normalized_value >= 0.d0) then
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

  if (do_warning) then
    print *,'Warning: normalized value is invalid, likely due to unstable simulation! process rank is ',myrank
    print *,'Warning: normalized_value = ',normalized_value,' will be set to zero'
  endif

!
!----  draw position of the sources and receivers
!
  if (DRAW_SOURCES_AND_RECEIVERS) then

! draw position of the sources with orange crosses
    do i = 1,NSOURCES

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
    do i = 1,nrec
! avoid edge effects for source or receiver symbols that can be partly outside of the image
      do iy = max(iy_image_color_receiver(i) - half_size_square,1), &
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

