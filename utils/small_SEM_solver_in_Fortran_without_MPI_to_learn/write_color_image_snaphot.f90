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

  subroutine write_color_image_snaphot(it,NX_IMAGE_color,NZ_IMAGE_color,NDIM,NGLOB,displ,veloc,accel, &
                  vector_field_display,image_color_data,iglob_image_color,ix_image_color_source,iy_image_color_source, &
                  ix_image_color_receiver,iy_image_color_receiver,isnapshot_number,NSOURCES,nrec,cp)

  implicit none

  include "precision.h"

  integer :: it,NX_IMAGE_color,NZ_IMAGE_color,NDIM,NGLOB,isnapshot_number,NSOURCES,nrec
  real(kind=CUSTOM_REAL) :: cp

  integer, dimension(NX_IMAGE_color,NZ_IMAGE_color) :: iglob_image_color
  double precision, dimension(NX_IMAGE_color,NZ_IMAGE_color) :: image_color_data
  integer, dimension(NSOURCES) :: ix_image_color_source,iy_image_color_source
  integer, dimension(nrec) :: ix_image_color_receiver,iy_image_color_receiver

! global displacement, velocity and acceleration vectors
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: displ,veloc,accel

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: vector_field_display

  ! local variables
  integer :: i,j

! parameter that can be changed if needed
! display 1=displ_Ux 2=displ_Uz 3=displ_norm 4=veloc_Vx 5=veloc_Vz 6=veloc_norm 7=accel_Ax 8=accel_Az 9=accel_norm
  integer, parameter :: imagetype_JPEG = 2

  write(*,*) 'Creating color image of size ',NX_IMAGE_color,' x ',NZ_IMAGE_color

  if (imagetype_JPEG >= 1 .and. imagetype_JPEG <= 3) then
    write(*,*) 'drawing scalar image of part of the displacement vector...'
    vector_field_display(:,:) = displ(:,:)

  else if (imagetype_JPEG >= 4 .and. imagetype_JPEG <= 6) then
    write(*,*) 'drawing scalar image of part of the velocity vector...'
    vector_field_display(:,:) = veloc(:,:)

  else if (imagetype_JPEG >= 7 .and. imagetype_JPEG <= 9) then
    write(*,*) 'drawing scalar image of part of the acceleration vector...'
    vector_field_display(:,:) = accel(:,:)

  else
    stop 'wrong type for JPEG snapshots'
  endif

  image_color_data(:,:) = 0.d0

  do j = 1, NZ_IMAGE_color
    do i = 1, NX_IMAGE_color

      ! P-SV waves, plot a component of vector, its norm, or else pressure
      if (iglob_image_color(i,j) /= -1) then
        if (imagetype_JPEG == 1 .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7) then
          ! draw the X component of the vector
          image_color_data(i,j) = vector_field_display(1,iglob_image_color(i,j))

        else if (imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8) then
          ! draw the Z component of the vector
          image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))

        else if (imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9) then
          ! draw the norm of the vector
          image_color_data(i,j) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2  &
                                     + vector_field_display(2,iglob_image_color(i,j))**2)

        else
          stop 'wrong type for JPEG snapshots'
        endif
      endif

    enddo
  enddo

  ! creates image
  call create_color_image(it,NSOURCES,nrec,NX_IMAGE_color,NZ_IMAGE_color,isnapshot_number, &
    image_color_data,iglob_image_color, &
    ix_image_color_source,iy_image_color_source,ix_image_color_receiver,iy_image_color_receiver,cp)

  ! user output
  write(*,*) 'Color image created'

  end subroutine write_color_image_snaphot

