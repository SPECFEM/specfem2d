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

! create a sorted version of the indirect addressing to reduce cache misses

  subroutine get_global()

  use constants, only: NGLLX,NGLLZ,MAX_STRING_LEN,IN_DATA_FILES

  use specfem_par, only: nspec,ibool,copy_ibool_ori,integer_mask_ibool,SAVE_MODEL,myrank

  implicit none

  ! local parameters
  integer :: inumber,ispec,i,j
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname

  ! initializes temporary arrays
  integer_mask_ibool(:) = -1
  copy_ibool_ori(:,:,:) = ibool(:,:,:)

  inumber = 0

  ! reduce cache misses in all the elements
  ! loop over spectral elements
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        if (integer_mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
          ! create a new point
          inumber = inumber + 1
          ibool(i,j,ispec) = inumber
          integer_mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
        else
          ! use an existing point created previously
          ibool(i,j,ispec) = integer_mask_ibool(copy_ibool_ori(i,j,ispec))
        endif
      enddo
    enddo
  enddo

  if (trim(SAVE_MODEL) == 'binary') then
    write(outputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_NSPEC_ibool.bin'

    open(888,file=trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening smoothed kernel file')
    write(888) nspec
    write(888) ibool
    close(888)
  endif

  end subroutine get_global

!
!-------------------------------------------------------------------------------------------------
!

! create a sorted version of the indirect addressing to reduce cache misses

  subroutine get_global_indirect_addressing(nspec,nglob,ibool,copy_ibool_ori,integer_mask_ibool)

  use constants, only: NGLLX,NGLLZ

  implicit none

  integer :: nspec,nglob

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool,copy_ibool_ori
  integer, dimension(nglob) :: integer_mask_ibool

  ! local parameters
  integer :: inumber,ispec,i,j

  ! initializes temporary arrays
  integer_mask_ibool(:) = -1
  copy_ibool_ori(:,:,:) = ibool(:,:,:)

  inumber = 0

  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        if (integer_mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
          ! create a new point
          inumber = inumber + 1
          ibool(i,j,ispec) = inumber
          integer_mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
        else
          ! use an existing point created previously
          ibool(i,j,ispec) = integer_mask_ibool(copy_ibool_ori(i,j,ispec))
        endif
      enddo
    enddo
  enddo

  end subroutine get_global_indirect_addressing

