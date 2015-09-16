
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

! create a sorted version of the indirect addressing to reduce cache misses

  subroutine get_global()

  use specfem_par, only : nspec,ibool,copy_ibool_ori,integer_mask_ibool,SAVE_MODEL,outputname,myrank,ier

  implicit none
  include "constants.h"

  ! local parameters
  integer :: inumber,ispec,i,j

  ! initializes temporary arrays
  integer_mask_ibool(:) = -1
  copy_ibool_ori(:,:,:) = ibool(:,:,:)

  inumber = 0

  ! reduce cache misses in all the elements
  ! loop over spectral elements
    do ispec = 1,nspec
      do j=1,NGLLZ
        do i=1,NGLLX
          if(integer_mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
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

if (save_model=='binary') then

  write(outputname,'(a,i6.6,a)') './DATA/proc',myrank,'_NSPEC_ibool.bin'

  open(888,file=trim(outputname),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error opening smoothed kernel file'
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

  implicit none
  include "constants.h"

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
    do j=1,NGLLZ
      do i=1,NGLLX
        if(integer_mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
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

