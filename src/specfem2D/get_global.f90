
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
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


  subroutine get_global(nspec_outer,nspec,nglob,ibool)

  implicit none
  include "constants.h"

  integer :: nspec_outer,nspec,nglob

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

  ! local parameters
  integer, dimension(:,:,:), allocatable :: copy_ibool_ori
  integer, dimension(:), allocatable :: mask_ibool
  integer :: inumber,ispec,i,j,ier

  ! allocates temporary arrays
  allocate(mask_ibool(nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocating mask_ibool'
  allocate(copy_ibool_ori(NGLLX,NGLLZ,nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating copy_ibool_ori'

  ! initializes temporary arrays
  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:) = ibool(:,:,:)

  inumber = 0

  if(.not. ACTUALLY_IMPLEMENT_PERM_WHOLE) then

  ! first reduce cache misses in outer elements, since they are taken first
  ! loop over spectral elements
    do ispec = 1,nspec_outer
      do j=1,NGLLZ
        do i=1,NGLLX
          if(mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
            ! create a new point
            inumber = inumber + 1
            ibool(i,j,ispec) = inumber
            mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
          else
            ! use an existing point created previously
            ibool(i,j,ispec) = mask_ibool(copy_ibool_ori(i,j,ispec))
          endif
        enddo
      enddo
    enddo

  ! then reduce cache misses in inner elements, since they are taken second
  ! loop over spectral elements
    do ispec = nspec_outer+1,nspec
      do j=1,NGLLZ
        do i=1,NGLLX
          if(mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
            ! create a new point
            inumber = inumber + 1
            ibool(i,j,ispec) = inumber
            mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
          else
            ! use an existing point created previously
            ibool(i,j,ispec) = mask_ibool(copy_ibool_ori(i,j,ispec))
          endif
        enddo
      enddo
    enddo

  else ! if ACTUALLY_IMPLEMENT_PERM_WHOLE

  ! reduce cache misses in all the elements
  ! loop over spectral elements
    do ispec = 1,nspec
      do j=1,NGLLZ
        do i=1,NGLLX
          if(mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
            ! create a new point
            inumber = inumber + 1
            ibool(i,j,ispec) = inumber
            mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
          else
            ! use an existing point created previously
            ibool(i,j,ispec) = mask_ibool(copy_ibool_ori(i,j,ispec))
          endif
        enddo
      enddo
    enddo

  endif

  deallocate(mask_ibool,copy_ibool_ori)

  end subroutine get_global


!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_global_indirect_addressing(nspec,nglob,ibool)

!- we can create a new indirect addressing to reduce cache misses

  implicit none
  include "constants.h"

  integer :: nspec,nglob
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

  ! local parameters
  integer, dimension(:,:,:), allocatable :: copy_ibool_ori
  integer, dimension(:), allocatable :: mask_ibool
  integer :: inumber,ispec,i,j,ier

  ! allocates temporary arrays
  allocate(mask_ibool(nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocating mask_ibool'
  allocate(copy_ibool_ori(NGLLX,NGLLZ,nspec),stat=ier)
  if( ier /= 0 ) stop 'error allocating copy_ibool_ori'

  ! initializes temporary arrays
  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:) = ibool(:,:,:)

  inumber = 0

  do ispec = 1,nspec
    do j=1,NGLLZ
      do i=1,NGLLX
        if(mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
          ! create a new point
          inumber = inumber + 1
          ibool(i,j,ispec) = inumber
          mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
        else
          ! use an existing point created previously
          ibool(i,j,ispec) = mask_ibool(copy_ibool_ori(i,j,ispec))
        endif
      enddo
    enddo
  enddo

  deallocate(mask_ibool,copy_ibool_ori)

  end subroutine get_global_indirect_addressing
