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

  subroutine add_manual_crack()

  use constants, only: &
    NDIM,IMAIN,ADD_A_SMALL_CRACK_IN_THE_MEDIUM,NB_POINTS_TO_ADD_TO_NPGEO,FAST_NUMBERING

  use specfem_par, only: &
    myrank,NPROC,nspec,kmato,knods,coorg,ngnod,npgeo

  implicit none

  ! local parameters
  integer :: ispec,ispec2,ignod,npgeo_ori,ier
  integer :: npoints_to_add_left,npoints_to_add_right,npoints_to_add
  integer :: current_last_point,original_value
  logical :: already_found_a_crack_element
  double precision,dimension(:,:),allocatable :: tmp_coorg

  ! checks if anything to do
  if (.not. ADD_A_SMALL_CRACK_IN_THE_MEDIUM) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'adding a crack manually'
    write(IMAIN,*) 'need to add ',NB_POINTS_TO_ADD_TO_NPGEO,' npgeo mesh points to do that'
    call flush_IMAIN()
  endif

  ! safety checks
  if (ngnod /= 4) &
    call stop_the_code('must currently have ngnod == 4 when adding a crack manually')
  if (FAST_NUMBERING) &
    call stop_the_code('must not have FAST_NUMBERING when adding a crack manually')

!! DK DK material number 2 indicates the spectral elements that form the left vertical side of the crack, and
!! DK DK material number 3 the right side

  ! all elements with material number 2 and 3 must be in the same slice
  npoints_to_add_left = count(kmato == 2)
  npoints_to_add_right = count(kmato == 3)

  ! both sides must have same number of elements/points to add
  if (npoints_to_add_left /= npoints_to_add_right) then
    print *, 'Error invalid partitioning to add crack:'
    print *, 'rank ',myrank,' has npoints to add left  = ',npoints_to_add_left
    print *, 'rank ',myrank,' has npoints to add right = ',npoints_to_add_right
    print *, 'Number of points must be the same!'
    if (NPROC > 1) print *,'Please try with serial simulation.'
    call stop_the_code('Invalid left and right number of points to add for manual crack')
  endif
  npoints_to_add = npoints_to_add_left

  ! checks if this slice contains the crack, then we must have the right number of points
  if (npoints_to_add > 0 .and. npoints_to_add /= NB_POINTS_TO_ADD_TO_NPGEO) &
      call stop_the_code('must have number of points to add == NB_POINTS_TO_ADD_TO_NPGEO when adding a crack manually')

  ! slice has no elements along the crack, nothing to do
  if (npoints_to_add == 0) return

  ! original points
  npgeo_ori = npgeo

  ! new number of points with crack
  npgeo = npgeo + npoints_to_add

  ! temporary coordinate array, will replace original coorg array
  allocate(tmp_coorg(NDIM,npgeo_ori),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating tmp_coorg array')

  ! initializes
  tmp_coorg(:,:) = coorg(:,:)

  ! replaces original coorg array with larger size one
  deallocate(coorg)

  ! re-sizes coorg
  allocate(coorg(NDIM,npgeo),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating new coorg array')
  coorg(:,:) = 0.d0

  ! copy over from temporary array
  coorg(:,1:npgeo_ori) = tmp_coorg(:,:)

  ! free temporary array
  deallocate(tmp_coorg)

  !! DK DK modify arrays "knods" and "coorg" to introduce the crack manually by duplicating and splitting the nodes
  already_found_a_crack_element = .false.
  current_last_point = npgeo_ori

  do ispec = 1,nspec-1
!! DK DK my convention is to introduce a vertical crack between two elements with material numbers 2 and 3
    if (kmato(ispec) == 2 .and. kmato(ispec+1) == 3) then

      if (myrank == 0) write(IMAIN,*) 'adding a crack between elements ',ispec,' and ',ispec+1

! knods indexing:
!
!  (1,ispec)     (2,ispec)
!   * ----------- *
!   |             |
!   |             |
!   |             |
!   * ----------- *
!  (4,ispec)     (3,ispec)
!
! debug
!if (current_last_point == npgeo_ori) then
!  print *,'knods 1: ',coorg(1,knods(1,ispec)),coorg(2,knods(1,ispec))
!  print *,'knods 2: ',coorg(1,knods(2,ispec)),coorg(2,knods(2,ispec))
!  print *,'knods 3: ',coorg(1,knods(3,ispec)),coorg(2,knods(3,ispec))
!  print *,'knods 4: ',coorg(1,knods(4,ispec)),coorg(2,knods(4,ispec))
!endif
! -> outputs knods 1: 9.18  8.03
!                  2: 9.5   8.03
!                  3: 9.5   8.35
!                  4: 9.18  8.35
!
! note: the following comments might confuse between upper/lower points? see the coordinates output, this indicates:
!       knods(2,..) - upper-right
!       knods(3,..) - lower-right


!! DK DK duplicate and split the lower-right corner of this element,
!! DK DK except if it is the first crack element found, because then it is the crack
!! DK DK tip and thus it should be assembled rather than split.
!! DK DK Lower-right corner of an element is local npgeo point #2
      if (already_found_a_crack_element .and. knods(2,ispec) <= npgeo_ori) then
        current_last_point = current_last_point + 1
        original_value = knods(2,ispec)

!! DK DK split this point number in all the elements in which it appears
        do ispec2 = 1,nspec
! do this only for elements that define the left vertical edge of the crack
          if (kmato(ispec2) /= 2) cycle
          do ignod = 1,ngnod
            if (knods(ignod,ispec2) == original_value) then
              knods(ignod,ispec2) = current_last_point
              coorg(:,current_last_point) = coorg(:,original_value)
            endif
          enddo
        enddo
      endif

!! DK DK duplicate and split the upper-right corner of this element
      already_found_a_crack_element = .true.

!! DK DK Upper-right corner of an element is local npgeo point #3
      if (knods(3,ispec) <= npgeo_ori) then

        current_last_point = current_last_point + 1
        original_value = knods(3,ispec)

!! DK DK split this point number in all the elements in which it appears
        do ispec2 = 1,nspec
! do this only for elements that define the left vertical edge of the crack
          if (kmato(ispec2) /= 2) cycle
          do ignod = 1,ngnod
            if (knods(ignod,ispec2) == original_value) then
              knods(ignod,ispec2) = current_last_point
              coorg(:,current_last_point) = coorg(:,original_value)
            endif
          enddo
        enddo
      endif

    endif ! of if (kmato(ispec) == 2 .and. kmato(ispec+1) == 3)

  enddo ! ispec

  if (current_last_point /= npgeo) then
    print *,'Error invalid number of points added for manual crack:'
    print *,'current_last_point = ',current_last_point
    print *,'npgeo_new = ',npgeo
    call stop_the_code('did not find the right total number of points, should have current_last_point == npgeo_new')
  endif

  end subroutine add_manual_crack


