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
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

  subroutine add_manual_crack(npgeo_ori)

  use constants, only: NB_POINTS_TO_ADD_TO_NPGEO,FAST_NUMBERING
  use specfem_par

  implicit none

  integer,intent(in) :: npgeo_ori

  ! local parameters
  integer :: ispec,ispec2,ignod
  integer :: check_nb_points_to_add_to_npgeo,current_last_point,original_value
  logical :: already_found_a_crack_element

  ! safety check
  if (NPROC > 1) stop 'currently only serial runs are handled when adding a crack manually'

!! DK DK material number 2 indicates the spectral elements that form the left vertical side of the crack
  check_nb_points_to_add_to_npgeo = count(kmato == 2)
  print *
  print *,'adding a crack manually'
  print *,'need to add ',nb_points_to_add_to_npgeo,' npgeo mesh points to do that'

  if (check_nb_points_to_add_to_npgeo /= NB_POINTS_TO_ADD_TO_NPGEO) &
    stop 'must have check_nb_points_to_add_to_npgeo == NB_POINTS_TO_ADD_TO_NPGEO when adding a crack manually'

  if (ngnod /= 4) stop 'must currently have ngnod == 4 when adding a crack manually'

  if (FAST_NUMBERING) stop 'must not have FAST_NUMBERING when adding a crack manually'

  !! DK DK modify arrays "knods" and "coorg" to introduce the crack manually by duplicating and splitting the nodes
  already_found_a_crack_element = .false.
  current_last_point = npgeo_ori

  do ispec = 1,nspec-1
!! DK DK my convention is to introduce a vertical crack between two elements with material numbers 2 and 3
    if (kmato(ispec) == 2 .and. kmato(ispec+1) == 3) then

      print *,'adding a crack between elements ',ispec,' and ',ispec+1

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
    print *,'current_last_point = ',current_last_point
    print *,'npgeo_new = ',npgeo
    stop 'did not find the right total number of points, should have current_last_point == npgeo_new'
  endif

  end subroutine add_manual_crack


