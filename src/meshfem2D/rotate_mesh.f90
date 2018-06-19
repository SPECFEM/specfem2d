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

  !-----------------------------------------------
  ! rotate_mesh_for_plane_wave.
  ! rotate mesh elements to make sure topological absorbing surface
  ! is aligned with physical absorbing surface, since this is the surface
  ! that we use to impose the plane wave and Bielak boundary conditions.
  !-----------------------------------------------

  subroutine rotate_mesh_for_plane_wave(ngnod)

  use part_unstruct_par, only: elmnts,nelmnts,nelemabs,abs_surface

  implicit none

  integer, intent(in)  :: ngnod

  ! local parameters
  integer :: i,j,ispec,i1,i2,inode,iswap
  logical :: found_this_point
  integer, dimension(:,:), allocatable :: ibool,ibool_rotated
!! DK DK beware here, "ibool" applies to the mesh corners (4 or 9 points) only,
!! DK DK not to the GLL points because there are no GLL points in the Gmsh mesh files
  integer :: index_rotation1,index_rotation2,index_rotation3,index_rotation4, &
             index_rotation5,index_rotation6,index_rotation7,index_rotation8,index_edge

  allocate(ibool(ngnod,nelmnts))
  allocate(ibool_rotated(ngnod,nelmnts))

 ! At the end of the loop, thank to ibool we can access to the global number of
 ! each node from the ispec of the element to which it belongs and from its
 ! geometrical position :
 !            4 . . 7 . . 3
 !            .           .
 !            .           .
 !            8     9     6
 !            .           .
 !            .           .
 !            1 . . 5 . . 2
 ! --> we just create a copy in an easier format for ease of use in this routine
  do ispec = 1, nelmnts
    if (ngnod == 4) then
      ibool(1,ispec) = elmnts((ispec-1)*ngnod)
      ibool(2,ispec) = elmnts((ispec-1)*ngnod+1)
      ibool(3,ispec) = elmnts((ispec-1)*ngnod+2)
      ibool(4,ispec) = elmnts((ispec-1)*ngnod+3)
    else if (ngnod == 9) then
      ibool(1,ispec) = elmnts((ispec-1)*ngnod)
      ibool(2,ispec) = elmnts((ispec-1)*ngnod+1)
      ibool(3,ispec) = elmnts((ispec-1)*ngnod+2)
      ibool(4,ispec) = elmnts((ispec-1)*ngnod+3)
      ibool(5,ispec) = elmnts((ispec-1)*ngnod+4)
      ibool(6,ispec) = elmnts((ispec-1)*ngnod+5)
      ibool(7,ispec) = elmnts((ispec-1)*ngnod+6)
      ibool(8,ispec) = elmnts((ispec-1)*ngnod+7)
      ibool(9,ispec) = elmnts((ispec-1)*ngnod+8)
    else
      call stop_the_code('error, ngnod should be either 4 or 9 for external meshes')
    endif
  enddo

  do j = 1, 4
    if (j == 1) then
      index_edge=3
      ibool_rotated(:,:) = ibool(:,:)
    else if (j == 2) then
      index_edge=1
      ibool(:,:) = ibool_rotated(:,:)
    else if (j == 3) then
      index_edge=4
      ibool(:,:) = ibool_rotated(:,:)
    else if (j == 4) then
      index_edge=2
      ibool(:,:) = ibool_rotated(:,:)
    else
      call stop_the_code('j should be >= 1 and <= 4')
    endif

    if (index_edge == 1) then
      ! bottom edge
      index_rotation1 = 1
      index_rotation2 = 2
      index_rotation3 = 2
      index_rotation4 = 3
      index_rotation5 = 3
      index_rotation6 = 4
      index_rotation7 = 1
      index_rotation8 = 4
      ! index_rotation9 does not need to exist because the center rotates on itself
    else if (index_edge == 2) then
      ! right edge
      index_rotation1 = 2
      index_rotation2 = 3
      index_rotation3 = 3
      index_rotation4 = 4
      index_rotation5 = 1
      index_rotation6 = 4
      index_rotation7 = 1
      index_rotation8 = 2
      ! index_rotation9 does not need to exist because the center rotates on itself
    else if (index_edge == 3) then
      ! top edge
      index_rotation1 = 3
      index_rotation2 = 4
      index_rotation3 = 1
      index_rotation4 = 4
      index_rotation5 = 1
      index_rotation6 = 2
      index_rotation7 = 2
      index_rotation8 = 3
      ! index_rotation9 does not need to exist because the center rotates on itself
    else if (index_edge == 4) then
      ! left edge
      index_rotation1 = 1
      index_rotation2 = 4
      index_rotation3 = 1
      index_rotation4 = 2
      index_rotation5 = 2
      index_rotation6 = 3
      index_rotation7 = 3
      index_rotation8 = 4
      ! index_rotation9 does not need to exist because the center rotates on itself
    else
      call stop_the_code('The edge on which abs_nodes is located should be defined')
    endif

    do i = 1,nelemabs
      if (index_edge == abs_surface(5,i)) then
        ispec = abs_surface(1,i) + 1  !!!! be careful: ispec from abs_surface(1,i) start at zero
        found_this_point = .false.
        do inode = 1,ngnod
          if (ibool(inode,ispec) == abs_surface(3,i)) then
            i1 = inode
            found_this_point = .true.
            exit
          endif
        enddo

        if (.not. found_this_point) call stop_the_code('point not found')

        found_this_point = .false.
        do inode = 1,4
          if (ibool(inode,ispec) == abs_surface(4,i)) then
            i2 = inode
            found_this_point = .true.
            exit
          endif
        enddo
        if (.not. found_this_point) call stop_the_code('point not found')

        ! swap points if needed for clarity, to avoid testing both cases each time below
        if (i1 > i2) then
          iswap = i1
          i1 = i2
          i2 = iswap
        endif

        ! test orientation
        if (i1 == index_rotation1 .and. i2 == index_rotation2) then
        ! print *,'orientation of element ',i,' is already good'

        else if (i1 == index_rotation3 .and. i2 == index_rotation4) then
          !for this one, remember that we have swapped, thus 41 is 14
          ! print *,'element ',i,' must be rotated 3 times'
          ibool_rotated(4,ispec) = ibool(1,ispec)
          ibool_rotated(1,ispec) = ibool(2,ispec)
          ibool_rotated(2,ispec) = ibool(3,ispec)
          ibool_rotated(3,ispec) = ibool(4,ispec)
          if (ngnod == 9) then
            ibool_rotated(8,ispec) = ibool(5,ispec)
            ibool_rotated(5,ispec) = ibool(6,ispec)
            ibool_rotated(6,ispec) = ibool(7,ispec)
            ibool_rotated(7,ispec) = ibool(8,ispec)
            ! 9th point is at the element center and thus never changes when we rotate an element
          endif

        else if (i1 == index_rotation5 .and. i2 == index_rotation6) then
          ! print *,'element ',i,ispec,' must be rotated 2 times top'
          ibool_rotated(3,ispec) = ibool(1,ispec)
          ibool_rotated(4,ispec) = ibool(2,ispec)
          ibool_rotated(1,ispec) = ibool(3,ispec)
          ibool_rotated(2,ispec) = ibool(4,ispec)
          if (ngnod == 9) then
            ibool_rotated(7,ispec) = ibool(5,ispec)
            ibool_rotated(8,ispec) = ibool(6,ispec)
            ibool_rotated(5,ispec) = ibool(7,ispec)
            ibool_rotated(6,ispec) = ibool(8,ispec)
            ! 9th point is at the element center and thus never changes when we rotate an element
          endif

        else if (i1 == index_rotation7 .and. i2 == index_rotation8) then
          ! print *,'element ',i,' must be rotated 1 time'
          ibool_rotated(2,ispec) = ibool(1,ispec)
          ibool_rotated(3,ispec) = ibool(2,ispec)
          ibool_rotated(4,ispec) = ibool(3,ispec)
          ibool_rotated(1,ispec) = ibool(4,ispec)
          if (ngnod == 9) then
            ibool_rotated(6,ispec) = ibool(5,ispec)
            ibool_rotated(7,ispec) = ibool(6,ispec)
            ibool_rotated(8,ispec) = ibool(7,ispec)
            ibool_rotated(5,ispec) = ibool(8,ispec)
            ! 9th point is at the element center and thus never changes when we rotate an element
          endif
        else
          call stop_the_code('problem in an element')
        endif
      endif
    enddo
  enddo

! here we put the result back in the not-so-easy to use format at the end of the routine
  do ispec = 1, nelmnts
    if (ngnod == 4) then
      elmnts((ispec-1)*ngnod)   = ibool_rotated(1,ispec)
      elmnts((ispec-1)*ngnod+1) = ibool_rotated(2,ispec)
      elmnts((ispec-1)*ngnod+2) = ibool_rotated(3,ispec)
      elmnts((ispec-1)*ngnod+3) = ibool_rotated(4,ispec)
    else if (ngnod == 9) then
      elmnts((ispec-1)*ngnod)   = ibool_rotated(1,ispec)
      elmnts((ispec-1)*ngnod+1) = ibool_rotated(2,ispec)
      elmnts((ispec-1)*ngnod+2) = ibool_rotated(3,ispec)
      elmnts((ispec-1)*ngnod+3) = ibool_rotated(4,ispec)
      elmnts((ispec-1)*ngnod+4) = ibool_rotated(5,ispec)
      elmnts((ispec-1)*ngnod+5) = ibool_rotated(6,ispec)
      elmnts((ispec-1)*ngnod+6) = ibool_rotated(7,ispec)
      elmnts((ispec-1)*ngnod+7) = ibool_rotated(8,ispec)
      elmnts((ispec-1)*ngnod+8) = ibool_rotated(9,ispec)
    else
      call stop_the_code('error, ngnod should be either 4 or 9 for external meshes')
    endif
  enddo

  end subroutine rotate_mesh_for_plane_wave

!
!---------------------------------------------------------------------------------------
!

  subroutine rotate_mesh_for_axisym(ngnod)

! This routine is almost the same than the one above with index_edge = 4 (left side of the model)
! It allows us to know that if i == 1 (first GLL point) we are on the axis.
! From axisymmetric elements mesh file with just know the ispec of the elements on the axis and the ibool of
! the two points that describe it :
! Ex :
!
!   x=0
!   ...
!    |
!    O ibool = 235   ---> We want that index to be saved in ibool(1,ispec) or ibool(4,ispec)
!    |
!    |
!    | ispec = 15
!    |
!    |
!    O ibool = 423   ---> We want that index to be saved in ibool(4,ispec) or ibool(1,ispec)
!    |
!    ...
!
! Indeed when we are using external mesher for axisymmetric simulations we do not control how the elements
! will be orientated. We could have code everything taking that into account but we have preferred
! to rotate the mesh. After that for each element on the axis we have got:
!
!           r=0
!            |
!            4 . . 7 . . 3
!            .           .
!            .           .
!            8     9     6
!            .           .
!            .           .
!            1 . . 5 . . 2
!            |
!           r=0

  use part_unstruct_par, only: elmnts,nelmnts,nelem_on_the_axis,ispec_of_axial_elements, &
    inode1_axial_elements,inode2_axial_elements

  implicit none

  integer, intent(in)  :: ngnod
  integer :: i,ispec,i1,i2,inode,iswap
  logical :: found_this_point
  integer, dimension(:,:), allocatable :: ibool,ibool_rotated
  !! DK DK be careful here, "ibool" applies to the mesh corners (4 or 9 points) only,

  allocate(ibool(ngnod,nelmnts))
  allocate(ibool_rotated(ngnod,nelmnts))

  do ispec = 1, nelmnts ! Loop on the elements
  ! At the end of the loop, thank to ibool we can access to the global number of
  ! each node from the ispec of the element to which it belongs and from its
  ! geometrical position :
  !            4 . . 7 . . 3
  !            .           .
  !            .           .
  !            8     9     6
  !            .           .
  !            .           .
  !            1 . . 5 . . 2
  ! --> we just create a copy in an easier format for ease of use in this routine
    if (ngnod == 4) then
      ibool(1,ispec) = elmnts((ispec-1)*ngnod)    ! Have to be zero if ispec is on the axis
      ibool(2,ispec) = elmnts((ispec-1)*ngnod+1)
      ibool(3,ispec) = elmnts((ispec-1)*ngnod+2)
      ibool(4,ispec) = elmnts((ispec-1)*ngnod+3)  ! Have to be zero if ispec is on the axis
    else if (ngnod == 9) then
      ibool(1,ispec) = elmnts((ispec-1)*ngnod)    ! Have to be zero if ispec is on the axis
      ibool(2,ispec) = elmnts((ispec-1)*ngnod+1)
      ibool(3,ispec) = elmnts((ispec-1)*ngnod+2)
      ibool(4,ispec) = elmnts((ispec-1)*ngnod+3)  ! Have to be zero if ispec is on the axis
      ibool(5,ispec) = elmnts((ispec-1)*ngnod+4)
      ibool(6,ispec) = elmnts((ispec-1)*ngnod+5)
      ibool(7,ispec) = elmnts((ispec-1)*ngnod+6)
      ibool(8,ispec) = elmnts((ispec-1)*ngnod+7)
      ibool(9,ispec) = elmnts((ispec-1)*ngnod+8)
    else
      call stop_the_code('rotate_mesh_for_axisym: error, ngnod should be either 4 or 9 for external meshes')
    endif
  enddo

  ibool_rotated(:,:) = ibool(:,:) ! We make a copy of ibool in ibool_rotated

  ! print *,"    Loop on the elements on the axis... : "
  do i = 1,nelem_on_the_axis ! Loop on the elements on the axis (read from the axisym file)
    ispec = ispec_of_axial_elements(i) + 1 ! ispec_of_axial_elements starts from 0
    found_this_point = .false.
    ! print *,"        Loop on the control points and look for ", inode1_axial_elements(i)
    do inode = 1,4 ! loop on the corners of axial element ispec_of_axial_elements(i) to see if we find inode1_axial_elements(i)
      if (ibool(inode,ispec) == inode1_axial_elements(i)) then
        i1 = inode
        found_this_point = .true.
        exit
      endif
    enddo
    if (.not. found_this_point) call stop_the_code('rotate_mesh_for_axisym: point not found 1')
    found_this_point = .false.
    do inode = 1,4 ! loop on the corners of axial element ispec_of_axial_elements(i) to see if we find inode2_axial_elements(i)
      if (ibool(inode,ispec) == inode2_axial_elements(i)) then
        i2 = inode
        found_this_point = .true.
        exit
      endif
    enddo
    if (.not. found_this_point) call stop_the_code('rotate_mesh_for_axisym: point not found 2')
    if (i1 > i2) then
      ! swap points if needed for clarity, to avoid testing both cases each time below
      ! Otherwise we would have done : if ((i1 == 1 .and. i2 == 4) .or. (i1 == 4 .and. i2 == 1)) then
      !                                   print *,'orientation of element ',i,' is already good'
      !                                 ...
      iswap = i1
      i1 = i2
      i2 = iswap
    endif

    ! print *,"point 1 (with coord 0):",inode1_axial_elements(i)
    ! print *,"point 2 (with coord 0):",inode2_axial_elements(i)
    ! print *,"ibool of this element :"
    ! print *,"  ibool(1,ispec) :",ibool(1,ispec)
    ! print *,"  ibool(2,ispec) :",ibool(2,ispec)
    ! print *,"  ibool(3,ispec) :",ibool(3,ispec)
    ! print *,"  ibool(4,ispec) :",ibool(4,ispec)
    ! print *,"    i1 -->",i1
    ! print *,"    i2 -->",i2

    ! test orientation
    if (i1 == 1 .and. i2 == 4) then ! Orientation of this element is already good'
      ! print *,'orientation of element ',i,' is already good'
    else if (i1 == 1 .and. i2 == 2) then ! Remember that we have swapped, thus no need to test i1 == 1 and i2 == 1
      ! Element must be rotated 3 times
      ! print *,'element ',i,' must be rotated 3 times'
      ibool_rotated(4,ispec) = ibool(1,ispec)
      ibool_rotated(1,ispec) = ibool(2,ispec)
      ibool_rotated(2,ispec) = ibool(3,ispec)
      ibool_rotated(3,ispec) = ibool(4,ispec)
      if (ngnod == 9) then
        ibool_rotated(8,ispec) = ibool(5,ispec)
        ibool_rotated(5,ispec) = ibool(6,ispec)
        ibool_rotated(6,ispec) = ibool(7,ispec)
        ibool_rotated(7,ispec) = ibool(8,ispec)
        ! 9th point is at the element center and thus never changes when we rotate an element
      endif
    else if (i1 == 2 .and. i2 == 3) then ! Element must be rotated 2 times
      ! print *,'element ',i,ispec,' must be rotated 2 times top'
      ibool_rotated(3,ispec) = ibool(1,ispec)
      ibool_rotated(4,ispec) = ibool(2,ispec)
      ibool_rotated(1,ispec) = ibool(3,ispec)
      ibool_rotated(2,ispec) = ibool(4,ispec)
      if (ngnod == 9) then
        ibool_rotated(7,ispec) = ibool(5,ispec)
        ibool_rotated(8,ispec) = ibool(6,ispec)
        ibool_rotated(5,ispec) = ibool(7,ispec)
        ibool_rotated(6,ispec) = ibool(8,ispec)
        ! 9th point is at the element center and thus never changes when we rotate an element
      endif
    else if (i1 == 3 .and. i2 == 4) then ! Element must be rotated 1 time
      ! print *,'element ',i,' must be rotated 1 time'
      ibool_rotated(2,ispec) = ibool(1,ispec)
      ibool_rotated(3,ispec) = ibool(2,ispec)
      ibool_rotated(4,ispec) = ibool(3,ispec)
      ibool_rotated(1,ispec) = ibool(4,ispec)
      if (ngnod == 9) then
        ibool_rotated(6,ispec) = ibool(5,ispec)
        ibool_rotated(7,ispec) = ibool(6,ispec)
        ibool_rotated(8,ispec) = ibool(7,ispec)
        ibool_rotated(5,ispec) = ibool(8,ispec)
        ! 9th point is at the element center and thus never changes when we rotate an element
      endif
    else
      call stop_the_code('rotate_mesh_for_axisym: problem in an element')
    endif
  enddo

  ! Here we put the result back in the not-so-easy to use format at the end of the routine
  do ispec = 1, nelmnts
    if (ngnod == 4) then
      elmnts((ispec-1)*ngnod)   = ibool_rotated(1,ispec)
      elmnts((ispec-1)*ngnod+1) = ibool_rotated(2,ispec)
      elmnts((ispec-1)*ngnod+2) = ibool_rotated(3,ispec)
      elmnts((ispec-1)*ngnod+3) = ibool_rotated(4,ispec)
    else if (ngnod == 9) then
      elmnts((ispec-1)*ngnod)   = ibool_rotated(1,ispec)
      elmnts((ispec-1)*ngnod+1) = ibool_rotated(2,ispec)
      elmnts((ispec-1)*ngnod+2) = ibool_rotated(3,ispec)
      elmnts((ispec-1)*ngnod+3) = ibool_rotated(4,ispec)
      elmnts((ispec-1)*ngnod+4) = ibool_rotated(5,ispec)
      elmnts((ispec-1)*ngnod+5) = ibool_rotated(6,ispec)
      elmnts((ispec-1)*ngnod+6) = ibool_rotated(7,ispec)
      elmnts((ispec-1)*ngnod+7) = ibool_rotated(8,ispec)
      elmnts((ispec-1)*ngnod+8) = ibool_rotated(9,ispec)
    else
      call stop_the_code('rotate_mesh_for_axisym: error, ngnod should be either 4 or 9 for external meshes')
    endif
  enddo

  end subroutine rotate_mesh_for_axisym

!
!---------------------------------------------------------------------------------------
!

  subroutine rotate_mesh_for_acoustic_forcing(ngnod)

! This routine is almost the same than the one at the top but for the acoustic forced elements
! Ex :
!
!   x=0
!   ...
!    |
!    O ibool = 235   ---> We want that index to be saved in ibool(1,ispec) or ibool(4,ispec)
!    |
!    |
!    | ispec = 15
!    |
!    |
!    O ibool = 423   ---> We want that index to be saved in ibool(4,ispec) or ibool(1,ispec)
!    |
!    ...
!
! Indeed when we are using external mesher we do not control how the elements
! will be orientated. We could have code everything taking that into account but we have preferred
! to rotate the mesh. After that for each element on the axis we have got:
!
!           r=0
!            |
!            4 . . 7 . . 3
!            .           .
!            .           .
!            8     9     6
!            .           .
!            .           .
!            1 . . 5 . . 2
!            |
!           r=0

  use part_unstruct_par, only: elmnts,nelmnts,nelemacforcing,acforcing_surface

  implicit none

  integer, intent(in)  :: ngnod

  ! local parameters
  integer :: i,j,ispec,i1,i2,inode,iswap
  logical :: found_this_point
  integer, dimension(:,:), allocatable :: ibool,ibool_rotated
!! DK DK be careful here, "ibool" applies to the mesh corners (4 or 9 points) only,
!! DK DK not to the GLL points because there are no GLL points in the Gmsh mesh files
  integer :: index_rotation1,index_rotation2,index_rotation3,index_rotation4, &
             index_rotation5,index_rotation6,index_rotation7,index_rotation8,index_edge

  allocate(ibool(ngnod,nelmnts))
  allocate(ibool_rotated(ngnod,nelmnts))

 ! At the end of the loop, thank to ibool we can access to the global number of
 ! each node from the ispec of the element to which it belongs and from its
 ! geometrical position :
 !            4 . . 7 . . 3
 !            .           .
 !            .           .
 !            8     9     6
 !            .           .
 !            .           .
 !            1 . . 5 . . 2
 ! --> we just create a copy in an easier format for ease of use in this routine
  do ispec = 1, nelmnts
    if (ngnod == 4) then
      ibool(1,ispec) = elmnts((ispec-1)*ngnod)
      ibool(2,ispec) = elmnts((ispec-1)*ngnod+1)
      ibool(3,ispec) = elmnts((ispec-1)*ngnod+2)
      ibool(4,ispec) = elmnts((ispec-1)*ngnod+3)
    else if (ngnod == 9) then
      ibool(1,ispec) = elmnts((ispec-1)*ngnod)
      ibool(2,ispec) = elmnts((ispec-1)*ngnod+1)
      ibool(3,ispec) = elmnts((ispec-1)*ngnod+2)
      ibool(4,ispec) = elmnts((ispec-1)*ngnod+3)
      ibool(5,ispec) = elmnts((ispec-1)*ngnod+4)
      ibool(6,ispec) = elmnts((ispec-1)*ngnod+5)
      ibool(7,ispec) = elmnts((ispec-1)*ngnod+6)
      ibool(8,ispec) = elmnts((ispec-1)*ngnod+7)
      ibool(9,ispec) = elmnts((ispec-1)*ngnod+8)
    else
      call stop_the_code('error, ngnod should be either 4 or 9 for external meshes')
    endif
  enddo

  do j = 1, 4
    if (j == 1) then
      index_edge=3
      ibool_rotated(:,:) = ibool(:,:)
    else if (j == 2) then
      index_edge=1
      ibool(:,:) = ibool_rotated(:,:)
    else if (j == 3) then
      index_edge=4
      ibool(:,:) = ibool_rotated(:,:)
    else if (j == 4) then
      index_edge=2
      ibool(:,:) = ibool_rotated(:,:)
    else
      call stop_the_code('j should be >= 1 and <= 4')
    endif

    if (index_edge == 1) then
      ! bottom edge
      index_rotation1 = 1
      index_rotation2 = 2
      index_rotation3 = 2
      index_rotation4 = 3
      index_rotation5 = 3
      index_rotation6 = 4
      index_rotation7 = 1
      index_rotation8 = 4
      ! index_rotation9 does not need to exist because the center rotates on itself
    else if (index_edge == 2) then
      ! right edge
      index_rotation1 = 2
      index_rotation2 = 3
      index_rotation3 = 3
      index_rotation4 = 4
      index_rotation5 = 1
      index_rotation6 = 4
      index_rotation7 = 1
      index_rotation8 = 2
      ! index_rotation9 does not need to exist because the center rotates on itself
    else if (index_edge == 3) then
      ! top edge
      index_rotation1 = 3
      index_rotation2 = 4
      index_rotation3 = 1
      index_rotation4 = 4
      index_rotation5 = 1
      index_rotation6 = 2
      index_rotation7 = 2
      index_rotation8 = 3
      ! index_rotation9 does not need to exist because the center rotates on itself
    else if (index_edge == 4) then
      ! left edge
      index_rotation1 = 1
      index_rotation2 = 4
      index_rotation3 = 1
      index_rotation4 = 2
      index_rotation5 = 2
      index_rotation6 = 3
      index_rotation7 = 3
      index_rotation8 = 4
      ! index_rotation9 does not need to exist because the center rotates on itself
    else
      call stop_the_code('The edge on which abs_nodes is located should be defined')
    endif

    do i = 1,nelemacforcing
      if (index_edge == acforcing_surface(5,i)) then
        ispec = acforcing_surface(1,i) + 1  !!!! be careful: ispec from acforcing_surface(1,i) start at zero
        found_this_point = .false.
        do inode = 1,ngnod
          if (ibool(inode,ispec) == acforcing_surface(3,i)) then
            i1 = inode
            found_this_point = .true.
            exit
          endif
        enddo

        if (.not. found_this_point) call stop_the_code('point not found')

        found_this_point = .false.
        do inode = 1,4
          if (ibool(inode,ispec) == acforcing_surface(4,i)) then
            i2 = inode
            found_this_point = .true.
            exit
          endif
        enddo
        if (.not. found_this_point) call stop_the_code('point not found')

        ! swap points if needed for clarity, to avoid testing both cases each time below
        if (i1 > i2) then
          iswap = i1
          i1 = i2
          i2 = iswap
        endif

        ! test orientation
        if (i1 == index_rotation1 .and. i2 == index_rotation2) then
        ! print *,'orientation of element ',i,' is already good'

        else if (i1 == index_rotation3 .and. i2 == index_rotation4) then
          !for this one, remember that we have swapped, thus 41 is 14
          ! print *,'element ',i,' must be rotated 3 times'
          ibool_rotated(4,ispec) = ibool(1,ispec)
          ibool_rotated(1,ispec) = ibool(2,ispec)
          ibool_rotated(2,ispec) = ibool(3,ispec)
          ibool_rotated(3,ispec) = ibool(4,ispec)
          if (ngnod == 9) then
            ibool_rotated(8,ispec) = ibool(5,ispec)
            ibool_rotated(5,ispec) = ibool(6,ispec)
            ibool_rotated(6,ispec) = ibool(7,ispec)
            ibool_rotated(7,ispec) = ibool(8,ispec)
            ! 9th point is at the element center and thus never changes when we rotate an element
          endif

        else if (i1 == index_rotation5 .and. i2 == index_rotation6) then
          ! print *,'element ',i,ispec,' must be rotated 2 times top'
          ibool_rotated(3,ispec) = ibool(1,ispec)
          ibool_rotated(4,ispec) = ibool(2,ispec)
          ibool_rotated(1,ispec) = ibool(3,ispec)
          ibool_rotated(2,ispec) = ibool(4,ispec)
          if (ngnod == 9) then
            ibool_rotated(7,ispec) = ibool(5,ispec)
            ibool_rotated(8,ispec) = ibool(6,ispec)
            ibool_rotated(5,ispec) = ibool(7,ispec)
            ibool_rotated(6,ispec) = ibool(8,ispec)
            ! 9th point is at the element center and thus never changes when we rotate an element
          endif

        else if (i1 == index_rotation7 .and. i2 == index_rotation8) then
          ! print *,'element ',i,' must be rotated 1 time'
          ibool_rotated(2,ispec) = ibool(1,ispec)
          ibool_rotated(3,ispec) = ibool(2,ispec)
          ibool_rotated(4,ispec) = ibool(3,ispec)
          ibool_rotated(1,ispec) = ibool(4,ispec)
          if (ngnod == 9) then
            ibool_rotated(6,ispec) = ibool(5,ispec)
            ibool_rotated(7,ispec) = ibool(6,ispec)
            ibool_rotated(8,ispec) = ibool(7,ispec)
            ibool_rotated(5,ispec) = ibool(8,ispec)
            ! 9th point is at the element center and thus never changes when we rotate an element
          endif
        else
          call stop_the_code('problem in an element')
        endif
      endif
    enddo
  enddo

! here we put the result back in the not-so-easy to use format at the end of the routine
  do ispec = 1, nelmnts
    if (ngnod == 4) then
      elmnts((ispec-1)*ngnod)   = ibool_rotated(1,ispec)
      elmnts((ispec-1)*ngnod+1) = ibool_rotated(2,ispec)
      elmnts((ispec-1)*ngnod+2) = ibool_rotated(3,ispec)
      elmnts((ispec-1)*ngnod+3) = ibool_rotated(4,ispec)
    else if (ngnod == 9) then
      elmnts((ispec-1)*ngnod)   = ibool_rotated(1,ispec)
      elmnts((ispec-1)*ngnod+1) = ibool_rotated(2,ispec)
      elmnts((ispec-1)*ngnod+2) = ibool_rotated(3,ispec)
      elmnts((ispec-1)*ngnod+3) = ibool_rotated(4,ispec)
      elmnts((ispec-1)*ngnod+4) = ibool_rotated(5,ispec)
      elmnts((ispec-1)*ngnod+5) = ibool_rotated(6,ispec)
      elmnts((ispec-1)*ngnod+6) = ibool_rotated(7,ispec)
      elmnts((ispec-1)*ngnod+7) = ibool_rotated(8,ispec)
      elmnts((ispec-1)*ngnod+8) = ibool_rotated(9,ispec)
    else
      call stop_the_code('error, ngnod should be either 4 or 9 for external meshes')
    endif
  enddo

  end subroutine rotate_mesh_for_acoustic_forcing

