
 program rotate_mesh

 implicit none

!! DK DK Dec 2011: rotate mesh elements to make sure topological absorbing surface is aligned with physical absorbing surface

integer, parameter :: ngnod = 9

integer :: nspec, nglob, nabs

integer :: i,ispec,iglob,idummy,i1,i2,inode,iswap

logical :: found_this_point

integer, dimension(:,:), allocatable :: ibool,ibool_rotated

integer, dimension(:), allocatable :: ielem,iglob1,iglob2

open(unit=12,file='Mesh_ht_test1',status='old')
read(12,*) nspec
allocate(ibool(ngnod,nspec))
allocate(ibool_rotated(ngnod,nspec))
do ispec = 1,nspec
  read(12,*) ibool(1,ispec), ibool(2,ispec), ibool(3,ispec), ibool(4,ispec), &
                  ibool(5,ispec), ibool(6,ispec), ibool(7,ispec), ibool(8,ispec), ibool(9,ispec)
enddo
close(12)

!!!!!!!!!!! rotate for TOP surface !!!!!!!!!!!!!

! first make a copy of ibool (for all elements that will be unchanged)
  ibool_rotated(:,:) = ibool(:,:)

open(unit=12,file='Surf_top_ht_test1',status='old')
read(12,*) nabs
allocate(ielem(nabs))
allocate(iglob1(nabs))
allocate(iglob2(nabs))
do i = 1,nabs
  read(12,*) ielem(i),idummy,iglob1(i),iglob2(i)
enddo
close(12)

! look for rotations needed
do i = 1,nabs
  ispec = ielem(i)

  found_this_point = .false.
  do inode = 1,4
    if(ibool(inode,ispec) == iglob1(i)) then
      i1 = inode
      found_this_point = .true.
      exit
    endif
  enddo
  if(.not. found_this_point) stop 'point not found'

  found_this_point = .false.
  do inode = 1,4
    if(ibool(inode,ispec) == iglob2(i)) then
      i2 = inode
      found_this_point = .true.
      exit
    endif
  enddo
  if(.not. found_this_point) stop 'point not found'

! swap points if needed for clarity, to avoid testing both cases each time below
  if(i1 > i2) then
    iswap = i1
    i1 = i2
    i2 = iswap
  endif

! test orientation
  if(i1 == 3 .and. i2 == 4) then
    print *,'orientation of element ',i,' is already good'

  else if (i1 == 1 .and. i2 == 4) then ! for this one, remember that we have swapped, thus 4 1 is 1 4
    print *,'element ',i,' must be rotated 3 times'
    ibool_rotated(4,ispec) = ibool(1,ispec)
    ibool_rotated(1,ispec) = ibool(2,ispec)
    ibool_rotated(2,ispec) = ibool(3,ispec)
    ibool_rotated(3,ispec) = ibool(4,ispec)
    ibool_rotated(8,ispec) = ibool(5,ispec)
    ibool_rotated(5,ispec) = ibool(6,ispec)
    ibool_rotated(6,ispec) = ibool(7,ispec)
    ibool_rotated(7,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element

  else if (i1 == 1 .and. i2 == 2) then
    print *,'element ',i,ispec,' must be rotated 2 times top'
    ibool_rotated(3,ispec) = ibool(1,ispec)
    ibool_rotated(4,ispec) = ibool(2,ispec)
    ibool_rotated(1,ispec) = ibool(3,ispec)
    ibool_rotated(2,ispec) = ibool(4,ispec)
    ibool_rotated(7,ispec) = ibool(5,ispec)
    ibool_rotated(8,ispec) = ibool(6,ispec)
    ibool_rotated(5,ispec) = ibool(7,ispec)
    ibool_rotated(6,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element

  else if (i1 == 2 .and. i2 == 3) then
    print *,'element ',i,' must be rotated 1 time'
    ibool_rotated(2,ispec) = ibool(1,ispec)
    ibool_rotated(3,ispec) = ibool(2,ispec)
    ibool_rotated(4,ispec) = ibool(3,ispec)
    ibool_rotated(1,ispec) = ibool(4,ispec)
    ibool_rotated(6,ispec) = ibool(5,ispec)
    ibool_rotated(7,ispec) = ibool(6,ispec)
    ibool_rotated(8,ispec) = ibool(7,ispec)
    ibool_rotated(5,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element

  else
    stop 'problem in an element'
  endif

enddo

deallocate(ielem)
deallocate(iglob1)
deallocate(iglob2)

!!!!!!!!!!! rotate for BOTTOM surface !!!!!!!!!!!!!

! copy the changes above back to ibool, which is used again below
  ibool(:,:) = ibool_rotated(:,:)

open(unit=12,file='Surf_bottom_ht_test1',status='old')
read(12,*) nabs
allocate(ielem(nabs))
allocate(iglob1(nabs))
allocate(iglob2(nabs))
do i = 1,nabs
  read(12,*) ielem(i),idummy,iglob1(i),iglob2(i)
enddo
close(12)

! look for rotations needed
do i = 1,nabs
  ispec = ielem(i)

  found_this_point = .false.
  do inode = 1,4
    if(ibool(inode,ispec) == iglob1(i)) then
      i1 = inode
      found_this_point = .true.
      exit
    endif
  enddo
  if(.not. found_this_point) stop 'point not found'

  found_this_point = .false.
  do inode = 1,4
    if(ibool(inode,ispec) == iglob2(i)) then
      i2 = inode
      found_this_point = .true.
      exit
    endif
  enddo
  if(.not. found_this_point) stop 'point not found'

! swap points if needed for clarity, to avoid testing both cases each time below
  if(i1 > i2) then
    iswap = i1
    i1 = i2
    i2 = iswap
  endif

! test orientation
  if(i1 == 1 .and. i2 == 2) then
    print *,'orientation of element ',i,' is already good'

  else if (i1 == 2 .and. i2 == 3) then
    print *,'element ',i,' must be rotated 3 times'
    ibool_rotated(4,ispec) = ibool(1,ispec)
    ibool_rotated(1,ispec) = ibool(2,ispec)
    ibool_rotated(2,ispec) = ibool(3,ispec)
    ibool_rotated(3,ispec) = ibool(4,ispec)
    ibool_rotated(8,ispec) = ibool(5,ispec)
    ibool_rotated(5,ispec) = ibool(6,ispec)
    ibool_rotated(6,ispec) = ibool(7,ispec)
    ibool_rotated(7,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element

  else if (i1 == 3 .and. i2 == 4) then
    print *,'element ',i,ispec,' must be rotated 2 times bottom'
    ibool_rotated(3,ispec) = ibool(1,ispec)
    ibool_rotated(4,ispec) = ibool(2,ispec)
    ibool_rotated(1,ispec) = ibool(3,ispec)
    ibool_rotated(2,ispec) = ibool(4,ispec)
    ibool_rotated(7,ispec) = ibool(5,ispec)
    ibool_rotated(8,ispec) = ibool(6,ispec)
    ibool_rotated(5,ispec) = ibool(7,ispec)
    ibool_rotated(6,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element

  else if (i1 == 1 .and. i2 == 4) then ! for this one, remember that we have swapped, thus 4 1 is 1 4
    print *,'element ',i,' must be rotated 1 time'
    ibool_rotated(2,ispec) = ibool(1,ispec)
    ibool_rotated(3,ispec) = ibool(2,ispec)
    ibool_rotated(4,ispec) = ibool(3,ispec)
    ibool_rotated(1,ispec) = ibool(4,ispec)
    ibool_rotated(6,ispec) = ibool(5,ispec)
    ibool_rotated(7,ispec) = ibool(6,ispec)
    ibool_rotated(8,ispec) = ibool(7,ispec)
    ibool_rotated(5,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element

  else
    stop 'problem in an element'
  endif

enddo

deallocate(ielem)
deallocate(iglob1)
deallocate(iglob2)

!!!!!!!!!!! rotate for LEFT surface !!!!!!!!!!!!!

! copy the changes above back to ibool, which is used again below
  ibool(:,:) = ibool_rotated(:,:)

open(unit=12,file='Surf_left_ht_test1',status='old')
read(12,*) nabs
allocate(ielem(nabs))
allocate(iglob1(nabs))
allocate(iglob2(nabs))
do i = 1,nabs
  read(12,*) ielem(i),idummy,iglob1(i),iglob2(i)
enddo
close(12)

! look for rotations needed
do i = 1,nabs
  ispec = ielem(i)

  found_this_point = .false.
  do inode = 1,4
    if(ibool(inode,ispec) == iglob1(i)) then
      i1 = inode
      found_this_point = .true.
      exit
    endif
  enddo
  if(.not. found_this_point) stop 'point not found'

  found_this_point = .false.
  do inode = 1,4
    if(ibool(inode,ispec) == iglob2(i)) then
      i2 = inode
      found_this_point = .true.
      exit
    endif
  enddo
  if(.not. found_this_point) stop 'point not found'

! swap points if needed for clarity, to avoid testing both cases each time below
  if(i1 > i2) then
    iswap = i1
    i1 = i2
    i2 = iswap
  endif

! test orientation
  if(i1 == 1 .and. i2 == 4) then ! for this one, remember that we have swapped, thus 4 1 is 1 4
    print *,'orientation of element ',i,' is already good'

  else if (i1 == 1 .and. i2 == 2) then
    print *,'element ',i,' must be rotated 3 times'
    ibool_rotated(4,ispec) = ibool(1,ispec)
    ibool_rotated(1,ispec) = ibool(2,ispec)
    ibool_rotated(2,ispec) = ibool(3,ispec)
    ibool_rotated(3,ispec) = ibool(4,ispec)
    ibool_rotated(8,ispec) = ibool(5,ispec)
    ibool_rotated(5,ispec) = ibool(6,ispec)
    ibool_rotated(6,ispec) = ibool(7,ispec)
    ibool_rotated(7,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element

  else if (i1 == 2 .and. i2 == 3) then
    print *,'element ',i,ispec,' must be rotated 2 times left'
    ibool_rotated(3,ispec) = ibool(1,ispec)
    ibool_rotated(4,ispec) = ibool(2,ispec)
    ibool_rotated(1,ispec) = ibool(3,ispec)
    ibool_rotated(2,ispec) = ibool(4,ispec)
    ibool_rotated(7,ispec) = ibool(5,ispec)
    ibool_rotated(8,ispec) = ibool(6,ispec)
    ibool_rotated(5,ispec) = ibool(7,ispec)
    ibool_rotated(6,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element

  else if (i1 == 3 .and. i2 == 4) then
    print *,'element ',i,' must be rotated 1 time'
    ibool_rotated(2,ispec) = ibool(1,ispec)
    ibool_rotated(3,ispec) = ibool(2,ispec)
    ibool_rotated(4,ispec) = ibool(3,ispec)
    ibool_rotated(1,ispec) = ibool(4,ispec)
    ibool_rotated(6,ispec) = ibool(5,ispec)
    ibool_rotated(7,ispec) = ibool(6,ispec)
    ibool_rotated(8,ispec) = ibool(7,ispec)
    ibool_rotated(5,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element

  else
    stop 'problem in an element'
  endif

enddo

deallocate(ielem)
deallocate(iglob1)
deallocate(iglob2)

!!!!!!!!!!! rotate for RIGHT surface !!!!!!!!!!!!!

! copy the changes above back to ibool, which is used again below
  ibool(:,:) = ibool_rotated(:,:)

open(unit=12,file='Surf_right_ht_test1',status='old')
read(12,*) nabs
allocate(ielem(nabs))
allocate(iglob1(nabs))
allocate(iglob2(nabs))
do i = 1,nabs
  read(12,*) ielem(i),idummy,iglob1(i),iglob2(i)
enddo
close(12)

! look for rotations needed
do i = 1,nabs
  ispec = ielem(i)

  found_this_point = .false.
  do inode = 1,4
    if(ibool(inode,ispec) == iglob1(i)) then
      i1 = inode
      found_this_point = .true.
      exit
    endif
  enddo
  if(.not. found_this_point) stop 'point not found'

  found_this_point = .false.
  do inode = 1,4
    if(ibool(inode,ispec) == iglob2(i)) then
      i2 = inode
      found_this_point = .true.
      exit
    endif
  enddo
  if(.not. found_this_point) stop 'point not found'

! swap points if needed for clarity, to avoid testing both cases each time below
  if(i1 > i2) then
    iswap = i1
    i1 = i2
    i2 = iswap
  endif

! test orientation
  if(i1 == 2 .and. i2 == 3) then
    print *,'orientation of element ',i,' is already good'

  else if (i1 == 3 .and. i2 == 4) then
    print *,'element ',i,' must be rotated 3 times'
    ibool_rotated(4,ispec) = ibool(1,ispec)
    ibool_rotated(1,ispec) = ibool(2,ispec)
    ibool_rotated(2,ispec) = ibool(3,ispec)
    ibool_rotated(3,ispec) = ibool(4,ispec)
    ibool_rotated(8,ispec) = ibool(5,ispec)
    ibool_rotated(5,ispec) = ibool(6,ispec)
    ibool_rotated(6,ispec) = ibool(7,ispec)
    ibool_rotated(7,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element

  else if (i1 == 1 .and. i2 == 4) then ! for this one, remember that we have swapped, thus 4 1 is 1 4
    print *,'element ',i,ispec,' must be rotated 2 times right'
    ibool_rotated(3,ispec) = ibool(1,ispec)
    ibool_rotated(4,ispec) = ibool(2,ispec)
    ibool_rotated(1,ispec) = ibool(3,ispec)
    ibool_rotated(2,ispec) = ibool(4,ispec)
    ibool_rotated(7,ispec) = ibool(5,ispec)
    ibool_rotated(8,ispec) = ibool(6,ispec)
    ibool_rotated(5,ispec) = ibool(7,ispec)
    ibool_rotated(6,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element

  else if (i1 == 1 .and. i2 == 2) then
    print *,'element ',i,' must be rotated 1 time'
    ibool_rotated(2,ispec) = ibool(1,ispec)
    ibool_rotated(3,ispec) = ibool(2,ispec)
    ibool_rotated(4,ispec) = ibool(3,ispec)
    ibool_rotated(1,ispec) = ibool(4,ispec)
    ibool_rotated(6,ispec) = ibool(5,ispec)
    ibool_rotated(7,ispec) = ibool(6,ispec)
    ibool_rotated(8,ispec) = ibool(7,ispec)
    ibool_rotated(5,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element

  else
    stop 'problem in an element'
  endif

enddo

deallocate(ielem)
deallocate(iglob1)
deallocate(iglob2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! save the rotated mesh file
open(unit=12,file='Mesh_ht_test1_rotated',status='unknown')
write(12,*) nspec
do ispec = 1,nspec
  write(12,*) ibool_rotated(1,ispec), ibool_rotated(2,ispec), ibool_rotated(3,ispec), ibool_rotated(4,ispec), &
            ibool_rotated(5,ispec), ibool_rotated(6,ispec), ibool_rotated(7,ispec), ibool_rotated(8,ispec), ibool_rotated(9,ispec)
enddo
close(12)

end

