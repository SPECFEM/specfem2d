
program add_bottom_PML

! Dimitri Komatitsch, CNRS Marseille, Nov 2012: add a PML layer to the existing CUBIT salt dome mesh

implicit none

integer, parameter :: npoin = 12715
integer, parameter :: nspec = 12487

integer, parameter :: nPML = 3 !!!! add three PML layers around
double precision, parameter :: thickness_PML = 50. !!! thickness of the PML layer
integer, parameter :: nspec_to_add = 127
integer, parameter :: npoin_to_add = nspec_to_add + 1

! the salt dome CUBIT mesh consists of 4-node elements
integer, parameter :: ngnod = 4

integer :: i,npoin_read,nspec_read,i2,i3,ispec,ispec2,inode
double precision :: xmin
logical :: element_found,first_point_found,second_point_found

double precision, dimension(npoin) :: x,z
double precision, dimension(npoin_to_add) :: x_bottom_edge,z_bottom_edge
double precision, dimension(npoin_to_add) :: x_bottom_edge_sorted,z_bottom_edge_sorted
integer, dimension(npoin_to_add) :: ibool_bottom_edge,ibool_bottom_edge_copy,ibool_bottom_edge_sorted
integer, dimension(nspec) :: imat
integer, dimension(ngnod,nspec) :: ibool

open(unit=27,file='modelY1_nodes_coords_file',status='old')
read(27,*) npoin_read
!!!!!!!if (npoin_read /= npoin) stop 'incorrect value of npoin'
do i = 1,npoin
 read(27,*) x(i),z(i)
enddo
close(27)

open(unit=27,file='modelY1_materials_file',status='old')
do ispec = 1,nspec
 read(27,*) imat(ispec)
enddo
close(27)

open(unit=27,file='modelY1_mesh_file',status='old')
read(27,*) nspec_read
!!!!!!!if (nspec_read /= nspec) stop 'incorrect value of nspec'
do ispec = 1,nspec
 read(27,*) ibool(1,ispec),ibool(2,ispec),ibool(3,ispec),ibool(4,ispec)
enddo
close(27)

open(unit=27,file='bottom_edge',status='old')
do i = 1,npoin_to_add
 read(27,*) ibool_bottom_edge(i)
 x_bottom_edge(i) = x(ibool_bottom_edge(i))
 z_bottom_edge(i) = z(ibool_bottom_edge(i))
enddo
close(27)

print *,'minval maxval x_bottom_edge = ',minval(x_bottom_edge),maxval(x_bottom_edge)
print *,'minval maxval z_bottom_edge = ',minval(z_bottom_edge),maxval(z_bottom_edge)

! loop on the x value from smallest to highest (inefficient bubble sort)
ibool_bottom_edge_copy(:) = ibool_bottom_edge(:)
do i2 = 1,npoin_to_add
xmin = + 100000000000000.d0
do i = 1,npoin_to_add
 if (ibool_bottom_edge_copy(i) /= -1) then ! not already sorted
   if (x_bottom_edge(i) < xmin) then
     xmin = x_bottom_edge(i)
     i3 = i
   endif
 endif
enddo
 x_bottom_edge_sorted(i2) = x_bottom_edge(i3)
 z_bottom_edge_sorted(i2) = z_bottom_edge(i3)
 ibool_bottom_edge_sorted(i2) = ibool_bottom_edge(i3)
 ibool_bottom_edge_copy(i3) = -1 !! mark as already used

print *,i2,x_bottom_edge_sorted(i2),z_bottom_edge_sorted(i2)
enddo

print *
print *,'in point list change npoin from ',npoin,' to ',npoin + npoin_to_add*nPML,' and add this to the list:'
do i = 1,NPML
do i2 = 1,npoin_to_add
print *,x_bottom_edge_sorted(i2),z_bottom_edge_sorted(i2) - thickness_PML*i
enddo
enddo
print *

print *
print *,'in element list change nspec from ',nspec,' to ',nspec + nspec_to_add*nPML,' and add this to the list:'
do i = 1,NPML
do ispec2 = 1,nspec_to_add
if (i == 1) then !! link first PML layer to the existing mesh
  print *,ibool_bottom_edge_sorted(ispec2), npoin + ispec2, npoin + ispec2 + 1,ibool_bottom_edge_sorted(ispec2+1)
else !! add nPML-1 more PML layers
  print *,npoin + ispec2 + (i-2)*npoin_to_add, npoin + ispec2 + (i-1)*npoin_to_add, &
          npoin + ispec2 + (i-1)*npoin_to_add + 1, npoin + ispec2 + (i-2)*npoin_to_add + 1
endif
enddo
enddo
print *

! flags for Stacey absorbing boundaries
! integer, parameter :: IBOTTOM = 1
! integer, parameter :: IRIGHT = 2
! integer, parameter :: ITOP = 3
! integer, parameter :: ILEFT = 4
print *
print *,'in absorbing edge list, replace all lines that end with '\0' 1'\0' with this:'
ispec = nspec
do i = 1,NPML
do ispec2 = 1,nspec_to_add
  ispec = ispec + 1
if (i == NPML) then
  print *,ispec,' 2 ',npoin + ispec2 + (i-1)*npoin_to_add, npoin + ispec2 + (i-1)*npoin_to_add + 1,' 1'
endif
enddo
enddo
print *

print *
print *,'in absorbing edge list, add these lines for new right absorbing edges:'
ispec = nspec
do i = 1,NPML
do ispec2 = 1,nspec_to_add
  ispec = ispec + 1
if (ispec2 == nspec_to_add) then
if (i == 1) then !! link first PML layer to the existing mesh
  print *,ispec,' 2 ', npoin + ispec2 + 1,ibool_bottom_edge_sorted(ispec2+1), ' 2'
else !! add nPML-1 more PML layers
  print *,ispec,' 2 ', npoin + ispec2 + (i-1)*npoin_to_add + 1, npoin + ispec2 + (i-2)*npoin_to_add + 1, ' 2'
endif
endif
enddo
enddo
print *

print *
print *,'in absorbing edge list, add these lines for new left absorbing edges:'
ispec = nspec
do i = 1,NPML
do ispec2 = 1,nspec_to_add
  ispec = ispec + 1
if (ispec2 == 1) then
if (i == 1) then !! link first PML layer to the existing mesh
  print *,ispec,' 2 ', ibool_bottom_edge_sorted(ispec2), npoin + ispec2, ' 4'
else !! add nPML-1 more PML layers
  print *,ispec,' 2 ', npoin + ispec2 + (i-2)*npoin_to_add, npoin + ispec2 + (i-1)*npoin_to_add, ' 4'
endif
endif
enddo
enddo
print *

print *
print *,'in list of material numbers, add this to the list:'
do i = 1,NPML
do ispec2 = 1,nspec_to_add

! detect to which spectral element this PML is added (inefficient loop on all the elements)
  element_found = .false.
  do ispec = 1,nspec
    first_point_found = .false.
    second_point_found = .false.
    do inode = 1,ngnod
      if (ibool(inode,ispec) == ibool_bottom_edge_sorted(ispec2)) first_point_found = .true.
    enddo
    do inode = 1,ngnod
      if (ibool(inode,ispec) == ibool_bottom_edge_sorted(ispec2+1)) second_point_found = .true.
    enddo
    if (first_point_found .and. second_point_found) goto 777
  enddo
  stop 'error: spectral element to which PML should be added never found'
777 continue
  print *,imat(ispec)
enddo
enddo
print *

! CPML Flag Meaning
! 1 element belongs to a left CPML layer only
! 2 element belongs to a right CPML layer only
! 3 element belongs to a bottom CPML layer
! 4 element belongs to a top CPML layer only
! 5 element belongs to a top-left CPML corner
! 6 element belongs to a top-right CPML corner
! 7 element belongs to a bottom-left CPML corner
! 8 element belongs to a bottom-right CPML corner
print *
print *,'in file Elements_CPML_list, add this:'
ispec = nspec
do i = 1,NPML
do ispec2 = 1,nspec_to_add
  ispec = ispec + 1
 if (ispec2 <= NPML) then  ! code for a bottom-left CPML corner
   print *,ispec,7
 else if (ispec2 >= nspec_to_add-NPML+1) then  ! code for a bottom-right CPML corner
   print *,ispec,8
 else
  print *,ispec,3  ! code for a pure CPML layer outside a corner
 endif
enddo
enddo
print *

end program add_bottom_PML

