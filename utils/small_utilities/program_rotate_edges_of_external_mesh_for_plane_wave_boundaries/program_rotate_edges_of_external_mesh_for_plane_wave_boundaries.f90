 program djjkfdkjdf

implicit none

integer, parameter :: ngnod = 4 !! DK DK for now, will switch to 9 later

double precision, parameter :: x_left_edge   = -1.d0 !! DK DK for now
double precision, parameter :: x_right_edge  = +1.d0 !! DK DK for now

double precision, parameter :: tol =  (x_right_edge - x_left_edge) / 1.d6

double precision, dimension(:), allocatable :: x,z
integer, dimension(:,:), allocatable :: ibool

integer :: npoin,nspec,ipoin,ispec,counter

! ========================== first process the left edge ==========================

! header
open(unit=12,file='Nodes_SqrCircMeta.msh',status='old')

read(12,*) npoin
allocate(x(npoin))
allocate(z(npoin))

do ipoin = 1,npoin
  read(12,*) x(ipoin),z(ipoin)
enddo
close(12)

! header
open(unit=12,file='Mesh_SqrCircMeta.msh',status='old')

read(12,*) nspec
allocate(ibool(ngnod,nspec))

do ispec = 1,nspec
  read(12,*) ibool(1,ispec),ibool(2,ispec),ibool(3,ispec),ibool(4,ispec)
enddo

close(12)

open(unit=12,file='Mesh_SqrCircMeta.msh',status='unknown')

write(12,*) nspec

! find left edge
counter = 0
do ispec = 1,nspec

  if (abs(x(ibool(1,ispec)) - x_left_edge) < tol .and. abs(x(ibool(2,ispec)) - x_left_edge) < tol) then
    print *,'elem ',ispec,' is on left edge with edge 1'
    counter = counter + 1
    write(12,*) ibool(2,ispec),ibool(3,ispec),ibool(4,ispec),ibool(1,ispec)

  else if (abs(x(ibool(2,ispec)) - x_left_edge) < tol .and. abs(x(ibool(3,ispec)) - x_left_edge) < tol) then
    print *,'elem ',ispec,' is on left edge with edge 2'
    counter = counter + 1
    write(12,*) ibool(3,ispec),ibool(4,ispec),ibool(1,ispec),ibool(2,ispec)

  else if (abs(x(ibool(3,ispec)) - x_left_edge) < tol .and. abs(x(ibool(4,ispec)) - x_left_edge) < tol) then
    print *,'elem ',ispec,' is on left edge with edge 3'
    counter = counter + 1
    write(12,*) ibool(4,ispec),ibool(1,ispec),ibool(2,ispec),ibool(3,ispec)

  else if (abs(x(ibool(4,ispec)) - x_left_edge) < tol .and. abs(x(ibool(1,ispec)) - x_left_edge) < tol) then
    print *,'elem ',ispec,' is on left edge with edge 4'
    counter = counter + 1
    write(12,*) ibool(1,ispec),ibool(2,ispec),ibool(3,ispec),ibool(4,ispec)

  else
    write(12,*) ibool(1,ispec),ibool(2,ispec),ibool(3,ispec),ibool(4,ispec)
  endif

enddo

print *
print *,'found ',counter,' elements on left edge'
print *

close(12)

deallocate(x,z,ibool)

! ========================== then process the right edge ==========================

! header
open(unit=12,file='Nodes_SqrCircMeta.msh',status='old')

read(12,*) npoin
allocate(x(npoin))
allocate(z(npoin))

do ipoin = 1,npoin
  read(12,*) x(ipoin),z(ipoin)
enddo
close(12)

! header
open(unit=12,file='Mesh_SqrCircMeta.msh',status='old')

read(12,*) nspec
allocate(ibool(ngnod,nspec))

do ispec = 1,nspec
  read(12,*) ibool(1,ispec),ibool(2,ispec),ibool(3,ispec),ibool(4,ispec)
enddo

close(12)

open(unit=12,file='Mesh_SqrCircMeta.msh',status='unknown')

write(12,*) nspec

! find right edge
counter = 0
do ispec = 1,nspec

  if (abs(x(ibool(1,ispec)) - x_right_edge) < tol .and. abs(x(ibool(2,ispec)) - x_right_edge) < tol) then
    print *,'elem ',ispec,' is on right edge with edge 1'
    counter = counter + 1
    write(12,*) ibool(4,ispec),ibool(1,ispec),ibool(2,ispec),ibool(3,ispec)

  else if (abs(x(ibool(2,ispec)) - x_right_edge) < tol .and. abs(x(ibool(3,ispec)) - x_right_edge) < tol) then
    print *,'elem ',ispec,' is on right edge with edge 2'
    counter = counter + 1
    write(12,*) ibool(1,ispec),ibool(2,ispec),ibool(3,ispec),ibool(4,ispec)

  else if (abs(x(ibool(3,ispec)) - x_right_edge) < tol .and. abs(x(ibool(4,ispec)) - x_right_edge) < tol) then
    print *,'elem ',ispec,' is on right edge with edge 3'
    counter = counter + 1
    write(12,*) ibool(2,ispec),ibool(3,ispec),ibool(4,ispec),ibool(1,ispec)

  else if (abs(x(ibool(4,ispec)) - x_right_edge) < tol .and. abs(x(ibool(1,ispec)) - x_right_edge) < tol) then
    print *,'elem ',ispec,' is on right edge with edge 4'
    counter = counter + 1
    write(12,*) ibool(3,ispec),ibool(4,ispec),ibool(1,ispec),ibool(2,ispec)

  else
    write(12,*) ibool(1,ispec),ibool(2,ispec),ibool(3,ispec),ibool(4,ispec)
  endif

enddo

print *
print *,'found ',counter,' elements on right edge'
print *

close(12)

deallocate(x,z,ibool)

end

