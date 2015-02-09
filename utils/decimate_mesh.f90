program subdivide_mesh

  implicit none

  include '../../constants.h'

  ! Number of nodes per elements.
  integer, parameter  :: ESIZE = 4
  ! Max number of neighbours per elements.
  integer,parameter  :: max_neighbour=30
  ! Max number of elements that can contain the same node.
  integer, parameter  :: nsize=20


  integer, parameter  :: NSUB=2

  integer :: nspec
  integer, dimension(:,:), allocatable  :: elmnts
  integer, dimension(:,:), allocatable  :: elmnts_new
  integer, dimension(:), allocatable  :: mat
  integer, dimension(:), allocatable  :: mat_new

  integer :: nnodes
  real, dimension(:,:), allocatable  :: nodes_coords
  real, dimension(:,:), allocatable  :: nodes_coords_new



  real, dimension(2,NSUB+1,NSUB+1)  :: temporary_nodes
  integer, dimension(NSUB+1,NSUB+1)  :: temporary_nodes_lookup

  integer, dimension(:), allocatable  :: xadj
  integer, dimension(:), allocatable  :: adjncy
  integer, dimension(:), allocatable  :: nnodes_elmnts
  integer, dimension(:), allocatable  :: nodes_elmnts


  integer  :: ispec, inode, ispec_neighbours, ispec_neighbours_new
  integer  :: num_nodes_new
  integer  :: i, j
  integer  :: ix, iz

  real :: xtol

  real :: xminval,yminval,xmaxval,ymaxval,xtypdist



  open(unit=98, file='./mesh_14', status='old', form='formatted')
  read(98,*) nspec
  allocate(elmnts(4,nspec))
   do ispec = 1, nspec
     read(98,*) elmnts(1,ispec), elmnts(2,ispec), elmnts(3,ispec), elmnts(4,ispec)
  enddo
  close(98)

  open(unit=98, file='./mat', status='old', form='formatted')
  allocate(mat(nspec))
   do ispec = 1, nspec
     read(98,*) mat(ispec)
  enddo
  close(98)

  open(unit=98, file='./nodes_coords', status='old', form='formatted')
  read(98,*) nnodes
  allocate(nodes_coords(2,nnodes))
  do inode = 1, nnodes
     read(98,*) nodes_coords(1,inode), nodes_coords(2,inode)
  enddo
  close(98)

! set up local geometric tolerances
  xtypdist=+HUGEVAL

  do ispec=1,nspec

  xminval=+HUGEVAL
  yminval=+HUGEVAL
  xmaxval=-HUGEVAL
  ymaxval=-HUGEVAL

  do inode = 1, 4
     xmaxval=max(nodes_coords(1,elmnts(inode,ispec)),xmaxval)
     xminval=min(nodes_coords(1,elmnts(inode,ispec)),xminval)
     ymaxval=max(nodes_coords(2,elmnts(inode,ispec)),ymaxval)
     yminval=min(nodes_coords(2,elmnts(inode,ispec)),yminval)
  enddo

! compute the minimum typical "size" of an element in the mesh
  xtypdist = min(xtypdist,xmaxval-xminval)
  xtypdist = min(xtypdist,ymaxval-yminval)

  enddo

! define a tolerance, small with respect to the minimum size
  xtol=smallvaltol*xtypdist

  print *, 'facteur de tolerance XTOL = ', xtol

  open(unit=98, file='./check', status='unknown', form='formatted')
  do ispec = 1, nspec
     write(98,*) nodes_coords(1,elmnts(1,ispec)), nodes_coords(2,elmnts(1,ispec))
     write(98,*) nodes_coords(1,elmnts(2,ispec)), nodes_coords(2,elmnts(2,ispec))
     write(98,*) nodes_coords(1,elmnts(3,ispec)), nodes_coords(2,elmnts(3,ispec))
     write(98,*) nodes_coords(1,elmnts(4,ispec)), nodes_coords(2,elmnts(4,ispec))
     write(98,*) nodes_coords(1,elmnts(1,ispec)), nodes_coords(2,elmnts(1,ispec))
     write(98,*) ' '
     write(98,*) ' '
  enddo
  close(98)


  allocate(elmnts_new(4,nspec*NSUB*NSUB))
  allocate(mat_new(nspec*NSUB*NSUB))
  allocate(nodes_coords_new(2,nspec*(NSUB+1)*(NSUB+1)))


  elmnts(:,:) = elmnts(:,:) - 1

  allocate(xadj(1:nspec+1))
  allocate(adjncy(1:max_neighbour*nspec))
  allocate(nnodes_elmnts(1:nnodes))
  allocate(nodes_elmnts(1:nsize*nnodes))


  call mesh2dual_ncommonnodes(nspec, nnodes, elmnts, xadj, adjncy, nnodes_elmnts, nodes_elmnts,2)

  elmnts(:,:) = elmnts(:,:) + 1
  adjncy(:) = adjncy(:) + 1
  xadj(:) = xadj(:) + 1

  num_nodes_new = 0

  do ispec = 1, nspec

     do ix = 1, NSUB+1

        temporary_nodes(1,ix,1) = nodes_coords(1,elmnts(1,ispec)) + &
             ( (nodes_coords(1,elmnts(2,ispec)) - nodes_coords(1,elmnts(1,ispec))) / real(NSUB))  * (ix-1)
        temporary_nodes(2,ix,1) = nodes_coords(2,elmnts(1,ispec)) + &
             ( (nodes_coords(2,elmnts(2,ispec)) - nodes_coords(2,elmnts(1,ispec))) / real(NSUB))  * (ix-1)

        temporary_nodes(1,ix,NSUB+1) = nodes_coords(1,elmnts(4,ispec)) + &
             ( (nodes_coords(1,elmnts(3,ispec)) - nodes_coords(1,elmnts(4,ispec))) / real(NSUB))  * (ix-1)
        temporary_nodes(2,ix,NSUB+1) = nodes_coords(2,elmnts(4,ispec)) + &
             ( (nodes_coords(2,elmnts(3,ispec)) - nodes_coords(2,elmnts(4,ispec))) / real(NSUB))  * (ix-1)


        do iz = 2, NSUB

           temporary_nodes(1,ix,iz) =   temporary_nodes(1,ix,1) + &
                (( temporary_nodes(1,ix,NSUB+1) - temporary_nodes(1,ix,1) ) / real(NSUB))  * (iz-1)
           temporary_nodes(2,ix,iz) =   temporary_nodes(2,ix,1) + &
                (( temporary_nodes(2,ix,NSUB+1) - temporary_nodes(2,ix,1) ) / real(NSUB))  * (iz-1)

        enddo

     enddo


     temporary_nodes_lookup(:,:) = 0


     do ispec_neighbours = xadj(ispec), xadj(ispec+1)-1

        if ( adjncy(ispec_neighbours) < ispec ) then
           do ispec_neighbours_new = (adjncy(ispec_neighbours)-1)*NSUB*NSUB + 1, adjncy(ispec_neighbours)*NSUB*NSUB

              do ix = 1, NSUB+1
                 do iz = 1, NSUB+1
                    do inode = 1, 4
                       if ( sqrt( (temporary_nodes(1,ix,iz)-nodes_coords_new(1,elmnts_new(inode,ispec_neighbours_new)))**2 + &
                            (temporary_nodes(2,ix,iz)-nodes_coords_new(2,elmnts_new(inode,ispec_neighbours_new)))**2 ) &
                            < xtol ) then
                          temporary_nodes_lookup(ix,iz) = elmnts_new(inode,ispec_neighbours_new)


                       endif


                    enddo

                 enddo
              enddo
           enddo
        endif
     enddo

  do ix = 1, NSUB+1
     do iz = 1, NSUB+1
        if (temporary_nodes_lookup(ix,iz) == 0 ) then
           num_nodes_new = num_nodes_new + 1
           temporary_nodes_lookup(ix,iz) = num_nodes_new
           nodes_coords_new(1,num_nodes_new) = temporary_nodes(1,ix,iz)
           nodes_coords_new(2,num_nodes_new) = temporary_nodes(2,ix,iz)
        endif
     enddo
  enddo

     do i = 1, NSUB
        do j = 1, NSUB
           elmnts_new(1,(ispec-1)*NSUB*NSUB+(i-1)*NSUB+j) = temporary_nodes_lookup(i,j)
           elmnts_new(2,(ispec-1)*NSUB*NSUB+(i-1)*NSUB+j) = temporary_nodes_lookup(i+1,j)
           elmnts_new(3,(ispec-1)*NSUB*NSUB+(i-1)*NSUB+j) = temporary_nodes_lookup(i+1,j+1)
           elmnts_new(4,(ispec-1)*NSUB*NSUB+(i-1)*NSUB+j) = temporary_nodes_lookup(i,j+1)
           mat_new((ispec-1)*NSUB*NSUB+(i-1)*NSUB+j) = mat(ispec)


        enddo
     enddo


  enddo


  open(unit=99, file='./mesh_new', status='unknown', form='formatted')
  write(99,*) nspec*NSUB*NSUB
  do ispec = 1, nspec*NSUB*NSUB
     write(99,*) elmnts_new(1,ispec), elmnts_new(2,ispec), elmnts_new(3,ispec), elmnts_new(4,ispec)
  enddo
  close(99)

  open(unit=99, file='./mat_new', status='unknown', form='formatted')
  do ispec = 1, nspec*NSUB*NSUB
     write(99,*) mat_new(ispec)
  enddo
  close(99)


  open(unit=99, file='./nodes_coords_new', status='unknown', form='formatted')
  write(99,*) num_nodes_new
  do inode = 1, num_nodes_new
     write(99,*) nodes_coords_new(1,inode), nodes_coords_new(2,inode)
  enddo
  close(99)

  open(unit=99, file='./check_new', status='unknown', form='formatted')
  do ispec = 1, nspec*NSUB*NSUB
     write(99,*) nodes_coords_new(1,elmnts_new(1,ispec)), nodes_coords_new(2,elmnts_new(1,ispec))
     write(99,*) nodes_coords_new(1,elmnts_new(2,ispec)), nodes_coords_new(2,elmnts_new(2,ispec))
     write(99,*) nodes_coords_new(1,elmnts_new(3,ispec)), nodes_coords_new(2,elmnts_new(3,ispec))
     write(99,*) nodes_coords_new(1,elmnts_new(4,ispec)), nodes_coords_new(2,elmnts_new(4,ispec))
  write(99,*) nodes_coords_new(1,elmnts_new(1,ispec)), nodes_coords_new(2,elmnts_new(1,ispec))
     write(99,*) ' '
     write(99,*) ' '
  enddo
  close(99)


end program subdivide_mesh




  !-----------------------------------------------
  ! Creating dual graph (adjacency is defined by 'ncommonnodes' between two elements).
  !-----------------------------------------------
  subroutine mesh2dual_ncommonnodes(nelmnts, nnodes, elmnts, xadj, adjncy, nnodes_elmnts, nodes_elmnts, ncommonnodes)

    ! Number of nodes per elements.
    integer, parameter  :: ESIZE = 4
    ! Max number of neighbours per elements.
    integer,parameter  :: max_neighbour=30
    ! Max number of elements that can contain the same node.
    integer, parameter  :: nsize=20


    integer, intent(in)  :: nelmnts
    integer, intent(in)  :: nnodes
    integer, dimension(0:esize*nelmnts-1), intent(in)  :: elmnts
    integer, dimension(0:nelmnts)  :: xadj
    integer, dimension(0:max_neighbour*nelmnts-1)  :: adjncy
    integer, dimension(0:nnodes-1)  :: nnodes_elmnts
    integer, dimension(0:nsize*nnodes-1)  :: nodes_elmnts
    integer, intent(in)  :: ncommonnodes

    integer  :: i, j, k, l, m, nb_edges
    logical  ::  is_neighbour
    integer  :: num_node, n
    integer  :: elem_base, elem_target
    integer  :: connectivity


    !allocate(xadj(0:nelmnts))
    xadj(:) = 0
    !allocate(adjncy(0:max_neighbour*nelmnts-1))
    adjncy(:) = 0
    !allocate(nnodes_elmnts(0:nnodes-1))
    nnodes_elmnts(:) = 0
    !allocate(nodes_elmnts(0:nsize*nnodes-1))
    nodes_elmnts(:) = 0

    nb_edges = 0


    ! list of elements per node
    do i = 0, esize*nelmnts-1
       nodes_elmnts(elmnts(i)*nsize+nnodes_elmnts(elmnts(i))) = i/esize
       nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1

    enddo

    print *, 'nnodes_elmnts'

    ! checking which elements are neighbours ('ncommonnodes' criteria)
    do j = 0, nnodes-1
       do k = 0, nnodes_elmnts(j)-1
          do l = k+1, nnodes_elmnts(j)-1

             connectivity = 0
             elem_base = nodes_elmnts(k+j*nsize)
             elem_target = nodes_elmnts(l+j*nsize)
             do n = 1, esize
                num_node = elmnts(esize*elem_base+n-1)
                do m = 0, nnodes_elmnts(num_node)-1
                   if ( nodes_elmnts(m+num_node*nsize) == elem_target ) then
                      connectivity = connectivity + 1
                   endif
                enddo
             enddo

             if ( connectivity >=  ncommonnodes) then

                is_neighbour = .false.

                do m = 0, xadj(nodes_elmnts(k+j*nsize))
                   if ( .not.is_neighbour ) then
                      if ( adjncy(nodes_elmnts(k+j*nsize)*max_neighbour+m) == nodes_elmnts(l+j*nsize) ) then
                         is_neighbour = .true.

                      endif
                   endif
                enddo
                if ( .not.is_neighbour ) then
                   adjncy(nodes_elmnts(k+j*nsize)*max_neighbour+xadj(nodes_elmnts(k+j*nsize))) = nodes_elmnts(l+j*nsize)
                   xadj(nodes_elmnts(k+j*nsize)) = xadj(nodes_elmnts(k+j*nsize)) + 1
                   adjncy(nodes_elmnts(l+j*nsize)*max_neighbour+xadj(nodes_elmnts(l+j*nsize))) = nodes_elmnts(k+j*nsize)
                   xadj(nodes_elmnts(l+j*nsize)) = xadj(nodes_elmnts(l+j*nsize)) + 1
                endif
             endif
          enddo
       enddo
    enddo

    ! making adjacency arrays compact (to be used for partitioning)
    do i = 0, nelmnts-1
       k = xadj(i)
       xadj(i) = nb_edges
       do j = 0, k-1
          adjncy(nb_edges) = adjncy(i*max_neighbour+j)
          nb_edges = nb_edges + 1
       enddo
    enddo

    xadj(nelmnts) = nb_edges


  end subroutine mesh2dual_ncommonnodes
