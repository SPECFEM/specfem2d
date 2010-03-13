
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

! implement reverse Cuthill-McKee (1969) ordering, introduced in
! E. Cuthill and J. McKee. Reducing the bandwidth of sparse symmetric matrices.
! In Proceedings of the 1969 24th national conference, pages 157-172,
! New-York, New-York, USA, 1969. ACM Press.
! see for instance http://en.wikipedia.org/wiki/Cuthill%E2%80%93McKee_algorithm

  subroutine get_perm(ibool,perm,limit,nspec,nglob)

  implicit none

  include "constants.h"

! input
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

! output
  integer, dimension(nspec) :: perm

! local variables
  integer nspec,nglob_GLL_full

! maximum number of neighbors of a spectral element (in principle, it could be any value)
  integer, parameter :: MAX_NUMBER_OF_NEIGHBORS = 50

! global corner numbers that need to be created
  integer, dimension(nglob) :: global_corner_number

  integer mn(nspec*NGNOD_QUADRANGLE),mp(nspec+1)
  integer, dimension(:), allocatable :: ne,np,adj
  integer xadj(nspec+1)

! arrays to store the permutation and inverse permutation of the Cuthill-McKee algorithm
  integer, dimension(nspec) :: invperm

  logical maskel(nspec)

  integer i,istart,istop,number_of_neighbors

  integer nglob_four_corners_only,nglob

! only count the total size of the array that will be created, or actually create it
  logical count_only
  integer total_size_ne,total_size_adj,limit

!
!-----------------------------------------------------------------------
!
  if(PERFORM_CUTHILL_MCKEE) then

  ! total number of points in the mesh
    nglob_GLL_full = nglob

  !---- call Charbel Farhat's routines
    call form_elt_connectivity_foelco(mn,mp,nspec,global_corner_number,nglob_GLL_full,ibool,nglob_four_corners_only)
    do i=1,nspec
        istart = mp(i)
        istop = mp(i+1) - 1
    enddo

    allocate(np(nglob_four_corners_only+1))
    count_only = .true.
    total_size_ne = 1
    allocate(ne(total_size_ne))
    call form_node_connectivity_fonoco(mn,mp,ne,np,nglob_four_corners_only,nspec,count_only,total_size_ne)
    deallocate(ne)
    allocate(ne(total_size_ne))
    count_only = .false.
    call form_node_connectivity_fonoco(mn,mp,ne,np,nglob_four_corners_only,nspec,count_only,total_size_ne)
    do i=1,nglob_four_corners_only
        istart = np(i)
        istop = np(i+1) - 1
    enddo

    count_only = .true.
    total_size_adj = 1
    allocate(adj(total_size_adj))
    call create_adjacency_table_adjncy(mn,mp,ne,np,adj,xadj,maskel,nspec,nglob_four_corners_only,&
    count_only,total_size_ne,total_size_adj)
    deallocate(adj)
    allocate(adj(total_size_adj))
    count_only = .false.
    call create_adjacency_table_adjncy(mn,mp,ne,np,adj,xadj,maskel,nspec,nglob_four_corners_only,&
    count_only,total_size_ne,total_size_adj)
    do i=1,nspec
        istart = xadj(i)
        istop = xadj(i+1) - 1
        number_of_neighbors = istop-istart+1
        if(number_of_neighbors < 1) stop 'error: your mesh seems to have at least one element not connected to any other'
        if(number_of_neighbors > MAX_NUMBER_OF_NEIGHBORS) stop 'error: your mesh seems to have an unlikely high valence'
    enddo
    deallocate(ne,np)

! call the Cuthill-McKee sorting algorithm
    call cuthill_mckee(adj,xadj,perm,invperm,nspec,total_size_adj,limit)
    deallocate(adj)
  else
! create identity permutation in order to do nothing
    do i=1,nspec
      perm(i) = i
    enddo
  endif

  end subroutine get_perm

!=======================================================================
!
!  Charbel Farhat's FEM topology routines
!
!  Dimitri Komatitsch, February 1996 - Code based on Farhat's original version
!  described in his technical report from 1987
!
!  modified and adapted by Dimitri Komatitsch, May 2006
!
!=======================================================================

  subroutine form_elt_connectivity_foelco(mn,mp,nspec,global_corner_number, &
                      nglob_GLL_full,ibool,nglob_four_corners_only)

!-----------------------------------------------------------------------
!
!   Forms the MN and MP arrays
!
!     Input :
!     -------
!           ibool    Array needed to build the element connectivity table
!           nspec    Number of elements in the domain
!           NGNOD_QUADRANGLE    number of nodes per hexahedron (brick with 8 corners)
!
!     Output :
!     --------
!           MN, MP   This is the element connectivity array pair.
!                    Array MN contains the list of the element
!                    connectivity, that is, the nodes contained in each
!                    element. They are stored in a stacked fashion.
!
!                    Pointer array MP stores the location of each
!                    element list. Its length is equal to the number
!                    of elements plus one.
!
!-----------------------------------------------------------------------

  implicit none

  include "constants.h"

  integer nspec,nglob_GLL_full

! arrays with mesh parameters per slice
  integer, intent(in), dimension(NGLLX,NGLLZ,nspec) :: ibool

! global corner numbers that need to be created
  integer, intent(out), dimension(nglob_GLL_full) :: global_corner_number
  integer, intent(out) :: mn(nspec*NGNOD_QUADRANGLE),mp(nspec+1)
  integer, intent(out) :: nglob_four_corners_only

  integer ninter,nsum,ispec,node,k,inumcorner,ix,iy

  ninter = 1
  nsum = 1
  mp(1) = 1

!---- define topology of the elements in the mesh
!---- we need to define adjacent numbers from the sub-mesh consisting of the corners only
  nglob_four_corners_only = 0
  global_corner_number(:) = -1

  do ispec=1,nspec

    inumcorner = 0
      do iy = 1,NGLLZ,NGLLZ-1
        do ix = 1,NGLLX,NGLLX-1

          inumcorner = inumcorner + 1
          if(inumcorner > NGNOD_QUADRANGLE) stop 'corner number too large'

! check if this point was already assigned a number previously, otherwise create one and store it
          if(global_corner_number(ibool(ix,iy,ispec)) == -1) then
            nglob_four_corners_only = nglob_four_corners_only + 1
            global_corner_number(ibool(ix,iy,ispec)) = nglob_four_corners_only
          endif

          node = global_corner_number(ibool(ix,iy,ispec))
            do k=nsum,ninter-1
              if(node == mn(k)) goto 200
            enddo

            mn(ninter) = node
            ninter = ninter + 1
  200 continue

      enddo
    enddo

      nsum = ninter
      mp(ispec + 1) = nsum

  enddo

  end subroutine form_elt_connectivity_foelco

!
!----------------------------------------------------
!

  subroutine form_node_connectivity_fonoco(mn,mp,ne,np,nglob_four_corners_only, &
                                nspec,count_only,total_size_ne)

!-----------------------------------------------------------------------
!
!   Forms the NE and NP arrays
!
!     Input :
!     -------
!           MN, MP, nspec
!           nglob_four_corners_only    Number of nodes in the domain
!
!     Output :
!     --------
!           NE, NP   This is the node-connected element array pair.
!                    Integer array NE contains a list of the
!                    elements connected to each node, stored in stacked fashion.
!
!                    Array NP is the pointer array for the
!                    location of a node's element list in the NE array.
!                    Its length is equal to the number of points plus one.
!
!-----------------------------------------------------------------------

  implicit none

  include "constants.h"

! only count the total size of the array that will be created, or actually create it
  logical count_only
  integer total_size_ne

  integer nglob_four_corners_only,nspec

  integer, intent(in) ::  mn(nspec*NGNOD_QUADRANGLE),mp(nspec+1)

  integer, intent(out) ::  ne(total_size_ne),np(nglob_four_corners_only+1)

  integer nsum,inode,ispec,j

  nsum = 1
  np(1) = 1

  do inode=1,nglob_four_corners_only
      do 200 ispec=1,nspec

            do j=mp(ispec),mp(ispec + 1) - 1
                  if (mn(j) == inode) then
                        if(count_only) then
                          total_size_ne = nsum
                        else
                          ne(nsum) = ispec
                        endif
                        nsum = nsum + 1
                        goto 200
                  endif
            enddo
  200 continue

      np(inode + 1) = nsum

  enddo

  end subroutine form_node_connectivity_fonoco

!
!----------------------------------------------------
!

  subroutine create_adjacency_table_adjncy(mn,mp,ne,np,adj,xadj,maskel,nspec, &
              nglob_four_corners_only,count_only,total_size_ne,total_size_adj)

!-----------------------------------------------------------------------
!
!   Establishes the element adjacency information of the mesh
!   Two elements are considered adjacent if they share a face.
!
!     Input :
!     -------
!           MN, MP, NE, NP, nspec
!           MASKEL    logical mask (length = nspec)
!
!     Output :
!     --------
!           ADJ, XADJ This is the element adjacency array pair. Array
!                     ADJ contains the list of the elements adjacent to
!                     element i. They are stored in a stacked fashion.
!                     Pointer array XADJ stores the location of each element list.
!
!-----------------------------------------------------------------------

  implicit none

  include "constants.h"

! only count the total size of the array that will be created, or actually create it
  logical count_only
  integer total_size_ne,total_size_adj

  integer nglob_four_corners_only

  integer, intent(in) :: mn(nspec*NGNOD_QUADRANGLE),mp(nspec+1),ne(total_size_ne),np(nglob_four_corners_only+1)

  integer, intent(out) :: adj(total_size_adj),xadj(nspec+1)

  logical maskel(nspec)
  integer countel(nspec)

  integer nspec,iad,ispec,istart,istop,ino,node,jstart,jstop,nelem,jel

  xadj(1) = 1
  iad = 1

  do ispec=1,nspec

! reset mask
  maskel(:) = .false.

! mask current element
  maskel(ispec) = .true.
  if (FACE) countel(:) = 0

  istart = mp(ispec)
  istop = mp(ispec+1) - 1
    do ino=istart,istop
      node = mn(ino)
      jstart = np(node)
      jstop = np(node + 1) - 1
        do 120 jel=jstart,jstop
            nelem = ne(jel)
            if(maskel(nelem)) goto 120
            if (FACE) then
!! DK DK this below implemented by David Michea in 3D, but not true anymore in 2D: should be
!! DK DK two corners instead of three. But does not matter because FACE is always .false.
!! DK DK and therefore this part of the routine is currently never used.
!! DK DK Let me add a stop statement just in case.
              stop 'FACE = .true. not implemented, check the above comment in the source code'
!! DK DK End of the stop statement added.
              ! if 2 elements share at least 3 corners, therefore they share a face
              countel(nelem) = countel(nelem) + 1
              if (countel(nelem)>=3) then
                if(count_only) then
                  total_size_adj = iad
                else
                  adj(iad) = nelem
                endif
                maskel(nelem) = .true.
                iad = iad + 1
              endif
            else
              if(count_only) then
                total_size_adj = iad
              else
                adj(iad) = nelem
              endif
              maskel(nelem) = .true.
              iad = iad + 1
            endif
  120   continue
    enddo

    xadj(ispec+1) = iad

  enddo

  end subroutine create_adjacency_table_adjncy

!
!----------------------------------------------------
!

  subroutine cuthill_mckee(adj,xadj,mask,invperm_all,nspec,total_size_adj,limit)

  implicit none
  include "constants.h"

  integer, intent(in) :: nspec,total_size_adj, limit
  integer, intent(in) :: adj(total_size_adj),xadj(nspec+1)

  integer, intent(out), dimension(nspec) :: mask,invperm_all
  integer, dimension(nspec) :: invperm_sub
  integer ispec,gsize,counter,nspec_sub,root,total_ordered_elts, next_root

! fill the mask with ones
  mask(:) = 1
  invperm_all(:) = 0
  counter = 0
  nspec_sub = limit
  root = 1
  total_ordered_elts = 0

  do while(total_ordered_elts < nspec)
    ! creation of a sublist of sorted elements which fit in the cache (the criterion of size is limit)
    ! limit = nb of element that can fit in the L2 cache
    call Cut_McK( root, nspec, total_size_adj, xadj, adj, mask, gsize, invperm_sub, limit, nspec_sub, next_root)
      ! add the sublist in the main permutation list
      invperm_all(total_ordered_elts+1:total_ordered_elts+nspec_sub) = invperm_sub(1:nspec_sub)
      total_ordered_elts = total_ordered_elts + nspec_sub
    ! seek for a new root to build the new sublist
    if (next_root > 0) then
      root = next_root
    else
      if (total_ordered_elts /= nspec) &
        call find_next_root(next_root,xadj,adj,total_size_adj,mask,invperm_all,total_ordered_elts,nspec)
      root = next_root
    endif
  enddo

  if (INVERSE) then
    do ispec=1,nspec
      mask(invperm_all(ispec)) = ispec
    enddo
  else
    mask(:) = invperm_all(:)
  endif

  end subroutine cuthill_mckee


!*******************************************************************************
! Objective: Cuthill-McKee ordering
!    The algorithm is:
!
!    X(1) = ROOT.
!    for ( I = 1 to N-1)
!      Find all unlabeled neighbors of X(I),
!      assign them the next available labels, in order of increasing degree.
!
!  Parameters:
!    root       the starting point for the cm ordering.
!    nbnodes    the number of nodes.
!    nnz        the number of adjacency entries.
!
!    xadj/adj   the graph
!    mask       only those nodes with nonzero mask are considered
!
!    gsize      the number of the connected component
!    invp       Inverse permutation (from new order to old order)
!*******************************************************************************

subroutine find_next_root(next_root,xadj,adj,total_size_adj,mask,invperm_all,total_ordered_elts,nspec)

  implicit none

  include "constants.h"

! input
  integer, intent(in) :: total_size_adj,total_ordered_elts,nspec
  integer, intent(in) :: adj(total_size_adj),xadj(nspec+1)
  integer, intent(in), dimension(nspec) :: mask,invperm_all
! output
  integer, intent(out) :: next_root
! variables
  integer :: cur_node,neighbor_node,i,j

  do i=total_ordered_elts, 1, -1
    cur_node = invperm_all(i)
    do j= xadj(cur_node), xadj(cur_node+1)-1
      neighbor_node = adj(j)
      if (mask(neighbor_node)/=0) then
        next_root=neighbor_node
        return
      endif
    enddo
  enddo

end subroutine find_next_root

!*******************************************************************************
! Objective: Cuthill-McKee ordering
!    The algorithm is:
!
!    X(1) = ROOT.
!    for ( I = 1 to N-1)
!      Find all unlabeled neighbors of X(I),
!      assign them the next available labels, in order of increasing degree.
!
!  Parameters:
!    root       the starting point for the cm ordering.
!    nbnodes    the number of nodes.
!    nnz        the number of adjacency entries.
!
!    xadj/adj   the graph
!    mask       only those nodes with nonzero mask are considered
!
!    gsize      the number of the connected component
!    invp       Inverse permutation (from new order to old order)
!*******************************************************************************

subroutine Cut_McK( root, nbnodes, nnz, xadj, adj, mask, gsize, invp, limit, nspec_sub, next_root)

  implicit none

  include "constants.h"

!--------------------------------------------------------------- Input Variables
  integer root, nnz, nbnodes, limit, nspec_sub, next_root

  integer xadj(nbnodes+1), adj(nnz), mask(nbnodes)

!-------------------------------------------------------------- Output Variables
  integer gsize
  integer invp(nbnodes)

!--------------------------------------------------------------- Local Variables
  integer i, j, k, l, lbegin, lnbr, linvp, lvlend, nbr, node, fnbr
  integer deg(nbnodes)

! Find the degrees of the nodes in the subgraph specified by mask and root
! Here invp is used to store a levelization of the subgraph
  invp(:)=0
  deg(:)=0
  call degree ( root, nbnodes, nnz, xadj, adj, mask, gsize, deg, invp)

  mask(root) = 0

  IF (gsize > 1) THEN
    !If there is at least 2 nodes in the subgraph
    lvlend = 0
    lnbr   = 1

    DO while (lvlend < lnbr)
      !lbegin/lvlend point to the begin/end of the present level
      lbegin = lvlend + 1
      lvlend = lnbr

      do i= lbegin, lvlend
        node = invp(i)

        !Find the unnumbered neighbours of node.
        !fnbr/lnbr point to the first/last neighbors of node
        fnbr = lnbr + 1
        do j= xadj(node), xadj(node+1)-1
          nbr = adj(j)

          if (mask(nbr) /= 0) then
            lnbr       = lnbr + 1
            mask(nbr)  = 0
            invp(lnbr) = nbr
          endif
        enddo

        !If no neighbors, go to next node in this level.
        IF (lnbr > fnbr) THEN
          !Sort the neighbors of NODE in increasing order by degree.
          !Linear insertion is used.
          k = fnbr
          do while (k < lnbr)
            l   = k
            k   = k + 1
            nbr = invp(k)

            DO WHILE (fnbr < l)
              linvp = invp(l)

              if (deg(linvp) <= deg(nbr)) then
                exit
              endif

              invp(l+1) = linvp
              l         = l-1
            ENDDO

            invp(l+1) = nbr
          enddo
        ENDIF
      enddo
    ENDDO

  ENDIF

  if (gsize > limit) then
    do i = limit + 1 , nbnodes
      node=invp(i)
      if (node /=0) mask(node) = 1
    enddo
    next_root = invp(limit +1)
    nspec_sub = limit
  else
    next_root = -1
    nspec_sub = gsize
  endif

END subroutine Cut_McK


!*******************************************************************************
! Objective: computes the degrees of the nodes in the connected graph
!
! Parameters:
!    root       the root node
!    nbnodes    the number of nodes in the graph
!    nnz        the graph size
!    xadj/adj   the whole graph
!    mask       Only nodes with mask == 0 are considered
!
!    gsize      the number of nodes in the connected graph
!    deg        degree for all the nodes in the connected graph
!    level      levelization of the connected graph
!
!*******************************************************************************

subroutine degree( root, nbnodes, nnz, xadj, adj, mask, gsize, deg, level )

  implicit none

!--------------------------------------------------------------- Input Variables
  integer root, nbnodes, nnz
  integer xadj(nbnodes+1), adj(nnz), mask(nbnodes)

!-------------------------------------------------------------- Output Variables
  integer gsize
  integer deg(nbnodes), level(nbnodes)

!--------------------------------------------------------------- Local Variables
  integer i, j, ideg, lbegin, lvlend, lvsize, nxt, nbr, node

! added a test to detect disconnected subsets in the mesh
! (in which case Cuthill-McKee fails and should be turned off)
  if(root > nbnodes+1) stop 'error: root > nbnodes+1 in Cuthill-McKee'
  if(root < 1) then
    print *,'error: root < 1 in Cuthill-McKee; you probably have a mesh composed of'
    print *,'two disconnected subsets of elements, in which case Cuthill-McKee fails and should be turned off.'
    print *,'please set PERFORM_CUTHILL_MCKEE = .false. in constants.h and recompile.'
    print *,'please also doublecheck that you indeed want to run two separate meshes simultaneously,'
    print *,'which is extremely unusual (but formally not incorrect).'
    stop 'fatal error in Cuthill-McKee'
  endif

! The sign of xadj(I) is used to indicate if node i has been considered
  xadj(root) = -xadj(root)
  level(1)   = root
  nxt        = 1
  lvlend     = 0
  lvsize     = 1

  DO WHILE (lvsize > 0)
    ! lbegin/lvlend points the begin/end of the present level
    lbegin = lvlend + 1
    lvlend = nxt

    ! Find the degrees of nodes in the present level and generate the next level
    DO i= lbegin, lvlend
      node  = level(i)
      ideg  = 0
      do j= ABS( xadj(node) ), ABS( xadj(node+1) )-1
        nbr = adj(j)

        if (mask(nbr) /= 0) then
          ideg = ideg + 1

          if (xadj(nbr) >= 0) then
            xadj(nbr)  = -xadj(nbr)
            nxt        = nxt  + 1
            level(nxt) = nbr
          endif
        endif
      enddo

      deg(node) = ideg
    ENDDO

    !Compute the level size of the next level
    lvsize = nxt - lvlend
  ENDDO

  !Reset xadj to its correct sign
  do i = 1, nxt
    node       = level(i)
    xadj(node) = -xadj(node)
  enddo

  gsize = nxt

END subroutine degree

!
!-----------------------------------------------------------------------
!

  subroutine permute_elements_real(array_to_permute,temp_array,perm,nspec)

  implicit none

  include "constants.h"

  integer, intent(in) :: nspec
  integer, intent(in), dimension(nspec) :: perm

  real(kind=CUSTOM_REAL), intent(inout), dimension(NGLLX,NGLLZ,nspec) :: array_to_permute,temp_array

  integer old_ispec,new_ispec

! copy the original array
  temp_array(:,:,:) = array_to_permute(:,:,:)

  do old_ispec = 1,nspec
    new_ispec = perm(old_ispec)
    array_to_permute(:,:,new_ispec) = temp_array(:,:,old_ispec)
  enddo

  end subroutine permute_elements_real

!
!-----------------------------------------------------------------------
!

! implement permutation of elements for arrays of integer type
  subroutine permute_elements_integer(array_to_permute,temp_array,perm,nspec)

  implicit none

  include "constants.h"

  integer, intent(in) :: nspec
  integer, intent(in), dimension(nspec) :: perm

  integer, intent(inout), dimension(NGLLX,NGLLZ,nspec) :: array_to_permute,temp_array

  integer old_ispec,new_ispec

! copy the original array
  temp_array(:,:,:) = array_to_permute(:,:,:)

  do old_ispec = 1,nspec
    new_ispec = perm(old_ispec)
    array_to_permute(:,:,new_ispec) = temp_array(:,:,old_ispec)
  enddo

  end subroutine permute_elements_integer

!
!-----------------------------------------------------------------------
!

! implement permutation of elements for arrays of double precision type
  subroutine permute_elements_dble(array_to_permute,temp_array,perm,nspec)

  implicit none

  include "constants.h"

  integer, intent(in) :: nspec
  integer, intent(in), dimension(nspec) :: perm

  double precision, intent(inout), dimension(NGLLX,NGLLZ,nspec) :: array_to_permute,temp_array

  integer old_ispec,new_ispec

! copy the original array
  temp_array(:,:,:) = array_to_permute(:,:,:)

  do old_ispec = 1,nspec
    new_ispec = perm(old_ispec)
    array_to_permute(:,:,new_ispec) = temp_array(:,:,old_ispec)
  enddo

  end subroutine permute_elements_dble

