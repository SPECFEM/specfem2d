
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, INRIA and University of Pau, France,
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

!
! This module contains subroutines related to unstructured meshes and partitioning of the
! corresponding graphs.
!

module part_unstruct

  implicit none

  integer :: nelmnts
  integer, dimension(:), pointer  :: elmnts
  integer, dimension(:), allocatable  :: elmnts_bis
  integer, dimension(:), allocatable  :: vwgt
  integer, dimension(:), allocatable  :: glob2loc_elmnts
  integer, dimension(:), allocatable  :: part

  integer :: nb_edges
  integer, dimension(:), allocatable  :: adjwgt

  integer, dimension(:), allocatable  :: xadj_g
  integer, dimension(:), allocatable  :: adjncy_g

  integer :: nnodes
  double precision, dimension(:,:), allocatable  :: nodes_coords
  integer, dimension(:), allocatable  :: nnodes_elmnts
  integer, dimension(:), allocatable  :: nodes_elmnts
  integer, dimension(:), allocatable  :: glob2loc_nodes_nparts
  integer, dimension(:), allocatable  :: glob2loc_nodes_parts
  integer, dimension(:), allocatable  :: glob2loc_nodes

  ! interface data
  integer :: ninterfaces
  integer, dimension(:), allocatable  :: tab_size_interfaces, tab_interfaces

  integer :: nelem_acoustic_surface
  integer, dimension(:,:), pointer  :: acoustic_surface
  integer :: nelem_acoustic_surface_loc

  integer :: nelemabs
  integer, dimension(:,:), allocatable  :: abs_surface
  logical, dimension(:,:), allocatable  :: abs_surface_char
  integer, dimension(:), allocatable  :: abs_surface_merge,abs_surface_type
  integer :: nelemabs_loc

  integer :: nelemabs_merge
  integer, dimension(:), allocatable  :: ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
       ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2

  ! for acoustic/elastic coupled elements
  integer :: nedges_coupled
  integer, dimension(:,:), pointer  :: edges_coupled

  ! for acoustic/poroelastic coupled elements
  integer :: nedges_acporo_coupled
  integer, dimension(:,:), pointer  :: edges_acporo_coupled

  ! for poroelastic/elastic coupled elements
  integer :: nedges_elporo_coupled
  integer, dimension(:,:), pointer  :: edges_elporo_coupled

contains

  !-----------------------------------------------
  ! Read the mesh and storing it in array 'elmnts' (which is allocated here).
  ! 'remove_min_to_start_at_zero' is used to have the numbering of the nodes starting at '0'.
  ! 'nelmnts' is the number of elements, 'nnodes' is the number of nodes in the mesh.
  !-----------------------------------------------
  subroutine read_external_mesh_file(filename, remove_min_to_start_at_zero, ngnod)

  implicit none
  !include "constants.h"

  character(len=256), intent(in)  :: filename
  integer, intent(out)  :: remove_min_to_start_at_zero
  integer, intent(in)  :: ngnod

  integer  :: i,ier

  open(unit=990, file=trim(filename), form='formatted' , status='old', action='read',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening file: ',trim(filename)
    stop 'error read external mesh file'
  endif

  read(990,*) nelmnts

  allocate(elmnts(0:ngnod*nelmnts-1))

  do i = 0, nelmnts-1
    if(ngnod == 4) then
      read(990,*) elmnts(i*ngnod), elmnts(i*ngnod+1), elmnts(i*ngnod+2), elmnts(i*ngnod+3)
    else if(ngnod == 9) then
      read(990,*) elmnts(i*ngnod), elmnts(i*ngnod+1), elmnts(i*ngnod+2), elmnts(i*ngnod+3), &
                  elmnts(i*ngnod+4), elmnts(i*ngnod+5), elmnts(i*ngnod+6), elmnts(i*ngnod+7), elmnts(i*ngnod+8)
    else
      stop 'error, ngnod should be either 4 or 9 for external meshes'
    endif
  enddo

  close(990)

  remove_min_to_start_at_zero = minval(elmnts)
  elmnts(:) = elmnts(:) - remove_min_to_start_at_zero
  nnodes = maxval(elmnts) + 1

  end subroutine read_external_mesh_file

  !-----------------------------------------------
  ! Read the nodes coordinates and storing it in array 'nodes_coords'
  !-----------------------------------------------
  subroutine read_nodes_coords(filename)

  implicit none

  character(len=256), intent(in)  :: filename

  integer  :: i,ier

  open(unit=991, file=trim(filename), form='formatted' , status='old', action='read', iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening file: ',trim(filename)
    stop 'error read external nodes coords file'
  endif

  read(991,*) nnodes
  allocate(nodes_coords(2,nnodes))
  do i = 1, nnodes
     read(991,*) nodes_coords(1,i), nodes_coords(2,i)
  enddo
  close(991)

  end subroutine read_nodes_coords


  !-----------------------------------------------
  ! Read the material for each element and storing it in array 'num_materials'
  !-----------------------------------------------
  subroutine read_mat(filename, num_material)

  implicit none

  character(len=256), intent(in)  :: filename
  integer, dimension(1:nelmnts), intent(out)  :: num_material

  integer  :: i,ier

  open(unit=992, file=trim(filename), form='formatted' , status='old', action='read',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening file: ',trim(filename)
    stop 'error read external mat file'
  endif

  do i = 1, nelmnts
     read(992,*) num_material(i)
  enddo
  close(992)

  end subroutine read_mat

  !-----------------------------------------------
  ! Read the position of pml element storing it in array 'region_pml_external_mesh'
  !-----------------------------------------------
  subroutine read_pml_element(filename, region_pml_external_mesh, nspec_cpml)

  implicit none

  character(len=256), intent(in)  :: filename
  integer, dimension(1:nelmnts), intent(out)  :: region_pml_external_mesh
  integer, intent(out)  :: nspec_cpml
!  integer, dimension(:,:), allocatable  :: local_pml  

  integer  :: i,j,k,ier

  open(unit=992, file=trim(filename), form='formatted' , status='old', action='read',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening file: ',trim(filename)
    stop 'error read external CPML_element_file'
  endif

  read(992,*) nspec_cpml

  do i = 1, nspec_cpml
     read(992,*) j, k
     region_pml_external_mesh(j) = k
  enddo

  close(992)

  end subroutine read_pml_element



  !-----------------------------------------------
  ! Read free surface.
  ! Edges from elastic elements are discarded.
  ! 'acoustic_surface' contains 1/ element number, 2/ number of nodes that form the free surface,
  ! 3/ first node on the free surface, 4/ second node on the free surface, if relevant (if 2/ is equal to 2)
  !-----------------------------------------------
  subroutine read_acoustic_surface(filename, num_material, &
                ANISOTROPIC_MATERIAL, nb_materials, icodemat, phi, remove_min_to_start_at_zero)

  implicit none

  !include "constants.h"

  character(len=256), intent(in)  :: filename
  integer, dimension(0:nelmnts-1)  :: num_material
  integer, intent(in)  :: ANISOTROPIC_MATERIAL
  integer, intent(in)  :: nb_materials
  integer, dimension(1:nb_materials), intent(in)  :: icodemat
  double precision, dimension(1:nb_materials), intent(in)  :: phi
  integer, intent(in)  :: remove_min_to_start_at_zero


  integer, dimension(:,:), allocatable  :: acoustic_surface_tmp
  integer  :: nelmnts_surface
  integer  :: i,ier
  integer  :: imaterial_number


  open(unit=993, file=trim(filename), form='formatted' , status='old', action='read', iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening file: ',trim(filename)
    stop 'error read acoustic surface file'
  endif

  read(993,*) nelmnts_surface

  allocate(acoustic_surface_tmp(4,nelmnts_surface))

  do i = 1, nelmnts_surface
     read(993,*) acoustic_surface_tmp(1,i), acoustic_surface_tmp(2,i), acoustic_surface_tmp(3,i), acoustic_surface_tmp(4,i)

  enddo

  close(993)
  acoustic_surface_tmp(1,:) = acoustic_surface_tmp(1,:) - remove_min_to_start_at_zero
  acoustic_surface_tmp(3,:) = acoustic_surface_tmp(3,:) - remove_min_to_start_at_zero
  acoustic_surface_tmp(4,:) = acoustic_surface_tmp(4,:) - remove_min_to_start_at_zero

  nelem_acoustic_surface = 0
  do i = 1, nelmnts_surface
     imaterial_number = num_material(acoustic_surface_tmp(1,i))
     if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. phi(imaterial_number) >= 1.d0 ) then
        nelem_acoustic_surface = nelem_acoustic_surface + 1

     endif
  enddo

  allocate(acoustic_surface(4,nelem_acoustic_surface))

  nelem_acoustic_surface = 0
  do i = 1, nelmnts_surface
     imaterial_number = num_material(acoustic_surface_tmp(1,i))
     if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. phi(imaterial_number) >= 1.d0 ) then
        nelem_acoustic_surface = nelem_acoustic_surface + 1
        acoustic_surface(:,nelem_acoustic_surface) = acoustic_surface_tmp(:,i)
     endif
  enddo

  end subroutine read_acoustic_surface


  !-----------------------------------------------
  ! Read absorbing surface.
  ! 'abs_surface' contains 1/ element number, 2/ number of nodes that form the absorbing edge
  ! (which currently must always be equal to two, see comment below),
  ! 3/ first node on the abs surface, 4/ second node on the abs surface
  !-----------------------------------------------
  subroutine read_abs_surface(filename, remove_min_to_start_at_zero)

  implicit none
  !include "constants.h"

  character(len=256), intent(in)  :: filename
  integer, intent(in)  :: remove_min_to_start_at_zero

  integer  :: i,ier

  open(unit=994, file=trim(filename), form='formatted' , status='old', action='read', iostat=ier)
  if( ier /= 0 ) then
    print *,'error opening file: ',trim(filename)
    stop 'error read absorbing surface file'
  endif

  read(994,*) nelemabs

  allocate(abs_surface(5,nelemabs))

  do i = 1, nelemabs
    read(994,*) abs_surface(1,i), abs_surface(2,i), abs_surface(3,i), abs_surface(4,i), abs_surface(5,i)

    if (abs_surface(2,i) /= 2) then
      print *,'Only two nodes per absorbing element can be listed.'
      print *,'If one of your elements has more than one edge along a given absorbing contour'
      print *,'(e.g., if that contour has a corner) then list it twice,'
      print *,'putting the first edge on the first line and the second edge on the second line.'
      print *,'If one of your elements has a single point along the absording contour rather than a full edge, do NOT list it'
      print *,'(it would have no weight in the contour integral anyway because it would consist of a single point).'
      print *,'If you use 9-node elements, list only the first and last points of the edge and not the intermediate point'
      print *,'located around the middle of the edge; the right 9-node curvature will be restored automatically by the code.'

      stop 'only two nodes per element should be listed for absorbing edges'
    endif

    if (abs_surface(5,i) < 1 .or. abs_surface(5,i) > 4) then
      stop 'absorbing element type must be between 1 (IBOTTOM) and 4 (ILEFT)'
    endif

  enddo

  close(994)

  abs_surface(1,:) = abs_surface(1,:) - remove_min_to_start_at_zero
  abs_surface(3,:) = abs_surface(3,:) - remove_min_to_start_at_zero
  abs_surface(4,:) = abs_surface(4,:) - remove_min_to_start_at_zero

  end subroutine read_abs_surface

  !-----------------------------------------------
  ! rotate_mesh_for_plane_wave.
  ! rotate mesh elements to make sure topological absorbing surface
  ! is aligned with physical absorbing surface, since this is the surface
  ! that we use to impose the plane wave and Bielak boundary conditions.
  !-----------------------------------------------

 subroutine rotate_mesh_for_plane_wave(ngnod)

 implicit none

 integer, intent(in)  :: ngnod
 integer :: i,j,ispec,i1,i2,inode,iswap
 logical :: found_this_point
 integer, dimension(:,:), allocatable :: ibool,ibool_rotated
!! DK DK be careful here, "ibool" applies to the mesh corners (4 or 9 points) only,
!! DK DK not to the GLL points because there are no GLL points in the Gmsh mesh files
  integer :: index_rotation1,index_rotation2,index_rotation3,index_rotation4,&
             index_rotation5,index_rotation6,index_rotation7,index_rotation8,index_edge

 allocate(ibool(ngnod,nelmnts))
 allocate(ibool_rotated(ngnod,nelmnts))

  do ispec = 1, nelmnts
    if(ngnod == 4) then
      ibool(1,ispec) = elmnts((ispec-1)*ngnod)
      ibool(2,ispec) = elmnts((ispec-1)*ngnod+1)
      ibool(3,ispec) = elmnts((ispec-1)*ngnod+2)
      ibool(4,ispec) = elmnts((ispec-1)*ngnod+3)
    else if(ngnod == 9) then
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
      stop 'error, ngnod should be either 4 or 9 for external meshes'
    endif
  enddo


do j = 1, 4
  if(j==1)then
    index_edge=3
    ibool_rotated(:,:) = ibool(:,:)
  else if(j==2)then
    index_edge=1
    ibool(:,:) = ibool_rotated(:,:)
  else if(j==3)then
    index_edge=4
    ibool(:,:) = ibool_rotated(:,:)
  else if(j==4)then
    index_edge=2
    ibool(:,:) = ibool_rotated(:,:)
  else
    stop 'j should be >=1 and <=4'
  endif

   if(index_edge==1)then
! bottome edge
      index_rotation1 = 1
      index_rotation2 = 2
      index_rotation3 = 2
      index_rotation4 = 3
      index_rotation5 = 3
      index_rotation6 = 4
      index_rotation7 = 1
      index_rotation8 = 4
   else if(index_edge==2)then
! right edge
      index_rotation1 = 2
      index_rotation2 = 3
      index_rotation3 = 3
      index_rotation4 = 4
      index_rotation5 = 1
      index_rotation6 = 4
      index_rotation7 = 1
      index_rotation8 = 2
   else if(index_edge==3)then
! top edge
      index_rotation1 = 3
      index_rotation2 = 4
      index_rotation3 = 1
      index_rotation4 = 4
      index_rotation5 = 1
      index_rotation6 = 2
      index_rotation7 = 2
      index_rotation8 = 3
   else if(index_edge==4)then
! left edge
      index_rotation1 = 1
      index_rotation2 = 4
      index_rotation3 = 1
      index_rotation4 = 2
      index_rotation5 = 2
      index_rotation6 = 3
      index_rotation7 = 3
      index_rotation8 = 4
   else
     stop 'The edge on which abs_nodes is located should be defined'
   endif

  do i = 1,nelemabs
   if(index_edge==abs_surface(5,i))then
   ispec = abs_surface(1,i) + 1  !!!! be careful: ispec from abs_surface(1,i) start at zero
   found_this_point = .false.
   do inode = 1,ngnod
      if(ibool(inode,ispec) == abs_surface(3,i)) then
       i1 = inode
       found_this_point = .true.
      exit
    endif
   enddo

  if(.not. found_this_point) stop 'point not found'

  found_this_point = .false.
  do inode = 1,4
    if(ibool(inode,ispec) == abs_surface(4,i)) then
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
   if(i1 == index_rotation1 .and. i2 == index_rotation2) then
!   print *,'orientation of element ',i,' is already good'

   else if (i1 == index_rotation3 .and. i2 == index_rotation4) then ! for this one, remember that we have swapped, thus 4 1 is 1 4
!   print *,'element ',i,' must be rotated 3 times'
    ibool_rotated(4,ispec) = ibool(1,ispec)
    ibool_rotated(1,ispec) = ibool(2,ispec)
    ibool_rotated(2,ispec) = ibool(3,ispec)
    ibool_rotated(3,ispec) = ibool(4,ispec)
    if(ngnod == 9) then
      ibool_rotated(8,ispec) = ibool(5,ispec)
      ibool_rotated(5,ispec) = ibool(6,ispec)
      ibool_rotated(6,ispec) = ibool(7,ispec)
      ibool_rotated(7,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element
    endif

   else if (i1 == index_rotation5 .and. i2 == index_rotation6) then
!   print *,'element ',i,ispec,' must be rotated 2 times top'
    ibool_rotated(3,ispec) = ibool(1,ispec)
    ibool_rotated(4,ispec) = ibool(2,ispec)
    ibool_rotated(1,ispec) = ibool(3,ispec)
    ibool_rotated(2,ispec) = ibool(4,ispec)
    if(ngnod == 9) then
      ibool_rotated(7,ispec) = ibool(5,ispec)
      ibool_rotated(8,ispec) = ibool(6,ispec)
      ibool_rotated(5,ispec) = ibool(7,ispec)
      ibool_rotated(6,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element
    endif

   else if (i1 == index_rotation7 .and. i2 == index_rotation8) then
!   print *,'element ',i,' must be rotated 1 time'
    ibool_rotated(2,ispec) = ibool(1,ispec)
    ibool_rotated(3,ispec) = ibool(2,ispec)
    ibool_rotated(4,ispec) = ibool(3,ispec)
    ibool_rotated(1,ispec) = ibool(4,ispec)
    if(ngnod == 9) then
      ibool_rotated(6,ispec) = ibool(5,ispec)
      ibool_rotated(7,ispec) = ibool(6,ispec)
      ibool_rotated(8,ispec) = ibool(7,ispec)
      ibool_rotated(5,ispec) = ibool(8,ispec)
! 9th point is at the element center and thus never changes when we rotate an element
    endif

   else
     stop 'problem in an element'
   endif

  endif

  enddo

 enddo

  do ispec = 1, nelmnts
    if(ngnod == 4) then
      elmnts((ispec-1)*ngnod)   = ibool_rotated(1,ispec)
      elmnts((ispec-1)*ngnod+1) = ibool_rotated(2,ispec)
      elmnts((ispec-1)*ngnod+2) = ibool_rotated(3,ispec)
      elmnts((ispec-1)*ngnod+3) = ibool_rotated(4,ispec)
    else if(ngnod == 9) then
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
      stop 'error, ngnod should be either 4 or 9 for external meshes'
    endif
  enddo

end  subroutine rotate_mesh_for_plane_wave

  !-----------------------------------------------
  ! Creating dual graph (adjacency is defined by 'ncommonnodes' between two elements).
  !-----------------------------------------------
  subroutine mesh2dual_ncommonnodes(elmnts_l,ncommonnodes,xadj,adjncy)

  implicit none
  include "constants.h"

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, intent(in)  :: ncommonnodes
  integer, dimension(0:nelmnts),intent(out)  :: xadj
  integer, dimension(0:MAX_NEIGHBORS*nelmnts-1),intent(out) :: adjncy

  ! local parameters
  integer  :: i, j, k, l, m, num_edges
  logical  ::  is_neighbour
  integer  :: num_node, n
  integer  :: elem_base, elem_target
  integer  :: connectivity

  ! allocates memory for arrays
  if( .not. allocated(nnodes_elmnts) ) allocate(nnodes_elmnts(0:nnodes-1))
  if( .not. allocated(nodes_elmnts) ) allocate(nodes_elmnts(0:nsize*nnodes-1))

  ! initializes
  xadj(:) = 0
  adjncy(:) = 0
  nnodes_elmnts(:) = 0
  nodes_elmnts(:) = 0
  num_edges = 0

  ! list of elements per node
  do i = 0, NCORNERS*nelmnts-1
    nodes_elmnts(elmnts_l(i)*nsize + nnodes_elmnts(elmnts_l(i))) = i/NCORNERS
    nnodes_elmnts(elmnts_l(i)) = nnodes_elmnts(elmnts_l(i)) + 1
  enddo

  ! checking which elements are neighbours ('ncommonnodes' criteria)
  do j = 0, nnodes-1
    do k = 0, nnodes_elmnts(j)-1
      do l = k+1, nnodes_elmnts(j)-1

        connectivity = 0
        elem_base = nodes_elmnts(k+j*nsize)
        elem_target = nodes_elmnts(l+j*nsize)
        do n = 1, NCORNERS
          num_node = elmnts_l(NCORNERS*elem_base+n-1)
          do m = 0, nnodes_elmnts(num_node)-1
            if ( nodes_elmnts(m+num_node*nsize) == elem_target ) then
              connectivity = connectivity + 1
            endif
          enddo
        enddo

        ! sets adjacency (adjncy) and number of neighbors (xadj)
        ! according to ncommonnodes criteria
        if ( connectivity >=  ncommonnodes) then

          is_neighbour = .false.

          do m = 0, xadj(nodes_elmnts(k+j*nsize))
            if ( .not.is_neighbour ) then
              if ( adjncy(nodes_elmnts(k+j*nsize)*MAX_NEIGHBORS+m) == nodes_elmnts(l+j*nsize) ) then
                is_neighbour = .true.
              endif
            endif
          enddo
          if ( .not.is_neighbour ) then
            adjncy(nodes_elmnts(k+j*nsize)*MAX_NEIGHBORS &
                   + xadj(nodes_elmnts(k+j*nsize))) = nodes_elmnts(l+j*nsize)

            xadj(nodes_elmnts(k+j*nsize)) = xadj(nodes_elmnts(k+j*nsize)) + 1
            if (xadj(nodes_elmnts(k+j*nsize)) > MAX_NEIGHBORS) &
              stop 'ERROR : too much neighbours per element, modify the mesh.'

            adjncy(nodes_elmnts(l+j*nsize)*MAX_NEIGHBORS &
                   + xadj(nodes_elmnts(l+j*nsize))) = nodes_elmnts(k+j*nsize)

            xadj(nodes_elmnts(l+j*nsize)) = xadj(nodes_elmnts(l+j*nsize)) + 1
            if (xadj(nodes_elmnts(l+j*nsize))>MAX_NEIGHBORS) &
              stop 'ERROR : too much neighbours per element, modify the mesh.'

          endif
        endif
      enddo
    enddo
  enddo

  ! making adjacency arrays compact (to be used for partitioning)
  do i = 0, nelmnts-1
    k = xadj(i)
    xadj(i) = num_edges
    do j = 0, k-1
      adjncy(num_edges) = adjncy(i*MAX_NEIGHBORS+j)
      num_edges = num_edges + 1
    enddo
  enddo

  xadj(nelmnts) = num_edges

  end subroutine mesh2dual_ncommonnodes


  !-----------------------------------------------
  ! Read the weight for each vertices and edges of the graph (not curretly used)
  !-----------------------------------------------
  subroutine read_weights()

  implicit none

  allocate(vwgt(0:nelmnts-1))
  allocate(adjwgt(0:nb_edges-1))

  vwgt(:) = 1
  adjwgt(:) = 1

  end subroutine read_weights


  !--------------------------------------------------
  ! construct local numbering for the elements in each partition
  !--------------------------------------------------
  subroutine Construct_glob2loc_elmnts(nparts)

  implicit none
  integer, intent(in)  :: nparts

  integer  :: num_glob, num_part
  integer, dimension(0:nparts-1)  :: num_loc


  allocate(glob2loc_elmnts(0:nelmnts-1))

  ! initializes number of local elements per partition
  do num_part = 0, nparts-1
    num_loc(num_part) = 0
  enddo

  ! local numbering
  do num_glob = 0, nelmnts-1
    num_part = part(num_glob)
    glob2loc_elmnts(num_glob) = num_loc(num_part)
    num_loc(num_part) = num_loc(num_part) + 1
  enddo

  end subroutine Construct_glob2loc_elmnts


  !--------------------------------------------------
  ! construct local numbering for the nodes in each partition
  !--------------------------------------------------
  subroutine Construct_glob2loc_nodes(nparts)

  implicit none
  include "constants.h"

  integer, intent(in)  :: nparts

  integer  :: num_node
  integer  :: el
  integer  ::  num_part
  integer  ::  size_glob2loc_nodes
  integer, dimension(0:nparts-1)  :: parts_node
  integer, dimension(0:nparts-1)  :: num_parts

  allocate(glob2loc_nodes_nparts(0:nnodes))

  size_glob2loc_nodes = 0

  parts_node(:) = 0


  do num_node = 0, nnodes-1
     glob2loc_nodes_nparts(num_node) = size_glob2loc_nodes
     do el = 0, nnodes_elmnts(num_node)-1
        parts_node(part(nodes_elmnts(el+nsize*num_node))) = 1
     enddo

     do num_part = 0, nparts-1
        if ( parts_node(num_part) == 1 ) then
           size_glob2loc_nodes = size_glob2loc_nodes + 1
           parts_node(num_part) = 0
        endif
     enddo

  enddo

  glob2loc_nodes_nparts(nnodes) = size_glob2loc_nodes

  allocate(glob2loc_nodes_parts(0:glob2loc_nodes_nparts(nnodes)-1))
  allocate(glob2loc_nodes(0:glob2loc_nodes_nparts(nnodes)-1))

  glob2loc_nodes(0) = 0

  parts_node(:) = 0
  num_parts(:) = 0
  size_glob2loc_nodes = 0


  do num_node = 0, nnodes-1
     do el = 0, nnodes_elmnts(num_node)-1
        parts_node(part(nodes_elmnts(el+nsize*num_node))) = 1
     enddo
     do num_part = 0, nparts-1

        if ( parts_node(num_part) == 1 ) then
           glob2loc_nodes_parts(size_glob2loc_nodes) = num_part
           glob2loc_nodes(size_glob2loc_nodes) = num_parts(num_part)
           size_glob2loc_nodes = size_glob2loc_nodes + 1
           num_parts(num_part) = num_parts(num_part) + 1
           parts_node(num_part) = 0
        endif

     enddo
  enddo

  end subroutine Construct_glob2loc_nodes


  !--------------------------------------------------
  ! Construct interfaces between each partitions.
  ! Two adjacent elements in distinct partitions make an entry in array tab_interfaces :
  ! 1/ first element, 2/ second element, 3/ number of common nodes, 4/ first node,
  ! 5/ second node, if relevant.
  ! No interface between acoustic, elastic, and poroelastic elements.
  !--------------------------------------------------
  subroutine Construct_interfaces(nparts, elmnts_l,  &
                                nb_materials, phi_material, num_material)

  implicit none
  include "constants.h"

  integer, intent(in)  :: nparts
  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, dimension(1:nelmnts), intent(in)  :: num_material
  integer, intent(in)  :: nb_materials
  double precision, dimension(1:nb_materials), intent(in)  :: phi_material

  integer  :: num_part, num_part_bis, el, el_adj, num_interface, num_edge, ncommon_nodes, &
       num_node, num_node_bis
  integer  :: i, j
  logical  :: is_acoustic_el, is_acoustic_el_adj, is_elastic_el, is_elastic_el_adj

  ninterfaces = 0
  do  i = 0, nparts-1
     do j = i+1, nparts-1
        ninterfaces = ninterfaces + 1
     enddo
  enddo

  allocate(tab_size_interfaces(0:ninterfaces))
  tab_size_interfaces(:) = 0

  num_interface = 0
  num_edge = 0

  do num_part = 0, nparts-1
     do num_part_bis = num_part+1, nparts-1
        do el = 0, nelmnts-1
           if ( part(el) == num_part ) then
              ! sets material flag
              if ( phi_material(num_material(el+1)) < TINYVAL) then
                ! elastic element
                is_acoustic_el = .false.
                is_elastic_el = .true.
              else if ( phi_material(num_material(el+1)) >= 1.d0) then
                ! acoustic element
                is_acoustic_el = .true.
                is_elastic_el = .false.
              else
                ! poroelastic element
                is_acoustic_el = .false.
                is_elastic_el = .false.
              endif

              ! looks at all neighbor elements
              do el_adj = xadj_g(el), xadj_g(el+1)-1
                ! sets neighbor material flag
                if ( phi_material(num_material(adjncy_g(el_adj)+1)) < TINYVAL) then
                  is_acoustic_el_adj = .false.
                  is_elastic_el_adj = .true.
                else if ( phi_material(num_material(adjncy_g(el_adj)+1)) >= 1.d0) then
                  is_acoustic_el_adj = .true.
                  is_elastic_el_adj = .false.
                else
                  is_acoustic_el_adj = .false.
                  is_elastic_el_adj = .false.
                endif
                ! adds element if neighbor element lies in next parition
                ! and belongs to same material
                if ( (part(adjncy_g(el_adj)) == num_part_bis) .and. &
                     (is_acoustic_el .eqv. is_acoustic_el_adj) .and. &
                     (is_elastic_el .eqv. is_elastic_el_adj) ) then
                    num_edge = num_edge + 1
                endif
              enddo
           endif
        enddo
        ! stores number of elements at interface
        tab_size_interfaces(num_interface+1) = tab_size_interfaces(num_interface) + num_edge
        num_edge = 0
        num_interface = num_interface + 1

     enddo
  enddo

  ! stores element indices for elements from above search at each interface
  num_interface = 0
  num_edge = 0

  allocate(tab_interfaces(0:(tab_size_interfaces(ninterfaces)*5-1)))
  tab_interfaces(:) = 0

  do num_part = 0, nparts-1
    do num_part_bis = num_part+1, nparts-1
      do el = 0, nelmnts-1
        if ( part(el) == num_part ) then
          if ( phi_material(num_material(el+1)) < TINYVAL) then
            is_acoustic_el = .false.
            is_elastic_el = .true.
          else if ( phi_material(num_material(el+1)) >= 1.d0) then
            is_acoustic_el = .true.
            is_elastic_el = .false.
          else
            is_acoustic_el = .false.
            is_elastic_el = .false.
          endif
          do el_adj = xadj_g(el), xadj_g(el+1)-1
            if ( phi_material(num_material(adjncy_g(el_adj)+1)) < TINYVAL) then
              is_acoustic_el_adj = .false.
              is_elastic_el_adj = .true.
            else if ( phi_material(num_material(adjncy_g(el_adj)+1)) >= 1.d0) then
              is_acoustic_el_adj = .true.
              is_elastic_el_adj = .false.
            else
              is_acoustic_el_adj = .false.
              is_elastic_el_adj = .false.
            endif
            if ( (part(adjncy_g(el_adj)) == num_part_bis) .and. &
                (is_acoustic_el .eqv. is_acoustic_el_adj) .and. &
                (is_elastic_el .eqv. is_elastic_el_adj) ) then
              tab_interfaces(tab_size_interfaces(num_interface)*5+num_edge*5+0) = el
              tab_interfaces(tab_size_interfaces(num_interface)*5+num_edge*5+1) = adjncy_g(el_adj)
              ncommon_nodes = 0
              do num_node = 0, 4-1
                do num_node_bis = 0, 4-1
                  if ( elmnts_l(el*NCORNERS+num_node) == &
                      elmnts_l(adjncy_g(el_adj)*NCORNERS+num_node_bis) ) then
                    tab_interfaces(tab_size_interfaces(num_interface)*5+num_edge*5+3+ncommon_nodes) &
                                = elmnts_l(el*NCORNERS+num_node)
                    ncommon_nodes = ncommon_nodes + 1
                  endif
                enddo
              enddo
              if ( ncommon_nodes > 0 ) then
                tab_interfaces(tab_size_interfaces(num_interface)*5+num_edge*5+2) = ncommon_nodes
              else
                print *, "Error while building interfaces!", ncommon_nodes
                stop 'fatal error'
              endif
              num_edge = num_edge + 1
            endif
          enddo
        endif

      enddo
      num_edge = 0
      num_interface = num_interface + 1
    enddo
  enddo

  end subroutine Construct_interfaces


  !--------------------------------------------------
  ! Write nodes (their coordinates) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_glob2loc_nodes_database(IIN_database, iproc, npgeo, num_phase)

  implicit none

  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: iproc, num_phase
  integer, intent(inout)  :: npgeo

  integer  :: i, j

  if ( num_phase == 1 ) then
     npgeo = 0

     do i = 0, nnodes-1
        do j = glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
           if ( glob2loc_nodes_parts(j) == iproc ) then
              npgeo = npgeo + 1
           endif
        enddo
     enddo
  else
     do i = 0, nnodes-1
        do j = glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
           if ( glob2loc_nodes_parts(j) == iproc ) then
              write(IIN_database,*) glob2loc_nodes(j)+1, nodes_coords(1,i+1), nodes_coords(2,i+1)
           endif
        enddo
     enddo
  endif

  end subroutine Write_glob2loc_nodes_database


  !--------------------------------------------------
  ! Write elements (their nodes) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_partition_database(IIN_database, iproc, nspec, &
                                      num_modele, num_pml, ngnod, num_phase)

  implicit none

  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: num_phase, iproc
  integer, intent(inout)  :: nspec
  integer, dimension(:)  :: num_modele
  integer, dimension(:)  :: num_pml
  integer, intent(in)  :: ngnod

  integer  :: i,j,k
  integer, dimension(0:ngnod-1)  :: loc_nodes

  if (num_phase == 1) then

     nspec = 0

     do i = 0, nelmnts-1
        if (part(i) == iproc) nspec = nspec + 1
     enddo

  else
     do i = 0, nelmnts-1
        if (part(i) == iproc) then

           do j = 0, ngnod-1
              do k = glob2loc_nodes_nparts(elmnts(i*ngnod+j)), glob2loc_nodes_nparts(elmnts(i*ngnod+j)+1)-1
                 if (glob2loc_nodes_parts(k) == iproc) loc_nodes(j) = glob2loc_nodes(k)
              enddo
           enddo
           write(IIN_database,*) glob2loc_elmnts(i)+1, num_modele(i+1), (loc_nodes(k)+1, k=0,ngnod-1), num_pml(i+1)
        endif
     enddo

  endif

  end subroutine write_partition_database


  !--------------------------------------------------
  ! Write interfaces (element and common nodes) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine Write_interfaces_database(IIN_database, nparts, iproc, &
                        my_ninterface, my_interfaces, my_nb_interfaces, num_phase)

  implicit none

  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: iproc
  integer, intent(in)  :: nparts
  integer, intent(inout)  :: my_ninterface
  integer, dimension(0:ninterfaces-1), intent(inout)  :: my_interfaces
  integer, dimension(0:ninterfaces-1), intent(inout)  :: my_nb_interfaces

  integer, dimension(2)  :: local_nodes
  integer  :: local_elmnt
  integer  :: num_phase

  integer  :: i, j, k, l
  integer  :: num_interface

  num_interface = 0

  if ( num_phase == 1 ) then

     my_interfaces(:) = 0
     my_nb_interfaces(:) = 0

     do i = 0, nparts-1
        do j = i+1, nparts-1
           if ( (tab_size_interfaces(num_interface) < tab_size_interfaces(num_interface+1)) .and. &
                (i == iproc .or. j == iproc) ) then
              my_interfaces(num_interface) = 1
              my_nb_interfaces(num_interface) = tab_size_interfaces(num_interface+1) &
                                              - tab_size_interfaces(num_interface)
           endif
           num_interface = num_interface + 1
        enddo
     enddo
     my_ninterface = sum(my_interfaces(:))

  else

    do i = 0, nparts-1
      do j = i+1, nparts-1
        if ( my_interfaces(num_interface) == 1 ) then
          if ( i == iproc ) then
            write(IIN_database,*) j, my_nb_interfaces(num_interface)
          else
            write(IIN_database,*) i, my_nb_interfaces(num_interface)
          endif

          do k = tab_size_interfaces(num_interface), tab_size_interfaces(num_interface+1)-1
            if ( i == iproc ) then
              local_elmnt = glob2loc_elmnts(tab_interfaces(k*5+0))+1
            else
              local_elmnt = glob2loc_elmnts(tab_interfaces(k*5+1))+1
            endif

            if ( tab_interfaces(k*5+2) == 1 ) then
              ! common node (single point)
              do l = glob2loc_nodes_nparts(tab_interfaces(k*5+3)), &
                        glob2loc_nodes_nparts(tab_interfaces(k*5+3)+1)-1
                if ( glob2loc_nodes_parts(l) == iproc ) then
                  local_nodes(1) = glob2loc_nodes(l)+1
                endif
              enddo

              write(IIN_database,*) local_elmnt, tab_interfaces(k*5+2), &
                                        local_nodes(1), -1
            else
              if ( tab_interfaces(k*5+2) == 2 ) then
                ! common edge (two nodes)
                ! first node
                do l = glob2loc_nodes_nparts(tab_interfaces(k*5+3)), &
                           glob2loc_nodes_nparts(tab_interfaces(k*5+3)+1)-1
                  if ( glob2loc_nodes_parts(l) == iproc ) then
                    local_nodes(1) = glob2loc_nodes(l)+1
                  endif
                enddo
                ! second node
                do l = glob2loc_nodes_nparts(tab_interfaces(k*5+4)), &
                         glob2loc_nodes_nparts(tab_interfaces(k*5+4)+1)-1
                  if ( glob2loc_nodes_parts(l) == iproc ) then
                    local_nodes(2) = glob2loc_nodes(l)+1
                  endif
                enddo

                write(IIN_database,*) local_elmnt, tab_interfaces(k*5+2), &
                                           local_nodes(1), local_nodes(2)
              else
                write(IIN_database,*) "erreur_write_interface_", tab_interfaces(k*5+2)
              endif
            endif
          enddo

        endif

        num_interface = num_interface + 1
      enddo
    enddo

  endif

  end subroutine Write_interfaces_database


  !--------------------------------------------------
  ! Write a surface (elements and nodes on the surface) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine Write_surface_database(IIN_database, nsurface, surface, &
                                nsurface_loc, iproc, num_phase)

  implicit none
  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: iproc
  integer  :: nsurface
  integer  :: nsurface_loc
  integer, dimension(:,:), pointer  :: surface

  integer, dimension(2)  :: local_nodes
  integer  :: local_elmnt
  integer  :: num_phase

  integer  :: i, l

  if ( num_phase == 1 ) then

    nsurface_loc = 0

    do i = 1, nsurface
      if ( part(surface(1,i)) == iproc ) then
        nsurface_loc = nsurface_loc + 1
      endif
    enddo

  else

    nsurface_loc = 0

    do i = 1, nsurface
      if ( part(surface(1,i)) == iproc ) then
        nsurface_loc = nsurface_loc + 1

        local_elmnt = glob2loc_elmnts(surface(1,i)) + 1

        if ( surface(2,i) == 1 ) then
          do l = glob2loc_nodes_nparts(surface(3,i)), &
                  glob2loc_nodes_nparts(surface(3,i)+1)-1
            if ( glob2loc_nodes_parts(l) == iproc ) then
              local_nodes(1) = glob2loc_nodes(l)+1
            endif
          enddo

          write(IIN_database,*) local_elmnt, surface(2,i), local_nodes(1), -1
        endif

        if ( surface(2,i) == 2 ) then
          do l = glob2loc_nodes_nparts(surface(3,i)), &
                  glob2loc_nodes_nparts(surface(3,i)+1)-1
            if ( glob2loc_nodes_parts(l) == iproc ) then
              local_nodes(1) = glob2loc_nodes(l)+1
            endif
          enddo
          do l = glob2loc_nodes_nparts(surface(4,i)), &
                  glob2loc_nodes_nparts(surface(4,i)+1)-1
            if ( glob2loc_nodes_parts(l) == iproc ) then
              local_nodes(2) = glob2loc_nodes(l)+1
            endif
          enddo

          write(IIN_database,*) local_elmnt, surface(2,i), local_nodes(1), local_nodes(2)
        endif

      endif

    enddo

  endif

  end subroutine Write_surface_database


  !--------------------------------------------------
  ! Set absorbing boundaries by elements instead of edges.
  ! Also excludes first or last GLL integration point
  ! if it has both absorbing condition and coupled fluid/solid relation
  ! (this is the reason why arrays ibegin_..., iend_... are included here).
  ! Under development : excluding points that have two different normals in two different elements.
  !--------------------------------------------------

  subroutine merge_abs_boundaries(nb_materials, phi_material, num_material, ngnod)

  implicit none
  include "constants.h"

  integer, intent(in)  :: ngnod
  integer  :: nb_materials
  double precision, dimension(nb_materials), intent(in)  :: phi_material
  integer, dimension(1:nelmnts), intent(in)  :: num_material

  logical, dimension(nb_materials)  :: is_acoustic
  integer  :: num_edge, nedge_bound
  integer  :: match
  integer  :: nb_elmnts_abs
  integer  :: i
  integer  :: temp
  integer  :: iedge, inode1, inode2

  allocate(abs_surface_char(4,nelemabs))
  allocate(abs_surface_merge(nelemabs))
  allocate(abs_surface_type(nelemabs))
  abs_surface_char(:,:) = .false.
  abs_surface_merge(:) = -1
  abs_surface_type(:) = -1

  nedge_bound = nelemabs
  nb_elmnts_abs = 0

  do num_edge = 1, nedge_bound

!! DK DK Sept 2012: in order to fix the rotated elements issue in external mesh
!! DK DK Sept 2012: we now use a type code and thus we must not merge elements that
!! DK DK Sept 2012: appear twice in the list any more because each occurrence now appears with a different type code
!! DK DK Sept 2012    match = 0
!! DK DK Sept 2012    do i = 1, nb_elmnts_abs
!! DK DK Sept 2012       if ( abs_surface(1,num_edge) == abs_surface_merge(i) ) then
!! DK DK Sept 2012          match = i
!! DK DK Sept 2012          exit
!! DK DK Sept 2012       endif
!! DK DK Sept 2012    enddo

!! DK DK Sept 2012    if ( match == 0 ) then
       nb_elmnts_abs = nb_elmnts_abs + 1
       match = nb_elmnts_abs
!! DK DK Sept 2012    endif

    abs_surface_merge(match) = abs_surface(1,num_edge)
!! DK DK Sept 2012 added the absorbing interface type for Stacey
    abs_surface_type(match) = abs_surface(5,num_edge)

    if ( (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
          abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+1)) ) then
       abs_surface_char(IEDGE1,match) = .true.
    endif

    if ( (abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
          abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+1)) ) then
       temp = abs_surface(4,num_edge)
       abs_surface(4,num_edge) = abs_surface(3,num_edge)
       abs_surface(3,num_edge) = temp
       abs_surface_char(IEDGE1,match) = .true.
    endif

    if ( (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
          abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+3)) ) then
       abs_surface_char(IEDGE4,match) = .true.
    endif

    if ( (abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
          abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+3)) ) then
       temp = abs_surface(4,num_edge)
       abs_surface(4,num_edge) = abs_surface(3,num_edge)
       abs_surface(3,num_edge) = temp
       abs_surface_char(IEDGE4,match) = .true.
    endif

    if ( (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+1) .and. &
          abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+2)) ) then
       abs_surface_char(IEDGE2,match) = .true.
    endif

    if ( (abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+1) .and. &
          abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+2)) ) then
       temp = abs_surface(4,num_edge)
       abs_surface(4,num_edge) = abs_surface(3,num_edge)
       abs_surface(3,num_edge) = temp
       abs_surface_char(IEDGE2,match) = .true.
    endif

    if ( (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+2) .and. &
          abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+3)) ) then
       temp = abs_surface(4,num_edge)
       abs_surface(4,num_edge) = abs_surface(3,num_edge)
       abs_surface(3,num_edge) = temp
       abs_surface_char(IEDGE3,match) = .true.
    endif

    if ( (abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+2) .and. &
          abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+3)) ) then
       abs_surface_char(IEDGE3,match) = .true.
    endif

  enddo

  nelemabs_merge = nb_elmnts_abs

! beware here and below that external meshes (for instance coming from CUBIT or Gmsh)
! may have rotated elements and thus edge 1 may not correspond to the bottom,
! edge 2 may not correspond to the right, edge 3 may not correspond to the top,
! and edge 4 may not correspond to the left.
  allocate(ibegin_edge1(nelemabs_merge))
  allocate(iend_edge1(nelemabs_merge))
  allocate(ibegin_edge2(nelemabs_merge))
  allocate(iend_edge2(nelemabs_merge))
  allocate(ibegin_edge3(nelemabs_merge))
  allocate(iend_edge3(nelemabs_merge))
  allocate(ibegin_edge4(nelemabs_merge))
  allocate(iend_edge4(nelemabs_merge))

  ibegin_edge1(:) = 1
  ibegin_edge2(:) = 1
  ibegin_edge3(:) = 1
  ibegin_edge4(:) = 1
  iend_edge1(:) = NGLLX
  iend_edge2(:) = NGLLZ
  iend_edge3(:) = NGLLX
  iend_edge4(:) = NGLLZ

  is_acoustic(:) = .false.

  do i = 1, nb_materials
     if (phi_material(i) >= 1.d0) then
        is_acoustic(i) = .true.
     endif
  enddo

  do num_edge = 1, nedge_bound

!! DK DK Sept 2012: in order to fix the rotated elements issue in external mesh
!! DK DK Sept 2012: we now use a type code and thus we must not merge elements that
!! DK DK Sept 2012: appear twice in the list any more because each occurrence now appears with a different type code.
!! DK DK Sept 2012: But here I think we must leave it because we are just fixing the fluid/solid elements in postprocessing.
  match = 0
  do i = 1, nelemabs_merge
    if ( abs_surface(1,num_edge) == abs_surface_merge(i) ) then
       match = i
       exit
    endif
  enddo

  if ( is_acoustic(num_material(abs_surface(1,num_edge)+1)) ) then

    do iedge = 1, nedges_coupled

      do inode1 = 0, 3
        if ( abs_surface(3,num_edge) == elmnts(ngnod*edges_coupled(1,iedge)+inode1) ) then
          do inode2 = 0, 3
            if ( abs_surface(3,num_edge) == elmnts(ngnod*edges_coupled(2,iedge)+inode2) ) then

              if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+0) .and. &
                   abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+1) )  then
                  ibegin_edge1(match) = 2
              endif

              if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+1) .and. &
                   abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+2) )  then
                  ibegin_edge2(match) = 2
              endif

              if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+3) .and. &
                   abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+2) )  then
                  ibegin_edge3(match) = 2
              endif

              if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+0) .and. &
                   abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+3) )  then
                  ibegin_edge4(match) = 2
              endif

            endif
          enddo

        endif

        if ( abs_surface(4,num_edge) == elmnts(ngnod*edges_coupled(1,iedge)+inode1) ) then
          do inode2 = 0, 3
            if ( abs_surface(4,num_edge) == elmnts(ngnod*edges_coupled(2,iedge)+inode2) ) then
              if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+0) .and. &
                   abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+1) )  then
                  iend_edge1(match) = NGLLX - 1
              endif

              if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+1) .and. &
                   abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+2) )  then
                  iend_edge2(match) = NGLLZ - 1
              endif

              if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+3) .and. &
                   abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+2) )  then
                  iend_edge3(match) = NGLLX - 1
              endif

              if ( abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+0) .and. &
                   abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+3) )  then
                  iend_edge4(match) = NGLLZ - 1
              endif

            endif
          enddo

        endif

      enddo


    enddo

  endif

  enddo

  end subroutine merge_abs_boundaries


  !--------------------------------------------------
  ! Write abs surface (elements and nodes on the surface) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------

  subroutine write_abs_merge_database(IIN_database, iproc, num_phase)

  implicit none

  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: iproc
  integer, intent(in)  :: num_phase

  integer  :: i

  if ( num_phase == 1 ) then
    nelemabs_loc = 0
    do i = 1, nelemabs_merge
       if ( part(abs_surface_merge(i)) == iproc ) then
          nelemabs_loc = nelemabs_loc + 1
       endif
    enddo
  else
    do i = 1, nelemabs_merge
       if ( part(abs_surface_merge(i)) == iproc ) then

! beware here and below that external meshes (for instance coming from CUBIT or Gmsh)
! may have rotated elements and thus edge 1 may not correspond to the bottom,
! edge 2 may not correspond to the right, edge 3 may not correspond to the top,
! and edge 4 may not correspond to the left.

! we add +1 to the "glob2loc_elmnts" number because internally in parts of our code numbers start at zero instead of one
          write(IIN_database,*) glob2loc_elmnts(abs_surface_merge(i))+1, abs_surface_char(1,i), &
!! DK DK sept 2012: added absorbing element type
!! DK DK sept 2012               abs_surface_char(2,i), abs_surface_char(3,i), abs_surface_char(4,i), &
               abs_surface_char(2,i), abs_surface_char(3,i), abs_surface_char(4,i), abs_surface(5,i), &
               ibegin_edge1(i), iend_edge1(i), &
               ibegin_edge2(i), iend_edge2(i), &
               ibegin_edge3(i), iend_edge3(i), &
               ibegin_edge4(i), iend_edge4(i)

       endif

    enddo
  endif

  end subroutine write_abs_merge_database


!! DK DK support for METIS now removed, we use SCOTCH instead
!#ifdef USE_METIS
! !--------------------------------------------------
! ! Partitioning using METIS
! !--------------------------------------------------
!    subroutine Part_metis(nelmnts, xadj, adjncy, vwgt, adjwgt, nparts, nb_edges, edgecut, part, metis_options)
!
!   include "constants.h"
!
!   integer, intent(in)  :: nelmnts, nparts, nb_edges
!   integer, intent(inout)  :: edgecut
!   integer, dimension(0:nelmnts), intent(in)  :: xadj
!   integer, dimension(0:MAX_NEIGHBORS*nelmnts-1), intent(in)  :: adjncy
!   integer, dimension(0:nelmnts-1), intent(in)  :: vwgt
!   integer, dimension(0:nb_edges-1), intent(in)  :: adjwgt
!   integer, dimension(:), pointer  :: part
!   integer, dimension(0:4)  :: metis_options
!
!   integer  :: wgtflag
!   integer  :: remove_min_to_start_at_zero
!
!   remove_min_to_start_at_zero = 0
!   wgtflag = 0
!
!   call METIS_PartGraphRecursive(nelmnts, xadj(0), adjncy(0), vwgt(0), adjwgt(0), wgtflag, remove_min_to_start_at_zero, nparts, &
!        metis_options, edgecut, part(0));
!   !call METIS_PartGraphVKway(nelmnts, xadj(0), adjncy(0), vwgt(0), adjwgt(0), wgtflag, remove_min_to_start_at_zero, nparts, &
!   !     options, edgecut, part(0));
!
! end subroutine Part_metis
!#endif


#ifdef USE_SCOTCH
  !--------------------------------------------------
  ! Partitioning using SCOTCH
  !--------------------------------------------------
  subroutine Part_scotch(nparts, edgecut)

  implicit none
  include "constants.h"

  include "scotchf.h"

  integer, intent(in)  :: nparts
  integer, intent(inout)  :: edgecut

  double precision, dimension(SCOTCH_GRAPHDIM)  :: SCOTCHGRAPH
  double precision, dimension(SCOTCH_STRATDIM)  :: SCOTCHSTRAT
  integer  :: IERR

  edgecut = vwgt(0)
  edgecut = 0

  ! we use the default strategy for partitioning
  ! thus no need to define an explicit strategy
  call scotchfstratinit (SCOTCHSTRAT(1), IERR)
   IF (IERR .NE. 0) THEN
     PRINT *, 'ERROR : MAIN : Cannot initialize strat'
     STOP
  ENDIF

  CALL SCOTCHFGRAPHINIT (SCOTCHGRAPH (1), IERR)
  IF (IERR .NE. 0) THEN
     PRINT *, 'ERROR : MAIN : Cannot initialize graph'
     STOP
  ENDIF

  ! fills graph structure : see user manual (scotch_user5.1.pdf, page 72/73)
  ! arguments: #(1) graph_structure       #(2) baseval(either 0/1)    #(3) number_of_vertices
  !                    #(4) adjacency_index_array         #(5) adjacency_end_index_array (optional)
  !                    #(6) vertex_load_array (optional) #(7) vertex_label_array
  !                    #(7) number_of_arcs                    #(8) adjacency_array
  !                    #(9) arc_load_array (optional)      #(10) ierror
  CALL SCOTCHFGRAPHBUILD (SCOTCHGRAPH (1), 0, nelmnts, &
                          xadj_g(0), xadj_g(0), &
                          xadj_g(0), xadj_g(0), &
                          nb_edges, &
                          adjncy_g(0), adjwgt (0), IERR)
  IF (IERR .NE. 0) THEN
     PRINT *, 'ERROR : MAIN : Cannot build graph'
     STOP
  ENDIF

  CALL SCOTCHFGRAPHCHECK (SCOTCHGRAPH (1), IERR)
  IF (IERR .NE. 0) THEN
     PRINT *, 'ERROR : MAIN : Invalid check'
     STOP
  ENDIF

  call scotchfgraphpart (SCOTCHGRAPH (1), nparts, SCOTCHSTRAT(1), part(0), IERR)
  IF (IERR .NE. 0) THEN
     PRINT *, 'ERROR : MAIN : Cannot part graph'
     STOP
  ENDIF

  CALL SCOTCHFGRAPHEXIT (SCOTCHGRAPH (1), IERR)
  IF (IERR .NE. 0) THEN
     PRINT *, 'ERROR : MAIN : Cannot destroy graph'
     STOP
  ENDIF

  call scotchfstratexit (SCOTCHSTRAT(1), IERR)
  IF (IERR .NE. 0) THEN
     PRINT *, 'ERROR : MAIN : Cannot destroy strat'
     STOP
  ENDIF

  end subroutine Part_scotch
#endif


  !--------------------------------------------------
  ! Repartitioning : two coupled acoustic/elastic elements are transfered to the same partition
  !--------------------------------------------------

  subroutine acoustic_elastic_repartitioning (elmnts_l, nb_materials, &
                                          phi_material, num_material, nproc)

  implicit none
  include "constants.h"

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, intent(in)  :: nproc, nb_materials
  double precision, dimension(nb_materials), intent(in)  :: phi_material
  integer, dimension(1:nelmnts), intent(in)  :: num_material

  ! local parameters
  integer, dimension(:), allocatable  :: xadj_l
  integer, dimension(:), allocatable  :: adjncy_l
  logical, dimension(nb_materials)  :: is_acoustic, is_elastic
  integer  :: i, num_edge
  integer  :: el, el_adj
  logical  :: is_repartitioned

  allocate(xadj_l(0:nelmnts))
  allocate(adjncy_l(0:MAX_NEIGHBORS*nelmnts-1))

  is_acoustic(:) = .false.
  is_elastic(:) = .false.

  do i = 1, nb_materials
     if (phi_material(i) >= 1.d0) then
        is_acoustic(i) = .true.
     endif
     if (phi_material(i) < TINYVAL) then
        is_elastic(i) = .true.
     endif
  enddo

  ! determines maximum neighbors based on 2 common nodes (common edge)
  call mesh2dual_ncommonnodes(elmnts_l, 2, xadj_l, adjncy_l)

  nedges_coupled = 0
  do el = 0, nelmnts-1
     if ( is_acoustic(num_material(el+1)) ) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if ( is_elastic(num_material(adjncy_l(el_adj)+1)) ) then
              nedges_coupled = nedges_coupled + 1
           endif
        enddo
     endif
  enddo

  allocate(edges_coupled(2,nedges_coupled))

  nedges_coupled = 0
  do el = 0, nelmnts-1
     if ( is_acoustic(num_material(el+1)) ) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if ( is_elastic(num_material(adjncy_l(el_adj)+1)) ) then
              nedges_coupled = nedges_coupled + 1
              edges_coupled(1,nedges_coupled) = el
              edges_coupled(2,nedges_coupled) = adjncy_l(el_adj)
           endif

        enddo
     endif
  enddo

  do i = 1, nedges_coupled*nproc
     is_repartitioned = .false.
     do num_edge = 1, nedges_coupled
        if ( part(edges_coupled(1,num_edge)) /= part(edges_coupled(2,num_edge)) ) then
           if ( part(edges_coupled(1,num_edge)) < part(edges_coupled(2,num_edge)) ) then
              part(edges_coupled(2,num_edge)) = part(edges_coupled(1,num_edge))
           else
              part(edges_coupled(1,num_edge)) = part(edges_coupled(2,num_edge))
           endif
           is_repartitioned = .true.
        endif

     enddo
     if ( .not. is_repartitioned ) then
        exit
     endif
  enddo

  deallocate(xadj_l,adjncy_l)

  end subroutine acoustic_elastic_repartitioning


  !--------------------------------------------------
  ! Repartitioning : two coupled acoustic/poroelastic elements are transfered to the same partition
  !--------------------------------------------------

  subroutine acoustic_poro_repartitioning (elmnts_l, nb_materials, &
                                        phi_material, num_material, nproc)

  implicit none
  include "constants.h"

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, intent(in)  :: nproc, nb_materials
  double precision, dimension(nb_materials), intent(in)  :: phi_material
  integer, dimension(1:nelmnts), intent(in)  :: num_material

  ! local parameters
  integer, dimension(:), allocatable  :: xadj_l
  integer, dimension(:), allocatable  :: adjncy_l
  logical, dimension(nb_materials)  :: is_acoustic,is_poroelastic
  integer  :: i, num_edge
  integer  :: el, el_adj
  logical  :: is_repartitioned

  allocate(xadj_l(0:nelmnts))
  allocate(adjncy_l(0:MAX_NEIGHBORS*nelmnts-1))

  is_acoustic(:) = .false.
  is_poroelastic(:) = .false.

  do i = 1, nb_materials
     if (phi_material(i) >=1.d0) then
        is_acoustic(i) = .true.
     endif
     if (phi_material(i) <1.d0 .and. phi_material(i) > TINYVAL) then
        is_poroelastic(i) = .true.
     endif
  enddo

  ! determines maximum neighbors based on 2 common nodes (common edge)
  call mesh2dual_ncommonnodes(elmnts_l, 2, xadj_l, adjncy_l)

  nedges_acporo_coupled = 0
  do el = 0, nelmnts-1
     if ( is_acoustic(num_material(el+1)) ) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if ( is_poroelastic(num_material(adjncy_l(el_adj)+1)) ) then
              nedges_acporo_coupled = nedges_acporo_coupled + 1
           endif

        enddo
     endif
  enddo

  print *, 'nedges_coupled (acoustic/poroelastic)', nedges_acporo_coupled

  allocate(edges_acporo_coupled(2,nedges_acporo_coupled))

  nedges_acporo_coupled = 0
  do el = 0, nelmnts-1
     if ( is_acoustic(num_material(el+1)) ) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if ( is_poroelastic(num_material(adjncy_l(el_adj)+1)) ) then
              nedges_acporo_coupled = nedges_acporo_coupled + 1
              edges_acporo_coupled(1,nedges_acporo_coupled) = el
              edges_acporo_coupled(2,nedges_acporo_coupled) = adjncy_l(el_adj)
           endif

        enddo
     endif
  enddo

  do i = 1, nedges_acporo_coupled*nproc
     is_repartitioned = .false.
     do num_edge = 1, nedges_acporo_coupled
        if ( part(edges_acporo_coupled(1,num_edge)) /= part(edges_acporo_coupled(2,num_edge)) ) then
           if ( part(edges_acporo_coupled(1,num_edge)) < part(edges_acporo_coupled(2,num_edge)) ) then
              part(edges_acporo_coupled(2,num_edge)) = part(edges_acporo_coupled(1,num_edge))
           else
              part(edges_acporo_coupled(1,num_edge)) = part(edges_acporo_coupled(2,num_edge))
           endif
           is_repartitioned = .true.
        endif

     enddo
     if ( .not. is_repartitioned ) then
        exit
     endif
  enddo

  deallocate(xadj_l,adjncy_l)

  end subroutine acoustic_poro_repartitioning


  !--------------------------------------------------
  ! Repartitioning : two coupled poroelastic/elastic elements are transfered to the same partition
  !--------------------------------------------------

  subroutine poro_elastic_repartitioning (elmnts_l, nb_materials, &
                                        phi_material, num_material, nproc)

  implicit none
  include "constants.h"

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, intent(in)  :: nproc, nb_materials
  double precision, dimension(nb_materials), intent(in)  :: phi_material
  integer, dimension(1:nelmnts), intent(in)  :: num_material

  ! local parameters
  integer, dimension(:), allocatable  :: xadj_l
  integer, dimension(:), allocatable  :: adjncy_l
  logical, dimension(nb_materials)  :: is_elastic,is_poroelastic
  integer  :: i, num_edge
  integer  :: el, el_adj
  logical  :: is_repartitioned

  allocate(xadj_l(0:nelmnts))
  allocate(adjncy_l(0:MAX_NEIGHBORS*nelmnts-1))

  is_elastic(:) = .false.
  is_poroelastic(:) = .false.

  do i = 1, nb_materials
     if (phi_material(i) < TINYVAL) then
        is_elastic(i) = .true.
     endif
     if (phi_material(i) <1.d0 .and. phi_material(i) > TINYVAL) then
        is_poroelastic(i) = .true.
     endif
  enddo

  ! determines maximum neighbors based on 2 common nodes (common edge)
  call mesh2dual_ncommonnodes(elmnts_l, 2, xadj_l, adjncy_l)

  nedges_elporo_coupled = 0
  do el = 0, nelmnts-1
     if ( is_poroelastic(num_material(el+1)) ) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if ( is_elastic(num_material(adjncy_l(el_adj)+1)) ) then
              nedges_elporo_coupled = nedges_elporo_coupled + 1
           endif

        enddo
     endif
  enddo

  print *, 'nedges_coupled (poroelastic/elastic)', nedges_elporo_coupled

  allocate(edges_elporo_coupled(2,nedges_elporo_coupled))

  nedges_elporo_coupled = 0
  do el = 0, nelmnts-1
     if ( is_poroelastic(num_material(el+1)) ) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if ( is_elastic(num_material(adjncy_l(el_adj)+1)) ) then
              nedges_elporo_coupled = nedges_elporo_coupled + 1
              edges_elporo_coupled(1,nedges_elporo_coupled) = el
              edges_elporo_coupled(2,nedges_elporo_coupled) = adjncy_l(el_adj)
           endif

        enddo
     endif
  enddo

  do i = 1, nedges_elporo_coupled*nproc
     is_repartitioned = .false.
     do num_edge = 1, nedges_elporo_coupled
        if ( part(edges_elporo_coupled(1,num_edge)) /= part(edges_elporo_coupled(2,num_edge)) ) then
           if ( part(edges_elporo_coupled(1,num_edge)) < part(edges_elporo_coupled(2,num_edge)) ) then
              part(edges_elporo_coupled(2,num_edge)) = part(edges_elporo_coupled(1,num_edge))
           else
              part(edges_elporo_coupled(1,num_edge)) = part(edges_elporo_coupled(2,num_edge))
           endif
           is_repartitioned = .true.
        endif

     enddo
     if ( .not. is_repartitioned ) then
        exit
     endif
  enddo

  deallocate(xadj_l,adjncy_l)

  end subroutine poro_elastic_repartitioning


  !--------------------------------------------------
  ! Write fluid/solid edges (fluid (or porous) elements and corresponding solid (or porous) elements)
  ! pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------

 subroutine write_fluidsolid_edges_database(IIN_database, nedges_coupled_bis, nedges_coupled_loc_bis, &
                                            edges_coupled_bis, iproc, num_phase)

  implicit none

  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: nedges_coupled_bis
  integer, intent(inout)  :: nedges_coupled_loc_bis
  integer, dimension(:,:), pointer  :: edges_coupled_bis
  integer, intent(in)  :: iproc
  integer, intent(in)  :: num_phase

  integer  :: i

  if ( num_phase == 1 ) then
     nedges_coupled_loc_bis = 0
     do i = 1, nedges_coupled_bis
        if ( part(edges_coupled_bis(1,i)) == iproc ) then
           nedges_coupled_loc_bis = nedges_coupled_loc_bis + 1
        endif
     enddo
  else
     do i = 1, nedges_coupled_bis
        if ( part(edges_coupled_bis(1,i)) == iproc ) then
           write(IIN_database,*) glob2loc_elmnts(edges_coupled_bis(1,i))+1, glob2loc_elmnts(edges_coupled_bis(2,i))+1
        endif
     enddo
  endif

  end subroutine write_fluidsolid_edges_database

end module part_unstruct
