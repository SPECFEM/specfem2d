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
  ! Creating dual graph (adjacency is defined by 'ncommonnodes' between two elements).
  !-----------------------------------------------
  subroutine mesh2dual_ncommonnodes(elmnts_l,ncommonnodes,xadj,adjncy)

  use constants, only: NCORNERS,MAX_NEIGHBORS,MAX_NSIZE_SHARED

  use part_unstruct_par, only: nelmnts,nnodes,nnodes_elmnts,nodes_elmnts

  implicit none

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, intent(in) :: ncommonnodes
  integer, dimension(0:nelmnts),intent(out)  :: xadj
  integer, dimension(0:MAX_NEIGHBORS*nelmnts-1),intent(out) :: adjncy

  ! local parameters
  integer  :: i, j, k, l, m, num_edges
  logical  ::  is_neighbor
  integer  :: num_node, n
  integer  :: elem_base, elem_target
  integer  :: connectivity

  ! allocates memory for arrays
  if (.not. allocated(nnodes_elmnts) ) allocate(nnodes_elmnts(0:nnodes-1))
  if (.not. allocated(nodes_elmnts) ) allocate(nodes_elmnts(0:MAX_NSIZE_SHARED*nnodes-1))

  ! initializes
  xadj(:) = 0
  adjncy(:) = 0
  nnodes_elmnts(:) = 0
  nodes_elmnts(:) = 0
  num_edges = 0

  ! list of elements per node
  do i = 0, NCORNERS*nelmnts-1
    nodes_elmnts(elmnts_l(i)*MAX_NSIZE_SHARED + nnodes_elmnts(elmnts_l(i))) = i/NCORNERS
    nnodes_elmnts(elmnts_l(i)) = nnodes_elmnts(elmnts_l(i)) + 1
  enddo

  ! checking which elements are neighbors ('ncommonnodes' criteria)
  do j = 0, nnodes-1
    do k = 0, nnodes_elmnts(j)-1
      do l = k+1, nnodes_elmnts(j)-1

        connectivity = 0
        elem_base = nodes_elmnts(k+j*MAX_NSIZE_SHARED)
        elem_target = nodes_elmnts(l+j*MAX_NSIZE_SHARED)
        do n = 1, NCORNERS
          num_node = elmnts_l(NCORNERS*elem_base+n-1)
          do m = 0, nnodes_elmnts(num_node)-1
            if (nodes_elmnts(m+num_node*MAX_NSIZE_SHARED) == elem_target) then
              connectivity = connectivity + 1
            endif
          enddo
        enddo

        ! sets adjacency (adjncy) and number of neighbors (xadj)
        ! according to ncommonnodes criteria
        if (connectivity >= ncommonnodes) then

          is_neighbor = .false.

          do m = 0, xadj(nodes_elmnts(k+j*MAX_NSIZE_SHARED))
            if (.not. is_neighbor) then
              if (adjncy(nodes_elmnts(k+j*MAX_NSIZE_SHARED)*MAX_NEIGHBORS+m) == nodes_elmnts(l+j*MAX_NSIZE_SHARED)) then
                is_neighbor = .true.
              endif
            endif
          enddo
          if (.not. is_neighbor) then
            adjncy(nodes_elmnts(k+j*MAX_NSIZE_SHARED)*MAX_NEIGHBORS &
                   + xadj(nodes_elmnts(k+j*MAX_NSIZE_SHARED))) = nodes_elmnts(l+j*MAX_NSIZE_SHARED)

            xadj(nodes_elmnts(k+j*MAX_NSIZE_SHARED)) = xadj(nodes_elmnts(k+j*MAX_NSIZE_SHARED)) + 1
            if (xadj(nodes_elmnts(k+j*MAX_NSIZE_SHARED)) > MAX_NEIGHBORS) &
              call stop_the_code('ERROR : too much neighbors per element, modify the mesh.')

            adjncy(nodes_elmnts(l+j*MAX_NSIZE_SHARED)*MAX_NEIGHBORS &
                   + xadj(nodes_elmnts(l+j*MAX_NSIZE_SHARED))) = nodes_elmnts(k+j*MAX_NSIZE_SHARED)

            xadj(nodes_elmnts(l+j*MAX_NSIZE_SHARED)) = xadj(nodes_elmnts(l+j*MAX_NSIZE_SHARED)) + 1
            if (xadj(nodes_elmnts(l+j*MAX_NSIZE_SHARED)) > MAX_NEIGHBORS) &
              call stop_the_code('ERROR : too much neighbors per element, modify the mesh.')

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

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! construct local numbering for the elements in each partition
  !--------------------------------------------------
  subroutine construct_glob2loc_elmnts(nparts)

  use part_unstruct_par, only: glob2loc_elmnts,nelmnts,part

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

  end subroutine construct_glob2loc_elmnts

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! construct local numbering for the nodes in each partition
  !--------------------------------------------------
  subroutine construct_glob2loc_nodes(nparts)

  use constants, only: MAX_NSIZE_SHARED
  use part_unstruct_par, only: nnodes,glob2loc_nodes_nparts,glob2loc_nodes_parts,glob2loc_nodes, &
    nodes_elmnts,nnodes_elmnts,part

  implicit none

  integer, intent(in)  :: nparts

  integer :: num_node
  integer :: el
  integer ::  num_part
  integer ::  size_glob2loc_nodes
  integer, dimension(0:nparts-1) :: parts_node
  integer, dimension(0:nparts-1) :: num_parts

  allocate(glob2loc_nodes_nparts(0:nnodes))

  size_glob2loc_nodes = 0

  parts_node(:) = 0


  do num_node = 0, nnodes-1
     glob2loc_nodes_nparts(num_node) = size_glob2loc_nodes
     do el = 0, nnodes_elmnts(num_node)-1
        parts_node(part(nodes_elmnts(el+MAX_NSIZE_SHARED*num_node))) = 1
     enddo

     do num_part = 0, nparts-1
        if (parts_node(num_part) == 1) then
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
        parts_node(part(nodes_elmnts(el+MAX_NSIZE_SHARED*num_node))) = 1
     enddo
     do num_part = 0, nparts-1

        if (parts_node(num_part) == 1) then
           glob2loc_nodes_parts(size_glob2loc_nodes) = num_part
           glob2loc_nodes(size_glob2loc_nodes) = num_parts(num_part)
           size_glob2loc_nodes = size_glob2loc_nodes + 1
           num_parts(num_part) = num_parts(num_part) + 1
           parts_node(num_part) = 0
        endif

     enddo
  enddo

  end subroutine construct_glob2loc_nodes

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! Construct interfaces between each partitions.
  ! Two adjacent elements in distinct partitions make an entry in array tab_interfaces :
  ! 1/ first element, 2/ second element, 3/ number of common nodes, 4/ first node,
  ! 5/ second node, if relevant.
  ! No interface between acoustic, elastic, and poroelastic elements.
  !--------------------------------------------------
  subroutine construct_interfaces(nparts, elmnts_l, &
                                  nbmodels, phi_material, num_material)

  use constants, only: NCORNERS,TINYVAL

  use part_unstruct_par, only: nelmnts,ninterfaces,tab_size_interfaces,tab_interfaces,part, &
    xadj_g,adjncy_g

  implicit none

  integer, intent(in)  :: nparts
  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, dimension(1:nelmnts), intent(in)  :: num_material
  integer, intent(in)  :: nbmodels
  double precision, dimension(1:nbmodels), intent(in)  :: phi_material

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
           if (part(el) == num_part) then
              ! sets material flag
              if (phi_material(num_material(el+1)) < TINYVAL) then
                ! elastic element
                is_acoustic_el = .false.
                is_elastic_el = .true.
              else if (phi_material(num_material(el+1)) >= 1.d0) then
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
                if (phi_material(num_material(adjncy_g(el_adj)+1)) < TINYVAL) then
                  is_acoustic_el_adj = .false.
                  is_elastic_el_adj = .true.
                else if (phi_material(num_material(adjncy_g(el_adj)+1)) >= 1.d0) then
                  is_acoustic_el_adj = .true.
                  is_elastic_el_adj = .false.
                else
                  is_acoustic_el_adj = .false.
                  is_elastic_el_adj = .false.
                endif
                ! adds element if neighbor element lies in next partition
                ! and belongs to same material
                if ((part(adjncy_g(el_adj)) == num_part_bis) .and. &
                     (is_acoustic_el .eqv. is_acoustic_el_adj) .and. &
                     (is_elastic_el .eqv. is_elastic_el_adj)) then
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
        if (part(el) == num_part) then
          if (phi_material(num_material(el+1)) < TINYVAL) then
            is_acoustic_el = .false.
            is_elastic_el = .true.
          else if (phi_material(num_material(el+1)) >= 1.d0) then
            is_acoustic_el = .true.
            is_elastic_el = .false.
          else
            is_acoustic_el = .false.
            is_elastic_el = .false.
          endif
          do el_adj = xadj_g(el), xadj_g(el+1)-1
            if (phi_material(num_material(adjncy_g(el_adj)+1)) < TINYVAL) then
              is_acoustic_el_adj = .false.
              is_elastic_el_adj = .true.
            else if (phi_material(num_material(adjncy_g(el_adj)+1)) >= 1.d0) then
              is_acoustic_el_adj = .true.
              is_elastic_el_adj = .false.
            else
              is_acoustic_el_adj = .false.
              is_elastic_el_adj = .false.
            endif
            if ((part(adjncy_g(el_adj)) == num_part_bis) .and. &
                (is_acoustic_el .eqv. is_acoustic_el_adj) .and. &
                (is_elastic_el .eqv. is_elastic_el_adj)) then
              tab_interfaces(tab_size_interfaces(num_interface)*5+num_edge*5+0) = el
              tab_interfaces(tab_size_interfaces(num_interface)*5+num_edge*5+1) = adjncy_g(el_adj)
              ncommon_nodes = 0
              do num_node = 0, 4-1
                do num_node_bis = 0, 4-1
                  if (elmnts_l(el*NCORNERS+num_node) == &
                      elmnts_l(adjncy_g(el_adj)*NCORNERS+num_node_bis)) then
                    tab_interfaces(tab_size_interfaces(num_interface)*5+num_edge*5+3+ncommon_nodes) &
                                = elmnts_l(el*NCORNERS+num_node)
                    ncommon_nodes = ncommon_nodes + 1
                  endif
                enddo
              enddo
              if (ncommon_nodes > 0) then
                tab_interfaces(tab_size_interfaces(num_interface)*5+num_edge*5+2) = ncommon_nodes
              else
                print *,'Error while building interfaces! invalid number of common nodes ', ncommon_nodes
                call stop_the_code('Error building interfaces')
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

  end subroutine construct_interfaces

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! Write nodes (their coordinates) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_glob2loc_nodes_database(IIN_database, iproc, npgeo, num_phase)

  use part_unstruct_par, only: nnodes,glob2loc_nodes_nparts,glob2loc_nodes_parts,glob2loc_nodes,nodes_coords

  implicit none

  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: iproc, num_phase
  integer, intent(inout)  :: npgeo

  integer  :: i, j

  if (num_phase == 1) then
     npgeo = 0

     do i = 0, nnodes-1
        do j = glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
           if (glob2loc_nodes_parts(j) == iproc) then
              npgeo = npgeo + 1
           endif
        enddo
     enddo
  else
     do i = 0, nnodes-1
        do j = glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
           if (glob2loc_nodes_parts(j) == iproc) then
              write(IIN_database) glob2loc_nodes(j)+1, nodes_coords(1,i+1), nodes_coords(2,i+1)
           endif
        enddo
     enddo
  endif

  end subroutine write_glob2loc_nodes_database

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! Write elements (their nodes) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_partition_database(IIN_database, iproc, nspec, &
                                      num_modele, num_pml, ngnod, num_phase)

  use part_unstruct_par, only: nelmnts,elmnts,part, &
    glob2loc_nodes_nparts,glob2loc_nodes_parts,glob2loc_nodes,glob2loc_elmnts

  implicit none

  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: num_phase, iproc
  integer, intent(inout)  :: nspec
  integer, dimension(1:nelmnts)  :: num_modele
  integer, dimension(1:nelmnts)  :: num_pml
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
           write(IIN_database) glob2loc_elmnts(i)+1, num_modele(i+1), (loc_nodes(k)+1, k=0,ngnod-1), num_pml(i+1)
        endif
     enddo

  endif

  end subroutine write_partition_database

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! Write interfaces (element and common nodes) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_interfaces_database(IIN_database, nparts, iproc, &
                        my_ninterface, my_interfaces, my_nb_interfaces, num_phase)

  use part_unstruct_par, only: ninterfaces,tab_size_interfaces,tab_interfaces, &
    glob2loc_elmnts,glob2loc_nodes_nparts,glob2loc_nodes_parts,glob2loc_nodes

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
  integer, parameter :: minus_one = -1

  num_interface = 0

  if (num_phase == 1) then

     my_interfaces(:) = 0
     my_nb_interfaces(:) = 0

     do i = 0, nparts-1
        do j = i+1, nparts-1
           if ((tab_size_interfaces(num_interface) < tab_size_interfaces(num_interface+1)) .and. &
                (i == iproc .or. j == iproc)) then
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
        if (my_interfaces(num_interface) == 1) then
          if (i == iproc) then
            write(IIN_database) j, my_nb_interfaces(num_interface)
          else
            write(IIN_database) i, my_nb_interfaces(num_interface)
          endif

          do k = tab_size_interfaces(num_interface), tab_size_interfaces(num_interface+1)-1
            if (i == iproc) then
              local_elmnt = glob2loc_elmnts(tab_interfaces(k*5+0))+1
            else
              local_elmnt = glob2loc_elmnts(tab_interfaces(k*5+1))+1
            endif

            if (tab_interfaces(k*5+2) == 1) then
              ! common node (single point)
              do l = glob2loc_nodes_nparts(tab_interfaces(k*5+3)), &
                        glob2loc_nodes_nparts(tab_interfaces(k*5+3)+1)-1
                if (glob2loc_nodes_parts(l) == iproc) then
                  local_nodes(1) = glob2loc_nodes(l)+1
                endif
              enddo

              write(IIN_database) local_elmnt, tab_interfaces(k*5+2), &
                                        local_nodes(1), minus_one
            else
              if (tab_interfaces(k*5+2) == 2) then
                ! common edge (two nodes)
                ! first node
                do l = glob2loc_nodes_nparts(tab_interfaces(k*5+3)), &
                           glob2loc_nodes_nparts(tab_interfaces(k*5+3)+1)-1
                  if (glob2loc_nodes_parts(l) == iproc) then
                    local_nodes(1) = glob2loc_nodes(l)+1
                  endif
                enddo
                ! second node
                do l = glob2loc_nodes_nparts(tab_interfaces(k*5+4)), &
                         glob2loc_nodes_nparts(tab_interfaces(k*5+4)+1)-1
                  if (glob2loc_nodes_parts(l) == iproc) then
                    local_nodes(2) = glob2loc_nodes(l)+1
                  endif
                enddo

                write(IIN_database) local_elmnt, tab_interfaces(k*5+2), &
                                           local_nodes(1), local_nodes(2)
              else
                print *,"Error: write_interfaces_database", tab_interfaces(k*5+2)
                call stop_the_code('Invalid interface found')
              endif
            endif
          enddo

        endif

        num_interface = num_interface + 1
      enddo
    enddo

  endif

  end subroutine write_interfaces_database

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! Write a surface (elements and nodes on the surface) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------
  subroutine write_surface_database(IIN_database, nsurface, surface, &
                                    nsurface_loc, iproc, num_phase)

  use part_unstruct_par, only: part,glob2loc_elmnts,glob2loc_nodes_nparts,glob2loc_nodes_parts,glob2loc_nodes

  implicit none
  integer, intent(in)  :: IIN_database
  integer, intent(in)  :: iproc
  integer  :: nsurface
  integer  :: nsurface_loc
  integer, dimension(4,nsurface)  :: surface

  integer, dimension(2)  :: local_nodes
  integer  :: local_elmnt
  integer  :: num_phase

  integer  :: i, l
  integer, parameter :: minus_one = -1

  if (num_phase == 1) then

    ! only counts surface elements
    nsurface_loc = 0

    do i = 1, nsurface
      if (part(surface(1,i)) == iproc) then
        nsurface_loc = nsurface_loc + 1
      endif
    enddo

  else

    nsurface_loc = 0

    do i = 1, nsurface
      if (part(surface(1,i)) == iproc) then
        nsurface_loc = nsurface_loc + 1

        local_elmnt = glob2loc_elmnts(surface(1,i)) + 1

        if (surface(2,i) == 1) then
          do l = glob2loc_nodes_nparts(surface(3,i)), &
                  glob2loc_nodes_nparts(surface(3,i)+1)-1
            if (glob2loc_nodes_parts(l) == iproc) then
              local_nodes(1) = glob2loc_nodes(l)+1
            endif
          enddo

          write(IIN_database) local_elmnt, surface(2,i), local_nodes(1), minus_one
        endif

        if (surface(2,i) == 2) then
          do l = glob2loc_nodes_nparts(surface(3,i)), &
                  glob2loc_nodes_nparts(surface(3,i)+1)-1
            if (glob2loc_nodes_parts(l) == iproc) then
              local_nodes(1) = glob2loc_nodes(l)+1
            endif
          enddo
          do l = glob2loc_nodes_nparts(surface(4,i)), &
                  glob2loc_nodes_nparts(surface(4,i)+1)-1
            if (glob2loc_nodes_parts(l) == iproc) then
              local_nodes(2) = glob2loc_nodes(l)+1
            endif
          enddo

          write(IIN_database) local_elmnt, surface(2,i), local_nodes(1), local_nodes(2)
        endif

      endif

    enddo

  endif

  end subroutine write_surface_database

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! Set absorbing boundaries by elements instead of edges.
  ! Also excludes first or last GLL integration point
  ! if it has both absorbing condition and coupled fluid/solid relation
  ! (this is the reason why arrays ibegin_..., iend_... are included here).
  ! Under development : excluding points that have two different normals in two different elements.
  !--------------------------------------------------

  subroutine merge_abs_boundaries(nbmodels, phi_material, num_material, ngnod)

  use constants, only: IEDGE1,IEDGE2,IEDGE3,IEDGE4,NGLLX,NGLLZ

  use part_unstruct_par, only: nelmnts,elmnts,nelemabs,nelemabs_merge,abs_surface, &
    abs_surface_char,abs_surface_merge,abs_surface_merge,abs_surface_type, &
    ibegin_edge1,iend_edge1,ibegin_edge2,iend_edge2,ibegin_edge3,iend_edge3,ibegin_edge4,iend_edge4, &
    nedges_coupled,edges_coupled

  implicit none

  integer, intent(in)  :: ngnod
  integer  :: nbmodels
  double precision, dimension(nbmodels), intent(in)  :: phi_material
  integer, dimension(1:nelmnts), intent(in)  :: num_material

  logical, dimension(nbmodels)  :: is_acoustic
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
!! DK DK Sept 2012       if (abs_surface(1,num_edge) == abs_surface_merge(i)) then
!! DK DK Sept 2012          match = i
!! DK DK Sept 2012          exit
!! DK DK Sept 2012       endif
!! DK DK Sept 2012    enddo

!! DK DK Sept 2012    if (match == 0) then
    nb_elmnts_abs = nb_elmnts_abs + 1
    match = nb_elmnts_abs
!! DK DK Sept 2012    endif

    abs_surface_merge(match) = abs_surface(1,num_edge)
!! DK DK Sept 2012 added the absorbing interface type for Stacey
    abs_surface_type(match) = abs_surface(5,num_edge)

    if ((abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
          abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+1))) then
       abs_surface_char(IEDGE1,match) = .true.
    endif

    if ((abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
          abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+1))) then
       temp = abs_surface(4,num_edge)
       abs_surface(4,num_edge) = abs_surface(3,num_edge)
       abs_surface(3,num_edge) = temp
       abs_surface_char(IEDGE1,match) = .true.
    endif

    if ((abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
          abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+3))) then
       abs_surface_char(IEDGE4,match) = .true.
    endif

    if ((abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+0) .and. &
          abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+3))) then
       temp = abs_surface(4,num_edge)
       abs_surface(4,num_edge) = abs_surface(3,num_edge)
       abs_surface(3,num_edge) = temp
       abs_surface_char(IEDGE4,match) = .true.
    endif

    if ((abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+1) .and. &
          abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+2))) then
       abs_surface_char(IEDGE2,match) = .true.
    endif

    if ((abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+1) .and. &
          abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+2))) then
       temp = abs_surface(4,num_edge)
       abs_surface(4,num_edge) = abs_surface(3,num_edge)
       abs_surface(3,num_edge) = temp
       abs_surface_char(IEDGE2,match) = .true.
    endif

    if ((abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+2) .and. &
          abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+3))) then
       temp = abs_surface(4,num_edge)
       abs_surface(4,num_edge) = abs_surface(3,num_edge)
       abs_surface(3,num_edge) = temp
       abs_surface_char(IEDGE3,match) = .true.
    endif

    if ((abs_surface(4,num_edge) == elmnts(ngnod*abs_surface_merge(match)+2) .and. &
          abs_surface(3,num_edge) == elmnts(ngnod*abs_surface_merge(match)+3))) then
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

  do i = 1, nbmodels
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
      if (abs_surface(1,num_edge) == abs_surface_merge(i)) then
        match = i
        exit
      endif
    enddo

    if (is_acoustic(num_material(abs_surface(1,num_edge)+1))) then

      do iedge = 1, nedges_coupled
        do inode1 = 0, 3

          if (abs_surface(3,num_edge) == elmnts(ngnod*edges_coupled(1,iedge)+inode1)) then
            do inode2 = 0, 3
              if (abs_surface(3,num_edge) == elmnts(ngnod*edges_coupled(2,iedge)+inode2)) then

                if (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+0) .and. &
                     abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+1)) then
                    ibegin_edge1(match) = 2
                endif

                if (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+1) .and. &
                     abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+2)) then
                    ibegin_edge2(match) = 2
                endif

                if (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+3) .and. &
                     abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+2)) then
                    ibegin_edge3(match) = 2
                endif

                if (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+0) .and. &
                     abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+3)) then
                    ibegin_edge4(match) = 2
                endif

              endif
            enddo

          endif

          if (abs_surface(4,num_edge) == elmnts(ngnod*edges_coupled(1,iedge)+inode1)) then
            do inode2 = 0, 3
              if (abs_surface(4,num_edge) == elmnts(ngnod*edges_coupled(2,iedge)+inode2)) then
                if (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+0) .and. &
                     abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+1)) then
                    iend_edge1(match) = NGLLX - 1
                endif

                if (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+1) .and. &
                     abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+2)) then
                    iend_edge2(match) = NGLLZ - 1
                endif

                if (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+3) .and. &
                     abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+2)) then
                    iend_edge3(match) = NGLLX - 1
                endif

                if (abs_surface(3,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+0) .and. &
                     abs_surface(4,num_edge) == elmnts(ngnod*abs_surface(1,num_edge)+3)) then
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

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! Write abs surface (elements and nodes on the surface) pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------

  subroutine write_abs_merge_database(IIN_database, iproc, num_phase)

  use part_unstruct_par, only: part,nelemabs_loc,nelemabs_merge,abs_surface_merge,abs_surface_char,abs_surface, &
    glob2loc_elmnts,ibegin_edge1,iend_edge1,ibegin_edge2,iend_edge2,ibegin_edge3,iend_edge3,ibegin_edge4,iend_edge4

  implicit none

  integer, intent(in) :: IIN_database
  integer, intent(in) :: iproc
  integer, intent(in) :: num_phase

  integer  :: i

  if (num_phase == 1) then
    ! only counts elements in this partition
    nelemabs_loc = 0
    do i = 1, nelemabs_merge
       if (part(abs_surface_merge(i)) == iproc) then
          nelemabs_loc = nelemabs_loc + 1
       endif
    enddo
  else
    do i = 1, nelemabs_merge
       if (part(abs_surface_merge(i)) == iproc) then

! beware here and below that external meshes (for instance coming from CUBIT or Gmsh)
! may have rotated elements and thus edge 1 may not correspond to the bottom,
! edge 2 may not correspond to the right, edge 3 may not correspond to the top,
! and edge 4 may not correspond to the left.

! we add +1 to the "glob2loc_elmnts" number because internally in parts of our code numbers start at zero instead of one
          write(IIN_database) glob2loc_elmnts(abs_surface_merge(i))+1, abs_surface_char(1,i), &
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

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! Set acoustic forcing boundaries by elements instead of edges.
  ! inspired by merge_abs_boundaries upper in this file
  !--------------------------------------------------

  subroutine merge_acoustic_forcing_boundaries(ngnod)

  use constants, only: IEDGE1,IEDGE2,IEDGE3,IEDGE4,NGLLX,NGLLZ

  use part_unstruct_par, only: elmnts,nelemacforcing,acforcing_surface, &
    acforcing_surface_char,acforcing_surface_merge,acforcing_surface_type, &
    ibegin_edge1_acforcing,iend_edge1_acforcing,ibegin_edge2_acforcing,iend_edge2_acforcing, &
    ibegin_edge3_acforcing,iend_edge3_acforcing,ibegin_edge4_acforcing,iend_edge4_acforcing, &
    nelemacforcing_merge

  implicit none

  integer, intent(in) :: ngnod

  ! local parameters
  integer  :: num_edge, nedge_bound
  integer  :: match
  integer  :: nb_elmnts_acforcing
  integer  :: temp

  allocate(acforcing_surface_char(4,nelemacforcing))
  allocate(acforcing_surface_merge(nelemacforcing))
  allocate(acforcing_surface_type(nelemacforcing))

  acforcing_surface_char(:,:) = .false.
  acforcing_surface_merge(:) = -1
  acforcing_surface_type(:) = -1

  nedge_bound = nelemacforcing
  nb_elmnts_acforcing = 0
  do num_edge = 1, nedge_bound

    nb_elmnts_acforcing = nb_elmnts_acforcing + 1
    match = nb_elmnts_acforcing

    acforcing_surface_merge(match) = acforcing_surface(1,num_edge)
    acforcing_surface_type(match) = acforcing_surface(5,num_edge)

    if ((acforcing_surface(3,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+0) .and. &
          acforcing_surface(4,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+1))) then
       acforcing_surface_char(IEDGE1,match) = .true.
    endif

    if ((acforcing_surface(4,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+0) .and. &
          acforcing_surface(3,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+1))) then
       temp = acforcing_surface(4,num_edge)
       acforcing_surface(4,num_edge) = acforcing_surface(3,num_edge)
       acforcing_surface(3,num_edge) = temp
       acforcing_surface_char(IEDGE1,match) = .true.
    endif

    if ((acforcing_surface(3,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+0) .and. &
          acforcing_surface(4,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+3))) then
       acforcing_surface_char(IEDGE4,match) = .true.
    endif

    if ((acforcing_surface(4,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+0) .and. &
          acforcing_surface(3,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+3))) then
       temp = acforcing_surface(4,num_edge)
       acforcing_surface(4,num_edge) = acforcing_surface(3,num_edge)
       acforcing_surface(3,num_edge) = temp
       acforcing_surface_char(IEDGE4,match) = .true.
    endif

    if ((acforcing_surface(3,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+1) .and. &
          acforcing_surface(4,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+2))) then
       acforcing_surface_char(IEDGE2,match) = .true.
    endif

    if ((acforcing_surface(4,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+1) .and. &
          acforcing_surface(3,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+2))) then
       temp = acforcing_surface(4,num_edge)
       acforcing_surface(4,num_edge) = acforcing_surface(3,num_edge)
       acforcing_surface(3,num_edge) = temp
       acforcing_surface_char(IEDGE2,match) = .true.
    endif

    if ((acforcing_surface(3,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+2) .and. &
          acforcing_surface(4,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+3))) then
       temp = acforcing_surface(4,num_edge)
       acforcing_surface(4,num_edge) = acforcing_surface(3,num_edge)
       acforcing_surface(3,num_edge) = temp
       acforcing_surface_char(IEDGE3,match) = .true.
    endif

    if ((acforcing_surface(4,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+2) .and. &
          acforcing_surface(3,num_edge) == elmnts(ngnod*acforcing_surface_merge(match)+3))) then
       acforcing_surface_char(IEDGE3,match) = .true.
    endif

  enddo

  nelemacforcing_merge = nb_elmnts_acforcing

! beware here and below that external meshes (for instance coming from CUBIT or Gmsh)
! may have rotated elements and thus edge 1 may not correspond to the bottom,
! edge 2 may not correspond to the right, edge 3 may not correspond to the top,
! and edge 4 may not correspond to the left.
  allocate(ibegin_edge1_acforcing(nelemacforcing_merge))
  allocate(iend_edge1_acforcing(nelemacforcing_merge))
  allocate(ibegin_edge2_acforcing(nelemacforcing_merge))
  allocate(iend_edge2_acforcing(nelemacforcing_merge))
  allocate(ibegin_edge3_acforcing(nelemacforcing_merge))
  allocate(iend_edge3_acforcing(nelemacforcing_merge))
  allocate(ibegin_edge4_acforcing(nelemacforcing_merge))
  allocate(iend_edge4_acforcing(nelemacforcing_merge))

  ibegin_edge1_acforcing(:) = 1
  ibegin_edge2_acforcing(:) = 1
  ibegin_edge3_acforcing(:) = 1
  ibegin_edge4_acforcing(:) = 1
  iend_edge1_acforcing(:) = NGLLX
  iend_edge2_acforcing(:) = NGLLZ
  iend_edge3_acforcing(:) = NGLLX
  iend_edge4_acforcing(:) = NGLLZ

  end subroutine merge_acoustic_forcing_boundaries

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! Write abs surface (elements and nodes on the surface) pertaining to iproc partition in the corresponding Database
  ! inspired by write_abs_merge_database upper in this file
  !--------------------------------------------------

  subroutine write_acoustic_forcing_merge_database(IIN_database, iproc, num_phase)

  use part_unstruct_par, only: part,nelemacforcing_loc,nelemacforcing_merge, &
    acforcing_surface_merge,acforcing_surface_char,acforcing_surface, &
    glob2loc_elmnts, &
    ibegin_edge1_acforcing,iend_edge1_acforcing,ibegin_edge2_acforcing,iend_edge2_acforcing, &
    ibegin_edge3_acforcing,iend_edge3_acforcing,ibegin_edge4_acforcing,iend_edge4_acforcing

  implicit none

  integer, intent(in) :: IIN_database
  integer, intent(in) :: iproc
  integer, intent(in) :: num_phase

  ! local parameters
  integer  :: i

  if (num_phase == 1) then
    ! only counts elements in this partition
    nelemacforcing_loc = 0
    do i = 1, nelemacforcing_merge
       if (part(acforcing_surface_merge(i)) == iproc) then
          nelemacforcing_loc = nelemacforcing_loc + 1
       endif
    enddo
  else
    do i = 1, nelemacforcing_merge
       if (part(acforcing_surface_merge(i)) == iproc) then

! beware here and below that external meshes (for instance coming from CUBIT or Gmsh)
! may have rotated elements and thus edge 1 may not correspond to the bottom,
! edge 2 may not correspond to the right, edge 3 may not correspond to the top,
! and edge 4 may not correspond to the left.

! we add +1 to the "glob2loc_elmnts" number because internally in parts of our code numbers start at zero instead of one
          write(IIN_database) glob2loc_elmnts(acforcing_surface_merge(i))+1, acforcing_surface_char(1,i), &

               acforcing_surface_char(2,i), acforcing_surface_char(3,i), acforcing_surface_char(4,i), acforcing_surface(5,i), &
               ibegin_edge1_acforcing(i), iend_edge1_acforcing(i), &
               ibegin_edge2_acforcing(i), iend_edge2_acforcing(i), &
               ibegin_edge3_acforcing(i), iend_edge3_acforcing(i), &
               ibegin_edge4_acforcing(i), iend_edge4_acforcing(i)
       endif
    enddo
  endif

  end subroutine write_acoustic_forcing_merge_database

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! Write fluid/solid edges (fluid (or porous) elements and corresponding solid (or porous) elements)
  ! pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------

  subroutine write_fluidsolid_edges_database(IIN_database, nedges_coupled_bis, edges_coupled_bis, &
                                             nedges_coupled_loc_bis, iproc, num_phase)

  use part_unstruct_par, only: part,glob2loc_elmnts

  implicit none

  integer, intent(in) :: IIN_database

  integer, intent(in) :: nedges_coupled_bis
  integer, dimension(2,nedges_coupled_bis) :: edges_coupled_bis

  integer, intent(inout)  :: nedges_coupled_loc_bis

  integer, intent(in) :: iproc
  integer, intent(in) :: num_phase

  ! local parameters
  integer  :: i

  if (num_phase == 1) then
    ! only counts elements in this partition
    nedges_coupled_loc_bis = 0
    do i = 1, nedges_coupled_bis
      if (part(edges_coupled_bis(1,i)) == iproc) then
        nedges_coupled_loc_bis = nedges_coupled_loc_bis + 1
      endif
    enddo
  else
    do i = 1, nedges_coupled_bis
      if (part(edges_coupled_bis(1,i)) == iproc) then
        write(IIN_database) glob2loc_elmnts(edges_coupled_bis(1,i))+1, glob2loc_elmnts(edges_coupled_bis(2,i))+1
      endif
    enddo
  endif

  end subroutine write_fluidsolid_edges_database

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! Write ispec of axial elements pertaining to iproc partition in the corresponding Database
  !--------------------------------------------------

  subroutine write_axial_elements_database(IIN_database, nelem_on_the_axis, ispec_of_axial_elements, &
                                           nelem_on_the_axis_loc, iproc, num_phase, remove_min_to_start_at_zero)

  use part_unstruct_par, only: part,glob2loc_elmnts

  implicit none

  integer, intent(in)  :: IIN_database

  integer, intent(in)  :: nelem_on_the_axis
  integer, dimension(nelem_on_the_axis)  :: ispec_of_axial_elements

  integer, intent(inout)  :: nelem_on_the_axis_loc

  integer, intent(in)  :: iproc
  integer, intent(in)  :: num_phase

  integer, intent(in)  :: remove_min_to_start_at_zero

  ! local parameters
  integer  :: i,ispec

  if (num_phase == 1) then
    ! only counts elements in this partition
    nelem_on_the_axis_loc = 0
    do i = 1, nelem_on_the_axis
      if (part(ispec_of_axial_elements(i)) == iproc) then
          nelem_on_the_axis_loc = nelem_on_the_axis_loc + 1
      endif
    enddo
  else
    do i = 1, nelem_on_the_axis
      ! if (part(ispec_of_axial_elements(i)) == 0 .and. iproc == 1) then
      !  print *,"ispec_of_axial_elements :",ispec_of_axial_elements(i)," -----> glob2loc_elmnts :", &
      !     glob2loc_elmnts(ispec_of_axial_elements(i))
      ! endif

      if (part(ispec_of_axial_elements(i)) == iproc) then
        ispec = glob2loc_elmnts(ispec_of_axial_elements(i)) + remove_min_to_start_at_zero
        write(IIN_database) ispec
      endif
    enddo
  endif

  end subroutine write_axial_elements_database

