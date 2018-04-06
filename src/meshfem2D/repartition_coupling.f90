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

  !--------------------------------------------------
  ! repartitioning: coupled acoustic/elastic elements are transferred to the same partition
  !--------------------------------------------------

  subroutine acoustic_elastic_repartitioning(elmnts_l, nbmodels, phi_material, num_material, nproc)

  use constants, only: IMAIN,NCORNERS,MAX_NEIGHBORS,TINYVAL
  use part_unstruct_par, only: nelmnts,edges_coupled,nedges_coupled

  implicit none

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, intent(in)  :: nproc, nbmodels
  double precision, dimension(nbmodels), intent(in)  :: phi_material
  integer, dimension(1:nelmnts), intent(in)  :: num_material

  ! local parameters
  integer, dimension(:), allocatable  :: xadj_l
  integer, dimension(:), allocatable  :: adjncy_l
  logical, dimension(nbmodels)  :: is_acoustic, is_elastic
  integer  :: i, ier
  integer  :: el, el_adj

  allocate(xadj_l(0:nelmnts))
  allocate(adjncy_l(0:MAX_NEIGHBORS*nelmnts-1))

  ! sets domain flags
  is_acoustic(:) = .false.
  is_elastic(:) = .false.

  do i = 1, nbmodels
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
    ! for acoustic element
    if (is_acoustic(num_material(el+1))) then
      ! loops over adjacent elements
      do el_adj = xadj_l(el), xadj_l(el+1) - 1
        ! adds its elastic neighbor
        if (is_elastic(num_material(adjncy_l(el_adj)+1))) then
          nedges_coupled = nedges_coupled + 1
        endif
      enddo
    endif
  enddo

  ! user output
  write(IMAIN,*) 'nedges_coupled (acoustic/elastic)     = ', nedges_coupled

  allocate(edges_coupled(2,nedges_coupled),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating array edges_coupled')
  edges_coupled(:,:) = 0

  ! repartitions elements
  if (nedges_coupled > 0) then
    call repartition_coupled_edges(nproc,nedges_coupled,edges_coupled, &
                                   num_material,nbmodels, &
                                   is_acoustic,is_elastic,xadj_l,adjncy_l)
  endif

  deallocate(xadj_l,adjncy_l)

  end subroutine acoustic_elastic_repartitioning


!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! repartitioning: coupled acoustic/poroelastic elements are transferred to the same partition
  !--------------------------------------------------

  subroutine acoustic_poro_repartitioning(elmnts_l, nbmodels, phi_material, num_material, nproc)

  use constants, only: IMAIN,NCORNERS,MAX_NEIGHBORS,TINYVAL
  use part_unstruct_par, only: nelmnts,edges_acporo_coupled,nedges_acporo_coupled

  implicit none

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, intent(in)  :: nproc, nbmodels
  double precision, dimension(nbmodels), intent(in)  :: phi_material
  integer, dimension(1:nelmnts), intent(in)  :: num_material

  ! local parameters
  integer, dimension(:), allocatable  :: xadj_l
  integer, dimension(:), allocatable  :: adjncy_l
  logical, dimension(nbmodels)  :: is_acoustic,is_poroelastic
  integer  :: i, ier
  integer  :: el, el_adj

  allocate(xadj_l(0:nelmnts))
  allocate(adjncy_l(0:MAX_NEIGHBORS*nelmnts-1))

  is_acoustic(:) = .false.
  is_poroelastic(:) = .false.

  do i = 1, nbmodels
     if (phi_material(i) >= 1.d0) then
        is_acoustic(i) = .true.
     endif
     if (phi_material(i) < 1.d0 .and. phi_material(i) > TINYVAL) then
        is_poroelastic(i) = .true.
     endif
  enddo

  ! determines maximum neighbors based on 2 common nodes (common edge)
  call mesh2dual_ncommonnodes(elmnts_l, 2, xadj_l, adjncy_l)

  nedges_acporo_coupled = 0
  do el = 0, nelmnts-1
     if (is_acoustic(num_material(el+1))) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if (is_poroelastic(num_material(adjncy_l(el_adj)+1))) then
              nedges_acporo_coupled = nedges_acporo_coupled + 1
           endif
        enddo
     endif
  enddo

  ! user output
  write(IMAIN,*) 'nedges_coupled (acoustic/poroelastic) = ', nedges_acporo_coupled

  allocate(edges_acporo_coupled(2,nedges_acporo_coupled),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating array edges_acporo_coupled')
  edges_acporo_coupled(:,:) = 0

  ! repartitions elements
  if (nedges_acporo_coupled > 0) then
    call repartition_coupled_edges(nproc,nedges_acporo_coupled,edges_acporo_coupled, &
                                   num_material,nbmodels, &
                                   is_acoustic,is_poroelastic,xadj_l,adjncy_l)
  endif

  deallocate(xadj_l,adjncy_l)

  end subroutine acoustic_poro_repartitioning

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! repartitioning: coupled poroelastic/elastic elements are transferred to the same partition
  !--------------------------------------------------

  subroutine poro_elastic_repartitioning(elmnts_l, nbmodels, phi_material, num_material, nproc)

  use constants, only: IMAIN,NCORNERS,MAX_NEIGHBORS,TINYVAL
  use part_unstruct_par, only: nelmnts,nedges_elporo_coupled,edges_elporo_coupled

  implicit none

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in)  :: elmnts_l
  integer, intent(in)  :: nproc, nbmodels
  double precision, dimension(nbmodels), intent(in)  :: phi_material
  integer, dimension(1:nelmnts), intent(in)  :: num_material

  ! local parameters
  integer, dimension(:), allocatable  :: xadj_l
  integer, dimension(:), allocatable  :: adjncy_l
  logical, dimension(nbmodels)  :: is_elastic,is_poroelastic
  integer  :: i, ier
  integer  :: el, el_adj

  allocate(xadj_l(0:nelmnts))
  allocate(adjncy_l(0:MAX_NEIGHBORS*nelmnts-1))

  is_elastic(:) = .false.
  is_poroelastic(:) = .false.

  do i = 1, nbmodels
     if (phi_material(i) < TINYVAL) then
        is_elastic(i) = .true.
     endif
     if (phi_material(i) < 1.d0 .and. phi_material(i) > TINYVAL) then
        is_poroelastic(i) = .true.
     endif
  enddo

  ! determines maximum neighbors based on 2 common nodes (common edge)
  call mesh2dual_ncommonnodes(elmnts_l, 2, xadj_l, adjncy_l)

  nedges_elporo_coupled = 0
  do el = 0, nelmnts-1
     if (is_poroelastic(num_material(el+1))) then
        do el_adj = xadj_l(el), xadj_l(el+1) - 1
           if (is_elastic(num_material(adjncy_l(el_adj)+1))) then
              nedges_elporo_coupled = nedges_elporo_coupled + 1
           endif
        enddo
     endif
  enddo

  ! user output
  write(IMAIN,*) 'nedges_coupled (poroelastic/elastic)  = ', nedges_elporo_coupled

  allocate(edges_elporo_coupled(2,nedges_elporo_coupled),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating array edges_elporo_coupled')
  edges_elporo_coupled(:,:) = 0

  ! repartitions elements
  if (nedges_elporo_coupled > 0) then
    call repartition_coupled_edges(nproc,nedges_elporo_coupled,edges_elporo_coupled, &
                                   num_material,nbmodels, &
                                   is_poroelastic,is_elastic,xadj_l,adjncy_l)
  endif

  deallocate(xadj_l,adjncy_l)

  end subroutine poro_elastic_repartitioning

!
!---------------------------------------------------------------------------------------
!

  subroutine repartition_coupled_edges(nproc,nedges_coupled,edges_coupled, &
                                       num_material,nbmodels, &
                                       is_domain_A,is_domain_B,xadj_l,adjncy_l)

  use constants, only: IMAIN,MAX_NEIGHBORS
  use part_unstruct_par, only: nelmnts,part

  implicit none

  integer, intent(in)  :: nproc, nedges_coupled
  integer, dimension(2,nedges_coupled),intent(inout) :: edges_coupled

  integer, dimension(1:nelmnts), intent(in)  :: num_material

  integer,intent(in) :: nbmodels
  logical, dimension(nbmodels),intent(in) :: is_domain_A, is_domain_B

  integer, dimension(0:nelmnts),intent(in) :: xadj_l
  integer, dimension(0:MAX_NEIGHBORS*nelmnts-1),intent(in) :: adjncy_l

  ! local parameters
  integer  :: i, iedge
  integer  :: el, el_adj
  logical  :: is_repartitioned

  ! sets edges
  iedge = 0
  do el = 0, nelmnts-1
    ! reference element
    if (is_domain_A(num_material(el+1))) then
      ! loops over adjacent elements
      do el_adj = xadj_l(el), xadj_l(el+1) - 1
        ! adds coupled edge
        if (is_domain_B(num_material(adjncy_l(el_adj)+1))) then
          iedge = iedge + 1
          edges_coupled(1,iedge) = el
          edges_coupled(2,iedge) = adjncy_l(el_adj)
        endif
      enddo
    endif
  enddo
  if (iedge /= nedges_coupled) call stop_the_code('Error in setting domain edges, number of edges invalid')

  ! only in case we have different partitions
  if (nproc > 1) then
    do i = 1, nedges_coupled * nproc
      is_repartitioned = .false.
      do iedge = 1, nedges_coupled
        ! puts coupled element in same partition
        if (part(edges_coupled(1,iedge)) /= part(edges_coupled(2,iedge))) then
          ! moves element into partition with smaller process id
          if (part(edges_coupled(1,iedge)) < part(edges_coupled(2,iedge))) then
            part(edges_coupled(2,iedge)) = part(edges_coupled(1,iedge))
          else
            part(edges_coupled(1,iedge)) = part(edges_coupled(2,iedge))
          endif
          is_repartitioned = .true.
        endif
      enddo
      ! check if there is still work to do
      if (.not. is_repartitioned) then
        exit
      endif
    enddo

    ! checks if initial coupled edges are repartitioned
    if (is_repartitioned) then
      ! checks count in case we need more
      i = 0
      do iedge = 1, nedges_coupled
        if (part(edges_coupled(1,iedge)) /= part(edges_coupled(2,iedge))) i = i + 1
      enddo
      write(IMAIN,*) '  repartitioning edges left = ',i
      if (i /= 0) then
        write(IMAIN,*) 'Error: repartitioning edges has still edges left = ',i
        call stop_the_code('Error: repartitioning coupled elements needs more iterations')
      else
        ! for user output
        i = nedges_coupled * nproc
      endif
    endif
    write(IMAIN,*) '  after iteration ',i,'repartitioning of all coupled elements done'
  endif

  end subroutine repartition_coupled_edges

!
!---------------------------------------------------------------------------------------
!

  !--------------------------------------------------
  ! repartitioning: coupled periodic elements are transferred to the same partition
  !--------------------------------------------------

  subroutine periodic_edges_repartitioning(elmnts_l,nnodes,nodes_coords,PERIODIC_HORIZ_DIST)

  use constants, only: IMAIN,NCORNERS
  use part_unstruct_par, only: nelmnts,part

  implicit none

  integer, dimension(0:NCORNERS*nelmnts-1), intent(in) :: elmnts_l

  integer :: nnodes
  double precision, dimension(2,nnodes) :: nodes_coords
  double precision :: PERIODIC_HORIZ_DIST

  ! local parameters
  logical, dimension(0:nelmnts-1) :: is_periodic

  integer :: el,el2,icorner,icorner2,num_node,num_node2,ifirst_partition_found

  double precision :: xtol,xtypdist
  double precision :: x,y,x2,y2

! set up a local geometric tolerance by computing the typical horizontal size of an element.
! the sqrt() assumes that the geometrical model is 'not too elongated' and thus 'not too far from a square'
! and thus contains more or less the same number of points along X and Y. If this is not the case i.e.
! if the model is very elongated then this trick will work anyway because we just want to have a rough idea
! of a typical length in the mesh, even if it is not very accurate it will work anyway.
  xtypdist = (maxval(nodes_coords(1,:)) - minval(nodes_coords(1,:))) / sqrt(dble(nnodes))

! define a tolerance, small with respect to the minimum size
  xtol = 1.d-4 * xtypdist

! detect the points that are on the same horizontal line (i.e. at the same height Z)
! and that have a value of the horizontal coordinate X that differs by exactly the periodicity length;
! if so, make them all have the same global number, which will then implement periodic boundary conditions automatically.
! We select the smallest value of iglob and assign it to all the points that are the same due to periodicity,
! this way the maximum value of the ibool() array will remain as small as possible.
!
! *** IMPORTANT: this simple algorithm will be slow for large meshes because it has a cost of NGLOB^2 / 2
! (where NGLOB is the number of points per MPI slice, not of the whole mesh though). This could be
! reduced to O(NGLOB log(NGLOB)) by using a quicksort algorithm on the coordinates of the points to detect the multiples
! (as implemented in routine createnum_fast() elsewhere in the code). This could be done one day if needed instead
! of the very simple double loop below.

  ! user output
  write(IMAIN,*) 'start detecting points for periodic boundary conditions &
                 &(the current algorithm can be slow and could be improved)...'

  is_periodic(:) = .false.

! loop on all the elements
  do el = 0, nelmnts-2 ! we call stop_the_code(one element before the end in order for the second loop to be OK in all cases)
    do el2 = el+1, nelmnts-1
      if (is_periodic(el2)) cycle
      ! it is sufficient to loop on the four corners to determine if this element has at least one periodic point
      do icorner = 0,NCORNERS-1
        num_node = elmnts_l(icorner + NCORNERS*el) + 1 ! the plus one is because elmnts_l() starts at zero
        x = nodes_coords(1,num_node)
        y = nodes_coords(2,num_node)
        do icorner2 = 0,NCORNERS-1
          num_node2 = elmnts_l(icorner2 + NCORNERS*el2) + 1 ! the plus one is because elmnts_l() starts at zero
          x2 = nodes_coords(1,num_node2)
          y2 = nodes_coords(2,num_node2)
          ! if the two points are at the same height Y
          if (abs(y2 - y) < xtol) then
            ! if in addition their X coordinates differ by exactly the periodicity distance
            if (abs(abs(x2 - x) - PERIODIC_HORIZ_DIST) < xtol) then
              ! then these two elements are in contact by a periodic edge
              is_periodic(el) = .true.
              is_periodic(el2) = .true.
              goto 100
            endif
          endif
        enddo
      enddo
 100  continue
    enddo
  enddo

  ! user output
  write(IMAIN,*) 'done detecting points for periodic boundary conditions.'
  write(IMAIN,*) 'number of periodic elements found and grouped in the same partition: ',count(is_periodic)
  call flush_IMAIN()

  ! loop on all the elements to find the first partition that contains a periodic element
  ifirst_partition_found = -1
  do el = 0, nelmnts-1
    if (is_periodic(el)) then
      ifirst_partition_found = part(el)
      exit
    endif
  enddo
  if (ifirst_partition_found < 0) call stop_the_code( &
'error: no periodic element found, even though ADD_PERIODIC_CONDITIONS is set')

  ! loop on all the elements to move all periodic elements to the first partition found
  do el = 0, nelmnts-1
    if (is_periodic(el)) part(el) = ifirst_partition_found
  enddo

  end subroutine periodic_edges_repartitioning

!
!---------------------------------------------------------------------------------------
!

  subroutine manual_crack_repartitioning(num_material,NPROC)

! puts elements along manual crack into the same partition
!
! note: For now, elements material numbers are hard-coded and must be 2 (for left side) and 3 (right side)
!       to indicate an element along the crack.
!       To split nodes, see routine in add_manual_crack.f90

  use constants, only: IMAIN,NCORNERS,ADD_A_SMALL_CRACK_IN_THE_MEDIUM
  use part_unstruct_par, only: nelmnts,part

  implicit none

  integer, intent(in) :: nproc
  integer, dimension(1:nelmnts), intent(in)  :: num_material

  ! local parameters
  logical, dimension(0:nelmnts-1) :: is_crack_element
  integer  :: ipartition_crack
  integer  :: el

  ! checks if anything to do
  if (.not. ADD_A_SMALL_CRACK_IN_THE_MEDIUM) return

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'detecting elements for manual crack.'

  ! sets flags for elements along the crack (material number 2 and 3)
  is_crack_element(:) = .false.
  ipartition_crack = nproc + 1
  do el = 0, nelmnts-1
    ! crack between material number 2 and 3
    if (num_material(el+1) == 2 .or. num_material(el+1) == 3) then
      is_crack_element(el) = .true.
      ! puts all crack elements into lowest partition possible
      if (part(el) < ipartition_crack) ipartition_crack = part(el)
    endif
  enddo

  ! user output
  write(IMAIN,*) 'number of crack elements ',count(is_crack_element),' found and grouped in the same partition ',ipartition_crack
  call flush_IMAIN()

  if (count(is_crack_element) == 0) &
    call stop_the_code('Error: no crack element found, even though ADD_A_SMALL_CRACK_IN_THE_MEDIUM is set')
  if (ipartition_crack > nproc) &
    call stop_the_code('Error: invalid partition number for crack elements')

  ! we will put all crack elements into the same partition
  do el = 0, nelmnts-1
    if (is_crack_element(el)) part(el) = ipartition_crack
  enddo

  end subroutine manual_crack_repartitioning

