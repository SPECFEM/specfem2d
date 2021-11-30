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


  subroutine decompose_mesh()

  use constants, only: IMAIN,MAX_NEIGHBORS,NCORNERS,MAX_NSIZE_SHARED,ADD_A_SMALL_CRACK_IN_THE_MEDIUM,myrank

  use shared_parameters, only: NPROC,ADD_PERIODIC_CONDITIONS,PERIODIC_HORIZ_DIST, &
    NGNOD,nbmodels,num_material,PARTITIONING_TYPE,phi_read

  use part_unstruct_par, only: part,nelmnts,xadj_g,adjncy_g,elmnts,elmnts_bis,nb_edges, &
    nnodes_elmnts,nnodes,nodes_elmnts,nodes_coords,ninterfaces

  use decompose_par
  use compute_elements_load_par

  implicit none

  ! local parameters
  integer :: i, iproc, ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  decomposing mesh using NPROC = ',NPROC
    call flush_IMAIN()
  endif

  ! allocates and initializes partitioning of elements
  allocate(part(0:nelmnts-1),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating partition array')
  part(:) = -1

  ! connectivity
  if (NPROC > 1) then
    allocate(xadj_g(0:nelmnts), &
             adjncy_g(0:MAX_NEIGHBORS*nelmnts-1),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating connectivity arrays')
    xadj_g(:) = 0
    adjncy_g(:) = -1
  endif

  ! construction of the graph
  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Graph construction:'
    call flush_IMAIN()
  endif

  ! if NGNOD == 9, we work on a subarray of elements that represents the elements with four nodes (four corners) only
  ! because the adjacency of the mesh elements can be entirely determined from the knowledge of the four corners only
  if (NGNOD == 9) then
    allocate(elmnts_bis(0:NCORNERS*nelmnts-1),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating elmnts_bis array')

    do i = 0, nelmnts-1
      elmnts_bis(i*NCORNERS:i*NCORNERS+NCORNERS-1) = elmnts(i*NGNOD:i*NGNOD+NCORNERS-1)
    enddo

    if (NPROC > 1) then
!! DK DK fixed problem in the previous implementation by Nicolas Le Goff:
!! DK DK (nxread+1)*(nzread+1) is OK for a regular internal mesh only, not for non structured external meshes
!! DK DK      call mesh2dual_ncommonnodes(nelmnts, (nxread+1)*(nzread+1), &
!! DK DK                                    elmnts_bis, xadj, adjncy, nnodes_elmnts, nodes_elmnts,1)
!! DK DK the subset of element corners is not renumbered therefore we must still use the nnodes computed for 9 nodes here
      ! determines maximum neighbors based on 1 common node
      call mesh2dual_ncommonnodes(elmnts_bis,1,xadj_g,adjncy_g)
    endif
  else
    if (NPROC > 1) then
      ! determines maximum neighbors based on 1 common node
      call mesh2dual_ncommonnodes(elmnts,1,xadj_g,adjncy_g)
    endif
  endif
  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  graph adjacency done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! partitions
  if (NPROC == 1) then
    part(:) = 0 ! single process has rank 0
  else
    ! number of common edges
    nb_edges = xadj_g(nelmnts)

    ! compute elements load for efficient partitioning (a PML element take more time to calculate than a fluid element for ex)
    call compute_elements_load()

    ! partitioning
    select case (PARTITIONING_TYPE)
    case(1)
      ! analytical
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'Partitioning type: analytical'
        write(IMAIN,*)
        call flush_IMAIN()
      endif

      do iproc = 0, NPROC-2
        part(iproc*floor(real(nelmnts)/real(NPROC)):(iproc+1)*floor(real(nelmnts)/real(NPROC))-1) = iproc
      enddo
      part(floor(real(nelmnts)/real(NPROC))*(NPROC-1):nelmnts-1) = NPROC - 1

    case(2)
      ! METIS
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'Partitioning type: METIS'
        write(IMAIN,*)
        call flush_IMAIN()
      endif

      call metis_partitioning()

    case(3)
      ! SCOTCH
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'Partitioning type: SCOTCH'
        write(IMAIN,*)
        call flush_IMAIN()
      endif

      call scotch_partitioning()

    case default
      call stop_the_code('Error invalid partitioning type value! must be 1, 2 or 3, please check your Par_file...')
    end select

  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Coupled interfaces:'
    call flush_IMAIN()
  endif

  ! fluid-solid edges: coupled elements are transferred to the same partition
  if (NGNOD == 9) then
     call acoustic_elastic_repartitioning(elmnts_bis, nbmodels, phi_read, num_material, NPROC)
  else
     call acoustic_elastic_repartitioning(elmnts, nbmodels, phi_read, num_material, NPROC)
  endif

  ! fluid-porous edges: coupled elements are transferred to the same partition
  if (NGNOD == 9) then
     call acoustic_poro_repartitioning(elmnts_bis, nbmodels, phi_read, num_material, NPROC)
  else
     call acoustic_poro_repartitioning(elmnts, nbmodels, phi_read, num_material, NPROC)
  endif

  ! porous-solid edges: coupled elements are transferred to the same partition
  if (NGNOD == 9) then
     call poro_elastic_repartitioning(elmnts_bis, nbmodels, phi_read, num_material, NPROC)
  else
     call poro_elastic_repartitioning(elmnts, nbmodels, phi_read, num_material, NPROC)
  endif

  ! periodic edges: coupled elements are transferred to the same partition
  if (ADD_PERIODIC_CONDITIONS .and. NPROC > 1) then
    if (NGNOD == 9) then
       call periodic_edges_repartitioning(elmnts_bis,nnodes,nodes_coords,PERIODIC_HORIZ_DIST)
    else
       call periodic_edges_repartitioning(elmnts,nnodes,nodes_coords,PERIODIC_HORIZ_DIST)
    endif
  endif

  ! manual crack elements
  if (ADD_A_SMALL_CRACK_IN_THE_MEDIUM .and. NPROC > 1) then
    ! safety check
    if (NGNOD /= 4) then
      call stop_the_code('must currently have NGNOD == 4 when adding a crack manually')
    else
      call manual_crack_repartitioning(num_material,NPROC)
    endif
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Local numbering:'
    write(IMAIN,*) '  NPROC: ', NPROC
    write(IMAIN,*) '  number of elements: ',nelmnts
    call flush_IMAIN()
  endif

  ! local number of each element for each partition
  call construct_glob2loc_elmnts(NPROC)

  if (NGNOD == 9) then
    if (allocated(nnodes_elmnts) ) deallocate(nnodes_elmnts)
    if (allocated(nodes_elmnts) ) deallocate(nodes_elmnts)

    allocate(nnodes_elmnts(0:nnodes-1))
    allocate(nodes_elmnts(0:MAX_NSIZE_SHARED*nnodes-1))

    nnodes_elmnts(:) = 0
    nodes_elmnts(:) = 0
    do i = 0, NGNOD*nelmnts-1
      nodes_elmnts(elmnts(i)*MAX_NSIZE_SHARED + nnodes_elmnts(elmnts(i))) = i/NGNOD
      nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1
    enddo
  else
    if (NPROC < 2) then
      if (.not. allocated(nnodes_elmnts) ) allocate(nnodes_elmnts(0:nnodes-1))
      if (.not. allocated(nodes_elmnts) ) allocate(nodes_elmnts(0:MAX_NSIZE_SHARED*nnodes-1))

      nnodes_elmnts(:) = 0
      nodes_elmnts(:) = 0

      do i = 0, NGNOD*nelmnts-1
        nodes_elmnts(elmnts(i)*MAX_NSIZE_SHARED+nnodes_elmnts(elmnts(i))) = i/NGNOD
        nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1
      enddo
    endif
  endif

  ! local number of each node for each partition
  call Construct_glob2loc_nodes(NPROC)

  ! construct the interfaces between partitions (used for MPI assembly)
  if (NPROC > 1) then
     if (NGNOD == 9) then
        call Construct_interfaces(NPROC, elmnts_bis, &
                                  nbmodels, phi_read, num_material)
     else
        call Construct_interfaces(NPROC, elmnts, &
                                  nbmodels, phi_read, num_material)
     endif
     allocate(my_interfaces(0:ninterfaces-1))
     allocate(my_nb_interfaces(0:ninterfaces-1))
  else
     ! dummy allocation
     ninterfaces=0
     allocate(my_interfaces(0:ninterfaces-1))
     allocate(my_nb_interfaces(0:ninterfaces-1))
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'decompose mesh all done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine decompose_mesh


