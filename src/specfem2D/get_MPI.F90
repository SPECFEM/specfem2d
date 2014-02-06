
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, Inria and University of Pau, France,
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
! license as circulated by CEA, CNRS and Inria at the following URL
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

#ifdef USE_MPI

  subroutine get_MPI(nspec,ibool,knods,ngnod,nglob,elastic,poroelastic, &
                    ninterface, max_interface_size, &
                    my_nelmnts_neighbours,my_interfaces,my_neighbours, &
                    ibool_interfaces_acoustic, ibool_interfaces_elastic, &
                    ibool_interfaces_poroelastic, &
                    nibool_interfaces_acoustic, nibool_interfaces_elastic, &
                    nibool_interfaces_poroelastic, &
                    inum_interfaces_acoustic, inum_interfaces_elastic, &
                    inum_interfaces_poroelastic, &
                    ninterface_acoustic, ninterface_elastic, ninterface_poroelastic, &
                    mask_ispec_inner_outer, &
                    myrank,ipass,coord)

! sets up the MPI interface for communication between partitions

  use :: mpi

  implicit none

  include "constants.h"

  integer, intent(in)  :: nspec, nglob, ngnod
  logical, dimension(nspec), intent(in)  :: elastic, poroelastic
  integer, dimension(ngnod,nspec), intent(in)  :: knods
  integer, dimension(NGLLX,NGLLZ,nspec), intent(in)  :: ibool

  integer  :: ninterface
  integer  :: max_interface_size
  integer, dimension(ninterface)  :: my_nelmnts_neighbours,my_neighbours
  integer, dimension(4,max_interface_size,ninterface)  :: my_interfaces

  integer, dimension(NGLLX*max_interface_size,ninterface)  :: &
       ibool_interfaces_acoustic,ibool_interfaces_elastic,ibool_interfaces_poroelastic
  integer, dimension(ninterface)  :: &
       nibool_interfaces_acoustic,nibool_interfaces_elastic,nibool_interfaces_poroelastic
  integer, dimension(ninterface), intent(out)  :: &
       inum_interfaces_acoustic, inum_interfaces_elastic, inum_interfaces_poroelastic
  integer, intent(out)  :: ninterface_acoustic, ninterface_elastic, ninterface_poroelastic

  logical, dimension(nspec), intent(inout)  :: mask_ispec_inner_outer

  integer :: myrank,ipass
  double precision, dimension(NDIM,nglob) :: coord

  !local parameters
  double precision, dimension(:), allocatable :: xp,zp
  double precision, dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: locval
  integer, dimension(:), allocatable :: nibool_interfaces_true
  ! for MPI buffers
  integer, dimension(:), allocatable :: reorder_interface,ind,ninseg,iwork
  integer, dimension(:), allocatable :: ibool_dummy
!  integer, dimension(:,:), allocatable :: ibool_interfaces_dummy
  logical, dimension(:), allocatable :: ifseg
  integer :: iinterface,ilocnum
  integer :: num_points1, num_points2
  ! assembly test
  integer :: i,j,ispec,iglob,count,inum,ier,idomain
  integer :: max_nibool_interfaces,num_nibool,num_interface
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: test_flag_cr
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_send_faces_vector_ac
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_recv_faces_vector_ac
  integer, dimension(:), allocatable  :: tab_requests_send_recv_acoustic

  ! gets global indices for points on MPI interfaces
  ! (defined by my_interfaces) between different partitions
  ! and stores them in ibool_interfaces*** & nibool_interfaces*** (number of total points)
  call prepare_assemble_MPI(nspec,ibool,knods, ngnod,nglob, elastic, poroelastic, &
                                ninterface, max_interface_size, &
                                my_nelmnts_neighbours, my_interfaces, &
                                ibool_interfaces_acoustic, ibool_interfaces_elastic, &
                                ibool_interfaces_poroelastic, &
                                nibool_interfaces_acoustic, nibool_interfaces_elastic, &
                                nibool_interfaces_poroelastic, &
                                inum_interfaces_acoustic, inum_interfaces_elastic, &
                                inum_interfaces_poroelastic, &
                                ninterface_acoustic, ninterface_elastic, ninterface_poroelastic, &
                                mask_ispec_inner_outer )


  ! sorts ibool comm buffers lexicographically for all MPI interfaces
  num_points1 = 0
  num_points2 = 0
  allocate(nibool_interfaces_true(ninterface))

  do idomain = 1,3

    ! checks number of interface in this domain
    num_interface = 0
    if( idomain == 1 ) then
      num_interface = ninterface_acoustic
    else if( idomain == 2 ) then
      num_interface = ninterface_elastic
    else if( idomain == 3 ) then
      num_interface = ninterface_poroelastic
    endif
    if( num_interface == 0 ) cycle

    ! loops over interfaces
    do iinterface = 1, ninterface

      ! number of global points in this interface
      num_nibool = 0
      if( idomain == 1 ) then
        num_nibool = nibool_interfaces_acoustic(iinterface)
      else if( idomain == 2 ) then
        num_nibool = nibool_interfaces_elastic(iinterface)
      else if( idomain == 3 ) then
        num_nibool = nibool_interfaces_poroelastic(iinterface)
      endif
      ! checks if anything to sort
      if( num_nibool == 0 ) cycle

      allocate(xp(num_nibool))
      allocate(zp(num_nibool))
      allocate(locval(num_nibool))
      allocate(ifseg(num_nibool))
      allocate(reorder_interface(num_nibool))
      allocate(ibool_dummy(num_nibool))
      allocate(ind(num_nibool))
      allocate(ninseg(num_nibool))
      allocate(iwork(num_nibool))
      allocate(work(num_nibool))

      ! works with a copy of ibool array
      if( idomain == 1 ) then
        ibool_dummy(:) = ibool_interfaces_acoustic(1:num_nibool,iinterface)
      else if( idomain == 2 ) then
        ibool_dummy(:) = ibool_interfaces_elastic(1:num_nibool,iinterface)
      else if( idomain == 3 ) then
        ibool_dummy(:) = ibool_interfaces_poroelastic(1:num_nibool,iinterface)
      endif

      ! gets x,y,z coordinates of global points on MPI interface
      do ilocnum = 1, num_nibool
        iglob = ibool_dummy(ilocnum)
        xp(ilocnum) = coord(1,iglob)
        zp(ilocnum) = coord(2,iglob)
      enddo

      ! sorts (lexicographically?) ibool_interfaces and updates value
      ! of total number of points nibool_interfaces_true(iinterface)
      call sort_array_coordinates(num_nibool,xp,zp, &
                                ibool_dummy, &
                                reorder_interface,locval,ifseg, &
                                nibool_interfaces_true(iinterface), &
                                ind,ninseg,iwork,work)

      ! checks that number of MPI points are still the same
      num_points1 = num_points1 + num_nibool
      num_points2 = num_points2 + nibool_interfaces_true(iinterface)
      if( num_points1 /= num_points2 ) then
        write(IOUT,*) 'error sorting MPI interface points:',myrank
        write(IOUT,*) '   domain:',idomain
        write(IOUT,*) '   interface:',iinterface,num_points1,num_points2
        call exit_MPI('error sorting MPI interface')
      endif

      ! stores new order of ibool array
      if( idomain == 1 ) then
        ibool_interfaces_acoustic(1:num_nibool,iinterface) = ibool_dummy(:)
      else if( idomain == 2 ) then
        ibool_interfaces_elastic(1:num_nibool,iinterface) = ibool_dummy(:)
      else if( idomain == 3 ) then
        ibool_interfaces_poroelastic(1:num_nibool,iinterface) = ibool_dummy(:)
      endif

      ! cleanup temporary arrays
      deallocate(xp)
      deallocate(zp)
      deallocate(locval)
      deallocate(ifseg)
      deallocate(reorder_interface)
      deallocate(ibool_dummy)
      deallocate(ind)
      deallocate(ninseg)
      deallocate(iwork)
      deallocate(work)
    enddo
  enddo

  ! cleanup
  deallocate(nibool_interfaces_true)

  ! outputs total number of MPI interface points
  call MPI_REDUCE(num_points2, num_points1, 1, MPI_INTEGER, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ier)
  if( myrank == 0 .and. ipass == 1 ) then
    write(IOUT,*) 'total MPI interface points: ',num_points1
  endif

  ! checks interfaces in acoustic domains
  inum = 0
  count = 0
  if ( ninterface_acoustic > 0) then

    ! checks with assembly of test fields
    allocate(test_flag_cr(nglob))
    test_flag_cr(:) = 0._CUSTOM_REAL
    count = 0
    do ispec = 1, nspec
      ! sets flags on global points
      do j = 1, NGLLZ
        do i = 1, NGLLX
          ! global index
          iglob = ibool(i,j,ispec)

          ! counts number of unique global points to set
          if( nint(test_flag_cr(iglob)) == 0 ) count = count+1

          ! sets identifier
          test_flag_cr(iglob) = myrank + 1.0
        enddo
      enddo
    enddo

    max_nibool_interfaces = maxval(nibool_interfaces_acoustic(:))

    allocate(tab_requests_send_recv_acoustic(ninterface_acoustic*2))
    allocate(buffer_send_faces_vector_ac(max_nibool_interfaces,ninterface_acoustic))
    allocate(buffer_recv_faces_vector_ac(max_nibool_interfaces,ninterface_acoustic))
    inum = 0
    do iinterface = 1, ninterface
      inum = inum + nibool_interfaces_acoustic(iinterface)
    enddo
  endif

  ! note: this mpi reduction awaits information from all processes.
  !          thus, avoid an mpi deadlock in case some of the paritions have no acoustic interface
  call MPI_REDUCE(inum, num_points1, 1, MPI_INTEGER, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ier)

  if( myrank == 0 .and. ipass == 1 ) then
    write(IOUT,*) '       acoustic interface points: ',num_points1
  endif

  ! checks if assembly works
  inum = 0
  if( ninterface_acoustic > 0 ) then
    ! adds contributions from different partitions to flag arrays
    ! custom_real arrays
    call assemble_MPI_vector_ac(test_flag_cr,nglob, &
                    ninterface, ninterface_acoustic,inum_interfaces_acoustic, &
                    max_interface_size, max_nibool_interfaces,&
                    ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
                    tab_requests_send_recv_acoustic,buffer_send_faces_vector_ac, &
                    buffer_recv_faces_vector_ac, my_neighbours)

    ! checks number of interface points
    inum = 0
    do iglob=1,nglob
      ! only counts flags with MPI contributions
      if( nint(test_flag_cr(iglob)) > myrank+1 ) inum = inum + 1
    enddo

    deallocate(tab_requests_send_recv_acoustic)
    deallocate(buffer_send_faces_vector_ac)
    deallocate(buffer_recv_faces_vector_ac)
    deallocate(test_flag_cr)
  endif

  ! note: this mpi reduction awaits information from all processes.
  call MPI_REDUCE(inum, num_points2, 1, MPI_INTEGER, &
                    MPI_SUM, 0, MPI_COMM_WORLD, ier)

  if( myrank == 0 ) then
    if( ipass == 1 ) then
      write(IOUT,*) '       assembly acoustic MPI interface points:',num_points2
    endif

    ! they don't need to fit, somehow..
    !if( num_points2 /= num_points1 ) then
    !  print*,'error acoustic assembly:' !,myrank
    !  print*,'  total = ',num_points1,' not equal to assembled ',num_points2
    !  call exit_MPI('error acoustic MPI assembly')
    !endif
  endif

  end subroutine get_MPI

#endif
