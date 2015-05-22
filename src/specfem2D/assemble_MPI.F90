
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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

!
! This file contains subroutines related to assembling (of the mass matrix, potential_dot_dot and
! accel_elastic, accels_poroelastic, accelw_poroelastic).
! These subroutines are for the most part not used in the sequential version.
!


#ifdef USE_MPI

!-----------------------------------------------
! Assembling the mass matrix.
!-----------------------------------------------
  subroutine assemble_MPI_scalar(array_val1,npoin_val1, &
                              array_val2,array_val5,npoin_val2, &
                              array_val3,array_val4,npoin_val3)

  use mpi
  use specfem_par, only:ninterface,max_ibool_interfaces_size_ac, &
                              max_ibool_interfaces_size_el, &
                              max_ibool_interfaces_size_po, &
                              ibool_interfaces_acoustic,ibool_interfaces_elastic, &
                              ibool_interfaces_poroelastic, &
                              nibool_interfaces_acoustic,nibool_interfaces_elastic, &
                              nibool_interfaces_poroelastic,my_neighbours,ier


  implicit none

  include 'constants.h'

  ! array to assemble
  ! acoustic
  integer :: npoin_val1
  real(kind=CUSTOM_REAL), dimension(npoin_val1), intent(inout) :: array_val1
  ! elastic
  integer :: npoin_val2
  real(kind=CUSTOM_REAL), dimension(npoin_val2), intent(inout) :: array_val2,array_val5
  ! poroelastic
  integer :: npoin_val3
  real(kind=CUSTOM_REAL), dimension(npoin_val3), intent(inout) :: array_val3,array_val4

  integer  :: ipoin, num_interface,i
! there are now two different mass matrices for the elastic case
! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
  double precision, dimension(max_ibool_interfaces_size_ac+2*max_ibool_interfaces_size_el+&
       2*max_ibool_interfaces_size_po, ninterface)  :: &
       buffer_send_faces_scalar, &
       buffer_recv_faces_scalar
  integer, dimension(ninterface)  :: msg_requests

  buffer_send_faces_scalar(:,:) = 0.d0
  buffer_recv_faces_scalar(:,:) = 0.d0

  do num_interface = 1, ninterface

     ipoin = 0
     do i = 1, nibool_interfaces_acoustic(num_interface)
        ipoin = ipoin + 1
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val1(ibool_interfaces_acoustic(i,num_interface))
     enddo

     do i = 1, nibool_interfaces_elastic(num_interface)
        ipoin = ipoin + 1
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val2(ibool_interfaces_elastic(i,num_interface))
     enddo

     do i = 1, nibool_interfaces_elastic(num_interface)
        ipoin = ipoin + 1
! there are now two different mass matrices for the elastic case
! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val5(ibool_interfaces_elastic(i,num_interface))
     enddo

     do i = 1, nibool_interfaces_poroelastic(num_interface)
        ipoin = ipoin + 1
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val3(ibool_interfaces_poroelastic(i,num_interface))
     enddo
     do i = 1, nibool_interfaces_poroelastic(num_interface)
        ipoin = ipoin + 1
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val4(ibool_interfaces_poroelastic(i,num_interface))
     enddo

     ! non-blocking send
     call MPI_ISEND( buffer_send_faces_scalar(1,num_interface), &
! there are now two different mass matrices for the elastic case
! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
          nibool_interfaces_acoustic(num_interface)+2*nibool_interfaces_elastic(num_interface)+&
          2*nibool_interfaces_poroelastic(num_interface), &
          MPI_DOUBLE_PRECISION, &
          my_neighbours(num_interface), 11, &
          MPI_COMM_WORLD, msg_requests(num_interface), ier)

  enddo

  do num_interface = 1, ninterface

     ! starts a blocking receive
     call MPI_RECV ( buffer_recv_faces_scalar(1,num_interface), &
! there are now two different mass matrices for the elastic case
! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
          nibool_interfaces_acoustic(num_interface)+2*nibool_interfaces_elastic(num_interface)+&
          2*nibool_interfaces_poroelastic(num_interface), &
          MPI_DOUBLE_PRECISION, &
          my_neighbours(num_interface), 11, &
          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

     ipoin = 0
     do i = 1, nibool_interfaces_acoustic(num_interface)
        ipoin = ipoin + 1
        array_val1(ibool_interfaces_acoustic(i,num_interface)) = &
            array_val1(ibool_interfaces_acoustic(i,num_interface))  &
             + buffer_recv_faces_scalar(ipoin,num_interface)
     enddo

     do i = 1, nibool_interfaces_elastic(num_interface)
        ipoin = ipoin + 1
        array_val2(ibool_interfaces_elastic(i,num_interface)) = &
            array_val2(ibool_interfaces_elastic(i,num_interface))  &
            + buffer_recv_faces_scalar(ipoin,num_interface)
     enddo

     do i = 1, nibool_interfaces_elastic(num_interface)
        ipoin = ipoin + 1
! there are now two different mass matrices for the elastic case
! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
        array_val5(ibool_interfaces_elastic(i,num_interface)) = &
            array_val5(ibool_interfaces_elastic(i,num_interface))  &
            + buffer_recv_faces_scalar(ipoin,num_interface)
     enddo

     do i = 1, nibool_interfaces_poroelastic(num_interface)
        ipoin = ipoin + 1
        array_val3(ibool_interfaces_poroelastic(i,num_interface)) = &
            array_val3(ibool_interfaces_poroelastic(i,num_interface))  &
            + buffer_recv_faces_scalar(ipoin,num_interface)
     enddo
     do i = 1, nibool_interfaces_poroelastic(num_interface)
        ipoin = ipoin + 1
        array_val4(ibool_interfaces_poroelastic(i,num_interface)) = &
            array_val4(ibool_interfaces_poroelastic(i,num_interface)) &
            + buffer_recv_faces_scalar(ipoin,num_interface)
     enddo

  enddo

  end subroutine assemble_MPI_scalar


!-----------------------------------------------
! Assembling potential_dot_dot for acoustic elements :
! the buffers are filled, the ISEND and IRECV are started here, then
! contributions are added.
! The previous version included communication overlap using persistent
! communication, but the merging of the outer and inner elements rendered
! overlap no longer possible, while persistent communications were removed
! because trace tool MPITrace does not yet instrument those.
! Particular care should be taken concerning possible optimisations of the
! communication scheme.
!-----------------------------------------------
  subroutine assemble_MPI_vector_ac(array_val1)

  use mpi

  use specfem_par, only: nglob,ninterface_acoustic, &
                         inum_interfaces_acoustic, &
                         ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
                         tab_requests_send_recv_acoustic, &
                         buffer_send_faces_vector_ac, &
                         buffer_recv_faces_vector_ac, &
                         my_neighbours,ier

  implicit none

  include 'constants.h'
  include 'precision.h'

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: array_val1

  ! local parameters
  integer  :: ipoin, num_interface,iinterface, iglob

  ! initializes buffers
  buffer_send_faces_vector_ac(:,:) = 0._CUSTOM_REAL
  buffer_recv_faces_vector_ac(:,:) = 0._CUSTOM_REAL
  tab_requests_send_recv_acoustic(:) = 0

  ! loops over acoustic interfaces only
  do iinterface = 1, ninterface_acoustic

    ! gets interface index in the range of all interfaces [1,ninterface]
    num_interface = inum_interfaces_acoustic(iinterface)

    ! loops over all interface points
    do ipoin = 1, nibool_interfaces_acoustic(num_interface)
      iglob = ibool_interfaces_acoustic(ipoin,num_interface)

      ! copies array values to buffer
      buffer_send_faces_vector_ac(ipoin,iinterface) = array_val1(iglob)
    enddo

  enddo

  do iinterface = 1, ninterface_acoustic

    ! gets global interface index
    num_interface = inum_interfaces_acoustic(iinterface)


    ! non-blocking send
    call MPI_ISEND( buffer_send_faces_vector_ac(1,iinterface), &
             nibool_interfaces_acoustic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_acoustic(iinterface), ier)

    if ( ier /= MPI_SUCCESS ) then
      call exit_mpi('MPI_ISEND unsuccessful in assemble_MPI_vector_start')
    endif

    ! starts a non-blocking receive
    call MPI_IRECV ( buffer_recv_faces_vector_ac(1,iinterface), &
             nibool_interfaces_acoustic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_acoustic(ninterface_acoustic+iinterface), ier)

    if ( ier /= MPI_SUCCESS ) then
      call exit_mpi('MPI_IRECV unsuccessful in assemble_MPI_vector')
    endif

  enddo


  ! waits for MPI requests to complete (recv)
  ! each wait returns once the specified MPI request completed
  do iinterface = 1, ninterface_acoustic
    call MPI_Wait(tab_requests_send_recv_acoustic(ninterface_acoustic+iinterface), &
                  MPI_STATUS_IGNORE, ier)
  enddo

  ! assembles the array values
  do iinterface = 1, ninterface_acoustic

    ! gets global interface index
    num_interface = inum_interfaces_acoustic(iinterface)

    ! loops over all interface points
    do ipoin = 1, nibool_interfaces_acoustic(num_interface)
      iglob = ibool_interfaces_acoustic(ipoin,num_interface)
      ! adds buffer contribution
      array_val1(iglob) = array_val1(iglob) + buffer_recv_faces_vector_ac(ipoin,iinterface)
    enddo

  enddo


  ! waits for MPI requests to complete (send)
  ! just to make sure that all sending is done
  do iinterface = 1, ninterface_acoustic
    call MPI_Wait(tab_requests_send_recv_acoustic(iinterface), MPI_STATUS_IGNORE, ier)
  enddo


  end subroutine assemble_MPI_vector_ac


!-----------------------------------------------
! Assembling accel_elastic for elastic elements :
! the buffers are filled, the ISEND and IRECV are started here, then
! contributions are added.
! The previous version included communication overlap using persistent
! communication, but the merging of the outer and inner elements rendered
! overlap no longer possible, while persistent communications were removed
! because trace tool MPITrace does not yet instrument those.
! Particular care should be taken concerning possible optimisations of the
! communication scheme.
!-----------------------------------------------
  subroutine assemble_MPI_vector_el(array_val2)

  use mpi

  use specfem_par, only: nglob,ninterface_elastic, &
                         inum_interfaces_elastic, &
                         ibool_interfaces_elastic, nibool_interfaces_elastic, &
                         tab_requests_send_recv_elastic, &
                         buffer_send_faces_vector_el, &
                         buffer_recv_faces_vector_el, &
                         my_neighbours

  implicit none

  include 'constants.h'
  include 'precision.h'

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(3,nglob), intent(inout) :: array_val2


  integer  :: ipoin, num_interface, iinterface, ier, i

  do iinterface = 1, ninterface_elastic

     num_interface = inum_interfaces_elastic(iinterface)

     ipoin = 0
     do i = 1, nibool_interfaces_elastic(num_interface)
        buffer_send_faces_vector_el(ipoin+1:ipoin+3,iinterface) = &
             array_val2(:,ibool_interfaces_elastic(i,num_interface))
        ipoin = ipoin + 3
     enddo

  enddo

  do iinterface = 1, ninterface_elastic

    num_interface = inum_interfaces_elastic(iinterface)

    call MPI_ISEND( buffer_send_faces_vector_el(1,iinterface), &
             3*nibool_interfaces_elastic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_elastic(iinterface), ier)

    if ( ier /= MPI_SUCCESS ) then
      call exit_mpi('MPI_ISEND unsuccessful in assemble_MPI_vector_el')
    endif

    call MPI_IRECV ( buffer_recv_faces_vector_el(1,iinterface), &
             3*nibool_interfaces_elastic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_elastic(ninterface_elastic+iinterface), ier)

    if ( ier /= MPI_SUCCESS ) then
      call exit_mpi('MPI_IRECV unsuccessful in assemble_MPI_vector_el')
    endif

  enddo

  do iinterface = 1, ninterface_elastic*2

    call MPI_Wait(tab_requests_send_recv_elastic(iinterface), MPI_STATUS_IGNORE, ier)

  enddo

  do iinterface = 1, ninterface_elastic

     num_interface = inum_interfaces_elastic(iinterface)

     ipoin = 0
     do i = 1, nibool_interfaces_elastic(num_interface)
        array_val2(:,ibool_interfaces_elastic(i,num_interface)) = &
            array_val2(:,ibool_interfaces_elastic(i,num_interface))  &
            + buffer_recv_faces_vector_el(ipoin+1:ipoin+3,iinterface)
        ipoin = ipoin + 3
     enddo

  enddo

  end subroutine assemble_MPI_vector_el


!-----------------------------------------------
! Assembling accel_elastic for poroelastic elements :
! the buffers are filled, the ISEND and IRECV are started here, then
! contributions are added.
! The previous version included communication overlap using persistent
! communication, but the merging of the outer and inner elements rendered
! overlap no longer possible, while persistent communications were removed
! because trace tool MPITrace does not yet instrument those.
! Particular care should be taken concerning possible optimisations of the
! communication scheme.
!-----------------------------------------------
  subroutine assemble_MPI_vector_po(array_val3,array_val4)

  use mpi

  use specfem_par


  implicit none

  include 'precision.h'

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob), intent(inout) :: array_val3,array_val4

  integer  :: ipoin, num_interface, iinterface,i

  do iinterface = 1, ninterface_poroelastic

     num_interface = inum_interfaces_poroelastic(iinterface)

     ipoin = 0
     do i = 1, nibool_interfaces_poroelastic(num_interface)
        buffer_send_faces_vector_pos(ipoin+1:ipoin+2,iinterface) = &
             array_val3(:,ibool_interfaces_poroelastic(i,num_interface))
        ipoin = ipoin + 2
     enddo

     ipoin = 0
     do i = 1, nibool_interfaces_poroelastic(num_interface)
        buffer_send_faces_vector_pow(ipoin+1:ipoin+2,iinterface) = &
             array_val4(:,ibool_interfaces_poroelastic(i,num_interface))
        ipoin = ipoin + 2
     enddo

  enddo

  do iinterface = 1, ninterface_poroelastic

    num_interface = inum_interfaces_poroelastic(iinterface)

    call MPI_ISEND( buffer_send_faces_vector_pos(1,iinterface), &
             NDIM*nibool_interfaces_poroelastic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_poro(iinterface), ier)

    if ( ier /= MPI_SUCCESS ) then
      call exit_mpi('MPI_ISEND unsuccessful in assemble_MPI_vector_pos')
    endif

    call MPI_IRECV ( buffer_recv_faces_vector_pos(1,iinterface), &
             NDIM*nibool_interfaces_poroelastic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_poro(ninterface_poroelastic+iinterface), ier)

    if ( ier /= MPI_SUCCESS ) then
      call exit_mpi('MPI_IRECV unsuccessful in assemble_MPI_vector_pos')
    endif

    call MPI_ISEND( buffer_send_faces_vector_pow(1,iinterface), &
             NDIM*nibool_interfaces_poroelastic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_poro(ninterface_poroelastic*2+iinterface), ier)

    if ( ier /= MPI_SUCCESS ) then
      call exit_mpi('MPI_ISEND unsuccessful in assemble_MPI_vector_pow')
    endif

    call MPI_IRECV ( buffer_recv_faces_vector_pow(1,iinterface), &
             NDIM*nibool_interfaces_poroelastic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_poro(ninterface_poroelastic*3+iinterface), ier)

    if ( ier /= MPI_SUCCESS ) then
      call exit_mpi('MPI_IRECV unsuccessful in assemble_MPI_vector_pow')
    endif

  enddo

  do iinterface = 1, ninterface_poroelastic*4

    call MPI_Wait (tab_requests_send_recv_poro(iinterface), MPI_STATUS_IGNORE, ier)

  enddo

  do iinterface = 1, ninterface_poroelastic

     num_interface = inum_interfaces_poroelastic(iinterface)

     ipoin = 0
     do i = 1, nibool_interfaces_poroelastic(num_interface)
        array_val3(:,ibool_interfaces_poroelastic(i,num_interface)) = &
             array_val3(:,ibool_interfaces_poroelastic(i,num_interface)) + &
             buffer_recv_faces_vector_pos(ipoin+1:ipoin+2,iinterface)
        ipoin = ipoin + 2
     enddo

     ipoin = 0
     do i = 1, nibool_interfaces_poroelastic(num_interface)
        array_val4(:,ibool_interfaces_poroelastic(i,num_interface)) = &
             array_val4(:,ibool_interfaces_poroelastic(i,num_interface)) + &
             buffer_recv_faces_vector_pow(ipoin+1:ipoin+2,iinterface)
        ipoin = ipoin + 2
     enddo

  enddo

  end subroutine assemble_MPI_vector_po


#endif


#ifdef USE_MPI
 subroutine assemble_MPI_scalar_send_cuda(NPROC, &
                                           buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           my_neighbours_ext_mesh, &
                                           tab_requests_send_recv_ext_mesh,ninterface_acoustic,inum_interfaces_acoustic)

! non-blocking MPI send

  ! sends data
  ! note: assembling data already filled into buffer_send_scalar_ext_mesh array

  use constants

  implicit none

  integer :: NPROC,ninterface_acoustic
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(2*num_interfaces_ext_mesh) :: tab_requests_send_recv_ext_mesh
integer, dimension(num_interfaces_ext_mesh), intent(in)  :: inum_interfaces_acoustic
  ! local parameters
  integer :: iinterface,num_interface

  tab_requests_send_recv_ext_mesh(:) = 0


  ! sends only if more than one partition
  if(NPROC > 1) then

    ! note: partition border copy into the buffer has already been done
    !          by routine transfer_boun_pot_from_device()

    ! send messages
    do iinterface = 1, ninterface_acoustic

     num_interface=inum_interfaces_acoustic(iinterface)

      call isend_cr(buffer_send_scalar_ext_mesh(1,num_interface), &
                    nibool_interfaces_ext_mesh(num_interface), &
                    my_neighbours_ext_mesh(num_interface), &
                    itag, &
                    tab_requests_send_recv_ext_mesh(num_interface) )

      ! receive request
      call irecv_cr(buffer_recv_scalar_ext_mesh(1,num_interface), &
                    nibool_interfaces_ext_mesh(num_interface), &
                    my_neighbours_ext_mesh(num_interface), &
                    itag, &
                    tab_requests_send_recv_ext_mesh(num_interface+num_interfaces_ext_mesh) )
    enddo

  endif

  end subroutine assemble_MPI_scalar_send_cuda
#endif

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_write_cuda(NPROC, &
                        Mesh_pointer, &
                        buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
                        tab_requests_send_recv_ext_mesh, &
                        FORWARD_OR_ADJOINT,ninterface_acoustic,inum_interfaces_acoustic)

! waits for send/receiver to be completed and assembles contributions

  use constants

  implicit none

  integer :: NPROC
  integer :: ninterface_acoustic
  integer(kind=8) :: Mesh_pointer

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_scalar_ext_mesh

  integer, dimension(2*num_interfaces_ext_mesh) :: tab_requests_send_recv_ext_mesh
  integer, dimension(num_interfaces_ext_mesh), intent(in)  :: inum_interfaces_acoustic
  integer :: FORWARD_OR_ADJOINT

  integer :: iinterface,num_interface ! ipoin


  ! assemble only if more than one partition
  if(NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, ninterface_acoustic
   num_interface=inum_interfaces_acoustic(iinterface)

#ifdef USE_MPI
      call wait_req(tab_requests_send_recv_ext_mesh(num_interface+num_interfaces_ext_mesh))
#else
!! DK DK this dummy statement just to avoid a compiler warning about an unused variable
      tab_requests_send_recv_ext_mesh(num_interface+num_interfaces_ext_mesh) = 0
#endif
    enddo

    ! adding contributions of neighbours
    call transfer_asmbl_pot_to_device(Mesh_pointer,buffer_recv_scalar_ext_mesh,FORWARD_OR_ADJOINT)

    ! wait for communications completion (send)
     do iinterface = 1, ninterface_acoustic
       num_interface=inum_interfaces_acoustic(iinterface)

#ifdef USE_MPI
      call wait_req(tab_requests_send_recv_ext_mesh(num_interface))
#endif
    enddo

  endif

  end subroutine assemble_MPI_scalar_write_cuda



!
!-------------------------------------------------------------------------------------------------
!


#ifdef USE_MPI
  subroutine assemble_MPI_vector_send_cuda(NPROC, &
                                          buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh, &
                                          my_neighbours_ext_mesh, &
                                          tab_requests_send_recv_vector,ninterface_elastic,inum_interfaces_elastic)

! sends data
! note: array to assemble already filled into buffer_send_vector_ext_mesh array

  use constants

  implicit none

  integer :: NPROC

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,ninterface_elastic

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(2*num_interfaces_ext_mesh) :: tab_requests_send_recv_vector
  integer, dimension(num_interfaces_ext_mesh), intent(in)  :: inum_interfaces_elastic
  ! local parameters
  integer :: iinterface,num_interface

  ! note: preparation of the contribution between partitions using MPI
  !          already done in transfer_boun_accel routine

  ! send only if more than one partition
  if(NPROC > 1) then
    ! send messages
    do iinterface = 1, ninterface_elastic
       num_interface=inum_interfaces_elastic(iinterface)
      call isend_cr(buffer_send_vector_ext_mesh(1,1,num_interface), &
                     NDIM*nibool_interfaces_ext_mesh(num_interface), &
                     my_neighbours_ext_mesh(num_interface), &
                     itag, &
                     tab_requests_send_recv_vector(num_interface) )
      call irecv_cr(buffer_recv_vector_ext_mesh(1,1,num_interface), &
                     NDIM*nibool_interfaces_ext_mesh(num_interface), &
                     my_neighbours_ext_mesh(num_interface), &
                     itag, &
                     tab_requests_send_recv_vector(num_interface + num_interfaces_ext_mesh ) )
    enddo

  endif

  end subroutine assemble_MPI_vector_send_cuda
#endif


!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_write_cuda(NPROC,Mesh_pointer, &
                                            buffer_recv_vector_ext_mesh, &
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            tab_requests_send_recv_vector, &
                                            FORWARD_OR_ADJOINT,ninterface_elastic,inum_interfaces_elastic )

! waits for data to receive and assembles

  use constants

  implicit none

  integer :: NPROC
  integer(kind=8) :: Mesh_pointer



  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,ninterface_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(2*num_interfaces_ext_mesh) :: tab_requests_send_recv_vector
  integer, dimension(num_interfaces_ext_mesh), intent(in)  :: inum_interfaces_elastic
  integer :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface, num_interface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if(NPROC > 1) then


    ! wait for communications completion (recv)
do iinterface = 1, ninterface_elastic
       num_interface=inum_interfaces_elastic(iinterface)
#ifdef USE_MPI
      call wait_req(tab_requests_send_recv_vector(num_interface + num_interfaces_ext_mesh))
#else
!! DK DK this dummy statement just to avoid a compiler warning about an unused variable
      tab_requests_send_recv_vector(num_interface + num_interfaces_ext_mesh) = 0
#endif
enddo



    ! adding contributions of neighbours
    call transfer_asmbl_accel_to_device(Mesh_pointer, &
                                        buffer_recv_vector_ext_mesh, &
                                        max_nibool_interfaces_ext_mesh, &
                                        nibool_interfaces_ext_mesh,&
                                        ibool_interfaces_ext_mesh,FORWARD_OR_ADJOINT)

    ! This step is done via previous function transfer_and_assemble...
    ! do iinterface = 1, num_interfaces_ext_mesh
    !   do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
    !     array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
    !          array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
    !   enddo
    ! enddo

    ! wait for communications completion (send)
    do iinterface = 1, ninterface_elastic
       num_interface=inum_interfaces_elastic(iinterface)
#ifdef USE_MPI
      call wait_req(tab_requests_send_recv_vector(num_interface))
#endif
    enddo

  endif

  end subroutine assemble_MPI_vector_write_cuda


!
!-------------------------------------------------------------------------------------------------
!

!
!-------------------------------------------------------------------------------------------------
!

! with cuda functions...

  subroutine transfer_boundary_to_device(NPROC, Mesh_pointer, &
                                            buffer_recv_vector_ext_mesh, &
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,&
                                            tab_requests_send_recv_vector,ninterface_elastic,inum_interfaces_elastic )

  use constants

  implicit none

  integer :: NPROC
  integer(kind=8) :: Mesh_pointer

  ! array to assemble
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,ninterface_elastic

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_ext_mesh

  integer, dimension(2*num_interfaces_ext_mesh) :: tab_requests_send_recv_vector
  integer, dimension(num_interfaces_ext_mesh), intent(in)  :: inum_interfaces_elastic
  ! local parameters
  integer :: iinterface,num_interface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if(NPROC > 1) then



    ! wait for communications completion (recv)
    !write(IMAIN,*) "sending MPI_wait"
    do iinterface = 1, ninterface_elastic
       num_interface=inum_interfaces_elastic(iinterface)
#ifdef USE_MPI
      call wait_req(tab_requests_send_recv_vector(num_interface+num_interfaces_ext_mesh))
#else
!! DK DK this dummy statement just to avoid a compiler warning about an unused variable
      tab_requests_send_recv_vector(num_interface+num_interfaces_ext_mesh) = 0
#endif
    enddo




    ! send contributions to GPU
    call transfer_boundary_to_device_a(Mesh_pointer, buffer_recv_vector_ext_mesh, max_nibool_interfaces_ext_mesh)
  endif

  ! This step is done via previous function transfer_and_assemble...
  ! do iinterface = 1, num_interfaces_ext_mesh
  !   do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
  !     array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
  !          array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
  !   enddo
  ! enddo

  end subroutine transfer_boundary_to_device

!
!-------------------------------------------------------------------------------------------------
!

