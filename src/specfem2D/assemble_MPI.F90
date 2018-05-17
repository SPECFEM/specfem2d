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

!
! This file contains subroutines related to assembling (of the mass matrix, potential_dot_dot and
! accel_elastic, accels_poroelastic, accelw_poroelastic).
!
! These subroutines are for the most part not used in the sequential version.


#ifdef USE_MPI
! only supported with parallel version...

!-----------------------------------------------
! Assembling the mass matrix.
!-----------------------------------------------
  subroutine assemble_MPI_scalar(array_val1,npoin_val1, &
                                 array_e1,n_sls_loc, &
                                 array_val2,npoin_val2, &
                                 array_val3,array_val4,npoin_val3)

  use mpi

  use constants, only: CUSTOM_REAL,NDIM,USE_A_STRONG_FORMULATION_FOR_E1

  use specfem_par, only: NPROC,ninterface,my_neighbors,N_SLS,ATTENUATION_VISCOACOUSTIC

  ! acoustic/elastic/poroelastic interfaces
  use specfem_par, only: max_ibool_interfaces_size_ac,max_ibool_interfaces_size_el,max_ibool_interfaces_size_po, &
    ibool_interfaces_acoustic,ibool_interfaces_elastic,ibool_interfaces_poroelastic, &
    nibool_interfaces_acoustic,nibool_interfaces_elastic,nibool_interfaces_poroelastic

  implicit none

  include "precision.h"

  ! arrays to assemble

  ! acoustic
  integer :: npoin_val1
  real(kind=CUSTOM_REAL), dimension(npoin_val1), intent(inout) :: array_val1

  ! viscoacoustic
  real(kind=CUSTOM_REAL), dimension(npoin_val1,N_SLS), intent(inout) :: array_e1
  integer :: i_sls,n_sls_loc

  ! elastic
  integer :: npoin_val2
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin_val2), intent(inout) :: array_val2

  ! poroelastic
  integer :: npoin_val3
  real(kind=CUSTOM_REAL), dimension(npoin_val3), intent(inout) :: array_val3,array_val4

  ! local parameters
  integer :: ipoin,iglob,i,idim,iinterface
  integer :: nbuffer_points

  ! there are now two different mass matrices for the elastic case
  ! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
  real(kind=CUSTOM_REAL), dimension(max_ibool_interfaces_size_ac + &
    NDIM * max_ibool_interfaces_size_el + &
    2 * max_ibool_interfaces_size_po + &
    n_sls_loc*max_ibool_interfaces_size_ac, ninterface)  :: buffer_send_faces_scalar, buffer_recv_faces_scalar

  integer, dimension(ninterface)  :: msg_send_requests,msg_recv_requests

  ! assemble only if more than one partition
  if (NPROC > 1) then

    buffer_send_faces_scalar(:,:) = 0._CUSTOM_REAL
    buffer_recv_faces_scalar(:,:) = 0._CUSTOM_REAL

    do iinterface = 1, ninterface
      ! acoustic
      ipoin = 0
      do i = 1, nibool_interfaces_acoustic(iinterface)
        ipoin = ipoin + 1
        iglob = ibool_interfaces_acoustic(i,iinterface)
        buffer_send_faces_scalar(ipoin,iinterface) = array_val1(iglob)
      enddo

      ! elastic
      ! there are now two different mass matrices for the elastic case
      ! in order to handle the C deltat / 2 contribution of the Stacey conditions to the mass matrix
      do idim = 1,NDIM
        do i = 1, nibool_interfaces_elastic(iinterface)
          ipoin = ipoin + 1
          iglob = ibool_interfaces_elastic(i,iinterface)
          buffer_send_faces_scalar(ipoin,iinterface) = array_val2(idim,iglob)
        enddo
      enddo

      ! poroelastic
      do i = 1, nibool_interfaces_poroelastic(iinterface)
        ipoin = ipoin + 1
        iglob = ibool_interfaces_poroelastic(i,iinterface)
        buffer_send_faces_scalar(ipoin,iinterface) = array_val3(iglob)
      enddo
      do i = 1, nibool_interfaces_poroelastic(iinterface)
        ipoin = ipoin + 1
        iglob = ibool_interfaces_poroelastic(i,iinterface)
        buffer_send_faces_scalar(ipoin,iinterface) = array_val4(iglob)
      enddo

      ! total number of points in buffer
      nbuffer_points = nibool_interfaces_acoustic(iinterface) &
                      + NDIM*nibool_interfaces_elastic(iinterface) &
                      + 2*nibool_interfaces_poroelastic(iinterface)

      if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1)) then
        ! viscoacoustic

        ! loop over relaxation mechanisms
        do i_sls = 1,N_SLS

          do i = 1, nibool_interfaces_acoustic(iinterface)
            ipoin = ipoin + 1
            iglob = ibool_interfaces_acoustic(i,iinterface)
            buffer_send_faces_scalar(ipoin,iinterface) = array_e1(iglob,i_sls)
          enddo

          ! update total number of points in buffer
          nbuffer_points = nbuffer_points + nibool_interfaces_acoustic(iinterface)

        enddo
      endif

      ! non-blocking send
      call isend_cr(buffer_send_faces_scalar(1,iinterface), nbuffer_points, my_neighbors(iinterface),11, &
                    msg_send_requests(iinterface))

      ! starts a blocking receive
      call irecv_cr(buffer_recv_faces_scalar(1,iinterface),nbuffer_points,my_neighbors(iinterface),11, &
                    msg_recv_requests(iinterface))

    enddo

    ! waits for MPI requests to complete (recv)
    ! each wait returns once the specified MPI request completed
    do iinterface = 1, ninterface
      call wait_req(msg_recv_requests(iinterface))
    enddo

    do iinterface = 1, ninterface
      ! acoustic
      ipoin = 0
      do i = 1, nibool_interfaces_acoustic(iinterface)
        ipoin = ipoin + 1
        iglob = ibool_interfaces_acoustic(i,iinterface)
        array_val1(iglob) = array_val1(iglob) + buffer_recv_faces_scalar(ipoin,iinterface)
      enddo

      ! elastic
      do idim = 1,NDIM
        do i = 1, nibool_interfaces_elastic(iinterface)
          ipoin = ipoin + 1
          iglob = ibool_interfaces_elastic(i,iinterface)
          array_val2(idim,iglob) = array_val2(idim,iglob) + buffer_recv_faces_scalar(ipoin,iinterface)
        enddo
      enddo

      ! poroelastic
      do i = 1, nibool_interfaces_poroelastic(iinterface)
        ipoin = ipoin + 1
        iglob = ibool_interfaces_poroelastic(i,iinterface)
        array_val3(iglob) = array_val3(iglob) + buffer_recv_faces_scalar(ipoin,iinterface)
      enddo
      do i = 1, nibool_interfaces_poroelastic(iinterface)
        ipoin = ipoin + 1
        iglob = ibool_interfaces_poroelastic(i,iinterface)
        array_val4(iglob) = array_val4(iglob) + buffer_recv_faces_scalar(ipoin,iinterface)
      enddo

      if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1)) then
        ! loop over relaxation mechanisms
        do i_sls = 1,N_SLS
          do i = 1, nibool_interfaces_acoustic(iinterface)
            ipoin = ipoin + 1
            iglob = ibool_interfaces_acoustic(i,iinterface)
            array_e1(iglob,i_sls) = array_e1(iglob,i_sls) + buffer_recv_faces_scalar(ipoin,iinterface)
          enddo
        enddo
      endif

    enddo

  endif ! NPROC > 1

  end subroutine assemble_MPI_scalar

#endif

!-------------------------------------------------------------------------------------------------
!
! acoustic domains
!
!-------------------------------------------------------------------------------------------------


!-----------------------------------------------
! Assembling potential_dot_dot for acoustic elements :
! the buffers are filled, the ISEND and IRECV are started here, then
! contributions are added.
! The previous version included communication overlap using persistent
! communication, but the merging of the outer and inner elements rendered
! overlap no longer possible, while persistent communications were removed
! because trace tool MPITrace does not yet instrument those.
! Particular care should be taken concerning possible optimizations of the
! communication scheme.
!-----------------------------------------------
  subroutine assemble_MPI_scalar_ac_blocking(array_val1)

  use constants, only: CUSTOM_REAL

  use specfem_par, only: NPROC,nglob,my_neighbors

  ! acoustic MPI interfaces
  use specfem_par, only: ninterface_acoustic,inum_interfaces_acoustic, &
    ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
    request_send_recv_acoustic,buffer_send_faces_vector_ac,buffer_recv_faces_vector_ac

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: array_val1

  ! local parameters
  integer  :: ipoin, num_interface,iinterface, iglob
  integer, parameter :: itag = 12

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! initializes buffers
    buffer_send_faces_vector_ac(:,:) = 0._CUSTOM_REAL
    buffer_recv_faces_vector_ac(:,:) = 0._CUSTOM_REAL
    request_send_recv_acoustic(:) = 0

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
      call isend_cr( buffer_send_faces_vector_ac(1,iinterface), &
                     nibool_interfaces_acoustic(num_interface), &
                     my_neighbors(num_interface), &
                     itag, &
                     request_send_recv_acoustic(iinterface) )

      ! starts a non-blocking receive
      call irecv_cr( buffer_recv_faces_vector_ac(1,iinterface), &
                     nibool_interfaces_acoustic(num_interface), &
                     my_neighbors(num_interface), &
                     itag, &
                     request_send_recv_acoustic(ninterface_acoustic+iinterface) )
    enddo

    ! waits for MPI requests to complete (recv)
    ! each wait returns once the specified MPI request completed
    do iinterface = 1, ninterface_acoustic
      call wait_req(request_send_recv_acoustic(ninterface_acoustic+iinterface))
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
      call wait_req(request_send_recv_acoustic(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_scalar_ac_blocking

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_ac_s(array_val1)

! sends MPI buffers (asynchronuous) - non-blocking routine

  use constants, only: CUSTOM_REAL

  use specfem_par, only: NPROC,nglob,my_neighbors

  ! acoustic MPI interfaces
  use specfem_par, only: ninterface_acoustic,inum_interfaces_acoustic, &
    ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
    request_send_recv_acoustic,buffer_send_faces_vector_ac,buffer_recv_faces_vector_ac

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: array_val1

  ! local parameters
  integer  :: ipoin, num_interface,iinterface, iglob
  integer, parameter :: itag = 12

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! initializes buffers
    buffer_send_faces_vector_ac(:,:) = 0._CUSTOM_REAL
    buffer_recv_faces_vector_ac(:,:) = 0._CUSTOM_REAL

    request_send_recv_acoustic(:) = 0

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

    ! send/request messages
    do iinterface = 1, ninterface_acoustic

      ! gets global interface index
      num_interface = inum_interfaces_acoustic(iinterface)

      call isend_cr( buffer_send_faces_vector_ac(1,iinterface), &
                     nibool_interfaces_acoustic(num_interface), &
                     my_neighbors(num_interface), &
                     itag, &
                     request_send_recv_acoustic(iinterface) )

      call irecv_cr( buffer_recv_faces_vector_ac(1,iinterface), &
                     nibool_interfaces_acoustic(num_interface), &
                     my_neighbors(num_interface), &
                     itag, &
                     request_send_recv_acoustic(ninterface_acoustic+iinterface) )

    enddo

  endif

  end subroutine assemble_MPI_scalar_ac_s

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_ac_w(array_val1)

! waits for data and assembles

  use constants, only: CUSTOM_REAL

  use specfem_par, only: NPROC,nglob

  ! acoustic MPI interfaces
  use specfem_par, only: ninterface_acoustic,inum_interfaces_acoustic, &
    ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
    request_send_recv_acoustic,buffer_recv_faces_vector_ac

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: array_val1

  ! local parameters
  integer  :: ipoin, num_interface,iinterface, iglob

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! waits for communication complection (all receive has finished)
    do iinterface = 1, ninterface_acoustic
      call wait_req(request_send_recv_acoustic(ninterface_acoustic+iinterface))
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

    ! waits for communication completion (all send has finished, we can clear the buffer...)
    do iinterface = 1, ninterface_acoustic
      call wait_req(request_send_recv_acoustic(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_scalar_ac_w

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_ac_s_e1(array_val1,array_val1_e1,n_sls_local)

! sends MPI buffers (asynchronously) - non-blocking routine

  use constants, only: CUSTOM_REAL,USE_A_STRONG_FORMULATION_FOR_E1

  use specfem_par, only: NPROC,nglob,my_neighbors

  ! acoustic MPI interfaces
  use specfem_par, only: ninterface_acoustic,inum_interfaces_acoustic, &
    ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
    request_send_recv_acoustic,buffer_send_faces_vector_ac,buffer_recv_faces_vector_ac, &
    N_SLS,ATTENUATION_VISCOACOUSTIC,nglob_att

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: array_val1
  real(kind=CUSTOM_REAL), dimension(nglob_att,N_SLS), intent(inout) :: array_val1_e1

  integer :: n_sls_local,n_sls_local_copy

  ! local parameters
  integer  :: ipoin, num_interface,iinterface, iglob, i
  integer, parameter :: itag = 12

  ! this is to avoid setting N_SLS to zero in the rest of the code
  ! when it is necessary to set it to zero temporarily here only
  n_sls_local_copy = n_sls_local

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! n_sls_local = 0 corresponds to a Newmark scheme, or to no viscoacousticity
    if (n_sls_local_copy > 0 .and. (.not. ATTENUATION_VISCOACOUSTIC .or. &
       (ATTENUATION_VISCOACOUSTIC .and. ( USE_A_STRONG_FORMULATION_FOR_E1)))) n_sls_local_copy = 0


    ! initializes buffers
    buffer_send_faces_vector_ac(:,:) = 0._CUSTOM_REAL
    buffer_recv_faces_vector_ac(:,:) = 0._CUSTOM_REAL

    request_send_recv_acoustic(:) = 0

    ! loops over acoustic interfaces only
    do iinterface = 1, ninterface_acoustic

      ! gets interface index in the range of all interfaces [1,ninterface]
      num_interface = inum_interfaces_acoustic(iinterface)

      ipoin = 0
      ! loops over all interface points
      do i = 1, nibool_interfaces_acoustic(num_interface)
        iglob = ibool_interfaces_acoustic(i,num_interface)

        ipoin = ipoin + 1
        ! copies array values to buffer
        buffer_send_faces_vector_ac(ipoin,iinterface) = array_val1(iglob)
        if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1) .and. n_sls_local_copy > 0) then
        buffer_send_faces_vector_ac(ipoin+1:ipoin+N_SLS,iinterface) = array_val1_e1(iglob,:)
        ipoin = ipoin + N_SLS
        endif
      enddo

    enddo

    ! send/request messages
    do iinterface = 1, ninterface_acoustic

      ! gets global interface index
      num_interface = inum_interfaces_acoustic(iinterface)

      call isend_cr( buffer_send_faces_vector_ac(1,iinterface), &
                     nibool_interfaces_acoustic(num_interface)*(n_sls_local_copy+1), &
                     my_neighbors(num_interface), &
                     itag, &
                     request_send_recv_acoustic(iinterface) )

      call irecv_cr( buffer_recv_faces_vector_ac(1,iinterface), &
                     nibool_interfaces_acoustic(num_interface)*(n_sls_local_copy+1), &
                     my_neighbors(num_interface), &
                     itag, &
                     request_send_recv_acoustic(ninterface_acoustic+iinterface) )

    enddo

  endif

  end subroutine assemble_MPI_scalar_ac_s_e1

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_ac_w_e1(array_val1,array_val1_e1,n_sls_local)

! waits for data and assembles

  use constants, only: CUSTOM_REAL,USE_A_STRONG_FORMULATION_FOR_E1

  use specfem_par, only: NPROC,nglob

  ! acoustic MPI interfaces
  use specfem_par, only: ninterface_acoustic,inum_interfaces_acoustic, &
    ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
    request_send_recv_acoustic,buffer_recv_faces_vector_ac, &
    N_SLS,ATTENUATION_VISCOACOUSTIC,nglob_att

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: array_val1
  real(kind=CUSTOM_REAL), dimension(nglob_att,N_SLS), intent(inout) :: array_val1_e1

  integer :: n_sls_local,n_sls_local_copy

  ! local parameters
  integer  :: ipoin, num_interface,iinterface, iglob, i

  ! this is to avoid setting N_SLS to zero in the rest of the code
  ! when it is necessary to set it to zero temporarily here only
  n_sls_local_copy = n_sls_local

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! n_sls_local = 0 corresponds to a Newmark scheme, or to no viscoacousticity
    if (n_sls_local_copy > 0 .and. ((.not. ATTENUATION_VISCOACOUSTIC) .or. &
       (ATTENUATION_VISCOACOUSTIC .and. ( USE_A_STRONG_FORMULATION_FOR_E1)))) n_sls_local_copy = 0

    ! waits for communication complection (all receive has finished)
    do iinterface = 1, ninterface_acoustic
      call wait_req(request_send_recv_acoustic(ninterface_acoustic+iinterface))
    enddo

    ! assembles the array values
    do iinterface = 1, ninterface_acoustic

      ! gets global interface index
      num_interface = inum_interfaces_acoustic(iinterface)

      ipoin = 0
      ! loops over all interface points
      do i = 1, nibool_interfaces_acoustic(num_interface)
        ipoin = ipoin + 1
        iglob = ibool_interfaces_acoustic(i,num_interface)
        ! adds buffer contribution
        array_val1(iglob) = array_val1(iglob) + buffer_recv_faces_vector_ac(ipoin,iinterface)

        if (ATTENUATION_VISCOACOUSTIC .and. (.not. USE_A_STRONG_FORMULATION_FOR_E1) .and. n_sls_local_copy > 0) then
        array_val1_e1(iglob,:) = array_val1_e1(iglob,:) &
                + buffer_recv_faces_vector_ac(ipoin+1:ipoin+N_SLS,iinterface)
        ipoin = ipoin + N_SLS
        endif
      enddo

    enddo

    ! waits for communication completion (all send has finished, we can clear the buffer...)
    do iinterface = 1, ninterface_acoustic
      call wait_req(request_send_recv_acoustic(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_scalar_ac_w_e1

!-------------------------------------------------------------------------------------------------
!
! elastic domains
!
!-------------------------------------------------------------------------------------------------


#ifdef USE_MPI
! only supported with parallel version...

!-----------------------------------------------
! Assembling accel_elastic for elastic elements :
! the buffers are filled, the ISEND and IRECV are started here, then
! contributions are added.
! The previous version included communication overlap using persistent
! communication, but the merging of the outer and inner elements rendered
! overlap no longer possible, while persistent communications were removed
! because trace tool MPITrace does not yet instrument those.
! Particular care should be taken concerning possible optimizations of the
! communication scheme.
!-----------------------------------------------
! note: although it uses asynchronuous isend/irecv, the routine is waiting to have all buffers transfered.
!       the scheme is therefore similar to a blocking MPI scheme.

  subroutine assemble_MPI_vector_el_blocking(array_val)

  use mpi

  use constants, only: CUSTOM_REAL,NDIM

  use specfem_par, only: NPROC, &
    nglob,ninterface_elastic, &
    inum_interfaces_elastic, &
    ibool_interfaces_elastic, nibool_interfaces_elastic, &
    request_send_recv_elastic, &
    buffer_send_faces_vector_el, &
    buffer_recv_faces_vector_el, &
    my_neighbors

  implicit none

  include "precision.h"

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob), intent(inout) :: array_val

  ! local parameters
  integer  :: ipoin, num_interface, iinterface, i

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! fills buffer
    do iinterface = 1, ninterface_elastic

       num_interface = inum_interfaces_elastic(iinterface)

       ipoin = 0
       do i = 1, nibool_interfaces_elastic(num_interface)
          buffer_send_faces_vector_el(ipoin+1:ipoin+NDIM,iinterface) = array_val(:,ibool_interfaces_elastic(i,num_interface))
          ipoin = ipoin + NDIM
       enddo
    enddo

    ! send/receive with neighbors
    do iinterface = 1, ninterface_elastic

      num_interface = inum_interfaces_elastic(iinterface)

      call isend_cr(buffer_send_faces_vector_el(1,iinterface),NDIM * nibool_interfaces_elastic(num_interface), &
                      my_neighbors(num_interface),12,request_send_recv_elastic(iinterface))

      call irecv_cr(buffer_recv_faces_vector_el(1,iinterface),NDIM * nibool_interfaces_elastic(num_interface), &
                    my_neighbors(num_interface), 12, request_send_recv_elastic(ninterface_elastic+iinterface))

    enddo

    ! waits for all send/receive has finished
    do iinterface = 1, ninterface_elastic*2
      call wait_req(request_send_recv_elastic(iinterface))
    enddo

    ! adds contributions
    do iinterface = 1, ninterface_elastic

       num_interface = inum_interfaces_elastic(iinterface)

       ipoin = 0
       do i = 1, nibool_interfaces_elastic(num_interface)
          array_val(:,ibool_interfaces_elastic(i,num_interface)) = &
              array_val(:,ibool_interfaces_elastic(i,num_interface)) + buffer_recv_faces_vector_el(ipoin+1:ipoin+NDIM,iinterface)
          ipoin = ipoin + NDIM
       enddo

    enddo

  endif

  end subroutine assemble_MPI_vector_el_blocking

#endif

!
!-------------------------------------------------------------------------------------------------
!


  subroutine assemble_MPI_vector_el_s(array_val)

! sends MPI buffers (asynchronuous) - non-blocking routine

  use constants, only: CUSTOM_REAL,NDIM

  use specfem_par, only: NPROC, &
    nglob,ninterface_elastic, &
    inum_interfaces_elastic, &
    ibool_interfaces_elastic, nibool_interfaces_elastic, &
    request_send_recv_elastic, &
    buffer_send_faces_vector_el, &
    buffer_recv_faces_vector_el, &
    my_neighbors

  implicit none

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob), intent(inout) :: array_val

  ! local parameters
  integer  :: ipoin, num_interface, iinterface, i
  integer, parameter :: itag = 12

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! fills buffer
    do iinterface = 1, ninterface_elastic

       num_interface = inum_interfaces_elastic(iinterface)

       ipoin = 0
       do i = 1, nibool_interfaces_elastic(num_interface)
          buffer_send_faces_vector_el(ipoin+1:ipoin+NDIM,iinterface) = array_val(:,ibool_interfaces_elastic(i,num_interface))
          ipoin = ipoin + NDIM
       enddo
    enddo

    ! send/request messages
    do iinterface = 1, ninterface_elastic

      num_interface = inum_interfaces_elastic(iinterface)

      call isend_cr( buffer_send_faces_vector_el(1,iinterface), &
                     NDIM * nibool_interfaces_elastic(num_interface), &
                     my_neighbors(num_interface), &
                     itag, &
                     request_send_recv_elastic(iinterface) )

      call irecv_cr( buffer_recv_faces_vector_el(1,iinterface), &
                     NDIM * nibool_interfaces_elastic(num_interface), &
                     my_neighbors(num_interface), &
                     itag, &
                     request_send_recv_elastic(ninterface_elastic+iinterface) )
    enddo

  endif

  end subroutine assemble_MPI_vector_el_s

!
!-------------------------------------------------------------------------------------------------
!


  subroutine assemble_MPI_vector_el_w(array_val)

! waits for data and assembles

  use constants, only: CUSTOM_REAL,NDIM

  use specfem_par, only: NPROC, &
    nglob,ninterface_elastic, &
    inum_interfaces_elastic, &
    ibool_interfaces_elastic, nibool_interfaces_elastic, &
    request_send_recv_elastic,buffer_recv_faces_vector_el

  implicit none

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob), intent(inout) :: array_val

  ! local parameters
  integer  :: ipoin, num_interface, iinterface, i

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! waits for communication complection (all receive has finished)
    do iinterface = 1, ninterface_elastic
      call wait_req(request_send_recv_elastic(ninterface_elastic+iinterface))
    enddo

    ! adds contributions
    do iinterface = 1, ninterface_elastic

       num_interface = inum_interfaces_elastic(iinterface)

       ipoin = 0
       do i = 1, nibool_interfaces_elastic(num_interface)
          array_val(:,ibool_interfaces_elastic(i,num_interface)) = &
              array_val(:,ibool_interfaces_elastic(i,num_interface)) + buffer_recv_faces_vector_el(ipoin+1:ipoin+NDIM,iinterface)
          ipoin = ipoin + NDIM
       enddo

    enddo

    ! waits for communication completion (all send has finished, we can clear the buffer...)
    do iinterface = 1, ninterface_elastic
      call wait_req(request_send_recv_elastic(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_el_w



!-------------------------------------------------------------------------------------------------
!
! poroelastic domains
!
!-------------------------------------------------------------------------------------------------

#ifdef USE_MPI
! only supported with parallel version...

!-----------------------------------------------
! Assembling accel_elastic for poroelastic elements :
! the buffers are filled, the ISEND and IRECV are started here, then
! contributions are added.
! The previous version included communication overlap using persistent
! communication, but the merging of the outer and inner elements rendered
! overlap no longer possible, while persistent communications were removed
! because trace tool MPITrace does not yet instrument those.
! Particular care should be taken concerning possible optimizations of the
! communication scheme.
!-----------------------------------------------
  subroutine assemble_MPI_vector_po_blocking(array_val3,array_val4)

  use mpi

  use constants, only: CUSTOM_REAL,NDIM

  use specfem_par, only: NPROC, &
    nglob,ninterface_poroelastic, &
    inum_interfaces_poroelastic, &
    ibool_interfaces_poroelastic, nibool_interfaces_poroelastic, &
    request_send_recv_poro, &
    buffer_send_faces_vector_pos,buffer_send_faces_vector_pow, &
    buffer_recv_faces_vector_pos,buffer_recv_faces_vector_pow, &
    my_neighbors

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob), intent(inout) :: array_val3,array_val4

  integer  :: ipoin, num_interface, iinterface, i

  ! assemble only if more than one partition
  if (NPROC > 1) then

    do iinterface = 1, ninterface_poroelastic

       num_interface = inum_interfaces_poroelastic(iinterface)

       ipoin = 0
       do i = 1, nibool_interfaces_poroelastic(num_interface)
          buffer_send_faces_vector_pos(ipoin+1:ipoin+NDIM,iinterface) = &
               array_val3(:,ibool_interfaces_poroelastic(i,num_interface))
          ipoin = ipoin + NDIM
       enddo

       ipoin = 0
       do i = 1, nibool_interfaces_poroelastic(num_interface)
          buffer_send_faces_vector_pow(ipoin+1:ipoin+NDIM,iinterface) = &
               array_val4(:,ibool_interfaces_poroelastic(i,num_interface))
          ipoin = ipoin + NDIM
       enddo

    enddo

    do iinterface = 1, ninterface_poroelastic

      num_interface = inum_interfaces_poroelastic(iinterface)

      call isend_cr(buffer_send_faces_vector_pos(1,iinterface),NDIM*nibool_interfaces_poroelastic(num_interface), &
                    my_neighbors(num_interface),12,request_send_recv_poro(iinterface))

      call irecv_cr(buffer_recv_faces_vector_pos(1,iinterface),NDIM*nibool_interfaces_poroelastic(num_interface), &
                       my_neighbors(num_interface),12,request_send_recv_poro(ninterface_poroelastic+iinterface))

      call isend_cr(buffer_send_faces_vector_pow(1,iinterface),NDIM*nibool_interfaces_poroelastic(num_interface), &
                    my_neighbors(num_interface),12,request_send_recv_poro(ninterface_poroelastic*2+iinterface))

      call irecv_cr(buffer_recv_faces_vector_pow(1,iinterface),NDIM*nibool_interfaces_poroelastic(num_interface), &
                       my_neighbors(num_interface),12,request_send_recv_poro(ninterface_poroelastic*3+iinterface))

    enddo

    do iinterface = 1, ninterface_poroelastic*4
      call wait_req(request_send_recv_poro(iinterface))
    enddo

    do iinterface = 1, ninterface_poroelastic

       num_interface = inum_interfaces_poroelastic(iinterface)

       ipoin = 0
       do i = 1, nibool_interfaces_poroelastic(num_interface)
          array_val3(:,ibool_interfaces_poroelastic(i,num_interface)) = &
               array_val3(:,ibool_interfaces_poroelastic(i,num_interface)) + &
               buffer_recv_faces_vector_pos(ipoin+1:ipoin+NDIM,iinterface)
          ipoin = ipoin + NDIM
       enddo

       ipoin = 0
       do i = 1, nibool_interfaces_poroelastic(num_interface)
          array_val4(:,ibool_interfaces_poroelastic(i,num_interface)) = &
               array_val4(:,ibool_interfaces_poroelastic(i,num_interface)) + &
               buffer_recv_faces_vector_pow(ipoin+1:ipoin+NDIM,iinterface)
          ipoin = ipoin + NDIM
       enddo

    enddo

  endif

  end subroutine assemble_MPI_vector_po_blocking

#endif

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_po_s(array_val3,array_val4)

! sends MPI buffers (asynchronuous) - non-blocking routine

  use constants, only: CUSTOM_REAL,NDIM

  use specfem_par, only: NPROC, &
    nglob,ninterface_poroelastic, &
    inum_interfaces_poroelastic, &
    ibool_interfaces_poroelastic, nibool_interfaces_poroelastic, &
    request_send_recv_poro, &
    buffer_send_faces_vector_pos,buffer_send_faces_vector_pow, &
    buffer_recv_faces_vector_pos,buffer_recv_faces_vector_pow, &
    my_neighbors

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob), intent(inout) :: array_val3,array_val4

  integer  :: ipoin, num_interface, iinterface, i
  integer, parameter :: itag = 12

  ! assemble only if more than one partition
  if (NPROC > 1) then

    do iinterface = 1, ninterface_poroelastic

       num_interface = inum_interfaces_poroelastic(iinterface)

       ipoin = 0
       do i = 1, nibool_interfaces_poroelastic(num_interface)
          buffer_send_faces_vector_pos(ipoin+1:ipoin+NDIM,iinterface) = &
               array_val3(:,ibool_interfaces_poroelastic(i,num_interface))
          ipoin = ipoin + NDIM
       enddo

       ipoin = 0
       do i = 1, nibool_interfaces_poroelastic(num_interface)
          buffer_send_faces_vector_pow(ipoin+1:ipoin+NDIM,iinterface) = &
               array_val4(:,ibool_interfaces_poroelastic(i,num_interface))
          ipoin = ipoin + NDIM
       enddo

    enddo

    ! send/request messages
    do iinterface = 1, ninterface_poroelastic

      num_interface = inum_interfaces_poroelastic(iinterface)

      ! solid
      call isend_cr( buffer_send_faces_vector_pos(1,iinterface), &
                     NDIM * nibool_interfaces_poroelastic(num_interface), &
                     my_neighbors(num_interface), &
                     itag, &
                     request_send_recv_poro(iinterface) )

      call irecv_cr( buffer_recv_faces_vector_pos(1,iinterface), &
                     NDIM * nibool_interfaces_poroelastic(num_interface), &
                     my_neighbors(num_interface), &
                     itag, &
                     request_send_recv_poro(ninterface_poroelastic+iinterface) )

      ! fluid
      call isend_cr( buffer_send_faces_vector_pow(1,iinterface), &
                     NDIM * nibool_interfaces_poroelastic(num_interface), &
                     my_neighbors(num_interface), &
                     itag, &
                     request_send_recv_poro(ninterface_poroelastic*2+iinterface) )


      call irecv_cr( buffer_recv_faces_vector_pow(1,iinterface), &
                     NDIM * nibool_interfaces_poroelastic(num_interface), &
                     my_neighbors(num_interface), &
                     itag, &
                     request_send_recv_poro(ninterface_poroelastic*3+iinterface) )

    enddo

  endif

  end subroutine assemble_MPI_vector_po_s

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_po_w(array_val3,array_val4)

! waits for data and assembles

  use constants, only: CUSTOM_REAL,NDIM

  use specfem_par, only: NPROC, &
    nglob,ninterface_poroelastic, &
    inum_interfaces_poroelastic, &
    ibool_interfaces_poroelastic, nibool_interfaces_poroelastic, &
    request_send_recv_poro, &
    buffer_recv_faces_vector_pos,buffer_recv_faces_vector_pow

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob), intent(inout) :: array_val3,array_val4

  integer  :: ipoin, num_interface, iinterface, i

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! waits for communication complection (all receive has finished)
    do iinterface = 1, ninterface_poroelastic
      call wait_req(request_send_recv_poro(ninterface_poroelastic+iinterface))
      call wait_req(request_send_recv_poro(ninterface_poroelastic*3+iinterface))
    enddo

    do iinterface = 1, ninterface_poroelastic

       num_interface = inum_interfaces_poroelastic(iinterface)

       ipoin = 0
       do i = 1, nibool_interfaces_poroelastic(num_interface)
          array_val3(:,ibool_interfaces_poroelastic(i,num_interface)) = &
               array_val3(:,ibool_interfaces_poroelastic(i,num_interface)) + &
               buffer_recv_faces_vector_pos(ipoin+1:ipoin+NDIM,iinterface)
          ipoin = ipoin + NDIM
       enddo

       ipoin = 0
       do i = 1, nibool_interfaces_poroelastic(num_interface)
          array_val4(:,ibool_interfaces_poroelastic(i,num_interface)) = &
               array_val4(:,ibool_interfaces_poroelastic(i,num_interface)) + &
               buffer_recv_faces_vector_pow(ipoin+1:ipoin+NDIM,iinterface)
          ipoin = ipoin + NDIM
       enddo

    enddo

    ! waits for communication completion (all send has finished, we can clear the buffer...)
    do iinterface = 1, ninterface_poroelastic
      call wait_req(request_send_recv_poro(iinterface))
      call wait_req(request_send_recv_poro(ninterface_poroelastic*2+iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_po_w



!-------------------------------------------------------------------------------------------------
!
! CUDA related
!
!-------------------------------------------------------------------------------------------------


 subroutine assemble_MPI_scalar_send_cuda(NPROC, &
                                           buffer_send_scalar_gpu,buffer_recv_scalar_gpu, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh, &
                                           my_neighbors_ext_mesh, &
                                           request_send_recv_gpu,ninterface_acoustic,inum_interfaces_acoustic)

! non-blocking MPI send

  ! sends data
  ! note: assembling data already filled into buffer_send_scalar_gpu array

  use constants

  implicit none

  integer :: NPROC,ninterface_acoustic
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_scalar_gpu,buffer_recv_scalar_gpu

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbors_ext_mesh
  integer, dimension(2*num_interfaces_ext_mesh) :: request_send_recv_gpu
integer, dimension(num_interfaces_ext_mesh), intent(in)  :: inum_interfaces_acoustic
  ! local parameters
  integer :: iinterface,num_interface

  request_send_recv_gpu(:) = 0


  ! sends only if more than one partition
  if (NPROC > 1) then

    ! note: partition border copy into the buffer has already been done
    !          by routine transfer_boun_pot_from_device()

    ! send messages
    do iinterface = 1, ninterface_acoustic
      num_interface = inum_interfaces_acoustic(iinterface)

      call isend_cr(buffer_send_scalar_gpu(1,num_interface), &
                    nibool_interfaces_ext_mesh(num_interface), &
                    my_neighbors_ext_mesh(num_interface), &
                    itag, &
                    request_send_recv_gpu(num_interface) )

      ! receive request
      call irecv_cr(buffer_recv_scalar_gpu(1,num_interface), &
                    nibool_interfaces_ext_mesh(num_interface), &
                    my_neighbors_ext_mesh(num_interface), &
                    itag, &
                    request_send_recv_gpu(num_interface+num_interfaces_ext_mesh) )
    enddo

  endif

  end subroutine assemble_MPI_scalar_send_cuda

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_write_cuda(NPROC, &
                        Mesh_pointer, &
                        buffer_recv_scalar_gpu,num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
                        request_send_recv_gpu, &
                        FORWARD_OR_ADJOINT,ninterface_acoustic,inum_interfaces_acoustic)

! waits for send/receiver to be completed and assembles contributions

  use constants

  implicit none

  integer :: NPROC
  integer :: ninterface_acoustic
  integer(kind=8) :: Mesh_pointer

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_scalar_gpu

  integer, dimension(2*num_interfaces_ext_mesh) :: request_send_recv_gpu
  integer, dimension(num_interfaces_ext_mesh), intent(in)  :: inum_interfaces_acoustic
  integer :: FORWARD_OR_ADJOINT

  integer :: iinterface,num_interface ! ipoin


  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, ninterface_acoustic
      num_interface = inum_interfaces_acoustic(iinterface)
      call wait_req(request_send_recv_gpu(num_interface + num_interfaces_ext_mesh))
    enddo

    ! adding contributions of neighbors
    call transfer_asmbl_pot_to_device(Mesh_pointer,buffer_recv_scalar_gpu,FORWARD_OR_ADJOINT)

    ! wait for communications completion (send)
    do iinterface = 1, ninterface_acoustic
      num_interface = inum_interfaces_acoustic(iinterface)
      call wait_req(request_send_recv_gpu(num_interface))
    enddo

  endif

  end subroutine assemble_MPI_scalar_write_cuda



!
!-------------------------------------------------------------------------------------------------
!


  subroutine assemble_MPI_vector_send_cuda(NPROC, &
                                          buffer_send_vector_gpu,buffer_recv_vector_gpu, &
                                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh, &
                                          my_neighbors_ext_mesh, &
                                          request_send_recv_vector_gpu,ninterface_elastic,inum_interfaces_elastic)

! sends data
! note: array to assemble already filled into buffer_send_vector_gpu array

  use constants

  implicit none

  integer :: NPROC

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,ninterface_elastic

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_vector_gpu,buffer_recv_vector_gpu

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbors_ext_mesh
  integer, dimension(2*num_interfaces_ext_mesh) :: request_send_recv_vector_gpu
  integer, dimension(num_interfaces_ext_mesh), intent(in)  :: inum_interfaces_elastic
  ! local parameters
  integer :: iinterface,num_interface

  ! note: preparation of the contribution between partitions using MPI
  !          already done in transfer_boun_accel routine

  ! send only if more than one partition
  if (NPROC > 1) then

    ! send messages
    do iinterface = 1, ninterface_elastic
      num_interface=inum_interfaces_elastic(iinterface)

      call isend_cr(buffer_send_vector_gpu(1,1,num_interface), &
                     NDIM*nibool_interfaces_ext_mesh(num_interface), &
                     my_neighbors_ext_mesh(num_interface), &
                     itag, &
                     request_send_recv_vector_gpu(num_interface) )

      call irecv_cr(buffer_recv_vector_gpu(1,1,num_interface), &
                     NDIM*nibool_interfaces_ext_mesh(num_interface), &
                     my_neighbors_ext_mesh(num_interface), &
                     itag, &
                     request_send_recv_vector_gpu(num_interface + num_interfaces_ext_mesh ) )
    enddo

  endif

  end subroutine assemble_MPI_vector_send_cuda


!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_write_cuda(NPROC,Mesh_pointer, &
                                            buffer_recv_vector_gpu, &
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                            nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                            request_send_recv_vector_gpu, &
                                            FORWARD_OR_ADJOINT,ninterface_elastic,inum_interfaces_elastic )

! waits for data to receive and assembles

  use constants

  implicit none

  integer :: NPROC
  integer(kind=8) :: Mesh_pointer



  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,ninterface_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_gpu

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(2*num_interfaces_ext_mesh) :: request_send_recv_vector_gpu
  integer, dimension(num_interfaces_ext_mesh), intent(in)  :: inum_interfaces_elastic
  integer :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface, num_interface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then


    ! wait for communications completion (recv)
    do iinterface = 1, ninterface_elastic
      num_interface=inum_interfaces_elastic(iinterface)
      call wait_req(request_send_recv_vector_gpu(num_interface + num_interfaces_ext_mesh))
    enddo

    ! adding contributions of neighbors
    call transfer_asmbl_accel_to_device(Mesh_pointer, &
                                        buffer_recv_vector_gpu, &
                                        max_nibool_interfaces_ext_mesh, &
                                        nibool_interfaces_ext_mesh, &
                                        ibool_interfaces_ext_mesh,FORWARD_OR_ADJOINT)

    ! This step is done via previous function transfer_and_assemble...
    ! do iinterface = 1, num_interfaces_ext_mesh
    !   do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
    !     array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
    !          array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_gpu(:,ipoin,iinterface)
    !   enddo
    ! enddo

    ! wait for communications completion (send)
    do iinterface = 1, ninterface_elastic
      num_interface=inum_interfaces_elastic(iinterface)
      call wait_req(request_send_recv_vector_gpu(num_interface))
    enddo

  endif

  end subroutine assemble_MPI_vector_write_cuda


!
!-------------------------------------------------------------------------------------------------
!

! with cuda functions...

  subroutine transfer_boundary_to_device(NPROC, Mesh_pointer, &
                                            buffer_recv_vector_gpu, &
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                            request_send_recv_vector_gpu,ninterface_elastic,inum_interfaces_elastic )

  use constants

  implicit none

  integer :: NPROC
  integer(kind=8) :: Mesh_pointer

  ! array to assemble
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh,ninterface_elastic

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_gpu

  integer, dimension(2*num_interfaces_ext_mesh) :: request_send_recv_vector_gpu
  integer, dimension(num_interfaces_ext_mesh), intent(in)  :: inum_interfaces_elastic
  ! local parameters
  integer :: iinterface,num_interface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! wait for communications completion (recv)
    !write(IMAIN,*) "sending MPI_wait"
    do iinterface = 1, ninterface_elastic
      num_interface=inum_interfaces_elastic(iinterface)
      call wait_req(request_send_recv_vector_gpu(num_interface + num_interfaces_ext_mesh))
    enddo

    ! send contributions to GPU
    call transfer_boundary_to_device_a(Mesh_pointer, buffer_recv_vector_gpu, max_nibool_interfaces_ext_mesh)
  endif

  ! This step is done via previous function transfer_and_assemble...
  ! do iinterface = 1, num_interfaces_ext_mesh
  !   do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
  !     array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
  !          array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_gpu(:,ipoin,iinterface)
  !   enddo
  ! enddo

  end subroutine transfer_boundary_to_device

