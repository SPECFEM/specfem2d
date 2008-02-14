
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic wave equation
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
! This file contains subroutines related to assembling (of the mass matrix, potential_dot_dot and
! accel_elastic).
! These subroutines are for the most part not used in the sequential version.
!

!-----------------------------------------------
! Determines the points that are on the interfaces with other partitions, to help
! build the communication buffers, and determines which elements are considered 'inner'
! (no points in common with other partitions) and 'outer' (at least one point in common
! with neighbouring partitions).
! We have both acoustic and elastic buffers, for coupling between acoustic and elastic elements
! led us to have two sets of communications.
!-----------------------------------------------
subroutine prepare_assemble_MPI (nspec,ibool, &
     knods, ngnod, &
     npoin, elastic, &
     ninterface, max_interface_size, &
     my_nelmnts_neighbours, my_interfaces, &
     ibool_interfaces_acoustic, ibool_interfaces_elastic, &
     nibool_interfaces_acoustic, nibool_interfaces_elastic, &
     inum_interfaces_acoustic, inum_interfaces_elastic, &
     ninterface_acoustic, ninterface_elastic, &
     mask_ispec_inner_outer &
     )

  implicit none

  include 'constants.h'

  integer, intent(in)  :: nspec, npoin, ngnod
  logical, dimension(nspec), intent(in)  :: elastic
  integer, dimension(ngnod,nspec), intent(in)  :: knods
  integer, dimension(NGLLX,NGLLZ,nspec), intent(in)  :: ibool

  integer  :: ninterface
  integer  :: max_interface_size
  integer, dimension(ninterface)  :: my_nelmnts_neighbours
  integer, dimension(4,max_interface_size,ninterface)  :: my_interfaces
  integer, dimension(NGLLX*max_interface_size,ninterface)  :: &
       ibool_interfaces_acoustic,ibool_interfaces_elastic
  integer, dimension(ninterface)  :: &
       nibool_interfaces_acoustic,nibool_interfaces_elastic

  integer, dimension(ninterface), intent(out)  :: &
       inum_interfaces_acoustic, inum_interfaces_elastic
  integer, intent(out)  :: ninterface_acoustic, ninterface_elastic

  integer  :: num_interface
  integer  :: ispec_interface

  logical, dimension(nspec), intent(inout)  :: mask_ispec_inner_outer

  logical, dimension(npoin)  :: mask_ibool_acoustic
  logical, dimension(npoin)  :: mask_ibool_elastic

  integer  :: ixmin, ixmax
  integer  :: izmin, izmax
  integer, dimension(ngnod)  :: n
  integer  :: e1, e2
  integer  :: type
  integer  :: ispec

  integer  :: k
  integer  :: npoin_interface_acoustic
  integer  :: npoin_interface_elastic

  integer  :: ix,iz

  integer  :: sens

  ibool_interfaces_acoustic(:,:) = 0
  nibool_interfaces_acoustic(:) = 0
  ibool_interfaces_elastic(:,:) = 0
  nibool_interfaces_elastic(:) = 0

  do num_interface = 1, ninterface
     npoin_interface_acoustic = 0
     npoin_interface_elastic = 0
     mask_ibool_acoustic(:) = .false.
     mask_ibool_elastic(:) = .false.

     do ispec_interface = 1, my_nelmnts_neighbours(num_interface)
        ispec = my_interfaces(1,ispec_interface,num_interface)
        type = my_interfaces(2,ispec_interface,num_interface)
        do k = 1, ngnod
           n(k) = knods(k,ispec)
        end do
        e1 = my_interfaces(3,ispec_interface,num_interface)
        e2 = my_interfaces(4,ispec_interface,num_interface)

        call get_edge(ngnod, n, type, e1, e2, ixmin, ixmax, izmin, izmax, sens)

        do iz = izmin, izmax, sens
           do ix = ixmin, ixmax, sens

              if ( elastic(ispec) ) then

                 if(.not. mask_ibool_elastic(ibool(ix,iz,ispec))) then
                    mask_ibool_elastic(ibool(ix,iz,ispec)) = .true.
                    npoin_interface_elastic = npoin_interface_elastic + 1
                    ibool_interfaces_elastic(npoin_interface_elastic,num_interface)=&
                         ibool(ix,iz,ispec)
                 end if
              else
                 if(.not. mask_ibool_acoustic(ibool(ix,iz,ispec))) then
                    mask_ibool_acoustic(ibool(ix,iz,ispec)) = .true.
                    npoin_interface_acoustic = npoin_interface_acoustic + 1
                    ibool_interfaces_acoustic(npoin_interface_acoustic,num_interface)=&
                         ibool(ix,iz,ispec)
                 end if
              end if
           end do
        end do

     end do
     nibool_interfaces_acoustic(num_interface) = npoin_interface_acoustic
     nibool_interfaces_elastic(num_interface) = npoin_interface_elastic

     do ispec = 1, nspec
       do iz = 1, NGLLZ
         do ix = 1, NGLLX
           if ( mask_ibool_acoustic(ibool(ix,iz,ispec)) &
            .or. mask_ibool_elastic(ibool(ix,iz,ispec)) ) then
               mask_ispec_inner_outer(ispec) = .true.
            endif

          enddo
        enddo
      enddo

  end do

  ninterface_acoustic = 0
  ninterface_elastic =  0
  do num_interface = 1, ninterface
     if ( nibool_interfaces_acoustic(num_interface) > 0 ) then
        ninterface_acoustic = ninterface_acoustic + 1
        inum_interfaces_acoustic(ninterface_acoustic) = num_interface
     end if
     if ( nibool_interfaces_elastic(num_interface) > 0 ) then
        ninterface_elastic = ninterface_elastic + 1
        inum_interfaces_elastic(ninterface_elastic) = num_interface
     end if
  end do

end subroutine prepare_assemble_MPI


!-----------------------------------------------
! Get the points (ixmin, ixmax, izmin and izmax) on an node/edge for one element.
! 'sens' is used to have DO loops with increment equal to 'sens' (-/+1).
!-----------------------------------------------
subroutine get_edge ( ngnod, n, type, e1, e2, ixmin, ixmax, izmin, izmax, sens )

  implicit none

  include "constants.h"

  integer, intent(in)  :: ngnod
  integer, dimension(ngnod), intent(in)  :: n
  integer, intent(in)  :: type, e1, e2
  integer, intent(out)  :: ixmin, ixmax, izmin, izmax
  integer, intent(out)  :: sens

   if ( type == 1 ) then
     if ( e1 == n(1) ) then
        ixmin = 1
        ixmax = 1
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(2) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        izmin = 1
        izmax = 1
     end if
     if ( e1 == n(3) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        izmin = NGLLZ
        izmax = NGLLZ
     end if
     if ( e1 == n(4) ) then
        ixmin = 1
        ixmax = 1
        izmin = NGLLZ
        izmax = NGLLZ
     end if
     sens = 1
  else
     if ( e1 ==  n(1) ) then
        ixmin = 1
        izmin = 1
        if ( e2 == n(2) ) then
           ixmax = NGLLX
           izmax = 1
           sens = 1
        end if
        if ( e2 == n(4) ) then
           ixmax = 1
           izmax = NGLLZ
           sens = 1
        end if
     end if
     if ( e1 == n(2) ) then
        ixmin = NGLLX
        izmin = 1
        if ( e2 == n(3) ) then
           ixmax = NGLLX
           izmax = NGLLZ
           sens = 1
        end if
        if ( e2 == n(1) ) then
           ixmax = 1
           izmax = 1
           sens = -1
        end if
     end if
     if ( e1 == n(3) ) then
        ixmin = NGLLX
        izmin = NGLLZ
        if ( e2 == n(4) ) then
           ixmax = 1
           izmax = NGLLZ
           sens = -1
        end if
        if ( e2 == n(2) ) then
           ixmax = NGLLX
           izmax = 1
           sens = -1
        end if
     end if
     if ( e1 == n(4) ) then
        ixmin = 1
        izmin = NGLLZ
        if ( e2 == n(1) ) then
           ixmax = 1
           izmax = 1
           sens = -1
        end if
        if ( e2 == n(3) ) then
           ixmax = NGLLX
           izmax = NGLLZ
           sens = 1
        end if
     end if
  end if

end subroutine get_edge


#ifdef USE_MPI


!-----------------------------------------------
! Creation of persistent communication requests (send and recv) for acoustic elements.
! Should be disposed of if using Paraver (with MPItrace), since it does not instrument persistent
! communications yet.
!-----------------------------------------------
subroutine create_MPI_req_SEND_RECV_ac( &
     ninterface, ninterface_acoustic, &
     nibool_interfaces_acoustic, &
     my_neighbours, &
     max_ibool_interfaces_size_ac, &
     buffer_send_faces_vector_ac, &
     buffer_recv_faces_vector_ac, &
     tab_requests_send_recv_acoustic, &
     inum_interfaces_acoustic &
     )

  implicit none

  include 'constants.h'
  include 'mpif.h'
  include 'precision_mpi.h'

  integer, intent(in)  :: ninterface, ninterface_acoustic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_acoustic
  integer, intent(in)  :: max_ibool_interfaces_size_ac
  real(kind=CUSTOM_REAL), dimension(max_ibool_interfaces_size_ac,ninterface_acoustic), intent(in)  :: &
       buffer_send_faces_vector_ac
  real(kind=CUSTOM_REAL), dimension(max_ibool_interfaces_size_ac,ninterface_acoustic), intent(in)  :: &
       buffer_recv_faces_vector_ac
  integer, dimension(ninterface_acoustic*2), intent(inout)  :: tab_requests_send_recv_acoustic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_acoustic
  integer, dimension(ninterface), intent(in) :: my_neighbours

  integer  :: inum_interface,num_interface
  integer  :: ier

  do inum_interface = 1, ninterface_acoustic

     num_interface = inum_interfaces_acoustic(inum_interface)

        call MPI_Send_init ( buffer_send_faces_vector_ac(1,inum_interface), &
             nibool_interfaces_acoustic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_acoustic(inum_interface), ier)
        call MPI_Recv_init ( buffer_recv_faces_vector_ac(1,inum_interface), &
             nibool_interfaces_acoustic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_acoustic(ninterface_acoustic+inum_interface), ier)
  end do

end subroutine create_MPI_req_SEND_RECV_ac


!-----------------------------------------------
! Creation of persistent communication requests (send and recv) for elastic elements.
! Should be disposed of if using Paraver (with MPItrace), since it does not instrument persistent
! communications yet.
!-----------------------------------------------
subroutine create_MPI_req_SEND_RECV_el( &
     ninterface, ninterface_elastic, &
     nibool_interfaces_elastic, &
     my_neighbours, &
     max_ibool_interfaces_size_el, &
     buffer_send_faces_vector_el, &
     buffer_recv_faces_vector_el, &
     tab_requests_send_recv_elastic, &
     inum_interfaces_elastic &
     )

  implicit none

  include 'constants.h'
  include 'mpif.h'
  include 'precision_mpi.h'


  integer, intent(in)  :: ninterface, ninterface_elastic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_elastic
  integer, intent(in)  :: max_ibool_interfaces_size_el
  real(kind=CUSTOM_REAL), dimension(max_ibool_interfaces_size_el,ninterface_elastic), intent(in)  :: &
       buffer_send_faces_vector_el
  real(kind=CUSTOM_REAL), dimension(max_ibool_interfaces_size_el,ninterface_elastic), intent(in)  :: &
       buffer_recv_faces_vector_el
  integer, dimension(ninterface_elastic*2), intent(inout)  :: tab_requests_send_recv_elastic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_elastic
  integer, dimension(ninterface), intent(in) :: my_neighbours

  integer  :: inum_interface,num_interface
  integer  :: ier

  do inum_interface = 1, ninterface_elastic

     num_interface = inum_interfaces_elastic(inum_interface)

        call MPI_Send_init ( buffer_send_faces_vector_el(1,inum_interface), &
             NDIM*nibool_interfaces_elastic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 13, MPI_COMM_WORLD, &
             tab_requests_send_recv_elastic(inum_interface), ier)
        call MPI_Recv_init ( buffer_recv_faces_vector_el(1,inum_interface), &
             NDIM*nibool_interfaces_elastic(num_interface), CUSTOM_MPI_TYPE, &
             my_neighbours(num_interface), 13, MPI_COMM_WORLD, &
             tab_requests_send_recv_elastic(ninterface_elastic+inum_interface), ier)
  end do

end subroutine create_MPI_req_SEND_RECV_el


!-----------------------------------------------
! Assembling the mass matrix.
!-----------------------------------------------
subroutine assemble_MPI_scalar(array_val1, array_val2,npoin, &
     ninterface, max_interface_size, max_ibool_interfaces_size_ac, max_ibool_interfaces_size_el, &
     ibool_interfaces_acoustic,ibool_interfaces_elastic, nibool_interfaces_acoustic,nibool_interfaces_elastic, my_neighbours)

  implicit none

  include 'constants.h'
  include 'mpif.h'

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(npoin), intent(inout) :: array_val1, array_val2

  integer, intent(in)  :: npoin
  integer, intent(in)  :: ninterface
  integer, intent(in)  :: max_interface_size
  integer, intent(in)  :: max_ibool_interfaces_size_ac, max_ibool_interfaces_size_el
  integer, dimension(NGLLX*max_interface_size,ninterface), intent(in)  :: &
       ibool_interfaces_acoustic,ibool_interfaces_elastic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_acoustic,nibool_interfaces_elastic
  integer, dimension(ninterface), intent(in)  :: my_neighbours

  double precision, dimension(max_ibool_interfaces_size_ac+max_ibool_interfaces_size_el, ninterface)  :: &
       buffer_send_faces_scalar, &
       buffer_recv_faces_scalar
  integer  :: msg_status(MPI_STATUS_SIZE)
  integer, dimension(ninterface)  :: msg_requests
  integer  :: ipoin, num_interface
  integer  :: ier

  integer  :: i

  do num_interface = 1, ninterface

     ipoin = 0
     do i = 1, nibool_interfaces_acoustic(num_interface)
        ipoin = ipoin + 1
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val1(ibool_interfaces_acoustic(i,num_interface))
     end do

     do i = 1, nibool_interfaces_elastic(num_interface)
        ipoin = ipoin + 1
        buffer_send_faces_scalar(ipoin,num_interface) = &
             array_val2(ibool_interfaces_elastic(i,num_interface))
     end do

     call MPI_isend ( buffer_send_faces_scalar(1,num_interface), &
          nibool_interfaces_acoustic(num_interface)+nibool_interfaces_elastic(num_interface), MPI_DOUBLE_PRECISION, &
          my_neighbours(num_interface), 11, &
          MPI_COMM_WORLD, msg_requests(num_interface), ier)

  end do

  do num_interface = 1, ninterface
     call MPI_recv ( buffer_recv_faces_scalar(1,num_interface), &
          nibool_interfaces_acoustic(num_interface)+nibool_interfaces_elastic(num_interface), MPI_DOUBLE_PRECISION, &
          my_neighbours(num_interface), 11, &
          MPI_COMM_WORLD, msg_status(1), ier)

     ipoin = 0
     do i = 1, nibool_interfaces_acoustic(num_interface)
        ipoin = ipoin + 1
        array_val1(ibool_interfaces_acoustic(i,num_interface)) = array_val1(ibool_interfaces_acoustic(i,num_interface)) + &
             buffer_recv_faces_scalar(ipoin,num_interface)
     end do

     do i = 1, nibool_interfaces_elastic(num_interface)
        ipoin = ipoin + 1
        array_val2(ibool_interfaces_elastic(i,num_interface)) = array_val2(ibool_interfaces_elastic(i,num_interface)) + &
             buffer_recv_faces_scalar(ipoin,num_interface)
     end do

  end do

  call MPI_BARRIER(mpi_comm_world,ier)

end subroutine assemble_MPI_scalar


!-----------------------------------------------
! Assembling potential_dot_dot for acoustic elements :
! the buffers are filled, and the send and recv are started here.
! We use MPI_Start (MPI_Startall is not used, since it causes problems in OpenMPI prior to v1.2).
!-----------------------------------------------
subroutine assemble_MPI_vector_ac_start(array_val1,npoin, &
     ninterface, ninterface_acoustic, &
     inum_interfaces_acoustic, &
     max_interface_size, max_ibool_interfaces_size_ac,&
     ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
     tab_requests_send_recv_acoustic, &
     buffer_send_faces_vector_ac &
     )

  implicit none

  include 'constants.h'
  include 'mpif.h'

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(npoin), intent(in) :: array_val1

  integer, intent(in)  :: npoin
  integer, intent(in)  :: ninterface, ninterface_acoustic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_acoustic
  integer, intent(in)  :: max_interface_size
  integer, intent(in)  :: max_ibool_interfaces_size_ac
  integer, dimension(NGLLX*max_interface_size,ninterface), intent(in)  :: ibool_interfaces_acoustic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_acoustic
  integer, dimension(ninterface_acoustic*2), intent(inout)  :: tab_requests_send_recv_acoustic
  real(kind=CUSTOM_REAL), dimension(max_ibool_interfaces_size_ac,ninterface_acoustic), intent(inout)  :: &
       buffer_send_faces_vector_ac

  integer  :: ipoin, num_interface, inum_interface
  integer  :: ier

  integer  :: i

  do inum_interface = 1, ninterface_acoustic

     num_interface = inum_interfaces_acoustic(inum_interface)

     ipoin = 0
     do i = 1, nibool_interfaces_acoustic(num_interface)
     ipoin = ipoin + 1
        buffer_send_faces_vector_ac(ipoin,inum_interface) = &
             array_val1(ibool_interfaces_acoustic(i,num_interface))
     end do

  end do

  do inum_interface = 1, ninterface_acoustic*2
     call MPI_START(tab_requests_send_recv_acoustic(inum_interface), ier)
     if ( ier /= MPI_SUCCESS ) then
        call exit_mpi('MPI_start unsuccessful in assemble_MPI_vector_start')
     end if
  end do

!call MPI_Startall ( ninterface*2, tab_requests_send_recv(1), ier )

end subroutine assemble_MPI_vector_ac_start


!-----------------------------------------------
! Assembling accel_elastic for elastic elements :
! the buffers are filled, and the send and recv are started here.
! We use MPI_Start (MPI_Startall is not used, since it causes problems in OpenMPI prior to v1.2).
!-----------------------------------------------
subroutine assemble_MPI_vector_el_start(array_val2,npoin, &
     ninterface, ninterface_elastic, &
     inum_interfaces_elastic, &
     max_interface_size, max_ibool_interfaces_size_el,&
     ibool_interfaces_elastic, nibool_interfaces_elastic, &
     tab_requests_send_recv_elastic, &
     buffer_send_faces_vector_el &
     )

  implicit none

  include 'constants.h'
  include 'mpif.h'

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin), intent(in) :: array_val2

  integer, intent(in)  :: npoin
  integer, intent(in)  :: ninterface, ninterface_elastic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_elastic
  integer, intent(in)  :: max_interface_size
  integer, intent(in)  :: max_ibool_interfaces_size_el
  integer, dimension(NGLLX*max_interface_size,ninterface), intent(in)  :: ibool_interfaces_elastic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_elastic
  integer, dimension(ninterface_elastic*2), intent(inout)  :: tab_requests_send_recv_elastic
  real(CUSTOM_REAL), dimension(max_ibool_interfaces_size_el,ninterface_elastic), intent(inout)  :: &
       buffer_send_faces_vector_el


  integer  :: ipoin, num_interface, inum_interface
  integer  :: ier

  integer  :: i


  do inum_interface = 1, ninterface_elastic

     num_interface = inum_interfaces_elastic(inum_interface)

     ipoin = 0
     do i = 1, nibool_interfaces_elastic(num_interface)
        buffer_send_faces_vector_el(ipoin+1:ipoin+2,inum_interface) = &
             array_val2(:,ibool_interfaces_elastic(i,num_interface))
        ipoin = ipoin + 2
     end do

  end do

  do inum_interface = 1, ninterface_elastic*2
     call MPI_START(tab_requests_send_recv_elastic(inum_interface), ier)
     if ( ier /= MPI_SUCCESS ) then
        call exit_mpi('MPI_start unsuccessful in assemble_MPI_vector_start')
     end if
  end do

!call MPI_Startall ( ninterface*2, tab_requests_send_recv(1), ier )

end subroutine assemble_MPI_vector_el_start


!-----------------------------------------------
! Assembling potential_dot_dot for acoustic elements :
! We wait for the completion of the communications, and add the contributions received
! for the points on the interfaces.
!-----------------------------------------------
subroutine assemble_MPI_vector_ac_wait(array_val1,npoin, &
     ninterface, ninterface_acoustic, &
     inum_interfaces_acoustic, &
     max_interface_size, max_ibool_interfaces_size_ac,&
     ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
     tab_requests_send_recv_acoustic, &
     buffer_recv_faces_vector_ac &
     )

  implicit none

  include 'constants.h'
  include 'mpif.h'

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(npoin), intent(inout) :: array_val1

  integer, intent(in)  :: npoin
  integer, intent(in)  :: ninterface, ninterface_acoustic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_acoustic
  integer, intent(in)  :: max_interface_size
  integer, intent(in)  :: max_ibool_interfaces_size_ac
  integer, dimension(NGLLX*max_interface_size,ninterface), intent(in)  :: ibool_interfaces_acoustic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_acoustic
  integer, dimension(ninterface_acoustic*2), intent(inout)  :: tab_requests_send_recv_acoustic
  real(kind=CUSTOM_REAL), dimension(max_ibool_interfaces_size_ac,ninterface_acoustic), intent(inout)  :: &
       buffer_recv_faces_vector_ac

  integer  :: ipoin, num_interface, inum_interface
  integer  :: ier
  integer, dimension(MPI_STATUS_SIZE,ninterface_acoustic*2)  :: tab_statuses_acoustic

  integer  :: i

  call MPI_Waitall ( ninterface_acoustic*2, tab_requests_send_recv_acoustic(1), tab_statuses_acoustic(1,1), ier )
  if ( ier /= MPI_SUCCESS ) then
     call exit_mpi('MPI_WAITALL unsuccessful in assemble_MPI_vector_wait')
  end if

  do inum_interface = 1, ninterface_acoustic

     num_interface = inum_interfaces_acoustic(inum_interface)

     ipoin = 0
     do i = 1, nibool_interfaces_acoustic(num_interface)
        ipoin = ipoin + 1
        array_val1(ibool_interfaces_acoustic(i,num_interface)) = array_val1(ibool_interfaces_acoustic(i,num_interface)) + &
             buffer_recv_faces_vector_ac(ipoin,inum_interface)
     end do

  end do

end subroutine assemble_MPI_vector_ac_wait


!-----------------------------------------------
! Assembling accel_elastic for elastic elements :
! We wait for the completion of the communications, and add the contributions received
! for the points on the interfaces.
!-----------------------------------------------
subroutine assemble_MPI_vector_el_wait(array_val2,npoin, &
     ninterface, ninterface_elastic, &
     inum_interfaces_elastic, &
     max_interface_size, max_ibool_interfaces_size_el,&
     ibool_interfaces_elastic, nibool_interfaces_elastic, &
     tab_requests_send_recv_elastic, &
     buffer_recv_faces_vector_el &
     )

  implicit none

  include 'constants.h'
  include 'mpif.h'

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin), intent(inout) :: array_val2

  integer, intent(in)  :: npoin
  integer, intent(in)  :: ninterface, ninterface_elastic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_elastic
  integer, intent(in)  :: max_interface_size
  integer, intent(in)  :: max_ibool_interfaces_size_el
  integer, dimension(NGLLX*max_interface_size,ninterface), intent(in)  :: ibool_interfaces_elastic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_elastic
  integer, dimension(ninterface_elastic*2), intent(inout)  :: tab_requests_send_recv_elastic
  real(kind=CUSTOM_REAL), dimension(max_ibool_interfaces_size_el,ninterface_elastic), intent(inout)  :: &
       buffer_recv_faces_vector_el

  integer  :: ipoin, num_interface, inum_interface
  integer  :: ier
  integer, dimension(MPI_STATUS_SIZE,ninterface_elastic*2)  :: tab_statuses_elastic

  integer  :: i

  call MPI_Waitall ( ninterface_elastic*2, tab_requests_send_recv_elastic(1), tab_statuses_elastic(1,1), ier )
  if ( ier /= MPI_SUCCESS ) then
     call exit_mpi('MPI_WAITALL unsuccessful in assemble_MPI_vector_wait')
  end if

  do inum_interface = 1, ninterface_elastic

     num_interface = inum_interfaces_elastic(inum_interface)

     ipoin = 0
     do i = 1, nibool_interfaces_elastic(num_interface)
        array_val2(:,ibool_interfaces_elastic(i,num_interface)) = array_val2(:,ibool_interfaces_elastic(i,num_interface)) + &
             buffer_recv_faces_vector_el(ipoin+1:ipoin+2,inum_interface)
        ipoin = ipoin + 2
     end do

  end do

end subroutine assemble_MPI_vector_el_wait

#endif


!-----------------------------------------------
! Dummy subroutine, to be able to stop the code whether sequential or parallel.
!-----------------------------------------------
subroutine exit_MPI(error_msg)

  implicit none

#ifdef USE_MPI
  ! standard include of the MPI library
  include 'mpif.h'
#endif

  ! identifier for error message file
  integer, parameter :: IERROR = 30

  character(len=*) error_msg

  integer ier

  ier = 0

  ! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI... proc '

  ! stop all the MPI processes, and exit
#ifdef USE_MPI
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
#endif

  stop 'error, program ended in exit_MPI'

end subroutine exit_MPI

