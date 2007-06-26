subroutine prepare_assemble_MPI (myrank,nspec,ibool, &
     knods, ngnod, & 
     npoin, elastic, &
     ninterface, max_interface_size, &
     my_neighbours, my_nelmnts_neighbours, my_interfaces, &
     ibool_interfaces_acoustic, ibool_interfaces_elastic, &
     nibool_interfaces_acoustic, nibool_interfaces_elastic, &
     inum_interfaces_acoustic, inum_interfaces_elastic, &
     ninterface_acoustic, ninterface_elastic &
     )


  implicit none
  include 'constants.h'
  

  integer, intent(in)  :: nspec, myrank, npoin, ngnod
  logical, dimension(nspec), intent(in)  :: elastic
  integer, dimension(ngnod,nspec), intent(in)  :: knods
  integer, dimension(NGLLX,NGLLZ,nspec), intent(in)  :: ibool

  !integer, dimension(nspec)  :: inner_to_glob_ispec
  !integer, dimension(nspec)  :: interface_to_glob_ispec

  !integer, intent(inout)  :: nspec_inner_known
  !integer, intent(inout)  :: nspec_interface_known


  integer  :: ninterface
  integer  :: max_interface_size
  integer, dimension(ninterface)  :: my_neighbours
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

  logical, dimension(nspec)  :: ispec_is_inner
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
  integer  :: ier

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


!!$  nspec_inner_known = 0
!!$  do ispec = 1, nspec
!!$     if ( ispec_is_inner(ispec) ) then
!!$        nspec_inner_known = nspec_inner_known + 1
!!$     end if
!!$  end do
 
  !allocate(inner_to_glob_ispec(nspec_inner_known))
  !allocate(interface_to_glob_ispec(nspec-nspec_inner_known))
  


!!$  nspec_inner_known = 0
!!$  nspec_interface_known = 0
!!$  do ispec = 1, nspec
!!$     if ( ispec_is_inner(ispec) ) then 
!!$        nspec_inner_known = nspec_inner_known + 1
!!$        inner_to_glob_ispec(nspec_inner_known) = ispec
!!$     else
!!$        nspec_interface_known = nspec_interface_known + 1
!!$        interface_to_glob_ispec(nspec_interface_known) = ispec
!!$     end if
!!$
!!$  end do

  
end subroutine prepare_assemble_MPI



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

subroutine create_MPI_requests_SEND_RECV_acoustic(myrank, &
     ninterface, ninterface_acoustic, &
     nibool_interfaces_acoustic, &
     my_neighbours, &
     max_ibool_interfaces_size_acoustic, &
     buffer_send_faces_vector_acoustic, &
     buffer_recv_faces_vector_acoustic, &
     tab_requests_send_recv_acoustic, &
     inum_interfaces_acoustic &
     )

  implicit none

  include 'constants.h'
  include 'mpif.h'
  

  integer, intent(in)  :: myrank
  integer, intent(in)  :: ninterface, ninterface_acoustic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_acoustic
  integer, intent(in)  :: max_ibool_interfaces_size_acoustic
  double precision, dimension(max_ibool_interfaces_size_acoustic,ninterface_acoustic), intent(in)  :: &
       buffer_send_faces_vector_acoustic
  double precision, dimension(max_ibool_interfaces_size_acoustic,ninterface_acoustic), intent(in)  :: &
       buffer_recv_faces_vector_acoustic
  integer, dimension(ninterface_acoustic*2), intent(inout)  :: tab_requests_send_recv_acoustic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_acoustic
  integer, dimension(ninterface), intent(in) :: my_neighbours


  integer  :: inum_interface,num_interface
  integer  :: ier
    

  do inum_interface = 1, ninterface_acoustic
     
     num_interface = inum_interfaces_acoustic(inum_interface)
     
        call MPI_Send_init ( buffer_send_faces_vector_acoustic(1,inum_interface), &
             nibool_interfaces_acoustic(num_interface), MPI_DOUBLE_PRECISION, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_acoustic(inum_interface), ier)
        call MPI_Recv_init ( buffer_recv_faces_vector_acoustic(1,inum_interface), &
             nibool_interfaces_acoustic(num_interface), MPI_DOUBLE_PRECISION, &
             my_neighbours(num_interface), 12, MPI_COMM_WORLD, &
             tab_requests_send_recv_acoustic(ninterface_acoustic+inum_interface), ier)
  end do


end subroutine create_MPI_requests_SEND_RECV_acoustic



subroutine create_MPI_requests_SEND_RECV_elastic(myrank, &
     ninterface, ninterface_elastic, &
     nibool_interfaces_elastic, &
     my_neighbours, &
     max_ibool_interfaces_size_elastic, &
     buffer_send_faces_vector_elastic, &
     buffer_recv_faces_vector_elastic, &
     tab_requests_send_recv_elastic, &
     inum_interfaces_elastic &
     )

  implicit none

  include 'constants.h'
  include 'mpif.h'
  

  integer, intent(in)  :: myrank
  integer, intent(in)  :: ninterface, ninterface_elastic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_elastic
  integer, intent(in)  :: max_ibool_interfaces_size_elastic
  double precision, dimension(max_ibool_interfaces_size_elastic,ninterface_elastic), intent(in)  :: &
       buffer_send_faces_vector_elastic
  double precision, dimension(max_ibool_interfaces_size_elastic,ninterface_elastic), intent(in)  :: &
       buffer_recv_faces_vector_elastic
  integer, dimension(ninterface_elastic*2), intent(inout)  :: tab_requests_send_recv_elastic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_elastic
  integer, dimension(ninterface), intent(in) :: my_neighbours

  integer  :: inum_interface,num_interface
  integer  :: ier
  
  

  do inum_interface = 1, ninterface_elastic
     
     num_interface = inum_interfaces_elastic(inum_interface)
     
        call MPI_Send_init ( buffer_send_faces_vector_elastic(1,inum_interface), &
             NDIM*nibool_interfaces_elastic(num_interface), MPI_DOUBLE_PRECISION, &
             my_neighbours(num_interface), 13, MPI_COMM_WORLD, &
             tab_requests_send_recv_elastic(inum_interface), ier)
        call MPI_Recv_init ( buffer_recv_faces_vector_elastic(1,inum_interface), &
             NDIM*nibool_interfaces_elastic(num_interface), MPI_DOUBLE_PRECISION, &
             my_neighbours(num_interface), 13, MPI_COMM_WORLD, &
             tab_requests_send_recv_elastic(ninterface_elastic+inum_interface), ier)
  end do


end subroutine create_MPI_requests_SEND_RECV_elastic



subroutine assemble_MPI_scalar(myrank,array_val1, array_val2,npoin, &
     ninterface, max_interface_size, max_ibool_interfaces_size_acoustic, max_ibool_interfaces_size_elastic, &
     ibool_interfaces_acoustic,ibool_interfaces_elastic, nibool_interfaces_acoustic,nibool_interfaces_elastic, my_neighbours)

  implicit none

  include 'constants.h'
  include 'mpif.h'


  ! array to assemble
  double precision, dimension(npoin), intent(inout) :: array_val1, array_val2

  integer, intent(in)  :: myrank
  integer, intent(in)  :: npoin
  integer, intent(in)  :: ninterface
  integer, intent(in)  :: max_interface_size
  integer, intent(in)  :: max_ibool_interfaces_size_acoustic, max_ibool_interfaces_size_elastic
  integer, dimension(NGLLX*max_interface_size,ninterface), intent(in)  :: &
       ibool_interfaces_acoustic,ibool_interfaces_elastic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_acoustic,nibool_interfaces_elastic
  integer, dimension(ninterface), intent(in)  :: my_neighbours


  double precision, dimension(max_ibool_interfaces_size_acoustic+max_ibool_interfaces_size_elastic, ninterface)  :: &
       buffer_send_faces_scalar, &
       buffer_recv_faces_scalar
  integer  :: msg_status(MPI_STATUS_SIZE)
  integer, dimension(ninterface)  :: msg_requests
  integer  :: ipoin, num_interface
  integer  :: ier
  
  integer  :: i


  do num_interface = 1, ninterface
   
     print *, 'QQQQQ', myrank,num_interface,nibool_interfaces_acoustic(num_interface), nibool_interfaces_elastic(num_interface)

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



subroutine assemble_MPI_vector_acoustic_start(myrank,array_val1,npoin, &
     ninterface, ninterface_acoustic, &
     inum_interfaces_acoustic, &
     max_interface_size, max_ibool_interfaces_size_acoustic,&
     ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
     tab_requests_send_recv_acoustic, &
     buffer_send_faces_vector_acoustic, &
     buffer_recv_faces_vector_acoustic &
     )

  implicit none

  include 'constants.h'
  include 'mpif.h'


  ! array to assemble
  double precision, dimension(npoin), intent(in) :: array_val1


  integer, intent(in)  :: myrank
  integer, intent(in)  :: npoin
  integer, intent(in)  :: ninterface, ninterface_acoustic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_acoustic
  integer, intent(in)  :: max_interface_size
  integer, intent(in)  :: max_ibool_interfaces_size_acoustic
  integer, dimension(NGLLX*max_interface_size,ninterface), intent(in)  :: ibool_interfaces_acoustic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_acoustic
  integer, dimension(ninterface_acoustic*2), intent(inout)  :: tab_requests_send_recv_acoustic
  double precision, dimension(max_ibool_interfaces_size_acoustic,ninterface_acoustic), intent(inout)  :: &
       buffer_send_faces_vector_acoustic, buffer_recv_faces_vector_acoustic

  integer  :: ipoin, num_interface, inum_interface
  integer  :: ier
  integer, dimension(MPI_STATUS_SIZE,ninterface_acoustic*2)  :: tab_statuses_acoustic

  integer  :: i


  do inum_interface = 1, ninterface_acoustic
     
     num_interface = inum_interfaces_acoustic(inum_interface)

     ipoin = 0
     do i = 1, nibool_interfaces_acoustic(num_interface)
     ipoin = ipoin + 1      
        buffer_send_faces_vector_acoustic(ipoin,inum_interface) = &
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


end subroutine assemble_MPI_vector_acoustic_start



subroutine assemble_MPI_vector_elastic_start(myrank,array_val2,npoin, &
     ninterface, ninterface_elastic, &
     inum_interfaces_elastic, &
     max_interface_size, max_ibool_interfaces_size_elastic,&
     ibool_interfaces_elastic, nibool_interfaces_elastic, &
     tab_requests_send_recv_elastic, &
     buffer_send_faces_vector_elastic, &
     buffer_recv_faces_vector_elastic &
     )

  implicit none

  include 'constants.h'
  include 'mpif.h'


  ! array to assemble
  double precision, dimension(NDIM,npoin), intent(in) :: array_val2


  integer, intent(in)  :: myrank
  integer, intent(in)  :: npoin
  integer, intent(in)  :: ninterface, ninterface_elastic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_elastic
  integer, intent(in)  :: max_interface_size
  integer, intent(in)  :: max_ibool_interfaces_size_elastic
  integer, dimension(NGLLX*max_interface_size,ninterface), intent(in)  :: ibool_interfaces_elastic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_elastic
  integer, dimension(ninterface_elastic*2), intent(inout)  :: tab_requests_send_recv_elastic
  double precision, dimension(max_ibool_interfaces_size_elastic,ninterface_elastic), intent(inout)  :: &
       buffer_send_faces_vector_elastic, buffer_recv_faces_vector_elastic
  


  integer  :: ipoin, num_interface, inum_interface
  integer  :: ier
  integer, dimension(MPI_STATUS_SIZE,ninterface_elastic)  :: tab_statuses_elastic

  integer  :: i


  do inum_interface = 1, ninterface_elastic
     
     num_interface = inum_interfaces_elastic(inum_interface)

     ipoin = 0
     do i = 1, nibool_interfaces_elastic(num_interface)
        buffer_send_faces_vector_elastic(ipoin+1:ipoin+2,inum_interface) = &
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


end subroutine assemble_MPI_vector_elastic_start



subroutine assemble_MPI_vector_acoustic_wait(myrank,array_val1,npoin, &
     ninterface, ninterface_acoustic, &
     inum_interfaces_acoustic, &
     max_interface_size, max_ibool_interfaces_size_acoustic,&
     ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
     tab_requests_send_recv_acoustic, &
     buffer_send_faces_vector_acoustic, &
     buffer_recv_faces_vector_acoustic &
     )

  implicit none

  include 'constants.h'
  include 'mpif.h'


  ! array to assemble
  double precision, dimension(npoin), intent(inout) :: array_val1


  integer, intent(in)  :: myrank
  integer, intent(in)  :: npoin
  integer, intent(in)  :: ninterface, ninterface_acoustic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_acoustic
  integer, intent(in)  :: max_interface_size
  integer, intent(in)  :: max_ibool_interfaces_size_acoustic
  integer, dimension(NGLLX*max_interface_size,ninterface), intent(in)  :: ibool_interfaces_acoustic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_acoustic
  integer, dimension(ninterface_acoustic*2), intent(inout)  :: tab_requests_send_recv_acoustic
  double precision, dimension(max_ibool_interfaces_size_acoustic,ninterface_acoustic), intent(inout)  :: &
       buffer_send_faces_vector_acoustic, buffer_recv_faces_vector_acoustic

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
             buffer_recv_faces_vector_acoustic(ipoin,inum_interface)
     end do
     
  end do


end subroutine assemble_MPI_vector_acoustic_wait



subroutine assemble_MPI_vector_elastic_wait(myrank,array_val2,npoin, &
     ninterface, ninterface_elastic, &
     inum_interfaces_elastic, &
     max_interface_size, max_ibool_interfaces_size_elastic,&
     ibool_interfaces_elastic, nibool_interfaces_elastic, &
     tab_requests_send_recv_elastic, &
     buffer_send_faces_vector_elastic, &
     buffer_recv_faces_vector_elastic &
     )

  implicit none

  include 'constants.h'
  include 'mpif.h'


  ! array to assemble
  double precision, dimension(NDIM,npoin), intent(inout) :: array_val2


  integer, intent(in)  :: myrank
  integer, intent(in)  :: npoin
  integer, intent(in)  :: ninterface, ninterface_elastic
  integer, dimension(ninterface), intent(in)  :: inum_interfaces_elastic
  integer, intent(in)  :: max_interface_size
  integer, intent(in)  :: max_ibool_interfaces_size_elastic
  integer, dimension(NGLLX*max_interface_size,ninterface), intent(in)  :: ibool_interfaces_elastic
  integer, dimension(ninterface), intent(in)  :: nibool_interfaces_elastic
  integer, dimension(ninterface_elastic*2), intent(inout)  :: tab_requests_send_recv_elastic
  double precision, dimension(max_ibool_interfaces_size_elastic,ninterface_elastic), intent(inout)  :: &
       buffer_send_faces_vector_elastic, buffer_recv_faces_vector_elastic

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
             buffer_recv_faces_vector_elastic(ipoin+1:ipoin+2,inum_interface)
        ipoin = ipoin + 2
     end do
     
  end do


end subroutine assemble_MPI_vector_elastic_wait

#endif



subroutine exit_MPI(error_msg)

  implicit none

#ifdef USE_MPI
  ! standard include of the MPI library
  include 'mpif.h'
#endif

  ! identifier for error message file
  integer, parameter :: IERROR = 30

  integer myrank
  character(len=*) error_msg

  integer ier
  character(len=80) outputname


  ! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI... proc '

  ! stop all the MPI processes, and exit
  ! on some machines, MPI_FINALIZE needs to be called before MPI_ABORT
#ifdef USE_MPI
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
  call MPI_FINALIZE(ier)
  
#endif
  stop 'error, program ended in exit_MPI'
      

end subroutine exit_MPI






