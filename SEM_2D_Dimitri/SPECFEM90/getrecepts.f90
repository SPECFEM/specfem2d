!=====================================================================
!
!                 S p e c f e m  V e r s i o n  4 . 2
!                 -----------------------------------
!
!                         Dimitri Komatitsch
!    Department of Earth and Planetary Sciences - Harvard University
!                         Jean-Pierre Vilotte
!                 Departement de Sismologie - IPGP - Paris
!                           (c) June 1998
!
!=====================================================================

  subroutine getrecepts(posrec,ndime,nrec)
!
!=======================================================================
!
!     "getrecepts" : lecture position recepteurs
!
!=======================================================================
!
  use iounit
  use infos

  implicit none

  integer ndime,nrec
  double precision posrec(ndime,nrec)

  double precision, dimension(:), allocatable :: posrecread

  integer i,j,irec
  character(len=80) datlin

!
!---- read receivers position
!
  irec = 0
  read(iin ,40) datlin
  allocate(posrecread(ndime))
  do i=1,nrec
   read(iin ,*) irec,(posrecread(j),j=1,ndime)
   if(irec<1 .or. irec>nrec) stop 'Wrong receiver number'
   posrec(:,irec) = posrecread
  enddo
  deallocate(posrecread)

  return

  40    format(a80)

  end subroutine getrecepts
