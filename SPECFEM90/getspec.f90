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

  subroutine getspec(coorg,npgeo,ndime)
!
!=======================================================================
!
!     "g e t s p e c" : Read spectral macroblocs nodal coordinates
!
!=======================================================================
!
  use iounit
  use infos
  use label1

  implicit none

  integer ndime,npgeo
  double precision coorg(ndime,npgeo)

  double precision, dimension(:), allocatable :: coorgread

  integer ip,ipoin,n,i,id
  character(len=80) datlin

!
!----
!
  ipoin = 0
  read(iin,40) datlin
  allocate(coorgread(ndime))
  do ip = 1,npgeo
   read(iin,*) ipoin,(coorgread(id),id =1,ndime)
   if(ipoin<1 .or. ipoin>npgeo) stop 'Wrong control point number'
   coorg(:,ipoin) = coorgread
  enddo
  deallocate(coorgread)

!
!----  check the input
!
  if(iecho == 2) then
   do n = 1,npgeo
      if(mod(n,50)  ==  1) write(iout,100) (labelc(i),i=1,ndime)
      write(iout,200) n, (coorg(i,n), i=1,ndime)
   enddo
  endif

  return
!
!---- formats
!
  40    format(a80)
  100   format(///' S p e c t r a l   c o n t r o l   p o i n t s'/1x, &
         45('=')///,4x,' node number ',10x,2(a5,12x))
  200   format(6x,i5,10x,3(1pe15.8,2x))
!
  end subroutine getspec
