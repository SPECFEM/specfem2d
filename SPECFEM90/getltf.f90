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

  subroutine getltf(gltfu,nltfl,initialfield)
!
!=======================================================================
!
!     "g e t l t f" : Read and store source functions
!
!=======================================================================
!
  use iounit
  use infos
  use defpi

  implicit none

  character(len=80) datlin
  character(len=20) funcname(10)

!
!-----------------------------------------------------------------------
!
  integer nltfl
  double precision gltfu(20,nltfl)
  logical initialfield

  integer n,isource,iexplo,k

  funcname(1) = ' '
  funcname(2) = ' '
  funcname(3) = ' '
  funcname(4) = 'Gaussian'
  funcname(5) = 'First derivative of a Gaussian'
  funcname(6) = 'Ricker'
  funcname(7) = 'Dirac'

!
!---- read load function parameters
!
  read(iin ,40) datlin

  do n = 1,nltfl
   read(iin ,*) (gltfu(k,n), k=1,9)
  enddo

!
!-----  check the input
!
 if(iecho  /=  0 .and. .not. initialfield) then
 do n = 1,nltfl
   if((nint(gltfu(1,n)) < 4).or.(nint(gltfu(1,n)) > 7)) &
        stop 'Wrong function number in getltf !'
   if(mod(n,50)  ==  1) write(iout,100) nltfl
   isource = nint(gltfu(1,n))
   iexplo = nint(gltfu(2,n))
   if (iexplo == 1) then
     write(iout,210) funcname(isource),(gltfu(k,n), k=3,8)
   else if(iexplo == 2) then
     write(iout,220) funcname(isource),(gltfu(k,n), k=3,7)
   else
     stop 'Unknown source type number !'
   endif
 enddo
 endif
!
!-----  convert angle from degrees to radians
!
  do n = 1,nltfl
   isource = nint(gltfu(1,n))
   iexplo = nint(gltfu(2,n))
   if(isource >= 4.and.isource <= 6.and.iexplo == 1) then
        gltfu(8,n) = gltfu(8,n) * pi / 180.d0
   endif
  enddo

  return
!
!---- formats
!
  40    format(a80)
  100   format(//,' S o u r c e  F u n c t i o n',/1x,28('='),//5x, &
         'Number of source functions. . . . . . . . (nltfl) =',i5)
  210   format(//,5x, &
  'Source Type. . . . . . . . . . . . . . = Collocated Force',/5x, &
     'Function Type. . . . . . . . . . . . . =',1x,a,/5x, &
     'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
     'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
     'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
     'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
     'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
     'Angle from vertical direction (deg). . =',1pe20.10,/5x)
  220   format(//,5x, &
     'Source Type. . . . . . . . . . . . . . = Explosion',/5x, &
     'Function Type. . . . . . . . . . . . . =',1x,a,/5x, &
     'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
     'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
     'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
     'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
     'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x)

  end subroutine getltf
