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

  subroutine datim (string1,string2,iout)
!
!=======================================================================
!
!     D a t i m : Get date and time using f90 portable routines
!     ---------
!
!=======================================================================
!
  implicit none

  character(len=*) string1
  character(len=50) string2
  character(len=8) datein
  character(len=10)  timein
  character(len=16) dateprint
  character(len=8)  timeprint

  integer iout

!-----------------------------------------------------------------------

  datein = ''
  timein = ''

  call date_and_time(datein,timein)

  dateprint = datein(7:8)//' - '//datein(5:6)//' - '//datein(1:4)
  timeprint = timein(1:2)//':'//timein(3:4)//':'//timein(5:6)

!
!-------------------------------------------------------------------
!
   write(iout,100) string1
   write(iout,101) string2
   write(iout,102) dateprint,timeprint

  return
!
!---- formats
!

  100   format(//1x,79('-')/1x,79('-')/1x,a)
  101   format(1x,79('-')/1x,79('-')/1x,a50)
  102   format(1x,79('-')/,1x,79('-')/' D a t e : ',a16, &
         30x,' T i m e  : ',a8/1x,79('-'),/1x,79('-'))

  end subroutine datim
