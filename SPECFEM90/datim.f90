
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) December 2004
!
!========================================================================

  subroutine datim(string_input)

! get date and time using f90 portable routines

  implicit none

  include "constants.h"

  character(len=50) string_input
  character(len=8) datein
  character(len=10) timein
  character(len=16) dateprint
  character(len=8) timeprint

  datein = ''
  timein = ''

  call date_and_time(datein,timein)

  dateprint = datein(7:8)//' - '//datein(5:6)//' - '//datein(1:4)
  timeprint = timein(1:2)//':'//timein(3:4)//':'//timein(5:6)

  write(iout,100)
  write(iout,101) string_input
  write(iout,102) dateprint,timeprint

!
!---- formats
!

 100 format(//1x,79('-')/1x,79('-')/1x,'Program SPECFEM2D: ')
 101 format(1x,79('-')/1x,79('-')/1x,a50)
 102 format(1x,79('-')/,1x,79('-')/' D a t e : ',a16,30x,' T i m e  : ',a8/1x,79('-'),/1x,79('-'))

  end subroutine datim

