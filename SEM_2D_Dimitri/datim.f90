
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!  Main authors: Dimitri Komatitsch, Nicolas Le Goff and Roland Martin
!                 University of Pau and CNRS, France
!
!        (c) University of Pau and CNRS, France, November 2007
!
!========================================================================

  subroutine datim(string_input)

! get date and time

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

  write(iout,"(//1x,79('-')/1x,79('-')/1x,'Program SPECFEM2D: ')")
  write(iout,"(1x,79('-')/1x,79('-')/1x,a50)") string_input
  write(iout,"(1x,79('-')/,1x,79('-')/' D a t e : ',a16,30x,' T i m e  : ',a8/1x,79('-'),/1x,79('-'))") dateprint,timeprint

  end subroutine datim

