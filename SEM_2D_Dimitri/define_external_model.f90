
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

  subroutine define_external_model(x,y,iflag_element,rho,vp,vs)

  implicit none

  include "constants.h"

! user can modify this routine to assign any different external Earth model (rho, vp, vs)
! based on the x and y coordinates of that grid point and the flag of the region it belongs to

  integer, intent(in) :: iflag_element

  double precision, intent(in) :: x,y

  double precision, intent(out) :: rho,vp,vs

! dummy routine here, just to demonstrate how the model can be assigned
  if(iflag_element == 1 .or. x < 1700.d0 .or. y >= 2300.d0) then
    rho = 2000.d0
    vp = 3000.d0
    vs = vp / sqrt(3.d0)
  else
    rho = 2500.d0
    vp = 3600.d0
    vs = vp / 2.d0
  endif

  end subroutine define_external_model

