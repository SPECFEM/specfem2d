
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

!
! From array 'surface' (element, type : node/edge, node(s) ) that describes the
! acoustic free surface, determines the points (ixmin, ixmax, izmin and izmax) on the surface
! for each element.
! We chose to have ixmin <= ixmax and izmin <= izmax, so as to be able to have DO loops on it with
! an increment of +1.
!
subroutine construct_acoustic_surface ( nspec, ngnod, knods, nsurface, surface, tab_surface )

  implicit none

  integer, intent(in)  :: nspec
    integer, intent(in)  :: ngnod
  integer, dimension(ngnod,nspec), intent(in)  :: knods
  integer, intent(in)  :: nsurface
  integer, dimension(4,nsurface), intent(in)  :: surface
  integer, dimension(5,nsurface), intent(out)  :: tab_surface

  integer  :: i, k
  integer  :: ixmin, ixmax
  integer  :: izmin, izmax
  integer, dimension(ngnod)  :: n
  integer  :: e1, e2
  integer  :: type

  do i = 1, nsurface
     tab_surface(1,i) = surface(1,i)
     type = surface(2,i)
     e1 = surface(3,i)
     e2 = surface(4,i)
     do k = 1, ngnod
        n(k) = knods(k,tab_surface(1,i))
     enddo

     call get_acoustic_edge ( ngnod, n, type, e1, e2, ixmin, ixmax, izmin, izmax )

     tab_surface(2,i) = ixmin
     tab_surface(3,i) = ixmax
     tab_surface(4,i) = izmin
     tab_surface(5,i) = izmax

  enddo

end subroutine construct_acoustic_surface


!-----------------------------------------------
! Get the points (ixmin, ixmax, izmin and izmax) on an node/edge for one element.
!-----------------------------------------------
subroutine get_acoustic_edge ( ngnod, n, type, e1, e2, ixmin, ixmax, izmin, izmax )

  implicit none
  include "constants.h"

  integer, intent(in)  :: ngnod
  integer, dimension(ngnod), intent(in)  :: n
  integer, intent(in)  :: type, e1, e2
  integer, intent(out)  :: ixmin, ixmax, izmin, izmax


  if ( type == 1 ) then
     if ( e1 == n(1) ) then
        ixmin = 1
        ixmax = 1
        izmin = 1
        izmax = 1
     endif
     if ( e1 == n(2) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        izmin = 1
        izmax = 1
     endif
     if ( e1 == n(3) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        izmin = NGLLZ
        izmax = NGLLZ
     endif
     if ( e1 == n(4) ) then
        ixmin = 1
        ixmax = 1
        izmin = NGLLZ
        izmax = NGLLZ
     endif

  else
     if ( e1 ==  n(1) ) then
        ixmin = 1
        izmin = 1
        if ( e2 == n(2) ) then
           ixmax = NGLLX
           izmax = 1

        endif
        if ( e2 == n(4) ) then
           ixmax = 1
           izmax = NGLLZ

        endif
     endif
     if ( e1 == n(2) ) then
        ixmin = NGLLX
        izmin = 1
        if ( e2 == n(3) ) then
           ixmax = NGLLX
           izmax = NGLLZ

        endif
        if ( e2 == n(1) ) then
           ixmax = ixmin
           ixmin = 1
           izmax = 1

        endif
     endif
     if ( e1 == n(3) ) then
        ixmin = NGLLX
        izmin = NGLLZ
        if ( e2 == n(4) ) then
           ixmax = ixmin
           ixmin = 1
           izmax = NGLLZ

        endif
        if ( e2 == n(2) ) then
           ixmax = NGLLX
           izmax = izmin
           izmin = 1

        endif
     endif
     if ( e1 == n(4) ) then
        ixmin = 1
        izmin = NGLLZ
        if ( e2 == n(1) ) then
           ixmax = 1
           izmax = izmin
           izmin = 1

        endif
        if ( e2 == n(3) ) then
           ixmax = NGLLX
           izmax = NGLLZ

        endif
     endif
  endif

end subroutine get_acoustic_edge

