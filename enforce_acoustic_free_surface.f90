
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!                         Dimitri Komatitsch
!                     University of Pau, France
!
!                          (c) April 2007
!
!========================================================================

  subroutine enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                potential_acoustic,ispecnum_acoustic_surface,iedgenum_acoustic_surface, &
                ibool,nelem_acoustic_surface,npoin,nspec)

! free surface for an acoustic medium
! if acoustic, the free surface condition is a Dirichlet condition for the potential,
! not Neumann, in order to impose zero pressure at the surface

  implicit none

  include "constants.h"

  integer :: nelem_acoustic_surface,npoin,nspec

  integer, dimension(nelem_acoustic_surface) :: ispecnum_acoustic_surface,iedgenum_acoustic_surface

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

  double precision, dimension(npoin) :: potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic

!---
!--- local variables
!---

  integer :: ispec_acoustic_surface,ispec,iedge,i,j,iglob

  do ispec_acoustic_surface = 1,nelem_acoustic_surface

    ispec = ispecnum_acoustic_surface(ispec_acoustic_surface)
    iedge = iedgenum_acoustic_surface(ispec_acoustic_surface)

    if(iedge == IBOTTOM .or. iedge == ITOP) then
      if(iedge == IBOTTOM) then
        j = 1
      else
        j = NGLLZ
      endif
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        potential_acoustic(iglob) = ZERO
        potential_dot_acoustic(iglob) = ZERO
        potential_dot_dot_acoustic(iglob) = ZERO
      enddo
    else
      if(iedge == ILEFT) then
        i = 1
      else
        i = NGLLX
      endif
      do j = 1,NGLLZ
        iglob = ibool(i,j,ispec)
        potential_acoustic(iglob) = ZERO
        potential_dot_acoustic(iglob) = ZERO
        potential_dot_dot_acoustic(iglob) = ZERO
      enddo
    endif

  enddo

  end subroutine enforce_acoustic_free_surface

