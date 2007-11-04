
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!  Main authors: Dimitri Komatitsch, Nicolas Le Goff and Roland Martin
!                     University of Pau, France
!
!                         (c) November 2007
!
!========================================================================

  subroutine enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                potential_acoustic,acoustic_surface, &
                ibool,nelem_acoustic_surface,npoin,nspec)

! free surface for an acoustic medium
! if acoustic, the free surface condition is a Dirichlet condition for the potential,
! not Neumann, in order to impose zero pressure at the surface

  implicit none

  include "constants.h"

  integer :: nelem_acoustic_surface,npoin,nspec

  integer, dimension(5,nelem_acoustic_surface) :: acoustic_surface

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

  real(kind=CUSTOM_REAL), dimension(npoin) :: potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic

!---
!--- local variables
!---

  integer :: ispec_acoustic_surface,ispec,i,j,iglob

  do ispec_acoustic_surface = 1, nelem_acoustic_surface

     ispec = acoustic_surface(1,ispec_acoustic_surface)
     do j = acoustic_surface(4,ispec_acoustic_surface), acoustic_surface(5,ispec_acoustic_surface)
        do i = acoustic_surface(2,ispec_acoustic_surface), acoustic_surface(3,ispec_acoustic_surface)
           iglob = ibool(i,j,ispec)
           potential_acoustic(iglob) = ZERO
           potential_dot_acoustic(iglob) = ZERO
           potential_dot_dot_acoustic(iglob) = ZERO

        enddo
     enddo

  enddo

  end subroutine enforce_acoustic_free_surface

