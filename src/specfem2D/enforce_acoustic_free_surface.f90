
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

  subroutine enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                                          potential_acoustic,acoustic_surface, &
                                          ibool,nelem_acoustic_surface,nglob,nspec)

! free surface for an acoustic medium
! if acoustic, the free surface condition is a Dirichlet condition for the potential,
! not Neumann, in order to impose zero pressure at the surface

  implicit none

  include "constants.h"

  integer :: nelem_acoustic_surface,nglob,nspec

  integer, dimension(5,nelem_acoustic_surface) :: acoustic_surface

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

  real(kind=CUSTOM_REAL), dimension(nglob) :: &
    potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic

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

