
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
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

! From array 'surface' (element, type : node/edge, node(s) ) that describes the
! acoustic free surface, determines the points (ixmin, ixmax, izmin and izmax) on the surface
! for each element.
! We chose to have ixmin <= ixmax and izmin <= izmax, so as to be able to have DO loops on it with
! an increment of +1.
!
subroutine construct_acoustic_surface ()

  use specfem_par, only : ngnod, knods, nelem_acoustic_surface, acoustic_edges, acoustic_surface

  implicit none

  integer  :: i, k
  integer  :: ixmin, ixmax
  integer  :: izmin, izmax
  integer, dimension(ngnod)  :: n
  integer  :: e1, e2
  integer  :: type_acoust

  do i = 1, nelem_acoustic_surface
     acoustic_surface(1,i) = acoustic_edges(1,i) ! Here we do a copy
     type_acoust = acoustic_edges(2,i)
     e1 = acoustic_edges(3,i)
     e2 = acoustic_edges(4,i)
     do k = 1, ngnod
        n(k) = knods(k,acoustic_surface(1,i))
     enddo

     call get_acoustic_edge ( ngnod, n, type_acoust, e1, e2, ixmin, ixmax, izmin, izmax )

     acoustic_surface(2,i) = ixmin
     acoustic_surface(3,i) = ixmax
     acoustic_surface(4,i) = izmin
     acoustic_surface(5,i) = izmax

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

