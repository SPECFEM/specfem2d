!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================


  subroutine determine_abs_surface()

! determines absorbing boundary elements

  use constants, only: IBOTTOM,IRIGHT,ITOP,ILEFT

  use part_unstruct_par, only: nelemabs,abs_surface,elmnts,nxread,nzread

  use shared_parameters, only: ngnod,absorbbottom,absorbleft,absorbright,absorbtop

  implicit none

  ! local parameters
  integer :: ix,iz
  integer :: inumelem

  !
  !--- definition of absorbing boundaries
  !
  nelemabs = 0
  if (absorbbottom) nelemabs = nelemabs + nxread
  if (absorbtop) nelemabs = nelemabs + nxread
  if (absorbleft) nelemabs = nelemabs + nzread
  if (absorbright) nelemabs = nelemabs + nzread

  allocate(abs_surface(5,nelemabs))

  ! generate the list of absorbing elements
  if (nelemabs > 0) then
    nelemabs = 0
    do iz = 1,nzread
       do ix = 1,nxread
          inumelem = (iz-1)*nxread + ix
          if (absorbbottom .and. iz == 1) then
             nelemabs = nelemabs + 1
             abs_surface(1,nelemabs) = inumelem-1
             abs_surface(2,nelemabs) = 2
             abs_surface(3,nelemabs) = elmnts(0+ngnod*(inumelem-1))
             abs_surface(4,nelemabs) = elmnts(1+ngnod*(inumelem-1))
             abs_surface(5,nelemabs) = IBOTTOM
             !is_abs_surf(inumelem) = .true.
          endif
          if (absorbright .and. ix == nxread) then
             nelemabs = nelemabs + 1
             abs_surface(1,nelemabs) = inumelem-1
             abs_surface(2,nelemabs) = 2
             abs_surface(3,nelemabs) = elmnts(1+ngnod*(inumelem-1))
             abs_surface(4,nelemabs) = elmnts(2+ngnod*(inumelem-1))
             abs_surface(5,nelemabs) = IRIGHT
             !is_abs_surf(inumelem) = .true.
          endif
          if (absorbtop .and. iz == nzread) then
             nelemabs = nelemabs + 1
             abs_surface(1,nelemabs) = inumelem-1
             abs_surface(2,nelemabs) = 2
             abs_surface(3,nelemabs) = elmnts(3+ngnod*(inumelem-1))
             abs_surface(4,nelemabs) = elmnts(2+ngnod*(inumelem-1))
             abs_surface(5,nelemabs) = ITOP
             !is_abs_surf(inumelem) = .true.
          endif
          if (absorbleft .and. ix == 1) then
             nelemabs = nelemabs + 1
             abs_surface(1,nelemabs) = inumelem-1
             abs_surface(2,nelemabs) = 2
             abs_surface(3,nelemabs) = elmnts(0+ngnod*(inumelem-1))
             abs_surface(4,nelemabs) = elmnts(3+ngnod*(inumelem-1))
             abs_surface(5,nelemabs) = ILEFT
             !is_abs_surf(inumelem) = .true.
          endif
       enddo
    enddo
  endif

  end subroutine determine_abs_surface

