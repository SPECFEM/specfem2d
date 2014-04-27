
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

  subroutine is_in_convex_quadrilateral(elmnt_coords, x_coord, z_coord, is_in)

  implicit none

  double precision, dimension(2,4)  :: elmnt_coords
  double precision, intent(in)  :: x_coord, z_coord
  logical, intent(out)  :: is_in

  real :: x1, x2, x3, x4, z1, z2, z3, z4
  real  :: normal1, normal2, normal3, normal4

  x1 = elmnt_coords(1,1)
  x2 = elmnt_coords(1,2)
  x3 = elmnt_coords(1,3)
  x4 = elmnt_coords(1,4)
  z1 = elmnt_coords(2,1)
  z2 = elmnt_coords(2,2)
  z3 = elmnt_coords(2,3)
  z4 = elmnt_coords(2,4)

  normal1 = (z_coord-z1) * (x2-x1) - (x_coord-x1) * (z2-z1)
  normal2 = (z_coord-z2) * (x3-x2) - (x_coord-x2) * (z3-z2)
  normal3 = (z_coord-z3) * (x4-x3) - (x_coord-x3) * (z4-z3)
  normal4 = (z_coord-z4) * (x1-x4) - (x_coord-x4) * (z1-z4)

  if ((normal1 < 0) .or. (normal2 < 0) .or. (normal3 < 0) .or. (normal4 < 0)) then
    is_in = .false.
  else
    is_in = .true.
  endif

  end subroutine is_in_convex_quadrilateral

