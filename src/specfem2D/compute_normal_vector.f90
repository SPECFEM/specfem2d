
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
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
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


  subroutine compute_normal_vector( angle, n1_x, n2_x, n3_x, n4_x, n1_z, n2_z, n3_z, n4_z )

  implicit none

  include 'constants.h'

  double precision :: angle
  double precision :: n1_x, n2_x, n3_x, n4_x, n1_z, n2_z, n3_z, n4_z

  double precision  :: theta1, theta2, theta3
  double precision  :: costheta1, costheta2, costheta3

  if ( abs(n2_z - n1_z) < TINYVAL ) then
     costheta1 = 0
  else
     costheta1 = (n2_z - n1_z) / sqrt((n2_x - n1_x)**2 + (n2_z - n1_z)**2)
  endif
  if ( abs(n3_z - n2_z) < TINYVAL ) then
     costheta2 = 0
  else
     costheta2 = (n3_z - n2_z) / sqrt((n3_x - n2_x)**2 + (n3_z - n2_z)**2)
  endif
  if ( abs(n4_z - n3_z) < TINYVAL ) then
     costheta3 = 0
  else
    costheta3 = (n4_z - n3_z) / sqrt((n4_x - n3_x)**2 + (n4_z - n3_z)**2)
  endif

  theta1 = - sign(1.d0,n2_x - n1_x) * acos(costheta1)
  theta2 = - sign(1.d0,n3_x - n2_x) * acos(costheta2)
  theta3 = - sign(1.d0,n4_x - n3_x) * acos(costheta3)

  ! a sum is needed here because in the case of a source force vector
  ! users can give an angle with respect to the normal to the topography surface,
  ! in which case we must compute the normal to the topography
  ! and add it the existing rotation angle
  angle = angle + (theta1 + theta2 + theta3) / 3.d0 + PI/2.d0

  end subroutine compute_normal_vector


!
!-------------------------------------------------------------------------------------------------
!

  subroutine tri_quad(n, n1, nnodes)

  implicit none

  integer  :: n1, nnodes
  integer, dimension(4)  :: n


  n(2) = n1

  if ( n1 == 1 ) then
     n(1) = nnodes
  else
     n(1) = n1-1
  endif

  if ( n1 == nnodes ) then
     n(3) = 1
  else
     n(3) = n1+1
  endif

  if ( n(3) == nnodes ) then
     n(4) = 1
  else
     n(4) = n(3)+1
  endif


  end subroutine tri_quad

