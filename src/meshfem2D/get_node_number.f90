
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


! *******************
! meshing subroutines
! *******************

!--- global node number

integer function num(i,j,nx)

  implicit none

  integer i,j,nx

  num = j*(nx+1) + i + 1

end function num


!---  global node number (when ngnod==4).
integer function num_4(i,j,nx)

  implicit none

  integer i,j,nx

  num_4 = j*(nx+1) + i + 1

end function num_4


!---  global node number (when ngnod==9).
integer function num_9(i,j,nx,nz)

  implicit none

  integer i,j,nx,nz


  if ( (mod(i,2) == 0) .and. (mod(j,2) == 0) ) then
     num_9 = j/2 * (nx+1) + i/2 + 1
  else
     if ( mod(j,2) == 0 ) then
        num_9 = (nx+1)*(nz+1) + j/2 * nx + ceiling(real(i)/real(2))
     else
        num_9 = (nx+1)*(nz+1) + nx*(nz+1) + floor(real(j)/real(2))*(nx*2+1) + i + 1

     endif
  endif

end function num_9
