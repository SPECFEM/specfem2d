!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently maNZ_IMAGE_color more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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

  subroutine write_postscript_snapshot()

  use constants, only: IMAIN

  use specfem_par, only: myrank,P_SV,it, &
                         potential_acoustic,displ_elastic,displs_poroelastic, &
                         potential_dot_acoustic,veloc_elastic,velocs_poroelastic, &
                         potential_dot_dot_acoustic,accel_elastic,accels_poroelastic

  use specfem_par_movie, only: imagetype_postscript

  implicit none

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Writing PostScript vector plot for time step ',it
    call flush_IMAIN()
  endif

  ! determines postscript output type
  if (P_SV) then
    select case (imagetype_postscript)
    case (1)
      ! displacement
      if (myrank == 0) write(IMAIN,*) 'drawing displacement vector as small arrows...'
      call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic)
    case (2)
      ! velocity
      if (myrank == 0) write(IMAIN,*) 'drawing velocity vector as small arrows...'
      call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic)
    case (3)
      ! acceleration
      if (myrank == 0) write(IMAIN,*) 'drawing acceleration vector as small arrows...'
      call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,accels_poroelastic)
    case default
      call exit_MPI(myrank,'wrong type for PostScript snapshots')
    end select

    ! postscript plotting
    call plot_post()

  else
    call exit_MPI(myrank,'cannot draw a SH scalar field as a vector plot, turn PostScript plots off')
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'PostScript file written'
    call flush_IMAIN()
  endif

  end subroutine write_postscript_snapshot

