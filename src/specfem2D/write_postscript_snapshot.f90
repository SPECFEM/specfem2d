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
  subroutine write_postscript_snapshot()

  use specfem_par, only: myrank,p_sv,it,imagetype_postscript, &
                         potential_acoustic,potential_gravitoacoustic, &
                         potential_gravito,displ_elastic,displs_poroelastic, &
                         potential_dot_acoustic,potential_dot_gravitoacoustic, &
                         potential_dot_gravito,veloc_elastic,velocs_poroelastic, &
                         potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
                         potential_dot_dot_gravito,accel_elastic,accels_poroelastic

  implicit none
  include "constants.h"

  if( myrank == 0 ) then
    write(IOUT,*)
    write(IOUT,*) 'Writing PostScript vector plot for time step ',it
  endif

  if( imagetype_postscript == 1 .and. p_sv ) then

    if( myrank == 0 ) write(IOUT,*) 'drawing displacement vector as small arrows...'
    call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                 potential_gravito,displ_elastic,displs_poroelastic)

    call plotpost()

  else if( imagetype_postscript == 2 .and. p_sv ) then

    if( myrank == 0 ) write(IOUT,*) 'drawing velocity vector as small arrows...'
    call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
                 potential_dot_gravito,veloc_elastic,velocs_poroelastic)

    call plotpost()

  else if( imagetype_postscript == 3 .and. p_sv ) then

    if( myrank == 0) write(IOUT,*) 'drawing acceleration vector as small arrows...'
    call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
                 potential_dot_dot_gravito,accel_elastic,accels_poroelastic)

    call plotpost()

  else if( .not. p_sv ) then
    call exit_MPI('cannot draw a SH scalar field as a vector plot, turn PostScript plots off')
  else
    call exit_MPI('wrong type for PostScript snapshots')
  endif

  if( myrank == 0 .and. imagetype_postscript /= 4 .and. p_sv ) write(IOUT,*) 'PostScript file written'

  end subroutine write_postscript_snapshot

