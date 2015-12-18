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
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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


! Partitioning using SCOTCH

  subroutine scotch_partitioning()

#ifdef USE_SCOTCH
  use part_unstruct_par,only: nb_edges,part,nelmnts,xadj_g,adjncy_g,adjwgt

  use parameter_file_par,only: nparts => NPROC
#endif

  implicit none

  include "constants.h"

#ifdef USE_SCOTCH
  include "scotchf.h"

  double precision, dimension(SCOTCH_GRAPHDIM)  :: SCOTCHGRAPH
  double precision, dimension(SCOTCH_STRATDIM)  :: SCOTCHSTRAT
  integer :: IERR

  ! we use the default strategy for partitioning
  ! thus no need to define an explicit strategy
  call scotchfstratinit (SCOTCHSTRAT(1), IERR)
   IF (IERR /= 0) THEN
     PRINT *, 'ERROR : MAIN : Cannot initialize strat'
     STOP
  endif

  CALL SCOTCHFGRAPHINIT (SCOTCHGRAPH (1), IERR)
  IF (IERR /= 0) THEN
     PRINT *, 'ERROR : MAIN : Cannot initialize graph'
     STOP
  endif

  ! fills graph structure : see user manual (scotch_user5.1.pdf, page 72/73)
  ! arguments: #(1) graph_structure       #(2) baseval(either 0/1)    #(3) number_of_vertices
  !                    #(4) adjacency_index_array         #(5) adjacency_end_index_array (optional)
  !                    #(6) vertex_load_array (optional) #(7) vertex_label_array
  !                    #(7) number_of_arcs                    #(8) adjacency_array
  !                    #(9) arc_load_array (optional)      #(10) ierror
  CALL SCOTCHFGRAPHBUILD (SCOTCHGRAPH (1), 0, nelmnts, &
                          xadj_g(0), xadj_g(0), &
                          xadj_g(0), xadj_g(0), &
                          nb_edges, &
                          adjncy_g(0), adjwgt (0), IERR)
  IF (IERR /= 0) THEN
     PRINT *, 'ERROR : MAIN : Cannot build graph'
     STOP
  endif

  CALL SCOTCHFGRAPHCHECK (SCOTCHGRAPH (1), IERR)
  IF (IERR /= 0) THEN
     PRINT *, 'ERROR : MAIN : Invalid check'
     STOP
  endif

  call scotchfgraphpart (SCOTCHGRAPH (1), nparts, SCOTCHSTRAT(1), part(0), IERR)
  IF (IERR /= 0) THEN
     PRINT *, 'ERROR : MAIN : Cannot part graph'
     STOP
  endif

  CALL SCOTCHFGRAPHEXIT (SCOTCHGRAPH (1), IERR)
  IF (IERR /= 0) THEN
     PRINT *, 'ERROR : MAIN : Cannot destroy graph'
     STOP
  endif

  call scotchfstratexit (SCOTCHSTRAT(1), IERR)
  IF (IERR /= 0) THEN
     PRINT *, 'ERROR : MAIN : Cannot destroy strat'
     STOP
  endif

#else
  ! safety stop
  print *, 'This version of SPECFEM was not compiled with support of SCOTCH.'
  print *, 'Please recompile with -DUSE_SCOTCH in order to enable use of SCOTCH.'
  stop 'SCOTCH partitioning not compiled'
#endif

  end subroutine scotch_partitioning

