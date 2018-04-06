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

! Partitioning using SCOTCH

  subroutine scotch_partitioning()

#ifdef USE_SCOTCH
  use part_unstruct_par, only: nb_edges,part,nelmnts,xadj_g,adjncy_g
  use compute_elements_load_par, only: elmnts_load,adjwgt

  use shared_parameters, only: nparts => NPROC
#endif

  implicit none

#ifdef USE_SCOTCH
  include "scotchf.h"

  double precision, dimension(SCOTCH_GRAPHDIM)  :: SCOTCHGRAPH
  double precision, dimension(SCOTCH_STRATDIM)  :: SCOTCHSTRAT
  integer :: ier

! workflow preferred by F. Pellegrini (SCOTCH):
!!
!!This comes from the fact that, in version 5.1.8, the name
!!for the "recursive bisection" method has changed from "b"
!!("bipartitioning") to "r" ("recursive").
!!
!!As a general rule, do not try to set up strategies by
!!yourself. The default strategy in Scotch will most probably
!!provide better results. To use it, just call:
!!
!!SCOTCHFstratInit (),
!!
!!and use this "empty" strategy in the mapping routine
!!(consequently, no call to SCOTCHFstratGraphMap () is
!!required).
!!
!!This will make you independent from further changes
!!(improvements) in the strategy syntax.
!!And you should see an improvement in performance, too,
!!as your hand-made strategy did not make use of the
!!multi-level framework.

  ! we use the default strategy for partitioning
  ! thus no need to define an explicit strategy
  call scotchfstratinit (SCOTCHSTRAT(1), ier)
   if (ier /= 0) then
     print *, 'ERROR : MAIN : Cannot initialize strategy'
     call stop_the_code('Error scotch init')
  endif

  ! resets SCOTCH random number generator to produce deterministic partitions
  call scotchfrandomReset()

  ! initializes graph
  call scotchfgraphinit (SCOTCHGRAPH (1), ier)
  if (ier /= 0) then
     print *, 'ERROR : MAIN : Cannot initialize graph'
     call stop_the_code('Error scotch graph')
  endif

  ! fills graph structure : see user manual (scotch_user5.1.pdf, page 72/73)
  ! arguments: #(1) graph_structure       #(2) baseval(either 0/1)    #(3) number_of_vertices
  !                    #(4) adjacency_index_array         #(5) adjacency_end_index_array (optional)
  !                    #(6) vertex_load_array (optional) #(7) vertex_label_array
  !                    #(7) number_of_arcs                    #(8) adjacency_array
  !                    #(9) arc_load_array (optional)      #(10) ierror
  call scotchfgraphbuild (SCOTCHGRAPH (1), 0, nelmnts, &
                          xadj_g(0), xadj_g(0), &
                          elmnts_load(0), xadj_g(0), &
                          nb_edges, &
                          adjncy_g(0), adjwgt (0), ier)
  if (ier /= 0) then
     print *, 'ERROR : MAIN : Cannot build graph'
     call stop_the_code('Error scotch graphbuild')
  endif

  call scotchfgraphcheck (SCOTCHGRAPH (1), ier)
  if (ier /= 0) then
     print *, 'ERROR : MAIN : Invalid check'
     call stop_the_code('Error scotch graphcheck')
  endif

  call scotchfgraphpart (SCOTCHGRAPH (1), nparts, SCOTCHSTRAT(1), part(0), ier)
  if (ier /= 0) then
     print *, 'ERROR : MAIN : Cannot part graph'
     call stop_the_code('Error scotch graphpart')
  endif

  call SCOTCHFGRAPHEXIT (SCOTCHGRAPH (1), ier)
  if (ier /= 0) then
     print *, 'ERROR : MAIN : Cannot destroy graph'
     call stop_the_code('Error scotch graphexit')
  endif

  call scotchfstratexit (SCOTCHSTRAT(1), ier)
  if (ier /= 0) then
     print *, 'ERROR : MAIN : Cannot destroy strat'
     call stop_the_code('Error scotch exit')
  endif

#else
  ! safety stop
  print *, 'This version of SPECFEM was not compiled with support of SCOTCH.'
  print *, 'Please recompile with -DUSE_SCOTCH in order to enable use of SCOTCH.'
  call stop_the_code('Error SCOTCH partitioning not compiled')
#endif

  end subroutine scotch_partitioning

