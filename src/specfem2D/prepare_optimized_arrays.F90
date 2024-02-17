!========================================================================
!
!                            S P E C F E M 2 D
!                            -----------------
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


! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJ, DO_LOOP_IJ, ENDDO_LOOP_IJ defined in config.fh
#include "config.fh"


  subroutine prepare_optimized_arrays()

! optimizes array memory layout to increase computational efficiency

  use constants, only: myrank,IMAIN

  implicit none

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing optimized arrays"
    call flush_IMAIN()
  endif

#ifdef USE_OPENMP
  ! prepares arrays for OpenMP
  call prepare_timerun_OpenMP()
#endif

  ! prepare fused array for computational kernel
  call prepare_fused_array()

  end subroutine prepare_optimized_arrays


!
!-------------------------------------------------------------------------------------------------
!

! OpenMP version uses "special" compute_forces_viscoelastic routine
! we need to set num_elem_colors_elastic arrays

#ifdef USE_OPENMP
  subroutine prepare_timerun_OpenMP()

  use constants, only: IMAIN
  use specfem_par

  implicit none

  ! local parameters
  integer :: ier
  integer :: max_threads
  integer,external :: OMP_GET_MAX_THREADS

  ! gets number of openMP threads
  ! will be determined by environment setting OMP_NUM_THREADS
  !
  ! for example, run executable with:
  ! OMP_NUM_THREADS=4 mpirun -np 2 ./bin/xspecfem3D
  !
  max_threads = OMP_GET_MAX_THREADS()

  ! output info
  if (myrank == 0) then
    write(IMAIN,*) '  OpenMP:'
    write(IMAIN,*) '    using',max_threads,' OpenMP threads'
    call flush_IMAIN()
  endif

  ! synchonizes
  call synchronize_all()

  end subroutine prepare_timerun_OpenMP
#endif


!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_fused_array()

! prepare fused array for computational kernel
!
! note: fusing separate arrays for xi/eta/gamma into a single one increases efficiency for hardware pre-fetching

  use constants, only: IMAIN
  use specfem_par

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLSQUARE
#endif

  implicit none

  ! local parameters
  integer :: ispec,ier

#ifdef FORCE_VECTORIZATION
  integer :: ij
#else
  integer :: i,j
#endif
  real(kind=CUSTOM_REAL) :: rhol

  ! allocates fused array
  allocate(deriv_mapping(6,NGLLX,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) call exit_mpi(myrank,'Error allocating array deriv_mapping')

  ! fused array of mapping matrix
  ! (mapping from reference element back to physical element)
  ! d(xi)/d(x)
  do ispec = 1,NSPEC

    DO_LOOP_IJ
      deriv_mapping(1,INDEX_IJ,ispec) = xix(INDEX_IJ,ispec)
      deriv_mapping(2,INDEX_IJ,ispec) = xiz(INDEX_IJ,ispec)

      deriv_mapping(3,INDEX_IJ,ispec) = gammax(INDEX_IJ,ispec)
      deriv_mapping(4,INDEX_IJ,ispec) = gammaz(INDEX_IJ,ispec)

      deriv_mapping(5,INDEX_IJ,ispec) = jacobian(INDEX_IJ,ispec)

      ! for acoustic elements
      ! density model
      rhol = rhostore(INDEX_IJ,ispec)
      deriv_mapping(6,INDEX_IJ,ispec) = jacobian(INDEX_IJ,ispec) / rhol
    ENDDO_LOOP_IJ

  enddo

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)"  fused array done"
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_fused_array
