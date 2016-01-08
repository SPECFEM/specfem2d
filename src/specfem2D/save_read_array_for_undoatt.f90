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
!=====================================================================

  subroutine save_forward_arrays_undoatt()

  use constants,only: IOUT_UNDO_ATT,MAX_STRING_LEN

  use specfem_par,only: myrank,iteration_on_subset, &
    any_acoustic,any_elastic,ATTENUATION_VISCOELASTIC_SOLID, &
    potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
    displ_elastic,veloc_elastic,accel_elastic, &
    e1,e11,e13


  implicit none

  ! local parameters
  integer :: iteration_on_subset_tmp
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname

  ! current subset iteration
  iteration_on_subset_tmp = iteration_on_subset

  ! saves frame of the forward simulation

  write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_save_frame_at',iteration_on_subset_tmp,'.bin'
  open(unit=IOUT_UNDO_ATT  ,file='OUTPUT_FILES/'//outputname, &
       status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_save_frame_at** for writing')

  if (any_acoustic) then
    write(IOUT_UNDO_ATT) potential_dot_dot_acoustic
    write(IOUT_UNDO_ATT) potential_dot_acoustic
    write(IOUT_UNDO_ATT) potential_acoustic
  endif

  if (any_elastic) then
    write(IOUT_UNDO_ATT) accel_elastic
    write(IOUT_UNDO_ATT) veloc_elastic
    write(IOUT_UNDO_ATT) displ_elastic

    if (ATTENUATION_VISCOELASTIC_SOLID) then
      write(IOUT_UNDO_ATT) e1
      write(IOUT_UNDO_ATT) e11
      write(IOUT_UNDO_ATT) e13
    endif
  endif

  close(IOUT_UNDO_ATT)

  end subroutine save_forward_arrays_undoatt

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays_undoatt()

! reads in saved wavefields

  use constants,only: IIN_UNDO_ATT,MAX_STRING_LEN

  use specfem_par,only: myrank, &
    iteration_on_subset,NSUBSET_ITERATIONS, &
    any_acoustic,any_elastic,ATTENUATION_VISCOELASTIC_SOLID, &
    b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
    b_displ_elastic,b_veloc_elastic,b_accel_elastic, &
    b_e1,b_e11,b_e13

  implicit none

  ! local parameters
  integer :: iteration_on_subset_tmp
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname

  ! current subset iteration
  iteration_on_subset_tmp = NSUBSET_ITERATIONS - iteration_on_subset + 1

  ! reads in saved wavefield
  write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_save_frame_at',iteration_on_subset_tmp,'.bin'

  ! opens corresponding snapshot file for reading
  open(unit=IIN_UNDO_ATT,file='OUTPUT_FILES/'//outputname, &
       status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_save_frame_at** for reading')

  if (any_acoustic) then
    read(IIN_UNDO_ATT) b_potential_dot_dot_acoustic
    read(IIN_UNDO_ATT) b_potential_dot_acoustic
    read(IIN_UNDO_ATT) b_potential_acoustic
  endif

  if (any_elastic) then
    read(IIN_UNDO_ATT) b_accel_elastic
    read(IIN_UNDO_ATT) b_veloc_elastic
    read(IIN_UNDO_ATT) b_displ_elastic

    if (ATTENUATION_VISCOELASTIC_SOLID) then
      read(IIN_UNDO_ATT) b_e1
      read(IIN_UNDO_ATT) b_e11
      read(IIN_UNDO_ATT) b_e13
    endif
  endif

  close(IIN_UNDO_ATT)

  end subroutine read_forward_arrays_undoatt

