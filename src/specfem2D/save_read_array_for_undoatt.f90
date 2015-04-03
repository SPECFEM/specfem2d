
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
!=====================================================================
  subroutine save_forward_arrays_undoatt(iteration_on_subset)

  use specfem_par

  implicit none

  ! local parameters
  integer :: iteration_on_subset,iteration_on_subset_tmp


  ! current subset iteration
  iteration_on_subset_tmp = iteration_on_subset

  ! saves frame of the forward simulation

  write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_save_frame_at',iteration_on_subset_tmp,'.bin'
  open(unit=IOUT_UNDO_ATT  ,file='OUTPUT_FILES/'//outputname, &
       status='unknown',form='unformatted',action='write',iostat=ier)
  if( ier /= 0 ) call exit_MPI('error opening file proc***_save_frame_at** for writing')

  if( any_acoustic ) then
    write( IOUT_UNDO_ATT ) potential_dot_dot_acoustic
    write( IOUT_UNDO_ATT ) potential_dot_acoustic
    write( IOUT_UNDO_ATT ) potential_acoustic
  endif

  if( any_elastic ) then
    write(IOUT_UNDO_ATT ) accel_elastic
    write(IOUT_UNDO_ATT ) veloc_elastic
    write(IOUT_UNDO_ATT ) displ_elastic

    if( ATTENUATION_VISCOELASTIC_SOLID ) then
      write(IOUT_UNDO_ATT ) e1
      write(IOUT_UNDO_ATT ) e11
      write(IOUT_UNDO_ATT ) e13
    endif
  endif

  close(IOUT_UNDO_ATT)

  end subroutine save_forward_arrays_undoatt

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays_undoatt(iteration_on_subset)

! reads in saved wavefields

  use specfem_par

  implicit none

  ! local parameters
  integer :: iteration_on_subset,iteration_on_subset_tmp

  ! current subset iteration
  iteration_on_subset_tmp = NSTEP/NT_DUMP_ATTENUATION - iteration_on_subset + 1

  ! reads in saved wavefield
  write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_save_frame_at',iteration_on_subset_tmp,'.bin'

  ! opens corresponding snapshot file for reading
  open(unit=IIN_UNDO_ATT,file='OUTPUT_FILES/'//outputname, &
       status='old',action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) call exit_MPI('error opening file proc***_save_frame_at** for reading')

  if( any_acoustic ) then
    read( IIN_UNDO_ATT ) b_potential_dot_dot_acoustic
    read( IIN_UNDO_ATT ) b_potential_dot_acoustic
    read( IIN_UNDO_ATT ) b_potential_acoustic
  endif

  if( any_elastic ) then
    read(IIN_UNDO_ATT ) b_accel_elastic
    read(IIN_UNDO_ATT ) b_veloc_elastic
    read(IIN_UNDO_ATT ) b_displ_elastic

    if( ATTENUATION_VISCOELASTIC_SOLID ) then
      read(IIN_UNDO_ATT ) b_e1
      read(IIN_UNDO_ATT ) b_e11
      read(IIN_UNDO_ATT ) b_e13
    endif
  endif

  close(IOUT_UNDO_ATT)

  end subroutine read_forward_arrays_undoatt

