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
!=====================================================================

  subroutine save_forward_arrays_undoatt()

  use constants, only: IOUT_UNDO_ATT,MAX_STRING_LEN,OUTPUT_FILES,NGLLX,NGLLZ

  use specfem_par, only: myrank,iteration_on_subset, &
    any_acoustic,any_elastic,ATTENUATION_VISCOACOUSTIC,ATTENUATION_VISCOELASTIC, &
    potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
    displ_elastic,veloc_elastic,accel_elastic, &
    e1,e11,e13,dux_dxl_old,duz_dzl_old,dux_dzl_plus_duz_dxl_old, &
    e1_acous_sf,sum_forces_old,GPU_MODE,nspec_ATT_ac,nglob

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! local parameters
  integer :: iteration_on_subset_tmp
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname

  ! current subset iteration
  iteration_on_subset_tmp = iteration_on_subset

  ! saves frame of the forward simulation

  write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_save_frame_at',iteration_on_subset_tmp,'.bin'
  open(unit=IOUT_UNDO_ATT  ,file=trim(OUTPUT_FILES)//outputname, &
       status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_save_frame_at** for writing')

  if (any_acoustic) then
    if (GPU_MODE) call transfer_fields_ac_from_device(nglob,potential_acoustic,potential_dot_acoustic, &
                                                      potential_dot_dot_acoustic,Mesh_pointer)
    write(IOUT_UNDO_ATT) potential_dot_dot_acoustic
    write(IOUT_UNDO_ATT) potential_dot_acoustic
    write(IOUT_UNDO_ATT) potential_acoustic

    if (ATTENUATION_VISCOACOUSTIC) then
      if (GPU_MODE) call transfer_viscoacoustic_var_from_device(NGLLX*NGLLZ*nspec_ATT_ac, &
                                                                e1_acous_sf,sum_forces_old,Mesh_pointer)
      write(IOUT_UNDO_ATT) e1_acous_sf
      write(IOUT_UNDO_ATT) sum_forces_old

    endif

  endif

  if (any_elastic) then
    write(IOUT_UNDO_ATT) accel_elastic
    write(IOUT_UNDO_ATT) veloc_elastic
    write(IOUT_UNDO_ATT) displ_elastic

    if (ATTENUATION_VISCOELASTIC) then
      write(IOUT_UNDO_ATT) e1
      write(IOUT_UNDO_ATT) e11
      write(IOUT_UNDO_ATT) e13
      write(IOUT_UNDO_ATT) dux_dxl_old
      write(IOUT_UNDO_ATT) duz_dzl_old
      write(IOUT_UNDO_ATT) dux_dzl_plus_duz_dxl_old
    endif
  endif

  close(IOUT_UNDO_ATT)

  end subroutine save_forward_arrays_undoatt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine save_forward_arrays_no_backward()

  use constants, only: IOUT_UNDO_ATT,MAX_STRING_LEN,OUTPUT_FILES,APPROXIMATE_HESS_KL

  use specfem_par, only: myrank,it,NSTEP,NSTEP_BETWEEN_COMPUTE_KERNELS, &
    any_acoustic,any_elastic,potential_acoustic, &
    displ_elastic,accel_elastic,nglob,no_backward_nframes, &!GPU_MODE
    no_backward_acoustic_buffer,no_backward_displ_buffer,no_backward_accel_buffer,no_backward_iframe

!  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname


  ! opens file to save at the first use of the routine
  if (it == 1) then
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_No_backward_reconstruction_database.bin'
    open(unit=IOUT_UNDO_ATT,file=trim(OUTPUT_FILES)//outputname, &
         status='unknown',asynchronous='yes',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_No_backward_reconstruction_database.bin for writing')
    no_backward_iframe = 0
  ! waits for previous asynchronous I/O
  else
    wait(IOUT_UNDO_ATT)
  endif

  if (NSTEP_BETWEEN_COMPUTE_KERNELS == 1 .or. it /= NSTEP) then
    if (any_acoustic) then
      no_backward_iframe = no_backward_iframe + 1
      !if (GPU_MODE) call transfer_fields_ac_from_device(nglob,potential_acoustic,potential_dot_acoustic,&
       !                                                 potential_dot_dot_acoustic,Mesh_pointer)
      no_backward_acoustic_buffer(:,2) = potential_acoustic
      write(IOUT_UNDO_ATT,asynchronous='yes') no_backward_acoustic_buffer(:,2)
    endif

    if (any_elastic) then
      write(IOUT_UNDO_ATT) accel_elastic
      write(IOUT_UNDO_ATT) displ_elastic
    endif
  endif

  ! this operation will automatically synchronize the remaining I/O to do
  if (it == NSTEP) then
   close(IOUT_UNDO_ATT)
  endif
  end subroutine save_forward_arrays_no_backward

