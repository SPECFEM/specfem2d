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

  use constants, only: IOUT_UNDO_ATT,MAX_STRING_LEN,OUTPUT_FILES,APPROXIMATE_HESS_KL,NDIM

  use specfem_par, only: myrank,it,NSTEP, &
    any_acoustic,any_elastic,potential_acoustic,displ_elastic,accel_elastic,GPU_MODE, &
    no_backward_acoustic_buffer,no_backward_displ_buffer,no_backward_accel_buffer, &
    no_backward_iframe,no_backward_Nframes,nglob

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! EB EB June 2018 : in this routine, in order to overlap both GPU = => RAM and RAM = => disk transfers, we transfer
  ! a wavefield in two iterations.
  ! At the first iteration, we transfer the wavefield from the GPU to the disk.
  ! At the second iteration, we write this wavefield on the disk from the RAM.
  ! In the text above, an iteration means NSTEP_BETWEEN_COMPUTE_KERNELS iterations of the timeloop.
  ! The buffer no_backward_acoustic_buffer is declared in only one dimension in
  ! order to allow the CUDA API to set it in pinned memory (HostRegister).
  ! To perform the async I/O, stream accesses are used for files, numerical
  ! experiences showed that it is the fastest way.

  ! local parameters
  integer :: ier,buffer_num_GPU_transfer,buffer_num_async_IO
  integer(KIND=8) :: offset
  character(len=MAX_STRING_LEN) :: outputname

  ! safety check
  if (GPU_MODE .and. any_elastic) call stop_the_code('No backward simulation is not available for elastic on GPU')

  ! increments counter of wavefield frames transfered
  no_backward_iframe = no_backward_iframe + 1

  ! opens file to save at the first use of the routine
  if (no_backward_iframe == 1) then
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_No_backward_reconstruction_database.bin'
    open(unit=IOUT_UNDO_ATT,file=trim(OUTPUT_FILES)//outputname, &
         status='unknown',asynchronous='yes',form='unformatted',action='write',access='stream',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_No_backward_reconstruction_database.bin for writing')
  else if (no_backward_iframe > 2) then
    ! waits for previous asynchronous I/O
    wait(IOUT_UNDO_ATT)
  endif

  buffer_num_GPU_transfer = mod(no_backward_iframe+2,3)
  buffer_num_async_IO = mod(no_backward_iframe,3)

  ! for the two first times, we only launch GPU = => RAM transfers
  if (no_backward_iframe < 3) then

    if (GPU_MODE) then
      call transfer_async_pot_ac_from_device(no_backward_acoustic_buffer(nglob*buffer_num_GPU_transfer+1),Mesh_pointer)
    else
      no_backward_acoustic_buffer(nglob*buffer_num_GPU_transfer+1:nglob*(buffer_num_GPU_transfer+1)) = potential_acoustic
    endif

  else if (no_backward_iframe <= no_backward_Nframes) then

    if (any_acoustic) then

      ! the offset is calculated in two steps in order to avoid integer overflow
      offset = no_backward_iframe - 3
      offset = offset * nglob * 4 + 1
      write(IOUT_UNDO_ATT,asynchronous='yes',pos=offset) &
                          no_backward_acoustic_buffer(nglob*buffer_num_async_IO+1:nglob*(buffer_num_async_IO+1))

      if (GPU_MODE) then
        call transfer_async_pot_ac_from_device(no_backward_acoustic_buffer(nglob*buffer_num_GPU_transfer+1),Mesh_pointer)
      else
        no_backward_acoustic_buffer(nglob*mod(no_backward_iframe+1,3)+1:nglob*(mod(no_backward_iframe+1,3)+1)) = potential_acoustic
      endif

      ! for the last transfer, we need to add a statement to wait for the last frame
      if (no_backward_iframe == no_backward_Nframes) then
        ! call to finalize disk writing
        wait(IOUT_UNDO_ATT)

        ! call to finalize GPU transfer, which also initiate a (dummy) transfer
        if (GPU_MODE) &
          call transfer_async_pot_ac_from_device(no_backward_acoustic_buffer(nglob*buffer_num_async_IO+1),Mesh_pointer)

        write(IOUT_UNDO_ATT,pos=offset+4*nglob) &
                            no_backward_acoustic_buffer(nglob*buffer_num_GPU_transfer+1:nglob*(buffer_num_GPU_transfer+1))
        write(IOUT_UNDO_ATT,asynchronous='yes',pos=offset+8*nglob) &
                            no_backward_acoustic_buffer(nglob*mod(no_backward_iframe+1,3)+1:nglob*(mod(no_backward_iframe+1,3)+1))

      endif
    endif

    if (any_elastic) then

      if (APPROXIMATE_HESS_KL) then
        offset = 4*2*(NDIM*nglob)*(no_backward_Nframes - no_backward_iframe) + 1
      else
        offset = 4*(NDIM*nglob)*(no_backward_Nframes - no_backward_iframe) + 1
      endif

      no_backward_displ_buffer(:,:) = displ_elastic(:,:)
      write(IOUT_UNDO_ATT,asynchronous='yes',pos=offset) no_backward_displ_buffer
      if (APPROXIMATE_HESS_KL) then
        no_backward_accel_buffer(:,:) = accel_elastic(:,:)
        write(IOUT_UNDO_ATT,asynchronous='yes',pos=offset) no_backward_accel_buffer
      endif

    endif

  endif

  ! this operation will automatically synchronize the remaining I/O to do
  if (it == NSTEP) close(IOUT_UNDO_ATT)

  end subroutine save_forward_arrays_no_backward

