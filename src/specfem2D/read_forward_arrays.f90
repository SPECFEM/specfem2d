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


  subroutine read_forward_arrays()

! restores last time snapshot saved for backward/reconstruction of wavefields
! note: this is done here after the Newmark time scheme, otherwise the indexing for sources
!          and adjoint sources will become more complicated
!          that is, index it for adjoint sources will match index NSTEP - 1 for backward/reconstructed wavefields

  use constants, only: OUTPUT_FILES
  use specfem_par
  use specfem_par_gpu

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname

  ! acoustic medium
  if (any_acoustic) then
    write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
    open(unit=55,file=trim(OUTPUT_FILES)//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'lastframe_acoustic**.bin')

    read(55) b_potential_acoustic
    read(55) b_potential_dot_acoustic
    read(55) b_potential_dot_dot_acoustic

    close(55)

    if (GPU_MODE) then
      ! transfers fields onto GPU
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic, &
                                          b_potential_dot_acoustic, &
                                          b_potential_dot_dot_acoustic, &
                                          Mesh_pointer)
    else
      ! free surface for an acoustic medium
      call enforce_acoustic_free_surface(b_potential_dot_dot_acoustic,b_potential_dot_acoustic,b_potential_acoustic)
    endif
  endif

  ! elastic medium
  if (any_elastic) then
    write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
    open(unit=55,file=trim(OUTPUT_FILES)//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'lastframe_elastic**.bin')

    read(55) b_displ_elastic
    read(55) b_veloc_elastic
    read(55) b_accel_elastic
    close(55)

    !SH (membrane) waves
    if (.not. P_SV) then
      ! only index array(1,:) contains SH wavefield
      b_displ_elastic(2,:) = 0._CUSTOM_REAL
      b_veloc_elastic(2,:) = 0._CUSTOM_REAL
      b_accel_elastic(2,:) = 0._CUSTOM_REAL
    endif

    if (GPU_MODE) then
      ! prepares wavefields for transfering
      if (P_SV) then
        tmp_displ_2D(1,:) = b_displ_elastic(1,:)
        tmp_displ_2D(2,:) = b_displ_elastic(2,:)
        tmp_veloc_2D(1,:) = b_veloc_elastic(1,:)
        tmp_veloc_2D(2,:) = b_veloc_elastic(2,:)
        tmp_accel_2D(1,:) = b_accel_elastic(1,:)
        tmp_accel_2D(2,:) = b_accel_elastic(2,:)
      else
        ! SH waves
        tmp_displ_2D(1,:) = b_displ_elastic(1,:)
        tmp_displ_2D(2,:) = 0._CUSTOM_REAL
        tmp_veloc_2D(1,:) = b_veloc_elastic(1,:)
        tmp_veloc_2D(2,:) = 0._CUSTOM_REAL
        tmp_accel_2D(1,:) = b_accel_elastic(1,:)
        tmp_accel_2D(2,:) = 0._CUSTOM_REAL
      endif
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,tmp_displ_2D,tmp_veloc_2D,tmp_accel_2D,Mesh_pointer)
    endif
  endif

  ! poroelastic medium
  if (any_poroelastic) then
    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_s',myrank,'.bin'
    open(unit=55,file=trim(OUTPUT_FILES)//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'lastframe_poroelastic_s**.bin')

    read(55) b_displs_poroelastic
    read(55) b_velocs_poroelastic
    read(55) b_accels_poroelastic
    close(55)

    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_w',myrank,'.bin'
    open(unit=56,file=trim(OUTPUT_FILES)//trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'lastframe_poroelastic_w**.bin')

    read(56) b_displw_poroelastic
    read(56) b_velocw_poroelastic
    read(56) b_accelw_poroelastic
    close(56)

    ! safety check
    if (GPU_MODE) then
      call stop_the_code('GPU_MODE error: sorry, reading lastframe from poroelastic simulation not implemented yet')
    endif
  endif

  end subroutine read_forward_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays_undoatt()

! reads in saved wavefields

  use constants, only: IIN_UNDO_ATT,MAX_STRING_LEN,OUTPUT_FILES,NGLLX,NGLLZ

  use specfem_par, only: myrank,iteration_on_subset,NSUBSET_ITERATIONS, &
    any_acoustic,any_elastic,ATTENUATION_VISCOACOUSTIC,ATTENUATION_VISCOELASTIC, &
    b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
    b_displ_elastic,b_veloc_elastic,b_accel_elastic,b_e1,b_e11,b_e13, &
    b_dux_dxl_old,b_duz_dzl_old,b_dux_dzl_plus_duz_dxl_old,b_e1_acous_sf,b_sum_forces_old, &
    GPU_MODE,nspec_ATT_ac,nglob

  use specfem_par_gpu, only: Mesh_pointer

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
  open(unit=IIN_UNDO_ATT,file=trim(OUTPUT_FILES)//outputname, &
       status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_save_frame_at** for reading')

  if (any_acoustic) then
    read(IIN_UNDO_ATT) b_potential_dot_dot_acoustic
    read(IIN_UNDO_ATT) b_potential_dot_acoustic
    read(IIN_UNDO_ATT) b_potential_acoustic
    if (GPU_MODE) call transfer_b_fields_ac_to_device(nglob,b_potential_acoustic,b_potential_dot_acoustic, &
                                                        b_potential_dot_dot_acoustic,Mesh_pointer)
    if (ATTENUATION_VISCOACOUSTIC) then
      read(IIN_UNDO_ATT) b_e1_acous_sf
      read(IIN_UNDO_ATT) b_sum_forces_old
      if (GPU_MODE) call transfer_viscoacoustic_b_var_to_device(NGLLX*NGLLZ*nspec_ATT_ac,b_e1_acous_sf, &
                                                                b_sum_forces_old,Mesh_pointer)
    endif
  endif

  if (any_elastic) then
    read(IIN_UNDO_ATT) b_accel_elastic
    read(IIN_UNDO_ATT) b_veloc_elastic
    read(IIN_UNDO_ATT) b_displ_elastic

    if (ATTENUATION_VISCOELASTIC) then
      read(IIN_UNDO_ATT) b_e1
      read(IIN_UNDO_ATT) b_e11
      read(IIN_UNDO_ATT) b_e13
      read(IIN_UNDO_ATT) b_dux_dxl_old
      read(IIN_UNDO_ATT) b_duz_dzl_old
      read(IIN_UNDO_ATT) b_dux_dzl_plus_duz_dxl_old
    endif
  endif

  close(IIN_UNDO_ATT)

  end subroutine read_forward_arrays_undoatt

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays_no_backward()

  use constants, only: IIN_UNDO_ATT,MAX_STRING_LEN,OUTPUT_FILES,APPROXIMATE_HESS_KL,NDIM

  use specfem_par, only: myrank,it,any_acoustic,any_elastic, &
    b_potential_acoustic,b_displ_elastic,b_accel_elastic, &
    nglob,no_backward_acoustic_buffer,no_backward_displ_buffer,no_backward_accel_buffer, &
    no_backward_iframe,no_backward_Nframes,GPU_MODE,NSTEP

  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! local parameters
  integer :: ier,buffer_num_async_IO,buffer_num_GPU_transfer
  integer(KIND=8) :: offset
  character(len=MAX_STRING_LEN) :: outputname

  ! EB EB June 2018 : in this routine, in order to overlap both disk = => RAM and RAM = => GPU transfers, we initiate the
  ! transfer of a wavefield two iterations before this wavefield is actually needed by the kernel computation.
  ! Two iterations before, the wavefield is read from the disk.
  ! One iteration before, this wavefield is transfered from the RAM to the GPU.
  ! In the text above, an iteration means NSTEP_BETWEEN_COMPUTE_KERNELS iterations of the timeloop.

  no_backward_iframe = no_backward_iframe + 1

  if (it == 1) then
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_No_backward_reconstruction_database.bin'
    ! opens corresponding file for reading
    open(unit=IIN_UNDO_ATT,asynchronous='yes',file=trim(OUTPUT_FILES)//outputname, &
       status='old',action='read',form='unformatted',access='stream',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_No_backward_reconstruction_database.bin for reading')
  else
    wait(IIN_UNDO_ATT)
  endif

  buffer_num_async_IO = mod(no_backward_iframe+2,3)
  buffer_num_GPU_transfer = mod(no_backward_iframe+1,3)

  if (any_acoustic) then

    ! offset is computed in two times to avoid integer overflow
    offset = 4*nglob
    offset = offset*(no_backward_Nframes - no_backward_iframe ) + 1
    if (no_backward_iframe <= no_backward_Nframes) &
      read(IIN_UNDO_ATT,asynchronous='yes',pos=offset) &
           no_backward_acoustic_buffer(nglob*buffer_num_async_IO+1:nglob*(buffer_num_async_IO+1))

    if (no_backward_iframe /= 1) then
      if (GPU_MODE) then
        ! this function ensures the previous async transfer is finished, and
        ! launches the transfer of the next wavefield
        call transfer_async_pot_ac_to_device(no_backward_acoustic_buffer(nglob*buffer_num_GPU_transfer+1),Mesh_pointer)
      else
        ! we get the wavefield from the previous iteration, because this RAM = =>
        ! RAM copy is blocking
        b_potential_acoustic(:) = no_backward_acoustic_buffer(nglob*mod(no_backward_iframe,3)+1: &
                                                                nglob*(mod(no_backward_iframe,3)+1))
      endif
    endif

    endif ! any_acoustic

    if (any_elastic) then

      b_displ_elastic(:,:) = no_backward_displ_buffer(:,:)
      if (APPROXIMATE_HESS_KL) b_accel_elastic(:,:) = no_backward_accel_buffer(:,:)

      if (APPROXIMATE_HESS_KL) then
        offset = 4*2*(NDIM*nglob)*(no_backward_Nframes - no_backward_iframe) + 1
      else
        offset = 4*(NDIM*nglob)*(no_backward_Nframes - no_backward_iframe) + 1
      endif

      read(IIN_UNDO_ATT,asynchronous='yes',pos=offset) no_backward_displ_buffer(:,:)
      if (APPROXIMATE_HESS_KL) read(IIN_UNDO_ATT,asynchronous='yes',pos=offset+8+4*nglob) no_backward_accel_buffer(:,:)
    endif

  if (it == NSTEP) close(IIN_UNDO_ATT)

  end subroutine read_forward_arrays_no_backward

