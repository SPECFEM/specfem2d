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

  subroutine write_movie_output(get_b_wavefield)

! outputs snapshots for movies

  use constants, only: MAX_STRING_LEN,CUSTOM_REAL,NDIM,NOISE_MOVIE_OUTPUT,OUTPUT_FILES

  use specfem_par, only: myrank,it,NSTEP,nspec,nglob,ibool, &
    potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
    b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
    displ_elastic,veloc_elastic,accel_elastic, &
    b_displ_elastic,rho_k,rho_kl, &
    any_acoustic,any_elastic,GPU_MODE,P_SV,UNDO_ATTENUATION_AND_OR_PML,SIMULATION_TYPE,NO_BACKWARD_RECONSTRUCTION

  use specfem_par_gpu, only: Mesh_pointer,tmp_displ_2D,tmp_veloc_2D,tmp_accel_2D,NGLOB_AB

  use specfem_par_noise, only: NOISE_TOMOGRAPHY,mask_noise, &
                              surface_movie_y_or_z_noise,noise_output_rhokl,noise_output_array,noise_output_ncol

  use specfem_par_movie, only: output_postscript_snapshot,output_color_image,output_wavefield_dumps, &
    NSTEP_BETWEEN_OUTPUT_IMAGES

  implicit none

  ! parameter useful for UNDO_ATTENUATION mode
  logical :: get_b_wavefield

  ! local parameters
  integer :: ier
  logical :: ex, is_opened,plot_b_wavefield_only
  character(len=MAX_STRING_LEN) :: noise_output_file

  ! checks if anything to do
  if (.not. (mod(it,NSTEP_BETWEEN_OUTPUT_IMAGES) == 0 .or. it == 5 .or. it == NSTEP)) return
  if ( (.not. output_postscript_snapshot) .and. (.not. NOISE_MOVIE_OUTPUT) .and. (.not. output_color_image) &
       .and. (.not. output_wavefield_dumps) ) return

  if (UNDO_ATTENUATION_AND_OR_PML .and. get_b_wavefield .and. .not. NO_BACKWARD_RECONSTRUCTION) then
    plot_b_wavefield_only = .true.
  else
    plot_b_wavefield_only = .false.
  endif

  ! transfers arrays from GPU to CPU
  if (GPU_MODE) then
    ! Fields transfer for imaging
    ! acoustic domains
    if (any_acoustic) then
      if (.not. plot_b_wavefield_only) &
        call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic,potential_dot_acoustic, &
                                            potential_dot_dot_acoustic,Mesh_pointer)
      if ((SIMULATION_TYPE == 3 .and. (.not. NO_BACKWARD_RECONSTRUCTION) .and. &
          (.not. UNDO_ATTENUATION_AND_OR_PML)) .or. plot_b_wavefield_only ) &
          call transfer_b_fields_ac_from_device(NGLOB_AB,b_potential_acoustic,b_potential_dot_acoustic, &
                                                b_potential_dot_dot_acoustic,Mesh_pointer)
    endif
    ! elastic domains
    if (any_elastic) then
      call transfer_fields_el_from_device(NDIM*NGLOB_AB,tmp_displ_2D,tmp_veloc_2D,tmp_accel_2D,Mesh_pointer)
      if (P_SV) then
        ! P-SV waves
        displ_elastic(1,:) = tmp_displ_2D(1,:)
        displ_elastic(2,:) = tmp_displ_2D(2,:)
        veloc_elastic(1,:) = tmp_veloc_2D(1,:)
        veloc_elastic(2,:) = tmp_veloc_2D(2,:)
        accel_elastic(1,:) = tmp_accel_2D(1,:)
        accel_elastic(2,:) = tmp_accel_2D(2,:)
      else
        ! SH waves
        displ_elastic(1,:) = tmp_displ_2D(1,:)
        veloc_elastic(1,:) = tmp_veloc_2D(1,:)
        accel_elastic(1,:) = tmp_accel_2D(1,:)
      endif
    endif
  endif

  ! noise simulations
  if (.not. GPU_MODE) then
    if (NOISE_TOMOGRAPHY == 3 .and. NOISE_MOVIE_OUTPUT) then
      ! load ensemble forward source
      inquire(unit=501,exist=ex,opened=is_opened)
      if (.not. is_opened) then
        open(unit=501,file=trim(OUTPUT_FILES)//'noise_eta.bin',access='direct', &
             recl=nglob*CUSTOM_REAL,action='read',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening noise eta file')
      endif
      ! for both P_SV/SH cases
      read(unit=501,rec=it) surface_movie_y_or_z_noise ! either y (SH) or z (P_SV) component

      ! load product of fwd, adj wavefields
      call spec2glob(nspec,nglob,ibool,rho_kl,noise_output_rhokl)

      ! prepares array
      ! noise distribution
      noise_output_array(1,:) = surface_movie_y_or_z_noise(:) * mask_noise(:)

      ! P_SV/SH-case
      noise_output_array(2,:) = b_displ_elastic(1,:)
      noise_output_array(3,:) = accel_elastic(1,:)

      ! rho kernel on global nodes
      noise_output_array(4,:) = rho_k(:)

      ! rho kernel on global nodes from local kernel (for comparison)
      noise_output_array(5,:) = noise_output_rhokl(:)

      ! writes out to text file
      write(noise_output_file,"(a,i6.6,a)") trim(OUTPUT_FILES)//'noise_snapshot_all_',it,'.txt'
      call snapshots_noise(noise_output_ncol,nglob,noise_output_file,noise_output_array)

      ! re-initializes noise array
      surface_movie_y_or_z_noise(:) = 0._CUSTOM_REAL
    endif
  endif

  ! output_postscript_snapshot
  if (output_postscript_snapshot) then
    call write_postscript_snapshot()
  endif

  ! display color image
  if (output_color_image) then
    call write_color_image_snaphot(plot_b_wavefield_only)
  endif

  ! dump the full (local) wavefield to a file
  ! note: in the case of MPI, in the future it would be more convenient to output a single file
  !       rather than one for each myrank
  if (output_wavefield_dumps) then
    call write_wavefield_dumps()
  endif

  end subroutine write_movie_output


