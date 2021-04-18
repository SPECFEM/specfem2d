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

  use constants, only: MAX_STRING_LEN,CUSTOM_REAL,NDIM,NOISE_MOVIE_OUTPUT

  use specfem_par, only: it,NSTEP, &
    potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
    b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
    displ_elastic,veloc_elastic,accel_elastic, &
    b_displ_elastic,b_veloc_elastic,b_accel_elastic, &
    any_acoustic,any_elastic, &
    GPU_MODE,UNDO_ATTENUATION_AND_OR_PML,SIMULATION_TYPE,NO_BACKWARD_RECONSTRUCTION

  use shared_parameters, only: output_postscript_snapshot,output_color_image,output_wavefield_dumps, &
    NSTEP_BETWEEN_OUTPUT_IMAGES

  use specfem_par_gpu, only: Mesh_pointer,NGLOB_AB

  implicit none

  ! parameter useful for UNDO_ATTENUATION mode
  logical :: get_b_wavefield

  ! local parameters
  logical :: plot_b_wavefield_only

  ! checks if anything to do
  if (.not. (mod(it,NSTEP_BETWEEN_OUTPUT_IMAGES) == 0 .or. it == 5 .or. it == NSTEP)) return

  if ( (.not. output_postscript_snapshot) .and. (.not. NOISE_MOVIE_OUTPUT) .and. (.not. output_color_image) &
       .and. (.not. output_wavefield_dumps) ) return

  ! checks plotting of backward wavefield
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
      call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ_elastic,veloc_elastic,accel_elastic,Mesh_pointer)
      ! backward/reconstructed wavefield
      if (SIMULATION_TYPE == 3) then
        call transfer_b_fields_from_device(NDIM*NGLOB_AB,b_displ_elastic,b_veloc_elastic,b_accel_elastic,Mesh_pointer)
      endif
    endif
  endif

  ! noise movie
  if (NOISE_MOVIE_OUTPUT) then
    call noise_save_movie_output()
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
  if (output_wavefield_dumps) then
    call write_wavefield_dumps()
  endif

  end subroutine write_movie_output

