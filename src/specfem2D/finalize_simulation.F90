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

  subroutine finalize_simulation()

  use constants, only: IMAIN,IOUT_ENERGY,ISTANDARD_OUTPUT

  use specfem_par
  use specfem_par_noise
  use specfem_par_gpu
  use specfem_par_movie, only: simulation_title,mask_ibool

  implicit none

  ! writes out kernel files
  if (SIMULATION_TYPE == 3) then
    call save_adjoint_kernels()
  endif

  ! saves model files
  if (trim(SAVE_MODEL) /= 'default' .and. trim(SAVE_MODEL) /= '.false.') then
    call save_model_files()
  endif

  ! For this mode, the forward model has been saved differently
  if ((.not. NO_BACKWARD_RECONSTRUCTION) ) then
    ! stores final wavefields into files
    if (SAVE_FORWARD .and. SIMULATION_TYPE == 1) then
      call save_forward_arrays()
    endif
  endif

  ! frees memory
  if (GPU_MODE) then
    ! frees arrays
    deallocate(request_send_recv_scalar_gpu,b_request_send_recv_scalar_gpu)
    deallocate(request_send_recv_vector_gpu,b_request_send_recv_vector_gpu)
    deallocate(buffer_send_scalar_gpu,b_buffer_send_scalar_gpu)
    deallocate(buffer_recv_scalar_gpu,b_buffer_recv_scalar_gpu)
    deallocate(buffer_send_vector_gpu,b_buffer_send_vector_gpu)
    deallocate(buffer_recv_vector_gpu,b_buffer_recv_vector_gpu)

    ! frees memory on GPU
    call prepare_cleanup_device(Mesh_pointer, &
                                any_acoustic,any_elastic,any_anisotropy, &
                                APPROXIMATE_HESS_KL, &
                                ATTENUATION_VISCOACOUSTIC, &
                                ATTENUATION_VISCOELASTIC, &
                                NO_BACKWARD_RECONSTRUCTION, &
                                no_backward_acoustic_buffer)
  endif

  if (output_wavefield_dumps) deallocate(mask_ibool)

  if (initialfield .and. over_critical_angle) then
    deallocate(v0x_left)
    deallocate(v0z_left)
    deallocate(t0x_left)
    deallocate(t0z_left)

    deallocate(v0x_right)
    deallocate(v0z_right)
    deallocate(t0x_right)
    deallocate(t0z_right)

    deallocate(v0x_bot)
    deallocate(v0z_bot)
    deallocate(t0x_bot)
    deallocate(t0z_bot)
  endif

  ! frees arrays (not complete, but at least a few...)
  deallocate(sisux,sisuz,siscurl)

  ! wavefields
  deallocate(displ_elastic,veloc_elastic,accel_elastic)
  deallocate(displ_elastic_old)
  deallocate(rmass_inverse_elastic)
  deallocate(b_displ_elastic,b_veloc_elastic,b_accel_elastic)
  deallocate(b_displ_elastic_old)
  deallocate(potential_acoustic,potential_acoustic_old)
  deallocate(potential_dot_acoustic,potential_dot_dot_acoustic)
  deallocate(rmass_inverse_acoustic)
  deallocate(b_potential_acoustic,b_potential_acoustic_old)
  if (.not. NO_BACKWARD_RECONSTRUCTION) then
    deallocate(b_potential_dot_acoustic,b_potential_dot_dot_acoustic)
  endif

  ! noise
  if (allocated(noise_sourcearray)) deallocate(noise_sourcearray)
  if (allocated(mask_noise)) deallocate(mask_noise)
  if (allocated(noise_surface_movie_y_or_z)) deallocate(noise_surface_movie_y_or_z)

  ! material
  deallocate(kappastore,mustore,rhostore,rho_vpstore,rho_vsstore)

  ! attenuation
  deallocate(e1,e11,e13)
  deallocate(dux_dxl_old,duz_dzl_old,dux_dzl_plus_duz_dxl_old)
  deallocate(A_newmark_nu1,B_newmark_nu1,A_newmark_nu2,B_newmark_nu2)
  deallocate(e1_acous_sf,sum_forces_old)
  if (SIMULATION_TYPE == 3) then
    if (any_acoustic) then
      deallocate(b_e1_acous_sf,b_sum_forces_old)
    endif
    if (any_elastic) then
      deallocate(b_e1,b_e11,b_e13,b_dux_dxl_old,b_duz_dzl_old,b_dux_dzl_plus_duz_dxl_old)
    endif
  endif

  ! close energy file
  if (OUTPUT_ENERGY .and. myrank == 0) close(IOUT_ENERGY)

  ! print exit banner
  if (myrank == 0) call datim(simulation_title)

  ! close output file
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

  end subroutine finalize_simulation
