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
!========================================================================

  subroutine finalize_simulation()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par
  use specfem_par_gpu
  use specfem_par_movie,only: simulation_title,output_wavefield_dumps,mask_ibool

  implicit none

  integer :: i,ispec,j,iglob
  integer :: ier
  real(kind=4),dimension(:,:,:),allocatable :: rho_save, vp_save, vs_save, kappa_save, x_save, z_save
  double precision :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic
  character(len=MAX_STRING_LEN) :: inputname,outputname

  ! writes out kernel files
  if (SIMULATION_TYPE == 3) then
    call save_adjoint_kernels()
  endif

  ! saves model files
  if (trim(SAVE_MODEL) /= 'default') then
    allocate(rho_save(NGLLX,NGLLZ,nspec))
    allocate(vp_save(NGLLX,NGLLZ,nspec))
    allocate(vs_save(NGLLX,NGLLZ,nspec))
    allocate(kappa_save(NGLLX,NGLLZ,nspec))
    allocate(x_save(NGLLX,NGLLZ,nspec))
    allocate(z_save(NGLLX,NGLLZ,nspec))
    do ispec= 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          rho_save(i,j,ispec) = density(1,kmato(ispec))
          lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
          mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
          if (AXISYM) then ! CHECK kappa
            kappa_save(i,j,ispec) = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic/3._CUSTOM_REAL
          else
            kappa_save(i,j,ispec) = lambdal_unrelaxed_elastic + mul_unrelaxed_elastic
          endif
          vp_save(i,j,ispec) = sqrt((kappa_save(i,j,ispec) + FOUR_THIRDS *mul_unrelaxed_elastic)/density(1,kmato(ispec)))
          vs_save(i,j,ispec) = sqrt(mul_unrelaxed_elastic/density(1,kmato(ispec)))

          iglob = ibool(i,j,ispec)
          x_save(i,j,ispec) = coord(1,iglob)
          z_save(i,j,ispec) = coord(2,iglob)
        enddo
      enddo
    enddo

    if (trim(SAVE_MODEL) == 'legacy') then
      write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_model_velocity.dat_input'
      open(unit=1001,file=inputname,status='unknown',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_model_velocity.dat_input')
      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            write(1001,'(6e15.5e4)') x_save(i,j,ispec), x_save(i,j,ispec),z_save(i,j,ispec),rho_save(i,j,ispec),&
                                     vp_save(i,j,ispec),vs_save(i,j,ispec)
          enddo
        enddo
      enddo
      close(1001)

    else if (trim(SAVE_MODEL) == 'ascii') then
      write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_rho_vp_vs.dat'
      open(unit=1001,file= inputname,status='unknown',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_rho_vp_vs.dat')
      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            write(1001,'(5e15.5e4)') x_save(i,j,ispec),z_save(i,j,ispec),rho_save(i,j,ispec),vp_save(i,j,ispec),vs_save(i,j,ispec)
          enddo
        enddo
      enddo
      close(1001)

    else if ((trim(SAVE_MODEL) == 'binary') .or. (trim(SAVE_MODEL) == 'gll')) then
      write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_rho.bin'
      open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_rho.bin')
      write(172) rho_save
      close(172)
      write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vp.bin'
      open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_vp.bin')
      write(172) vp_save
      close(172)
      write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vs.bin'
      open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_vs.bin')
      write(172) vs_save
      close(172)
      write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_x.bin'
      open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_x.bin')
      write(172) x_save
      close(172)
      write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_z.bin'
      open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_z.bin')
      write(172) z_save
      close(172)

      write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_jacobian.bin'
      open(unit=172,file=outputname,status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_jacobian.bin')
      write(172) jacobian
      close(172)

    else
      stop 'Save Model not implemented for external and tomo'
    endif !Type of model
  endif !save model

  ! frees memory
  if (GPU_MODE) then
    ! frees temporary arrays
    if (any_elastic) deallocate(tmp_displ_2D,tmp_veloc_2D,tmp_accel_2D)

    deallocate(tab_requests_send_recv_scalar,b_tab_requests_send_recv_scalar)
    deallocate(tab_requests_send_recv_vector,b_tab_requests_send_recv_vector)
    deallocate(buffer_send_scalar_ext_mesh,b_buffer_send_scalar_ext_mesh)
    deallocate(buffer_recv_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh)
    deallocate(buffer_send_vector_ext_mesh,b_buffer_send_vector_ext_mesh)
    deallocate(buffer_recv_vector_ext_mesh,b_buffer_recv_vector_ext_mesh)

    ! frees memory on GPU
    call prepare_cleanup_device(Mesh_pointer,any_acoustic,any_elastic, &
                                STACEY_BOUNDARY_CONDITIONS, &
                                ANISOTROPY, &
                                APPROXIMATE_HESS_KL)
  endif

  if (output_wavefield_dumps) deallocate(mask_ibool)

  ! stores absorbing boundary contributions into files
  if (anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. STACEY_BOUNDARY_CONDITIONS) then

    if (any_acoustic) then
      !--- left absorbing boundary
      if (nspec_left > 0) write(65) b_absorb_acoustic_left
      !--- right absorbing boundary
      if (nspec_right > 0) write(66) b_absorb_acoustic_right
      !--- bottom absorbing boundary
      if (nspec_bottom > 0) write(67) b_absorb_acoustic_bottom
      !--- top absorbing boundary
      if (nspec_top > 0) write(68) b_absorb_acoustic_top
      ! closes absorbing files
      close(65)
      close(66)
      close(67)
      close(68)
    endif !any acoustic

    if (any_elastic) then
      !--- left absorbing boundary
      if (nspec_left > 0) write(35) b_absorb_elastic_left
      !--- right absorbing boundary
      if (nspec_right > 0)  write(36) b_absorb_elastic_right
      !--- bottom absorbing boundary
      if (nspec_bottom > 0)  write(37) b_absorb_elastic_bottom
      !--- top absorbing boundary
      if (nspec_top > 0) write(38) b_absorb_elastic_top
      ! closes absorbing files
      close(35)
      close(36)
      close(37)
      close(38)
    endif !any elastic

    if (any_poroelastic) then
      !--- left absorbing boundary
      if(nspec_left > 0) then
        write(45) b_absorb_poro_s_left
        write(25) b_absorb_poro_w_left
      endif
      !--- right absorbing boundary
      if(nspec_right > 0) then
        write(46) b_absorb_poro_s_right
        write(26) b_absorb_poro_w_right
      endif
      !--- bottom absorbing boundary
      if(nspec_bottom > 0)  then
        write(47) b_absorb_poro_s_bottom
        write(29) b_absorb_poro_w_bottom
      endif
      !--- top absorbing boundary
      if(nspec_top > 0) then
        write(48) b_absorb_poro_s_top
        write(28) b_absorb_poro_w_top
      endif
      ! closes absorbing files
      close(45)
      close(46)
      close(47)
      close(48)
      close(25)
      close(26)
      close(29)
      close(28)
    endif
  endif

  ! PML
  if (anyabs_glob .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. PML_BOUNDARY_CONDITIONS) then
    if (any_elastic .and. nglob_interface > 0) close(71)
    if (any_acoustic .and. nglob_interface > 0) close(72)
  endif

  ! save last frame
  if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. any_elastic) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Saving elastic last frame...'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_elastic**.bin')

    write(55) displ_elastic
    write(55) veloc_elastic
    write(55) accel_elastic

    close(55)
  endif

  if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. any_poroelastic) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Saving poroelastic last frame...'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_s',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_poroelastic_s**.bin')

    write(55) displs_poroelastic
    write(55) velocs_poroelastic
    write(55) accels_poroelastic
    close(55)

    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_w',myrank,'.bin'
    open(unit=56,file='OUTPUT_FILES/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_poroelastic_w**.bin')

    write(56) displw_poroelastic
    write(56) velocw_poroelastic
    write(56) accelw_poroelastic
    close(56)
  endif

  if (SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. any_acoustic) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Saving acoustic last frame...'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_acoustic**.bin')

    write(55) potential_acoustic
    write(55) potential_dot_acoustic
    write(55) potential_dot_dot_acoustic
    close(55)
  endif

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

  ! close energy file
  if (output_energy .and. myrank == 0) close(IOUT_ENERGY)

  ! print exit banner
  if (myrank == 0) call datim(simulation_title)

  ! close output file
  if (IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

  end subroutine finalize_simulation
