
subroutine finalize_simulation()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

integer i,ispec,j,iglob

#ifdef USE_MPI
  include "precision.h"
#endif


  real(kind=4),dimension(:,:,:),allocatable :: rho_save, vp_save, vs_save, kappa_save, x_save, z_save

if ( trim(SAVE_MODEL) /= 'default' ) then
   allocate(rho_save(NGLLX,NGLLZ,nspec))
   allocate(vp_save(NGLLX,NGLLZ,nspec))
   allocate(vs_save(NGLLX,NGLLZ,nspec))
   allocate(kappa_save(NGLLX,NGLLZ,nspec))
   allocate(x_save(NGLLX,NGLLZ,nspec))
   allocate(z_save(NGLLX,NGLLZ,nspec))
do ispec=1,nspec
          do j = 1,NGLLZ
              do i = 1,NGLLX

              rho_save(i,j,ispec)            = density(1,kmato(ispec))
              lambdal_unrelaxed_elastic      = poroelastcoef(1,1,kmato(ispec))
              mul_unrelaxed_elastic          = poroelastcoef(2,1,kmato(ispec))
              kappa_save(i,j,ispec)          = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic/3._CUSTOM_REAL
              vp_save(i,j,ispec)             = sqrt((kappa_save(i,j,ispec) + &
                                                4._CUSTOM_REAL*mul_unrelaxed_elastic/ &
                                                3._CUSTOM_REAL)/density(1,kmato(ispec)))
              vs_save(i,j,ispec)             = sqrt(mul_unrelaxed_elastic/density(1,kmato(ispec)))
              iglob = ibool(i,j,ispec)
              x_save(i,j,ispec)              = coord(1,iglob)
              z_save(i,j,ispec)              = coord(2,iglob)
              enddo
        enddo
enddo


  if(trim(SAVE_MODEL) == 'legacy') then

    write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_model_velocity.dat_input'
    open(unit=1001,file=inputname,status='unknown')
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


  else if(trim(SAVE_MODEL)=='ascii') then

    write(inputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_rho_vp_vs.dat'
    open(unit=1001,file= inputname,status='unknown')
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          write(1001,'(5e15.5e4)') x_save(i,j,ispec),z_save(i,j,ispec),rho_save(i,j,ispec),vp_save(i,j,ispec),vs_save(i,j,ispec)
        enddo
      enddo
    enddo
    close(1001)

  else if((trim(SAVE_MODEL) == 'binary') .or. (trim(SAVE_MODEL) == 'gll')) then

          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_rho.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) rho_save
          close(172)
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vp.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) vp_save
          close(172)
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_vs.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) vs_save
          close(172)
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_x.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) x_save
          close(172)
          write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_z.bin'
          open(unit=172,file=outputname,status='unknown',form='unformatted')
          write(172) z_save
          close(172)

  else
       stop 'Save Model not implemented for external and tomo'

  endif !Type of model


endif !save model




if (GPU_MODE) call prepare_cleanup_device(Mesh_pointer, &
                              any_acoustic,any_elastic, &
                              STACEY_BOUNDARY_CONDITIONS, &
                              ANISOTROPY, &
                              APPROXIMATE_HESS_KL)


  if(output_wavefield_dumps) deallocate(mask_ibool)


!!!! Displacement Etienne GPU

! stores absorbing boundary contributions into files
      if(anyabs .and. SAVE_FORWARD .and. SIMULATION_TYPE == 1 .and. (.not. PML_BOUNDARY_CONDITIONS)) then

      if (any_acoustic) then

        !--- left absorbing boundary
        if(nspec_left >0) write(65) b_absorb_acoustic_left
        !--- right absorbing boundary
        if(nspec_right >0) write(66) b_absorb_acoustic_right
        !--- bottom absorbing boundary
        if(nspec_bottom >0) write(67) b_absorb_acoustic_bottom
        !--- top absorbing boundary
        if(nspec_top >0) write(68) b_absorb_acoustic_top

      endif !any acoustic

      close(65)
      close(66)
      close(67)
      close(68)
      close(72)

 if(any_elastic) then

        !--- left absorbing boundary
        if(nspec_left >0) write(35) b_absorb_elastic_left
        !--- right absorbing boundary
        if(nspec_right >0)  write(36) b_absorb_elastic_right
        !--- bottom absorbing boundary
        if(nspec_bottom >0)  write(37) b_absorb_elastic_bottom
        !--- top absorbing boundary
        if(nspec_top >0) write(38) b_absorb_elastic_top

   endif !any elastic

      close(35)
      close(36)
      close(37)
      close(38)
      close(71)


    if(any_poroelastic) then
      close(25)
      close(45)
      close(26)
      close(46)
      close(29)
      close(47)
      close(28)
      close(48)
    endif

  endif

!
!--- save last frame
!
  if(SAVE_FORWARD .and. SIMULATION_TYPE ==1 .and. any_elastic) then
    if ( myrank == 0 ) then
      write(IOUT,*)
      write(IOUT,*) 'Saving elastic last frame...'
      write(IOUT,*)
    endif
    write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')

        write(55) displ_elastic
        write(55) veloc_elastic
        write(55) accel_elastic

    close(55)
  endif

  if(SAVE_FORWARD .and. SIMULATION_TYPE ==1 .and. any_poroelastic) then
    if ( myrank == 0 ) then
      write(IOUT,*)
      write(IOUT,*) 'Saving poroelastic last frame...'
      write(IOUT,*)
    endif
    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_s',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_w',myrank,'.bin'
    open(unit=56,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')

      write(55) displs_poroelastic
      write(55) velocs_poroelastic
      write(55) accels_poroelastic
      write(56) displw_poroelastic
      write(56) velocw_poroelastic
      write(56) accelw_poroelastic

    close(55)
    close(56)
  endif

  if(SAVE_FORWARD .and. SIMULATION_TYPE ==1 .and. any_acoustic) then
    if ( myrank == 0 ) then
      write(IOUT,*)
      write(IOUT,*) 'Saving acoustic last frame...'
      write(IOUT,*)
    endif
    write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
    open(unit=55,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
      write(55) potential_acoustic
      write(55) potential_dot_acoustic
      write(55) potential_dot_dot_acoustic
    close(55)
  endif


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

!----  close energy file
  if(output_energy .and. myrank == 0) close(IOUT_ENERGY)

  if (OUTPUT_MODEL_VELOCITY_FILE .and. .not. any_poroelastic) then
    write(outputname,'(a,i6.6,a)') 'DATA/proc',myrank,'_rho_vp_vs.dat_output'
    open(unit=1001,file=outputname,status='unknown')
    if ( .NOT. assign_external_model) then
      allocate(rho_local(ngllx,ngllz,nspec)); rho_local=0.
      allocate(vp_local(ngllx,ngllz,nspec)); vp_local=0.
      allocate(vs_local(ngllx,ngllz,nspec)); vs_local=0.
      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            rho_local(i,j,ispec) = density(1,kmato(ispec))
            vp_local(i,j,ispec) = sqrt(poroelastcoef(3,1,kmato(ispec))/density(1,kmato(ispec)))
            vs_local(i,j,ispec) = sqrt(poroelastcoef(2,1,kmato(ispec))/density(1,kmato(ispec)))
            write(1001,'(I10, 5F13.4)') iglob, coord(1,iglob),coord(2,iglob),&
                                      rho_local(i,j,ispec),vp_local(i,j,ispec),vs_local(i,j,ispec)
          enddo
        enddo
      enddo
    else
      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            write(1001,'(I10,5F13.4)') iglob, coord(1,iglob),coord(2,iglob),&
                                       rhoext(i,j,ispec),vpext(i,j,ispec),vsext(i,j,ispec)
          enddo
        enddo
      enddo
    endif
    close(1001)
  endif

! print exit banner
  if (myrank == 0) call datim(simulation_title)

!
!----  close output file
!
  if(IOUT /= ISTANDARD_OUTPUT) close(IOUT)

!
!----  end MPI
!
#ifdef USE_MPI
  call MPI_FINALIZE(ier)
#endif


end subroutine finalize_simulation
