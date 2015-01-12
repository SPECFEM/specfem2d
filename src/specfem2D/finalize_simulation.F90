
subroutine finalize_simulation()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

#ifdef USE_MPI
  include "precision.h"
#endif

  if(output_wavefield_dumps) deallocate(mask_ibool)

  if((SAVE_FORWARD .and. SIMULATION_TYPE==1) .or. SIMULATION_TYPE == 3) then
    if(any_acoustic) then
      close(65)
      close(66)
      close(67)
      close(68)
      close(72)
    endif
    if(any_elastic) then
      close(35)
      close(36)
      close(37)
      close(38)
      close(71)

    endif
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
    if(p_sv)then !P-SV waves
      do j=1,nglob
        write(55) displ_elastic(1,j), displ_elastic(3,j), &
                  veloc_elastic(1,j), veloc_elastic(3,j), &
                  accel_elastic(1,j), accel_elastic(3,j)
      enddo
    else !SH (membrane) waves
      do j=1,nglob
        write(55) displ_elastic(2,j), &
                  veloc_elastic(2,j), &
                  accel_elastic(2,j)
      enddo
    endif
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
       do j=1,nglob
      write(55) (displs_poroelastic(i,j), i=1,NDIM), &
                  (velocs_poroelastic(i,j), i=1,NDIM), &
                  (accels_poroelastic(i,j), i=1,NDIM)
      write(56) (displw_poroelastic(i,j), i=1,NDIM), &
                  (velocw_poroelastic(i,j), i=1,NDIM), &
                  (accelw_poroelastic(i,j), i=1,NDIM)
       enddo
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
       do j=1,nglob
      write(55) potential_acoustic(j),&
               potential_dot_acoustic(j),&
               potential_dot_dot_acoustic(j)
       enddo
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
    open(unit=1001,file='DATA/model_velocity.dat_output',status='unknown')
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
