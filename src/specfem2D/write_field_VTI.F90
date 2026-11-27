!========================================================================
!
!                            S P E C F E M 2 D
!                            -----------------
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

  subroutine write_field_VTI(plot_b_wavefield_only,field_id)

  use constants, only: IMAIN,NGLLX,NGLLZ,REMOVE_PMLS_FROM_JPEG_IMAGES

  use specfem_par, only: myrank,nspec,it,it_begin,deltat,NSTEP,NPROC, &
                         ibool,rhostore,P_SV, &
                         potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                         displ_elastic,veloc_elastic,accel_elastic, &
                         displs_poroelastic,velocs_poroelastic,accels_poroelastic, &
                         displ_electromagnetic, &
                         b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                         b_displ_elastic,b_veloc_elastic,b_accel_elastic, &
                         b_displs_poroelastic,b_velocs_poroelastic,b_accels_poroelastic,SIMULATION_TYPE, &
                         UNDO_ATTENUATION_AND_OR_PML,NO_BACKWARD_RECONSTRUCTION

  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML
  use shared_parameters, only: NTSTEP_BETWEEN_OUTPUT_IMAGES

  use specfem_par_movie, only: vector_field_display,image_color_data,iglob_image_color, &
    nb_pixel_loc, &
    NX_IMAGE_color,NZ_IMAGE_color, &
    num_pixel_loc, &
    size_pixel_horizontal,size_pixel_vertical, &
    field_name

#ifdef WITH_MPI
  use specfem_par_movie, only: data_pixel_recv,data_pixel_send,nb_pixel_per_proc,num_pixel_recv
#endif

  ! helper adjoint wavefields
  use specfem_par, only: coupled_acoustic_elastic
  use specfem_par_movie, only: potential_adj_acoustic,potential_adj_int_acoustic,potential_adj_int_int_acoustic

  implicit none

  ! parameter useful for UNDO_ATTENUATION
  logical,intent(in) :: plot_b_wavefield_only
  integer,intent(in) :: field_id
  ! local variables
  integer :: i,j,k,ispec,iglob,iproc,i_field,n_fields
  double precision :: rhol

  if (SIMULATION_TYPE == 1 .or. &
      SIMULATION_TYPE == 2 .or. &
      UNDO_ATTENUATION_AND_OR_PML .or. &
      NO_BACKWARD_RECONSTRUCTION) then
    ! single field
    n_fields = 1
  else
    ! forward & backward field
    n_fields = 2
  endif

  ! loops over each wavefield
  do i_field = 1, n_fields

    if (myrank == 0) then
      write(IMAIN,*) 'Writing VTI of size ',NX_IMAGE_color,' x ',NZ_IMAGE_color,' for time step ',it
      call flush_IMAIN()
    endif

    if (field_id >= 1 .and. field_id <= 3) then
      if (i_field == 1 .and. .not. plot_b_wavefield_only) then
        if (myrank == 0 .and. SIMULATION_TYPE == 1 ) then
          write(IMAIN,*) 'drawing scalar image of the forward wavefield displacement...'
        else if (myrank == 0) then
          write(IMAIN,*) 'drawing scalar image of the adjoint wavefield displacement...'
        endif
        if (coupled_acoustic_elastic .and. SIMULATION_TYPE == 3) then
          ! uses integrated acoustic adjoint potential
          call compute_vector_whole_medium(potential_adj_int_int_acoustic,displ_elastic,displs_poroelastic)
        else
          ! forward/adjoint displacement wavefields
          call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic)
        endif
      else
        if (myrank == 0) write(IMAIN,*) 'drawing scalar image of the reconstructed forward wavefield displacement...'
        call compute_vector_whole_medium(b_potential_acoustic,b_displ_elastic,b_displs_poroelastic)
      endif

    else if (field_id >= 4 .and. field_id <= 6) then

      if (i_field == 1 .and. .not. plot_b_wavefield_only) then
        if (myrank == 0 .and. SIMULATION_TYPE == 1 ) then
          write(IMAIN,*) 'drawing scalar image of the forward wavefield velocity...'
        else if (myrank == 0) then
          write(IMAIN,*) 'drawing scalar image of the adjoint wavefield velocity...'
        endif
        if (coupled_acoustic_elastic .and. SIMULATION_TYPE == 3) then
          ! uses integrated acoustic adjoint potential
          call compute_vector_whole_medium(potential_adj_int_acoustic,veloc_elastic,velocs_poroelastic)
        else
          ! forward/adjoint velocity wavefields
          call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic)
        endif
      else
        if (myrank == 0) write(IMAIN,*) 'drawing scalar image of the reconstructed forward wavefield velocity...'
        call compute_vector_whole_medium(b_potential_dot_acoustic,b_veloc_elastic,b_velocs_poroelastic)
      endif

    else if (field_id >= 7 .and. field_id <= 9) then

      if (i_field == 1 .and. .not. plot_b_wavefield_only) then
        if (myrank == 0 .and. SIMULATION_TYPE == 1 ) then
          write(IMAIN,*) 'drawing scalar image of the forward wavefield acceleration...'
        else if (myrank == 0) then
          write(IMAIN,*) 'drawing scalar image of the adjoint wavefield acceleration...'
        endif
        if (coupled_acoustic_elastic .and. SIMULATION_TYPE == 3) then
          ! uses integrated acoustic adjoint potential
          call compute_vector_whole_medium(potential_adj_acoustic,accel_elastic,accels_poroelastic)
        else
          ! forward/adjoint acceleration wavefields
          call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,accels_poroelastic)
        endif
      else
        if (myrank == 0) write(IMAIN,*) 'drawing scalar image of the reconstructed forward wavefield acceleration...'
        call compute_vector_whole_medium(b_potential_dot_dot_acoustic,b_accel_elastic,b_accels_poroelastic)
      endif

    else if (field_id >= 11 .and. field_id <= 13) then
      ! allocation for normalized representation in JPEG image
      ! for an atmosphere model
      if (myrank == 0) write(IMAIN,*) 'drawing scalar image of part of normalized displacement vector...'
      ! todo: for coupled simulations, displacement in the adjoint acoustic wavefield would need to be adapted...
      call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic)

      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            rhol = rhostore(i,j,ispec)
            iglob = ibool(i,j,ispec)

            vector_field_display(1,iglob) = sqrt(rhol) * vector_field_display(1,iglob)
            vector_field_display(2,iglob) = sqrt(rhol) * vector_field_display(2,iglob)
          enddo
        enddo
      enddo

    else if (field_id >= 14 .and. field_id <= 16) then
      ! allocation for normalized representation in JPEG image
      ! for an atmosphere model
      ! todo: for coupled simulations, velocity in the adjoint acoustic wavefield would need to be adapted...
      call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic)

      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            rhol = rhostore(i,j,ispec)
            iglob = ibool(i,j,ispec)

            vector_field_display(1,iglob) = sqrt(rhol) * vector_field_display(1,iglob)
            vector_field_display(2,iglob) = sqrt(rhol) * vector_field_display(2,iglob)
          enddo
        enddo
      enddo

    else if (field_id >= 17 .and. field_id <= 19) then
    ! added entries for Ex, Ez, E norm
      if (myrank == 0) write(IMAIN,*) 'drawing scalar image of part of the Electric vector...'
      do ispec = 1,nspec
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            vector_field_display(1,iglob) = displ_electromagnetic(1,iglob)
            vector_field_display(2,iglob) = displ_electromagnetic(2,iglob)
          enddo
        enddo
      enddo

    else if (field_id == 10 .and. P_SV) then

      if (i_field == 1 .and. .not. plot_b_wavefield_only) then
        if (myrank == 0 .and. SIMULATION_TYPE == 1 ) then
          write(IMAIN,*) 'drawing scalar image of the forward wavefield pressure...'
        else
          if (myrank == 0) write(IMAIN,*) 'drawing scalar image of the adjoint wavefield pressure...'
        endif
        ! todo: for coupled simulations, pressure in the adjoint acoustic wavefield would need to be adapted...
        call compute_pressure_whole_medium(1)
      else
        if (myrank == 0) write(IMAIN,*) 'drawing scalar image of the reconstructed forward wavefield pressure...'
        call compute_pressure_whole_medium(2)
      endif

    else if (field_id == 10 .and. .not. P_SV) then
      call exit_MPI(myrank,'cannot draw pressure field for SH (membrane) waves')

    else
      call exit_MPI(myrank,'wrong type for JPEG snapshots')
    endif

    !! DK DK quick hack to remove the PMLs from JPEG images if needed: set the vector field to zero there
    if (PML_BOUNDARY_CONDITIONS .and. REMOVE_PMLS_FROM_JPEG_IMAGES) then
      do ispec = 1,nspec
        if (ispec_is_PML(ispec)) then
          do j = 1,NGLLZ
            do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              vector_field_display(1,iglob) = 0.d0
              vector_field_display(2,iglob) = 0.d0
            enddo
          enddo
        endif
      enddo
    endif
    !! DK DK quick hack to remove the PMLs from JPEG images if needed

    image_color_data(:,:) = 0.d0

    do k = 1, nb_pixel_loc
      j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
      i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

      ! avoid edge effects
      if (i < 1 ) i = 1
      if (j < 1 ) j = 1

      if (i > NX_IMAGE_color ) i = NX_IMAGE_color
      if (j > NZ_IMAGE_color ) j = NZ_IMAGE_color

      if (P_SV) then
        ! P-SV waves, plot a component of vector, its norm, or else pressure
        if (iglob_image_color(i,j) /= -1) then
          if (field_id == 1 .or. field_id == 4 .or. field_id == 7 .or. &
              field_id == 11 .or. field_id == 14 .or. field_id == 17) then
            ! draw the X component of the vector
            image_color_data(i,j) = vector_field_display(1,iglob_image_color(i,j))

          else if (field_id == 2 .or. field_id == 5 .or. field_id == 8 .or. &
                  field_id == 12 .or. field_id == 15 .or.  field_id == 18) then
            ! draw the Z component of the vector
            image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))

          else if (field_id == 3 .or. field_id == 6 .or. field_id == 9 .or. &
                  field_id == 13 .or. field_id == 16 .or.  field_id == 19) then
            ! draw the norm of the vector
            image_color_data(i,j) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2  &
                                       + vector_field_display(2,iglob_image_color(i,j))**2)

          else if (field_id == 10) then
            ! by convention we have stored pressure in the 2. component of the array
            image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))

          else
            call exit_MPI(myrank,'wrong type for JPEG snapshots')
          endif
        endif

      else
        ! SH (membrane) waves, plot y-component
        if (iglob_image_color(i,j) /= -1) image_color_data(i,j) = vector_field_display(1,iglob_image_color(i,j))
      endif
    enddo

    ! assembling array image_color_data on process zero for color output
#ifdef WITH_MPI
    if (NPROC > 1) then
      if (myrank == 0) then
        do iproc = 1, NPROC-1
          call recv_dp(data_pixel_recv(1), nb_pixel_per_proc(iproc), iproc, 43)

          do k = 1, nb_pixel_per_proc(iproc)
            j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
            i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color

            ! avoid edge effects
            if (i < 1) i = 1
            if (j < 1) j = 1

            if (i > NX_IMAGE_color) i = NX_IMAGE_color
            if (j > NZ_IMAGE_color) j = NZ_IMAGE_color

            image_color_data(i,j) = data_pixel_recv(k)
          enddo
        enddo
      else
        do k = 1, nb_pixel_loc
          j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
          i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

          ! avoid edge effects
          if (i < 1) i = 1
          if (j < 1) j = 1

          if (i > NX_IMAGE_color) i = NX_IMAGE_color
          if (j > NZ_IMAGE_color) j = NZ_IMAGE_color

          if (P_SV) then ! P-SV waves, plot a component of vector, its norm, or else pressure

            if (field_id == 1 .or. field_id == 4 .or. field_id == 7 .or. &
                field_id == 11 .or. field_id == 14 .or. field_id == 17) then
               data_pixel_send(k) = vector_field_display(1,iglob_image_color(i,j))  ! draw the X component of the vector

            else if (field_id == 2 .or. field_id == 5 .or. field_id == 8 .or. &
                    field_id == 12 .or. field_id == 15 .or. field_id == 18) then
               data_pixel_send(k) = vector_field_display(2,iglob_image_color(i,j))  ! draw the Z component of the vector

            else if (field_id == 3 .or. field_id == 6 .or. field_id == 9 .or. &
                    field_id == 13 .or. field_id == 16 .or. field_id == 19) then
              data_pixel_send(k) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2 + &
                                        vector_field_display(2,iglob_image_color(i,j))**2)  ! draw the norm of the vector

            else if (field_id == 10) then
              ! by convention we have stored pressure in the 2. component of the array
              data_pixel_send(k) = vector_field_display(2,iglob_image_color(i,j))

            else
              call exit_MPI(myrank,'wrong field type for VTI snapshots')
            endif

          else ! SH (membrane) waves, plot y-component
            if (iglob_image_color(i,j) /= -1) data_pixel_send(k) = vector_field_display(1,iglob_image_color(i,j))
          endif
        enddo

        call send_dp(data_pixel_send(1), nb_pixel_loc, 0, 43)

      endif
    endif
    call synchronize_all()
#else
    ! dummy to avoid compiler warning
    iproc = NPROC
#endif

    call synchronize_all()

    ! creates image
    if (myrank == 0) then
      if(it==it_begin)then
        call write_vtk_imagec(NTSTEP_BETWEEN_OUTPUT_IMAGES,NX_IMAGE_color,     &
        NZ_IMAGE_color,real(size_pixel_horizontal,kind=4),                     &
        real(size_pixel_vertical,kind=4),0,real(deltat,kind=4),nstep,          &
        real(0.d0*image_color_data,kind=4),trim(field_name(field_id))//char(0))
      end if
      call write_vtk_imagec(NTSTEP_BETWEEN_OUTPUT_IMAGES,NX_IMAGE_color,       &
      NZ_IMAGE_color,real(size_pixel_horizontal,kind=4),                       &
      real(size_pixel_vertical,kind=4),it,real(deltat,kind=4),nstep,           &
      real(image_color_data),trim(field_name(field_id))//char(0))
      ! user output
      write(IMAIN,*) 'VTI data written'
      call flush_IMAIN()
    endif
    call synchronize_all()

  enddo !loop over wavefields (forward and adjoint if necessary)

  ! output
  if (myrank == 0) then
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine write_field_VTI

