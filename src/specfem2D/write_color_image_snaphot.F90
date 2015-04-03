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
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================
  subroutine write_color_image_snaphot()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

  !local variables
  integer :: i,j,k,ispec,iglob

  if( myrank == 0 ) then
    write(IOUT,*)
    write(IOUT,*) 'Creating color image of size ',NX_IMAGE_color,' x ',NZ_IMAGE_color,' for time step ',it
  endif

  if( imagetype_JPEG >= 1 .and. imagetype_JPEG <= 3 ) then
    if( myrank == 0 ) write(IOUT,*) 'drawing scalar image of part of the displacement vector...'
    call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                                     potential_gravito,displ_elastic,displs_poroelastic)

  else if( imagetype_JPEG >= 4 .and. imagetype_JPEG <= 6 ) then
    if( myrank == 0 ) write(IOUT,*) 'drawing scalar image of part of the velocity vector...'
    call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
              potential_dot_gravito,veloc_elastic,velocs_poroelastic)

  else if( imagetype_JPEG >= 7 .and. imagetype_JPEG <= 9 ) then
    if( myrank == 0 ) write(IOUT,*) 'drawing scalar image of part of the acceleration vector...'
    call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic, &
              potential_dot_dot_gravito,accel_elastic,accels_poroelastic)

  else if( imagetype_JPEG >= 11 .and. imagetype_JPEG <= 13 ) then
    ! allocation for normalized representation in JPEG image
    ! for an atmosphere model
    if( myrank == 0 ) write(IOUT,*) 'drawing scalar image of part of normalized displacement vector...'
    call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
              potential_gravito,displ_elastic,displs_poroelastic)

    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if(  assign_external_model ) then
            rhol = rhoext(i,j,ispec)
          else
            rhol = density(1,kmato(ispec))
          endif
          iglob = ibool(i,j,ispec)
          vector_field_display(1,iglob) = sqrt(rhol) * vector_field_display(1,iglob)
          vector_field_display(2,iglob) = sqrt(rhol) * vector_field_display(2,iglob)
          vector_field_display(3,iglob) = sqrt(rhol) * vector_field_display(3,iglob)
        enddo
      enddo
    enddo

  else if( imagetype_JPEG >= 14 .and. imagetype_JPEG <= 16 ) then
    ! allocation for normalized representation in JPEG image
    ! for an atmosphere model
    call compute_vector_whole_medium(potential_dot_acoustic,potential_dot_gravitoacoustic, &
              potential_dot_gravito,veloc_elastic,velocs_poroelastic)

    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if(  assign_external_model ) then
            rhol = rhoext(i,j,ispec)
          else
            rhol = density(1,kmato(ispec))
          endif
          iglob = ibool(i,j,ispec)
          vector_field_display(1,iglob) = sqrt(rhol) * vector_field_display(1,iglob)
          vector_field_display(2,iglob) = sqrt(rhol) * vector_field_display(2,iglob)
          vector_field_display(3,iglob) = sqrt(rhol) * vector_field_display(3,iglob)
        enddo
      enddo
    enddo

  else if( imagetype_JPEG == 10 .and. p_sv ) then
    if( myrank == 0 ) write(IOUT,*) 'drawing image of pressure field...'
    call compute_pressure_whole_medium()

  else if( imagetype_JPEG == 10 .and. .not. p_sv ) then
    call exit_MPI('cannot draw pressure field for SH (membrane) waves')

  else
    call exit_MPI('wrong type for JPEG snapshots')
  endif

!! DK DK quick hack to remove the PMLs from JPEG images if needed: set the vector field to zero there
  if( PML_BOUNDARY_CONDITIONS .and. REMOVE_PMLS_FROM_JPEG_IMAGES ) then
    do ispec = 1,nspec
      if( is_PML(ispec) ) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            vector_field_display(1,iglob) = 0.d0
            vector_field_display(2,iglob) = 0.d0
            vector_field_display(3,iglob) = 0.d0
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
    if( i < 1 ) i = 1
    if( j < 1 ) j = 1

    if( i > NX_IMAGE_color ) i = NX_IMAGE_color
    if( j > NZ_IMAGE_color ) j = NZ_IMAGE_color

    if( p_sv ) then ! P-SH waves, plot a component of vector, its norm, or else pressure
      if( iglob_image_color(i,j) /= -1 ) then
        if( imagetype_JPEG == 1  .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7 .or. &
            imagetype_JPEG == 11 .or. imagetype_JPEG == 14 ) then
          image_color_data(i,j) = vector_field_display(1,iglob_image_color(i,j))  ! draw the X component of the vector

        else if( imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8 .or. &
                imagetype_JPEG == 12 .or. imagetype_JPEG == 15 ) then
          image_color_data(i,j) = vector_field_display(3,iglob_image_color(i,j))  ! draw the Z component of the vector

        else if( imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9 .or. &
                imagetype_JPEG == 13 .or. imagetype_JPEG == 16 ) then
          image_color_data(i,j) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2 + &
                           vector_field_display(3,iglob_image_color(i,j))**2)  ! draw the norm of the vector

        else if( imagetype_JPEG == 10 ) then
! by convention we have stored pressure in the third component of the array
          image_color_data(i,j) = vector_field_display(3,iglob_image_color(i,j))

        else
          call exit_MPI('wrong type for JPEG snapshots')
        endif
      endif

    else ! SH (membrane) waves, plot y-component
      if( iglob_image_color(i,j) /= -1) image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))
    endif
  enddo

! assembling array image_color_data on process zero for color output
#ifdef USE_MPI

  if( nproc > 1 ) then
    if( myrank == 0 ) then
      do iproc = 1, nproc-1
        call MPI_RECV(data_pixel_recv(1),nb_pixel_per_proc(iproc+1), MPI_DOUBLE_PRECISION, &
                      iproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        do k = 1, nb_pixel_per_proc(iproc+1)
          j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
          i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color

          ! avoid edge effects
          if( i < 1) i = 1
          if( j < 1) j = 1

          if( i > NX_IMAGE_color) i = NX_IMAGE_color
          if( j > NZ_IMAGE_color) j = NZ_IMAGE_color

          image_color_data(i,j) = data_pixel_recv(k)
        enddo
      enddo
    else
      do k = 1, nb_pixel_loc
        j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
        i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

        ! avoid edge effects
        if( i < 1) i = 1
        if( j < 1) j = 1

        if( i > NX_IMAGE_color) i = NX_IMAGE_color
        if( j > NZ_IMAGE_color) j = NZ_IMAGE_color

        if( p_sv ) then ! P-SH waves, plot a component of vector, its norm, or else pressure

          if( imagetype_JPEG == 1 .or. imagetype_JPEG == 4 .or. imagetype_JPEG == 7 .or. &
              imagetype_JPEG == 11 .or. imagetype_JPEG == 14 ) then
             data_pixel_send(k) = vector_field_display(1,iglob_image_color(i,j))  ! draw the X component of the vector

          else if( imagetype_JPEG == 2 .or. imagetype_JPEG == 5 .or. imagetype_JPEG == 8 .or. &
                  imagetype_JPEG == 12 .or. imagetype_JPEG == 15 ) then
             data_pixel_send(k) = vector_field_display(3,iglob_image_color(i,j))  ! draw the Z component of the vector

          else if( imagetype_JPEG == 3 .or. imagetype_JPEG == 6 .or. imagetype_JPEG == 9 .or. &
                  imagetype_JPEG == 13 .or. imagetype_JPEG == 16 ) then
            data_pixel_send(k) = sqrt(vector_field_display(1,iglob_image_color(i,j))**2 + &
                                      vector_field_display(3,iglob_image_color(i,j))**2)  ! draw the norm of the vector

          else if( imagetype_JPEG == 10 ) then
            ! by convention we have stored pressure in the third component of the array
            data_pixel_send(k) = vector_field_display(3,iglob_image_color(i,j))

          else
            call exit_MPI('wrong type for JPEG snapshots')
          endif

        else ! SH (membrane) waves, plot y-component
          if( iglob_image_color(i,j) /= -1) data_pixel_send(k) = vector_field_display(2,iglob_image_color(i,j))
        endif
      enddo
      call MPI_SEND(data_pixel_send(1),nb_pixel_loc,MPI_DOUBLE_PRECISION, 0, 43, MPI_COMM_WORLD, ier)
    endif
  endif
#endif

  if( myrank == 0 ) then
    call create_color_image()
    write(IOUT,*) 'Color image created'
  endif

  end subroutine write_color_image_snaphot

