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

  subroutine prepare_timerun()

  use constants, only: USE_ENFORCE_FIELDS,IOUT_ENERGY,IMAIN
  use specfem_par
  use specfem_par_movie
  use specfem_par_noise, only: NOISE_TOMOGRAPHY

  implicit none

  ! Test compatibility with axisymmetric formulation
  if (AXISYM) call check_compatibility_axisym()

  ! prepares constant factors for time scheme and seismograms
  call prepare_timerun_constants()

  ! wavefield array initialization
  call prepare_wavefields()

  ! PML preparation
  call prepare_PML()

  ! prepares mass matrices
  call prepare_timerun_mass_matrix()

  ! check which GLL will be forced or not
  if (USE_ENFORCE_FIELDS) call build_forced()

  ! postscript images for grids and snapshots
  call prepare_timerun_postscripts()

  ! jpeg images
  call prepare_timerun_image_coloring()

  ! for adjoint kernel runs
  call prepare_timerun_adjoint()

  ! reads initial fields from external file if needed
  call prepare_timerun_initialfield

  ! compute the source time function and stores it in a text file
  call prepare_source_time_function()

  ! prepares noise simulations
  if (NOISE_TOMOGRAPHY /= 0) call prepare_timerun_noise()

  ! attenuation
  call prepare_timerun_attenuation()

  ! prepares GPU arrays
  if (GPU_MODE) call prepare_GPU()

  !-------------------------------------------------------------

  ! creates a Gnuplot script to display the energy curve in log scale
  if (output_energy .and. myrank == 0) then
    close(IOUT_ENERGY)
    open(unit=IOUT_ENERGY,file='OUTPUT_FILES/plot_energy.gnu',status='unknown',action='write')
    write(IOUT_ENERGY,*) 'set term wxt'
    write(IOUT_ENERGY,*) '#set term postscript landscape color solid "Helvetica" 22'
    write(IOUT_ENERGY,*) '#set output "energy.ps"'
    write(IOUT_ENERGY,*) '# set xrange [0:60]'
    write(IOUT_ENERGY,*) 'set logscale y'
    write(IOUT_ENERGY,*) 'set xlabel "Time (s)"'
    write(IOUT_ENERGY,*) 'set ylabel "Energy (J)"'
    write(IOUT_ENERGY,*) 'set loadpath "./OUTPUT_FILES"'
    write(IOUT_ENERGY,'(A)') &
      'plot "energy.dat" us 1:4 t ''Total Energy'' w l lc 1, "energy.dat" us 1:3 t ''Potential Energy'' w l lc 2'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    close(IOUT_ENERGY)
  endif

  ! open the file in which we will store the energy curve
  if (output_energy .and. myrank == 0) open(unit=IOUT_ENERGY,file='OUTPUT_FILES/energy.dat',status='unknown',action='write')

  ! synchronizes all processes
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "done, preparation successful"
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun


!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_constants()

  use constants, only: HALF,ZERO
  use specfem_par

  implicit none

  ! local parameters
  integer :: ier

  ! defines coefficients of the Newmark time scheme
  deltatover2 = HALF * deltat
  deltatsquareover2 = HALF * deltat * deltat

  !  define coefficients of the Newmark time scheme for the backward wavefield
  if (SIMULATION_TYPE == 3) then
    if (UNDO_ATTENUATION) then
      ! moves forward
      b_deltat = deltat
      b_deltatover2 = deltatover2
      b_deltatsquareover2 = deltatsquareover2
    else
      ! reconstructed wavefield moves backward in time from last snapshot
      b_deltat = - deltat
      b_deltatover2 = HALF * b_deltat
      b_deltatsquareover2 = HALF * b_deltat * b_deltat
    endif
  else
    ! will not be used, but initialized
    b_deltat = 0._CUSTOM_REAL
    b_deltatover2 = 0._CUSTOM_REAL
    b_deltatsquareover2 = 0._CUSTOM_REAL
  endif

  ! seismograms
  ! allocate seismogram arrays
  allocate(sisux(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrecloc), &
           sisuz(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrecloc), &
           siscurl(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrecloc),stat=ier)
  if (ier /= 0) stop 'Error allocating seismogram arrays'

  sisux(:,:) = ZERO ! double precision zero
  sisuz(:,:) = ZERO
  siscurl(:,:) = ZERO

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_constants


!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_mass_matrix()

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IMAIN
  use specfem_par

  implicit none

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Preparing mass matrices'
    call flush_IMAIN()
  endif

  ! builds the global mass matrix
  call invert_mass_matrix_init()

#ifdef USE_MPI
  ! assembling the mass matrix of shared nodes on MPI partition interfaces
  call assemble_MPI_scalar(rmass_inverse_acoustic,nglob_acoustic, &
                           rmass_inverse_elastic,nglob_elastic, &
                           rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic,nglob_poroelastic)
#endif

  ! inverts mass matrix
  ! note: we will divide the acceleration terms by this diagonal mass matrix (which is trivial to invert)
  call invert_mass_matrix()

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_mass_matrix

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_postscripts()

  use constants, only: IMAIN
  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: i,j,ier

  ! for Lagrange interpolants
  double precision, external :: hgll, hglj

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Preparing image coloring'
    call flush_IMAIN()
  endif

  ! arrays for display images
  allocate(shape2D_display(ngnod,pointsdisp,pointsdisp), &
           dershape2D_display(NDIM,ngnod,pointsdisp,pointsdisp),stat=ier)
  if (ier /= 0) stop 'Error allocating shape arrays for display'

  ! computes shape functions and their derivatives for regular interpolated display grid
  do j = 1,pointsdisp
    do i = 1,pointsdisp
      xirec  = 2.d0*dble(i-1)/dble(pointsdisp-1) - 1.d0
      gammarec  = 2.d0*dble(j-1)/dble(pointsdisp-1) - 1.d0
      call define_shape_functions(shape2D_display(:,i,j),dershape2D_display(:,:,i,j),xirec,gammarec,ngnod)
    enddo
  enddo

  ! for postscript snapshots
  ! arrays for display images as snapshot postscript images
  allocate(flagrange(NGLLX,pointsdisp))
  if (AXISYM) then
    allocate(flagrange_GLJ(NGLJ,pointsdisp))
  else
    allocate(flagrange_GLJ(1,1))
  endif

  ! compute Lagrange interpolants on a regular interpolated grid in (xi,gamma)
  ! for display (assumes NGLLX = NGLLZ)
  do j = 1,NGLLX
    do i = 1,pointsdisp
      xirec  = 2.d0*dble(i-1)/dble(pointsdisp-1) - 1.d0
      flagrange(j,i) = hgll(j-1,xirec,xigll,NGLLX)
      if (AXISYM) flagrange_GLJ(j,i) = hglj(j-1,xirec,xiglj,NGLJ)
    enddo
  enddo

  allocate(xinterp(pointsdisp,pointsdisp))
  allocate(zinterp(pointsdisp,pointsdisp))
  allocate(Uxinterp(pointsdisp,pointsdisp))
  allocate(Uzinterp(pointsdisp,pointsdisp))

  ! to display the whole vector field (it needs to be computed from the potential in acoustic elements,
  ! thus it does not exist as a whole in case of simulations that contain some acoustic elements
  ! and it thus needs to be computed specifically for display purposes)
  allocate(vector_field_display(NDIM,nglob))

  ! when periodic boundary conditions are on, some global degrees of freedom are going to be removed,
  ! thus we need to set this array to zero otherwise some of its locations may contain random values
  ! if the memory is not cleaned
  vector_field_display(:,:) = 0.d0

  ! check the mesh, stability and number of points per wavelength
  call check_grid()

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_postscripts


!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_image_coloring()

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IMAIN
#ifdef USE_MPI
  use constants, only: DISPLAY_COLORS,DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT
#endif
  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: i,j,k,iproc,ipixel
  integer :: ier

  ! postscripts
  integer :: d1_coorg_send_ps_velocity_model,d2_coorg_send_ps_velocity_model, &
             d1_coorg_recv_ps_velocity_model,d2_coorg_recv_ps_velocity_model, &
             d1_RGB_send_ps_velocity_model,d2_RGB_send_ps_velocity_model, &
             d1_RGB_recv_ps_velocity_model,d2_RGB_recv_ps_velocity_model

  integer :: d1_coorg_send_ps_element_mesh,d2_coorg_send_ps_element_mesh, &
             d1_coorg_recv_ps_element_mesh,d2_coorg_recv_ps_element_mesh, &
             d1_color_send_ps_element_mesh, &
             d1_color_recv_ps_element_mesh

  integer :: d1_coorg_send_ps_abs, d2_coorg_send_ps_abs, &
             d1_coorg_recv_ps_abs, d2_coorg_recv_ps_abs

  integer :: d1_coorg_send_ps_free_surface, d2_coorg_send_ps_free_surface, &
             d1_coorg_recv_ps_free_surface, d2_coorg_recv_ps_free_surface

  integer :: d1_coorg_send_ps_vector_field, d2_coorg_send_ps_vector_field, &
             d1_coorg_recv_ps_vector_field, d2_coorg_recv_ps_vector_field

  ! checks if anything to do
  if (.not. (output_color_image .or. output_postscript_snapshot)) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Preparing image coloring'
    call flush_IMAIN()
  endif

  ! for color images
  if (output_color_image) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  allocating color image arrays'
      call flush_IMAIN()
    endif

    ! prepares dimension of image
    call prepare_color_image_init()

    ! allocate an array for image data
    allocate(image_color_data(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
    if (ier /= 0) stop 'error in an allocate statement 1'
    allocate(image_color_vp_display(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
    if (ier /= 0) stop 'error in an allocate statement 2'

    ! allocate an array for the grid point that corresponds to a given image data point
    allocate(iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
    if (ier /= 0) stop 'error in an allocate statement 3'
    allocate(copy_iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
    if (ier /= 0) stop 'error in an allocate statement 4'

    !remember which image are going to produce
    if (USE_SNAPSHOT_NUMBER_IN_FILENAME) then
      ! initializes counter
      isnapshot_number = 0
    endif

    ! creates pixels indexing
    call prepare_color_image_pixels()

    ! creating and filling array num_pixel_loc with the positions of each colored
    ! pixel owned by the local process (useful for parallel jobs)
    allocate(num_pixel_loc(nb_pixel_loc))

    ipixel = 0
    do i = 1, NX_IMAGE_color
       do j = 1, NZ_IMAGE_color
          if (iglob_image_color(i,j) /= -1) then
             ipixel = ipixel + 1
             num_pixel_loc(ipixel) = (j-1)*NX_IMAGE_color + i
          endif
       enddo
    enddo
    ! checks
    if (ipixel /= nb_pixel_loc) then
      print *,'Error: pixel count ',ipixel,nb_pixel_loc
      stop 'Error invalid pixel count for color image'
    endif
    ! filling array iglob_image_color, containing info on which process owns which pixels.
    iproc = 0
    k = 0
#ifdef USE_MPI
    allocate(nb_pixel_per_proc(0:NPROC-1))
    nb_pixel_per_proc(:) = 0
    call gather_all_singlei(nb_pixel_loc,nb_pixel_per_proc,NPROC)

    if (myrank == 0) then
      allocate(num_pixel_recv(maxval(nb_pixel_per_proc(:)),NPROC))
      allocate(data_pixel_recv(maxval(nb_pixel_per_proc(:))))
    endif

    allocate(data_pixel_send(nb_pixel_loc))
    if (NPROC > 1) then
      if (myrank == 0) then
        ! master collects
        do iproc = 1, NPROC-1
          call MPI_RECV(num_pixel_recv(1,iproc+1),nb_pixel_per_proc(iproc), MPI_INTEGER, &
                        iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
          do k = 1, nb_pixel_per_proc(iproc)
            j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
            i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color

            ! avoid edge effects
            if (i < 1) i = 1
            if (j < 1) j = 1

            if (i > NX_IMAGE_color) i = NX_IMAGE_color
            if (j > NZ_IMAGE_color) j = NZ_IMAGE_color

            iglob_image_color(i,j) = iproc

          enddo
        enddo

      else
        call MPI_SEND(num_pixel_loc(1),nb_pixel_loc,MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)
      endif
    endif
    call synchronize_all()
#endif

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  done locating all the pixels of color images'
      call flush_IMAIN()
    endif

    ! prepares image background
    call prepare_color_image_vp()

  endif ! output_color_image

  ! synchronizes all processes
  call synchronize_all()

  ! postscript output
  if (output_postscript_snapshot) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  allocating postscript image arrays'
      call flush_IMAIN()
    endif

    ! allocate arrays for postscript output
#ifdef USE_MPI
    if (modelvect) then
      d1_coorg_recv_ps_velocity_model=2
      call mpi_allreduce(nspec,d2_coorg_recv_ps_velocity_model,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
      d2_coorg_recv_ps_velocity_model=d2_coorg_recv_ps_velocity_model*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
         ((NGLLX-subsamp_postscript)/subsamp_postscript)*4
      d1_RGB_recv_ps_velocity_model=1
      call mpi_allreduce(nspec,d2_RGB_recv_ps_velocity_model,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
      d2_RGB_recv_ps_velocity_model=d2_RGB_recv_ps_velocity_model*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
         ((NGLLX-subsamp_postscript)/subsamp_postscript)*4
    else
      d1_coorg_recv_ps_velocity_model=1
      d2_coorg_recv_ps_velocity_model=1
      d1_RGB_recv_ps_velocity_model=1
      d2_RGB_recv_ps_velocity_model=1
    endif

    d1_coorg_send_ps_element_mesh=2
    if (ngnod == 4) then
      if (DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT == 1) then
        d2_coorg_send_ps_element_mesh=nspec*5
        if (DISPLAY_COLORS == 1) then
          d1_color_send_ps_element_mesh=2*nspec
        else
          d1_color_send_ps_element_mesh=1*nspec
        endif
      else
        d2_coorg_send_ps_element_mesh=nspec*6
        if (DISPLAY_COLORS == 1) then
          d1_color_send_ps_element_mesh=1*nspec
        endif
      endif
    else
      if (DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT == 1) then
        d2_coorg_send_ps_element_mesh=nspec*((pointsdisp-1)*3+max(0,pointsdisp-2)+1+1)
        if (DISPLAY_COLORS == 1) then
          d1_color_send_ps_element_mesh=2*nspec
        else
          d1_color_send_ps_element_mesh=1*nspec
        endif
      else
        d2_coorg_send_ps_element_mesh=nspec*((pointsdisp-1)*3+max(0,pointsdisp-2)+1)
        if (DISPLAY_COLORS == 1) then
          d1_color_send_ps_element_mesh=1*nspec
        endif
      endif
    endif

    call mpi_allreduce(d1_coorg_send_ps_element_mesh,d1_coorg_recv_ps_element_mesh,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
    call mpi_allreduce(d2_coorg_send_ps_element_mesh,d2_coorg_recv_ps_element_mesh,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
    call mpi_allreduce(d1_color_send_ps_element_mesh,d1_color_recv_ps_element_mesh,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

    d1_coorg_send_ps_abs=4
    d2_coorg_send_ps_abs=4*nelemabs
    call mpi_allreduce(d1_coorg_send_ps_abs,d1_coorg_recv_ps_abs,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
    call mpi_allreduce(d2_coorg_send_ps_abs,d2_coorg_recv_ps_abs,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

    d1_coorg_send_ps_free_surface=4
    d2_coorg_send_ps_free_surface=4*nelem_acoustic_surface
    call mpi_allreduce(d1_coorg_send_ps_free_surface,d1_coorg_recv_ps_free_surface,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
    call mpi_allreduce(d2_coorg_send_ps_free_surface,d2_coorg_recv_ps_free_surface,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

    d1_coorg_send_ps_vector_field=8
    if (interpol) then
      if (plot_lowerleft_corner_only) then
        d2_coorg_send_ps_vector_field=nspec*1*1
      else
        d2_coorg_send_ps_vector_field=nspec*pointsdisp*pointsdisp
      endif
    else
      d2_coorg_send_ps_vector_field=nglob
    endif
    call mpi_allreduce(d1_coorg_send_ps_vector_field,d1_coorg_recv_ps_vector_field,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
    call mpi_allreduce(d2_coorg_send_ps_vector_field,d2_coorg_recv_ps_vector_field,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

#else
    ! dummy values
    d1_coorg_recv_ps_velocity_model=1
    d2_coorg_recv_ps_velocity_model=1
    d1_RGB_recv_ps_velocity_model=1
    d2_RGB_recv_ps_velocity_model=1

    d1_coorg_send_ps_element_mesh=1
    d2_coorg_send_ps_element_mesh=1
    d1_coorg_recv_ps_element_mesh=1
    d2_coorg_recv_ps_element_mesh=1
    d1_color_send_ps_element_mesh=1
    d1_color_recv_ps_element_mesh=1

    d1_coorg_send_ps_abs=1
    d2_coorg_send_ps_abs=1
    d1_coorg_recv_ps_abs=1
    d2_coorg_recv_ps_abs=1
    d1_coorg_send_ps_free_surface=1
    d2_coorg_send_ps_free_surface=1
    d1_coorg_recv_ps_free_surface=1
    d2_coorg_recv_ps_free_surface=1

    d1_coorg_send_ps_vector_field=1
    d2_coorg_send_ps_vector_field=1
    d1_coorg_recv_ps_vector_field=1
    d2_coorg_recv_ps_vector_field=1
#endif

    d1_coorg_send_ps_velocity_model=2
    d2_coorg_send_ps_velocity_model=nspec*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
                                          ((NGLLX-subsamp_postscript)/subsamp_postscript)*4
    d1_RGB_send_ps_velocity_model=1
    d2_RGB_send_ps_velocity_model=nspec*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
                                        ((NGLLX-subsamp_postscript)/subsamp_postscript)

    allocate(coorg_send_ps_velocity_model(d1_coorg_send_ps_velocity_model,d2_coorg_send_ps_velocity_model))
    allocate(RGB_send_ps_velocity_model(d1_RGB_send_ps_velocity_model,d2_RGB_send_ps_velocity_model))

    allocate(coorg_recv_ps_velocity_model(d1_coorg_recv_ps_velocity_model,d2_coorg_recv_ps_velocity_model))
    allocate(RGB_recv_ps_velocity_model(d1_RGB_recv_ps_velocity_model,d2_RGB_recv_ps_velocity_model))

    allocate(coorg_send_ps_element_mesh(d1_coorg_send_ps_element_mesh,d2_coorg_send_ps_element_mesh))
    allocate(coorg_recv_ps_element_mesh(d1_coorg_recv_ps_element_mesh,d2_coorg_recv_ps_element_mesh))
    allocate(color_send_ps_element_mesh(d1_color_send_ps_element_mesh))
    allocate(color_recv_ps_element_mesh(d1_color_recv_ps_element_mesh))

    allocate(coorg_send_ps_abs(d1_coorg_send_ps_abs,d2_coorg_send_ps_abs))
    allocate(coorg_recv_ps_abs(d1_coorg_recv_ps_abs,d2_coorg_recv_ps_abs))

    allocate(coorg_send_ps_free_surface(d1_coorg_send_ps_free_surface,d2_coorg_send_ps_free_surface))
    allocate(coorg_recv_ps_free_surface(d1_coorg_recv_ps_free_surface,d2_coorg_recv_ps_free_surface))

    allocate(coorg_send_ps_vector_field(d1_coorg_send_ps_vector_field,d2_coorg_send_ps_vector_field))
    allocate(coorg_recv_ps_vector_field(d1_coorg_recv_ps_vector_field,d2_coorg_recv_ps_vector_field),stat=ier)
    if (ier /= 0) stop 'Error allocating postscript arrays'

    ! to dump the wave field
    this_is_the_first_time_we_dump = .true.

  endif ! postscript

  end subroutine prepare_timerun_image_coloring

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_adjoint()

! prepares adjoint runs

  use constants, only: IMAIN,USE_PORO_VISCOUS_DAMPING
  use specfem_par

  implicit none

  ! local parameters
  integer :: reclen,ier
  character(len=MAX_STRING_LEN) :: outputname,outputname2

  ! checks if anything to do
  if (.not. (SAVE_FORWARD .or. SIMULATION_TYPE == 3)) return


  ! prepares kernels
  if (SIMULATION_TYPE == 3) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Preparing adjoint simulation'
      call flush_IMAIN()
    endif

    call prepare_timerun_kernels()

    ! Absorbing boundaries
    if (STACEY_ABSORBING_CONDITIONS) then
      ! Reads last frame for forward wavefield reconstruction
      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  reading Stacey boundary arrays'
        call flush_IMAIN()
      endif

      ! reads in absorbing boundary data
      if (any_acoustic) call prepare_absorb_read_acoustic()
      if (any_elastic) call prepare_absorb_read_elastic()
      if (any_poroelastic) call prepare_absorb_read_poroelastic()
    endif
  endif

  ! Files where viscous damping are saved during forward wavefield calculation
  if (POROELASTIC_SIMULATION .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3)) then
    if (USE_PORO_VISCOUS_DAMPING) then
      ! user output
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'Preparing save forward/adjoint simulation for poroelasticity'
        write(IMAIN,*) '  using viscous damping arrays for poroelastic domain'
        call flush_IMAIN()
      endif

      ! allocate only when this slice contains poroelastic elements
      if (any_poroelastic) then
        allocate(b_viscodampx(NGLLX,NGLLZ,nspec), &
                 b_viscodampz(NGLLX,NGLLZ,nspec),stat=ier)
        if (ier /= 0) stop 'Error allocating b_viscodamp arrays'
        b_viscodampx(:,:,:) = 0._CUSTOM_REAL
        b_viscodampz(:,:,:) = 0._CUSTOM_REAL

        ! file i/o
        ! array size
        reclen = CUSTOM_REAL * NGLLX * NGLLZ * nspec

        write(outputname,'(a,i6.6,a)') 'viscodampingx',myrank,'.bin'
        write(outputname2,'(a,i6.6,a)') 'viscodampingz',myrank,'.bin'

        if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
          open(unit=23,file='OUTPUT_FILES/'//outputname,status='unknown', &
                form='unformatted',access='direct',recl=reclen,iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening file OUTPUT_FILES/viscodampingx**.bin')

          open(unit=24,file='OUTPUT_FILES/'//outputname2,status='unknown', &
                form='unformatted',access='direct',recl=reclen,iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening file OUTPUT_FILES/viscodampingz**.bin')

        else if (SIMULATION_TYPE == 3) then
          open(unit=23,file='OUTPUT_FILES/'//outputname,status='old', &
                action='read',form='unformatted',access='direct',recl=reclen,iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening file OUTPUT_FILES/viscodampingx**.bin')

          open(unit=24,file='OUTPUT_FILES/'//outputname2,status='old', &
                action='read',form='unformatted',access='direct',recl=reclen,iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening file OUTPUT_FILES/viscodampingz**.bin')
        endif
      endif
    endif
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_adjoint


!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_kernels()


#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IMAIN,APPROXIMATE_HESS_KL
  use specfem_par

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  initializing adjoint sensitivity kernels'
    call flush_IMAIN()
  endif

  ! Allocates sensitivity kernel arrays
  ! elastic domains
  if (any_elastic) then

    if (save_ASCII_kernels) then
      ! ascii format
      if (count(ispec_is_anisotropic(:) .eqv. .true.) >= 1) then ! anisotropic
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_cijkl_kernel.dat'
        open(unit = 97, file='OUTPUT_FILES/'//outputname,status='unknown',iostat=ier)
        if (ier /= 0) stop 'Error writing kernel file to disk'
      else
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_mu_kernel.dat'
        open(unit = 97, file = 'OUTPUT_FILES/'//outputname,status='unknown',iostat=ier)
        if (ier /= 0) stop 'Error writing kernel file to disk'

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_alpha_beta_kernel.dat'
        open(unit = 98, file = 'OUTPUT_FILES/'//outputname,status='unknown',iostat=ier)
        if (ier /= 0) stop 'Error writing kernel file to disk'
      endif
    else
      ! binary format
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kernel.bin'
      open(unit = 204, file = 'OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error writing kernel file to disk'

      if (count(ispec_is_anisotropic(:) .eqv. .true.) >= 1) then ! anisotropic
         write(outputname,'(a,i6.6,a)')'proc',myrank,'_c11_kernel.bin'
         open(unit = 205,file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
         if (ier /= 0) stop 'Error writing kernel file to disk'

         write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c13_kernel.bin'
         open(unit = 206,file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
         if (ier /= 0) stop 'Error writing kernel file to disk'

         write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c15_kernel.bin'
         open(unit = 207,file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
         if (ier /= 0) stop 'Error writing kernel file to disk'

         write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c33_kernel.bin'
         open(unit = 208,file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
         if (ier /= 0) stop 'Error writing kernel file to disk'

         write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c35_kernel.bin'
         open(unit = 209,file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
         if (ier /= 0) stop 'Error writing kernel file to disk'

         write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c55_kernel.bin'
         open(unit = 210,file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
         if (ier /= 0) stop 'Error writing kernel file to disk'
      else
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_kappa_kernel.bin'
        open(unit = 205, file ='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
        if (ier /= 0) stop 'Error writing kernel file to disk'

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mu_kernel.bin'
        open(unit = 206, file ='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
        if (ier /= 0) stop 'Error writing kernel file to disk'

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_kernel.bin'
        open(unit = 207, file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
        if (ier /= 0) stop 'Error writing kernel file to disk'

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_alpha_kernel.bin'
        open(unit = 208, file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
        if (ier /= 0) stop 'Error writing kernel file to disk'

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_beta_kernel.bin'
        open(unit = 209, file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
        if (ier /= 0) stop 'Error writing kernel file to disk'

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_bulk_c_kernel.bin'
        open(unit = 210,file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
        if (ier /= 0) stop 'Error writing kernel file to disk'

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_bulk_beta_kernel.bin'
        open(unit = 211,file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
        if (ier /= 0) stop 'Error writing kernel file to disk'

        if (APPROXIMATE_HESS_KL) then
          write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian1_kernel.bin'
          open(unit =214,file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
          if (ier /= 0) stop 'Error writing kernel file to disk'

          write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian2_kernel.bin'
          open(unit=215,file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
          if (ier /= 0) stop 'Error writing kernel file to disk'
        endif
      endif

    endif

    rho_kl(:,:,:) = 0._CUSTOM_REAL
    mu_kl(:,:,:) = 0._CUSTOM_REAL
    kappa_kl(:,:,:) = 0._CUSTOM_REAL

    rhop_kl(:,:,:) = 0._CUSTOM_REAL
    beta_kl(:,:,:) = 0._CUSTOM_REAL
    alpha_kl(:,:,:) = 0._CUSTOM_REAL
    bulk_c_kl(:,:,:) = 0._CUSTOM_REAL
    bulk_beta_kl(:,:,:) = 0._CUSTOM_REAL

    if (APPROXIMATE_HESS_KL) then
      rhorho_el_Hessian_final2(:,:,:) = 0._CUSTOM_REAL
      rhorho_el_Hessian_final1(:,:,:) = 0._CUSTOM_REAL
    endif
  endif

  ! poro-elastic domains
  if (any_poroelastic) then

    if (.not. save_ASCII_kernels) stop 'poroelastic simulations must use save_ASCII_kernels'

    ! Primary kernels
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mu_B_C_kernel.dat'
    open(unit = 144, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) stop 'Error writing kernel file to disk'
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_M_rho_rhof_kernel.dat'
    open(unit = 155, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) stop 'Error writing kernel file to disk'
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_m_eta_kernel.dat'
    open(unit = 16, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) stop 'Error writing kernel file to disk'
    ! Wavespeed kernels
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_cpI_cpII_cs_kernel.dat'
    open(unit = 20, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) stop 'Error writing kernel file to disk'
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhobb_rhofbb_ratio_kernel.dat'
    open(unit = 21, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) stop 'Error writing kernel file to disk'
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_phib_eta_kernel.dat'
    open(unit = 22, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) stop 'Error writing kernel file to disk'
    ! Density normalized kernels
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mub_Bb_Cb_kernel.dat'
    open(unit = 17, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) stop 'Error writing kernel file to disk'
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Mb_rhob_rhofb_kernel.dat'
    open(unit = 18, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) stop 'Error writing kernel file to disk'
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mb_etab_kernel.dat'
    open(unit = 19, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) stop 'Error writing kernel file to disk'

    rhot_kl(:,:,:) = 0._CUSTOM_REAL
    rhof_kl(:,:,:) = 0._CUSTOM_REAL
    eta_kl(:,:,:) = 0._CUSTOM_REAL
    sm_kl(:,:,:) = 0._CUSTOM_REAL
    mufr_kl(:,:,:) = 0._CUSTOM_REAL
    B_kl(:,:,:) = 0._CUSTOM_REAL
    C_kl(:,:,:) = 0._CUSTOM_REAL
    M_kl(:,:,:) = 0._CUSTOM_REAL

    rhob_kl(:,:,:) = 0._CUSTOM_REAL
    rhofb_kl(:,:,:) = 0._CUSTOM_REAL
    phi_kl(:,:,:) = 0._CUSTOM_REAL
    mufrb_kl(:,:,:) = 0._CUSTOM_REAL

    rhobb_kl(:,:,:) = 0._CUSTOM_REAL
    rhofbb_kl(:,:,:) = 0._CUSTOM_REAL
    phib_kl(:,:,:) = 0._CUSTOM_REAL
    cs_kl(:,:,:) = 0._CUSTOM_REAL
    cpI_kl(:,:,:) = 0._CUSTOM_REAL
    cpII_kl(:,:,:) = 0._CUSTOM_REAL
    ratio_kl(:,:,:) = 0._CUSTOM_REAL
  endif

  ! acoustic domains
  if (any_acoustic) then

    if (save_ASCII_kernels) then
      ! ascii format
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_kernel.dat'
      open(unit = 95, file = 'OUTPUT_FILES/'//outputname,status ='unknown',iostat=ier)
      if (ier /= 0) stop 'Error writing kernel file to disk'

      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_c_kernel.dat'
      open(unit = 96, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ier)
      if (ier /= 0) stop 'Error writing kernel file to disk'

    else
      ! binary format
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_acoustic_kernel.bin'
      open(unit = 200, file = 'OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error writing kernel file to disk'

      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_kappa_acoustic_kernel.bin'
      open(unit = 201, file = 'OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error writing kernel file to disk'

      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_acoustic_kernel.bin'
      open(unit = 202, file = 'OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error writing kernel file to disk'

      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c_acoustic_kernel.bin'
      open(unit = 203, file = 'OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) stop 'Error writing kernel file to disk'

      if (APPROXIMATE_HESS_KL) then
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian1_acoustic_kernel.bin'
        open(unit=212,file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
        if (ier /= 0) stop 'Error writing kernel file to disk'

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian2_acoustic_kernel.bin'
        open(unit=213,file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
        if (ier /= 0) stop 'Error writing kernel file to disk'
      endif
    endif

    rho_ac_kl(:,:,:) = 0._CUSTOM_REAL
    kappa_ac_kl(:,:,:) = 0._CUSTOM_REAL

    rhop_ac_kl(:,:,:) = 0._CUSTOM_REAL
    alpha_ac_kl(:,:,:) = 0._CUSTOM_REAL

    if (APPROXIMATE_HESS_KL) then
      rhorho_ac_Hessian_final2(:,:,:) = 0._CUSTOM_REAL
      rhorho_ac_Hessian_final1(:,:,:) = 0._CUSTOM_REAL
    endif
  endif


  end subroutine prepare_timerun_kernels

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_initialfield()

! reads initial fields from external file if needed

  use constants, only: IMAIN
  use specfem_par

  implicit none

  ! local parameters
  double precision :: cploc,csloc

  ! if we are looking a plane wave beyond critical angle we use other method
  over_critical_angle = .false.

  ! checks if anything to do
  if (.not. initialfield) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Preparing initial field for plane wave source'
    call flush_IMAIN()
  endif

  ! safety checks
  if (.not. any_elastic) &
    stop 'Sorry, initial field (plane wave source) only implemented for elastic simulations so far...'
  if (any_acoustic .or. any_poroelastic) &
    stop 'Initial field currently implemented for purely elastic simulation only'

  ! Calculation of the initial field for a plane wave
  if (any_elastic) then
    call prepare_initial_field(cploc,csloc)

    if (over_critical_angle) then

      allocate(left_bound(nelemabs*NGLLX))
      allocate(right_bound(nelemabs*NGLLX))
      allocate(bot_bound(nelemabs*NGLLZ))

      call prepare_initial_field_paco()

      allocate(v0x_left(count_left,NSTEP))
      allocate(v0z_left(count_left,NSTEP))
      allocate(t0x_left(count_left,NSTEP))
      allocate(t0z_left(count_left,NSTEP))

      allocate(v0x_right(count_right,NSTEP))
      allocate(v0z_right(count_right,NSTEP))
      allocate(t0x_right(count_right,NSTEP))
      allocate(t0z_right(count_right,NSTEP))

      allocate(v0x_bot(count_bottom,NSTEP))
      allocate(v0z_bot(count_bottom,NSTEP))
      allocate(t0x_bot(count_bottom,NSTEP))
      allocate(t0z_bot(count_bottom,NSTEP))

      ! call Paco's routine to compute in frequency and convert to time by Fourier transform
      call paco_beyond_critical(anglesource(1), &
                                f0_source(1),QKappa_attenuation(1),source_type(1),left_bound(1:count_left), &
                                right_bound(1:count_right),bot_bound(1:count_bottom), &
                                count_left,count_right,count_bottom,x_source(1),cploc,csloc)

      deallocate(left_bound)
      deallocate(right_bound)
      deallocate(bot_bound)

      if (myrank == 0) then
        write(IMAIN,*)  '***********'
        write(IMAIN,*)  'done calculating the initial wave field'
        write(IMAIN,*)  '***********'
        call flush_IMAIN()
      endif

    endif ! beyond critical angle

    if (myrank == 0) then
      write(IMAIN,*) 'Max norm of initial elastic displacement = ', &
                      maxval(sqrt(displ_elastic(1,:)**2 + displ_elastic(2,:)**2))
      call flush_IMAIN()
    endif

  endif

  end subroutine prepare_timerun_initialfield


!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_noise()

! for noise simulations

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: NGLLX,NGLLZ,NDIM,IMAIN,NOISE_MOVIE_OUTPUT,TWO_THIRDS

  use specfem_par, only: myrank,AXISYM,NSTEP,nglob,nspec,ibool,coord, &
                         rhoext,vpext,vsext,density,poroelastcoef,kmato,assign_external_model
  use specfem_par_noise

  implicit none

  integer :: i,j,iglob,ispec,ier

  ! checks if anything to do
  if (NOISE_TOMOGRAPHY <= 0) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Preparing noise simulation'
    call flush_IMAIN()
  endif

  ! allocates arrays for noise tomography
  allocate(time_function_noise(NSTEP), &
           source_array_noise(NDIM,NGLLX,NGLLZ,NSTEP), &
           mask_noise(nglob), &
           surface_movie_y_or_z_noise(nglob),stat=ier)
  if (ier /= 0) stop 'Error allocating noise arrays'

  ! initializes
  source_array_noise(:,:,:,:) = 0._CUSTOM_REAL
  surface_movie_y_or_z_noise(:) = 0._CUSTOM_REAL

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  reading noise parameters'
    call flush_IMAIN()
  endif

  !read in parameters for noise tomography
  call read_parameters_noise()

  if (NOISE_TOMOGRAPHY == 1) then
    call compute_source_array_noise()

    ! write out coordinates of mesh
    open(unit=504,file='OUTPUT_FILES/mesh_spec',status='unknown',action='write')
      do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            write(504,'(1pe11.3,1pe11.3,2i3,i7)') coord(1,iglob), coord(2,iglob), i, j, ispec
         enddo
        enddo
      enddo
    close(504)

    open(unit=504,file='OUTPUT_FILES/mesh_glob',status='unknown',action='write')
      do iglob = 1, nglob
        write(504,'(1pe11.3,1pe11.3,i7)') coord(1,iglob), coord(2,iglob), iglob
      enddo
    close(504)

    ! write out spatial distribution of noise sources
    call create_mask_noise()
    open(unit=504,file='OUTPUT_FILES/mask_noise',status='unknown',action='write')
      do iglob = 1, nglob
            write(504,'(1pe11.3,1pe11.3,1pe11.3)') coord(1,iglob), coord(2,iglob), mask_noise(iglob)
      enddo
    close(504)

    ! write out velocity model
    if (assign_external_model) then
      open(unit=504,file='OUTPUT_FILES/model_rho_vp_vs',status='unknown',action='write')
        do ispec = 1, nspec
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob = ibool(i,j,ispec)
              write(504,'(1pe11.3,1pe11.3,1pe11.3,1pe11.3,1pe11.3)') &
                coord(1,iglob), coord(2,iglob), rhoext(i,j,ispec), vpext(i,j,ispec), vsext(i,j,ispec)
            enddo
          enddo
        enddo
      close(504)
    else
      open(unit=504,file='OUTPUT_FILES/model_rho_kappa_mu',status='unknown',action='write')
        do ispec = 1, nspec
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob = ibool(i,j,ispec)
              if (AXISYM) then ! CHECK kappa
                write(504,'(1pe11.3,1pe11.3,1pe11.3,1pe11.3,1pe11.3)') &
                  coord(1,iglob), coord(2,iglob), density(1,kmato(ispec)), &
                  poroelastcoef(1,1,kmato(ispec)) + TWO_THIRDS * poroelastcoef(2,1,kmato(ispec)), &
                  poroelastcoef(2,1,kmato(ispec))
              else
                write(504,'(1pe11.3,1pe11.3,1pe11.3,1pe11.3,1pe11.3)') &
                  coord(1,iglob), coord(2,iglob), density(1,kmato(ispec)), &
                  poroelastcoef(1,1,kmato(ispec)) + poroelastcoef(2,1,kmato(ispec)), &
                  poroelastcoef(2,1,kmato(ispec))
              endif
            enddo
          enddo
        enddo
      close(504)
    endif

  else if (NOISE_TOMOGRAPHY == 2) then
    call create_mask_noise()

  else if (NOISE_TOMOGRAPHY == 3) then
    ! noise movie
    if (NOISE_MOVIE_OUTPUT) then
      call create_mask_noise()
      ! prepare array that will hold wavefield snapshots
      noise_output_ncol = 5
      allocate(noise_output_array(noise_output_ncol,nglob), &
               noise_output_rhokl(nglob))
    endif

  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_noise


!
!-------------------------------------------------------------------------------------
!


  subroutine prepare_timerun_attenuation()

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IMAIN,TWO,PI,FOUR_THIRDS,TWO_THIRDS
  use specfem_par

  implicit none

  ! local parameters
  integer :: i,j,ispec,n,ier
  ! for shifting of velocities if needed in the case of viscoelasticity
  double precision :: vp,vs,rhol,mul,lambdal,kappal
  double precision :: qkappal,qmul
  ! attenuation factors
  real(kind=CUSTOM_REAL) :: Mu_nu1_sent,Mu_nu2_sent
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: tau_epsilon_nu1_sent,tau_epsilon_nu2_sent
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: inv_tau_sigma_nu1_sent,inv_tau_sigma_nu2_sent, &
                                                       phi_nu1_sent,phi_nu2_sent

  ! attenuation
  if (ATTENUATION_VISCOELASTIC) then
    nspec_ATT = nspec
  else
    nspec_ATT = 1
  endif

  ! allocate memory variables for attenuation
  allocate(e1(NGLLX,NGLLZ,nspec_ATT,N_SLS), &
           e11(NGLLX,NGLLZ,nspec_ATT,N_SLS), &
           e13(NGLLX,NGLLZ,nspec_ATT,N_SLS),stat=ier)
  if (ier /= 0) stop 'Error allocating attenuation arrays'
  e1(:,:,:,:) = 0._CUSTOM_REAL
  e11(:,:,:,:) = 0._CUSTOM_REAL
  e13(:,:,:,:) = 0._CUSTOM_REAL

  if (time_stepping_scheme == 2) then
    allocate(e1_LDDRK(NGLLX,NGLLZ,nspec_ATT,N_SLS))
    allocate(e11_LDDRK(NGLLX,NGLLZ,nspec_ATT,N_SLS))
    allocate(e13_LDDRK(NGLLX,NGLLZ,nspec_ATT,N_SLS))
  else
    allocate(e1_LDDRK(1,1,1,1))
    allocate(e11_LDDRK(1,1,1,1))
    allocate(e13_LDDRK(1,1,1,1))
  endif
  e1_LDDRK(:,:,:,:) = 0._CUSTOM_REAL
  e11_LDDRK(:,:,:,:) = 0._CUSTOM_REAL
  e13_LDDRK(:,:,:,:) = 0._CUSTOM_REAL

  if (time_stepping_scheme == 3) then
    allocate(e1_initial_rk(NGLLX,NGLLZ,nspec_ATT,N_SLS))
    allocate(e11_initial_rk(NGLLX,NGLLZ,nspec_ATT,N_SLS))
    allocate(e13_initial_rk(NGLLX,NGLLZ,nspec_ATT,N_SLS))
    allocate(e1_force_rk(NGLLX,NGLLZ,nspec_ATT,N_SLS,stage_time_scheme))
    allocate(e11_force_rk(NGLLX,NGLLZ,nspec_ATT,N_SLS,stage_time_scheme))
    allocate(e13_force_rk(NGLLX,NGLLZ,nspec_ATT,N_SLS,stage_time_scheme))
  else
    allocate(e1_initial_rk(1,1,1,1))
    allocate(e11_initial_rk(1,1,1,1))
    allocate(e13_initial_rk(1,1,1,1))
    allocate(e1_force_rk(1,1,1,1,1))
    allocate(e11_force_rk(1,1,1,1,1))
    allocate(e13_force_rk(1,1,1,1,1))
  endif
  e1_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
  e11_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
  e13_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
  e1_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL
  e11_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL
  e13_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL

  ! attenuation arrays
  if (.not. assign_external_model) then
    allocate(already_shifted_velocity(numat),stat=ier)
    if (ier /= 0) stop 'Error allocating attenuation Qkappa,Qmu,.. arrays'
    already_shifted_velocity(:) = .false.
  endif

  allocate(inv_tau_sigma_nu1(NGLLX,NGLLZ,nspec_ATT,N_SLS), &
           inv_tau_sigma_nu2(NGLLX,NGLLZ,nspec_ATT,N_SLS), &
           phi_nu1(NGLLX,NGLLZ,nspec_ATT,N_SLS), &
           phi_nu2(NGLLX,NGLLZ,nspec_ATT,N_SLS), &
           Mu_nu1(NGLLX,NGLLZ,nspec_ATT), &
           Mu_nu2(NGLLX,NGLLZ,nspec_ATT),stat=ier)
  if (ier /= 0) stop 'Error allocating attenuation arrays'

  ! temporary arrays for function argument
  allocate(tau_epsilon_nu1_sent(N_SLS), &
           tau_epsilon_nu2_sent(N_SLS), &
           inv_tau_sigma_nu1_sent(N_SLS), &
           inv_tau_sigma_nu2_sent(N_SLS), &
           phi_nu1_sent(N_SLS), &
           phi_nu2_sent(N_SLS),stat=ier)
  if (ier /= 0) stop 'Error allocating attenuation coefficient arrays'

  ! initialize to dummy values
  ! convention to indicate that Q = 9999 in that element i.e. that there is no viscoelasticity in that element
  inv_tau_sigma_nu1(:,:,:,:) = -1._CUSTOM_REAL
  phi_nu1(:,:,:,:) = -1._CUSTOM_REAL

  inv_tau_sigma_nu2(:,:,:,:) = -1._CUSTOM_REAL
  phi_nu2(:,:,:,:) = -1._CUSTOM_REAL

  Mu_nu1(:,:,:) = -1._CUSTOM_REAL
  Mu_nu2(:,:,:) = -1._CUSTOM_REAL

  ! if source is not a Dirac or Heavyside then f0_attenuation is f0 of the first source
  if (.not. (time_function_type(1) == 4 .or. time_function_type(1) == 5)) then
    f0_attenuation = f0_source(1)
  endif

  ! setup attenuation
  if (ATTENUATION_VISCOELASTIC) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Preparing attenuation'
      write(IMAIN,*) '  using external model for Qkappa & Qmu: ',assign_external_model
      write(IMAIN,*) '  reading velocity at f0               : ',READ_VELOCITIES_AT_f0
      write(IMAIN,*) '  using f0 attenuation = ',f0_attenuation,'Hz'
      call flush_IMAIN()
    endif

    ! define the attenuation quality factors.
!! DK DK if needed in the future, here the quality factor could be different for each point
    do ispec = 1,nspec

      ! attenuation is not implemented in acoustic (i.e. fluid) media for now, only in viscoelastic (i.e. solid) media
      if (ispec_is_acoustic(ispec)) cycle

      ! get values for internal meshes
      if (.not. assign_external_model) then
        qkappal = QKappa_attenuation(kmato(ispec))
        qmul = Qmu_attenuation(kmato(ispec))

        ! check that attenuation values entered by the user make sense
        if ((qkappal <= 9998.999d0 .and. qmul > 9998.999d0) .or. &
            (qkappal > 9998.999d0 .and. qmul <= 9998.999d0)) &
           stop 'need to have Qkappa and Qmu both above or both below 9999 for a given material; &
                &trick: use 9998 if you want to turn off one'

        ! if no attenuation in that elastic element
        if (qkappal > 9998.999d0) cycle

        ! determines attenuation factors
        call attenuation_model(qkappal,qmul,f0_attenuation,N_SLS, &
                               tau_epsilon_nu1_sent,inv_tau_sigma_nu1_sent,phi_nu1_sent,Mu_nu1_sent, &
                               tau_epsilon_nu2_sent,inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu2_sent)
      endif

      do j = 1,NGLLZ
        do i = 1,NGLLX

          ! get values for external meshes
          if (assign_external_model) then
            qkappal = QKappa_attenuationext(i,j,ispec)
            qmul = Qmu_attenuationext(i,j,ispec)

            ! check that attenuation values entered by the user make sense
            if ((qkappal <= 9998.999d0 .and. qmul > 9998.999d0) .or. &
                (qkappal > 9998.999d0 .and. qmul <= 9998.999d0)) &
               stop 'need to have Qkappa and Qmu both above or both below 9999 for a given material; &
                    &trick: use 9998 if you want to turn off one'

            ! if no attenuation in that elastic element
            if (qkappal > 9998.999d0) cycle

            ! determines attenuation factors
            call attenuation_model(qkappal,qmul,f0_attenuation,N_SLS, &
                                 tau_epsilon_nu1_sent,inv_tau_sigma_nu1_sent,phi_nu1_sent,Mu_nu1_sent, &
                                 tau_epsilon_nu2_sent,inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu2_sent)
          endif

          ! stores attenuation values
          inv_tau_sigma_nu1(i,j,ispec,:) = inv_tau_sigma_nu1_sent(:)
          phi_nu1(i,j,ispec,:) = phi_nu1_sent(:)

          inv_tau_sigma_nu2(i,j,ispec,:) = inv_tau_sigma_nu2_sent(:)
          phi_nu2(i,j,ispec,:) = phi_nu2_sent(:)

          Mu_nu1(i,j,ispec) = Mu_nu1_sent
          Mu_nu2(i,j,ispec) = Mu_nu2_sent

          ! shifts velocities
          if (READ_VELOCITIES_AT_f0) then
            ! safety check
            if (ispec_is_anisotropic(ispec) .or. ispec_is_poroelastic(ispec) .or. ispec_is_gravitoacoustic(ispec)) &
              stop 'READ_VELOCITIES_AT_f0 only implemented for non anisotropic, &
                    &non poroelastic, non gravitoacoustic materials for now'

            if (assign_external_model) then
              ! external mesh model
              rhol = dble(rhoext(i,j,ispec))
              vp = dble(vpext(i,j,ispec))
              vs = dble(vsext(i,j,ispec))

              ! shifts vp & vs (according to f0 and attenuation band)
              call shift_velocities_from_f0(vp,vs,rhol, &
                                    f0_attenuation,N_SLS, &
                                    tau_epsilon_nu1_sent,tau_epsilon_nu2_sent,inv_tau_sigma_nu1_sent,inv_tau_sigma_nu2_sent)

              ! stores shifted values
              vpext(i,j,ispec) = vp
              vsext(i,j,ispec) = vs
            else
              ! internal mesh
              n = kmato(ispec)
              if (.not. already_shifted_velocity(n)) then
                rhol = density(1,n)
                lambdal = poroelastcoef(1,1,n)
                mul = poroelastcoef(2,1,n)

                vp = sqrt((lambdal + TWO * mul) / rhol)
                vs = sqrt(mul/rhol)

                ! shifts vp & vs
                call shift_velocities_from_f0(vp,vs,rhol, &
                                      f0_attenuation,N_SLS, &
                                      tau_epsilon_nu1_sent,tau_epsilon_nu2_sent,inv_tau_sigma_nu1_sent,inv_tau_sigma_nu2_sent)

                ! stores shifted mu,lambda
                mul = rhol * vs*vs
                lambdal = rhol * vp*vp - TWO * mul

                poroelastcoef(1,1,n) = lambdal
                poroelastcoef(2,1,n) = mul
                poroelastcoef(3,1,n) = lambdal + TWO * mul

                already_shifted_velocity(n) = .true.
              endif
            endif
          endif
        enddo
      enddo
    enddo
  endif ! ATTENUATION_VISCOELASTIC

  ! sets new material properties
  ! note: velocities might have been shifted by attenuation
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Preparing material arrays'
    call flush_IMAIN()
  endif

  ! allocates material arrays
  allocate(kappastore(NGLLX,NGLLZ,nspec), &
           mustore(NGLLX,NGLLZ,nspec), &
           rhostore(NGLLX,NGLLZ,nspec), &
           rho_vp(NGLLX,NGLLZ,nspec), &
           rho_vs(NGLLX,NGLLZ,nspec),stat=ier)
  if (ier /= 0) stop 'Error allocating material arrays'

  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        if (assign_external_model) then
          ! external model
          rhol = rhoext(i,j,ispec)
          vp = vpext(i,j,ispec)
          vs = vsext(i,j,ispec)
          ! determins mu & kappa
          mul = rhol * vs * vs
          if (AXISYM) then ! CHECK kappa
            kappal = rhol * vp * vp - FOUR_THIRDS * mul
          else
            kappal = rhol * vp * vp - mul
          endif
        else
          ! internal mesh
          rhol = density(1,kmato(ispec))
          lambdal = poroelastcoef(1,1,kmato(ispec))  ! lambdal_unrelaxed_elastic
          mul = poroelastcoef(2,1,kmato(ispec)) ! mul_unrelaxed_elastic
          if (AXISYM) then ! CHECK kappa
            kappal = lambdal + TWO_THIRDS * mul
            vp = sqrt((kappal + FOUR_THIRDS * mul)/rhol)
          else
            kappal = lambdal + mul
            vp = sqrt((kappal + mul)/rhol)
          endif
        endif
        ! stores moduli
        rhostore(i,j,ispec) = rhol
        mustore(i,j,ispec) = mul
        kappastore(i,j,ispec) = kappal

        ! stores density times vp & vs

        vs = sqrt(mul/rhol)

        rho_vp(i,j,ispec) = rhol * vp
        rho_vs(i,j,ispec) = rhol * vs
      enddo
    enddo
  enddo

  ! allocate memory variables for viscous attenuation (poroelastic media)
  if (ATTENUATION_PORO_FLUID_PART) then
    allocate(rx_viscous(NGLLX,NGLLZ,nspec))
    allocate(rz_viscous(NGLLX,NGLLZ,nspec))
    allocate(viscox(NGLLX,NGLLZ,nspec))
    allocate(viscoz(NGLLX,NGLLZ,nspec))
    ! initialize memory variables for attenuation
    rx_viscous(:,:,:) = 0.d0
    rz_viscous(:,:,:) = 0.d0
    viscox(:,:,:) = 0.d0
    viscoz(:,:,:) = 0.d0

    if (time_stepping_scheme == 2) then
      allocate(rx_viscous_LDDRK(NGLLX,NGLLZ,nspec))
      allocate(rz_viscous_LDDRK(NGLLX,NGLLZ,nspec))
      rx_viscous_LDDRK(:,:,:) = 0.d0
      rz_viscous_LDDRK(:,:,:) = 0.d0
    endif

    if (time_stepping_scheme == 3) then
      allocate(rx_viscous_initial_rk(NGLLX,NGLLZ,nspec))
      allocate(rz_viscous_initial_rk(NGLLX,NGLLZ,nspec))
      allocate(rx_viscous_force_RK(NGLLX,NGLLZ,nspec,stage_time_scheme))
      allocate(rz_viscous_force_RK(NGLLX,NGLLZ,nspec,stage_time_scheme))
      rx_viscous_initial_rk(:,:,:) = 0.d0
      rz_viscous_initial_rk(:,:,:) = 0.d0
      rx_viscous_force_RK(:,:,:,:) = 0.d0
      rz_viscous_force_RK(:,:,:,:) = 0.d0
    endif

    ! Precompute Runge Kutta coefficients if viscous attenuation
    ! viscous attenuation is implemented following the memory variable formulation of
    ! J. M. Carcione Wave fields in real media: wave propagation in anisotropic,
    ! anelastic and porous media, Elsevier, p. 304-305, 2007
    theta_e = (sqrt(Q0**2+1.d0) +1.d0)/(2.d0*pi*freq0*Q0)
    theta_s = (sqrt(Q0**2+1.d0) -1.d0)/(2.d0*pi*freq0*Q0)

    thetainv = - 1.d0 / theta_s
    alphaval = 1.d0 + deltat*thetainv + deltat**2*thetainv**2 / 2.d0 &
                    + deltat**3*thetainv**3 / 6.d0 + deltat**4*thetainv**4 / 24.d0
    betaval = deltat / 2.d0 + deltat**2*thetainv / 3.d0 + deltat**3*thetainv**2 / 8.d0 + deltat**4*thetainv**3 / 24.d0
    gammaval = deltat / 2.d0 + deltat**2*thetainv / 6.d0 + deltat**3*thetainv**2 / 24.d0
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_attenuation

