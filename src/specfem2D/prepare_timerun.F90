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

  subroutine prepare_timerun()

  use constants, only: USE_ENFORCE_FIELDS,IOUT_ENERGY,IMAIN,OUTPUT_FILES
  use specfem_par
  use specfem_par_movie
  use specfem_par_noise, only: NOISE_TOMOGRAPHY

  implicit none

  if (setup_with_binary_database == 2 ) then
   call setup_mesh_external_models()
   call read_sources_receivers()
   if (SIMULATION_TYPE == 3) call setup_adjoint_sources()
   call read_binary_database_part1()
  endif

  ! Test compatibility with axisymmetric formulation
  if (AXISYM) call check_compatibility_axisym()

  ! prepares constant factors for time scheme and seismograms
  call prepare_timerun_constants()

  ! wavefield array initialization
  call prepare_wavefields()

  ! PML preparation
  call prepare_PML()

  if (setup_with_binary_database /= 2) then

    ! attenuation
    !! DK DK moved preparation of attenuation to before preparation of mass matrix
    !! DK DK because now that we have support for viscoacoustic fluids, we need to use
    !! DK DK the unrelaxed Kappa modulus in the fluid mass matrix when ATTENUATION_VISCOACOUSTIC is on
    !! DK DK and thus we need to prepare attenuation before preparing the mass matrix
    call prepare_timerun_attenuation()

    ! prepares mass matrices
    call prepare_timerun_mass_matrix()

    ! check which GLL will be forced or not
    if (USE_ENFORCE_FIELDS) call build_forced()

    ! postscript images for grids and snapshots
    call prepare_timerun_postscripts()

    ! for adjoint kernel runs
    call prepare_timerun_adjoint()

    ! reads initial fields from external file if needed
    call prepare_timerun_initialfield

    ! compute the source time function and stores it in a text file
    call prepare_source_time_function()

    ! prepares noise simulations
    if (NOISE_TOMOGRAPHY /= 0) call prepare_timerun_noise()

    if (setup_with_binary_database == 1 ) call save_binary_database()

  else

    call read_binary_database_part2()

  endif

  ! compute material property arrays
  call prepare_material_arrays()

  ! jpeg images
  call prepare_timerun_image_coloring()

  ! init specific to NO_BACKWARD_RECONSTRUCTION option
  call prepare_timerun_no_backward_reconstruction()

  ! prepares GPU arrays
  if (GPU_MODE) call prepare_GPU()

  !-------------------------------------------------------------

  ! creates a Gnuplot script to display the energy curve in log scale
  if (OUTPUT_ENERGY .and. myrank == 0) then
    close(IOUT_ENERGY)
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'plot_energy.gnu',status='unknown',action='write')
    write(IOUT_ENERGY,*) 'set term wxt'
    write(IOUT_ENERGY,*) '#set term postscript landscape color solid "Helvetica" 22'
    write(IOUT_ENERGY,*) '#set output "energy.ps"'
    write(IOUT_ENERGY,*) '# set xrange [0:60]'
    write(IOUT_ENERGY,*) 'set logscale y'
    write(IOUT_ENERGY,*) 'set xlabel "Time (s)"'
    write(IOUT_ENERGY,*) 'set ylabel "Energy (J)"'
    write(IOUT_ENERGY,*) 'set loadpath "'//trim(OUTPUT_FILES)//'"'
    write(IOUT_ENERGY,'(a)') &
      'plot "energy.dat" us 1:4 t "Total Energy" w l lc 1, "energy.dat" us 1:3 t "Potential Energy" w l lc 2'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    close(IOUT_ENERGY)
  endif

  ! open the file in which we will store the energy curve
  if (OUTPUT_ENERGY .and. myrank == 0) open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'energy.dat',status='unknown',action='write')

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
    if (UNDO_ATTENUATION_AND_OR_PML) then
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
  if (ier /= 0) call stop_the_code('Error allocating seismogram arrays')

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

  ! local variable
#ifdef USE_MPI
  integer :: n_sls_loc
#endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Preparing mass matrices'
    call flush_IMAIN()
  endif

  ! builds the global mass matrix
  call invert_mass_matrix_init()

#ifdef USE_MPI
  n_sls_loc = 0
  if (ATTENUATION_VISCOACOUSTIC) n_sls_loc = N_SLS
  ! assembling the mass matrix of shared nodes on MPI partition interfaces
  call assemble_MPI_scalar(rmass_inverse_acoustic,nglob_acoustic, &
                           rmass_inverse_e1,n_sls_loc, &
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
  if (ier /= 0) call stop_the_code('Error allocating shape arrays for display')

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

  ! wavefield dump
  integer :: d1_dump_send, d2_dump_send, d1_dump_recv, d2_dump_recv

  ! checks if anything to do
  if (.not. (output_color_image .or. output_postscript_snapshot .or. output_wavefield_dumps)) return

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
    if (ier /= 0) call stop_the_code('error in an allocate statement 1')
    allocate(image_color_vp_display(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
    if (ier /= 0) call stop_the_code('error in an allocate statement 2')

    ! allocate an array for the grid point that corresponds to a given image data point
    allocate(iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
    if (ier /= 0) call stop_the_code('error in an allocate statement 3')
    allocate(copy_iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
    if (ier /= 0) call stop_the_code('error in an allocate statement 4')

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
      call stop_the_code('Error invalid pixel count for color image')
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
          call recv_i(num_pixel_recv(1,iproc+1), nb_pixel_per_proc(iproc), iproc, 42)

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
        call send_i(num_pixel_loc(1), nb_pixel_loc, 0, 42)
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
      d1_coorg_recv_ps_velocity_model = 2
      call max_all_all_i(nspec,d2_coorg_recv_ps_velocity_model)
      d2_coorg_recv_ps_velocity_model = d2_coorg_recv_ps_velocity_model*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
         ((NGLLX-subsamp_postscript)/subsamp_postscript)*4
      d1_RGB_recv_ps_velocity_model = 1

      call max_all_all_i(nspec,d2_RGB_recv_ps_velocity_model)
      d2_RGB_recv_ps_velocity_model = d2_RGB_recv_ps_velocity_model*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
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

    call max_all_all_i(d1_coorg_send_ps_element_mesh,d1_coorg_recv_ps_element_mesh)
    call max_all_all_i(d2_coorg_send_ps_element_mesh,d2_coorg_recv_ps_element_mesh)
    call max_all_all_i(d1_color_send_ps_element_mesh,d1_color_recv_ps_element_mesh)

    d1_coorg_send_ps_abs=5
    d2_coorg_send_ps_abs=4*nelemabs
    call max_all_all_i(d1_coorg_send_ps_abs,d1_coorg_recv_ps_abs)
    call max_all_all_i(d2_coorg_send_ps_abs,d2_coorg_recv_ps_abs)

    d1_coorg_send_ps_free_surface=4
    d2_coorg_send_ps_free_surface=4*nelem_acoustic_surface
    call max_all_all_i(d1_coorg_send_ps_free_surface,d1_coorg_recv_ps_free_surface)
    call max_all_all_i(d2_coorg_send_ps_free_surface,d2_coorg_recv_ps_free_surface)

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
    call max_all_all_i(d1_coorg_send_ps_vector_field,d1_coorg_recv_ps_vector_field)
    call max_all_all_i(d2_coorg_send_ps_vector_field,d2_coorg_recv_ps_vector_field)

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
    if (ier /= 0) call stop_the_code('Error allocating postscript arrays')

    ! to dump the wave field
    this_is_the_first_time_we_dump = .true.

  endif ! postscript

  if (output_wavefield_dumps) then
    ! allocate arrays for wavefield dump
    d1_dump_recv = 2
    d1_dump_send = 2

    d2_dump_send = nspec*NGLLX*NGLLZ

#ifdef USE_MPI
    call mpi_allreduce(d2_dump_send,d2_dump_recv,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
#else
    d2_dump_recv = d2_dump_send
#endif

    allocate(dump_send(d1_dump_send, d2_dump_send))
    allocate(dump_recv(d1_dump_recv, d2_dump_recv))

    allocate(dump_duplicate_send(d2_dump_send))
    allocate(dump_duplicate_recv(d2_dump_recv))

    allocate(dump_recv_counts(0:NPROC-1))

    this_is_the_first_time_we_dump = .true.

  endif ! wavefield dump

  end subroutine prepare_timerun_image_coloring

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_adjoint()

! prepares adjoint runs

  use constants, only: IMAIN,USE_PORO_VISCOUS_DAMPING,OUTPUT_FILES
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
        if (ier /= 0) call stop_the_code('Error allocating b_viscodamp arrays')
        b_viscodampx(:,:,:) = 0._CUSTOM_REAL
        b_viscodampz(:,:,:) = 0._CUSTOM_REAL

        ! file i/o
        ! array size
        reclen = CUSTOM_REAL * NGLLX * NGLLZ * nspec

        write(outputname,'(a,i6.6,a)') 'viscodampingx',myrank,'.bin'
        write(outputname2,'(a,i6.6,a)') 'viscodampingz',myrank,'.bin'

        if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
          open(unit=23,file=trim(OUTPUT_FILES)//outputname,status='unknown', &
                form='unformatted',access='direct',recl=reclen,iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'viscodampingx**.bin')

          open(unit=24,file=trim(OUTPUT_FILES)//outputname2,status='unknown', &
                form='unformatted',access='direct',recl=reclen,iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'viscodampingz**.bin')

        else if (SIMULATION_TYPE == 3) then
          open(unit=23,file=trim(OUTPUT_FILES)//outputname,status='old', &
                action='read',form='unformatted',access='direct',recl=reclen,iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'viscodampingx**.bin')

          open(unit=24,file=trim(OUTPUT_FILES)//outputname2,status='old', &
                action='read',form='unformatted',access='direct',recl=reclen,iostat=ier)
          if (ier /= 0) call exit_MPI(myrank,'Error opening file '//trim(OUTPUT_FILES)//'viscodampingz**.bin')
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

  use constants, only: IMAIN,APPROXIMATE_HESS_KL,OUTPUT_FILES
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
        open(unit = 97, file=trim(OUTPUT_FILES)//outputname,status='unknown',iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      else
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_mu_kernel.dat'
        open(unit = 97, file = trim(OUTPUT_FILES)//outputname,status='unknown',iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_alpha_beta_kernel.dat'
        open(unit = 98, file = trim(OUTPUT_FILES)//outputname,status='unknown',iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      endif
    else
      ! binary format
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kernel.bin'
      open(unit = 204, file = trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

      if (count(ispec_is_anisotropic(:) .eqv. .true.) >= 1) then ! anisotropic
         write(outputname,'(a,i6.6,a)')'proc',myrank,'_c11_kernel.bin'
         open(unit = 205,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
         if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

         write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c13_kernel.bin'
         open(unit = 206,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
         if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

         write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c15_kernel.bin'
         open(unit = 207,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
         if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

         write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c33_kernel.bin'
         open(unit = 208,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
         if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

         write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c35_kernel.bin'
         open(unit = 209,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
         if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

         write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c55_kernel.bin'
         open(unit = 210,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
         if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      else
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_kappa_kernel.bin'
        open(unit = 205, file =trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mu_kernel.bin'
        open(unit = 206, file =trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_kernel.bin'
        open(unit = 207, file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_alpha_kernel.bin'
        open(unit = 208, file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_beta_kernel.bin'
        open(unit = 209, file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_bulk_c_kernel.bin'
        open(unit = 210,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_bulk_beta_kernel.bin'
        open(unit = 211,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

        if (APPROXIMATE_HESS_KL) then
          write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian1_kernel.bin'
          open(unit =214,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
          if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

          write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian2_kernel.bin'
          open(unit=215,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
          if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
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

    if (.not. save_ASCII_kernels) call stop_the_code('poroelastic simulations must use save_ASCII_kernels')

    ! Primary kernels
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mu_B_C_kernel.dat'
    open(unit = 144, file = trim(OUTPUT_FILES)//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_M_rho_rhof_kernel.dat'
    open(unit = 155, file = trim(OUTPUT_FILES)//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_m_eta_kernel.dat'
    open(unit = 16, file = trim(OUTPUT_FILES)//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    ! Wavespeed kernels
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_cpI_cpII_cs_kernel.dat'
    open(unit = 20, file = trim(OUTPUT_FILES)//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhobb_rhofbb_ratio_kernel.dat'
    open(unit = 21, file = trim(OUTPUT_FILES)//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_phib_eta_kernel.dat'
    open(unit = 22, file = trim(OUTPUT_FILES)//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    ! Density normalized kernels
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mub_Bb_Cb_kernel.dat'
    open(unit = 17, file = trim(OUTPUT_FILES)//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Mb_rhob_rhofb_kernel.dat'
    open(unit = 18, file = trim(OUTPUT_FILES)//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mb_etab_kernel.dat'
    open(unit = 19, file = trim(OUTPUT_FILES)//outputname,status = 'unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

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
      open(unit = 95, file = trim(OUTPUT_FILES)//outputname,status ='unknown',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_c_kernel.dat'
      open(unit = 96, file = trim(OUTPUT_FILES)//outputname,status = 'unknown',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

    else
      ! binary format
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_acoustic_kernel.bin'
      open(unit = 200, file = trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_kappa_acoustic_kernel.bin'
      open(unit = 201, file = trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_acoustic_kernel.bin'
      open(unit = 202, file = trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c_acoustic_kernel.bin'
      open(unit = 203, file = trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

      if (APPROXIMATE_HESS_KL) then
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian1_acoustic_kernel.bin'
        open(unit=212,file = trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing kernel file to disk')

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian2_acoustic_kernel.bin'
        open(unit=213,file = trim(OUTPUT_FILES)//outputname,status='unknown',action='write',form='unformatted',iostat=ier)
        if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
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
    call stop_the_code('Sorry, initial field (plane wave source) only implemented for elastic simulations so far...')
  if (any_acoustic .or. any_poroelastic) &
    call stop_the_code('Initial field currently implemented for purely elastic simulation only')

  ! Calculation of the initial field for a plane wave
  if (any_elastic) then

    ! calculates initial plane wave coefficients
    call prepare_initial_field(cploc,csloc)

    ! special case for Rayleigh waves, SV waves above critical angle
    if (over_critical_angle) then

      ! allocates boundaries
      allocate(left_bound(nelemabs*NGLLX))
      allocate(right_bound(nelemabs*NGLLX))
      allocate(bot_bound(nelemabs*NGLLZ))

      ! sets up boundary points
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
      call paco_beyond_critical(anglesource(1),f0_source(1), &
                                QKappa_attenuation(1),source_type(1), &
                                left_bound(1:count_left),right_bound(1:count_right),bot_bound(1:count_bottom), &
                                count_left,count_right,count_bottom, &
                                x_source(1),cploc,csloc)

      ! frees memory
      deallocate(left_bound)
      deallocate(right_bound)
      deallocate(bot_bound)

      ! user output
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*)  '***********'
        write(IMAIN,*)  'done calculating the initial wave field'
        write(IMAIN,*)  '***********'
        write(IMAIN,*)
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

  use constants, only: NGLLX,NGLLZ,NDIM,IMAIN,NOISE_MOVIE_OUTPUT,TWO_THIRDS,OUTPUT_FILES

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
  if (ier /= 0) call stop_the_code('Error allocating noise arrays')

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
    open(unit=504,file=trim(OUTPUT_FILES)//'mesh_spec',status='unknown',action='write')
      do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            write(504,'(1pe11.3,1pe11.3,2i3,i7)') coord(1,iglob), coord(2,iglob), i, j, ispec
         enddo
        enddo
      enddo
    close(504)

    open(unit=504,file=trim(OUTPUT_FILES)//'mesh_glob',status='unknown',action='write')
      do iglob = 1, nglob
        write(504,'(1pe11.3,1pe11.3,i7)') coord(1,iglob), coord(2,iglob), iglob
      enddo
    close(504)

    ! write out spatial distribution of noise sources
    call create_mask_noise()
    open(unit=504,file=trim(OUTPUT_FILES)//'mask_noise',status='unknown',action='write')
      do iglob = 1, nglob
            write(504,'(1pe11.3,1pe11.3,1pe11.3)') coord(1,iglob), coord(2,iglob), mask_noise(iglob)
      enddo
    close(504)

    ! write out velocity model
    if (assign_external_model) then
      open(unit=504,file=trim(OUTPUT_FILES)//'model_rho_vp_vs',status='unknown',action='write')
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
      open(unit=504,file=trim(OUTPUT_FILES)//'model_rho_kappa_mu',status='unknown',action='write')
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

  use constants, only: IMAIN,TWO,PI,FOUR_THIRDS,TWO_THIRDS,USE_A_STRONG_FORMULATION_FOR_E1
  use specfem_par

  implicit none

  ! local parameters
  integer :: i,j,ispec,n,ier
  ! for shifting of velocities if needed in the case of viscoelasticity
  double precision :: vp,vs,rhol,mul,lambdal
  double precision :: qkappal,qmul
  ! attenuation factors
  real(kind=CUSTOM_REAL) :: Mu_nu1_sent,Mu_nu2_sent
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: tau_epsilon_nu1_sent,tau_epsilon_nu2_sent
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: inv_tau_sigma_nu1_sent,inv_tau_sigma_nu2_sent, &
                                                       phi_nu1_sent,phi_nu2_sent
  real(kind=CUSTOM_REAL), dimension(N_SLS) ::  phinu,tauinvnu,temp,coef

  ! attenuation
  if (ATTENUATION_VISCOELASTIC) then
    nspec_ATT_el = nspec
  else
    nspec_ATT_el = 1
  endif

  ! allocate memory variables for attenuation
  allocate(e1(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
           e11(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
           e13(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
           dux_dxl_old(NGLLX,NGLLZ,nspec_ATT_el), &
           duz_dzl_old(NGLLX,NGLLZ,nspec_ATT_el), &
           dux_dzl_plus_duz_dxl_old(NGLLX,NGLLZ,nspec_ATT_el), &
           A_newmark_nu1(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
           B_newmark_nu1(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
           A_newmark_nu2(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), &
           B_newmark_nu2(N_SLS,NGLLX,NGLLZ,nspec_ATT_el), stat=ier)

  ! attenuation
  if (ATTENUATION_VISCOACOUSTIC) then
    nglob_ATT = nglob
    nspec_ATT_ac = nspec
  else
    nglob_ATT = 1
    nspec_ATT_ac = 1
  endif

  allocate(e1_acous_sf(N_SLS,NGLLX,NGLLZ,nspec_ATT_ac), &
           sum_forces_old(NGLLX,NGLLZ,nspec_ATT_ac), stat=ier)

  if (ATTENUATION_VISCOACOUSTIC .and. .not. USE_A_STRONG_FORMULATION_FOR_E1) then

    allocate(e1_acous(nglob_acoustic,N_SLS), &
             dot_e1(nglob_acoustic,N_SLS), &
             A_newmark_e1_sf(1,1,1,1), &
             B_newmark_e1_sf(1,1,1,1),stat=ier)
    if (time_stepping_scheme == 1) then
      allocate(dot_e1_old(nglob_acoustic,N_SLS), &
               A_newmark_e1(nglob_acoustic,N_SLS), &
               B_newmark_e1(nglob_acoustic,N_SLS),stat=ier)
    else
      allocate(dot_e1_old(1,N_SLS), &
               A_newmark_e1(1,N_SLS), &
               B_newmark_e1(1,N_SLS),stat=ier)
    endif

    if (time_stepping_scheme == 2) then
      allocate(e1_acous_temp(nglob_acoustic,N_SLS),stat=ier)
    else
      allocate(e1_acous_temp(1,N_SLS),stat=ier)
    endif

  else if (ATTENUATION_VISCOACOUSTIC .and. USE_A_STRONG_FORMULATION_FOR_E1) then

    allocate(e1_acous(1,N_SLS), &
             e1_acous_temp(1,N_SLS), &
             dot_e1(1,N_SLS), &
             dot_e1_old(1,N_SLS), &
             A_newmark_e1(1,N_SLS), &
             B_newmark_e1(1,N_SLS),stat=ier)
    allocate(A_newmark_e1_sf(N_SLS,NGLLX,NGLLZ,nspec), &
             B_newmark_e1_sf(N_SLS,NGLLX,NGLLZ,nspec),stat=ier)

  else ! no ATTENUATION_VISCOACOUSTIC

    allocate(e1_acous(1,N_SLS), &
             e1_acous_temp(1,N_SLS), &
             dot_e1(1,N_SLS), &
             dot_e1_old(1,N_SLS), &
             A_newmark_e1(1,N_SLS), &
             B_newmark_e1(1,N_SLS), &
             A_newmark_e1_sf(1,1,1,1), &
             B_newmark_e1_sf(1,1,1,1),stat=ier)
  endif

  if (ier /= 0) call stop_the_code('Error allocating attenuation arrays')
  e1(:,:,:,:) = 0._CUSTOM_REAL
  e11(:,:,:,:) = 0._CUSTOM_REAL
  e13(:,:,:,:) = 0._CUSTOM_REAL
  dux_dxl_old(:,:,:) = 0._CUSTOM_REAL
  duz_dzl_old(:,:,:) = 0._CUSTOM_REAL
  dux_dzl_plus_duz_dxl_old(:,:,:) = 0._CUSTOM_REAL

  e1_acous(:,:) = 0._CUSTOM_REAL
  e1_acous_sf(:,:,:,:) = 0._CUSTOM_REAL

  dot_e1_old = 0._CUSTOM_REAL
  dot_e1     = 0._CUSTOM_REAL
  sum_forces_old = 0._CUSTOM_REAL

  if (time_stepping_scheme == 2) then
    if (ATTENUATION_VISCOELASTIC) then
      allocate(e1_LDDRK(NGLLX,NGLLZ,nspec_ATT_el,N_SLS))
      allocate(e11_LDDRK(NGLLX,NGLLZ,nspec_ATT_el,N_SLS))
      allocate(e13_LDDRK(NGLLX,NGLLZ,nspec_ATT_el,N_SLS))
    else
      allocate(e1_LDDRK(1,1,1,1))
      allocate(e11_LDDRK(1,1,1,1))
      allocate(e13_LDDRK(1,1,1,1))
    endif

    if (ATTENUATION_VISCOACOUSTIC) then
        allocate(e1_LDDRK_acous(nglob_att,N_SLS))
    else
        allocate(e1_LDDRK_acous(1,1))
    endif
  else
    allocate(e1_LDDRK(1,1,1,1))
    allocate(e11_LDDRK(1,1,1,1))
    allocate(e13_LDDRK(1,1,1,1))

    allocate(e1_LDDRK_acous(1,1))
  endif
  e1_LDDRK(:,:,:,:) = 0._CUSTOM_REAL
  e11_LDDRK(:,:,:,:) = 0._CUSTOM_REAL
  e13_LDDRK(:,:,:,:) = 0._CUSTOM_REAL

  e1_LDDRK_acous(:,:) = 0._CUSTOM_REAL

  if (time_stepping_scheme == 3) then
    allocate(e1_initial_rk(NGLLX,NGLLZ,nspec_ATT_el,N_SLS))
    allocate(e11_initial_rk(NGLLX,NGLLZ,nspec_ATT_el,N_SLS))
    allocate(e13_initial_rk(NGLLX,NGLLZ,nspec_ATT_el,N_SLS))
    allocate(e1_force_rk(NGLLX,NGLLZ,nspec_ATT_el,N_SLS,stage_time_scheme))
    allocate(e11_force_rk(NGLLX,NGLLZ,nspec_ATT_el,N_SLS,stage_time_scheme))
    allocate(e13_force_rk(NGLLX,NGLLZ,nspec_ATT_el,N_SLS,stage_time_scheme))

    if (ATTENUATION_VISCOACOUSTIC) then
            allocate(e1_initial_rk_acous(nglob_att,N_SLS))
            allocate(e1_force_rk_acous(nglob_att,N_SLS,stage_time_scheme))
    else
            allocate(e1_initial_rk_acous(1,1))
            allocate(e1_force_rk_acous(1,1,1))
    endif
  else
    allocate(e1_initial_rk(1,1,1,1))
    allocate(e11_initial_rk(1,1,1,1))
    allocate(e13_initial_rk(1,1,1,1))
    allocate(e1_force_rk(1,1,1,1,1))
    allocate(e11_force_rk(1,1,1,1,1))
    allocate(e13_force_rk(1,1,1,1,1))

    allocate(e1_initial_rk_acous(1,1))
    allocate(e1_force_rk_acous(1,1,1))
  endif
  e1_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
  e11_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
  e13_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
  e1_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL
  e11_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL
  e13_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL

  e1_initial_rk_acous(:,:) = 0._CUSTOM_REAL
  e1_force_rk_acous(:,:,:) = 0._CUSTOM_REAL

  ! attenuation arrays
  if (.not. assign_external_model) then
    allocate(already_shifted_velocity(numat),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating attenuation Qkappa,Qmu,.. arrays')
    already_shifted_velocity(:) = .false.
  endif

  allocate(inv_tau_sigma_nu1(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac),N_SLS), &
           inv_tau_sigma_nu2(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac),N_SLS), &
           phi_nu1(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac),N_SLS), &
           phi_nu2(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac),N_SLS), &
           Mu_nu1(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac)), &
           Mu_nu2(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac)), &
! ZX ZX needed for further optimization with nspec_ATT_el replaced with nspec_PML
           tau_epsilon_nu1(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac),N_SLS), &
           tau_epsilon_nu2(NGLLX,NGLLZ,max(nspec_ATT_el,nspec_ATT_ac),N_SLS),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating attenuation arrays')

  ! temporary arrays for function argument
  allocate(tau_epsilon_nu1_sent(N_SLS), &
           tau_epsilon_nu2_sent(N_SLS), &
           inv_tau_sigma_nu1_sent(N_SLS), &
           inv_tau_sigma_nu2_sent(N_SLS), &
           phi_nu1_sent(N_SLS), &
           phi_nu2_sent(N_SLS),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating attenuation coefficient arrays')

  ! initialize to dummy values
  ! convention to indicate that Q = 9999 in that element i.e. that there is no viscoelasticity in that element
  inv_tau_sigma_nu1(:,:,:,:) = -1._CUSTOM_REAL
  inv_tau_sigma_nu2(:,:,:,:) = -1._CUSTOM_REAL

  tau_epsilon_nu1(:,:,:,:) = -1._CUSTOM_REAL
  tau_epsilon_nu2(:,:,:,:) = -1._CUSTOM_REAL

  phi_nu1(:,:,:,:) = 0._CUSTOM_REAL
  phi_nu2(:,:,:,:) = 0._CUSTOM_REAL

  ! do not change this, in the case of a viscoacoustic medium the mass matrix is multiplied by this,
  ! and thus the factor needs to be equal to +1 when QKappa = 9999 i.e. when viscoacousticity is turned off in parts of the medium
  Mu_nu1(:,:,:) = +1._CUSTOM_REAL
  Mu_nu2(:,:,:) = +1._CUSTOM_REAL

  ! if source is not a Dirac or Heavyside then ATTENUATION_f0_REFERENCE is f0 of the first source
  if (.not. (time_function_type(1) == 4 .or. time_function_type(1) == 5)) then
    ATTENUATION_f0_REFERENCE = f0_source(1)
  endif

  ! setup attenuation
  if (ATTENUATION_VISCOELASTIC .or. ATTENUATION_VISCOACOUSTIC) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Preparing attenuation in viscoelastic or viscoacoustic parts of the model:'
      write(IMAIN,*) '  using external model for Qkappa and Qmu: ',assign_external_model
      write(IMAIN,*) '  reading velocity at f0                 : ',READ_VELOCITIES_AT_f0
      write(IMAIN,*) '  using an attenuation reference frequency of ',ATTENUATION_f0_REFERENCE,'Hz'
      call flush_IMAIN()
    endif

    ! define the attenuation quality factors.
    do ispec = 1,nspec

      ! get values for internal meshes
      if (.not. assign_external_model) then
        qkappal = QKappa_attenuation(kmato(ispec))
        qmul = Qmu_attenuation(kmato(ispec))
        ! if no attenuation in that elastic element
        if (qkappal > 9998.999d0 .and. qmul > 9998.999d0) cycle

        ! determines attenuation factors
        call attenuation_model(qkappal,qmul,ATTENUATION_f0_REFERENCE,N_SLS, &
                               tau_epsilon_nu1_sent,inv_tau_sigma_nu1_sent,phi_nu1_sent,Mu_nu1_sent, &
                               tau_epsilon_nu2_sent,inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu2_sent)
      endif

      do j = 1,NGLLZ
        do i = 1,NGLLX

          ! get values for external meshes
          if (assign_external_model) then

            qkappal = QKappa_attenuationext(i,j,ispec)
            qmul = Qmu_attenuationext(i,j,ispec)

            ! if no attenuation in that elastic element
            if (qkappal > 9998.999d0 .and. qmul > 9998.999d0) cycle

            ! determines attenuation factors
            call attenuation_model(qkappal,qmul,ATTENUATION_f0_REFERENCE,N_SLS, &
                                   tau_epsilon_nu1_sent,inv_tau_sigma_nu1_sent,phi_nu1_sent,Mu_nu1_sent, &
                                   tau_epsilon_nu2_sent,inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu2_sent)
          endif

          ! stores attenuation values
          inv_tau_sigma_nu1(i,j,ispec,:) = inv_tau_sigma_nu1_sent(:)
          tau_epsilon_nu1(i,j,ispec,:) = tau_epsilon_nu1_sent(:)
          phi_nu1(i,j,ispec,:) = phi_nu1_sent(:)

          inv_tau_sigma_nu2(i,j,ispec,:) = inv_tau_sigma_nu2_sent(:)
          tau_epsilon_nu2(i,j,ispec,:) = tau_epsilon_nu2_sent(:)
          phi_nu2(i,j,ispec,:) = phi_nu2_sent(:)

          Mu_nu1(i,j,ispec) = Mu_nu1_sent
          Mu_nu2(i,j,ispec) = Mu_nu2_sent

          if (ATTENUATION_VISCOACOUSTIC .and. USE_A_STRONG_FORMULATION_FOR_E1 .and. time_stepping_scheme == 1 ) then
            phinu(:)    = phi_nu1(i,j,ispec,:)
            tauinvnu(:) = inv_tau_sigma_nu1(i,j,ispec,:)
            temp(:)      = exp(- 0.5d0 * tauinvnu(:) * deltat)
            coef(:)     = (1.d0 - temp(:)) / tauinvnu(:)
            A_newmark_e1_sf(:,i,j,ispec) = temp(:)
            B_newmark_e1_sf(:,i,j,ispec) = phinu(:) * coef(:)
          endif

          if (ATTENUATION_VISCOELASTIC .and. time_stepping_scheme == 1 ) then
            phinu(:)    = phi_nu1(i,j,ispec,:)
            tauinvnu(:) = inv_tau_sigma_nu1(i,j,ispec,:)
            temp(:)      = exp(- 0.5d0 * tauinvnu(:) * deltat)
            coef(:)     = (1.d0 - temp(:)) / tauinvnu(:)
            A_newmark_nu1(:,i,j,ispec) = temp(:)
            B_newmark_nu1(:,i,j,ispec) = phinu(:) * coef(:)

            phinu(:)    = phi_nu2(i,j,ispec,:)
            tauinvnu(:) = inv_tau_sigma_nu2(i,j,ispec,:)
            temp(:)      = exp(- 0.5d0 * tauinvnu(:) * deltat)
            coef(:)     = (1.d0 - temp(:)) / tauinvnu(:)
            A_newmark_nu2(:,i,j,ispec) = temp(:)
            B_newmark_nu2(:,i,j,ispec) = phinu(:) * coef(:)
          endif

          ! shifts velocities
          if (READ_VELOCITIES_AT_f0) then

            ! safety check
            if (ispec_is_anisotropic(ispec) .or. ispec_is_poroelastic(ispec)) &
              call stop_the_code('READ_VELOCITIES_AT_f0 only implemented for non anisotropic, non poroelastic materials for now')

            if (ispec_is_acoustic(ispec)) then
              do n = 1,100
                print *,'WARNING: READ_VELOCITIES_AT_f0 in viscoacoustic elements may imply having to rebuild the mass matrix &
                   &with the shifted velocities, since the fluid mass matrix contains Kappa; not implemented yet, BEWARE!!'
              enddo
            endif

            if (assign_external_model) then
              ! external mesh model
              rhol = dble(rhoext(i,j,ispec))
              vp = dble(vpext(i,j,ispec))
              vs = dble(vsext(i,j,ispec))

              ! shifts vp and vs (according to f0 and attenuation band)
              call shift_velocities_from_f0(vp,vs,rhol, &
                                    ATTENUATION_f0_REFERENCE,N_SLS, &
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

                ! shifts vp and vs
                call shift_velocities_from_f0(vp,vs,rhol, &
                                      ATTENUATION_f0_REFERENCE,N_SLS, &
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

    if (PML_BOUNDARY_CONDITIONS) call prepare_timerun_attenuation_with_PML()

  endif ! of if (ATTENUATION_VISCOELASTIC .or. ATTENUATION_VISCOACOUSTIC)

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

    ! precompute Runge Kutta coefficients if viscous attenuation
    ! viscous attenuation is implemented following the memory variable formulation of
    ! J. M. Carcione Wave fields in real media: wave propagation in anisotropic,
    ! anelastic and porous media, Elsevier, p. 304-305, 2007
    theta_e = (sqrt(Q0_poroelastic**2+1.d0) +1.d0)/(2.d0*pi*freq0_poroelastic*Q0_poroelastic)
    theta_s = (sqrt(Q0_poroelastic**2+1.d0) -1.d0)/(2.d0*pi*freq0_poroelastic*Q0_poroelastic)

    thetainv = - 1.d0 / theta_s
    alphaval = 1.d0 + deltat*thetainv + deltat**2*thetainv**2 / 2.d0 &
                    + deltat**3*thetainv**3 / 6.d0 + deltat**4*thetainv**4 / 24.d0
    betaval = deltat / 2.d0 + deltat**2*thetainv / 3.d0 + deltat**3*thetainv**2 / 8.d0 + deltat**4*thetainv**3 / 24.d0
    gammaval = deltat / 2.d0 + deltat**2*thetainv / 6.d0 + deltat**3*thetainv**2 / 24.d0
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_attenuation

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_material_arrays()

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IMAIN,FOUR_THIRDS,TWO_THIRDS
  use specfem_par

  implicit none

  ! local parameters
  integer :: i,j,ispec,ier
  ! for shifting of velocities if needed in the case of viscoelasticity
  double precision :: vp,vs,rhol,mul,lambdal,kappal

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
  if (ier /= 0) call stop_the_code('Error allocating material arrays')

  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        if (assign_external_model) then
          ! external model
          rhol = rhoext(i,j,ispec)
          vp = vpext(i,j,ispec)
          vs = vsext(i,j,ispec)
          ! determins mu and kappa
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

        ! stores density times vp and vs
        vs = sqrt(mul/rhol)

        rho_vp(i,j,ispec) = rhol * vp
        rho_vs(i,j,ispec) = rhol * vs
      enddo
    enddo
  enddo

  end subroutine prepare_material_arrays

!
!-------------------------------------------------------------------------------------
!
  subroutine prepare_timerun_attenuation_with_PML()

  use constants
  use specfem_par

  double precision, dimension(2*N_SLS) :: tauinvnu
  double precision, dimension(2*N_SLS + 1) :: tauinvplusone
  double precision :: qkappal,qmul
  double precision :: d_x, d_z, K_x, K_z, alpha_x, alpha_z, beta_x, beta_z
  double precision :: const_for_separation_two

  integer :: i,j,ispec,i_sls,ispec_PML

  const_for_separation_two = min_distance_between_CPML_parameter * 2.d0
  do ispec = 1,nspec
    ispec_PML = spec_to_PML(ispec)
    if (ispec_is_PML(ispec)) then
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if (.not. assign_external_model) then
            qkappal = QKappa_attenuation(kmato(ispec))
            qmul = Qmu_attenuation(kmato(ispec))
          else
            qkappal = QKappa_attenuationext(i,j,ispec)
            qmul = Qmu_attenuationext(i,j,ispec)
          endif
          if (qkappal > 9998.999d0 .and. qmul > 9998.999d0) cycle
          K_x = K_x_store(i,j,ispec_PML)
          d_x = d_x_store(i,j,ispec_PML)
          alpha_x = alpha_x_store(i,j,ispec_PML)

          K_z = K_z_store(i,j,ispec_PML)
          d_z = d_z_store(i,j,ispec_PML)
          alpha_z = alpha_z_store(i,j,ispec_PML)

          if (ATTENUATION_VISCOELASTIC) then
            if (region_CPML(ispec) == CPML_XZ) then
              if (ispec_is_elastic(ispec)) then
                do i_sls = 1,N_SLS
                  tauinvnu(i_sls) = inv_tau_sigma_nu1(i,j,ispec,i_sls)
                  tauinvnu(i_sls+N_SLS) = inv_tau_sigma_nu2(i,j,ispec,i_sls)
                enddo

                call tauinvnu_arange_from_lt_to_gt(2*N_SLS,tauinvnu)

                do i_sls = 1,2*N_SLS
                  if (abs(alpha_x - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                    alpha_x = tauinvnu(i_sls) + const_for_separation_two
                  endif
                  if (abs(alpha_x - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                    call stop_the_code('error in separation of alpha_x, tauinvnu')
                  endif
                enddo
                do i_sls = 1,2*N_SLS
                  tauinvplusone(i_sls) = tauinvnu(i_sls)
                enddo
                tauinvplusone(2*N_SLS + 1) = alpha_x

                call tauinvnu_arange_from_lt_to_gt(2 * N_SLS + 1,tauinvplusone)

                do i_sls = 1,2*N_SLS + 1
                  if (abs(alpha_z - tauinvplusone(i_sls)) < min_distance_between_CPML_parameter) then
                    alpha_z = tauinvplusone(i_sls) + const_for_separation_two
                  endif
                  if (abs(alpha_z - tauinvplusone(i_sls)) < min_distance_between_CPML_parameter) then
                    call stop_the_code('error in separation of alpha_z, alpha_x,tauinvnu')
                  endif
                enddo

                beta_z = alpha_z + d_z / K_z

                do i_sls = 1,2*N_SLS + 1
                  if (abs(beta_z - tauinvplusone(i_sls)) < min_distance_between_CPML_parameter) then
                    beta_z = tauinvplusone(i_sls) + const_for_separation_two
                  endif
                  if (abs(beta_z - tauinvplusone(i_sls)) < min_distance_between_CPML_parameter) then
                    call stop_the_code('error in separation of beta_z, alpha_x,tauinvnu')
                  endif
                enddo

                do i_sls = 1,2*N_SLS
                  tauinvplusone(i_sls) = tauinvnu(i_sls)
                enddo
                tauinvplusone(2*N_SLS + 1) = alpha_z

                call tauinvnu_arange_from_lt_to_gt(2 * N_SLS + 1,tauinvplusone)
                beta_x = alpha_x + d_x / K_x
                do i_sls = 1,2*N_SLS + 1
                  if (abs(beta_x - tauinvplusone(i_sls)) < min_distance_between_CPML_parameter) then
                    beta_x = tauinvplusone(i_sls) + const_for_separation_two
                  endif
                  if (abs(beta_x - tauinvplusone(i_sls)) < min_distance_between_CPML_parameter) then
                    call stop_the_code('error in separation of beta_x, alpha_z,tauinvnu')
                  endif
                enddo

                d_x = (beta_x - alpha_x) * K_x
                d_z = (beta_z - alpha_z) * K_z

                d_x_store(i,j,ispec_PML) = d_x
                alpha_x_store(i,j,ispec_PML) = alpha_x
                d_z_store(i,j,ispec_PML) = d_z
                alpha_z_store(i,j,ispec_PML) = alpha_z

              endif
            endif

            if (region_CPML(ispec) == CPML_X_ONLY) then
              do i_sls = 1,2*N_SLS
                if (abs(alpha_x - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  alpha_x = tauinvnu(i_sls) + const_for_separation_two
                endif
                if (abs(alpha_x - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  call stop_the_code('error in separation of alpha_x, tauinvnu')
                endif
              enddo
              beta_x = alpha_x + d_x / K_x
              do i_sls = 1,2*N_SLS
                if (abs(beta_x - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  beta_x = tauinvnu(i_sls) + const_for_separation_two
                endif
                if (abs(beta_x - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  call stop_the_code('error in separation of beta_x, tauinvnu')
                endif
              enddo

              d_x = (beta_x - alpha_x) * K_x
              d_x_store(i,j,ispec_PML) = d_x
              alpha_x_store(i,j,ispec_PML) = alpha_x

            endif

            if (region_CPML(ispec) == CPML_Z_ONLY) then
              do i_sls = 1,2*N_SLS
                if (abs(alpha_z - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  alpha_z = tauinvnu(i_sls) + const_for_separation_two
                endif
                if (abs(alpha_z - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  call stop_the_code('error in separation of alpha_z, tauinvnu')
                endif
              enddo
              beta_z = alpha_z + d_z / K_z
              do i_sls = 1,2*N_SLS
                if (abs(beta_z - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  beta_z = tauinvnu(i_sls) + const_for_separation_two
                endif
                if (abs(beta_z - tauinvnu(i_sls)) < min_distance_between_CPML_parameter) then
                  call stop_the_code('error in separation of beta_z, tauinvnu')
                endif
              enddo

              d_z = (beta_z - alpha_z) * K_z
              d_z_store(i,j,ispec_PML) = d_z
              alpha_z_store(i,j,ispec_PML) = alpha_z

            endif
          endif
        enddo
      enddo
    endif
  enddo

  end subroutine prepare_timerun_attenuation_with_PML
!
!-------------------------------------------------------------------------------------
!
  subroutine tauinvnu_arange_from_lt_to_gt(siz,tauinvnu)

  implicit none
  integer, intent(in) :: siz
  double precision, dimension(siz), intent(inout) :: tauinvnu

  !local parameters
  integer :: i,j
  double precision :: temp

  do i= 1,siz
    do j= i+1,siz
      if (tauinvnu(i) > tauinvnu(j)) then
        temp = tauinvnu(i)
        tauinvnu(i) = tauinvnu(j)
        tauinvnu(j) =  temp
      endif
    enddo
  enddo

  end subroutine tauinvnu_arange_from_lt_to_gt
!
!-------------------------------------------------------------------------------------
!
  subroutine prepare_timerun_no_backward_reconstruction()

  use constants
  use specfem_par
  implicit none

  integer :: ier
  double precision :: sizeval

  no_backward_iframe = 0

  !checks if anything to do
  if (.not. NO_BACKWARD_RECONSTRUCTION) then
    allocate(no_backward_acoustic_buffer(1),no_backward_displ_buffer(1,1),no_backward_accel_buffer(1,1))
    return
  endif
  !safety checks
  if (time_stepping_scheme /= 1) call exit_MPI(myrank,'for NO_BACKWARD_RECONSTRUCTION, only Newmark scheme has implemented ')
  if (UNDO_ATTENUATION_AND_OR_PML) call exit_MPI(myrank, &
                                      'NO_BACKWARD_RECONSTRUCTION is not compatible with UNDO_ATTENUATION_AND_OR_PML')

  ! gets the number of frames to store/read in the NO BACKWARD RECONSTRUCTION
  ! database
  no_backward_nframes = 0
  do while (no_backward_nframes * NSTEP_BETWEEN_COMPUTE_KERNELS <= NSTEP )
    no_backward_nframes = no_backward_nframes + 1
  enddo
  no_backward_nframes = no_backward_nframes -1

  ! user output
  if (SAVE_FORWARD .or. SIMULATION_TYPE == 3) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Preparing NO_BACKWARD_RECONSTRUCTION :'
      write(IMAIN,*) '  number of frames to save :',no_backward_nframes
      if (any_acoustic) then
        sizeval = dble(nglob) * dble(no_backward_nframes) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
        write(IMAIN,*) '  size of No_backward_reconstruction file per slice = ', sngl(sizeval),'MB'
      endif
      if (any_elastic) then
        sizeval = dble(NDIM) * dble(nglob) * dble(no_backward_nframes) * &
                  dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
        if (APPROXIMATE_HESS_KL) sizeval = 2 * sizeval
        write(IMAIN,*) '  size of No_backward_reconstruction file per slice= ', sngl(sizeval),'MB'
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif ! myrank
  endif ! SAVE_FORWARD .or. SIMULATION_TYPE == 3

  ! deallocates arrays that won't be used with this mode
  if (SIMULATION_TYPE == 3) then
    if (any_acoustic) deallocate(b_potential_dot_dot_acoustic,b_potential_dot_acoustic,stat=ier)
    if (any_elastic) deallocate(b_veloc_elastic,b_accel_elastic,stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'error deallocating b_accel_elastic etc')
  endif

  ! allocates buffers for I/O
  if (SAVE_FORWARD .or. SIMULATION_TYPE == 3) then
    if (any_acoustic) then
      allocate(no_backward_acoustic_buffer(3*nglob),stat=ier)
    else
      allocate(no_backward_acoustic_buffer(1),stat=ier)
    endif
    if (any_elastic) then
      allocate(no_backward_displ_buffer(NDIM,nglob),stat=ier)
      if (APPROXIMATE_HESS_KL) then
        allocate(no_backward_accel_buffer(NDIM,nglob),stat=ier)
      else
        allocate(no_backward_accel_buffer(1,1),stat=ier)
      endif
    else
      allocate(no_backward_displ_buffer(1,1),no_backward_accel_buffer(1,1),stat=ier)
    endif
    if (ier /= 0 ) call exit_MPI(myrank,'error allocating no_backward_***_buffer')
  endif

  end subroutine prepare_timerun_no_backward_reconstruction
