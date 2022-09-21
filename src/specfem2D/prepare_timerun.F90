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

  implicit none

  ! local parameters
  double precision :: tCPU,tstart
  double precision, external :: wtime

  ! user output
  call synchronize_all()

  ! get MPI starting time
  tstart = wtime()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "Preparing timerun:"
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! reads binary database for setup
  if (setup_with_binary_database == 2) then
    ! updates external model and material properties in case
    call setup_mesh_material_properties()

    ! updates source/receivers
    call read_sources_receivers()
    if (SIMULATION_TYPE == 3) call setup_adjoint_sources()

    ! reads binary database
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

  ! Stacey preparation
  call prepare_timerun_Stacey()

  ! attenuation
  !! DK DK moved preparation of attenuation to before preparation of mass matrix
  !! DK DK because now that we have support for viscoacoustic fluids, we need to use
  !! DK DK the unrelaxed Kappa modulus in the fluid mass matrix when ATTENUATION_VISCOACOUSTIC is on
  !! DK DK and thus we need to prepare attenuation before preparing the mass matrix
  call prepare_attenuation()

  ! default preparation
  if (setup_with_binary_database /= 2) then

    ! prepares mass matrices
    call prepare_timerun_mass_matrix()

    ! check which GLL will be forced or not
    if (USE_ENFORCE_FIELDS) call build_forced()

    ! for adjoint kernel runs
    call prepare_timerun_adjoint()

    ! reads initial fields from external file if needed
    call prepare_timerun_initialfield

    ! compute the source time function and stores it in a text file
    call prepare_source_time_function()

    ! prepares noise simulations
    if (NOISE_TOMOGRAPHY /= 0) call prepare_timerun_noise()

    ! saves setup as binary database to skip these steps for repeated runs (mode == 2 used reading database)
    !
    ! daniel todo: note that most of the mesh setup steps should go into the mesher as is done in 3D versions.
    !              the 2D version historically does only basic mesh setups in the mesher, but then assigns most
    !              of the mesh population tasks (like xgenerate_databases in 3D Cartesian) to the solver.
    !              in future, one could move the time consuming meshing parts to the mesher and avoid
    !              this binary_database parameter handling.
    !
    if (setup_with_binary_database == 1) call save_binary_database()

  else
    ! reads binary setup database
    call read_binary_database_part2()
  endif

  ! movie simulations
  call prepare_timerun_movies()

  ! init specific to NO_BACKWARD_RECONSTRUCTION option
  call prepare_timerun_no_backward_reconstruction()

  ! prepares GPU arrays
  if (GPU_MODE) call prepare_GPU()

  !-------------------------------------------------------------

  ! check the mesh, stability and number of points per wavelength
  call check_grid()

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
  if (OUTPUT_ENERGY) then
    if (myrank == 0) then
      open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'energy.dat',status='unknown',action='write')
      ! file header
      write(IOUT_ENERGY,"('# Energy')")
      write(IOUT_ENERGY,"('#')")
      write(IOUT_ENERGY,"('# simulation setup: ')")
      if (ACOUSTIC_SIMULATION) then
        write(IOUT_ENERGY,"('#   using    acoustic wavefield')")
      endif
      if (ELASTIC_SIMULATION) then
        if (P_SV) then
          write(IOUT_ENERGY,"('#   using     elastic wavefield  -  P-SV waves')")
        else
          write(IOUT_ENERGY,"('#   using     elastic wavefield  -  SH waves')")
        endif
      endif
      if (POROELASTIC_SIMULATION) then
        write(IOUT_ENERGY,"('#   using poroelastic wavefield')")
      endif
      write(IOUT_ENERGY,"('#')")
      write(IOUT_ENERGY,"('#   NTSTEP_BETWEEN_OUTPUT_ENERGY : ',i7)") NTSTEP_BETWEEN_OUTPUT_ENERGY
      write(IOUT_ENERGY,"('#')")
      write(IOUT_ENERGY,"('# format:')")
      write(IOUT_ENERGY,"('#time  #E_kin(kinetic energy)  #E_pot(potential energy)  #E_tot(total energy)')")
    endif
  endif

  ! synchronizes all processes
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    ! elapsed time since beginning of preparation
    tCPU = wtime() - tstart
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for preparing timerun in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) '************'
    write(IMAIN,*) ' time loop'
    write(IMAIN,*) '************'
    select case(time_stepping_scheme)
    case (1)
      ! Newmark time scheme
      write(IMAIN,*) '              scheme:         Newmark'
    case (2)
      ! LDDRK
      write(IMAIN,*) '              scheme:         LDDRK      with',NSTAGE_TIME_SCHEME,'stages'
    case (3)
      ! RK4
      write(IMAIN,*) '              scheme:         RK4        with',NSTAGE_TIME_SCHEME,'stages'
    case (4)
      ! symplectic PEFRL
      write(IMAIN,*) '              scheme:         symplectic with',NSTAGE_TIME_SCHEME,'stages'
    case default
      call stop_the_code('Error time scheme not implemente yet')
    end select
    write(IMAIN,*)
    write(IMAIN,*) '           time step: ',sngl(DT),' s'
    write(IMAIN,*) 'number of time steps: ',NSTEP
    write(IMAIN,*) 'total simulated time: ',sngl(NSTEP*DT),' seconds'
    write(IMAIN,*) 'start time:',sngl(-t0),' seconds'
    write(IMAIN,*)

    ! flushes file buffer for main output file (IMAIN)
    call flush_IMAIN()
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_constants()

  use constants, only: IMAIN,HALF,CUSTOM_REAL
  use specfem_par

  implicit none

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Preparing timerun constants'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! defines coefficients of the Newmark time scheme
  !
  ! note: whenever possible, we will use deltat,.. values with CUSTOM_REAL precision to avoid implicit conversions
  !       when multiplying with CUSTOM_REAL wavefields etc.
  !
  !       DT will be kept in double precision and used when higher accuracy is needed or for non-critical computation
  !       (e.g., determining time for seismogram trace outputs)
  deltat = real(DT,kind=CUSTOM_REAL)
  deltatover2 = real(HALF * deltat,kind=CUSTOM_REAL)
  deltatsquareover2 = real(HALF * deltat * deltat,kind=CUSTOM_REAL)

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
      b_deltatover2 = real(HALF * b_deltat,kind=CUSTOM_REAL)
      b_deltatsquareover2 = real(HALF * b_deltat * b_deltat,kind=CUSTOM_REAL)
    endif
  else
    ! will not be used, but initialized
    b_deltat = 0._CUSTOM_REAL
    b_deltatover2 = 0._CUSTOM_REAL
    b_deltatsquareover2 = 0._CUSTOM_REAL
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_constants

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_mass_matrix()

  use constants, only: IMAIN
  use specfem_par

  implicit none

  ! local variable
#ifdef WITH_MPI
  integer :: n_sls_loc
#endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Preparing mass matrices'
    call flush_IMAIN()
  endif

  ! builds the global mass matrix
  call invert_mass_matrix_init()

#ifdef WITH_MPI
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

  subroutine prepare_timerun_movies()

  use constants, only: IMAIN,NOISE_MOVIE_OUTPUT
  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: ier

  ! sets movie flag
  if (output_postscript_snapshot .or. &
      output_color_image .or. &
      output_wavefield_dumps .or. &
      NOISE_MOVIE_OUTPUT) then
    MOVIE_SIMULATION = .true.
  else
    MOVIE_SIMULATION = .false.
  endif

  if (MOVIE_SIMULATION) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'Movie simulation:'
      if (output_postscript_snapshot) then
        write(IMAIN,*) '  postscript snapshots - PS image type  : ',imagetype_postscript
      endif
      if (output_color_image) then
        write(IMAIN,*) '  color images         - JPEG image type: ',imagetype_JPEG
      endif
      if (output_wavefield_dumps) then
        write(IMAIN,*) '  wavefield dumps      - image type     : ',imagetype_wavefield_dumps
      endif
      if (NOISE_MOVIE_OUTPUT) then
        write(IMAIN,*) '  noise movie'
      endif
      write(IMAIN,*) '  number of steps between outputs       : ',NTSTEP_BETWEEN_OUTPUT_IMAGES
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! to display the whole vector field (it needs to be computed from the potential in acoustic elements,
    ! thus it does not exist as a whole in case of simulations that contain some acoustic elements
    ! and it thus needs to be computed specifically for display purposes)
    allocate(vector_field_display(NDIM,nglob),stat=ier)
    ! when periodic boundary conditions are on, some global degrees of freedom are going to be removed,
    ! thus we need to set this array to zero otherwise some of its locations may contain random values
    ! if the memory is not cleaned
    vector_field_display(:,:) = 0.d0

    ! postscript images for grids and snapshots
    if (output_postscript_snapshot) then
      call prepare_timerun_postscripts()
    endif

    ! jpeg images
    if (output_color_image) then
      call prepare_timerun_image_coloring()
    endif

    ! wavefield dumps
    if (output_wavefield_dumps) then
      call prepare_timerun_wavefield_dumps()
    endif
  endif

  end subroutine prepare_timerun_movies

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_postscripts()

  use constants, only: IMAIN
  use specfem_par
  use specfem_par_movie
#ifdef WITH_MPI
  use constants, only: DISPLAY_COLORS,DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT
#endif

  implicit none

  ! local parameters
  integer :: i,j,ier

  ! for Lagrange interpolants
  double precision, external :: hgll, hglj

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
  if (.not. output_postscript_snapshot) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Preparing image coloring: postscripts'
    call flush_IMAIN()
  endif

  ! arrays for display images
  allocate(shape2D_display(NGNOD,pointsdisp,pointsdisp), &
           dershape2D_display(NDIM,NGNOD,pointsdisp,pointsdisp),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating shape arrays for display')
  shape2D_display(:,:,:) = 0.d0; dershape2D_display(:,:,:,:) = 0.d0

  ! computes shape functions and their derivatives for regular interpolated display grid
  do j = 1,pointsdisp
    do i = 1,pointsdisp
      xirec  = 2.d0*dble(i-1)/dble(pointsdisp-1) - 1.d0
      gammarec  = 2.d0*dble(j-1)/dble(pointsdisp-1) - 1.d0
      call define_shape_functions(shape2D_display(:,i,j),dershape2D_display(:,:,i,j),xirec,gammarec,NGNOD)
    enddo
  enddo

  ! for postscript snapshots
  ! arrays for display images as snapshot postscript images
  allocate(flagrange(NGLLX,pointsdisp))
  flagrange(:,:) = 0.d0

  if (AXISYM) then
    allocate(flagrange_GLJ(NGLJ,pointsdisp))
  else
    allocate(flagrange_GLJ(1,1))
  endif
  flagrange_GLJ(:,:) = 0.d0

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
  xinterp(:,:) = 0.d0; zinterp(:,:) = 0.d0
  Uxinterp(:,:) = 0.d0; Uzinterp(:,:) = 0.d0

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  allocating postscript image arrays'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocate arrays for postscript output
#ifdef WITH_MPI
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
  if (NGNOD == 4) then
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
  d2_coorg_send_ps_abs=4*num_abs_boundary_faces
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

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_postscripts

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_image_coloring()

  use constants, only: IMAIN
  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: i,j,k,iproc,ipixel
  integer :: ier

  ! checks if anything to do
  if (.not. output_color_image) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Preparing image coloring: jpeg'
    call flush_IMAIN()
  endif

  ! prepares dimension of image
  call prepare_color_image_init()

  ! for color images
  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  allocating color image arrays'
    call flush_IMAIN()
  endif

  ! allocate an array for image data
  allocate(image_color_data(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
  if (ier /= 0) call stop_the_code('error in an allocate statement 1')
  allocate(image_color_vp_display(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
  if (ier /= 0) call stop_the_code('error in an allocate statement 2')
  image_color_data(:,:) = 0.d0; image_color_vp_display(:,:) = 0.d0

  ! allocate an array for the grid point that corresponds to a given image data point
  allocate(iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
  if (ier /= 0) call stop_the_code('error in an allocate statement 3')
  allocate(copy_iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier)
  if (ier /= 0) call stop_the_code('error in an allocate statement 4')
  iglob_image_color(:,:) = 0; copy_iglob_image_color(:,:) = 0

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
  num_pixel_loc(:) = 0
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
#ifdef WITH_MPI
  allocate(nb_pixel_per_proc(0:NPROC-1))
  nb_pixel_per_proc(:) = 0
  call gather_all_singlei(nb_pixel_loc,nb_pixel_per_proc,NPROC)

  if (myrank == 0) then
    allocate(num_pixel_recv(maxval(nb_pixel_per_proc(:)),NPROC))
    allocate(data_pixel_recv(maxval(nb_pixel_per_proc(:))))
    num_pixel_recv(:,:) = 0; data_pixel_recv(:) = 0.d0
  endif
  allocate(data_pixel_send(nb_pixel_loc))
  data_pixel_send(:) = 0.d0

  if (NPROC > 1) then
    if (myrank == 0) then
      ! main collects
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

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_image_coloring

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_wavefield_dumps()

  use constants, only: IMAIN
  use specfem_par
  use specfem_par_movie

  implicit none

  ! wavefield dump
  integer :: d1_dump_send, d2_dump_send, d1_dump_recv, d2_dump_recv
  integer :: ier

  ! checks if anything to do
  if (.not. output_wavefield_dumps) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Preparing wavefield dumps:'
    write(IMAIN,*) '  allocating dump arrays'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocate arrays for wavefield dump
  d1_dump_recv = 2
  d1_dump_send = 2

  d2_dump_send = nspec*NGLLX*NGLLZ
  call max_all_all_i(d2_dump_send,d2_dump_recv)

  ! allocates temporary arrays for wavefield outputs
  allocate(dump_send(d1_dump_send, d2_dump_send))
  allocate(dump_recv(d1_dump_recv, d2_dump_recv))
  dump_send(:,:) = 0.d0
  dump_recv(:,:) = 0.d0

  allocate(dump_duplicate_send(d2_dump_send))
  allocate(dump_duplicate_recv(d2_dump_recv))
  dump_duplicate_send(:) = .false.
  dump_duplicate_recv(:) = .false.

  allocate(dump_recv_counts(0:NPROC-1),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating wavefield_dumps arrays')
  dump_recv_counts(:) = 0

  this_is_the_first_time_we_dump = .true.

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_wavefield_dumps

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

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    if (SAVE_FORWARD) then
      write(IMAIN,*) 'Preparing save forward simulation:'
    else
      write(IMAIN,*) 'Preparing kernel simulation:'
    endif
    write(IMAIN,*) '  estimated minimum period resolved by mesh        : ',sngl(mesh_T_min)
    write(IMAIN,*) '  estimated number of time steps for minimum period: ',int(mesh_T_min / DT)
    ! note: to reconstruct kernels, the stepping to approximate kernel depends on both source and adjoint source frequency,
    !       as well as the minimum period resolved by the mesh.
    !       by comparison, we see that close-enough reconstructions (within ~90% accuracy) need about
    !       a stepping of the size of a quarter of the minimum period resolved.
    !       this lower estimate can then be compared against the NTSTEP_BETWEEN_COMPUTE_KERNELS setting.
    write(IMAIN,*) '  estimated steps between compute kernels (for a fair reconstruction): ',int(mesh_T_min / DT / 4.0)
    write(IMAIN,*)
    write(IMAIN,*) '  number of steps between compute kernels: ',NTSTEP_BETWEEN_COMPUTE_KERNELS
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! prepares kernels
  if (SIMULATION_TYPE == 3) then
    ! user output
    if (myrank == 0) then
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


  use constants, only: IMAIN
  use specfem_par

  implicit none

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  initializing adjoint sensitivity kernels'
    call flush_IMAIN()
  endif

  ! initalizes sensitivity kernel arrays
  ! elastic domains
  if (any_elastic) then
    ! initializes
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
    ! safety checks
    if (.not. save_ASCII_kernels) call stop_the_code('poroelastic simulations must use save_ASCII_kernels')
    if (GPU_MODE) call stop_the_code('poroelastic kernel output not implemented on GPUs yet')

    ! initializes
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
    ! initializes
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
  if (any_poroelastic) &
    call stop_the_code('Initial field currently implemented for purely elastic simulation only')

  ! Calculation of the initial field for a plane wave
  ! calculates initial plane wave coefficients
  call prepare_initial_field(cploc,csloc)

  ! special case for Rayleigh waves, SV waves above critical angle
  if (over_critical_angle) then

    ! allocates boundaries
    allocate(left_bound(num_abs_boundary_faces*NGLLX))
    allocate(right_bound(num_abs_boundary_faces*NGLLX))
    allocate(bot_bound(num_abs_boundary_faces*NGLLZ))

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
                              QKappa_attenuationcoef(1),source_type(1), &
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
    if (any_elastic) then
      write(IMAIN,*) 'Max norm of initial elastic displacement = ', &
                      maxval(sqrt(displ_elastic(1,:)**2 + displ_elastic(2,:)**2))
    endif
    if (any_acoustic) then
      write(IMAIN,*) 'Max norm of initial acoustic displacement = ', &
                      maxval(potential_acoustic(:))
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine prepare_timerun_initialfield

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_noise()

! for noise simulations

  use constants, only: NGLLX,NGLLZ,NDIM,IMAIN,NOISE_MOVIE_OUTPUT,TWO_THIRDS,OUTPUT_FILES

  use specfem_par, only: myrank,NSTEP,nglob,nspec,ibool,coord, &
                         rhostore,rho_vpstore,rho_vsstore, &
                         NOISE_TOMOGRAPHY

  use specfem_par_noise

  implicit none

  integer :: i,j,iglob,ispec,ier
  double precision :: rhol,vs,vp

  ! checks if anything to do
  if (NOISE_TOMOGRAPHY <= 0) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Preparing noise simulation'
    call flush_IMAIN()
  endif

  ! allocates arrays for noise tomography
  allocate(noise_sourcearray(NDIM,NGLLX,NGLLZ,NSTEP), &
           mask_noise(nglob), &
           noise_surface_movie_y_or_z(nglob),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating noise arrays')
  noise_sourcearray(:,:,:,:) = 0._CUSTOM_REAL
  mask_noise(:) = 0._CUSTOM_REAL
  noise_surface_movie_y_or_z(:) = 0._CUSTOM_REAL

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  reading noise parameters'
    call flush_IMAIN()
  endif

  !read in parameters for noise tomography
  call read_parameters_noise()

  if (NOISE_TOMOGRAPHY == 1) then
    ! creates generating noise source
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
    open(unit=504,file=trim(OUTPUT_FILES)//'model_rho_vp_vs',status='unknown',action='write')
    do ispec = 1, nspec
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          rhol = dble(rhostore(i,j,ispec))
          vs = dble(rho_vsstore(i,j,ispec)/rhol)
          vp = dble(rho_vpstore(i,j,ispec)/rhol)

          write(504,'(1pe11.3,1pe11.3,1pe11.3,1pe11.3,1pe11.3)') coord(1,iglob), coord(2,iglob), rhol, vp, vs
        enddo
      enddo
    enddo
    close(504)

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
      noise_output_array(:,:) = 0._CUSTOM_REAL
      noise_output_rhokl(:) = 0._CUSTOM_REAL
    endif

  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_timerun_noise

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
    ! dummy
    allocate(no_backward_acoustic_buffer(1),no_backward_displ_buffer(1,1),no_backward_accel_buffer(1,1))
    return
  endif

  !safety checks
  if (time_stepping_scheme /= 1) &
    call exit_MPI(myrank,'for NO_BACKWARD_RECONSTRUCTION, only Newmark scheme has implemented ')
  if (UNDO_ATTENUATION_AND_OR_PML) &
    call exit_MPI(myrank,'NO_BACKWARD_RECONSTRUCTION is not compatible with UNDO_ATTENUATION_AND_OR_PML')

  ! gets the number of frames to store/read in the NO BACKWARD RECONSTRUCTION
  ! database
  no_backward_nframes = 0
  do while (no_backward_nframes * NTSTEP_BETWEEN_COMPUTE_KERNELS <= NSTEP )
    no_backward_nframes = no_backward_nframes + 1
  enddo
  no_backward_nframes = no_backward_nframes -1

  ! user output
  if (SAVE_FORWARD .or. SIMULATION_TYPE == 3) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Preparing NO_BACKWARD_RECONSTRUCTION :'
      write(IMAIN,*) '  number of steps between compute kernels: ',NTSTEP_BETWEEN_COMPUTE_KERNELS
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

    ! checks integer overflow
    if (any_acoustic) then
      !offset = CUSTOM_REAL * nglob * (no_backward_Nframes - no_backward_iframe) + 1
      ! checks
      if (dble(nglob) * dble(no_backward_Nframes - 1) > 2147483646.d0 / dble(CUSTOM_REAL) ) &
        call exit_MPI(myrank,'Error no_backward buffer offset might exceed integer limit')
    endif
    if (any_elastic) then
      if (APPROXIMATE_HESS_KL) then
        !offset = 2 * CUSTOM_REAL * (NDIM*nglob) * (no_backward_Nframes - no_backward_iframe) + 1
        ! checks
        if (2.d0 * dble(NDIM) * dble(nglob) * dble(no_backward_Nframes - 1) > 2147483646.d0 / dble(CUSTOM_REAL) ) &
          call exit_MPI(myrank,'Error no_backward buffer offset might exceed integer limit')
      else
        !offset = CUSTOM_REAL * (NDIM*nglob) * (no_backward_Nframes - no_backward_iframe) + 1
        ! checks
        if (dble(NDIM) * dble(nglob) * dble(no_backward_Nframes - 1) > 2147483646.d0 / dble(CUSTOM_REAL) ) &
          call exit_MPI(myrank,'Error no_backward buffer offset might exceed integer limit')
      endif
    endif
  endif ! SAVE_FORWARD .or. SIMULATION_TYPE == 3

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

!
!-------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_Stacey()

  use constants, only: IEDGE1,IEDGE2,IEDGE3,IEDGE4,IMAIN
  use specfem_par

  implicit none

  integer :: i,j,ispec,ispecabs,ier,iedge,num_all
  real(kind=CUSTOM_REAL) :: xxi,zxi,xgamma,zgamma,jacobian1D
  logical :: has_abs_edge

  ! sets up arrays for stacey boundary routines
  if (STACEY_ABSORBING_CONDITIONS) then
    ! gets total number of boundary faces/edges
    call sum_all_i(num_abs_boundary_faces,num_all)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Preparing Stacey boundaries'
      write(IMAIN,*) '  total number of absorbing boundary faces/edges: ',num_all
      call flush_IMAIN()
    endif

    ! stacey boundaries
    if (num_abs_boundary_faces > 0) then
      allocate(abs_boundary_ij(2,NGLLX,num_abs_boundary_faces), &
               abs_boundary_jacobian1Dw(NGLLX,num_abs_boundary_faces), &
               abs_boundary_normal(NDIM,NGLLX,num_abs_boundary_faces),stat=ier)
      if (ier /= 0) stop 'error allocating array abs_boundary_ij etc.'
    else
      ! dummy
      allocate(abs_boundary_ij(1,1,1), &
               abs_boundary_jacobian1Dw(1,1), &
               abs_boundary_normal(1,1,1),stat=ier)
      if (ier /= 0) stop 'error allocating array abs_boundary_ij etc.'
    endif
    abs_boundary_ij(:,:,:) = 0; abs_boundary_jacobian1Dw(:,:) = 0.0_CUSTOM_REAL
    abs_boundary_normal(:,:,:) = 0.0_CUSTOM_REAL

    ! needed for gpu boundary array storage
    if (num_abs_boundary_faces > 0) then
      allocate(edge_abs(num_abs_boundary_faces),stat=ier)
      if (ier /= 0) stop 'error allocating array edge_abs etc.'
    else
      ! dummy
      allocate(edge_abs(1),stat=ier)
      if (ier /= 0) stop 'error allocating array edge_abs etc.'
    endif
    edge_abs(:) = 0

    do ispecabs = 1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(ispecabs)

      ! checks if edge has correct absorbing flags
      has_abs_edge = .false.
      do iedge = 1,4
        if (codeabs(iedge,ispecabs)) then
          if (has_abs_edge) then
            ! boundary face has already been set, should not occur
            print *,'Error: absorbing boundary element ',ispec,'has edge ',ispecabs, &
                    ' with multiple flags l/r/bottom/top:',codeabs(:,ispecabs)
            stop 'Invalid absorbing edge found'
          else
            has_abs_edge = .true.
          endif
        endif
      enddo
      ! must have at least one absorbing flag set
      if (.not. has_abs_edge) then
        print *,'Error: absorbing boundary element ',ispec,'has edge ',ispecabs, &
                ' with multiple flags l/r/bottom/top:',codeabs(:,ispecabs)
        stop 'Invalid absorbing edge with no flags found'
      endif

      !--- left absorbing boundary
      if (codeabs(IEDGE4,ispecabs)) then
        i = 1
        do j = 1,NGLLZ
          abs_boundary_ij(1,j,ispecabs) = i
          abs_boundary_ij(2,j,ispecabs) = j

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)

          abs_boundary_normal(1,j,ispecabs) = - zgamma / jacobian1D
          abs_boundary_normal(2,j,ispecabs) = + xgamma / jacobian1D

          abs_boundary_jacobian1Dw(j,ispecabs) = jacobian1D * wzgll(j)

          edge_abs(ispecabs) = 4
        enddo
        ! note: for elastic elements, edges are looped over the full range, i.e. from 1 to NGLLX/NGLLY/NGLLZ.
        !       for acoustic & poroelastic elements however, common points on coupling interfaces are removed
        !       by looping only over ibegin_edge* to iend_edge* ranges.
        !       here, we assign a value of NGLL* + 1 to indicate such duplicate points on edges.
        !       this will avoid the need to check if the edge is in an acoustic/elastic/poroelastic element.
        !
        ! uses NGLLZ+1 to indicate points which duplicate contributions and can be left out
        if (ispec_is_acoustic(ispec)) then
          if (ibegin_edge4(ispecabs) == 2) abs_boundary_ij(2,1,ispecabs) = NGLLZ+1
          if (iend_edge4(ispecabs) == NGLLZ-1) abs_boundary_ij(2,NGLLZ,ispecabs) = NGLLZ+1
        endif
        if (ispec_is_poroelastic(ispec)) then
          if (ibegin_edge4_poro(ispecabs) == 2) abs_boundary_ij(2,1,ispecabs) = NGLLZ+1
          if (iend_edge4_poro(ispecabs) == NGLLZ-1) abs_boundary_ij(2,NGLLZ,ispecabs) = NGLLZ+1
        endif
      endif

      !--- right absorbing boundary
      if (codeabs(IEDGE2,ispecabs)) then
        i = NGLLX
        do j = 1,NGLLZ
          abs_boundary_ij(1,j,ispecabs) = i
          abs_boundary_ij(2,j,ispecabs) = j

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)

          abs_boundary_normal(1,j,ispecabs) = + zgamma / jacobian1D
          abs_boundary_normal(2,j,ispecabs) = - xgamma / jacobian1D

          abs_boundary_jacobian1Dw(j,ispecabs) = jacobian1D * wzgll(j)

          edge_abs(ispecabs) = 2
        enddo
        ! uses NGLLZ+1 to indicate points which duplicate contributions and can be left out
        if (ispec_is_acoustic(ispec)) then
          if (ibegin_edge2(ispecabs) == 2) abs_boundary_ij(2,1,ispecabs) = NGLLZ+1
          if (iend_edge2(ispecabs) == NGLLZ-1) abs_boundary_ij(2,NGLLZ,ispecabs) = NGLLZ+1
        endif
        if (ispec_is_poroelastic(ispec)) then
          if (ibegin_edge2_poro(ispecabs) == 2) abs_boundary_ij(2,1,ispecabs) = NGLLZ+1
          if (iend_edge2_poro(ispecabs) == NGLLZ-1) abs_boundary_ij(2,NGLLZ,ispecabs) = NGLLZ+1
        endif
      endif

      !--- bottom absorbing boundary
      if (codeabs(IEDGE1,ispecabs)) then
        j = 1
        do i = 1,NGLLX
          abs_boundary_ij(1,i,ispecabs) = i
          abs_boundary_ij(2,i,ispecabs) = j

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)

          abs_boundary_normal(1,i,ispecabs) = + zxi / jacobian1D
          abs_boundary_normal(2,i,ispecabs) = - xxi / jacobian1D

          abs_boundary_jacobian1Dw(i,ispecabs) = jacobian1D * wxgll(i)

          edge_abs(ispecabs) = 1
        enddo
        ! uses NGLLX+1 to indicate points which duplicate contributions and can be left out
        if (ispec_is_acoustic(ispec)) then
          if (ibegin_edge1(ispecabs) == 2) abs_boundary_ij(1,1,ispecabs) = NGLLX+1
          if (iend_edge1(ispecabs) == NGLLX-1) abs_boundary_ij(1,NGLLX,ispecabs) = NGLLX+1
        endif
        if (ispec_is_poroelastic(ispec)) then
          if (ibegin_edge1_poro(ispecabs) == 2) abs_boundary_ij(1,1,ispecabs) = NGLLX+1
          if (iend_edge1_poro(ispecabs) == NGLLX-1) abs_boundary_ij(1,NGLLX,ispecabs) = NGLLX+1
        endif
        ! exclude duplicate points on corners (elements with left/right and bottom/top edge on boundary)
        if (codeabs_corner(1,ispecabs)) abs_boundary_ij(1,1,ispecabs) = NGLLX+1
        if (codeabs_corner(2,ispecabs)) abs_boundary_ij(1,NGLLX,ispecabs) = NGLLX+1
      endif

      !--- top absorbing boundary
      if (codeabs(IEDGE3,ispecabs)) then
        j = NGLLZ
        do i = 1,NGLLX
          abs_boundary_ij(1,i,ispecabs) = i
          abs_boundary_ij(2,i,ispecabs) = j

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)

          abs_boundary_normal(1,i,ispecabs) = - zxi / jacobian1D
          abs_boundary_normal(2,i,ispecabs) = + xxi / jacobian1D

          abs_boundary_jacobian1Dw(i,ispecabs) = jacobian1D * wxgll(i)

          edge_abs(ispecabs) = 3
        enddo
        ! uses NGLLX+1 to indicate points which duplicate contributions and can be left out
        if (ispec_is_acoustic(ispec)) then
          if (ibegin_edge3(ispecabs) == 2) abs_boundary_ij(1,1,ispecabs) = NGLLX+1
          if (iend_edge3(ispecabs) == NGLLX-1) abs_boundary_ij(1,NGLLX,ispecabs) = NGLLX+1
        endif
        if (ispec_is_poroelastic(ispec)) then
          if (ibegin_edge3_poro(ispecabs) == 2) abs_boundary_ij(1,1,ispecabs) = NGLLX+1
          if (iend_edge3_poro(ispecabs) == NGLLX-1) abs_boundary_ij(1,NGLLX,ispecabs) = NGLLX+1
        endif
        ! exclude duplicate points on corners (elements with left/right and bottom/top edge on boundary)
        if (codeabs_corner(3,ispecabs)) abs_boundary_ij(1,1,ispecabs) = NGLLX+1
        if (codeabs_corner(4,ispecabs)) abs_boundary_ij(1,NGLLX,ispecabs) = NGLLX+1
      endif
    enddo
  else
    ! dummy allocation
    allocate(abs_boundary_ij(1,1,1), &
             abs_boundary_jacobian1Dw(1,1), &
             abs_boundary_normal(1,1,1), &
             edge_abs(1),stat=ier)
    if (ier /= 0) stop 'error allocating dummy array abs_boundary_ij etc.'
  endif ! STACEY_ABSORBING_CONDITIONS

  ! allocates arrays for backward field
  if (anyabs .and. STACEY_ABSORBING_CONDITIONS .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3)) then
    ! files to save absorbed waves needed to reconstruct backward wavefield for adjoint method
    ! elastic domains
    if (any_elastic) then
      allocate(b_absorb_elastic_left(NDIM,NGLLZ,nspec_left,NSTEP))
      allocate(b_absorb_elastic_right(NDIM,NGLLZ,nspec_right,NSTEP))
      allocate(b_absorb_elastic_bottom(NDIM,NGLLX,nspec_bottom,NSTEP))
      allocate(b_absorb_elastic_top(NDIM,NGLLX,nspec_top,NSTEP))
      b_absorb_elastic_left(:,:,:,:) = 0.0_CUSTOM_REAL; b_absorb_elastic_right(:,:,:,:) = 0.0_CUSTOM_REAL
      b_absorb_elastic_bottom(:,:,:,:) = 0.0_CUSTOM_REAL; b_absorb_elastic_top(:,:,:,:) = 0.0_CUSTOM_REAL
    endif
    ! poroelastic domains
    if (any_poroelastic) then
      allocate(b_absorb_poro_s_left(NDIM,NGLLZ,nspec_left,NSTEP))
      allocate(b_absorb_poro_s_right(NDIM,NGLLZ,nspec_right,NSTEP))
      allocate(b_absorb_poro_s_bottom(NDIM,NGLLX,nspec_bottom,NSTEP))
      allocate(b_absorb_poro_s_top(NDIM,NGLLX,nspec_top,NSTEP))
      allocate(b_absorb_poro_w_left(NDIM,NGLLZ,nspec_left,NSTEP))
      allocate(b_absorb_poro_w_right(NDIM,NGLLZ,nspec_right,NSTEP))
      allocate(b_absorb_poro_w_bottom(NDIM,NGLLX,nspec_bottom,NSTEP))
      allocate(b_absorb_poro_w_top(NDIM,NGLLX,nspec_top,NSTEP))
      b_absorb_poro_s_left(:,:,:,:) = 0.0_CUSTOM_REAL; b_absorb_poro_s_right(:,:,:,:) = 0.0_CUSTOM_REAL
      b_absorb_poro_s_bottom(:,:,:,:) = 0.0_CUSTOM_REAL; b_absorb_poro_s_top(:,:,:,:) = 0.0_CUSTOM_REAL
      b_absorb_poro_w_left(:,:,:,:) = 0.0_CUSTOM_REAL; b_absorb_poro_w_right(:,:,:,:) = 0.0_CUSTOM_REAL
      b_absorb_poro_w_bottom(:,:,:,:) = 0.0_CUSTOM_REAL; b_absorb_poro_w_top(:,:,:,:) = 0.0_CUSTOM_REAL
    endif
    ! acoustic domains
    if (any_acoustic) then
      allocate(b_absorb_acoustic_left(NGLLZ,nspec_left,NSTEP))
      allocate(b_absorb_acoustic_right(NGLLZ,nspec_right,NSTEP))
      allocate(b_absorb_acoustic_bottom(NGLLX,nspec_bottom,NSTEP))
      allocate(b_absorb_acoustic_top(NGLLX,nspec_top,NSTEP))
      b_absorb_acoustic_left(:,:,:) = 0.0_CUSTOM_REAL; b_absorb_acoustic_right(:,:,:) = 0.0_CUSTOM_REAL
      b_absorb_acoustic_bottom(:,:,:) = 0.0_CUSTOM_REAL; b_absorb_acoustic_top(:,:,:) = 0.0_CUSTOM_REAL
    endif
  else
    ! dummy arrays
    ! elastic domains
    if (.not. allocated(b_absorb_elastic_left)) then
      allocate(b_absorb_elastic_left(1,1,1,1))
      allocate(b_absorb_elastic_right(1,1,1,1))
      allocate(b_absorb_elastic_bottom(1,1,1,1))
      allocate(b_absorb_elastic_top(1,1,1,1))
    endif
    ! poroelastic domains
    if (.not. allocated(b_absorb_poro_s_left)) then
      allocate(b_absorb_poro_s_left(1,1,1,1))
      allocate(b_absorb_poro_s_right(1,1,1,1))
      allocate(b_absorb_poro_s_bottom(1,1,1,1))
      allocate(b_absorb_poro_s_top(1,1,1,1))
      allocate(b_absorb_poro_w_left(1,1,1,1))
      allocate(b_absorb_poro_w_right(1,1,1,1))
      allocate(b_absorb_poro_w_bottom(1,1,1,1))
      allocate(b_absorb_poro_w_top(1,1,1,1))
    endif
    ! acoustic domains
    if (.not. allocated(b_absorb_acoustic_left)) then
      allocate(b_absorb_acoustic_left(1,1,1))
      allocate(b_absorb_acoustic_right(1,1,1))
      allocate(b_absorb_acoustic_bottom(1,1,1))
      allocate(b_absorb_acoustic_top(1,1,1))
    endif
  endif

  ! done
  call synchronize_all()

  end subroutine prepare_timerun_Stacey
