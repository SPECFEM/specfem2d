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

module moving_sources_par
  ! This module has been created because it allows for the use of optional arguments
  ! (only inside modules)
  ! Also it is more tidy this way
  ! Read the comments if you are interested in how moving sources work

  use constants, only: NDIM,NGLLX,NGLLZ,IMAIN,HUGEVAL,TINYVAL,NUM_ITER, &
                       USE_BEST_LOCATION_FOR_SOURCE, &
                       IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC

  use specfem_par, only: SOURCE_IS_MOVING, AXISYM, is_on_the_axis, xiglj, nspec, &
                         nglob, ibool, coord, xigll, zigll, NGNOD, npgeo, knods, coorg, &
                         ispec_is_acoustic, ispec_is_elastic, ispec_is_poroelastic, &
                         CUSTOM_REAL

  implicit none

  ! This type is convenient when dealing with arrays of variable size
  type :: custom_array_pointer_type
      integer, pointer, dimension(:) :: array
      logical :: allocated = .false.
      integer :: size = 0
  end type custom_array_pointer_type

  integer, dimension(:), allocatable :: nsources_local_moving
  integer, dimension(:,:), allocatable :: ispec_selected_source_moving

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: sourcearrays_moving
  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: source_time_function_moving

  integer, parameter :: IMOV = 1684  ! File unit for moving sources databases

  contains

  !
  !-----------------------------------------------------------------------------------------
  !

  subroutine init_moving_sources_GPU()
  !----
  !---- Compute the sourcearrays for all timesteps and send them to the device
  !---- It is really similar to subroutine init_source but
  !---- in the case of moving sources this routine if a little bit more complicated
  !---- because we have to loop over all the time steps which can be very costly if
  !---- we do it dumbly.
  !---- Hopefully, the source not moving too quickly, we know that it will stay in the same
  !---- element (or the ones around) for some time so we can avoid looping over
  !---- all the elements of the mesh too often
  !---- Also, in the case of GPU we have to compute all the positions
  !---- and source arrays in advance (for all the sources and time steps)
  !---- to avoid exchanging with the device afterward
  !----
  !---- See old version in compute_gpu_acoutic.f90
  !---- (subroutine compute_add_sources_acoustic_GPU_moving_sources_old)
  !----
  !---- Please let the comments here, they can be very useful
  !----

  use constants, only: NGLLX, NGLLZ, TINYVAL, IMAIN, IN_DATA_FILES, &
                       USE_BEST_LOCATION_FOR_SOURCE, HUGEVAL

  use specfem_par
  use specfem_par_gpu

  implicit none

  ! Local variables
  double precision :: xsrc, zsrc, time_val, t_used
  integer :: i_source, it_l, i_stage_loc, i_adjacent, ier
  integer :: ispec_source_first_guess,ispec_selected_source_local, ispec_adjacent
  integer, dimension(NSOURCES) :: ispec_first_guess_vec
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  ! For debugging
  !double precision :: start_time, stop_time
  !integer :: iproc

  logical :: file_exists, is_force_source
  logical :: src_in_first_guess_element, src_in_adjacent_element
  logical, dimension(NSOURCES) :: source_belonged_to_this_rank_all
  logical, dimension(NSOURCES) :: reset_source,any_reset_source

  character(len=MAX_STRING_LEN) :: pathToMovingDatabaseDir, outputname, pathToMovingDatabase

  ! To store list of elements (of variable size):
  type(custom_array_pointer_type), dimension(NSOURCES) :: pt_first_guess

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Moving sources on GPU:'
    write(IMAIN,*) '  write moving database = ',write_moving_sources_database
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! Allocations
  allocate(nsources_local_moving(NSTEP), &
           sourcearrays_moving(NDIM,NGLLX,NGLLZ,NSOURCES,NSTEP), &
           source_time_function_moving(NSOURCES,NSTEP), &
           ispec_selected_source_moving(NSOURCES,NSTEP),stat=ier)
  if (ier /= 0) stop 'Error allocating moving sources arrays'

  ! Initializations
  nsources_local_moving(:) = 0
  sourcearrays_moving(:,:,:,:,:) = 0.0_CUSTOM_REAL
  source_time_function_moving(:,:) = 0.0_CUSTOM_REAL
  ispec_selected_source_moving(:,:) = 0  ! zero indicates not-a-valid ispec
  ispec_first_guess_vec(:) = 0

  ! initializes slice which hold source
  islice_selected_source(:) = -1
  iglob_source(:) = 0
  ! Initialize all the flags
  file_exists = .false.

  source_belonged_to_this_rank_all(:) = .true.
  reset_source(:) = .false.
  any_reset_source(:) = .false.

  ! checks
  if (setup_with_binary_database /= 0 ) &
    call exit_MPI(myrank,'setup_with_binary_database not available with moving sources yet')

  ! Warnings and checks
  if (myrank == 0) then
    ! citation
    write(IMAIN,*)
    write(IMAIN,*) '****************************************************************************************'
    write(IMAIN,*) 'Your are using moving source capabilities. Please cite:'
    write(IMAIN,*) 'Bottero (2018) Full-wave numerical simulation of T-waves and of moving acoustic sources'
    write(IMAIN,*) 'PhD thesis'
    write(IMAIN,*) 'https://tel.archives-ouvertes.fr/tel-01893011'
    write(IMAIN,*) '****************************************************************************************'
    write(IMAIN,*)
    do i_source = 1,NSOURCES
      ! timing warning
      if ((abs(tshift_src(i_source)) > 0.0d0) .or. (abs(t0) > 0.0d0)) then
        write(IMAIN,*) 'Source #',i_source
        write(IMAIN,*) ' !! BEWARE !! Parameters tshift and/or t0 are used with moving source !'
        write(IMAIN,*) ' The time step for the moving source is: '
        write(IMAIN,*) '    t_used = (it_l-1)*DT-t0-tshift_src(i_source)'
        write(IMAIN,*) ' And the source position is calculated like:'
        write(IMAIN,*) '  xsrc = x_source + vx_source*t_used'
        write(IMAIN,*)
      endif
    enddo
  endif

  ! Check that the sources are all force sources
  do i_source = 1,NSOURCES
    if (source_type(i_source) /= 1) &
      call exit_MPI(myrank,'Only force sources are implemented for moving sources on GPUs')
  enddo
  is_force_source = .true.

  if (write_moving_sources_database) then
    ! Saves data in a binary file
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_movingDatabase.bin'
    ! pathToMovingDatabaseDir = '/directory/you/want/'  !!! Must end with '/' !!!
    pathToMovingDatabaseDir = trim(IN_DATA_FILES) ! OR '.'
    pathToMovingDatabase = trim(pathToMovingDatabaseDir)//trim(outputname)
    inquire(file=pathToMovingDatabase, exist=file_exists)
  else
    pathToMovingDatabase = './this_file_does_not_exist_for_sure_and_this_should_never_happen_but_we_never_know'
  endif

  ! Generate the moving source databases (can be expensive)
  if ((.not. file_exists) .or. (.not. write_moving_sources_database)) then
    if (myrank == 0) then
      write(IMAIN,*) '    Generating moving source databases ...'
      if (.not. write_moving_sources_database) then
        write(IMAIN,*) '      If this step takes too much time in your case, you may want to turn in the Par_file:'
        write(IMAIN,*) '         write_moving_sources_database = .true.'
        write(IMAIN,*) '      In this case you may also want to goto subroutine init_moving_sources and read the comments.'
        ! There:
        ! Turn write_moving_sources_database to .true. in the case of GPU computing
        ! if the generation of moving source databases takes a long time.
        ! The simulation is done in two steps:
        !  1. you run the code and it writes the databases to file (in DATA folder by default).
        !  2. then you rerun the code and it will read the databases in there directly saving a lot of time.
      endif
      call flush_IMAIN()
    endif

    ! Loop over the time steps to precompute all source positions
    do it_l = 1, NSTEP
      ! user output
      if (myrank == 0) then
        if (mod(it_l-1, int(NSTEP/50)) == 0) &
          write(IMAIN,*) "        Computing source(s) position(s) for timestep",it_l-1,"/",NSTEP,"..."
          call flush_IMAIN()
      endif

      ! ! For debugging:
      ! print *,"debug: init_moving_sources_GPU: source(s) position(s) for timestep",it_l," out of ",NSTEP
      ! call cpu_time(stop_time)
      ! call synchronize_all()
      ! if (myrank == 0) then
      ! call sleep(1)
      !   print *
      !   print *
      !   print *,myrank," Duration of this time step: ",stop_time - start_time, "seconds"
      ! endif
      ! call synchronize_all()
      ! if (myrank == 1) then
      !   print *,myrank," Duration of this time step: ",stop_time - start_time, "seconds"
      ! endif
      ! call synchronize_all()
      ! if (myrank == 2) then
      !   print *,myrank," Duration of this time step: ",stop_time - start_time, "seconds"
      ! endif
      ! call synchronize_all()
      ! if (myrank == 3) then
      !   print *,myrank," Duration of this time step: ",stop_time - start_time, "seconds"
      ! endif
      ! call synchronize_all()
      ! call cpu_time(start_time)

      ! note: the following routine and arrays assume there is only a single stage for the time scheme.
      if (NSTAGE_TIME_SCHEME /= 1) stop 'Invalid number of stages of time scheme for moving sources'

      ! counts sources in this process slice for this time step
      nsources_local = 0

      do i_stage_loc = 1, NSTAGE_TIME_SCHEME  ! For now NSTAGE_TIME_SCHEME == 1 only because only Newmark time scheme are supported

        ! current time
        if (time_stepping_scheme == 1) then
          ! Newmark
          time_val = (it_l-1)*DT
        else
          call exit_MPI(myrank,'Only Newmark time scheme is implemented for moving sources (1)')
        endif

        ! checks if source needs to be re-located in all slices
        ! Loop on the sources
        do i_source = 1,NSOURCES
          ! No need to do the following if amplitude is zero
          if (abs(source_time_function(i_source,it_l,i_stage_loc)) <= TINYVAL) cycle

          ! updated source position
          ! time
          t_used = (time_val - t0 - tshift_src(i_source))
          ! Moves and re-locates the sources along x and z-axis
          xsrc = x_source(i_source) + vx_source(i_source) * t_used
          zsrc = z_source(i_source) + vz_source(i_source) * t_used

          ! if (it_l > 51) then
          !   call synchronize_all()
          !   if (myrank == 0) then
          !     print *
          !     print *,"it:",it_l," xsrc :", xsrc
          !     print *
          !     print *,myrank," source_belonged_to_this_rank :", source_belonged_to_this_rank_all(i_source), &
          !             "reset:", reset_source(i_source), &
          !             "ispec_source_first_guess:",ispec_source_first_guess
          !   endif
          !   call synchronize_all()
          !   if (myrank == 1) then
          !     print *,myrank," source_belonged_to_this_rank :", source_belonged_to_this_rank_all(i_source), &
          !             "reset:", reset_source(i_source), &
          !             "ispec_source_first_guess:",ispec_source_first_guess
          !   endif
          !   call synchronize_all()
          !   if (myrank == 2) then
          !     print *,myrank," source_belonged_to_this_rank :", source_belonged_to_this_rank_all(i_source), &
          !             "reset:", reset_source(i_source), &
          !             "ispec_source_first_guess:",ispec_source_first_guess
          !   endif
          !   call synchronize_all()
          !   if (myrank == 3) then
          !     print *,myrank," source_belonged_to_this_rank :", source_belonged_to_this_rank_all(i_source), &
          !             "reset:", reset_source(i_source), &
          !             "ispec_source_first_guess:",ispec_source_first_guess
          !     print *
          !   endif
          !   call synchronize_all()
          ! endif

          ispec_source_first_guess = ispec_first_guess_vec(i_source)

          !  if (it_l > 51) then
          !    call synchronize_all()
          !    if (myrank == 0) then
          !      print *
          !      print *,myrank," ispec_source_first_guess :", ispec_source_first_guess, "xmin:", &
          !              coord(1,ibool(1,1,ispec_source_first_guess)), &
          !              "xmax:", coord(1,ibool(NGLLX,1,ispec_source_first_guess))
          !    endif
          !    call synchronize_all()
          !    if (myrank == 1) then
          !      print *,myrank," ispec_source_first_guess :", ispec_source_first_guess, "xmin:", &
          !              coord(1,ibool(1,1,ispec_source_first_guess)), &
          !              "xmax:", coord(1,ibool(NGLLX,1,ispec_source_first_guess))
          !    endif
          !    call synchronize_all()
          !    if (myrank == 2) then
          !      print *,myrank," ispec_source_first_guess :", ispec_source_first_guess, "xmin:", &
          !              coord(1,ibool(1,1,ispec_source_first_guess)), &
          !              "xmax:", coord(1,ibool(NGLLX,1,ispec_source_first_guess))
          !    endif
          !    call synchronize_all()
          !    if (myrank == 3) then
          !      print *,myrank," ispec_source_first_guess :", ispec_source_first_guess, "xmin:", &
          !              coord(1,ibool(1,1,ispec_source_first_guess)), &
          !              "xmax:", coord(1,ibool(NGLLX,1,ispec_source_first_guess))
          !      print *
          !    endif
          !    call synchronize_all()
          ! endif

          if (ispec_source_first_guess /= 0) then
            src_in_first_guess_element = in_element(xsrc, zsrc, ispec_source_first_guess)

            if ((.not. src_in_first_guess_element) .and. source_belonged_to_this_rank_all(i_source)) then
              ! All the procs have to recompute the source position
              reset_source(i_source) = .true.
              ! Added this also to check if the source is in adjacent elements
              ! (in this case we do not need to reset either because these are scanned later)
              if (pt_first_guess(i_source)%allocated) then
                do i_adjacent = 1,pt_first_guess(i_source)%size
                  ispec_adjacent = pt_first_guess(i_source)%array(i_adjacent)
                  src_in_adjacent_element = in_element(xsrc, zsrc, ispec_adjacent)
                  if (src_in_adjacent_element) then
                    reset_source(i_source) = .false.
                    exit
                  endif
                enddo
              else
                call stop_the_code('pt_first_guess should have been allocated if ispec_source_first_guess /= 0')
              endif
            else
              ! The source is in first guess element, no need to scan the others
              reset_source(i_source) = .false.
            endif
          endif
        enddo

        ! synchronizes between slices
        ! (moves synchronization out of do-loop above to avoid dead-locks)
#ifdef WITH_MPI
        if (NPROC > 1) then
          do i_source = 1,NSOURCES
            call any_all_l(reset_source(i_source), any_reset_source(i_source))
          enddo
          call synchronize_all()
        endif
#endif

        ! determines best source positions
        do i_source = 1,NSOURCES
          ! time
          t_used = (time_val - t0 - tshift_src(i_source))
          ! Moves and re-locates the sources along x and z-axis
          xsrc = x_source(i_source) + vx_source(i_source) * t_used
          zsrc = z_source(i_source) + vz_source(i_source) * t_used

          ispec_source_first_guess = ispec_first_guess_vec(i_source)

          !debug
          ! if (it_l > 51) then
          !   call synchronize_all()
          !   if (myrank == 0) then
          !     print *,myrank,"reset :", reset_source(i_source)," any_reset :", any_reset_source(i_source), &
          !            " src_in_first_guess_element :", src_in_first_guess_element
          !   endif
          !   call synchronize_all()
          !   if (myrank == 1) then
          !     print *,myrank,"reset :", reset_source(i_source)," any_reset :", any_reset_source(i_source), &
          !            " src_in_first_guess_element :", src_in_first_guess_element
          !   endif
          !   call synchronize_all()
          !   if (myrank == 2) then
          !     print *,myrank,"reset :", reset_source(i_source)," any_reset :", any_reset_source(i_source), &
          !            " src_in_first_guess_element :", src_in_first_guess_element
          !   endif
          !   call synchronize_all()
          !   if (myrank == 3) then
          !     print *,myrank,"reset :", reset_source(i_source)," any_reset :", any_reset_source(i_source), &
          !            " src_in_first_guess_element :", src_in_first_guess_element
          !     call sleep(1)
          !   endif
          !   call synchronize_all()
          !
          !   if (myrank == 0) then
          !     print *
          !     print *,myrank,"NOW LOCATING THE SOURCE ? ", &
          !             (source_belonged_to_this_rank_all(i_source) .or. any_reset_source(i_source))
          !   endif
          !   call synchronize_all()
          !   if (myrank == 1) then
          !     print *,myrank,"NOW LOCATING THE SOURCE ? ", &
          !             (source_belonged_to_this_rank_all(i_source) .or. any_reset_source(i_source))
          !   endif
          !   call synchronize_all()
          !   if (myrank == 2) then
          !     print *,myrank,"NOW LOCATING THE SOURCE ? ", &
          !             (source_belonged_to_this_rank_all(i_source) .or. any_reset_source(i_source))
          !   endif
          !   call synchronize_all()
          !   if (myrank == 3) then
          !     print *,myrank,"NOW LOCATING THE SOURCE ? ", &
          !             (source_belonged_to_this_rank_all(i_source) .or. any_reset_source(i_source))
          !   endif
          !   call synchronize_all()
          ! endif

          ! This is done only by the slice carrying the source (without the expensive loops)
          ! Except if a reset has been sent by any of the procs
          if (source_belonged_to_this_rank_all(i_source) .or. any_reset_source(i_source)) then
            ! initializes
            ispec_selected_source_local = 1
            ! locate force source (most expensive step)
            call locate_source_moving(xsrc,zsrc, &
                                      ispec_selected_source_local,islice_selected_source(i_source), &
                                      NPROC,myrank,xi_source(i_source),gamma_source(i_source),is_force_source, &
                                      source_belonged_to_this_rank=source_belonged_to_this_rank_all(i_source), &
                                      ispec_first_guess=ispec_source_first_guess, &
                                      pt_first_guess=pt_first_guess(i_source), &
                                      reset=any_reset_source(i_source))

            ispec_first_guess_vec(i_source) = ispec_selected_source_local

            if (myrank == islice_selected_source(i_source)) then
              ! Compute the sourcearray for this source and this time step
              ! (for interpolation: the source does not fall exactly at a GLL point)
              ! based on the source position in the reference element (xi_source, gamma_source)
              ! do some checking also
              call setup_source_interpolation_moving(xi_source(i_source), gamma_source(i_source), &
                                                     sourcearray, &
                                                     ispec_selected_source_local, &
                                                     source_type(i_source))

              ! counts sources in this process slice for this time step
              ! (will be used to fill sourcearrays_moving for local sources only)
              nsources_local = nsources_local + 1

              ! stores best position (non-zero) element ispec
              ispec_selected_source_moving(nsources_local, it_l) = ispec_selected_source_local
              ! stores sourcearray
              sourcearrays_moving(:,:,:,nsources_local,it_l) = sourcearray(:,:,:)
              ! STF
              source_time_function_moving(nsources_local,it_l) = source_time_function(i_source,it_l,i_stage_loc)

            endif  ! myrank == islice_selected_source(i_source)
          endif  ! source_belonged_to_this_rank .or. any_reset

          ! note: sourcearrays_moving(:,:,:,isource_local,it_l) is only non-zero for process
          !       which holds (best) source position for source i_source, sorted with local sources only.
          !       nsources_local will indicate that the process has local sources at time it.
          !
          !       when adding source contributions on the GPU, we assume that the arrays
          !       sourcearrays_moving(..) and ispec_selected_source_moving(..) are sorted with only local source entries.
          !       this avoids the need of storing islice_selected_source(NSOURCES,NSTEP) for all sources and time steps.
          !debug
          !if (.false.) then
          !  call synchronize_all()
          !  do iproc = 0,NPROC-1
          !    if (myrank == iproc .and. myrank == islice_selected_source(i_source)) then
          !      print *,it_l,"rank",myrank,"source ",i_source,"source local",nsources_local, &
          !              " ispec_selected_source_moving :", ispec_selected_source_moving(nsources_local, it_l), &
          !              " islice_selected_source(i_source) :",islice_selected_source(i_source), &
          !              "reset :", reset_source(i_source),"sourcearray",sourcearrays_moving(1,2,2,nsources_local,it_l)
          !    endif
          !  enddo
          !  call synchronize_all()
          !endif

        enddo  ! NSOURCES
      enddo ! NSTAGE_TIME_SCHEME

      ! check
      if (nsources_local > NSOURCES) stop 'Invalid nsources_local for moving sources'

      ! counts sources in this process slice for this time step
      nsources_local_moving(it_l) = nsources_local

    enddo ! NSTEP

    ! Write moving databases to a binary file
    if (write_moving_sources_database) then
      call write_moving_databases(pathToMovingDatabase)
    endif
  else
    ! Read the binary file file_exists .and. write_moving_sources_database
    call read_moving_databases(pathToMovingDatabase)
  endif

  ! Beware:
  ! When the source is moving we don't know where it is going: all the slices
  ! need to know the source_time_function
  ! If the source is not moving only the slice containing the source knows the source_time_function

  ! Send variables to device, once and for all:
  call prepare_moving_sources_cuda(Mesh_pointer, nsources_local_moving, NSOURCES, &
                                   sourcearrays_moving, ispec_selected_source_moving, &
                                   NSTEP, source_time_function_moving)

  ! Free memory
  do i_source = 1,NSOURCES
    if (pt_first_guess(i_source)%allocated) then
      deallocate(pt_first_guess(i_source)%array)
      pt_first_guess(i_source)%allocated = .false.
    else
      nullify(pt_first_guess(i_source)%array)
    endif
    pt_first_guess(i_source)%size = 0
  enddo

  deallocate(ispec_selected_source_moving)
  deallocate(sourcearrays_moving)
  deallocate(source_time_function_moving)

  ! synchronizes processes
  call synchronize_all()

end subroutine init_moving_sources_GPU

!
!-----------------------------------------------------------------------------------------
!

  subroutine locate_source_moving(x_source,z_source, &
                                  ispec_selected_source,islice_selected_source, &
                                  NPROC,myrank, xi_source,gamma_source, is_force_source, &
                                  source_belonged_to_this_rank, &
                                  ispec_first_guess, pt_first_guess, reset)

   !----
   !---- subroutine locate_source_moving below finds the correct position of
   !---- the (point force/momen-tensor) source in the case of moving sources.
   !---- It is really similar to subroutine locate_source but
   !---- in the case of moving sources this routine if a little bit more complicated
   !---- because we have to loop over it quite a lot which can be very costly if
   !---- we do it dumbly.
   !---- Hopefully, the source not moving too quickly, we know that it will stay in the same
   !---- element (or the ones around) for some time so we can avoid looping over
   !---- all the elements of the mesh too often
   !----

  ! Inputs
  double precision, intent(in) :: x_source, z_source
  integer, intent(in)  :: NPROC, myrank
  logical, intent(in) :: is_force_source
  integer, intent(in), optional :: ispec_first_guess
  logical, intent(in), optional :: reset

  ! source information, outputs
  integer, intent(inout) :: ispec_selected_source, islice_selected_source
  ! We use a custom pointer here because this can change size
  type(custom_array_pointer_type),intent(inout), optional :: pt_first_guess
  logical, intent(inout), optional :: source_belonged_to_this_rank
  double precision, intent(inout) :: xi_source, gamma_source

  ! local variables
  !! DK DK dec 2017: also loop on all the elements in contact with the initial guess
  !! element to improve accuracy of estimate
  logical, dimension(nglob) :: flag_topological
  integer :: number_of_mesh_elements_for_the_initial_guess
  integer, dimension(:), allocatable :: array_of_all_elements_of_ispec_selected_source

  ! other local parameters
  integer :: i, j, ispec, iglob, iter_loop, ix_initial_guess, iz_initial_guess, number_of_iterations
  integer :: imin, imax, jmin, jmax, ispecmin, ispecmax
  integer :: idomain, iglob_source, ispec_guess, ispec_test
  integer :: is_proc_source
  integer, dimension(1:NPROC) :: allgather_is_proc_source
  integer, dimension(1) :: locate_is_proc_source
  double precision :: dist_squared
  double precision :: xi, gamma, dx, dz, dxi, dgamma
  double precision :: x, z, xix, xiz, gammax, gammaz, jacobian
  double precision :: distmin_squared, dist_glob_squared
  double precision :: final_distance, final_distance_this_element
  logical :: do_next_loops, return_array_of_adj_elem, point_in_first_guess_element, rst
  logical :: ldummy

  ! check
  if (.not. SOURCE_IS_MOVING) then
    call exit_MPI(myrank,'subroutine locate_source_moving should be used when the source is moving')
  endif

  rst = .false.
  if (present(reset)) then ! We define another variable to avoid repeating if (present) ...
    rst = reset
  endif

  ! initializes slice holding the source
  is_proc_source = 0
  ! initialize closest global point
  iglob_source = 0
  ! initialize the first guess element for the source
  ispec_guess = -1
  ! initialize ispec for the loop
  ispec_test = -1

  ldummy = is_force_source ! unused yet, dummy to avoid compiler warning

  ! These have to be initialized otherwise the compiler can complain
  ! In the case without MPI (inside preprocessor directive)
  allgather_is_proc_source(:) = 0
  dist_glob_squared = 0.0d0
  locate_is_proc_source(:) = 0

  ! set distance to huge initial value
  distmin_squared = HUGEVAL

  ! search range in element (not used here, but I let it here if someday we need
  ! to get back to a scan in only the inner elements)
  imin = 1
  imax = NGLLX
  jmin = 1
  jmax = NGLLZ

  ! The second expensive loops to do (we will try to avoid doing it too much)
  do_next_loops = .true.

  ! Do we return the an array containing adjacent element indices ? Then we can
  ! look for the source only in these
  return_array_of_adj_elem = .false.

  ! Search range in the mesh (we will try to reduce it if possible!)
  ispecmin = 1
  ispecmax = nspec

  ! In which condition do we skip the second loop or reduce the search range ?
  if (present(ispec_first_guess)) then  ! If a first guess has been supplied
    ispec_guess = ispec_first_guess  ! To avoid having to repeat : if (present(..))...
    if (present(pt_first_guess)) then   ! If an array containing adjacent elements has been supplied
      return_array_of_adj_elem = .true.
      if (pt_first_guess%allocated) then
        do_next_loops = .false.
      endif
    else
      call stop_the_code('Need pt_first_guess with ispec_first_guess (even if not allocated)')
    endif

    if (.not. present(source_belonged_to_this_rank)) &
      call stop_the_code('Need source_belonged_to_this_rank with ispec_first_guess')

    if (ispec_guess /= 0) then  ! If a first guess has been supplied
      ! Is the source in this first guess element ?
      point_in_first_guess_element = in_element(x_source, z_source, ispec_guess)
      if (point_in_first_guess_element) then ! If yes, no need to do the first loop
        ispecmin = ispec_guess  ! The first loop will be a loop over one element
        ispecmax = ispec_guess  ! The first loop will be a loop over one element
      else ! If the source is not in the first guess element we have to perform the loops
        do_next_loops = .true.
      endif
    endif
  endif

  if (rst) then  ! If reset has been given all the loops are performed
    do_next_loops = .true.
    ispecmin = 1
    ispecmax = nspec
  endif

  ! loops over all elements, except if a first guess is given:
  ! in this case we only check that one
  do ispec = ispecmin, ispecmax
    do j = jmin,jmax
      do i = imin,imax
        iglob = ibool(i,j,ispec)

        ! we compare squared distances instead of distances themselves to significantly speed up calculations
        dist_squared = (x_source-dble(coord(1,iglob)))**2 + (z_source-dble(coord(2,iglob)))**2

        ! keep this point if it is closer to the source
        if (dist_squared < distmin_squared) then
          iglob_source = iglob
          distmin_squared = dist_squared
          ispec_test = ispec
          ix_initial_guess = i
          iz_initial_guess = j
          ! determines domain for outputting element type
          ! not used here but could be in the future
          if (ispec_is_acoustic(ispec)) then
            idomain = IDOMAIN_ACOUSTIC
          else if (ispec_is_elastic(ispec)) then
            idomain = IDOMAIN_ELASTIC
          else if (ispec_is_poroelastic(ispec)) then
            idomain = IDOMAIN_POROELASTIC
          else
            call stop_the_code('Invalid element type in locating source found!')
          endif
        endif
      enddo
    enddo
  enddo

  if (do_next_loops) then
    !! DK DK dec 2017: also loop on all the elements in contact with the initial guess element to improve accuracy of estimate
    flag_topological(:) = .false.

    ! mark the four corners of the initial guess element
    flag_topological(ibool(1,1,ispec_test)) = .true.
    flag_topological(ibool(NGLLX,1,ispec_test)) = .true.
    flag_topological(ibool(NGLLX,NGLLZ,ispec_test)) = .true.
    flag_topological(ibool(1,NGLLZ,ispec_test)) = .true.

    ! loop on all the elements to count how many are shared with the initial guess
    number_of_mesh_elements_for_the_initial_guess = 1
    do ispec = 1,nspec
      if (ispec == ispec_test) cycle  ! skip the initial guess element
      ! loop on the four corners only, no need to loop on the rest since we just want to detect adjacency
      do j = 1,NGLLZ,NGLLZ-1
        do i = 1,NGLLX,NGLLX-1
          if (flag_topological(ibool(i,j,ispec))) then
            ! this element is in contact with the initial guess
            number_of_mesh_elements_for_the_initial_guess = number_of_mesh_elements_for_the_initial_guess + 1
            ! let us not count it more than once, it may have a full edge in contact with it and would then be counted twice
            goto 700  ! Exit the loop over the four corners: we know this element is adjacent
          endif
        enddo
      enddo
      700 continue
    enddo

    ! now that we know the number of adjacent elements, we can allocate the list of elements and create it
    allocate(array_of_all_elements_of_ispec_selected_source(number_of_mesh_elements_for_the_initial_guess))

    ! first store the initial guess itself
    number_of_mesh_elements_for_the_initial_guess = 1
    array_of_all_elements_of_ispec_selected_source(number_of_mesh_elements_for_the_initial_guess) = ispec_test

    ! then store all the others
    do ispec = 1,nspec
      if (ispec == ispec_test) cycle
      ! loop on the four corners only, no need to loop on the rest since we just want to detect adjacency
      do j = 1,NGLLZ,NGLLZ-1
        do i = 1,NGLLX,NGLLX-1
          if (flag_topological(ibool(i,j,ispec))) then
            ! this element is in contact with the initial guess
            number_of_mesh_elements_for_the_initial_guess = number_of_mesh_elements_for_the_initial_guess + 1
            array_of_all_elements_of_ispec_selected_source(number_of_mesh_elements_for_the_initial_guess) = ispec
            ! let us not count it more than once, it may have a full edge in contact with it and would then be counted twice
            goto 800
          endif
        enddo
      enddo
      800 continue
    enddo
    ! The following returns array_of_all_elements_of_ispec_selected_source
    ! so it can be used in future calls (next time iterations) without having
    ! to recompute the previous loops over all the elements
    if (return_array_of_adj_elem) then
      ! deallocate pt_first_guess
      if (pt_first_guess%allocated) then
        deallocate(pt_first_guess%array)
        pt_first_guess%allocated = .false.
        pt_first_guess%size = 0
      else
        nullify(pt_first_guess%array)
      endif
      ! Reallocate with new values
      pt_first_guess%size = number_of_mesh_elements_for_the_initial_guess
      allocate(pt_first_guess%array(number_of_mesh_elements_for_the_initial_guess))
      pt_first_guess%allocated = .true.
      pt_first_guess%array = array_of_all_elements_of_ispec_selected_source
    endif
  else  ! If the loop is not done, we use the array supplied
    if (present(pt_first_guess)) then
      number_of_mesh_elements_for_the_initial_guess = pt_first_guess%size
      allocate(array_of_all_elements_of_ispec_selected_source(number_of_mesh_elements_for_the_initial_guess))
      array_of_all_elements_of_ispec_selected_source = pt_first_guess%array
    else
      call stop_the_code('do_next_loops == .false. needs pt_first_guess!')
    endif
  endif

  ! ************************************************************************************
  ! find the best (xi,gamma) for each source
  ! Given the coordinate (x,y,z) of the source that we know is in a spectral element
  ! How do we get the best corresponding coordinates in the reference element [-1,1]^2 ?
  ! We apply a small non linear algorithm, at each iteration we compute the distance
  ! between the point proposed and the source coordinate and we try to minimize this
  ! following the gradient
  ! ************************************************************************************

  final_distance = HUGEVAL

  do i = 1,number_of_mesh_elements_for_the_initial_guess ! Loop on the neighbors elements to ispec_selected_source

    !! DK DK dec 2017 set initial guess in the middle of the element,
    !! DK DK dec 2017 since we computed the true one only for the true initial guess
    !! DK DK dec 2017 the nonlinear process below will converge anyway
    if (i > 1) then
      ix_initial_guess = int(NGLLX / 2.0)
      iz_initial_guess = int(NGLLZ / 2.0)
    endif

    ispec = array_of_all_elements_of_ispec_selected_source(i)

    ! use initial guess in xi and gamma
    xi = xigll(ix_initial_guess)  ! Even with AXISYM this is OK! this is only an initial guess
                                  ! The algorithm will converge anyway
    gamma = zigll(iz_initial_guess)

    ! iterate to solve the nonlinear system
    if (USE_BEST_LOCATION_FOR_SOURCE) then
      number_of_iterations = NUM_ITER
    else
      number_of_iterations = 0 ! this means that the loop below will not be executed, i.e. we will not iterate
    endif

    do iter_loop = 1,number_of_iterations

      ! recompute jacobian for the new point
      call recompute_jacobian_with_negative_stop(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                                                 coorg,knods,ispec,NGNOD,nspec,npgeo,.true.)

      ! compute distance to target location
      dx = - (x - x_source)
      dz = - (z - z_source)

      ! compute increments (in direction of the distance gradient)
      dxi  = xix*dx + xiz*dz
      dgamma = gammax*dx + gammaz*dz

      ! update values
      xi = xi + dxi
      gamma = gamma + dgamma

      ! impose that we stay in that element
      ! (useful if user gives a source outside the mesh for instance)
      ! here we can't go outside the [1,1] segment even if with finite elements
      ! the polynomial solution is defined everywhere because then it becomes
      ! a mess to determine if the point considered is inside or outside the
      ! first guess element.
      ! Maybe this can poses a problem for convergence of itertive scheme
      ! with distorted elements
      if (xi > 1.00d0) xi = 1.00d0   ! if (xi > 1.01d0) xi = 1.01d0
      if (xi < -1.00d0) xi = -1.00d0
      if (gamma > 1.00d0) gamma = 1.00d0
      if (gamma < -1.00d0) gamma = -1.00d0

    ! end of nonlinear iterations
    enddo

    ! compute final coordinates of point found
    call recompute_jacobian_with_negative_stop(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                                               coorg,knods,ispec,NGNOD,nspec,npgeo,.true.)

    ! compute final distance between asked and found (using square value since taking the root is expensive)
    final_distance_this_element = (x_source-x)**2 + (z_source-z)**2

    ! if we have found an element that gives a shorter distance
    if (final_distance_this_element < final_distance) then
      ! store element number found
      ispec_selected_source = ispec
      ! store xi,gamma of point found
      xi_source = xi
      gamma_source = gamma

      ! store final distance between asked and found
      final_distance = final_distance_this_element
    endif
  enddo

#ifdef WITH_MPI
  ! for MPI version, gather information from all the nodes
  if (NPROC > 1) then
    ! MPI commands are necessary! The processes does not know the same part of the
    ! mesh. But we do it only when necessary to avoid exchanging too much in the
    ! case of moving sources
    if (rst .or. (ispec_guess == 0)) then

      ! global minimum distance computed over all processes
      call min_all_all_dp(final_distance, dist_glob_squared)

      ! check if this process contains the source
      if (abs(dist_glob_squared - final_distance) < TINYVAL ) is_proc_source = 1

      if (present(source_belonged_to_this_rank)) then
        if (is_proc_source == 1) then ! ADDED TODO
          source_belonged_to_this_rank = .true.
        else
          source_belonged_to_this_rank = .false.
        endif
      endif

      ! main collects info
      call gather_all_singlei(is_proc_source,allgather_is_proc_source,NPROC)

      if (myrank == 0) then
        ! select slice with maximum rank which contains source
        locate_is_proc_source = maxloc(allgather_is_proc_source) - 1
        islice_selected_source = locate_is_proc_source(1)
      endif

      ! selects slice which holds source
      call bcast_all_singlei(islice_selected_source)

      ! if (islice_selected_source /= 0) then
      !   ! source is in another slice than the main process
      !   if (myrank == islice_selected_source) then
      !     ! send information from process holding source
      !     call send_singlei(ispec_selected_source,0,0)
      !     call send_singledp(xi_source,0,2)
      !     call send_singledp(gamma_source,0,3)
      !   else if (myrank == 0) then
      !     ! main collects
      !     call recv_singlei(ispec_selected_source,islice_selected_source,0)
      !     call recv_singledp(xi_source,islice_selected_source,2)
      !     call recv_singledp(gamma_source,islice_selected_source,3)
      !   endif
      ! endif
    endif
  else
    ! single MPI process
    islice_selected_source = 0
  endif
#else
  ! no MPI, single process
  islice_selected_source = 0
#endif

  ! user output
  if (myrank == 0) then
    ! check position
    if (final_distance == HUGEVAL) call exit_MPI(myrank,'Error locating source')
  endif

  deallocate(array_of_all_elements_of_ispec_selected_source)

end subroutine locate_source_moving

!
!-----------------------------------------------------------------------------------------
!

subroutine write_moving_databases(pathToMovingDatabase)
  !---
  !--- Write databases for moving source position to binary file so it can be read quickly
  !---

  use constants, only: MAX_STRING_LEN, IMAIN

  use specfem_par, only: myrank, NSOURCES, NSTEP, DT, nglob, nglob_elastic, nglob_acoustic, &
                         x_source, z_source, vx_source, vz_source

  implicit none

  character(len=MAX_STRING_LEN), intent(in) :: pathToMovingDatabase

  ! Local variables
  integer :: ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Moving source databases will be written in',pathToMovingDatabase,' ...'
    call flush_IMAIN()
  endif

  open(unit = IMOV, file = pathToMovingDatabase,status='unknown',action='write',form='unformatted',iostat=ier)
  if (ier /= 0) call stop_the_code('  Error writing moving source data file to disk !')

  write(IMOV) NSTEP
  write(IMOV) DT
  write(IMOV) nglob
  write(IMOV) nglob_elastic
  write(IMOV) nglob_acoustic
  write(IMOV) NSOURCES
  write(IMOV) x_source
  write(IMOV) z_source
  write(IMOV) vx_source
  write(IMOV) vz_source
  write(IMOV) nsources_local_moving
  write(IMOV) sourcearrays_moving
  write(IMOV) ispec_selected_source_moving
  write(IMOV) source_time_function_moving

  close(IMOV)

  ! user output
  write(IMAIN,*) '  Moving source databases have been written in ', pathToMovingDatabase
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  Just rerun the code now, they will be read there'
    call flush_IMAIN()
  endif
  call exit_MPI(myrank, '  Terminating ...')

end subroutine write_moving_databases

!
!-----------------------------------------------------------------------------------------
!

subroutine read_moving_databases(pathToMovingDatabase)
  !---
  !--- Read databases for moving source position from file, make some checking
  !---

  use constants, only: MAX_STRING_LEN, IMAIN

  use specfem_par, only: myrank, NSOURCES, NSTEP, DT, nglob, nglob_elastic, nglob_acoustic, &
                         x_source, z_source, vx_source, vz_source

  implicit none

  character(len=MAX_STRING_LEN), intent(in) :: pathToMovingDatabase

  ! Local variables
  integer :: ier, i_source
  integer :: NSTEP_read, nglob_read, nglob_elastic_read, nglob_acoustic_read, NSOURCES_read
  double precision :: DT_read
  double precision, dimension(:), allocatable :: x_source_read, z_source_read
  double precision, dimension(:), allocatable :: vx_source_read, vz_source_read
  character(len=MAX_STRING_LEN) :: str

  ! Read data from a binary file
  write(IMAIN,*) '  Reading moving source databases from file: ', pathToMovingDatabase
  open(unit = IMOV,file = pathToMovingDatabase,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_data.bin')

  read(IMOV) NSTEP_read
  read(IMOV) DT_read
  read(IMOV) nglob_read
  read(IMOV) nglob_elastic_read
  read(IMOV) nglob_acoustic_read
  read(IMOV) NSOURCES_read
  str = ''
  if (NSTEP_read /= NSTEP) then
    write(str,*) 'Error: Database read does not match the Par_file (different NSTEP:', &
                 NSTEP,'instead of',NSTEP_read,')'
    call stop_the_code(NEW_LINE('') // trim(str) // NEW_LINE('') // '  - > Terminating' // NEW_LINE(''))
  endif
  if (abs(DT_read - DT) > TINYVAL) then
    write(str,*) 'Error: Database read does not match the Par_file (different DT:', &
                  DT,'instead of',DT_read,')'
    call stop_the_code(NEW_LINE('') // trim(str) // NEW_LINE('') // '  - > Terminating' // NEW_LINE(''))
  endif
  if (nglob_read /= nglob) then
    write(str,*) 'Error: Database read does not match the current configuration', &
                   '(different number of GLL points:',nglob,'instead of',nglob_read,')'
    call stop_the_code(NEW_LINE('') // trim(str) // NEW_LINE('') // '  - > Terminating' // NEW_LINE(''))
  endif
  if (nglob_elastic_read /= nglob_elastic) then
    write(str,*) 'Error: Database read does not match the current configuration (different number of', &
                   'elastic elements',nglob_elastic,'instead of',nglob_elastic_read,')'
    call stop_the_code(NEW_LINE('') // trim(str) // NEW_LINE('') // '  - > Terminating' // NEW_LINE(''))
  endif
  if (nglob_acoustic_read /= nglob_acoustic) then
    write(str,*) 'Error: Database read does not match the current configuration (different number of', &
                   'acoustic elements',nglob_acoustic,'instead of',nglob_acoustic_read,')'
    call stop_the_code(NEW_LINE('') // trim(str) // NEW_LINE('') // '  - > Terminating' // NEW_LINE(''))
  endif
  if (NSOURCES_read /= NSOURCES) then
    write(str,*) 'Error: Database read does not match the SOURCES file (different number of sources:', &
                  NSOURCES,'instead of',NSOURCES_read,')'
    call stop_the_code(NEW_LINE('') // trim(str) // NEW_LINE('') // '  - > Terminating' // NEW_LINE(''))
  endif
  allocate(x_source_read(NSOURCES), z_source_read(NSOURCES), vx_source_read(NSOURCES), &
           vz_source_read(NSOURCES))
  read(IMOV) x_source_read
  read(IMOV) z_source_read
  read(IMOV) vx_source_read
  read(IMOV) vz_source_read
  read(IMOV) nsources_local_moving
  do i_source = 1,NSOURCES
    if (abs(x_source_read(i_source) - x_source(i_source)) > TINYVAL) then
      write(str,*) 'Error: Database read does not match the current SOURCES file', &
                     '(different sources x coord:',x_source,'instead of',x_source_read,')'
      call stop_the_code(trim(str) // '  Terminating')
    endif
    if (abs(z_source_read(i_source) - z_source(i_source)) > TINYVAL) then
      write(str,*) 'Error: Database read does not match the current SOURCES file', &
                     '(different sources z coord:',z_source,'instead of',z_source_read,')'
      call stop_the_code(NEW_LINE('') // trim(str) // NEW_LINE('') // '  - > Terminating' // NEW_LINE(''))
    endif
    if (abs(vx_source_read(i_source) - vx_source(i_source)) > TINYVAL) then
      write(str,*) 'Error: Database read does not match the current SOURCES file', &
                     '(different sources x velocity:',vx_source,'instead of',vx_source_read,')'
      call stop_the_code(NEW_LINE('') // trim(str) // NEW_LINE('') // '  - > Terminating' // NEW_LINE(''))
    endif
    if (abs(vz_source_read(i_source) - vz_source(i_source)) > TINYVAL) then
      write(str,*) 'Error: Database read does not match the current SOURCES file', &
                     '(different sources z velocity:',vz_source,'instead of',vz_source_read,')'
      call stop_the_code(NEW_LINE('') // trim(str) // NEW_LINE('') // '  - > Terminating' // NEW_LINE(''))
    endif
  enddo

  ! re-allocates arrays
  deallocate(ispec_selected_source_moving)
  deallocate(sourcearrays_moving)
  deallocate(source_time_function_moving)

  allocate(ispec_selected_source_moving(NSOURCES, NSTEP), &
           source_time_function_moving(NSOURCES,NSTEP), &
           sourcearrays_moving(NDIM,NGLLX,NGLLX,NSOURCES,NSTEP),stat=ier)
  if (ier /= 0) stop 'Error allocating moving sources arrays'

  ! reads in arrays from database
  read(IMOV) sourcearrays_moving
  read(IMOV) ispec_selected_source_moving
  read(IMOV) source_time_function_moving
  close(IMOV)

end subroutine read_moving_databases

!
!-----------------------------------------------------------------------------------------
!

subroutine setup_source_interpolation_moving(xi_source, gamma_source, sourcearray, &
                                             ispec_source, this_source_type)
  !----
  !---- Compute sourcearray (for interpolation: the source does not fall at a GLL point)
  !---- based on the source position in the reference element (xi_source, gamma_source)
  !---- Do some checking also
  !----

  use constants, only: NDIM,NGLLX,NGLLZ,ZERO,CUSTOM_REAL

  use specfem_par, only: myrank, &
    xigll, zigll, hxis, hpxis, hgammas, hpgammas

  implicit none

  ! single source array
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ), intent(out) :: sourcearray
  double precision, intent(in) :: xi_source, gamma_source
  ! element containing the source
  integer, intent(in) :: ispec_source, this_source_type

  ! local parameters
  integer :: i,j
  double precision :: hlagrange

  ! (re)initializes (it has been done already in setup_sources_receivers.f90)
  hxis(:) = ZERO
  hgammas(:) = ZERO
  hpxis(:) = ZERO
  hpgammas(:) = ZERO
  sourcearray(:,:,:) = 0._CUSTOM_REAL

  ! Compute Lagrange interpolators for the source considered
  call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
  call lagrange_any(gamma_source,NGLLZ,zigll,hgammas,hpgammas)

  ! computes source arrays
  if (this_source_type == 1) then
    ! collocated force source
    do j = 1,NGLLZ
      do i = 1,NGLLX
        hlagrange = hxis(i) * hgammas(j)
        ! source element is acoustic
        if (ispec_is_acoustic(ispec_source)) then
          ! In the acoustic case both dimensions are equal (pressure)
          sourcearray(:,i,j) = hlagrange
        else
          call exit_MPI(myrank,'Moving source not in acoustic element (GPU, not implemented)')
        endif
      enddo
    enddo
  else
    call exit_MPI(myrank,'Moving source with source_type != 1 (GPU, not implemented)')
  endif

end subroutine setup_source_interpolation_moving

!
!-----------------------------------------------------------------------------------------
!

logical function in_element(x_coord, z_coord, ispec, debug)
  !----
  !---- This function in_element returns true if (x_coord, z_coord) belongs to element ispec
  !---- geometrically (taking the four corners into account)
  !----

  use specfem_par, only: ibool, coord, nspec

  implicit none

  integer,intent(in) :: ispec
  double precision,intent(in) :: x_coord, z_coord

  ! local variables
  integer :: iglob_1, iglob_2, iglob_3, iglob_4
  real (kind = 8) :: x1, x2, x3, x4, z1, z2, z3, z4
  integer :: in_or_out
  real (kind = 8), dimension(4) :: x, z
  logical, optional :: debug

  if ((ispec > nspec) .or. (ispec <= 0)) then
    ! Deal with pathologic cases to avoid segmentation faults
    in_element = .false.
    return
  endif

  iglob_1 = ibool(1,1,ispec)
  iglob_2 = ibool(NGLLX,1,ispec)
  iglob_3 = ibool(NGLLX,NGLLZ,ispec)
  iglob_4 = ibool(1,NGLLZ,ispec)
  x1 = real(coord(1,iglob_1), 8)
  z1 = real(coord(2,iglob_1), 8)
  x2 = real(coord(1,iglob_2), 8)
  z2 = real(coord(2,iglob_2), 8)
  x3 = real(coord(1,iglob_3), 8)
  z3 = real(coord(2,iglob_3), 8)
  x4 = real(coord(1,iglob_4), 8)
  z4 = real(coord(2,iglob_4), 8)
  x = (/ x1, x2, x3, x4 /)
  z = (/ z1, z2, z3, z4 /)

  if (present(debug)) then
    if (debug) then
      print *, x_coord, z_coord, " in "
      print *, "  ", x1, ",", z1
      print *, "  ", x2, ",", z2
      print *, "  ", x3, ",", z3
      print *, "  ", x4, ",", z4
    endif
  endif

  call pnpoly(4, x, z, x_coord, z_coord, in_or_out)
  in_element = (in_or_out == 0) .or. (in_or_out == 1)

  ! Other option:
  ! in_element = point_in_polygon(4, x, z, x_coord, z_coord)
  ! But not possible to deal directly with the case where the point falls on the edges

  return

end function in_element

!
!-----------------------------------------------------------------------------------------
!

function point_in_polygon ( n, x, y, x0, y0 )
  !----
  !---- function point_in_polygon determines if a point is inside a polygon
  !----
  !----  Discussion:
  !----
  !----    If the points ( x(i), y(i) ) ( i = 1, 2, ..., n ) are,
  !----    in this cyclic order, the vertices of a simple closed polygon and
  !----    (x0,y0) is a point not on any side of the polygon, then the
  !----    procedure determines, by setting "point_in_polygon" to TRUE or FALSE,
  !----    whether (x0,y0) lies in the interior of the polygon.
  !----
  !----  Licensing:
  !----
  !----    This code is distributed under the GNU LGPL license.
  !----
  !----  Modified:
  !----
  !----    07 November 2016
  !----
  !----  Author:
  !----
  !----    John Burkardt
  !----
  !----  Reference:
  !----
  !----    Moshe Shimrat,
  !----    ACM Algorithm 112,
  !----    Position of Point Relative to Polygon,
  !----    communications of the ACM,
  !----    Volume 5, Number 8, page 434, August 1962.
  !----
  !----    Richard Hacker,
  !----    Certification of Algorithm 112,
  !----    communications of the ACM,
  !----    Volume 5, Number 12, page  606, December 1962.
  !----
  !----  Parameters:
  !----
  !----    Input, integer ( kind = 4 ) N, the number of nodes or vertices in
  !----    the polygon.  N must be at least 3.
  !----
  !----    Input, real ( kind = 8 ) V(2,N), the vertices of the polygon.
  !----
  !----    Input, real ( kind = 8 ) P(2), the coordinates of the point to be tested.
  !----
  !----    Output, logical ( kind = 4 ) INSIDE, is TRUE if the point is
  !----    inside the polygon.
  !----

    implicit none

    integer (kind = 4) :: n
    integer (kind = 4) :: i, ip1
    logical :: b, point_in_polygon
    real (kind = 8) :: t, x0, y0
    real (kind = 8), dimension(n) :: x, y

    b = .false.

    do i = 1, n
      ip1 = mod ( i, n ) + 1
      !
      !   if ( ( y(ip1) < y0 .and. y0 <= y(i)   ) .or. &
      !        ( y(i) < y0 .and. y0 <= y(ip1) ) ) then
      !
      if ( y(ip1) < y0 .eqv. y0 <= y(i) ) then
        t = x0 - x(i) - ( y0 - y(i) ) * ( x(ip1) - x(i) ) / ( y(ip1) - y(i) )
        print *,"t:",t
        if ( t <= 0.0D+00 ) then
          b = .not. b
        endif
      endif
    enddo

    point_in_polygon = b

    return
end function point_in_polygon

!
!-----------------------------------------------------------------------------------------
!

function eor_condition (ix, iy)
  !----
  !---- Small logical function used in subroutine pnpoly
  !----

    implicit none

    logical :: ix , iy
    logical :: eor_condition

    eor_condition = (ix .or. iy) .and. .not. (ix .and. iy)

    return

end function eor_condition


!
!-----------------------------------------------------------------------------------------
!

subroutine pnpoly(n, x, y, x0, y0, in_or_out)
    !----
    !----  Courtesy: Jay Sandhu
    !----               email: jsandhu@esri.com
    !----
    !----
    !----  Come from David H. Douglas, COLLECTED ALGORITHMS, Cambridge MA:
    !----  Harvard Laboratory for Computer Graphics, 1974
    !----
    !----        Purpose
    !----           to determine whether a point is inside a polygon
    !----
    !----        Usage
    !----           call pnpoly (x0, y0, x, y, n, inout )
    !----
    !----        Description of the parameters
    !----           x0      - x-coordinate of point in question.
    !----           y0      - y-coordinate of point in question.
    !----           x       - n long vector containing x-coordinates of
    !----                     vertices of polygon.
    !----           y       - n long vector containing y-coordinates of
    !----                     vertices of polygon.
    !----           n       - number of vertices in the polygon.
    !----           in_or_out   - the signal returned:
    !----                     -1 if the point is outside of the polygon,
    !----                      0 if the point is on an edge or at a vertex,
    !----                      1 if the point is inside of the polygon.
    !----
    !----        Remarks
    !----           the vertices may be listed in clockwise or anticlockwise
    !----           order.  for this subroutine a point is considered inside
    !----           the polygon if it is located in the enclosed area defined
    !----           by the line forming the polygon.
    !----           the input polygon may be a compound polygon consisting
    !----           of several separate subpolygons. if so, the first vertex
    !----           of each subpolygon must be repeated, and when calculating
    !----           n, these first vertices must be counted twice.
    !----           inout is the only parameter whose value is changed.
    !----           pnpoly can handle any number of vertices in the polygon.
    !----           written by randolph franklin, university of ottawa, 6/72.
    !----
    !----        subroutines and function subprograms required
    !----           none
    !----
    !----        Method
    !----           a vertical semi-infinite line is drawn up from the point
    !----           in question. if it crosses the polygon an odd number of
    !----           times, the point is inside the polygon.
    !----

    implicit none

    ! Inputs, outputs
    integer, intent(in) :: n
    real (kind = 8), dimension(n), intent(in) :: x, y
    real (kind = 8), intent(in) :: x0 , y0
    integer, intent(out) :: in_or_out

    ! Local variables
    integer :: i , j
    real (kind = 8) :: xi , yi , xj , yj
    logical :: ix , iy , jx , jy

    in_or_out = -1

    do i = 1 , n
       xi = x(i) - x0
       yi = y(i) - y0
       ! check whether the point in question is at this vertex.
       if ( xi == 0.0 .and. yi == 0.0 ) then
          in_or_out = 0
          return
       endif
       ! j is next vertex number of polygon.
       j = 1 + mod(i,n)
       xj = x(j) - x0
       yj = y(j) - y0
       ! is this line of 0 length ?
       if ( xi == xj .and. yi == yj ) goto 100
       ix = xi >= 0.0
       iy = yi >= 0.0
       jx = xj >= 0.0
       jy = yj >= 0.0
       ! check whether (x0,y0) is on vertical side of polygon.
       if ( xi == 0.0 .and. xj == 0.0 .and. eor_condition(iy,jy) ) then
         in_or_out = 0
         return
       endif
       ! check whether (x0,y0) is on horizontal side of polygon.
       if ( yi == 0.0 .and. yj == 0.0 .and. eor_condition(ix,jx) ) then
         in_or_out = 0
         return
       endif
       ! check whether both ends of this side are completely 1) to right
       ! of, 2) to left of, or 3) below (x0,y0).
       if (.not. ((iy .or. jy) .and. eor_condition(ix,jx)) ) goto 100
       ! does this side obviously cross line rising vertically from (x0,y0)
       if (.not. (iy .and. jy .and. eor_condition(ix,jx)) ) then
         if ( (yi*xj-xi*yj)/(xj-xi) < 0.0 ) then
           goto 100
         else if ( (yi*xj-xi*yj)/(xj-xi) == 0.0 ) then
           in_or_out = 0
           return
         else
           in_or_out = -in_or_out
         endif
       else
         in_or_out = -in_or_out
       endif
100   enddo

    continue
    end subroutine pnpoly

end module moving_sources_par
