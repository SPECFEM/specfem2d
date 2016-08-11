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
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
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

  subroutine setup_sources_receivers()

  implicit none

  ! locates sources and determines simulation start time t0
  call setup_sources()

  ! reads in stations file and locates receivers
  call setup_receivers()

  ! reads in adjoint sources
  call setup_adjoint_sources()

  ! tangential components
  call setup_source_receiver_tangentials()

  ! pre-compute lagrangians and sourcearrays for sources
  call setup_source_interpolation()

  ! pre-compute lagrangians for receivers
  call setup_receiver_interpolation()

  ! synchronizes processes
  call synchronize_all()

  end subroutine setup_sources_receivers

!
!----------------------------------------------------------------------------
!

  subroutine setup_sources()

  use constants, only: NGLLX,NGLLZ,NDIM,IMAIN,IIN,MAX_STRING_LEN

  use specfem_par, only: NSOURCES,initialfield,source_type, &
                         coord,ibool,nglob,nspec,nelem_acoustic_surface,acoustic_surface, &
                         ispec_is_elastic,ispec_is_poroelastic, &
                         x_source,z_source,ispec_selected_source, &
                         islice_selected_source, &
                         xigll,zigll,npgeo, &
                         NPROC,myrank,xi_source,gamma_source,coorg,knods,ngnod, &
                         iglob_source
  implicit none

  ! Local variables
  integer :: ispec_acoustic_surface
  integer :: ixmin, ixmax, izmin, izmax,i_source,ispec

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'sources:'
    call flush_IMAIN()
  endif

  ! locates sources
  do i_source = 1,NSOURCES

    if (source_type(i_source) == 1) then
      ! collocated force source
      call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                         x_source(i_source),z_source(i_source), &
                         ispec_selected_source(i_source),islice_selected_source(i_source), &
                         NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo, &
                         iglob_source(i_source),.true.) ! flag .true. indicates force source

      ! check
      if (myrank == islice_selected_source(i_source)) then
        ! checks that acoustic source is not exactly on the free surface because pressure is zero there
        do ispec_acoustic_surface = 1,nelem_acoustic_surface
          ispec = acoustic_surface(1,ispec_acoustic_surface)
          ixmin = acoustic_surface(2,ispec_acoustic_surface)
          ixmax = acoustic_surface(3,ispec_acoustic_surface)
          izmin = acoustic_surface(4,ispec_acoustic_surface)
          izmax = acoustic_surface(5,ispec_acoustic_surface)
          if (.not. ispec_is_elastic(ispec) .and. .not. ispec_is_poroelastic(ispec) .and. &
            ispec == ispec_selected_source(i_source)) then
            if ((izmin == 1 .and. izmax == 1 .and. ixmin == 1 .and. ixmax == NGLLX .and. &
                gamma_source(i_source) < -0.99d0) .or. &
                (izmin == NGLLZ .and. izmax == NGLLZ .and. ixmin == 1 .and. ixmax == NGLLX .and. &
                gamma_source(i_source) > 0.99d0) .or. &
                (izmin == 1 .and. izmax == NGLLZ .and. ixmin == 1 .and. ixmax == 1 .and. &
                xi_source(i_source) < -0.99d0) .or. &
                (izmin == 1 .and. izmax == NGLLZ .and. ixmin == NGLLX .and. ixmax == NGLLX .and. &
                xi_source(i_source) > 0.99d0) .or. &
                (izmin == 1 .and. izmax == 1 .and. ixmin == 1 .and. ixmax == 1 .and. &
                gamma_source(i_source) < -0.99d0 .and. xi_source(i_source) < -0.99d0) .or. &
                (izmin == 1 .and. izmax == 1 .and. ixmin == NGLLX .and. ixmax == NGLLX .and. &
                gamma_source(i_source) < -0.99d0 .and. xi_source(i_source) > 0.99d0) .or. &
                (izmin == NGLLZ .and. izmax == NGLLZ .and. ixmin == 1 .and. ixmax == 1 .and. &
                gamma_source(i_source) > 0.99d0 .and. xi_source(i_source) < -0.99d0) .or. &
                (izmin == NGLLZ .and. izmax == NGLLZ .and. ixmin == NGLLX .and. ixmax == NGLLX .and. &
                gamma_source(i_source) > 0.99d0 .and. xi_source(i_source) > 0.99d0)) then
              call exit_MPI(myrank,'an acoustic source cannot be located exactly '// &
                            'on the free surface because pressure is zero there')
            endif
          endif
        enddo

        ! daniel debug
        !open(unit=1234,file='OUTPUT_FILES/debug_source_contribution.txt',status='unknown')
      endif

    else if (source_type(i_source) == 2) then
      ! moment-tensor source
      ! note: iglob_source is not really needed for moment-tensor sources, but left as argument since it's already allocated
      !       and might help for future routines...
      call locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                         x_source(i_source),z_source(i_source), &
                         ispec_selected_source(i_source),islice_selected_source(i_source), &
                         NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo, &
                         iglob_source(i_source),.false.) ! flag .false. indicates moment-tensor source

    else if (.not. initialfield) then

      call exit_MPI(myrank,'incorrect source type')

    endif

  enddo ! do i_source= 1,NSOURCES

!! DK DK this below not supported in the case of MPI yet, we should do a MPI_GATHER() of the values
!! DK DK and use "if (myrank == islice_selected_rec(irec)) then" to display the right sources
!! DK DK and receivers carried by each mesh slice, and not fictitious values coming from other slices
#ifndef USE_MPI
  if (myrank == 0) then
     ! write actual source locations to file
     ! note that these may differ from input values, especially if source_surf = .true. in SOURCE
     ! note that the exact source locations are determined from (ispec,xi,gamma) values
     open(unit=14,file='OUTPUT_FILES/for_information_SOURCE_actually_used',status='unknown')
     do i_source= 1,NSOURCES
        write(14,*) x_source(i_source), z_source(i_source)
     enddo
     close(14)
  endif
#endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_sources

!
!----------------------------------------------------------------------------
!

  subroutine setup_receivers()

  use constants, only: NGLLX,NGLLZ,NDIM,IMAIN,IIN,MAX_STRING_LEN
#ifndef USE_MPI
  use constants, only: IOUT
#endif

  use specfem_par, only: coord,ibool,nglob,nspec, &
                         ispec_selected_rec, &
                         NPROC,myrank,coorg,knods,ngnod, &
                         xigll,zigll,npgeo, &
                         nrec,nrecloc,recloc,islice_selected_rec,st_xval,st_zval, &
                         xi_receiver,gamma_receiver,station_name,network_name, &
                         x_final_receiver,z_final_receiver, &
                         x_source,z_source

  implicit none

  ! Local variables
  integer :: irec,nrec_tot_found
  integer :: ier

  character(len=MAX_STRING_LEN) :: dummystring

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'receivers:'
    call flush_IMAIN()
  endif

  ! get number of stations from receiver file
  open(unit=IIN,file='DATA/STATIONS',status='old',action='read',iostat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error opening DATA/STATIONS file')
  nrec = 0
  do while(ier == 0)
    read(IIN,"(a)",iostat=ier) dummystring
    if (ier == 0) nrec = nrec + 1
  enddo
  close(IIN)

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of receivers = ',nrec
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (nrec < 1) call exit_MPI(myrank,'need at least one receiver')

  ! receiver information
  allocate(ispec_selected_rec(nrec), &
           st_xval(nrec), &
           st_zval(nrec), &
           xi_receiver(nrec), &
           gamma_receiver(nrec), &
           station_name(nrec), &
           network_name(nrec), &
           recloc(nrec), &
           islice_selected_rec(nrec), &
           x_final_receiver(nrec), &
           z_final_receiver(nrec),stat=ier)
  if (ier /= 0) stop 'Error allocating receiver arrays'

  ! locate receivers in the mesh
  call locate_receivers(ibool,coord,nspec,nglob,xigll,zigll, &
                        nrec,nrecloc,recloc,islice_selected_rec,NPROC,myrank, &
                        st_xval,st_zval,ispec_selected_rec, &
                        xi_receiver,gamma_receiver,station_name,network_name, &
                        x_source(1),z_source(1), &
                        coorg,knods,ngnod,npgeo, &
                        x_final_receiver,z_final_receiver)

!! DK DK this below not supported in the case of MPI yet, we should do a MPI_GATHER() of the values
!! DK DK and use "if (myrank == islice_selected_rec(irec)) then" to display the right sources
!! DK DK and receivers carried by each mesh slice, and not fictitious values coming from other slices
  irec = 0
#ifndef USE_MPI
  if (myrank == 0) then
     ! write out actual station locations (compare with STATIONS from meshfem2D)
     ! NOTE: this will be written out even if use_existing_STATIONS = .true.
     open(unit=IOUT,file='OUTPUT_FILES/for_information_STATIONS_actually_used',status='unknown')
     do irec = 1,nrec
        write(IOUT,"('S',i4.4,'    AA ',f20.7,1x,f20.7,'       0.0         0.0')") &
             irec,x_final_receiver(irec),z_final_receiver(irec)
     enddo
     close(IOUT)
  endif
#endif

  ! synchronizes all processes
  call synchronize_all()

  ! checks that the sum of the number of receivers in each slice is nrec
  call sum_all_i(nrecloc,nrec_tot_found)
  if (myrank == 0) then
    ! checks total
    if (nrec_tot_found /= nrec) then
      call exit_MPI(myrank,'problem when dispatching the receivers')
    endif

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'found a total of ',nrec_tot_found,' receivers, this is okay'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! synchronizes all processes
  call synchronize_all()

  ! checks if acoustic receiver is exactly on the free surface because pressure is zero there
  call setup_receivers_check_acoustic()

  end subroutine setup_receivers


!
!-----------------------------------------------------------------------------------------
!

  subroutine setup_receivers_check_acoustic()

! checks if acoustic receiver is exactly on the free surface because pressure is zero there

  use constants, only: NGLLX,NGLLZ,IMAIN

  use specfem_par, only: myrank,ispec_selected_rec,nrecloc,recloc, &
                         xi_receiver,gamma_receiver,seismotype, &
                         nelem_acoustic_surface,acoustic_surface,ispec_is_acoustic
  implicit none

  ! Local variables
  integer :: irec,irecloc
  integer :: ispec,ispec_acoustic_surface
  integer :: ixmin, ixmax, izmin, izmax

  ! check if acoustic receiver is exactly on the free surface because pressure is zero there
  do ispec_acoustic_surface = 1,nelem_acoustic_surface
    ispec = acoustic_surface(1,ispec_acoustic_surface)
    ixmin = acoustic_surface(2,ispec_acoustic_surface)
    ixmax = acoustic_surface(3,ispec_acoustic_surface)
    izmin = acoustic_surface(4,ispec_acoustic_surface)
    izmax = acoustic_surface(5,ispec_acoustic_surface)
    do irecloc = 1,nrecloc
      irec = recloc(irecloc)
      if (ispec_is_acoustic(ispec) .and. ispec == ispec_selected_rec(irec)) then
        if ((izmin == 1 .and. izmax == 1 .and. ixmin == 1 .and. ixmax == NGLLX .and. &
        gamma_receiver(irec) < -0.99d0) .or. &
        (izmin == NGLLZ .and. izmax == NGLLZ .and. ixmin == 1 .and. ixmax == NGLLX .and. &
        gamma_receiver(irec) > 0.99d0) .or. &
        (izmin == 1 .and. izmax == NGLLZ .and. ixmin == 1 .and. ixmax == 1 .and. &
        xi_receiver(irec) < -0.99d0) .or. &
        (izmin == 1 .and. izmax == NGLLZ .and. ixmin == NGLLX .and. ixmax == NGLLX .and. &
        xi_receiver(irec) > 0.99d0) .or. &
        (izmin == 1 .and. izmax == 1 .and. ixmin == 1 .and. ixmax == 1 .and. &
        gamma_receiver(irec) < -0.99d0 .and. xi_receiver(irec) < -0.99d0) .or. &
        (izmin == 1 .and. izmax == 1 .and. ixmin == NGLLX .and. ixmax == NGLLX .and. &
        gamma_receiver(irec) < -0.99d0 .and. xi_receiver(irec) > 0.99d0) .or. &
        (izmin == NGLLZ .and. izmax == NGLLZ .and. ixmin == 1 .and. ixmax == 1 .and. &
        gamma_receiver(irec) > 0.99d0 .and. xi_receiver(irec) < -0.99d0) .or. &
        (izmin == NGLLZ .and. izmax == NGLLZ .and. ixmin == NGLLX .and. ixmax == NGLLX .and. &
        gamma_receiver(irec) > 0.99d0 .and. xi_receiver(irec) > 0.99d0)) then
          ! checks
          if (seismotype == 4) then
            call exit_MPI(myrank,'an acoustic pressure receiver cannot be located exactly '// &
                            'on the free surface because pressure is zero there')
          else
            write(IMAIN,*) '**********************************************************************'
            write(IMAIN,*) '*** Warning: acoustic receiver located exactly on the free surface ***'
            write(IMAIN,*) '*** Warning: tangential component will be zero there               ***'
            write(IMAIN,*) '**********************************************************************'
            write(IMAIN,*)
          endif
        endif
      endif
    enddo
  enddo

  end subroutine setup_receivers_check_acoustic


!
!-----------------------------------------------------------------------------------------
!


  subroutine setup_adjoint_sources

! compute source array for adjoint source

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,MAX_STRING_LEN,IMAIN

  use specfem_par, only: nadj_rec_local,nrec,nrecloc,NSTEP,NPROC,SIMULATION_TYPE,SU_FORMAT, &
                        adj_sourcearrays, &
                        myrank,islice_selected_rec,seismotype, &
                        xi_receiver,gamma_receiver, &
                        network_name,station_name,GPU_MODE

  use specfem_par_gpu, only: source_adjointe

  implicit none

  ! local parameters
  integer :: irec,irec_local
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: adj_sourcearray
  character(len=MAX_STRING_LEN) :: adj_source_file

  ! number of adjoint receivers in this slice
  nadj_rec_local = 0

  ! note: SIMULATION_TYPE == 2 for "pure" adjoint simulations is not fully supported yet,
  !       the run would stop after initialization.
  !       it is left here for compatibility with 3D versions however and in case we plan to support it in future...

  ! adjoint calculation
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'adjoint sources:'
      call flush_IMAIN()
    endif

    ! for GPU-version
    if (GPU_MODE) then
      allocate(source_adjointe(nrecloc,NSTEP,2))
    endif

    ! counts number of adjoint sources in this slice
    do irec = 1,nrec
      ! counts local adjoint receiver stations
      if (myrank == islice_selected_rec(irec)) then
        ! check that the source proc number is okay
        if (islice_selected_rec(irec) < 0 .or. islice_selected_rec(irec) > NPROC-1) then
          call exit_MPI(myrank,'something is wrong with the source proc number in adjoint simulation')
        endif
        ! counter
        nadj_rec_local = nadj_rec_local + 1
      endif
    enddo

    ! array for all adjoint sources
    if (nadj_rec_local > 0) then
      allocate(adj_sourcearrays(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLZ))
    else
      allocate(adj_sourcearrays(1,1,1,1,1))
    endif
    adj_sourcearrays(:,:,:,:,:) = 0._CUSTOM_REAL

    ! reads in adjoint source files
    if (.not. SU_FORMAT) then
      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  reading ASCII adjoint source files'
        call flush_IMAIN()
      endif

      ! temporary array
      allocate(adj_sourcearray(NSTEP,NDIM,NGLLX,NGLLZ))

      ! reads in ascii adjoint source files **.adj
      irec_local = 0
      do irec = 1, nrec
        ! compute only adjoint source arrays in the local proc
        if (myrank == islice_selected_rec(irec)) then
          irec_local = irec_local + 1
          adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))

          call compute_arrays_adj_source(xi_receiver(irec),gamma_receiver(irec),irec_local,adj_source_file,adj_sourcearray)

          adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
        endif
      enddo
      ! checks
      if (irec_local /= nadj_rec_local) stop 'Error invalid number of local adjoint sources found'
      ! frees temporary array
      deallocate(adj_sourcearray)
    else
      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  reading SU-format adjoint source files'
        call flush_IMAIN()
      endif

      ! (SU_FORMAT)
      call compute_arrays_adj_source_SU(seismotype)
    endif

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  number of adjoint sources = ',nrec
      call flush_IMAIN()
    endif
  else
    ! dummy allocation
    allocate(adj_sourcearrays(1,1,1,1,1))
  endif ! SIMULATION_TYPE == 3

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_adjoint_sources

!
!-----------------------------------------------------------------------------------------
!

  subroutine setup_source_receiver_tangentials()

! tangential computation

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: PI,HUGEVAL
  use specfem_par

  implicit none

  ! local parameters
  integer :: i,irec,i_source
  integer :: ier,nrec_alloc
  integer :: irecloc
  double precision :: x_final_receiver_dummy, z_final_receiver_dummy

  ! receiver arrays
  if (nrecloc > 0) then
    nrec_alloc = nrecloc
  else
    ! dummy allocation
    nrec_alloc = 1
  endif

  allocate(anglerec_irec(nrecloc), &
           cosrot_irec(nrecloc), &
           sinrot_irec(nrecloc), &
           rec_tangential_detection_curve(nrecloc),stat=ier)
  if (ier /= 0) stop 'Error allocating tangential arrays'

  ! checks angle
  if (rec_normal_to_surface .and. abs(anglerec) > 1.d-6) &
    stop 'anglerec should be zero when receivers are normal to the topography'

  ! convert receiver angle to radians
  anglerec = anglerec * pi / 180.d0

  anglerec_irec(:) = anglerec
  cosrot_irec(:) = cos(anglerec_irec(:))
  sinrot_irec(:) = sin(anglerec_irec(:))

  ! tangential computation
  ! for receivers
  if (rec_normal_to_surface) then
    irecloc = 0
    do irec = 1, nrec
      if (myrank == islice_selected_rec(irec)) then
        irecloc = irecloc + 1
        distmin = HUGEVAL
        do i = 1, nnodes_tangential_curve
          dist_current = sqrt((x_final_receiver(irec)-nodes_tangential_curve(1,i))**2 + &
             (z_final_receiver(irec)-nodes_tangential_curve(2,i))**2)
          if (dist_current < distmin) then
            n1_tangential_detection_curve = i
            distmin = dist_current
          endif
       enddo

       rec_tangential_detection_curve(irecloc) = n1_tangential_detection_curve
       call tri_quad(n_tangential_detection_curve, n1_tangential_detection_curve, &
                     nnodes_tangential_curve)

       call compute_normal_vector( anglerec_irec(irecloc), &
                                   nodes_tangential_curve(1,n_tangential_detection_curve(1)), &
                                   nodes_tangential_curve(1,n_tangential_detection_curve(2)), &
                                   nodes_tangential_curve(1,n_tangential_detection_curve(3)), &
                                   nodes_tangential_curve(1,n_tangential_detection_curve(4)), &
                                   nodes_tangential_curve(2,n_tangential_detection_curve(1)), &
                                   nodes_tangential_curve(2,n_tangential_detection_curve(2)), &
                                   nodes_tangential_curve(2,n_tangential_detection_curve(3)), &
                                   nodes_tangential_curve(2,n_tangential_detection_curve(4)) )
      endif

    enddo
    cosrot_irec(:) = cos(anglerec_irec(:))
    sinrot_irec(:) = sin(anglerec_irec(:))
  endif

  ! for the source
  if (force_normal_to_surface) then
    do i_source = 1,NSOURCES
      if (myrank == islice_selected_source(i_source)) then
        distmin = HUGEVAL
        do i = 1, nnodes_tangential_curve
          dist_current = sqrt((coord(1,iglob_source(i_source))-nodes_tangential_curve(1,i))**2 + &
                              (coord(2,iglob_source(i_source))-nodes_tangential_curve(2,i))**2)
          if (dist_current < distmin) then
            n1_tangential_detection_curve = i
            distmin = dist_current
          endif
        enddo

        call tri_quad(n_tangential_detection_curve, n1_tangential_detection_curve, &
                     nnodes_tangential_curve)

        ! in the case of a source force vector
        ! users can give an angle with respect to the normal to the topography surface,
        ! in which case we must compute the normal to the topography
        ! and add it the existing rotation angle
        call compute_normal_vector( anglesource(i_source), &
                                    nodes_tangential_curve(1,n_tangential_detection_curve(1)), &
                                    nodes_tangential_curve(1,n_tangential_detection_curve(2)), &
                                    nodes_tangential_curve(1,n_tangential_detection_curve(3)), &
                                    nodes_tangential_curve(1,n_tangential_detection_curve(4)), &
                                    nodes_tangential_curve(2,n_tangential_detection_curve(1)), &
                                    nodes_tangential_curve(2,n_tangential_detection_curve(2)), &
                                    nodes_tangential_curve(2,n_tangential_detection_curve(3)), &
                                    nodes_tangential_curve(2,n_tangential_detection_curve(4)) )

        source_courbe_eros(i_source) = n1_tangential_detection_curve
        if (myrank == 0 .and. myrank == islice_selected_source(i_source)) then
          source_courbe_eros(i_source) = n1_tangential_detection_curve
          anglesource_recv = anglesource(i_source)
#ifdef USE_MPI
        else if (myrank == 0) then
          call MPI_recv(source_courbe_eros(i_source),1,MPI_INTEGER, &
                        MPI_ANY_SOURCE,42,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
          call MPI_recv(anglesource_recv,1,MPI_DOUBLE_PRECISION, &
                        MPI_ANY_SOURCE,43,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

        else if (myrank == islice_selected_source(i_source)) then
          call MPI_send(n1_tangential_detection_curve,1,MPI_INTEGER,0,42,MPI_COMM_WORLD,ier)
          call MPI_send(anglesource(i_source),1,MPI_DOUBLE_PRECISION,0,43,MPI_COMM_WORLD,ier)
#endif
        endif

#ifdef USE_MPI
        call bcast_all_singledp(anglesource_recv)
        anglesource(i_source) = anglesource_recv
#endif


      endif
    enddo ! do i_source= 1,NSOURCES
  endif !  if (force_normal_to_surface)

! CHRIS --- how to deal with multiple source. Use first source now. ---
! compute distance from source to receivers following the curve
  if (force_normal_to_surface .and. rec_normal_to_surface) then
    dist_tangential_detection_curve(source_courbe_eros(1)) = 0
    do i = source_courbe_eros(1)+1, nnodes_tangential_curve
      dist_tangential_detection_curve(i) = dist_tangential_detection_curve(i-1) + &
          sqrt((nodes_tangential_curve(1,i)-nodes_tangential_curve(1,i-1))**2 + &
          (nodes_tangential_curve(2,i)-nodes_tangential_curve(2,i-1))**2)
    enddo
    dist_tangential_detection_curve(1) = dist_tangential_detection_curve(nnodes_tangential_curve) + &
         sqrt((nodes_tangential_curve(1,1)-nodes_tangential_curve(1,nnodes_tangential_curve))**2 + &
         (nodes_tangential_curve(2,1)-nodes_tangential_curve(2,nnodes_tangential_curve))**2)
    do i = 2, source_courbe_eros(1)-1
      dist_tangential_detection_curve(i) = dist_tangential_detection_curve(i-1) + &
          sqrt((nodes_tangential_curve(1,i)-nodes_tangential_curve(1,i-1))**2 + &
          (nodes_tangential_curve(2,i)-nodes_tangential_curve(2,i-1))**2)
    enddo
    do i = source_courbe_eros(1)-1, 1, -1
      dist_current = dist_tangential_detection_curve(i+1) + &
          sqrt((nodes_tangential_curve(1,i)-nodes_tangential_curve(1,i+1))**2 + &
          (nodes_tangential_curve(2,i)-nodes_tangential_curve(2,i+1))**2)
      if (dist_current < dist_tangential_detection_curve(i)) then
        dist_tangential_detection_curve(i) = dist_current
      endif
    enddo
    dist_current = dist_tangential_detection_curve(1) + &
       sqrt((nodes_tangential_curve(1,1)-nodes_tangential_curve(1,nnodes_tangential_curve))**2 + &
       (nodes_tangential_curve(2,1)-nodes_tangential_curve(2,nnodes_tangential_curve))**2)
    if (dist_current < dist_tangential_detection_curve(nnodes_tangential_curve)) then
      dist_tangential_detection_curve(nnodes_tangential_curve) = dist_current
    endif
    do i = nnodes_tangential_curve-1, source_courbe_eros(1)+1, -1
      dist_current = dist_tangential_detection_curve(i+1) + &
          sqrt((nodes_tangential_curve(1,i)-nodes_tangential_curve(1,i+1))**2 + &
          (nodes_tangential_curve(2,i)-nodes_tangential_curve(2,i+1))**2)
      if (dist_current < dist_tangential_detection_curve(i)) then
        dist_tangential_detection_curve(i) = dist_current
      endif
    enddo

    if (myrank == 0) then
      open(unit=11,file='OUTPUT_FILES/dist_rec_tangential_detection_curve', &
            form='formatted', status='unknown')
    endif

    irecloc = 0
    do irec = 1,nrec

      if (myrank == 0) then
        if (myrank == islice_selected_rec(irec)) then
          irecloc = irecloc + 1
          n1_tangential_detection_curve = rec_tangential_detection_curve(irecloc)
          x_final_receiver_dummy = x_final_receiver(irec)
          z_final_receiver_dummy = z_final_receiver(irec)
#ifdef USE_MPI
        else

          call MPI_RECV(n1_tangential_detection_curve,1,MPI_INTEGER, &
             islice_selected_rec(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
          call MPI_RECV(x_final_receiver_dummy,1,MPI_DOUBLE_PRECISION, &
             islice_selected_rec(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
          call MPI_RECV(z_final_receiver_dummy,1,MPI_DOUBLE_PRECISION, &
             islice_selected_rec(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

#endif
        endif

#ifdef USE_MPI
      else
        if (myrank == islice_selected_rec(irec)) then
          irecloc = irecloc + 1
          call MPI_SEND(rec_tangential_detection_curve(irecloc),1,MPI_INTEGER,0,irec,MPI_COMM_WORLD,ier)
          call MPI_SEND(x_final_receiver(irec),1,MPI_DOUBLE_PRECISION,0,irec,MPI_COMM_WORLD,ier)
          call MPI_SEND(z_final_receiver(irec),1,MPI_DOUBLE_PRECISION,0,irec,MPI_COMM_WORLD,ier)
        endif
#endif

      endif
      if (myrank == 0) then
        write(11,*) dist_tangential_detection_curve(n1_tangential_detection_curve)
        write(12,*) x_final_receiver_dummy
        write(13,*) z_final_receiver_dummy
      endif
    enddo

    if (myrank == 0) then
      close(11)
      close(12)
      close(13)
    endif

  endif ! force_normal_to_surface

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_source_receiver_tangentials

!
!-----------------------------------------------------------------------------------------
!

  subroutine setup_source_interpolation()

  use constants, only: NDIM,NGLLX,NGLLZ,NGLJ,ZERO,CUSTOM_REAL

  use specfem_par, only: myrank,nspec,NSOURCES,source_type,anglesource,P_SV, &
    sourcearrays,Mxx,Mxz,Mzz, &
    ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic, &
    ispec_selected_source,islice_selected_source, &
    xi_source,gamma_source, &
    xix,xiz,gammax,gammaz,xigll,zigll, &
    hxis_store,hgammas_store,hxis,hpxis,hgammas,hpgammas, &
    AXISYM,is_on_the_axis,xiglj

  implicit none

  ! local parameters
  integer :: i_source,ispec,i,j,ier
  double precision :: hlagrange
  ! single source array
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  ! allocates Lagrange interpolators for sources
  allocate(hxis_store(NSOURCES,NGLLX), &
           hgammas_store(NSOURCES,NGLLZ),stat=ier)
  if (ier /= 0) stop 'Error allocating source h**_store arrays'

  ! initializes
  hxis_store(:,:) = ZERO
  hgammas_store(:,:) = ZERO

  sourcearrays(:,:,:,:) = 0._CUSTOM_REAL

  ! define and store Lagrange interpolators at all the sources
  do i_source = 1,NSOURCES

    if (myrank == islice_selected_source(i_source)) then

      ! element containing source
      ispec = ispec_selected_source(i_source)

      ! Lagrange interpolators
      if (AXISYM) then
        if (is_on_the_axis(ispec)) then
          call lagrange_any(xi_source(i_source),NGLJ,xiglj,hxis,hpxis)
          !do j = 1,NGLJ ! ABAB same result with that loop, this is good
          !  hxis(j) = hglj(j-1,xi_source(i_source),xiglj,NGLJ)
          !enddo
        else
          call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
        endif
      else
        call lagrange_any(xi_source(i_source),NGLLX,xigll,hxis,hpxis)
      endif
      call lagrange_any(gamma_source(i_source),NGLLZ,zigll,hgammas,hpgammas)

      ! stores Lagrangians for source
      hxis_store(i_source,:) = hxis(:)
      hgammas_store(i_source,:) = hgammas(:)

      sourcearray(:,:,:) = 0._CUSTOM_REAL

      ! computes source arrays
      select case (source_type(i_source))
      case (1)
        ! collocated force source
        do j = 1,NGLLZ
          do i = 1,NGLLX
            hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)

            ! source element is acoustic
            if (ispec_is_acoustic(ispec)) then
              sourcearray(:,i,j) = hlagrange
            endif

            ! source element is elastic
            if (ispec_is_elastic(ispec)) then
              if (P_SV) then
                ! P_SV case
                sourcearray(1,i,j) = - sin(anglesource(i_source)) * hlagrange
                sourcearray(2,i,j) =   cos(anglesource(i_source)) * hlagrange
              else
                ! SH case (membrane)
                sourcearray(:,i,j) = hlagrange
              endif
            endif

            ! source element is poroelastic
            if (ispec_is_poroelastic(ispec)) then
              sourcearray(1,i,j) = - sin(anglesource(i_source)) * hlagrange
              sourcearray(2,i,j) =   cos(anglesource(i_source)) * hlagrange
            endif

          enddo
        enddo

      case (2)
        ! moment-tensor source
        call compute_arrays_source(ispec,xi_source(i_source),gamma_source(i_source),sourcearray, &
                                   Mxx(i_source),Mzz(i_source),Mxz(i_source),xix,xiz,gammax,gammaz,xigll,zigll,nspec)
        ! checks source
        if (ispec_is_acoustic(ispec)) then
          call exit_MPI(myrank,'cannot have moment tensor source in acoustic element')
        endif

        ! checks wave type
        if (ispec_is_elastic(ispec)) then
          if (.not. P_SV ) call exit_MPI(myrank,'cannot have moment tensor source in SH (membrane) waves calculation')
        endif

      end select

      ! stores sourcearray for all sources
      sourcearrays(i_source,:,:,:) = sourcearray(:,:,:)

    endif
  enddo

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_source_interpolation


!
!-----------------------------------------------------------------------------------------
!

  subroutine setup_receiver_interpolation()

  use constants, only: NGLLX,NGLLZ,NGLJ

  use specfem_par, only: myrank,nrec,nrecloc, &
    ispec_selected_rec,islice_selected_rec, &
    xigll,zigll, &
    hxir_store,hgammar_store,xi_receiver,gamma_receiver,hxir,hpxir,hgammar,hpgammar, &
    AXISYM,is_on_the_axis,xiglj

  use specfem_par_gpu, only: xir_store_loc,gammar_store_loc

  implicit none

  ! local parameters
  integer :: irec,irec_local,ier

  ! allocate Lagrange interpolators for receivers
  allocate(hxir_store(nrec,NGLLX), &
           hgammar_store(nrec,NGLLZ),stat=ier)
  if (ier /= 0) stop 'Error allocating receiver h**_store arrays'

  allocate(xir_store_loc(nrecloc,NGLLX), &
           gammar_store_loc(nrecloc,NGLLZ),stat=ier)
  if (ier /= 0) stop 'Error allocating local receiver h**_store arrays'

  ! define and store Lagrange interpolators at all the receivers
  irec_local = 0
  do irec = 1,nrec
    if (AXISYM) then
      if (is_on_the_axis(ispec_selected_rec(irec)) .and. myrank == islice_selected_rec(irec)) then
        call lagrange_any(xi_receiver(irec),NGLJ,xiglj,hxir,hpxir)
        !do j = 1,NGLJ ! ABAB Exactly same result with that loop. This is good
        !  hxir2(j) = hglj(j-1,xi_receiver(irec),xiglj,NGLJ)
        !  print *,hxir(j),' = ',hxir2(j)
        !enddo
      else
        call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
      endif
    else
      call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
      !do j = 1,NGLLX !do j = 1,NGLJ ! ABAB Exactly same result with that loop. This is good
      !  hxir2(j) = hgll(j-1,xi_receiver(irec),xigll,NGLLX)
      !  print *,hxir(j),' = ',hxir2(j)
      !enddo
      ! Defined:
      !  double precision, dimension(NGLL) :: hxir2
      !  double precision, external :: hglj,hgll
    endif
    call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)

    ! stores Lagrangians for receiver
    hxir_store(irec,:) = hxir(:)
    hgammar_store(irec,:) = hgammar(:)

    ! local receivers in this slice
    if (myrank == islice_selected_rec(irec)) then
      irec_local = irec_local + 1
      xir_store_loc(irec_local,:) = hxir(:)
      gammar_store_loc(irec_local,:) = hgammar(:)
    endif
  enddo

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_receiver_interpolation

