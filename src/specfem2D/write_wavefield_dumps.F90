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

  subroutine write_wavefield_dumps()

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IMAIN,SIZE_REAL,NGLLX,NGLLZ,OUTPUT_FILES

  use specfem_par, only: myrank,nglob,nspec, &
                         ibool,coord,P_SV,it,SIMULATION_TYPE, &
                         potential_acoustic,displ_elastic,displs_poroelastic, &
                         potential_dot_acoustic,veloc_elastic,velocs_poroelastic, &
                         potential_dot_dot_acoustic,accel_elastic,accels_poroelastic, &
                         NPROC

  use specfem_par_movie, only: this_is_the_first_time_we_dump, &
                               mask_ibool,imagetype_wavefield_dumps, &
                               use_binary_for_wavefield_dumps,vector_field_display, &
                               dump_recv_counts, dump_send, dump_recv, dump_gather, dump_write, &
                               mask_duplicate, dump_duplicate_send, dump_duplicate_recv, dump_duplicate_gather

  use specfem_par, only: ninterface, ibool_interfaces_ext_mesh

  implicit none

  !local variables
  integer :: ispec,iglob,icounter,nb_of_values_to_save
  integer :: ier
  ! name of wavefield snapshot file
  character(len=150) :: wavefield_file

  integer :: iproc              ! counter over MPI process
  integer :: gcounter           ! index counter for variable send/receive counts
  integer :: ii,jj,kk              ! Loop counters

  integer :: dummy

#ifndef USE_MPI
! To avoid warnings by the compiler about unused variables in case of a serial code.

  iproc = 0
  dump_recv_counts = 0
  gcounter = 0
  ier = 0

  allocate(dump_gather(1,1))
  deallocate(dump_gather)

  allocate(dump_duplicate_gather(1))
  deallocate(dump_duplicate_gather)

  allocate(mask_duplicate(1))
  deallocate(mask_duplicate)
  dummy = NPROC
#else
  dummy = 0
#endif



  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Dumping the wave field to a file for time step ',it
    call flush_IMAIN()
  endif


  if (this_is_the_first_time_we_dump) then

    if (.not. allocated(mask_ibool)) allocate(mask_ibool(nglob))
    ! These has been allocated in prepare_timerun.F90
    dump_duplicate_send = .false.
    dump_duplicate_recv = .false.
    ! Loop over all elements and mask inner duplicate points.
    icounter = 0
    mask_ibool(:) = .false.
    do ispec = 1,nspec
      do jj = 1,NGLLZ
        do ii = 1,NGLLX
          iglob = ibool(ii,jj,ispec)
          if (.not. mask_ibool(iglob)) then
            icounter = icounter + 1
            mask_ibool(iglob) = .true.

            ! Storing directly in recv buffer for proc 0
            if (myrank == 0) then
              dump_recv(1,icounter) = coord(1,iglob)
              dump_recv(2,icounter) = coord(2,iglob)
              ! Loop through interfaces and if current node is marked as being on an interface, mark it as duplicate.
              do kk = 1, ninterface
                if (any(iglob == ibool_interfaces_ext_mesh(:,kk))) dump_duplicate_recv(icounter) = .true.
              enddo
            else
              ! Send data to master
              dump_send(1,icounter) = coord(1,iglob)
              dump_send(2,icounter) = coord(2,iglob)
              do kk = 1, ninterface
                if (any(iglob == ibool_interfaces_ext_mesh(:,kk))) dump_duplicate_send(icounter) = .true.
              enddo
            endif
          endif
        enddo
      enddo
    enddo

#ifdef USE_MPI
    if (NPROC > 1) then
      if (myrank == 0) then
        ! Master collects.
        ! Collect first receive counts to allocate gather array.
        dump_recv_counts(0) = icounter
        do iproc = 1, NPROC-1
          call MPI_RECV (dump_recv_counts(iproc), 1, MPI_INTEGER, iproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        enddo
        ! As this is first time, dump_gather is not yet allocated.
        allocate(dump_gather(2,sum(dump_recv_counts)))
        dump_gather = 0.0
        allocate(dump_duplicate_gather(sum(dump_recv_counts)))
        dump_duplicate_gather = .false.

        ! Start gathering with proc 0 data
        dump_gather(:,1:dump_recv_counts(0)) = dump_recv(:,1:dump_recv_counts(0))
        dump_duplicate_gather(1:dump_recv_counts(0)) = dump_duplicate_recv(1:dump_recv_counts(0))
        gcounter = dump_recv_counts(0)
        ! Collect from other process.
        do iproc = 1, NPROC-1
          call MPI_RECV (dump_recv(1,1), 2*dump_recv_counts(iproc), &
               MPI_DOUBLE_PRECISION, iproc, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

          dump_gather(:,gcounter+1:gcounter+dump_recv_counts(iproc)) = dump_recv(:,1:dump_recv_counts(iproc))

          call MPI_RECV (dump_duplicate_recv(1), dump_recv_counts(iproc), &
               MPI_LOGICAL, iproc, 45, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

          dump_duplicate_gather(gcounter+1:gcounter+dump_recv_counts(iproc)) = dump_duplicate_recv(1:dump_recv_counts(iproc))
          ! Update index
          gcounter = gcounter + dump_recv_counts(iproc)
        enddo

        call mask_duplicates()

        ! As array mask_duplicate is now filled, this array is no longer needed.
        deallocate(dump_duplicate_gather)

        call mask_write_matrix()

      else
        call MPI_SEND (icounter, 1, MPI_INTEGER, 0, 43, MPI_COMM_WORLD, ier)
        call MPI_SEND (dump_send(1,1), 2*icounter, &
             MPI_DOUBLE_PRECISION, 0, 44, MPI_COMM_WORLD, ier)
        call MPI_SEND (dump_duplicate_send(1), icounter, &
             MPI_LOGICAL, 0, 45, MPI_COMM_WORLD, ier)

      endif ! if (myrank == 0)

    else
      ! In this case, all data is located in dump_recv (as myrank = 0) as run in single process.
      ! Rather ineffecient to copy data into another variable just to write, but helas.
      ! As dump_write gets allocated in mask_write_matrix, which is not called in this case, the array must be allocated here.
      allocate(dump_write(2, size(dump_recv,2)))
      dump_write = dump_recv(:,1:icounter)
    endif ! if (NPROC > 1)
#else
    ! Array must be allocated first time we dump.
    allocate(dump_write(2, icounter))
    dump_write = dump_recv(:,1:icounter)
#endif

    this_is_the_first_time_we_dump = .false.

    ! Proc 0 handles mesh coordinate file write. Done only first time we dump.
    if (myrank == 0) then
      call write_file_grid()
    endif

  endif ! if (this_is_the_first_time_we_dump)

  ! Prepare the requested data for dump.
  if (imagetype_wavefield_dumps == 1) then
    if (myrank == 0) write(IMAIN,*) 'Dumping the displacement vector...'
      call compute_vector_whole_medium(potential_acoustic,displ_elastic,displs_poroelastic)

  else if (imagetype_wavefield_dumps == 2) then
    if (myrank == 0) write(IMAIN,*) 'Dumping the velocity vector...'
    call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,velocs_poroelastic)

  else if (imagetype_wavefield_dumps == 3) then
    if (myrank == 0) write(IMAIN,*) 'Dumping the acceleration vector...'
    call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,accels_poroelastic)

  else if (imagetype_wavefield_dumps == 4 .and. P_SV) then
    if (myrank == 0) write(IMAIN,*) 'Dumping the pressure field...'
    call compute_pressure_whole_medium(1)

  else if (imagetype_wavefield_dumps == 4 .and. .not. P_SV) then
    call exit_MPI(myrank,'cannot dump the pressure field for SH (membrane) waves')

  else
    call exit_MPI(myrank,'wrong type of flag for wavefield dumping')
  endif

  ! Copy the requested data to send arrays
  icounter = 0
  mask_ibool(:) = .false.
  do ispec = 1,nspec
    do jj = 1,NGLLZ
      do ii = 1,NGLLX
        iglob = ibool(ii,jj,ispec)
        if (.not. mask_ibool(iglob)) then
          icounter = icounter + 1
          mask_ibool(iglob) = .true.

          if (myrank == 0) then
            dump_recv(1,icounter) = vector_field_display(1,iglob)
            dump_recv(2,icounter) = vector_field_display(2,iglob)
          else
            dump_send(1,icounter) = vector_field_display(1,iglob)
            dump_send(2,icounter) = vector_field_display(2,iglob)
          endif

        endif
      enddo ! do ii
    enddo ! do jj
  enddo ! do ispec

#ifdef USE_MPI
  if (NPROC > 1) then
    if (myrank == 0) then
      ! Master collects
      ! Start gathering with proc 0 data
      dump_gather(:,1:dump_recv_counts(0)) = dump_recv(:,1:dump_recv_counts(0))
      gcounter = dump_recv_counts(0)
      ! The receiver counts dump_recv_counts has already been found and stored and these does not change from first to subsequent dumps.
      do iproc = 1, NPROC-1
        call MPI_RECV (dump_recv(1,1), 2*dump_recv_counts(iproc), &
             MPI_DOUBLE_PRECISION, iproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        dump_gather(:,gcounter+1:gcounter+dump_recv_counts(iproc)) = dump_recv(:,1:dump_recv_counts(iproc))
        gcounter = gcounter + dump_recv_counts(iproc)
      enddo

      ! Mask array mask_duplicate is generated first time we dump.
      call mask_write_matrix()

    else
      ! Send collected vector field to master
      call MPI_SEND (dump_send(1,1), 2*icounter, &
           MPI_DOUBLE_PRECISION, 0, 43, MPI_COMM_WORLD, ier)
    endif ! if (myrank == 0)
  else
    dump_write = dump_recv(:,1:icounter)
  endif ! if (NPROC > 1)
#else
  ! Array dump_write has already been allocated first time we dump.
  dump_write = dump_recv(:,1:icounter)
#endif

  ! Proc 0 handles file write.
  if (myrank == 0) then
    call write_file_dump()
  endif


  if (myrank == 0) then
    write(IMAIN,*) 'Wave field dumped'
    call flush_IMAIN()
  endif

  call synchronize_all()

  return

  contains
  !
  ! ------------------------------------------------------------
  !
  subroutine write_file_grid()

  integer :: ii

  ! Open file handle
  if (use_binary_for_wavefield_dumps) then
    wavefield_file = 'OUTPUT_FILES/wavefield_grid_for_dumps.bin'
    open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='replace', &
         action='write',recl=2*SIZE_REAL,iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield_grid_for_dumps,bin')
  else
    wavefield_file = 'OUTPUT_FILES/wavefield_grid_for_dumps.txt'
    open(unit=27,file=wavefield_file,status='replace',action='write',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield_grid_for_dumps.txt')
  endif

  ! Write file content

  do ii = 1, size(dump_write,2)
    if (use_binary_for_wavefield_dumps) then
      write(27,rec=ii) sngl(dump_write(1,ii)), sngl(dump_write(2,ii))
    else
      write(27,'(2e16.6)') dump_write(1,ii), dump_write(2,ii)
    endif
  enddo

  ! Close file
  close(27)

  return
  end subroutine write_file_grid
  !
  ! ------------------------------------------------------------
  !
  subroutine write_file_dump()

  integer :: ii

  ! Open file handle
  if (use_binary_for_wavefield_dumps) then
    if (P_SV .and. .not. imagetype_wavefield_dumps == 4) then
      nb_of_values_to_save = 2
    else
      nb_of_values_to_save = 1
    endif
    write(wavefield_file,"(a,i7.7,'_',i2.2,a)") trim(OUTPUT_FILES)//'wavefield',it,SIMULATION_TYPE,'.bin'
    open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='replace', &
         action='write',recl=nb_of_values_to_save*SIZE_REAL,iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield**.bin')
  else
    write(wavefield_file,"(a,i7.7,'_',i2.2,a)") trim(OUTPUT_FILES)//'wavefield',it,SIMULATION_TYPE,'.txt'
    open(unit=27,file=wavefield_file,status='replace',action='write',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield**.txt')
  endif

  ! Write file content

  do ii=1, size(dump_write, 2)

    if (use_binary_for_wavefield_dumps) then
      if (P_SV) then
        ! P-SV waves
        if (imagetype_wavefield_dumps == 4) then
          ! by convention we use the 2. component of the array to store the pressure above
          write(27,rec=ii) sngl(dump_write(2,ii))
        else
          write(27,rec=ii) sngl(dump_write(1,ii)), sngl(dump_write(2,ii))
        endif
      else
        ! SH case
        write(27,rec=ii) sngl(dump_write(1,ii))
      endif
    else
      if (P_SV) then
        ! P-SV waves
        if (imagetype_wavefield_dumps == 4) then
          ! by convention we use the 2. component of the array to store the pressure above
          write(27,*) sngl(dump_write(2,ii))
        else
          write(27,*) sngl(dump_write(1,ii)), sngl(dump_write(2,ii))
        endif
      else
        ! SH case
        write(27,*) sngl(dump_write(1,ii))
      endif
    endif
  enddo

  ! Close file
  close(27)

  return
  end subroutine write_file_dump
  !
  ! ------------------------------------------------------------
  !
  subroutine mask_duplicates()
    !

  implicit none

  integer :: ii,jj
  integer :: nduplicate
  integer, dimension(:), allocatable :: duplicate_index, duplicate_index_pack
  logical, dimension(:), allocatable :: duplicate_index_mask

  ! This variable is only generated once
  allocate(mask_duplicate(size(dump_gather,2)))

  ! Initialise all elements to be included. Duplicate elements will later be set to false.
  mask_duplicate(:) = .true.

  ! Nothing to do then
  if (count(dump_duplicate_gather) == 0) return

  ! Counter for duplicate removal cycles.
  ii = 0
  do while (count(dump_duplicate_gather) > 0)
    ii = ii + 1
    ! Indices of duplicates in gather array
    allocate(duplicate_index(count(dump_duplicate_gather)))
    ! Mask vector for duplicates indices
    allocate(duplicate_index_mask(count(dump_duplicate_gather)))
    ! Might be inefficient to reform the entire index array every time, but avoids storage.
    duplicate_index = pack([(jj, jj=1, size(dump_gather, 2))], dump_duplicate_gather)
    ! Search for duplicates of first entry still marked as duplicates.
    duplicate_index_mask = dump_gather(1,duplicate_index(1)) == dump_gather(1,duplicate_index) .and. &
                           dump_gather(2,duplicate_index(1)) == dump_gather(2,duplicate_index)
    nduplicate = count(duplicate_index_mask)
    ! The reduced (masked) duplicate indices (i.e. excluding the first duplicate that is being searched for)
    allocate(duplicate_index_pack(nduplicate))
    ! Indices of duplicates set
    duplicate_index_pack = pack(duplicate_index, duplicate_index_mask)
    ! Set these entries as being handled
    dump_duplicate_gather(duplicate_index_pack) = .false.
    ! The current duplicated handled is first entry, so mask this for writing to file.
    mask_duplicate(duplicate_index_pack(1)) = .true.
    ! Remove the rest from writing
    mask_duplicate(duplicate_index_pack(2:nduplicate)) = .false.
    ! Finally deallocate for generation of new indices in next cycle
    deallocate(duplicate_index)
    deallocate(duplicate_index_mask)
    deallocate(duplicate_index_pack)
    ! Final print out message
    if (count(dump_duplicate_gather) == 0) print *, ' Number of duplicates found: ', ii
  enddo

  return
  end subroutine mask_duplicates
  !
  ! ------------------------------------------------------------
  !
  subroutine mask_write_matrix()
  ! Masks the gathered data to exclude duplicates and packs to the reduced size write arrays.
  implicit none

  if (.not. allocated(dump_write)) then
    allocate(dump_write(2,count(mask_duplicate)))
  endif

  dump_write = 0.0
  dump_write(1,:) = pack(dump_gather(1,:), mask_duplicate)
  dump_write(2,:) = pack(dump_gather(2,:), mask_duplicate)

  return
  end subroutine mask_write_matrix
end subroutine write_wavefield_dumps
