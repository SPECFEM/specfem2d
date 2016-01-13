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

  ! pre-compute lagrangians
  call setup_source_receiver_interpolation()

  ! synchronizes processes
  call synchronize_all()

  end subroutine setup_sources_receivers

!
!----------------------------------------------------------------------------
!

  subroutine setup_sources()

  use constants,only: NGLLX,NGLLZ,NDIM,IMAIN,IIN,MAX_STRING_LEN

  use specfem_par, only: NSOURCES,initialfield,source_type, &
                         coord,ibool,nglob,nspec,nelem_acoustic_surface,acoustic_surface, &
                         ispec_is_elastic,ispec_is_poroelastic, &
                         x_source,z_source,ispec_selected_source, &
                         is_proc_source,nb_proc_source, &
                         sourcearray,Mxx,Mzz,Mxz, &
                         xix,xiz,gammax,gammaz,xigll,zigll,npgeo, &
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

  do i_source= 1,NSOURCES

    if (source_type(i_source) == 1) then

      ! collocated force source
      call locate_source_force(ibool,coord,nspec,nglob,xigll,zigll,x_source(i_source),z_source(i_source), &
                               ispec_selected_source(i_source),is_proc_source(i_source),nb_proc_source(i_source), &
                               NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo, &
                               iglob_source(i_source))

      ! check that acoustic source is not exactly on the free surface because pressure is zero there
      if (is_proc_source(i_source) == 1) then
        do ispec_acoustic_surface = 1,nelem_acoustic_surface
          ispec = acoustic_surface(1,ispec_acoustic_surface)
          ixmin = acoustic_surface(2,ispec_acoustic_surface)
          ixmax = acoustic_surface(3,ispec_acoustic_surface)
          izmin = acoustic_surface(4,ispec_acoustic_surface)
          izmax = acoustic_surface(5,ispec_acoustic_surface)
          if (.not. ispec_is_elastic(ispec) .and. .not. ispec_is_poroelastic(ispec) .and. &
            ispec == ispec_selected_source(i_source)) then
            if ((izmin==1 .and. izmax==1 .and. ixmin==1 .and. ixmax==NGLLX .and. &
                gamma_source(i_source) < -0.99d0) .or.&
                (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==NGLLX .and. &
                gamma_source(i_source) > 0.99d0) .or.&
                (izmin==1 .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==1 .and. &
                xi_source(i_source) < -0.99d0) .or.&
                (izmin==1 .and. izmax==NGLLZ .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
                xi_source(i_source) > 0.99d0) .or.&
                (izmin==1 .and. izmax==1 .and. ixmin==1 .and. ixmax==1 .and. &
                gamma_source(i_source) < -0.99d0 .and. xi_source(i_source) < -0.99d0) .or.&
                (izmin==1 .and. izmax==1 .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
                gamma_source(i_source) < -0.99d0 .and. xi_source(i_source) > 0.99d0) .or.&
                (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==1 .and. &
                gamma_source(i_source) > 0.99d0 .and. xi_source(i_source) < -0.99d0) .or.&
                (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
                gamma_source(i_source) > 0.99d0 .and. xi_source(i_source) > 0.99d0)) then
              call exit_MPI(myrank,'an acoustic source cannot be located exactly '// &
                            'on the free surface because pressure is zero there')
            endif
          endif
        enddo
      endif

    else if (source_type(i_source) == 2) then
      ! moment-tensor source
      call locate_source_moment_tensor(ibool,coord,nspec,nglob,xigll,zigll,x_source(i_source),z_source(i_source), &
             ispec_selected_source(i_source),is_proc_source(i_source),nb_proc_source(i_source), &
             NPROC,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo)

      ! compute source array for moment-tensor source
      call compute_arrays_source(ispec_selected_source(i_source),xi_source(i_source),gamma_source(i_source),&
             sourcearray(i_source,1,1,1), &
             Mxx(i_source),Mzz(i_source),Mxz(i_source),xix,xiz,gammax,gammaz,xigll,zigll,nspec)

    else if (.not.initialfield) then

      call exit_MPI(myrank,'incorrect source type')

    endif

  enddo ! do i_source= 1,NSOURCES

!! DK DK this below not supported in the case of MPI yet, we should do a MPI_GATHER() of the values
!! DK DK and use "if (myrank == which_proc_receiver(irec)) then" to display the right sources
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

  use constants,only: NGLLX,NGLLZ,NDIM,IMAIN,IIN,MAX_STRING_LEN

  use specfem_par, only: coord,ibool,nglob,nspec, &
                         ispec_selected_rec, &
                         NPROC,myrank,coorg,knods,ngnod, &
                         xigll,zigll,npgeo, &
                         nrec,nrecloc,recloc,which_proc_receiver,st_xval,st_zval, &
                         xi_receiver,gamma_receiver,station_name,network_name, &
                         x_final_receiver,z_final_receiver, &
                         x_source,z_source

  implicit none

  ! Local variables
  integer :: irec,nrec_tot_found
  integer :: ier

  character(len=MAX_STRING_LEN) :: dummystring

  ! user output
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
           which_proc_receiver(nrec), &
           x_final_receiver(nrec), &
           z_final_receiver(nrec),stat=ier)
  if (ier /= 0) stop 'Error allocating receiver arrays'

  ! locate receivers in the mesh
  call locate_receivers(ibool,coord,nspec,nglob,xigll,zigll, &
                        nrec,nrecloc,recloc,which_proc_receiver,NPROC,myrank, &
                        st_xval,st_zval,ispec_selected_rec, &
                        xi_receiver,gamma_receiver,station_name,network_name, &
                        x_source(1),z_source(1), &
                        coorg,knods,ngnod,npgeo, &
                        x_final_receiver,z_final_receiver)

!! DK DK this below not supported in the case of MPI yet, we should do a MPI_GATHER() of the values
!! DK DK and use "if (myrank == which_proc_receiver(irec)) then" to display the right sources
!! DK DK and receivers carried by each mesh slice, and not fictitious values coming from other slices
  irec = 0
#ifndef USE_MPI
  if (myrank == 0) then
     ! write out actual station locations (compare with STATIONS from meshfem2D)
     ! NOTE: this will be written out even if use_existing_STATIONS = .true.
     open(unit=15,file='OUTPUT_FILES/for_information_STATIONS_actually_used',status='unknown')
     do irec = 1,nrec
        write(15,"('S',i4.4,'    AA ',f20.7,1x,f20.7,'       0.0         0.0')") &
             irec,x_final_receiver(irec),z_final_receiver(irec)
     enddo
     close(15)
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

  use constants,only: NGLLX,NGLLZ,IMAIN

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
        if ((izmin==1 .and. izmax==1 .and. ixmin==1 .and. ixmax==NGLLX .and. &
        gamma_receiver(irec) < -0.99d0) .or.&
        (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==NGLLX .and. &
        gamma_receiver(irec) > 0.99d0) .or.&
        (izmin==1 .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==1 .and. &
        xi_receiver(irec) < -0.99d0) .or.&
        (izmin==1 .and. izmax==NGLLZ .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
        xi_receiver(irec) > 0.99d0) .or.&
        (izmin==1 .and. izmax==1 .and. ixmin==1 .and. ixmax==1 .and. &
        gamma_receiver(irec) < -0.99d0 .and. xi_receiver(irec) < -0.99d0) .or.&
        (izmin==1 .and. izmax==1 .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
        gamma_receiver(irec) < -0.99d0 .and. xi_receiver(irec) > 0.99d0) .or.&
        (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==1 .and. &
        gamma_receiver(irec) > 0.99d0 .and. xi_receiver(irec) < -0.99d0) .or.&
        (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
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

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,MAX_STRING_LEN

  use specfem_par,only: nadj_rec_local,nrec,nrecloc,NSTEP,NPROC,SIMULATION_TYPE,SU_FORMAT, &
                        adj_sourcearrays, &
                        myrank,which_proc_receiver,seismotype, &
                        xi_receiver,gamma_receiver, &
                        network_name,station_name

  use specfem_par_gpu,only: source_adjointe

  implicit none

  ! local parameters
  integer :: irec,irec_local
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: adj_sourcearray
  character(len=MAX_STRING_LEN) :: adj_source_file

  ! number of adjoint receivers in this slice
  nadj_rec_local = 0

  ! adjoint calculation
  if (SIMULATION_TYPE == 3) then

    allocate(source_adjointe(nrecloc,NSTEP,2))

    ! counts number of adjoint sources in this slice
    do irec = 1,nrec
      ! counts local adjoint receiver stations
      if (myrank == which_proc_receiver(irec)) then
        ! check that the source proc number is okay
        if (which_proc_receiver(irec) < 0 .or. which_proc_receiver(irec) > NPROC-1) then
          call exit_MPI(myrank,'something is wrong with the source proc number in adjoint simulation')
        endif
        ! counter
        nadj_rec_local = nadj_rec_local + 1
      endif
    enddo

    ! array for all adjoint sources
    if (nadj_rec_local > 0)  then
      allocate(adj_sourcearrays(nadj_rec_local,NSTEP,NDIM,NGLLX,NGLLZ))
    else
      allocate(adj_sourcearrays(1,1,1,1,1))
    endif
    adj_sourcearrays(:,:,:,:,:) = 0._CUSTOM_REAL

    ! reads in adjoint source files
    if (.not. SU_FORMAT) then
      ! temporary array
      allocate(adj_sourcearray(NSTEP,NDIM,NGLLX,NGLLZ))

      ! reads in ascii adjoint source files **.adj
      if (seismotype == 1 .or. seismotype == 2 .or. seismotype == 3) then
         irec_local = 0
         do irec = 1, nrec
           ! compute only adjoint source arrays in the local proc
           if (myrank == which_proc_receiver(irec)) then
             irec_local = irec_local + 1
             adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))

             call compute_arrays_adj_source(xi_receiver(irec),gamma_receiver(irec),irec_local,adj_source_file,adj_sourcearray)

             adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
           endif
         enddo
      else if (seismotype == 4 .or. seismotype == 6) then
        ! single component
        irec_local = 0
        do irec = 1, nrec
          ! compute only adjoint source arrays in the local proc
          if (myrank == which_proc_receiver(irec)) then
            irec_local = irec_local + 1
            adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))

            call compute_arrays_adj_source(xi_receiver(irec),gamma_receiver(irec),irec_local,adj_source_file,adj_sourcearray)

            adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
          endif
        enddo
      endif
      ! frees temporary array
      deallocate(adj_sourcearray)
    else
      ! (SU_FORMAT)
      call add_adjoint_sources_SU(seismotype)
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


  subroutine add_adjoint_sources_SU(seismotype)

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLZ,NGLJ,NDIM

  use specfem_par, only: AXISYM,xiglj,is_on_the_axis, &
                         myrank, NSTEP, nrec, &
                         xi_receiver, gamma_receiver, which_proc_receiver, &
                         xigll,zigll,hxir,hgammar,hpxir,hpgammar, &
                         adj_sourcearrays, &
                         GPU_MODE, ispec_selected_rec, &
                         P_SV

  use specfem_par_gpu,only: source_adjointe

  implicit none

  integer,intent(in) :: seismotype

  ! local parameters
  integer :: i, j, irec, irec_local, ier
  character(len=MAX_STRING_LEN) :: filename

  ! SU
  integer(kind=4) :: r4head(60)
  integer(kind=2) :: header2(2)
  real(kind=4),dimension(:,:),allocatable :: adj_src_s

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: adj_sourcearray

  ! opens adjoint source files
  if (seismotype == 4 .or. seismotype == 6) then
    write(filename, "('./SEM/Up_file_single.su.adj')")
    open(111,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')
  else
    write(filename, "('./SEM/Ux_file_single.su.adj')")
    open(111,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

    write(filename, "('./SEM/Uy_file_single.su.adj')")
    open(112,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

    write(filename, "('./SEM/Uz_file_single.su.adj')")
    open(113,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
    if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')
  endif

  ! allocates temporary array
  allocate(adj_src_s(NSTEP,3))
  adj_src_s(:,:) = 0.0

  ! temporary array
  allocate(adj_sourcearray(NSTEP,NDIM,NGLLX,NGLLZ))

  irec_local = 0
  do irec = 1, nrec
    if (myrank == which_proc_receiver(irec)) then
      irec_local = irec_local + 1

      adj_sourcearray(:,:,:,:) = 0._CUSTOM_REAL

      if (seismotype == 4 .or. seismotype == 6) then
        read(111,rec=irec,iostat=ier) r4head, adj_src_s(:,1)
        if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' read error')
      else
        read(111,rec=irec,iostat=ier) r4head, adj_src_s(:,1)
        if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' read error')

        read(112,rec=irec,iostat=ier) r4head, adj_src_s(:,2)
        if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' read error')

        read(113,rec=irec,iostat=ier) r4head, adj_src_s(:,3)
        if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' read error')
      endif

      header2 = int(r4head(29), kind=2)
      if (irec==1) print *, r4head(1),r4head(19),r4head(20),r4head(21),r4head(22),header2(2)

      if (AXISYM) then
        if (is_on_the_axis(ispec_selected_rec(irec))) then ! TODO verify ispec_selected_rec and not ispec_selected_source
          call lagrange_any(xi_receiver(irec),NGLJ,xiglj,hxir,hpxir)
        else
          call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
        endif
      else
        call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
      endif

      call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)

      if (P_SV) then
        ! P_SV-case
        source_adjointe(irec_local,:,1) = adj_src_s(:,1)
        source_adjointe(irec_local,:,2) = adj_src_s(:,3)
      else
        ! SH-case
        source_adjointe(irec_local,:,1) = adj_src_s(:,2)
      endif

      if (.not. GPU_MODE) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            if (P_SV) then
              ! P_SV-case
              adj_sourcearray(:,1,i,j) = hxir(i) * hgammar(j) * adj_src_s(:,1)
              adj_sourcearray(:,2,i,j) = hxir(i) * hgammar(j) * adj_src_s(:,3)
            else
              ! SH-case
              adj_sourcearray(:,1,i,j) = hxir(i) * hgammar(j) * adj_src_s(:,2)
            endif
          enddo
        enddo
        adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
      endif

    endif !  if (myrank == which_proc_receiver(irec))
  enddo ! irec

  ! closes files
  if (seismotype == 4) then
    close(111)
  else
    close(111)
    close(112)
    close(113)
  endif

  ! frees memory
  deallocate(adj_src_s)
  deallocate(adj_sourcearray)

  end subroutine add_adjoint_sources_SU

!
!-----------------------------------------------------------------------------------------
!

  subroutine setup_source_receiver_tangentials()

! tangential computation

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

  ! local parameters
  integer :: i,irec,i_source
  integer :: ier,nrec_alloc
  integer :: irecloc
  double precision :: x_final_receiver_dummy, z_final_receiver_dummy

  ! for Lagrange interpolants
  double precision, external :: hgll, hglj

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
      if (which_proc_receiver(irec) == myrank) then
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
      if (is_proc_source(i_source) == 1) then
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
        if (myrank == 0 .and. is_proc_source(i_source) == 1 .and. nb_proc_source(i_source) == 1) then
          source_courbe_eros(i_source) = n1_tangential_detection_curve
          anglesource_recv = anglesource(i_source)
#ifdef USE_MPI
        else if (myrank == 0) then
          do i = 1, nb_proc_source(i_source) - is_proc_source(i_source)
            call MPI_recv(source_courbe_eros(i_source),1,MPI_INTEGER, &
                        MPI_ANY_SOURCE,42,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
            call MPI_recv(anglesource_recv,1,MPI_DOUBLE_PRECISION, &
                        MPI_ANY_SOURCE,43,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
          enddo
        else if (is_proc_source(i_source) == 1) then
          call MPI_send(n1_tangential_detection_curve,1,MPI_INTEGER,0,42,MPI_COMM_WORLD,ier)
          call MPI_send(anglesource(i_source),1,MPI_DOUBLE_PRECISION,0,43,MPI_COMM_WORLD,ier)
#endif
        endif

#ifdef USE_MPI
        call MPI_bcast(anglesource_recv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
        anglesource(i_source) = anglesource_recv
#endif
      endif !  if (is_proc_source(i_source) == 1)
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
        if (which_proc_receiver(irec) == myrank) then
          irecloc = irecloc + 1
          n1_tangential_detection_curve = rec_tangential_detection_curve(irecloc)
          x_final_receiver_dummy = x_final_receiver(irec)
          z_final_receiver_dummy = z_final_receiver(irec)
#ifdef USE_MPI
        else

          call MPI_RECV(n1_tangential_detection_curve,1,MPI_INTEGER,&
             which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
          call MPI_RECV(x_final_receiver_dummy,1,MPI_DOUBLE_PRECISION,&
             which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
          call MPI_RECV(z_final_receiver_dummy,1,MPI_DOUBLE_PRECISION,&
             which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

#endif
        endif

#ifdef USE_MPI
      else
        if (which_proc_receiver(irec) == myrank) then
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

  subroutine setup_source_receiver_interpolation()

  use constants,only: NGLLX,NGLLZ,NGLJ

  use specfem_par,only: myrank,nrec,nrecloc,NSOURCES, &
    ispec_selected_rec,ispec_selected_source,which_proc_receiver,is_proc_source, &
    xigll,zigll, &
    hxir_store,hgammar_store,xi_receiver,gamma_receiver,hxir,hpxir,hgammar,hpgammar, &
    hxis_store,hgammas_store,xi_source,gamma_source,hxis,hpxis,hgammas,hpgammas, &
    AXISYM,is_on_the_axis,xiglj

  use specfem_par_gpu,only: xir_store_loc,gammar_store_loc

  implicit none

  ! local parameters
  integer :: i,irec,irec_local,ier

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
      if (is_on_the_axis(ispec_selected_rec(irec)) .and. myrank == which_proc_receiver(irec)) then
        call lagrange_any(xi_receiver(irec),NGLJ,xiglj,hxir,hpxir)
        !do j = 1,NGLJ ! AB AB Same result with that loop
        !  hxir(j) = hglj(j-1,xi_receiver(irec),xiglj,NGLJ)
        !enddo
      else
        call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
      endif
    else
      call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
    endif
    call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)

    ! stores Lagrangians for receiver
    hxir_store(irec,:) = hxir(:)
    hgammar_store(irec,:) = hgammar(:)

    ! local receivers in this slice
    if (myrank == which_proc_receiver(irec)) then
      irec_local = irec_local + 1
      xir_store_loc(irec_local,:) = hxir(:)
      gammar_store_loc(irec_local,:) = hgammar(:)
    endif
  enddo

  ! allocate Lagrange interpolators for sources
  allocate(hxis_store(NSOURCES,NGLLX), &
           hgammas_store(NSOURCES,NGLLZ),stat=ier)
  if (ier /= 0) stop 'Error allocating source h**_store arrays'

  ! define and store Lagrange interpolators at all the sources
  do i = 1,NSOURCES
    if (AXISYM) then
      if (is_on_the_axis(ispec_selected_source(i)) .and. is_proc_source(i) == 1) then
        call lagrange_any(xi_source(i),NGLJ,xiglj,hxis,hpxis)
        !do j = 1,NGLJ ! AB AB same result with that loop
        !  hxis(j) = hglj(j-1,xi_source(i),xiglj,NGLJ)
        !enddo
      else
        call lagrange_any(xi_source(i),NGLLX,xigll,hxis,hpxis)
      endif
    else
      call lagrange_any(xi_source(i),NGLLX,xigll,hxis,hpxis)
    endif
    call lagrange_any(gamma_source(i),NGLLZ,zigll,hgammas,hpgammas)

    ! stores Lagrangians for source
    hxis_store(i,:) = hxis(:)
    hgammas_store(i,:) = hgammas(:)
  enddo

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_source_receiver_interpolation

