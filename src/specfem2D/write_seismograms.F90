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

! write seismograms to text files

  subroutine write_seismograms()

  use constants, only: ZERO,IMAIN
  use specfem_par
  use specfem_par_gpu, only: Mesh_pointer

  implicit none

  ! local parameters
  integer :: i, j, iglob, irecloc, ispec, i_sig
  integer :: seismotype_l
  double precision :: valux,valuz,valcurl

  ! vector field in an element
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: vector_field_element
  ! pressure in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: pressure_element
  ! curl in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: curl_element

  ! timing
  double precision, external :: wtime
  double precision :: write_time_begin,write_time

  ! flag to indicate that traces for kernel runs are taken from backward/reconstructed wavefields instead of adjoint wavefields;
  ! useful for debugging.
  ! default (.false.) is to output adjoint wavefield traces
  logical, parameter :: OUTPUT_BACKWARD_WAVEFIELD = .false.

  ! checks subsampling recurrence
  if (mod(it-1,NTSTEP_BETWEEN_OUTPUT_SAMPLE) == 0) then

    ! update position in seismograms
    seismo_current = seismo_current + 1

    ! check for edge effects
    if (seismo_current < 1 .or. seismo_current > nlength_seismogram) then
      print *,'Error: seismo_current ',seismo_current,' should be between 1 and ',nlength_seismogram
      call stop_the_code('Error: seismo_current out of bounds in recording of seismograms')
    endif

    do i_sig = 1,NSIGTYPE
      seismotype_l = seismotypeVec(i_sig)

      ! updates local receiver records
      if (nrecloc > 0) then

        ! computes seismogram entry for all local receivers
        if (.not. GPU_MODE) then
          ! on CPU
          do irecloc = 1,nrecloc

            ispec = ispec_selected_rec_loc(irecloc)

            ! initializes local element arrays
            vector_field_element(:,:,:) = 0._CUSTOM_REAL
            pressure_element(:,:) = 0._CUSTOM_REAL
            curl_element(:,:) = 0._CUSTOM_REAL

            ! compute seismo/pressure/curl in this element
            select case (seismotype_l)
            case (1)
              ! displacement
              if (OUTPUT_BACKWARD_WAVEFIELD .and. SIMULATION_TYPE == 3) then
                call compute_vector_one_element(b_potential_acoustic,b_displ_elastic,b_displs_poroelastic,ispec, &
                                                vector_field_element)
              else
                call compute_vector_one_element(potential_acoustic,displ_elastic,displs_poroelastic,ispec, &
                                                vector_field_element)
              endif
            case (2)
              ! velocity
              if (OUTPUT_BACKWARD_WAVEFIELD .and. SIMULATION_TYPE == 3) then
                call compute_vector_one_element(b_potential_dot_acoustic,b_veloc_elastic,b_velocs_poroelastic,ispec, &
                                                vector_field_element)
              else
                call compute_vector_one_element(potential_dot_acoustic,veloc_elastic,velocs_poroelastic,ispec, &
                                                vector_field_element)
              endif
            case (3)
              ! acceleration
              if (OUTPUT_BACKWARD_WAVEFIELD .and. SIMULATION_TYPE == 3) then
                call compute_vector_one_element(b_potential_dot_dot_acoustic,b_accel_elastic,b_accels_poroelastic,ispec, &
                                                vector_field_element)
              else
                call compute_vector_one_element(potential_dot_dot_acoustic,accel_elastic,accels_poroelastic,ispec, &
                                                vector_field_element)
              endif
            case (4)
              ! pressure
              if (OUTPUT_BACKWARD_WAVEFIELD .and. SIMULATION_TYPE == 3) then
                call compute_pressure_one_element(ispec,pressure_element, &
                                                  b_displ_elastic,b_displs_poroelastic,b_displw_poroelastic, &
                                                  b_potential_dot_dot_acoustic,b_potential_acoustic)
              else
                call compute_pressure_one_element(ispec,pressure_element, &
                                                  displ_elastic,displs_poroelastic,displw_poroelastic, &
                                                  potential_dot_dot_acoustic,potential_acoustic)
              endif
            case (5)
              ! displacement
              if (OUTPUT_BACKWARD_WAVEFIELD .and. SIMULATION_TYPE == 3) then
                call compute_vector_one_element(b_potential_acoustic,b_displ_elastic,b_displs_poroelastic,ispec, &
                                                vector_field_element)
                ! curl of displacement
                call compute_curl_one_element(b_displ_elastic,b_displs_poroelastic,ispec,curl_element)
              else
                call compute_vector_one_element(potential_acoustic,displ_elastic,displs_poroelastic,ispec, &
                                                vector_field_element)
                ! curl of displacement
                call compute_curl_one_element(displ_elastic,displs_poroelastic,ispec,curl_element)
              endif

            case (6)
              ! fluid potential
              ! uses pressure_element to store local element values
              if (ispec_is_acoustic(ispec)) then
                if (OUTPUT_BACKWARD_WAVEFIELD .and. SIMULATION_TYPE == 3) then
                  ! backward
                  do j = 1,NGLLZ
                    do i = 1,NGLLX
                      iglob = ibool(i,j,ispec)
                      pressure_element(i,j) = b_potential_acoustic(iglob)
                    enddo
                  enddo
                else
                  ! forward/backward wavefield
                  do j = 1,NGLLZ
                    do i = 1,NGLLX
                      iglob = ibool(i,j,ispec)
                      pressure_element(i,j) = potential_acoustic(iglob)
                    enddo
                  enddo
                endif
              endif

            case default
              call stop_the_code('Invalid seismotype for writing seismograms')
            end select

            ! perform the general interpolation using Lagrange polynomials
            call compute_interpolated_dva(irecloc,ispec,vector_field_element,pressure_element,curl_element, &
                                          valux,valuz,valcurl,seismotype_l)

            ! rotate seismogram components if needed, except if recording pressure, which is a scalar
            if (seismotype_l == 4 .or. seismotype_l == 6) then
              ! pressure/potential type has only single component
              sisux(seismo_current,irecloc,i_sig) = valux
              sisuz(seismo_current,irecloc,i_sig) = ZERO
            else
              if (P_SV) then
                sisux(seismo_current,irecloc,i_sig) =   cosrot_irec(irecloc)*valux + sinrot_irec(irecloc)*valuz
                sisuz(seismo_current,irecloc,i_sig) = - sinrot_irec(irecloc)*valux + cosrot_irec(irecloc)*valuz
              else
                sisux(seismo_current,irecloc,i_sig) = valux
                sisuz(seismo_current,irecloc,i_sig) = ZERO
              endif
              ! additional curl case
              if (seismotype_l == 5) siscurl(seismo_current,irecloc,i_sig) = valcurl
            endif
          enddo ! irecloc

        else
          ! on GPU
          ! Simulating seismograms
          call compute_seismograms_cuda(Mesh_pointer,i_sig,sisux(:,:,i_sig),sisuz(:,:,i_sig), &
                                        seismo_current,nlength_seismogram, &
                                        ELASTIC_SIMULATION,ACOUSTIC_SIMULATION,USE_TRICK_FOR_BETTER_PRESSURE, &
                                        ATTENUATION_VISCOELASTIC, &
                                        it,NSTEP)
          ! note: curl not implemented yet
        endif ! GPU_MODE
      endif ! nrecloc
    enddo  ! loop on signal types (seismotype)

  endif  ! NTSTEP_BETWEEN_OUTPUT_SAMPLE

  ! save temporary or final seismograms
  if (mod(it,NTSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! timing
    write_time_begin = wtime()

    ! Loop on signal types
    do i_sig = 1,NSIGTYPE
      seismotype_l = seismotypeVec(i_sig)

      call write_seismograms_to_file(sisux(:,:,i_sig),sisuz(:,:,i_sig),siscurl(:,:,i_sig),seismotype_l,seismo_current, &
                                     seismo_offset)

    enddo ! loop on signal types (seismotype)

    ! updates current seismogram offsets
    seismo_offset = seismo_offset + seismo_current
    seismo_current = 0

    ! user output
    if (myrank == 0) then
      ! timing
      write_time = wtime() - write_time_begin
      ! output
      write(IMAIN,*)
      write(IMAIN,*) 'Total number of time steps done: ', it-it_begin+1
      if (WRITE_SEISMOGRAMS_BY_MAIN) then
        write(IMAIN,*) 'Writing the seismograms by main proc alone took ',sngl(write_time),' seconds'
      else
        write(IMAIN,*) 'Writing the seismograms in parallel took ',sngl(write_time),' seconds'
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  end subroutine write_seismograms

!================================================================

  subroutine write_seismograms_to_file(sisux_l,sisuz_l,siscurl_l,seismotype_l,seismo_current_l,seismo_offset_l)

  use constants, only: NDIM,MAX_LENGTH_NETWORK_NAME,MAX_LENGTH_STATION_NAME,OUTPUT_FILES, RegInt_K

  use specfem_par, only: station_name,network_name,NSTEP,islice_selected_rec,nrec,myrank,DT,t0, &
                         NTSTEP_BETWEEN_OUTPUT_SAMPLE,nrecloc,nlength_seismogram, &
                         P_SV,SU_FORMAT,save_ASCII_seismograms, &
                         save_binary_seismograms_single,save_binary_seismograms_double,x_source,z_source, &
                         WRITE_SEISMOGRAMS_BY_MAIN

  implicit none

  double precision,dimension(nlength_seismogram,nrecloc),intent(in) :: sisux_l,sisuz_l,siscurl_l
  integer,intent(in) :: seismotype_l,seismo_current_l,seismo_offset_l

  ! local parameters
  logical :: save_binary_seismograms
  integer :: irec,length_station_name,length_network_name,iorientation,isample,number_of_components

  character(len=4) :: channel
  character(len=1) :: component
  character(len=4) :: suffix
  character(len=150) :: sisname,sismo_statut

  ! to write seismograms in single precision SEP and double precision binary format
  double precision, dimension(:,:,:), allocatable :: buffer_binary
  double precision :: time_t
  real, dimension(:), allocatable :: single_precision_seismo

  integer :: irecloc,ier

  integer (kind=RegInt_K) :: ioffset

  logical :: file_unit_12_has_been_opened,file_unit_13_has_been_opened,file_unit_14_has_been_opened
  logical :: file_unit_15_has_been_opened,file_unit_16_has_been_opened,file_unit_17_has_been_opened

  file_unit_12_has_been_opened = .false.
  file_unit_13_has_been_opened = .false.
  file_unit_14_has_been_opened = .false.
  file_unit_15_has_been_opened = .false.
  file_unit_16_has_been_opened = .false.
  file_unit_17_has_been_opened = .false.

  ! safety stop
  if (.not. WRITE_SEISMOGRAMS_BY_MAIN) &
    stop 'Error writing seismograms in parallel not supported yet!'

! write seismograms

  ! save displacement, velocity, acceleration or pressure
  if (seismotype_l == 1) then
    component = 'd'
  else if (seismotype_l == 2) then
    component = 'v'
  else if (seismotype_l == 3) then
    component = 'a'
  else if (seismotype_l == 4) then
    component = 'p'
  else if (seismotype_l == 6) then
    component = 'x'
  else if (seismotype_l == 5) then
    component = 'c'
  else
    call exit_MPI(myrank,'wrong component to save for seismograms')
  endif

  ! only one seismogram if pressures or SH (membrane) waves
  if (seismotype_l == 4 .or. seismotype_l == 6 .or. .not. P_SV) then
    number_of_components = 1
  else if (seismotype_l == 5) then
    ! adds curl
    number_of_components = NDIM + 1
  else
    number_of_components = NDIM
  endif

  allocate(buffer_binary(seismo_current_l,nrec,number_of_components),stat=ier)
  allocate(single_precision_seismo(seismo_current_l),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating array buffer_binary')
  buffer_binary(:,:,:) = 0.d0

  ! see if we need to save any seismogram in binary format
  save_binary_seismograms = ( save_binary_seismograms_single .or. save_binary_seismograms_double )

  ! binary file output
  if ( (save_binary_seismograms .or. SU_FORMAT) .and. myrank == 0) then
    ! filename extension
    if (.not. SU_FORMAT) then
      suffix = '.bin'
    else
      suffix = '.su'
    endif

    ! deletes old files if first pass
    if (seismo_offset_l == 0) then
      write(sismo_statut,'(a)') 'replace'
    else
      write(sismo_statut,'(a)') 'old'
    endif

    ! write the new files
    if (save_binary_seismograms_single .or. SU_FORMAT) then
      if (seismotype_l == 4 .or. seismotype_l == 6) then
        open(unit=12,file=trim(OUTPUT_FILES)//'Up_file_single_'//component//suffix,status=sismo_statut,access='stream')
        file_unit_12_has_been_opened = .true.
      else if (.not. P_SV) then
        open(unit=12,file=trim(OUTPUT_FILES)//'Uy_file_single_'//component//suffix,status=sismo_statut,access='stream')
        file_unit_12_has_been_opened = .true.
      else
        open(unit=12,file=trim(OUTPUT_FILES)//'Ux_file_single_'//component//suffix,status=sismo_statut,access='stream')
        file_unit_12_has_been_opened = .true.
      endif
    endif

    if (save_binary_seismograms_double) then
      if (seismotype_l == 4 .or. seismotype_l == 6) then
        open(unit=13,file=trim(OUTPUT_FILES)//'Up_file_double_'//component//suffix,status=sismo_statut,access='stream')
        file_unit_13_has_been_opened = .true.
      else if (.not. P_SV) then
        open(unit=13,file=trim(OUTPUT_FILES)//'Uy_file_double_'//component//suffix,status=sismo_statut,access='stream')
        file_unit_13_has_been_opened = .true.
      else
        open(unit=13,file=trim(OUTPUT_FILES)//'Ux_file_double_'//component//suffix,status=sismo_statut,access='stream')
        file_unit_13_has_been_opened = .true.
      endif
    endif

    ! no Z component seismogram if pressure
    if (seismotype_l /= 4 .and. seismotype_l /= 6 .and. P_SV) then
      if (save_binary_seismograms_single .or. SU_FORMAT) &
        open(unit=14,file=trim(OUTPUT_FILES)//'Uz_file_single_'//component//suffix,status=sismo_statut,access='stream')
        file_unit_14_has_been_opened = .true.
      if (save_binary_seismograms_double) &
        open(unit=15,file=trim(OUTPUT_FILES)//'Uz_file_double_'//component//suffix,status=sismo_statut,access='stream')
        file_unit_15_has_been_opened = .true.
    endif

    ! curl output
    if (seismotype_l == 5) then
      if (save_binary_seismograms_single .or. SU_FORMAT) &
        open(unit=16,file=trim(OUTPUT_FILES)//'Uc_file_single_'//component//suffix,status=sismo_statut,access='stream')
        file_unit_16_has_been_opened = .true.
      if (save_binary_seismograms_double) &
        open(unit=17,file=trim(OUTPUT_FILES)//'Uc_file_double_'//component//suffix,status=sismo_statut,access='stream')
        file_unit_17_has_been_opened = .true.
    endif
  endif ! save_binary_seismograms

  ! sets seismogram values
  irecloc = 0
  do irec = 1,nrec
    if (myrank == 0) then

      if (myrank == islice_selected_rec(irec)) then
        irecloc = irecloc + 1

        ! fills buffer
        buffer_binary(:,irec,1) = sisux_l(1:seismo_current_l,irecloc)
        if (number_of_components == 2) then
          buffer_binary(:,irec,2) = sisuz_l(1:seismo_current_l,irecloc)
        else if (number_of_components == 3) then
          ! adds curl trace
          buffer_binary(:,irec,2) = sisuz_l(1:seismo_current_l,irecloc)
          buffer_binary(:,irec,3) = siscurl_l(1:seismo_current_l,irecloc)
        endif

#ifdef WITH_MPI
      else
        ! gets from other processes
        if (WRITE_SEISMOGRAMS_BY_MAIN) then
          ! collects seismogram components on main
          call recv_dp(buffer_binary(1,irec,1), seismo_current_l, islice_selected_rec(irec), irec)

          if (number_of_components == 2) then
            call recv_dp(buffer_binary(1,irec,2), seismo_current_l, islice_selected_rec(irec), irec)
          endif

          if (number_of_components == 3) then
            call recv_dp(buffer_binary(1,irec,2), seismo_current_l, islice_selected_rec(irec), irec)
            call recv_dp(buffer_binary(1,irec,3), seismo_current_l, islice_selected_rec(irec), irec)
          endif
        endif

#endif
      endif ! islice_selected_rec

#ifdef WITH_MPI
    else
      ! secondary processes (myrank > 0)
      if (WRITE_SEISMOGRAMS_BY_MAIN) then
        ! sends seismogram values to main
        if (myrank == islice_selected_rec(irec)) then
          irecloc = irecloc + 1
          call send_dp(sisux_l(1,irecloc), seismo_current_l, 0, irec)

          if (number_of_components >= 2) then
            call send_dp(sisuz_l(1,irecloc), seismo_current_l, 0, irec)
          endif

          if (number_of_components == 3) then
            call send_dp(siscurl_l(1,irecloc), seismo_current_l, 0, irec)
          endif
        endif
      endif
#endif

    endif ! myrank == 0
  enddo

  ! outputs to file
  if (myrank == 0) then
    irecloc = 0
    do irec = 1,nrec

      if (.not. SU_FORMAT) then

        if (save_ASCII_seismograms) then

          ! write trace
          do iorientation = 1,number_of_components

            if (iorientation == 1) then
              channel = 'BXX'
            else if (iorientation == 2) then
              channel = 'BXZ'
            else if (iorientation == 3) then
              channel = 'cur'
            else
              call exit_MPI(myrank,'incorrect channel value')
            endif

            ! in case of pressure or potential, use different abbreviation
            if (seismotype_l == 4) channel = 'PRE'
            if (seismotype_l == 6) channel = 'POT'

            ! in case of SH (membrane) waves, use different abbreviation
            if (.not. P_SV) channel = 'BXY'

            ! create the name of the seismogram file for each slice
            ! file name includes the name of the station, the network and the component
            length_station_name = len_trim(station_name(irec))
            length_network_name = len_trim(network_name(irec))

            ! check that length conforms to standard
            if (length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) then
              call exit_MPI(myrank,'wrong length of station name')
            endif
            if (length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) then
              call exit_MPI(myrank,'wrong length of network name')
            endif

            write(sisname,"(a,a,'.',a,'.',a3,'.sem',a1)") &
            trim(OUTPUT_FILES),network_name(irec)(1:length_network_name),station_name(irec)(1:length_station_name),channel,component

            ! deletes old seismogram file when starting to write output
            if (seismo_offset_l == 0) then
              open(unit=11,file=sisname(1:len_trim(sisname)),status='unknown')
              close(11,status='delete')
            endif

            open(unit=11,file=sisname(1:len_trim(sisname)),status='unknown',position='append')

            ! make sure we never write more than the maximum number of time steps
            ! subtract offset of the source to make sure travel time is correct
            do isample = 1,seismo_current_l

              ! forward time
              time_t = dble(seismo_offset_l + isample - 1) * DT * NTSTEP_BETWEEN_OUTPUT_SAMPLE - t0

              write(11,*) time_t,' ',buffer_binary(isample,irec,iorientation)
            enddo

            close(11)
          enddo ! iorientation

        endif ! save_ASCII_seismograms

        ! write binary seismogram
        if (save_binary_seismograms) then

          ioffset = (irec-1) * NSTEP / NTSTEP_BETWEEN_OUTPUT_SAMPLE + seismo_offset_l

          if (save_binary_seismograms_single) then
            do isample = 1, seismo_current_l
              single_precision_seismo(isample) = sngl(buffer_binary(isample,irec,1))
            enddo
            write(12,pos=4*ioffset+1) single_precision_seismo
          endif
          if (save_binary_seismograms_double) write(13,pos=8*ioffset+1) buffer_binary(:,irec,1)

          if (seismotype_l /= 4 .and. seismotype_l /= 6 .and. P_SV) then
            if (save_binary_seismograms_single) then
              do isample = 1, seismo_current_l
                single_precision_seismo(isample) = sngl(buffer_binary(isample,irec,2))
              enddo
              write(14,pos=4*ioffset+1) single_precision_seismo
            endif
            if (save_binary_seismograms_double) write(15,pos=8*ioffset+1) buffer_binary(:,irec,2)
          endif

          if (seismotype_l == 5) then
            if (save_binary_seismograms_single) then
              do isample = 1, seismo_current_l
                single_precision_seismo(isample) = sngl(buffer_binary(isample,irec,3))
              enddo
              write(16,pos=4*ioffset+1) single_precision_seismo
            endif
            if (save_binary_seismograms_double) write(17,pos=8*ioffset+1) buffer_binary(:,irec,3)
          endif

        endif

      else
        ! if SU_FORMAT
        call write_output_SU(x_source(1),z_source(1),irec,buffer_binary,number_of_components, &
                             seismo_offset_l,seismo_current_l,seismotype_l)
      endif
    enddo
  endif ! myrank

  ! close files
  if (file_unit_12_has_been_opened) close(12)
  if (file_unit_13_has_been_opened) close(13)
  if (file_unit_14_has_been_opened) close(14)
  if (file_unit_15_has_been_opened) close(15)
  if (file_unit_16_has_been_opened) close(16)
  if (file_unit_17_has_been_opened) close(17)

  ! frees temporary buffer
  if (allocated(buffer_binary)) deallocate(buffer_binary)
  if (allocated(single_precision_seismo)) deallocate(single_precision_seismo)

  end subroutine write_seismograms_to_file

