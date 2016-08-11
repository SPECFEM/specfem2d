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

  subroutine write_wavefield_dumps()

  use constants, only: IMAIN,SIZE_REAL,NGLLX,NGLLZ

  use specfem_par, only: myrank,nglob,nspec, &
                         ibool,coord,P_SV,it,SIMULATION_TYPE, &
                         potential_acoustic,potential_gravitoacoustic, &
                         potential_gravito,displ_elastic,displs_poroelastic, &
                         potential_dot_acoustic,veloc_elastic,velocs_poroelastic, &
                         potential_dot_dot_acoustic,accel_elastic,accels_poroelastic

  use specfem_par_movie, only: this_is_the_first_time_we_dump,mask_ibool,imagetype_wavefield_dumps, &
    use_binary_for_wavefield_dumps,vector_field_display

  implicit none

  !local variables
  integer :: i,j,ispec,iglob,icounter,nb_of_values_to_save
  integer :: ier
  ! name of wavefield snapshot file
  character(len=150) :: wavefield_file

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Dumping the wave field to a file for time step ',it
    call flush_IMAIN()
  endif

  if (this_is_the_first_time_we_dump) then

    if (.not. allocated(mask_ibool)) allocate(mask_ibool(nglob))

! save the grid separately once and for all
    if (use_binary_for_wavefield_dumps) then
      write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.bin')") myrank
      open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='unknown', &
           action='write',recl=2*SIZE_REAL,iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield_grid_for_dumps_**.bin')
    else
      write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.txt')") myrank
      open(unit=27,file=wavefield_file,status='unknown',action='write',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield_grid_for_dumps_**.txt')
    endif

    icounter = 0
    mask_ibool(:) = .false.
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
           iglob = ibool(i,j,ispec)
           if (.not. mask_ibool(iglob)) then
             icounter = icounter + 1
             mask_ibool(iglob) = .true.
             if (use_binary_for_wavefield_dumps) then
               write(27,rec=icounter) sngl(coord(1,iglob)),sngl(coord(2,iglob))
             else
               write(27,'(2e16.6)') coord(1,iglob),coord(2,iglob)
             endif
           endif
        enddo
      enddo
    enddo

    close(27)

    ! save nglob to a file once and for all
    write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_value_of_nglob_',i3.3,'.txt')") myrank
    open(unit=27,file=wavefield_file,status='unknown',action='write',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield_grid_value_of_nglob_**.txt')
    write(27,*) icounter
    close(27)
    if (icounter /= nglob) stop 'error: should have icounter == nglob in wavefield dumps'

    this_is_the_first_time_we_dump = .false.

  endif

  if (imagetype_wavefield_dumps == 1) then
    if (myrank == 0) write(IMAIN,*) 'dumping the displacement vector...'
      call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                                       potential_gravito,displ_elastic,displs_poroelastic)
  else if (imagetype_wavefield_dumps == 2) then
    if (myrank == 0) write(IMAIN,*) 'dumping the velocity vector...'
    call compute_vector_whole_medium(potential_dot_acoustic,potential_gravitoacoustic, &
                                     potential_gravito,veloc_elastic,velocs_poroelastic)

  else if (imagetype_wavefield_dumps == 3) then
    if (myrank == 0) write(IMAIN,*) 'dumping the acceleration vector...'
    call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_gravitoacoustic, &
                                     potential_gravito,accel_elastic,accels_poroelastic)

  else if (imagetype_wavefield_dumps == 4 .and. P_SV) then
    if (myrank == 0) write(IMAIN,*) 'dumping the pressure field...'
    call compute_pressure_whole_medium()

  else if (imagetype_wavefield_dumps == 4 .and. .not. P_SV) then
    call exit_MPI(myrank,'cannot dump the pressure field for SH (membrane) waves')

  else
    call exit_MPI(myrank,'wrong type of flag for wavefield dumping')
  endif

  if (use_binary_for_wavefield_dumps) then
    if (P_SV .and. .not. imagetype_wavefield_dumps == 4) then
      nb_of_values_to_save = 2
    else
      nb_of_values_to_save = 1
    endif
    write(wavefield_file,"('OUTPUT_FILES/wavefield',i7.7,'_',i2.2,'_',i3.3,'.bin')") it,SIMULATION_TYPE,myrank
    open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='unknown', &
             action='write',recl=nb_of_values_to_save*SIZE_REAL,iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield**.bin')
  else
    write(wavefield_file,"('OUTPUT_FILES/wavefield',i7.7,'_',i2.2,'_',i3.3,'.txt')") it,SIMULATION_TYPE,myrank
    open(unit=27,file=wavefield_file,status='unknown',action='write',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file wavefield**.txt')
  endif

  icounter = 0
  mask_ibool(:) = .false.
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        if (.not. mask_ibool(iglob)) then
          icounter = icounter + 1
          mask_ibool(iglob) = .true.
          if (use_binary_for_wavefield_dumps) then
            if (P_SV) then
              ! P-SV waves
              if (imagetype_wavefield_dumps == 4) then
                ! by convention we use the 2. component of the array to store the pressure above
                write(27,rec=icounter) sngl(vector_field_display(2,iglob))
              else
                write(27,rec=icounter) sngl(vector_field_display(1,iglob)),sngl(vector_field_display(2,iglob))
              endif
            else
              ! SH case
              write(27,rec=icounter) sngl(vector_field_display(1,iglob))
            endif
          else
            if (P_SV) then
              ! P-SV waves
              if (imagetype_wavefield_dumps == 4) then
                ! by convention we use the 2. component of the array to store the pressure above
                write(27,*) sngl(vector_field_display(2,iglob))
              else
                write(27,*) sngl(vector_field_display(1,iglob)),sngl(vector_field_display(2,iglob))
              endif
            else
              ! SH case
              write(27,*) sngl(vector_field_display(1,iglob))
            endif
          endif
        endif
      enddo
    enddo
  enddo

  close(27)

  if (myrank == 0) then
    write(IMAIN,*) 'Wave field dumped'
    call flush_IMAIN()
  endif

  end subroutine write_wavefield_dumps

