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
  subroutine write_wavefield_dumps()

  use specfem_par, only: this_is_the_first_time_we_dump,myrank,nglob,nspec,mask_ibool, &
                         ibool,coord,imagetype_wavefield_dumps,p_sv,it,SIMULATION_TYPE, &
                         use_binary_for_wavefield_dumps,wavefield_file, &
                         vector_field_display,potential_acoustic,potential_gravitoacoustic, &
                         potential_gravito,displ_elastic,displs_poroelastic, &
                         potential_dot_acoustic,veloc_elastic,velocs_poroelastic,  &
                         potential_dot_dot_acoustic,accel_elastic,accels_poroelastic

  implicit none
  include "constants.h"

  !local variables
  integer :: i,j,ispec,iglob,icounter,nb_of_values_to_save

  if( myrank == 0 ) then
    write(IOUT,*)
    write(IOUT,*) 'Dumping the wave field to a file for time step ',it
  endif

  if( this_is_the_first_time_we_dump ) then

    if (.not. allocated(mask_ibool)) allocate(mask_ibool(nglob))

! save the grid separately once and for all
    if( use_binary_for_wavefield_dumps ) then
      write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.bin')") myrank
      open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='unknown', &
           action='write',recl=2*SIZE_REAL)
    else
      write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.txt')") myrank
      open(unit=27,file=wavefield_file,status='unknown',action='write')
    endif

    icounter = 0
    mask_ibool(:) = .false.
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
           iglob = ibool(i,j,ispec)
           if( .not. mask_ibool(iglob) ) then
             icounter = icounter + 1
             mask_ibool(iglob) = .true.
             if( use_binary_for_wavefield_dumps ) then
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
    open(unit=27,file=wavefield_file,status='unknown',action='write')
    write(27,*) icounter
    close(27)
    if( icounter /= nglob) stop 'error: should have icounter == nglob in wavefield dumps'

    this_is_the_first_time_we_dump = .false.

  endif

  if( imagetype_wavefield_dumps == 1 ) then
    if( myrank == 0 ) write(IOUT,*) 'dumping the displacement vector...'
      call compute_vector_whole_medium(potential_acoustic,potential_gravitoacoustic, &
                                       potential_gravito,displ_elastic,displs_poroelastic)
  else if( imagetype_wavefield_dumps == 2 ) then
    if( myrank == 0 ) write(IOUT,*) 'dumping the velocity vector...'
    call compute_vector_whole_medium(potential_dot_acoustic,potential_gravitoacoustic, &
                                     potential_gravito,veloc_elastic,velocs_poroelastic)

  else if( imagetype_wavefield_dumps == 3 ) then
    if( myrank == 0 ) write(IOUT,*) 'dumping the acceleration vector...'
    call compute_vector_whole_medium(potential_dot_dot_acoustic,potential_gravitoacoustic, &
                                     potential_gravito,accel_elastic,accels_poroelastic)

  else if( imagetype_wavefield_dumps == 4 .and. p_sv ) then
    if( myrank == 0 ) write(IOUT,*) 'dumping the pressure field...'
    call compute_pressure_whole_medium()

  else if( imagetype_wavefield_dumps == 4 .and. .not. p_sv ) then
    call exit_MPI('cannot dump the pressure field for SH (membrane) waves')

  else
    call exit_MPI('wrong type of flag for wavefield dumping')
  endif

  if( use_binary_for_wavefield_dumps ) then
    if( p_sv .and. .not. imagetype_wavefield_dumps == 4 ) then
      nb_of_values_to_save = 2
    else
      nb_of_values_to_save = 1
    endif
    write(wavefield_file,"('OUTPUT_FILES/wavefield',i7.7,'_',i2.2,'_',i3.3,'.bin')") it,SIMULATION_TYPE,myrank
    open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='unknown', &
             action='write',recl=nb_of_values_to_save*SIZE_REAL)
  else
    write(wavefield_file,"('OUTPUT_FILES/wavefield',i7.7,'_',i2.2,'_',i3.3,'.txt')") it,SIMULATION_TYPE,myrank
    open(unit=27,file=wavefield_file,status='unknown',action='write')
  endif

  icounter = 0
  mask_ibool(:) = .false.
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        if( .not. mask_ibool(iglob) ) then
          icounter = icounter + 1
          mask_ibool(iglob) = .true.
          if( use_binary_for_wavefield_dumps ) then
            if( p_sv .and. .not. imagetype_wavefield_dumps == 4 ) then
              write(27,rec=icounter) sngl(vector_field_display(1,iglob)),sngl(vector_field_display(3,iglob))
            else if( p_sv .and. imagetype_wavefield_dumps == 4 ) then
              ! by convention we use the third component of the array to store the pressure above
              write(27,rec=icounter) sngl(vector_field_display(3,iglob))
            else ! SH case
              write(27,rec=icounter) sngl(vector_field_display(2,iglob))
            endif
          else
            if( p_sv .and. .not. imagetype_wavefield_dumps == 4 ) then
              write(27,*) sngl(vector_field_display(1,iglob)),sngl(vector_field_display(3,iglob))
            else if( p_sv .and. imagetype_wavefield_dumps == 4 ) then
              ! by convention we use the third component of the array to store the pressure above
              write(27,*) sngl(vector_field_display(3,iglob))
            else ! SH case
              write(27,*) sngl(vector_field_display(2,iglob))
            endif
          endif
        endif
      enddo
    enddo
  enddo

  close(27)
  if( myrank ==0 ) write(IOUT,*) 'Wave field dumped'

  end subroutine write_wavefield_dumps

