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


  subroutine save_forward_arrays()

  use constants, only: IMAIN,OUTPUT_FILES

  use specfem_par

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname,outputname2

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "Saving forward arrays:"
    call flush_IMAIN()
  endif

  ! absorbing boundaries
  if (STACEY_ABSORBING_CONDITIONS) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'Saving Stacey absorbing boundary contributions...'
      call flush_IMAIN()
    endif

    if (any_acoustic) then
      !--- left absorbing boundary
      if (nspec_left > 0) then
        ! opens file
        write(outputname,'(a,i6.6,a)') 'absorb_acoustic_left',myrank,'.bin'
        open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        ! writes boundary contributions
        write(35) b_absorb_acoustic_left
        close(35)
      endif
      !--- right absorbing boundary
      if (nspec_right > 0) then
        write(outputname,'(a,i6.6,a)') 'absorb_acoustic_right',myrank,'.bin'
        open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        ! writes boundary contributions
        write(35) b_absorb_acoustic_right
        close(35)
      endif
      !--- bottom absorbing boundary
      if (nspec_bottom > 0) then
        write(outputname,'(a,i6.6,a)') 'absorb_acoustic_bottom',myrank,'.bin'
        open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        ! writes boundary contributions
        write(35) b_absorb_acoustic_bottom
        close(35)
      endif
      !--- top absorbing boundary
      if (nspec_top > 0) then
        write(outputname,'(a,i6.6,a)') 'absorb_acoustic_top',myrank,'.bin'
        open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        ! writes boundary contributions
        write(35) b_absorb_acoustic_top
        close(35)
      endif
    endif !any acoustic

    if (any_elastic) then
      !--- left absorbing boundary
      if (nspec_left > 0) then
        write(outputname,'(a,i6.6,a)') 'absorb_elastic_left',myrank,'.bin'
        open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        ! writes boundary contributions
        write(35) b_absorb_elastic_left
        close(35)
      endif
      !--- right absorbing boundary
      if (nspec_right > 0) then
        write(outputname,'(a,i6.6,a)') 'absorb_elastic_right',myrank,'.bin'
        open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        ! writes boundary contributions
        write(35) b_absorb_elastic_right
        close(35)
      endif
      !--- bottom absorbing boundary
      if (nspec_bottom > 0) then
        write(outputname,'(a,i6.6,a)') 'absorb_elastic_bottom',myrank,'.bin'
        open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        ! writes boundary contributions
        write(35) b_absorb_elastic_bottom
        close(35)
      endif
      !--- top absorbing boundary
      if (nspec_top > 0) then
        write(outputname,'(a,i6.6,a)') 'absorb_elastic_top',myrank,'.bin'
        open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        ! writes boundary contributions
        write(35) b_absorb_elastic_top
        close(35)
      endif
    endif !any elastic

    if (any_poroelastic) then
      !--- left absorbing boundary
      if (nspec_left > 0) then
        write(outputname,'(a,i6.6,a)') 'absorb_poro_s_left',myrank,'.bin'
        write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_left',myrank,'.bin'
        open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        open(unit=36,file=trim(OUTPUT_FILES)//outputname2,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        ! writes boundary contributions
        write(35) b_absorb_poro_s_left
        write(36) b_absorb_poro_w_left
        close(35)
        close(36)
      endif
      !--- right absorbing boundary
      if (nspec_right > 0) then
        write(outputname,'(a,i6.6,a)') 'absorb_poro_s_right',myrank,'.bin'
        write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_right',myrank,'.bin'
        open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        open(unit=36,file=trim(OUTPUT_FILES)//outputname2,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        ! writes boundary contributions
        write(35) b_absorb_poro_s_right
        write(36) b_absorb_poro_w_right
        close(35)
        close(36)
      endif
      !--- bottom absorbing boundary
      if (nspec_bottom > 0) then
        write(outputname,'(a,i6.6,a)') 'absorb_poro_s_bottom',myrank,'.bin'
        write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_bottom',myrank,'.bin'
        open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        open(unit=36,file=trim(OUTPUT_FILES)//outputname2,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        ! writes boundary contributions
        write(35) b_absorb_poro_s_bottom
        write(36) b_absorb_poro_w_bottom
        close(35)
        close(36)
      endif
      !--- top absorbing boundary
      if (nspec_top > 0) then
        write(outputname,'(a,i6.6,a)') 'absorb_poro_s_top',myrank,'.bin'
        write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_top',myrank,'.bin'
        open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        open(unit=36,file=trim(OUTPUT_FILES)//outputname2,status='unknown',form='unformatted',iostat=ier)
        if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
        ! writes boundary contributions
        write(35) b_absorb_poro_s_top
        write(36) b_absorb_poro_w_top
        close(35)
        close(36)
      endif
    endif ! poroelastic
  endif ! STACEY_ABSORBING_CONDITIONS

  ! PML
  if (anyabs_glob .and. PML_BOUNDARY_CONDITIONS) then
    if (any_elastic .and. nglob_interface > 0) close(71)
    if (any_acoustic .and. nglob_interface > 0) close(72)
  endif

  ! save last frame
  if (any_elastic) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'Saving elastic last frame...'
      call flush_IMAIN()
    endif

    write(outputname,'(a,i6.6,a)') 'lastframe_elastic',myrank,'.bin'
    open(unit=55,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_elastic**.bin')

    write(55) displ_elastic
    write(55) veloc_elastic
    write(55) accel_elastic
    close(55)
  endif ! elastic

  if (any_poroelastic) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'Saving poroelastic last frame...'
      call flush_IMAIN()
    endif

    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_s',myrank,'.bin'
    open(unit=55,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_poroelastic_s**.bin')

    write(55) displs_poroelastic
    write(55) velocs_poroelastic
    write(55) accels_poroelastic
    close(55)

    write(outputname,'(a,i6.6,a)') 'lastframe_poroelastic_w',myrank,'.bin'
    open(unit=56,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_poroelastic_w**.bin')

    write(56) displw_poroelastic
    write(56) velocw_poroelastic
    write(56) accelw_poroelastic
    close(56)
  endif ! poroelastic

  if (any_acoustic) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'Saving acoustic last frame...'
      call flush_IMAIN()
    endif

    write(outputname,'(a,i6.6,a)') 'lastframe_acoustic',myrank,'.bin'
    open(unit=55,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening file lastframe_acoustic**.bin')

    write(55) potential_acoustic
    write(55) potential_dot_acoustic
    write(55) potential_dot_dot_acoustic
    close(55)
  endif ! acoustic

  end subroutine save_forward_arrays
