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


  subroutine save_model_files()

  use constants, only: FOUR_THIRDS,TWO_THIRDS,IMAIN,IN_DATA_FILES

  use specfem_par

  implicit none

  ! local parameters
  integer :: i,ispec,j,iglob
  integer :: ier
  real(kind=4),dimension(:,:,:),allocatable :: rho_save, vp_save, vs_save, kappa_save, x_save, z_save, Qkappa_save,Qmu_save
  double precision :: rhol
  character(len=MAX_STRING_LEN) :: filename

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Saving model files to directory: ',trim(IN_DATA_FILES)
    write(IMAIN,*) '  SAVE_MODEL type: ',trim(SAVE_MODEL)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! allocates temporary arrays for file storage
  allocate(rho_save(NGLLX,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 01')
  rho_save(:,:,:) = 0.0

  allocate(vp_save(NGLLX,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 02')
  vp_save(:,:,:) = 0.0

  allocate(vs_save(NGLLX,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 03')
  vs_save(:,:,:) = 0.0

  allocate(kappa_save(NGLLX,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 04')
  kappa_save(:,:,:) = 0.0

  allocate(x_save(NGLLX,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 05')
  x_save(:,:,:) = 0.0

  allocate(z_save(NGLLX,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 06')
  z_save(:,:,:) = 0.0

  if (ATTENUATION_VISCOACOUSTIC) then
    allocate(Qkappa_save(NGLLX,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 07')
    Qkappa_save(:,:,:) = 0.0
  endif

  if (ATTENUATION_VISCOELASTIC) then
    if (ATTENUATION_VISCOACOUSTIC) &
      call exit_MPI(myrank,'Not possible yet to save model with both acoustic and elastic attenuation')

    allocate(Qkappa_save(NGLLX,NGLLZ,nspec),Qmu_save(NGLLX,NGLLZ,nspec),stat=ier)
    if (ier /= 0) call exit_MPI(myrank, 'error allocating save model arrays 08')
    Qkappa_save(:,:,:) = 0.0
  endif

  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        rhol = rhostore(i,j,ispec)
        vp_save(i,j,ispec) = rho_vpstore(i,j,ispec) / rhol
        vs_save(i,j,ispec) = rho_vsstore(i,j,ispec) / rhol
        rho_save(i,j,ispec) = rhol

        iglob = ibool(i,j,ispec)
        x_save(i,j,ispec) = coord(1,iglob)
        z_save(i,j,ispec) = coord(2,iglob)

        ! attenuation arrays
        if (ATTENUATION_VISCOACOUSTIC) then
          Qkappa_save(i,j,ispec) = qkappa_attenuation_store(i,j,ispec)
        endif

        if (ATTENUATION_VISCOELASTIC) then
          Qkappa_save(i,j,ispec) = qkappa_attenuation_store(i,j,ispec)
          Qmu_save(i,j,ispec) = qmu_attenuation_store(i,j,ispec)
        endif
      enddo
    enddo
  enddo

  ! SMNSR For compatibility with NUMBER_OF_SIMULTANEOUS_RUNS we have to change the lines trim(IN_DATA_FILES)//'proc'

  ! outputs model files
  if (trim(SAVE_MODEL) == 'legacy') then
    ! legacy format
    write(filename,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_model_velocity.dat_input'
    open(unit=1001,file=trim(filename),status='unknown',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_model_velocity.dat_input')
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          ! format: #iglob  #x  #z  #rho  #vp  #vs
          write(1001,'(I10,5e15.5e4)') iglob, x_save(i,j,ispec),z_save(i,j,ispec),rho_save(i,j,ispec), &
                                       vp_save(i,j,ispec),vs_save(i,j,ispec)
        enddo
      enddo
    enddo
    close(1001)
    if (myrank == 0) write(IMAIN,*) '  written to file(s): ',trim(filename)

  else if (trim(SAVE_MODEL) == 'ascii') then
    ! ascii format
    write(filename,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_rho_vp_vs.dat'
    open(unit=1001,file=trim(filename),status='unknown',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_rho_vp_vs.dat')
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! format: #x  #z  #rho  #vp  #vs
          write(1001,'(5e15.5e4)') x_save(i,j,ispec),z_save(i,j,ispec),rho_save(i,j,ispec), &
                                   vp_save(i,j,ispec),vs_save(i,j,ispec)
        enddo
      enddo
    enddo
    close(1001)
    if (myrank == 0) write(IMAIN,*) '  written to file(s): ',trim(filename)

  else if ((trim(SAVE_MODEL) == 'binary') .or. (trim(SAVE_MODEL) == 'gll')) then
    ! binary and GLL format
    write(filename,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_rho.bin'
    open(unit=172,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_rho.bin')
    write(172) rho_save
    close(172)
    if (myrank == 0) write(IMAIN,*) '  written to file(s): ',trim(filename)

    write(filename,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_vp.bin'
    open(unit=172,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_vp.bin')
    write(172) vp_save
    close(172)
    if (myrank == 0) write(IMAIN,*) '  written to file(s): ',trim(filename)

    write(filename,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_vs.bin'
    open(unit=172,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_vs.bin')
    write(172) vs_save
    close(172)
    if (myrank == 0) write(IMAIN,*) '  written to file(s): ',trim(filename)

    write(filename,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_x.bin'
    open(unit=172,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_x.bin')
    write(172) x_save
    close(172)
    if (myrank == 0) write(IMAIN,*) '  written to file(s): ',trim(filename)

    write(filename,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_z.bin'
    open(unit=172,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_z.bin')
    write(172) z_save
    close(172)
    if (myrank == 0) write(IMAIN,*) '  written to file(s): ',trim(filename)

    write(filename,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_jacobian.bin'
    open(unit=172,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_jacobian.bin')
    write(172) jacobian
    close(172)
    if (myrank == 0) write(IMAIN,*) '  written to file(s): ',trim(filename)

    if (ATTENUATION_VISCOACOUSTIC) then
      write(filename,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_Qkappa.bin'
      open(unit=172,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_Qkappa.bin')
      write(172) Qkappa_save
      close(172)
      if (myrank == 0) write(IMAIN,*) '  written to file(s): ',trim(filename)
    endif

    if (ATTENUATION_VISCOELASTIC) then
      write(filename,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_Qkappa.bin'
      open(unit=172,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_Qkappa.bin')
      write(172) Qkappa_save
      close(172)
      if (myrank == 0) write(IMAIN,*) '  written to file(s): ',trim(filename)

      write(filename,'(a,i6.6,a)')trim(IN_DATA_FILES)//'proc',myrank,'_Qmu.bin'
      open(unit=172,file=trim(filename),status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening model file proc**_Qmu.bin')
      write(172) Qmu_save
      close(172)
      if (myrank == 0) write(IMAIN,*) '  written to file(s): ',trim(filename)
    endif

  else
    call stop_the_code('Save Model not implemented for external and tomo')
  endif !Type of model

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  done'
    call flush_IMAIN()
  endif

  end subroutine save_model_files

