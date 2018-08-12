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

  subroutine read_external_model()

! reads in external model files

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,TINYVAL,IMAIN,IN_DATA_FILES,IIN,MAX_STRING_LEN

  use specfem_par, only: nspec,nglob,ibool,ispec_is_elastic,ispec_is_anisotropic, &
    coord,kmato,MODEL,myrank,setup_with_binary_database

  ! external model parameters
  use specfem_par, only: rhoext,vpext,vsext, &
    QKappa_attenuationext,Qmu_attenuationext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext, &
    ATTENUATION_VISCOELASTIC,ATTENUATION_VISCOACOUSTIC

  implicit none

  ! Local variables
  integer :: i,j,ispec,itmp
  integer :: ier
  real(kind=CUSTOM_REAL) :: tmp1,tmp2
  double precision :: vs_val,vp_val,rho_val
  character(len=MAX_STRING_LEN) :: inputname, line
  logical :: read_next_line

  ! note: we read in external models once the basic mesh with its geometry and GLL points has been setup.
  !       External models define new velocity/material parameters which need to be defined on all GLL points.
  if (myrank == 0) write(IMAIN,*) '  model: ',trim(MODEL)

  select case (trim(MODEL))
  case ('legacy')
    ! old model format
    write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_model_velocity.dat_input'
    if (myrank == 0) write(IMAIN,*) '  reading external files: ','DATA/proc*****_model_velocity.dat_input'

    open(unit=IIN,file=inputname,status='old',action='read',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_model_velocity.dat_input file.')
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! reads next data line
          read_next_line = .true.
          do while (read_next_line)
            ! format: #unused #unused #unused #rho #vp #vs
            read(IIN,'(a256)',iostat=ier) line
            if (ier /= 0) call stop_the_code('Error reading file model_velocity.dat_input')

            ! left adjust
            line = adjustl(line)

            ! trim white space
            line = trim(line)

            ! skip comment lines
            if (line(1:1) == '#') then
              read_next_line = .true.
            else
              read_next_line = .false.
            endif
          enddo

          ! reads in values
          read(line,*) itmp,tmp1,tmp2,rho_val,vp_val,vs_val
          rhoext(i,j,ispec) = rho_val
          vpext(i,j,ispec) = vp_val
          vsext(i,j,ispec) = vs_val
        enddo
      enddo
    enddo
    close(IIN)

    ! default no attenuation
    QKappa_attenuationext(:,:,:) = 9999.d0
    Qmu_attenuationext(:,:,:) = 9999.d0

  case ('ascii')
    ! ascii model format
    write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_rho_vp_vs.dat'
    if (myrank == 0) write(IMAIN,*) '  reading external files: ','DATA/proc*****_rho_vp_vs.dat'

    open(unit=IIN,file=inputname,status='old',action='read',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_rho_vp_vs.dat file.')
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! format: #unused #unused #rho #vp #vs
          read(IIN,*) tmp1,tmp2,rho_val,vp_val,vs_val
          rhoext(i,j,ispec) = rho_val
          vpext(i,j,ispec) = vp_val
          vsext(i,j,ispec) = vs_val
        enddo
      enddo
    enddo
    close(IIN)
    ! default no attenuation
    QKappa_attenuationext(:,:,:) = 9999.d0
    Qmu_attenuationext(:,:,:) = 9999.d0

  case ('binary','gll')
    ! binary formats
    write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_rho.bin'
    if (myrank == 0) write(IMAIN,*) '  reading external files: ','DATA/proc*****_rho.bin, .._vp.bin, .._vs.bin'

    open(unit = IIN, file = inputname, status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_rho.bin file.')

    read(IIN) rhoext
    close(IIN)
    print *, 'rho', minval(rhoext), maxval(rhoext)

    write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_vp.bin'
    open(unit = IIN, file = inputname, status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_vp.bin file.')

    read(IIN) vpext
    close(IIN)

    write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_vs.bin'
    open(unit = IIN, file = inputname, status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_vs.bin file.')

    read(IIN) vsext
    close(IIN)

    if (ATTENUATION_VISCOACOUSTIC) then
      write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_Qkappa.bin'
      open(unit = IIN, file = inputname, status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_Qkappa.bin file.')

      read(IIN) QKappa_attenuationext
      close(IIN)
      Qmu_attenuationext(:,:,:) = 9999.d0
    ! for the moment we don't do external model with both viscoacoustics and
    ! viscoelastics
    else if (ATTENUATION_VISCOELASTIC) then
      write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_Qkappa.bin'
      open(unit = IIN, file = inputname,status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_Qkappa.bin file.')

      read(IIN) QKappa_attenuationext
      close(IIN)

      write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_Qmu.bin'
      open(unit = IIN, file = inputname,status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_Qmu.bin file.')

      read(IIN) Qmu_attenuationext
      close(IIN)

    else
      ! default no attenuation
      QKappa_attenuationext(:,:,:) = 9999.d0
      Qmu_attenuationext(:,:,:) = 9999.d0
    endif

  case ('binary_voigt')
    ! Voigt model
    write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_rho.bin'
    if (myrank == 0) write(IMAIN,*) '  reading external files: ','DATA/proc*****_rho.bin, .._c11.bin, .._c55.bin'

    open(unit = IIN, file = inputname,status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_rho.bin file.')

    read(IIN) rhoext
    close(IIN)
    print *, 'rho', minval(rhoext), maxval(rhoext)

    write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_c11.bin'
    open(unit = IIN, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_c11.bin file.')

    read(IIN) c11ext
    close(IIN)
    print *, 'c11ext', minval(c11ext), maxval(c11ext)

    write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_c13.bin'
    open(unit = IIN, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_c13.bin file.')

    read(IIN) c13ext
    close(IIN)
    print *, 'c13ext', minval(c13ext), maxval(c13ext)

    write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_c15.bin'
    open(unit = IIN, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_c15.bin file.')

    read(IIN) c15ext
    close(IIN)
    print *, 'c15ext', minval(c15ext), maxval(c15ext)

    write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_c33.bin'
    open(unit = IIN, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_c33.bin file.')

    read(IIN) c33ext
    close(IIN)
    print *, 'c33ext', minval(c33ext), maxval(c33ext)

    write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_c35.bin'
    open(unit = IIN, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_c35.bin file.')

    read(IIN) c35ext
    close(IIN)
    print *, 'c35ext', minval(c35ext), maxval(c35ext)

    write(inputname,'(a,i6.6,a)') trim(IN_DATA_FILES)//'proc',myrank,'_c55.bin'
    open(unit = IIN, file = inputname,status='old',action='read',form='unformatted', iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening DATA/proc*****_c55.bin file.')

    read(IIN) c55ext
    close(IIN)
    print *, 'c55ext', minval(c55ext), maxval(c55ext)

    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if (c55ext(i,j,ispec) < TINYVAL) then
             c33ext(i,j,ispec) = 0.0_CUSTOM_REAL
             c12ext(i,j,ispec) = 0.0_CUSTOM_REAL
             c23ext(i,j,ispec) = 0.0_CUSTOM_REAL
             vpext(i,j,ispec) = 1500.0
             vsext(i,j,ispec) = 0.0_CUSTOM_REAL
          else
             c12ext(i,j,ispec) = 1.d-6
             c23ext(i,j,ispec) = 1.d-6
             vpext(i,j,ispec) = sqrt(c33ext(i,j,ispec)/rhoext(i,j,ispec))
             vsext(i,j,ispec) = sqrt(c55ext(i,j,ispec)/rhoext(i,j,ispec))
          endif
          c25ext(i,j,ispec) = 0.0_CUSTOM_REAL
        enddo
      enddo
    enddo

    ! default no attenuation
    QKappa_attenuationext(:,:,:) = 9999.d0
    Qmu_attenuationext(:,:,:) = 9999.d0

  case ('external')
    ! generic model defined in external files
    call define_external_model(coord,kmato,ibool,rhoext,vpext,vsext, &
                               QKappa_attenuationext,Qmu_attenuationext, &
                               c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,nspec,nglob)
  case ('marmousi')
    ! marmousi type model
    call define_external_model_from_marmousi(coord,ibool,rhoext,vpext,vsext, &
                                             QKappa_attenuationext,Qmu_attenuationext, &
                                             c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,nspec,nglob)

  case ('tomo')
    ! tomographic file
    call define_external_model_from_tomo_file()

  case default
    print *,"Error: unrecognized model = ",trim(MODEL)
    print *,"Invalid MODEL chosen, please check your Par_file settings..."
    call stop_the_code('Invalid MODEL parameter')

  end select

  ! check that the external model that has just been defined makes sense
  if (trim(MODEL) == 'external' .or. trim(MODEL) == 'tomo' .or. trim(MODEL) == 'marmousi') then
    ! checks velocities for each element
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! anisotropic elastic
          if (c11ext(i,j,ispec) > TINYVAL .or. c13ext(i,j,ispec) > TINYVAL .or. c15ext(i,j,ispec) > TINYVAL .or. &
              c33ext(i,j,ispec) > TINYVAL .or. c35ext(i,j,ispec) > TINYVAL .or. c55ext(i,j,ispec) > TINYVAL) then
            ! vp, vs : assign dummy values, trick to avoid floating point errors in the case of an anisotropic medium
            vpext(i,j,ispec) = 20.d0
            vsext(i,j,ispec) = 10.d0
          endif

          ! check that the element type is not redefined compared to what is defined initially in DATA/Par_file
          if ((c11ext(i,j,ispec) > TINYVAL .or. c13ext(i,j,ispec) > TINYVAL .or. c15ext(i,j,ispec) > TINYVAL .or. &
               c33ext(i,j,ispec) > TINYVAL .or. c35ext(i,j,ispec) > TINYVAL .or. c55ext(i,j,ispec) > TINYVAL) &
              .and. .not. ispec_is_anisotropic(ispec)) &
            call stop_the_code('Error: non anisotropic material in DATA/Par_file or &
                 &external mesh redefined as anisotropic in define_external_model()')

          ! acoustic element
          if (vsext(i,j,ispec) < TINYVAL .and. (ispec_is_elastic(ispec) .or. ispec_is_anisotropic(ispec))) &
            call stop_the_code('Error: non acoustic material in DATA/Par_file or &
                 &external mesh redefined as acoustic in define_external_model()')

          ! elastic element
          if (vsext(i,j,ispec) > TINYVAL .and. .not. ispec_is_elastic(ispec)) &
            call stop_the_code('Error: acoustic material in DATA/Par_file or &
                 &external mesh redefined as non acoustic in define_external_model()')
        enddo
      enddo
    enddo
  endif

  ! resets domain flags
  if (setup_with_binary_database /= 2) then
    call get_simulation_domains_from_external_models()
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_external_model

