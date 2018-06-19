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

  subroutine read_source_file(NSOURCES)

  ! reads in source file DATA/SOURCE

  use constants, only: IMAIN,IGNORE_JUNK,NLINES_PER_SOURCE,TINYVAL,PI,mygroup,IN_DATA_FILES
  use source_file_par
  use shared_parameters, only: DT

  implicit none

  integer,intent(in) :: NSOURCES

  ! local parameters
  integer :: ier,icounter,i_source,num_sources
  character(len=256) string_read
  character(len=MAX_STRING_LEN) :: source_filename,path_to_add
  integer, parameter :: IIN_SOURCE = 22

  ! allocates memory arrays
  allocate(source_surf(NSOURCES), &
           xs(NSOURCES), &
           zs(NSOURCES), &
           source_type(NSOURCES), &
           time_function_type(NSOURCES), &
           name_of_source_file(NSOURCES), &
           burst_band_width(NSOURCES), &
           f0_source(NSOURCES), &
           tshift_src(NSOURCES), &
           anglesource(NSOURCES), &
           Mxx(NSOURCES), &
           Mxz(NSOURCES), &
           Mzz(NSOURCES), &
           factor(NSOURCES),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating source arrays')

  ! initializes
  xs(:) = 0.d0
  zs(:) = 0.d0

  source_type(:) = 0
  time_function_type(:) = 0

  f0_source(:) = 0.d0
  tshift_src(:) = 0.d0
  anglesource(:) = 0.d0
  Mxx(:) = 0.d0
  Mxz(:) = 0.d0
  Mzz(:) = 0.d0

  source_filename = trim(IN_DATA_FILES)//'SOURCE'

  ! mygroup has been initialized with negative value. It is positive just in the case NUMBER_OF_SIMULTANEOUS_RUNS > 1
  if (mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    source_filename = path_to_add(1:len_trim(path_to_add))//source_filename(1:len_trim(source_filename))
  endif

  ! counts lines
  open(unit=IIN_SOURCE,file=trim(source_filename),status='old',action='read',iostat=ier)
  if (ier /= 0) call stop_the_code('Error opening source file, please make sure file exists...')

  ! counts number of lines
  icounter = 0
  do while(ier == 0)
    read(IIN_SOURCE,"(a)",iostat=ier) string_read

    if (ier == 0) then
      ! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
      if (index(string_read,achar(13)) > 0) string_read = string_read(1:index(string_read,achar(13))-1)

      ! suppress leading and trailing white spaces, if any
      string_read = adjustl(string_read)
      string_read = string_read(1:len_trim(string_read))

      ! if the line is not empty and is not a comment, count it
      if (len_trim(string_read) > 0 .and. (index(string_read,'#') == 0 .or. index(string_read,'#') > 1)) icounter = icounter + 1
    endif
  enddo
  close(IIN_SOURCE)

  ! checks counter
  if (mod(icounter,NLINES_PER_SOURCE) /= 0) &
    call stop_the_code('total number of non blank and non comment lines in SOURCE file should be a multiple of NLINES_PER_SOURCE')

  ! total number of sources
  num_sources = icounter / NLINES_PER_SOURCE

  ! checks number of sources
  if (num_sources < 1) call stop_the_code('need at least one source in SOURCE file')
  if (num_sources /= NSOURCES) then
       print *,'Error invalid num_sources :',num_sources
       print *,'NSOURCES :',NSOURCES
       call stop_the_code( &
'Error: Total number of sources in DATA/SOURCE is different from that declared in the Par_file, please check...')
  endif

  ! reads in source parameters
  open(unit=IIN_SOURCE,file=trim(source_filename),status='old',action='read',iostat=ier)
  if (ier /= 0) call stop_the_code('Error opening source file, please make sure file exists...')

  ! reads in all source informations
  do  i_source= 1,NSOURCES

    ! source set to surface
    call read_value_logical(IIN_SOURCE,IGNORE_JUNK,source_surf(i_source))

    ! x/z location
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,xs(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,zs(i_source))

    ! source and source time function type
    call read_value_integer(IIN_SOURCE,IGNORE_JUNK,source_type(i_source))
    call read_value_integer(IIN_SOURCE,IGNORE_JUNK,time_function_type(i_source))

    ! external source time function file (sft type == 8)
    name_of_source_file(i_source) = ''
    call read_value_string(IIN_SOURCE,IGNORE_JUNK,name_of_source_file(i_source))

    ! burst (stf type == 9)
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,burst_band_width(i_source))

    ! dominant frequency
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,f0_source(i_source))

    ! time shift
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,tshift_src(i_source))

    ! force source angle
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,anglesource(i_source))

    ! moment tensor
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,Mxx(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,Mzz(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,Mxz(i_source))

    ! amplification factor
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,factor(i_source))


    ! Dirac/Heaviside
    ! if Dirac source time function, use a very thin Gaussian instead
    ! if Heaviside source time function, use a very thin error function instead
    if (time_function_type(i_source) == 4 .or. time_function_type(i_source) == 5) then
      f0_source(i_source) = 1.d0 / (10.d0 * DT)
    endif

    ! checks source frequency
    if (abs(f0_source(i_source)) < TINYVAL) then
      call exit_MPI(0,'Error source frequency is zero')
    endif

    ! convert angle from degrees to radians
    anglesource(i_source) = anglesource(i_source) * PI / 180.d0

    ! user output
    ! note: we will further process source info in solver,
    !         here we just read in the given specifics and show them
    write(IMAIN,*) 'Source', i_source
    write(IMAIN,*) '  Position xs, zs = ',xs(i_source),zs(i_source)
    write(IMAIN,*)

    ! source type
    write(IMAIN,*) '  Source type (1=force, 2=moment tensor, 3=Rayleigh wave, 4=plane incident P, &
                   &5=plane incident S): ',source_type(i_source)
    select case (source_type(i_source))
    case (1)
      ! force
      write(IMAIN,*) '  Force source:'
      write(IMAIN,*) '  Angle of the source (deg) = ',anglesource(i_source)
      write(IMAIN,*) '  Multiplying factor  = ',factor(i_source)
    case (2)
      ! moment tensor
      write(IMAIN,*) '  Moment tensor source:'
      write(IMAIN,*) '  Mxx of the source = ',Mxx(i_source)
      write(IMAIN,*) '  Mzz of the source = ',Mzz(i_source)
      write(IMAIN,*) '  Mxz of the source = ',Mxz(i_source)
    case (3)
      ! Rayleigh wave
      write(IMAIN,*) '  Rayleigh wave source:'
    case (4)
      ! plane P wave without converted/refracted phases
      write(IMAIN,*) '  Plane P-wave source without converted/refracted phases:'
      write(IMAIN,*) '  Angle of the incident wave (deg) = ',anglesource(i_source)
    case (5)
      ! plane S wave without converted/refracted phases
      write(IMAIN,*) '  Plane S-wave source without converted/refracted phases:'
      write(IMAIN,*) '  Angle of the incident wave (deg) = ',anglesource(i_source)
    case default
      ! not supported yet
      call stop_the_code('Error invalid source type! must be 1, 2, 3, 4 or 5, exiting...')
    end select
    write(IMAIN,*)

    ! STF
    write(IMAIN,*) '  Time function type (1=Ricker, 2=First derivative, 3=Gaussian, 4=Dirac, 5=Heaviside, &
                   &8=Read from file, 9=burst):',time_function_type(i_source)
    select case (time_function_type(i_source))
    case (1)
      write(IMAIN,*) '  Ricker wavelet (second-derivative):'
      write(IMAIN,*) '  Frequency, delay = ',f0_source(i_source),tshift_src(i_source)
    case (2)
      write(IMAIN,*) '  Ricker wavelet (first-derivative):'
      write(IMAIN,*) '  Frequency, delay = ',f0_source(i_source),tshift_src(i_source)
    case (3)
      write(IMAIN,*) '  Gaussian:'
      write(IMAIN,*) '  Frequency, delay = ',f0_source(i_source),tshift_src(i_source)
    case (4)
      write(IMAIN,*) '  Dirac:'
      write(IMAIN,*) '  Frequency, delay = ',f0_source(i_source),tshift_src(i_source)
    case (5)
      write(IMAIN,*) '  Heaviside:'
      write(IMAIN,*) '  Frequency, delay = ',f0_source(i_source),tshift_src(i_source)
    case (6)
      write(IMAIN,*) '  Ocean acoustics (type I):'
      write(IMAIN,*) '  Frequency, delay = ',f0_source(i_source),tshift_src(i_source)
    case (7)
      write(IMAIN,*) '  Ocean acoustics (type II):'
      write(IMAIN,*) '  Frequency, delay = ',f0_source(i_source),tshift_src(i_source)
    case (8)
      write(IMAIN,*) '  External source time function file:'
      write(IMAIN,*) '  Source read from file:',trim(name_of_source_file(i_source))
    case (9)
      write(IMAIN,*) '  Burst wavelet:'
      write(IMAIN,*) '  Burst band width: ',burst_band_width(i_source)
    case (10)
      write(IMAIN,*) '  Sinus source time function:'
      write(IMAIN,*) '  Frequency, delay = ',f0_source(i_source),tshift_src(i_source)
    case (11)
      write(IMAIN,*) '  Ormsby source time function:'
      write(IMAIN,*) '  Frequency, delay = ',f0_source(i_source),tshift_src(i_source)
    case default
      call stop_the_code('Error invalid source time function type! must be between 1 and 9, exiting...')
    end select
    write(IMAIN,*)

  enddo ! do i_source= 1,NSOURCES
  close(IIN_SOURCE)

  end subroutine read_source_file

