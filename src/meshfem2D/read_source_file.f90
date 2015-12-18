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


  subroutine read_source_file(NSOURCES)

  ! reads in source file DATA/SOURCE

  use source_file_par

  implicit none

  include "constants.h"

  integer,intent(in) :: NSOURCES

  ! local parameters
  integer :: ios,icounter,i_source,num_sources
  character(len=150) string_read
  integer, parameter :: IIN_SOURCE = 22

  ! allocates memory arrays
  allocate(source_surf(NSOURCES))
  allocate(xs(NSOURCES))
  allocate(zs(NSOURCES))
  allocate(source_type(NSOURCES))
  allocate(time_function_type(NSOURCES))
  allocate(name_of_source_file(NSOURCES))
  allocate(burst_band_width(NSOURCES))
  allocate(f0_source(NSOURCES))
  allocate(tshift_src(NSOURCES))
  allocate(anglesource(NSOURCES))
  allocate(Mxx(NSOURCES))
  allocate(Mxz(NSOURCES))
  allocate(Mzz(NSOURCES))
  allocate(factor(NSOURCES))

  ! counts lines
  open(unit=IIN_SOURCE,file='DATA/SOURCE',iostat=ios,status='old',action='read')
  if (ios /= 0) stop 'error opening DATA/SOURCE file'

  icounter = 0
  do while(ios == 0)
     read(IIN_SOURCE,"(a)",iostat=ios) string_read

     if (ios == 0) then

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
    stop 'total number of non blank and non comment lines in SOURCE file should be a multiple of NLINES_PER_SOURCE'

  ! total number of sources
  num_sources = icounter / NLINES_PER_SOURCE

  if (num_sources < 1) stop 'need at least one source in SOURCE file'
  if (num_sources /= NSOURCES) then
       print *,'num_sources :',num_sources
       print *,'NSOURCES :',NSOURCES
       stop 'error: Total number of sources read is different from that declared in the Par_file'
  endif

  ! reads in source parameters
  open(unit=IIN_SOURCE,file='DATA/SOURCE',status='old',action='read')
  do  i_source= 1,NSOURCES
    call read_value_logical(IIN_SOURCE,IGNORE_JUNK,source_surf(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,xs(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,zs(i_source))
    call read_value_integer(IIN_SOURCE,IGNORE_JUNK,source_type(i_source))
    call read_value_integer(IIN_SOURCE,IGNORE_JUNK,time_function_type(i_source))
    name_of_source_file(i_source)=''
    call read_value_string(IIN_SOURCE,IGNORE_JUNK,name_of_source_file(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,burst_band_width(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,f0_source(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,tshift_src(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,anglesource(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,Mxx(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,Mzz(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,Mxz(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,factor(i_source))

    ! note: we will further process source info in solver,
    !         here we just read in the given specifics and show them

    print *
    print *,'Source', i_source
    print *,'Position xs, zs = ',xs(i_source),zs(i_source)
    print *,'Source type (1=force, 2=explosion): ',source_type(i_source)
    print *,'Angle of the source if force = ',anglesource(i_source)
    print *,'Multiplying factor = ',factor(i_source)
    if (time_function_type(i_source) == 8) then
      print *,"Source read from file:", name_of_source_file(i_source)
    else
      if (time_function_type(i_source) == 9) then
        print *,'Burst wavelet'
        print *,'Burst band width: ',burst_band_width(i_source)
      else
        print *,'Frequency, delay = ',f0_source(i_source),tshift_src(i_source)
        print *,'Time function type (1=Ricker, 2=First derivative, 3=Gaussian, 4=Dirac, 5=Heaviside, 8=Read from file, 9=burst):'&
               ,time_function_type(i_source)
        print *,'Mxx of the source if moment tensor = ',Mxx(i_source)
        print *,'Mzz of the source if moment tensor = ',Mzz(i_source)
        print *,'Mxz of the source if moment tensor = ',Mxz(i_source)
      endif
    endif

  enddo ! do i_source= 1,NSOURCES
  close(IIN_SOURCE)

  end subroutine read_source_file

