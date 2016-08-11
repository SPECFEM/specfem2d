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

! read values from parameter file, ignoring white lines and comments

  subroutine read_value_integer(iin,ignore_junk,value_to_read)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer iin
  logical ignore_junk
  integer value_to_read
  character(len=MAX_STRING_LEN) string_read

  call read_next_line(iin,ignore_junk,string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_integer

!--------------------

  subroutine read_value_double_precision(iin,ignore_junk,value_to_read)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer iin
  logical ignore_junk
  double precision value_to_read
  character(len=MAX_STRING_LEN) string_read

  call read_next_line(iin,ignore_junk,string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_double_precision

!--------------------

  subroutine read_value_logical(iin,ignore_junk,value_to_read)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer iin
  logical ignore_junk
  logical value_to_read
  character(len=MAX_STRING_LEN) string_read

  call read_next_line(iin,ignore_junk,string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_logical

!--------------------

  subroutine read_value_string(iin,ignore_junk,value_to_read)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer iin
  logical ignore_junk
  character(len=MAX_STRING_LEN) value_to_read
  character(len=MAX_STRING_LEN) string_read

  call read_next_line(iin,ignore_junk,string_read)
  read(string_read,'(a)') value_to_read

  end subroutine read_value_string

!--------------------

  subroutine read_two_interface_points(iin,ignore_junk,value_to_read_1,value_to_read_2)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer iin
  logical ignore_junk
  double precision value_to_read_1,value_to_read_2
  character(len=MAX_STRING_LEN) string_read

  call read_next_line(iin,ignore_junk,string_read)
  read(string_read,*) value_to_read_1,value_to_read_2

  end subroutine read_two_interface_points

!--------------------

  subroutine read_next_line(iin,ignore_junk,string_read)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical ignore_junk
  character(len=MAX_STRING_LEN) string_read

  integer ios,iin,index_equal_sign

  do
    ! daniel: actually MAX_STRING_LEN set to 512...
    read(unit=iin,fmt="(a256)",iostat=ios) string_read
    if (ios /= 0) stop 'error while reading input file'

! suppress leading white spaces, if any
    string_read = adjustl(string_read)

! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
    if (index(string_read,achar(13)) > 0) string_read = string_read(1:index(string_read,achar(13))-1)

! exit loop when we find the first line that is not a comment or a white line
    if (len_trim(string_read) == 0) cycle
    if (string_read(1:1) /= '#') exit

  enddo

! suppress trailing white spaces, if any
  string_read = string_read(1:len_trim(string_read))

! suppress trailing comments, if any
  if (index(string_read,'#') > 0) string_read = string_read(1:index(string_read,'#')-1)

! suppress leading junk (up to the first equal sign, included) if needed
  if (ignore_junk) then
    index_equal_sign = index(string_read,'=')
    if (index_equal_sign <= 1 .or. index_equal_sign == len_trim(string_read)) stop 'incorrect syntax detected in DATA/Par_file'
    string_read = string_read(index_equal_sign + 1:len_trim(string_read))
  endif

! suppress leading and trailing white spaces again, if any, after having suppressed the leading junk
  string_read = adjustl(string_read)
  string_read = string_read(1:len_trim(string_read))

  end subroutine read_next_line

!--------------------




!--------------------
!--------------------
! uses param_reader.c functions
!--------------------
!--------------------


  subroutine read_value_integer_p(value_to_read, name)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer value_to_read
  character(len=*) name
  character(len=MAX_STRING_LEN) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  read(string_read,*) value_to_read

  end subroutine read_value_integer_p

!--------------------

  subroutine read_value_double_precision_p(value_to_read, name)

  use constants, only: MAX_STRING_LEN

  implicit none

  double precision value_to_read
  character(len=*) name
  character(len=MAX_STRING_LEN) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  read(string_read,*) value_to_read

  end subroutine read_value_double_precision_p

!--------------------

  subroutine read_value_logical_p(value_to_read, name)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical value_to_read
  character(len=*) name
  character(len=MAX_STRING_LEN) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  read(string_read,*) value_to_read

  end subroutine read_value_logical_p

!--------------------

  subroutine read_value_string_p(value_to_read, name)

  use constants, only: MAX_STRING_LEN

  implicit none

  character(len=*) value_to_read
  character(len=*) name
  character(len=MAX_STRING_LEN) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  value_to_read = string_read

  end subroutine read_value_string_p

!--------------------

  subroutine read_value_integer_next_p(value_to_read, name)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer value_to_read
  character(len=*) name
  character(len=MAX_STRING_LEN) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read_nextparam(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  read(string_read,*) value_to_read

  end subroutine read_value_integer_next_p

!--------------------

  subroutine read_value_double_prec_next_p(value_to_read, name)

  use constants, only: MAX_STRING_LEN

  implicit none

  double precision value_to_read
  character(len=*) name
  character(len=MAX_STRING_LEN) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read_nextparam(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  read(string_read,*) value_to_read

  end subroutine read_value_double_prec_next_p

!--------------------

  subroutine read_value_logical_next_p(value_to_read, name)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical value_to_read
  character(len=*) name
  character(len=MAX_STRING_LEN) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read_nextparam(string_read, len(string_read), name, len(name), ierr)
  if (ierr /= 0) return
  read(string_read,*) value_to_read

  end subroutine read_value_logical_next_p


!--------------------

  subroutine read_material_parameters_p(i,icodematread,val0read,val1read,val2read,val3read, &
                         val4read,val5read,val6read,val7read,val8read,val9read,val10read, &
                         val11read,val12read)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer i,icodematread
  double precision val0read,val1read,val2read,val3read,val4read,val5read,val6read,val7read, &
                   val8read,val9read,val10read,val11read,val12read

  character(len=MAX_STRING_LEN) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read_nextline(string_read, len(string_read), ierr)
  if (ierr /= 0) stop 'error reading material parameter'

  !print *,trim(string_read)

  read(string_read,*,iostat=ierr) i,icodematread,val0read,val1read,val2read,val3read,val4read,val5read, &
                      val6read,val7read,val8read,val9read,val10read,val11read,val12read

  if (ierr /= 0) stop 'error reading material parameters line'

  end subroutine read_material_parameters_p

!--------------------

  subroutine read_region_coordinates_p(value_to_read_1,value_to_read_2, &
                          value_to_read_3,value_to_read_4,value_to_read_5)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer value_to_read_1,value_to_read_2,value_to_read_3,value_to_read_4,value_to_read_5
  character(len=MAX_STRING_LEN) string_read
  integer ierr
  common /param_err_common/ ierr

  call param_read_nextline(string_read, len(string_read), ierr)
  if (ierr /= 0) stop 'error reading region coordinates'

  !print *,string_read

  read(string_read,*,iostat=ierr) value_to_read_1,value_to_read_2,value_to_read_3,value_to_read_4,value_to_read_5

  if (ierr /= 0) stop 'error reading region coordinates line'

  end subroutine read_region_coordinates_p


!--------------------


  subroutine open_parameter_file()

  use constants, only: MAX_STRING_LEN

  implicit none
  integer ierr
  common /param_err_common/ ierr
  character(len=MAX_STRING_LEN) filename

  ! to use c routines
  filename = 'DATA/Par_file'

  call param_open(filename, len_trim(filename), ierr)
  if (ierr /= 0 ) stop 'error opening DATA/Par_file file'

  end subroutine open_parameter_file

!--------------------

  subroutine close_parameter_file()

  implicit none

  ! to use c routines
  call param_close()

  end subroutine close_parameter_file

!--------------------

  integer function err_occurred()

  implicit none

  integer ierr
  common /param_err_common/ ierr

  err_occurred = ierr

  end function err_occurred

