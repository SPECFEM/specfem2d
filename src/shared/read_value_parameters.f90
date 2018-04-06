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
    if (ios /= 0) call stop_the_code('error while reading input file')

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
    if (index_equal_sign <= 1 .or. index_equal_sign == len_trim(string_read)) call stop_the_code( &
'incorrect syntax detected in DATA/Par_file')
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
  if (ierr /= 0) call stop_the_code('error reading material parameter')

  !print *,trim(string_read)

  read(string_read,*,iostat=ierr) i,icodematread,val0read,val1read,val2read,val3read,val4read,val5read, &
                      val6read,val7read,val8read,val9read,val10read,val11read,val12read

  if (ierr /= 0) call stop_the_code('error reading material parameters line')

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
  if (ierr /= 0) call stop_the_code('error reading region coordinates')

  !print *,string_read

  read(string_read,*,iostat=ierr) value_to_read_1,value_to_read_2,value_to_read_3,value_to_read_4,value_to_read_5

  if (ierr /= 0) call stop_the_code('error reading region coordinates line')

  end subroutine read_region_coordinates_p

!--------------------

subroutine open_parameter_file_from_master_only()

  use constants, only: MAX_STRING_LEN,IN_DATA_FILES

  implicit none

  character(len=MAX_STRING_LEN) :: filename_main,filename_run0001
  logical :: exists_main_Par_file,exists_run0001_Par_file
  integer :: ier

  filename_main = IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'Par_file'

! also see if we are running several independent runs in parallel
! to do so, add the right directory for that run for the master process only here
  filename_run0001 = 'run0001/'//filename_main(1:len_trim(filename_main))
  call param_open(filename_main, len(filename_main), ier)
  if (ier == 0) then
    exists_main_Par_file = .true.
    call close_parameter_file()
  else
    exists_main_Par_file    = .false.
  endif
  call param_open(filename_run0001, len(filename_run0001), ier)
  if (ier == 0) then
    exists_run0001_Par_file = .true.
    call close_parameter_file()
  else
    exists_run0001_Par_file = .false.
  endif

  !if (exists_main_Par_file .and. exists_run0001_Par_file) then ! TODO why is it like that in the 3D version??
  !  print *
  !  print *,'cannot have both DATA/Par_file and run0001/DATA/Par_file present, please remove one of them'
  !  stop 'error: two different copies of the Par_file'
  !endif

  call param_open(filename_main, len(filename_main), ier)
  if (ier /= 0) then
    call param_open(filename_run0001, len(filename_run0001), ier)
    if (ier /= 0) then
      print *
      print *,'opening file failed, please check your file path and run-directory.'
      call stop_the_code('error opening Par_file')
    endif
  endif

  end subroutine open_parameter_file_from_master_only

!--------------------


  subroutine open_parameter_file()

  use constants, only: MAX_STRING_LEN,mygroup,IN_DATA_FILES
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

  integer ierr
  common /param_err_common/ ierr
  character(len=MAX_STRING_LEN) filename,path_to_add

  ! to use c routines
  filename = trim(IN_DATA_FILES)//'Par_file'

  ! see if we are running several independent runs in parallel
  ! if so, add the right directory for that run
  ! (group numbers start at zero, but directory names start at run0001, thus we add one)
  ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    filename = path_to_add(1:len_trim(path_to_add))//filename(1:len_trim(filename))
  endif

  call param_open(filename, len_trim(filename), ierr)
  if (ierr /= 0) then
    print *
    print *,'opening file failed, please check your file path and run-directory.'
    call stop_the_code('error opening Par_file')
  endif

  end subroutine open_parameter_file

!--------------------

  subroutine close_parameter_file()

  implicit none

  ! to use C routines
  call param_close()

  end subroutine close_parameter_file

!--------------------

  integer function err_occurred()

  implicit none

  integer ierr
  common /param_err_common/ ierr

  err_occurred = ierr

  end function err_occurred

!--------------------
!--------------------
!--------------------

  subroutine dummy_routine()

! dummy routine that does nothing, it is there just to fix an Intel ifort compiler bug
! with some releases of that compiler, in file src/specfem2D/locate_receivers.F90

  implicit none

  integer :: i,j,k

  i = 12
  j = 14
  k = i + j

  end subroutine dummy_routine

