
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!                         Dimitri Komatitsch
!                     University of Pau, France
!
!                          (c) April 2007
!
!========================================================================

! read values from parameter file, ignoring white lines and comments

  subroutine read_value_integer(iin,ignore_junk,value_to_read)

  implicit none

  integer iin
  logical ignore_junk
  integer value_to_read
  character(len=100) string_read

  call read_next_line(iin,ignore_junk,string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_integer

!--------------------

  subroutine read_value_double_precision(iin,ignore_junk,value_to_read)

  implicit none

  integer iin
  logical ignore_junk
  double precision value_to_read
  character(len=100) string_read

  call read_next_line(iin,ignore_junk,string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_double_precision

!--------------------

  subroutine read_value_logical(iin,ignore_junk,value_to_read)

  implicit none

  integer iin
  logical ignore_junk
  logical value_to_read
  character(len=100) string_read

  call read_next_line(iin,ignore_junk,string_read)
  read(string_read,*) value_to_read

  end subroutine read_value_logical

!--------------------

  subroutine read_value_string(iin,ignore_junk,value_to_read)

  implicit none

  integer iin
  logical ignore_junk
  character(len=*) value_to_read
  character(len=100) string_read

  call read_next_line(iin,ignore_junk,string_read)
  value_to_read = string_read

  end subroutine read_value_string

!--------------------

  subroutine read_two_interface_points(iin,ignore_junk,value_to_read_1,value_to_read_2)

  implicit none

  integer iin
  logical ignore_junk
  double precision value_to_read_1,value_to_read_2
  character(len=100) string_read

  call read_next_line(iin,ignore_junk,string_read)
  read(string_read,*) value_to_read_1,value_to_read_2

  end subroutine read_two_interface_points

!--------------------

  subroutine read_region_coordinates(iin,ignore_junk,value_to_read_1,value_to_read_2, &
                          value_to_read_3,value_to_read_4,value_to_read_5)

  implicit none

  integer iin
  logical ignore_junk
  integer value_to_read_1,value_to_read_2,value_to_read_3,value_to_read_4,value_to_read_5
  character(len=100) string_read

  call read_next_line(iin,ignore_junk,string_read)
  read(string_read,*) value_to_read_1,value_to_read_2,value_to_read_3,value_to_read_4,value_to_read_5

  end subroutine read_region_coordinates

!--------------------

  subroutine read_material_parameters(iin,ignore_junk,i,icodematread,rhoread,cpread,csread,aniso3read,aniso4read)

  implicit none

  integer iin
  logical ignore_junk
  integer i,icodematread
  double precision rhoread,cpread,csread,aniso3read,aniso4read
  character(len=100) string_read

  call read_next_line(iin,ignore_junk,string_read)
  read(string_read,*) i,icodematread,rhoread,cpread,csread,aniso3read,aniso4read

  end subroutine read_material_parameters

!--------------------

  subroutine read_next_line(iin,ignore_junk,string_read)

  implicit none

  logical ignore_junk
  character(len=100) string_read

  integer ios,iin,index_equal_sign

  do
    read(unit=iin,fmt="(a100)",iostat=ios) string_read
    if(ios /= 0) stop 'error while reading input file'

! suppress leading white spaces, if any
    string_read = adjustl(string_read)

! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
    if(index(string_read,achar(13)) > 0) string_read = string_read(1:index(string_read,achar(13))-1)

! exit loop when we find the first line that is not a comment or a white line
    if(len_trim(string_read) == 0) cycle
    if(string_read(1:1) /= '#') exit

  enddo

! suppress trailing white spaces, if any
  string_read = string_read(1:len_trim(string_read))

! suppress trailing comments, if any
  if(index(string_read,'#') > 0) string_read = string_read(1:index(string_read,'#')-1)

! suppress leading junk (up to the first equal sign, included) if needed
  if(ignore_junk) then
    index_equal_sign = index(string_read,'=')
    if(index_equal_sign <= 1 .or. index_equal_sign == len_trim(string_read)) stop 'incorrect syntax detected in DATA/Par_file'
    string_read = string_read(index_equal_sign + 1:len_trim(string_read))
  endif

! suppress leading and trailing white spaces again, if any, after having suppressed the leading junk
  string_read = adjustl(string_read)
  string_read = string_read(1:len_trim(string_read))

  end subroutine read_next_line

