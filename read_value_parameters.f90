
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) January 2005
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
  read(string_read,"(a)") value_to_read

  end subroutine read_value_string

!--------------------

  subroutine read_next_line(iin,ignore_junk,string_read)

  implicit none

  logical ignore_junk
  character(len=100) string_read

  integer ios,iin,index_equal_sign

  do
    read(unit=iin,fmt=200,iostat=ios) string_read
    if(ios /= 0) stop 'error while reading input file'

! suppress leading white spaces, if any
    string_read = adjustl(string_read)

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

! format
 200 format(a100)

  end subroutine read_next_line

