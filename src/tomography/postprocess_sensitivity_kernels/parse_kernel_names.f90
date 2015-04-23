
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


subroutine parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)

  use postprocess_par,only: MAX_STRING_LEN, MAX_KERNEL_NAMES

  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: kernel_names_comma_delimited
  character(len=MAX_STRING_LEN),intent(inout) :: kernel_names(MAX_KERNEL_NAMES)
  integer,intent(out) :: nker

  ! local parameters
  integer :: iker
  character(len=MAX_STRING_LEN) :: tmp
  character,parameter :: delimiter = ','

  ! gets first name/token
  iker = 1
  call strtok(kernel_names_comma_delimited, delimiter, kernel_names(iker))

  ! null-string as first argument for successive strtok-calls
  tmp(1:1) = char(0)

  ! gets next names/tokens (strtok will return null-terminated token when finished)
  do while (kernel_names(iker)(1:1) /= char(0))
    ! increases name/token number
    iker = iker + 1
    ! gets next successive token (with null-terminated string as first argument)
    call strtok(tmp, delimiter, kernel_names(iker))
  enddo

  ! number of kernel names
  nker = iker-1

  ! checks name lengths (e.g. if kernel_name argument is "vsv,vsh," we will have a 3. kernel name with empty string)
  do iker = 1,nker
    if (len_trim(kernel_names(iker)) == 0) then
      print *,'Error encountered kernel name with zero length: kernel name number ',iker,' out of ',nker,' is empty'
      print *,'Please check your kernel_names argument...'
      stop 'Error kernel name with zero length'
    endif
  enddo

end subroutine parse_kernel_names


!
!-------------------------------------------------------------------------------------------------
!

! The following utility function was modified from http://fortranwiki.org/fortran/show/strtok
!
subroutine strtok (source_string, delimiter, token)

!     @(#) Tokenize a string in a similar manner to C routine strtok(3c).
!
!     Usage:  First call STRTOK() with the string to tokenize as SOURCE_STRING,
!             and the delimiter list used to tokenize SOURCE_STRING in delimiter.
!
!             then, if the returned value is not equal to char(0), keep calling until it is
!             with SOURCE_STRING set to char(0).
!
!            STRTOK will return a token on each call until the entire line is processed,
!            which it signals by returning char(0).
!
!     Input:  source_string =   Source string to tokenize.
!             delimiter    =   delimiter string.  Used to determine the beginning/end of each token in a string.
!
!     Output: strtok()
!
!     LIMITATIONS:
!     can not be called with a different string until current string is totally processed, even from different procedures

  use postprocess_par,only: MAX_STRING_LEN

  !     PARAMETERS:
  character(len=MAX_STRING_LEN), intent(in)  :: source_string
  character(len=1), intent(in)  :: delimiter
  character(len=MAX_STRING_LEN), intent(out) :: token

  !     SAVED VALUES:
  character(len=MAX_STRING_LEN),save :: saved_string
  integer,save :: isaved_start  ! points to beginning of unprocessed data
  integer,save :: isource_len   ! length of original input string

  !     LOCAL VALUES:
  integer :: ibegin        ! beginning of token to return
  integer :: ifinish       ! end of token to return

  ! initialize stored copy of input string and pointer into input string on first call
  if (source_string(1:1) /= char(0)) then
    isaved_start = 1                 ! beginning of unprocessed data
    saved_string = source_string     ! save input string from first call in series
    isource_len = LEN(saved_string)  ! length of input string from first call
  endif

  token = ''
  ibegin = isaved_start

  ! sets first index ibegin to beginning of (next) token
  do while (.true.)
    if ( (ibegin <= isource_len) .and. (index(delimiter,saved_string(ibegin:ibegin)) /= 0)) then
      ! delimiter is encountered, starts with next index (next token)
      ibegin = ibegin + 1
    else
      ! exits do-loop
      exit
    endif
  enddo

  if (ibegin > isource_len) then
    token = char(0)
    return
  endif

  ! sets second index ifinish to end of token (including delimiter)
  ifinish = ibegin

  do while (.true.)
    if ((ifinish <= isource_len) .and.  (index(delimiter,saved_string(ifinish:ifinish)) == 0)) then
      ! delimiter is not encountered yet, increases finish index
      ifinish = ifinish + 1
    else
      ! exits do-loop
      exit
    endif
  enddo

  ! sets token string
  !strtok = "["//saved_string(ibegin:ifinish-1)//"]"
  token = saved_string(ibegin:ifinish-1)
  isaved_start = ifinish

end subroutine strtok

