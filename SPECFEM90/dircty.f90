!=====================================================================
!
!                 S p e c f e m  V e r s i o n  4 . 2
!                 -----------------------------------
!
!                         Dimitri Komatitsch
!    Department of Earth and Planetary Sciences - Harvard University
!                         Jean-Pierre Vilotte
!                 Departement de Sismologie - IPGP - Paris
!                           (c) June 1998
!
!=====================================================================

  subroutine dircty
!
!=======================================================================
!
!     Dynamic storage allocation :
!     --------------------------
!
!     Print a directory listing of all dynamically allocated arrays
!       and their properties
!
!=======================================================================

  use iounit
  use arraydir

  implicit none

  integer itotsize,iarray
  character(len=7) label(3)
  integer isizevars(3)

! ici codage en dur des tailles des variables en octets
  isizevars(1) = 4  ! integer
  isizevars(2) = 4  ! single precision
  isizevars(3) = 8  ! double precision

  label(1) = 'Integer'
  label(2) = 'Real   '
  label(3) = 'Double '

! compute total size in bytes
  itotsize = 0
  do iarray = 1,nbarrays
    itotsize = itotsize + arraysizes(iarray)*isizevars(arraytypes(iarray))
  enddo

  write(iout,100) nbarrays,dble(itotsize)/dble(1024*1024),itotsize, &
                      itotsize/isizevars(3)

  do iarray = 1,nbarrays
    write(iout,110) iarray,arraysizes(iarray),arraynames(iarray), &
        label(arraytypes(iarray))
  enddo

  100   format(//1x,41('=')/ &
  ' =  D i r e c t o r y     l i s t i n g  ='/1x,41('=')// &
  ' Total number of allocated arrays. . . . . . . . . .',i11/ &
  ' Total size of arrays in megabytes . . . . . . . . .',f11.3/ &
  ' Total size of arrays in bytes . . . . . . . . . . .',i11/ &
  ' Total size of arrays in double precision words. . .',i11/// &
  '  Array nb    Size         Name        Type'/1x,47('=')/)
  110   format(i6,3x,i10,5x,a12,2x,a7)

  return
  end subroutine dircty
