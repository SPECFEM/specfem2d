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

  subroutine storearray(name,isize,itype)
!
!=======================================================================
!
!     Dynamic storage : store the array properties
!     ----------------
!
!=======================================================================

  use iounit
  use arraydir

  implicit none

  character(len=*) name
  integer isize,itype

  if(itype /= iinteg .and. itype /= isngl .and. itype /= idouble) &
    stop 'Wrong array type in dynamic allocation'

  if(isize <= 0) &
    stop 'Incoherent array size in dynamic allocation'

  nbarrays = nbarrays + 1
  if(nbarrays > maxnbarrays) stop 'Maximum number of arrays reached'

  arraysizes(nbarrays) = isize
  arraytypes(nbarrays) = itype
  arraynames(nbarrays) = name

  return
  end subroutine storearray
