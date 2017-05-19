----------------------------------------------------------------------
README
----------------------------------------------------------------------

*********************************************************************************************************
 This example requires setting:

! select fast (Paul Fischer) or slow (topology only) global numbering algorithm
  logical, parameter :: FAST_NUMBERING = .false.

and

! add a small crack (discontinuity) in the medium manually
  logical, parameter :: ADD_A_SMALL_CRACK_IN_THE_MEDIUM = .true.

 in the source code file setup/constants.h
*********************************************************************************************************

