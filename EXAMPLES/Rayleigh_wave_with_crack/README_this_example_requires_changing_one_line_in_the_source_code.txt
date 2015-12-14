----------------------------------------------------------------------
README
----------------------------------------------------------------------

##################
##################  IMPORTANT: for this example (only) do not change NPROC = 1 in DATA/Par_file to any other value because
##################  IMPORTANT: this example with a crack in the mesh can ony run in serial mode (I use a trick
##################  IMPORTANT: inside the code to create a mathematical crack, and that trick is not ported to MPI yet)
##################

*********************************************************************************************************
 This example requires setting:
! add a small crack (discontinuity) in the medium manually
  logical, parameter :: ADD_A_SMALL_CRACK_IN_THE_MEDIUM = .true.
 in the source code file src/specfem2D/specfem2D.F90
*********************************************************************************************************

