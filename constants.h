
! size of real and double precision numbers in bytes
  integer, parameter :: SIZE_REAL = 4
  integer, parameter :: SIZE_DOUBLE = 8

! set to SIZE_REAL to run in single precision
! set to SIZE_DOUBLE to run in double precision (increases memory size by 2)
!
! DO NOT forget to change precision_mpi.h accordingly
!
  integer, parameter :: CUSTOM_REAL = SIZE_DOUBLE
! integer, parameter :: CUSTOM_REAL = SIZE_REAL

! polynomial degree
  integer, parameter :: NGLLX = 5
! the code does NOT work if NGLLZ /= NGLLX because it then cannot handle a non-structured mesh
  integer, parameter :: NGLLZ = NGLLX

! further reduce cache misses inner/outer in two passes in the case of an MPI simulation
! this flag is ignored in the case of a serial simulation
  logical, parameter :: FURTHER_REDUCE_CACHE_MISSES = .true.

! for inverse Cuthill-McKee (1969) permutation
  logical, parameter :: PERFORM_CUTHILL_MCKEE = .true.
  logical, parameter :: INVERSE = .true.
  logical, parameter :: FACE = .false.
  integer, parameter :: NGNOD_QUADRANGLE = 4
! perform classical or multi-level Cuthill-McKee ordering
  logical, parameter :: CMcK_MULTI = .false.
! maximum size if multi-level Cuthill-McKee ordering
  integer, parameter :: LIMIT_MULTI_CUTHILL = 50

! implement Cuthill-McKee or replace with identity permutation
  logical, parameter :: ACTUALLY_IMPLEMENT_PERM_OUT = .false.
  logical, parameter :: ACTUALLY_IMPLEMENT_PERM_INN = .false.
  logical, parameter :: ACTUALLY_IMPLEMENT_PERM_WHOLE = .true.

! add MPI barriers and suppress seismograms if we generate traces of the run for analysis with "ParaVer"
  logical, parameter :: GENERATE_PARAVER_TRACES = .false.

! option to display only part of the mesh and not the whole mesh,
! for instance to analyze Cuthill-McKee mesh partitioning etc.
! Possible values are:
!  1: display all the elements (i.e., the whole mesh)
!  2: display inner elements only
!  3: display outer elements only
!  4: display a fixed number of elements (in each partition) only
  integer, parameter :: DISPLAY_SUBSET_OPTION = 1
! number of spectral elements to display in each subset when a fixed subset size is used (option 4 above)
  integer, parameter :: NSPEC_DISPLAY_SUBSET = 2300

!--- beginning of Nicolas Le Goff's constants for an unstructured CUBIT/METIS/SCOTCH mesh

! number of nodes per element
  integer, parameter :: ESIZE = 4

! maximum number of neighbors per element
  integer, parameter :: max_neighbor = 30

! maximum number of elements that can contain the same node
  integer, parameter :: nsize = 20

!--- end of Nicolas Le Goff's constants for an unstructured CUBIT/METIS/SCOTCH mesh

! compute and output acoustic and elastic energy (slows down the code significantly)
  logical, parameter :: OUTPUT_ENERGY = .false.

! output file for energy
  integer, parameter :: IENERGY = 77

! select fast (Paul Fischer) or slow (topology only) global numbering algorithm
  logical, parameter :: FAST_NUMBERING = .true.

! mesh tolerance for fast global numbering
  double precision, parameter :: SMALLVALTOL = 0.00001d0

! displacement threshold above which we consider the code became unstable
  double precision, parameter :: STABILITY_THRESHOLD = 1.d+25

! input and output files
  integer, parameter :: IIN  = 40
  integer, parameter :: ISTANDARD_OUTPUT = 6
! uncomment this to write to standard output
  integer, parameter :: IOUT = ISTANDARD_OUTPUT
! uncomment this to write to file instead
! integer, parameter :: IOUT = 41

! number of lines per source in SOURCE file
  integer, parameter :: NLINES_PER_SOURCE = 13

! flags for absorbing boundaries
  integer, parameter :: IBOTTOM = 1
  integer, parameter :: IRIGHT = 2
  integer, parameter :: ITOP = 3
  integer, parameter :: ILEFT = 4

! number of edges of each element
  integer, parameter :: NEDGES = 4

! a few useful constants
  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0
  double precision, parameter :: HALF = 0.5d0,TWO = 2.d0,QUART = 0.25d0

! pi
  double precision, parameter :: PI = 3.141592653589793d0

! 4/3
  double precision, parameter :: FOUR_THIRDS = 4.d0/3.d0

! 1/24
  double precision, parameter :: ONE_OVER_24 = 1.d0 / 24.d0

! parameters to define the Gauss-Lobatto-Legendre points
  double precision, parameter :: GAUSSALPHA = ZERO,GAUSSBETA = ZERO

! very large and very small values
  double precision, parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9

! number of spatial dimensions
  integer, parameter :: NDIM = 2

! maximum length of station and network name for receivers
  integer, parameter :: MAX_LENGTH_STATION_NAME = 32
  integer, parameter :: MAX_LENGTH_NETWORK_NAME = 8

! number of iterations to solve the system for xi and eta
  integer, parameter :: NUM_ITER = 4

! we mimic a triangle of half duration equal to half_duration_triangle
! using a Gaussian having a very close shape, as explained in Figure 4.2
! of the manual. This source decay rate to mimic an equivalent triangle
! was found by trial and error
  double precision, parameter :: SOURCE_DECAY_MIMIC_TRIANGLE = 1.628d0

! non linear display to enhance small amplitudes in color images
  double precision, parameter :: POWER_DISPLAY_COLOR = 0.30d0

! US letter paper or European A4
  logical, parameter :: US_LETTER = .false.

! X and Z axis origin of PostScript plot in centimeters
  double precision, parameter :: ORIG_X = 2.4d0
  double precision, parameter :: ORIG_Z = 2.9d0

! dot to centimeter conversion for PostScript
  double precision, parameter :: CENTIM = 28.5d0

! parameters for arrows for PostScript snapshot
  double precision, parameter :: ARROW_ANGLE = 20.d0
  double precision, parameter :: ARROW_RATIO = 0.40d0

! size of frame used for Postscript display in percentage of the size of the page
  double precision, parameter :: RPERCENTX = 70.0d0,RPERCENTZ = 77.0d0

!-----------------------------------------------------------------------

!
! anisotropic copper crystal (cubic symmetry)
!

! regular c_ijkl with no rotation
  double precision, parameter :: c11val = 169.d9
  double precision, parameter :: c12val = 122.d9
  double precision, parameter :: c13val = c12val
  double precision, parameter :: c14val = 0.d0
  double precision, parameter :: c15val = 0.d0
  double precision, parameter :: c16val = 0.d0

  double precision, parameter :: c22val = c11val
  double precision, parameter :: c23val = c12val
  double precision, parameter :: c24val = 0.d0
  double precision, parameter :: c25val = 0.d0
  double precision, parameter :: c26val = 0.d0

  double precision, parameter :: c33val = c11val
  double precision, parameter :: c34val = 0.d0
  double precision, parameter :: c35val = 0.d0
  double precision, parameter :: c36val = 0.d0

  double precision, parameter :: c44val = 75.3d9
  double precision, parameter :: c45val = 0.d0
  double precision, parameter :: c46val = 0.d0

  double precision, parameter :: c55val = c44val
  double precision, parameter :: c56val = 0.d0

  double precision, parameter :: c66val = c44val

