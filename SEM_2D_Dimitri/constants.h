
! polynomial degree
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLZ = NGLLX

! select fast (Paul Fischer) or slow (topology only) global numbering algorithm
  logical, parameter :: FAST_NUMBERING = .true.

! mesh tolerance for fast global numbering
  double precision, parameter :: SMALLVALTOL = 0.000001d0

! displacement threshold above which we consider the code became unstable
  double precision, parameter :: STABILITY_THRESHOLD = 1.d+25

! input and output files
  integer, parameter :: IIN  = 40

! uncomment this to write to standard output
  integer, parameter :: IOUT = 6
! uncomment this to write to file instead
! integer, parameter :: IOUT = 41

! flags for absorbing boundaries
  integer, parameter :: ITOP = 1
  integer, parameter :: IBOTTOM = 2
  integer, parameter :: ILEFT = 3
  integer, parameter :: IRIGHT = 4

! a few useful constants
  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0
  double precision, parameter :: HALF = 0.5d0,TWO = 2.0d0,QUART = 0.25d0

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

! error function source decay rate for Heaviside
  double precision, parameter :: SOURCE_DECAY_RATE = 1.628d0

! display non lineaire pour rehausser les faibles amplitudes sur les images PNM
  double precision, parameter :: POWER_DISPLAY_PNM = 0.30d0

! X and Z scaling du display pour PostScript
  double precision, parameter :: SCALEX = 1.d0
  double precision, parameter :: SCALEZ = 1.d0

! US letter paper or European A4
  logical, parameter :: US_LETTER = .false.

! write symbols on PostScript display
  logical, parameter :: ISYMBOLS = .true.

! X and Z axis origin of PostScript plot in centimeters
  double precision, parameter :: ORIG_X = 2.4d0
  double precision, parameter :: ORIG_Z = 2.9d0

! dot to centimeter conversion for PostScript
  double precision, parameter :: CENTIM = 28.5d0

! parameters for arrows for PostScript snapshot
  double precision, parameter :: ANGLE = 20.d0
  double precision, parameter :: RAPPORT = 0.40d0

! ecrire legendes ou non in PostScript display
  logical, parameter :: LEGENDES = .true.

! limite pour afficher des points a la place des recepteurs
  integer, parameter :: NDOTS = 30

! taille de la fenetre de display Postscript en pourcentage de la feuille
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

