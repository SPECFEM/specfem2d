
! polynomial degree
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = NGLLX

! select fast (Paul Fischer) or slow (topology only) global numbering algorithm
  logical, parameter :: FAST_NUMBERING = .true.

! mesh tolerance for fast global numbering
  double precision, parameter :: SMALLVALTOL = 0.000001d0

! input and output files
  integer, parameter :: IIN  = 40

! uncomment this to write to standard output
  integer, parameter :: IOUT = 6
! uncomment this to write to file instead
! integer, parameter :: IOUT = 41

! flags for absorbing boundaries
  integer, parameter :: IHAUT   = 1
  integer, parameter :: IBAS    = 2
  integer, parameter :: IGAUCHE = 3
  integer, parameter :: IDROITE = 4

  integer, parameter :: IARETEBAS    = 1
  integer, parameter :: IARETEDROITE = 2
  integer, parameter :: IARETEHAUT   = 3
  integer, parameter :: IARETEGAUCHE = 4

! a few useful constants
  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0
  double precision, parameter :: HALF = 0.5d0,TWO = 2.0d0
  double precision, parameter :: PI = 3.141592653589793d0

! parameters to define the Gauss-Lobatto-Legendre points
  double precision, parameter :: GAUSSALPHA = ZERO,GAUSSBETA = ZERO

! large value for maximum
  double precision, parameter :: HUGEVAL = 1.d+30

! number of spatial dimensions
  integer, parameter :: NDIME = 2

! X and Z scaling du display pour PostScript
  double precision, parameter :: SCALEX = 1.d0
  double precision, parameter :: SCALEZ = 1.d0

! taille de la plus grande fleche en centimetres
  double precision, parameter :: SIZEMAX = 1.d0

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

