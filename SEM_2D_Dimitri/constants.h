
! polynomial degree
  integer, parameter :: NGLLX = 5
! the code does NOT work if NGLLZ /= NGLLX because it then cannot handle a non-structured mesh
  integer, parameter :: NGLLZ = NGLLX

! select fast (Paul Fischer) or slow (topology only) global numbering algorithm
  logical, parameter :: FAST_NUMBERING = .true.

! mesh tolerance for fast global numbering
  double precision, parameter :: SMALLVALTOL = 0.000001d0

! displacement threshold above which we consider the code became unstable
  double precision, parameter :: STABILITY_THRESHOLD = 1.d+25

! input and output files
  integer, parameter :: IIN  = 40
  integer, parameter :: ISTANDARD_OUTPUT = 6
! uncomment this to write to standard output
  integer, parameter :: IOUT = ISTANDARD_OUTPUT
! uncomment this to write to file instead
! integer, parameter :: IOUT = 41

! flags for absorbing boundaries
  integer, parameter :: IBOTTOM = 1
  integer, parameter :: IRIGHT = 2
  integer, parameter :: ITOP = 3
  integer, parameter :: ILEFT = 4

! number of edges of each element
  integer, parameter :: NEDGES = 4

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

! attenuation constants from J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).
! Beware: these values implement specific values of the quality factors:
! Qp approximately equal to 13 and Qs approximately equal to 10,
! which means very high attenuation, see that paper for details.
! double precision, parameter :: tau_epsilon_nu1_mech1 = 0.0334d0
! double precision, parameter :: tau_sigma_nu1_mech1   = 0.0303d0
! double precision, parameter :: tau_epsilon_nu2_mech1 = 0.0352d0
! double precision, parameter :: tau_sigma_nu2_mech1   = 0.0287d0

! double precision, parameter :: tau_epsilon_nu1_mech2 = 0.0028d0
! double precision, parameter :: tau_sigma_nu1_mech2   = 0.0025d0
! double precision, parameter :: tau_epsilon_nu2_mech2 = 0.0029d0
! double precision, parameter :: tau_sigma_nu2_mech2   = 0.0024d0

! attenuation constants from J. M. Carcione, D. Kosloff and R. Kosloff,
! Wave propagation simulation in a linear viscoelastic medium, Geophysical Journal International,
! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).
! Beware: these values implement specific values of the quality factors:
! Qp approximately equal to 27 and Qs approximately equal to 20,
! which means very high attenuation, see that paper for details.
  double precision, parameter :: tau_epsilon_nu1_mech1 = 0.0325305d0
  double precision, parameter :: tau_sigma_nu1_mech1   = 0.0311465d0
  double precision, parameter :: tau_epsilon_nu2_mech1 = 0.0332577d0
  double precision, parameter :: tau_sigma_nu2_mech1   = 0.0304655d0

  double precision, parameter :: tau_epsilon_nu1_mech2 = 0.0032530d0
  double precision, parameter :: tau_sigma_nu1_mech2   = 0.0031146d0
  double precision, parameter :: tau_epsilon_nu2_mech2 = 0.0033257d0
  double precision, parameter :: tau_sigma_nu2_mech2   = 0.0030465d0

  double precision, parameter :: inv_tau_sigma_nu1_mech1 = ONE / tau_sigma_nu1_mech1
  double precision, parameter :: inv_tau_sigma_nu2_mech1 = ONE / tau_sigma_nu2_mech1
  double precision, parameter :: inv_tau_sigma_nu1_mech2 = ONE / tau_sigma_nu1_mech2
  double precision, parameter :: inv_tau_sigma_nu2_mech2 = ONE / tau_sigma_nu2_mech2

  double precision, parameter :: phi_nu1_mech1 = (ONE - tau_epsilon_nu1_mech1/tau_sigma_nu1_mech1) / tau_sigma_nu1_mech1
  double precision, parameter :: phi_nu2_mech1 = (ONE - tau_epsilon_nu2_mech1/tau_sigma_nu2_mech1) / tau_sigma_nu2_mech1
  double precision, parameter :: phi_nu1_mech2 = (ONE - tau_epsilon_nu1_mech2/tau_sigma_nu1_mech2) / tau_sigma_nu1_mech2
  double precision, parameter :: phi_nu2_mech2 = (ONE - tau_epsilon_nu2_mech2/tau_sigma_nu2_mech2) / tau_sigma_nu2_mech2

  double precision, parameter :: Mu_nu1 = ONE - (ONE - tau_epsilon_nu1_mech1/tau_sigma_nu1_mech1) &
                                              - (ONE - tau_epsilon_nu1_mech2/tau_sigma_nu1_mech2)
  double precision, parameter :: Mu_nu2 = ONE - (ONE - tau_epsilon_nu2_mech1/tau_sigma_nu2_mech1) &
                                              - (ONE - tau_epsilon_nu2_mech2/tau_sigma_nu2_mech2)

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

