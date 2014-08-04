
! select single precision or double precision for all the calculations
! since this 2D example is cheap to run it is better to use double precision
! integer, parameter :: CUSTOM_REAL = 4
  integer, parameter :: CUSTOM_REAL = 8

! number of spatial dimensions of the problem
! here we define a 2D formulation
  integer, parameter :: NDIM = 2

! number of GLL points in each of the two directions of a given spectral element
! they MUST be equal in the case of a non-structured mesh because the two directions
! can correspond if rotated elements are in contact
  integer, parameter :: NGLLX = 5, NGLLZ = NGLLX

! radius of the spherical earth (in meters)
  double precision, parameter :: R_EARTH = 6371000.d0

! we define the geometrical mesh using geometrical elements that have 9 control nodes instead of 4
! in order to be able to take the curvature of the Earth into account, i.e. we need to have elements
! whose edges can be curved rather than being straight lines
  integer, parameter :: NGNOD = 9

! number of iterations to solve the system for xi and gamma
  integer, parameter :: NUM_ITER = 4

! a few useful constants
  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0
  double precision, parameter :: HALF = 0.5d0,TWO = 2.d0,QUARTER = 0.25d0

! pi
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: TWO_PI = 2.d0 * PI

! very large and very small values
  double precision, parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9

! to convert angles from degrees to radians
  double precision, parameter :: DEGREES_TO_RADIANS = PI / 180.d0
! to convert angles from radians to degrees
  double precision, parameter :: RADIANS_TO_DEGREES = 180.d0 / PI

