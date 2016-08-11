module small_specfem_par

  implicit none

  ! Parameter for simulation

  logical :: is_acoustic=.true.
  ! if is_acoustic=.true.  an Qkappa > 9999 then it is an acoustic simulation
  ! if is_acoustic=.true.  an Qkappa < 9999 then it is an viscoacoustic simulation
  ! if is_acoustic=.false. an Qkappa and Qmu are >9999 then is is an elastic simulation
  ! if is_acoustic=.false. an Qkappa and Qmu are both < 9999 then is is an viscoelastic simulation
  ! logical ::  is_viscoelastic, is_viscoacoustic

  ! density, P wave velocity and S wave velocity of the medium under study
  real(kind=4), parameter :: rho = 2700.
  real(kind=4), parameter :: cp = 3000.
  real(kind=4), parameter :: cs = cp / 1.732
  real(kind=4), parameter :: kappal = rho * cp**2 ! kappal =lambda+2mu, kappal= lambda for acoustic media first Lame coefficient
  real(kind=8), parameter :: QKappa=9999, Qmu=9999

  ! 2-D simulation
  integer, parameter :: NDIM = 2
  ! constant value of the time step in the main time loop, and total number of time steps to simulate
  real, parameter :: deltat = 1.1d-3
  integer, parameter :: NSTEP = 500
  integer it ! time step
  real(kind=4), parameter :: deltatover2 = 0.5*deltat, deltatsqover2 = 0.5*deltat*deltat

  ! to create the mesh
  ! geometry of the model (origin lower-left corner = 0,0) and mesh description
  double precision, parameter :: xmin = 0.d0           ! abscissa of left side of the model
  double precision, parameter :: xmax = 4000.d0        ! abscissa of right side of the model
  double precision, parameter :: zmin = 0.d0           ! abscissa of bottom side of the model
  double precision, parameter :: zmax = 3000.d0        ! abscissa of top side of the model
  integer, parameter :: nx = 80             ! number of spectral elements along X
  integer, parameter :: nz = 60             ! number of spectral elements along Z

  ! number of GLL integration points in each direction of an element (degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLZ = NGLLX
  ! use 4-node elements to describe the geometry of the mesh
  integer, parameter :: ngnod = 4
  ! number of spectral elements and of unique grid points of the mesh
  integer, parameter :: NSPEC = nx * nz
  integer, parameter :: NGLOB = (nx*(NGLLX-1) + 1) * (nz*(NGLLZ-1) + 1)
  ! arrays with the mesh in double precision
  double precision, dimension(NDIM,NGLOB) :: coord
  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLZ,NSPEC) :: ibool

  ! parameters and arrays needed for the simple mesh creation routine
  integer, parameter :: npgeo = (nx+1)*(nz+1) ! total number of geometrical points that describe the geometry
  integer, dimension(ngnod,NSPEC) :: knods ! numbering of the four corners of each mesh element
  double precision, dimension(0:nx,0:nz) :: xgrid,zgrid ! coordinates of all the corners of the mesh elements in a first format
  double precision, dimension(NDIM,npgeo) :: coorg ! coordinates of all the corners of the mesh elements in another format
  integer, external :: num ! function that numbers the mesh points with a unique number


  ! we locate the source and the receiver in arbitrary elements here for this demo code
  integer, parameter :: NSPEC_SOURCE = NSPEC/2 - nx/2, IGLL_SOURCE = 2, JGLL_SOURCE = 2
  integer, parameter :: NSPEC_RECEIVER = 2721, IGLL_RECEIVER = 1, JGLL_RECEIVER = 1
  integer :: iglob_source
  double precision :: x_source,z_source
  double precision :: x_receiver,z_receiver

  ! for the source time function
  double precision, parameter :: f0 = 10.d0
  real, parameter :: t0 = 1.2 / f0
  real, parameter :: factor_amplitude = 1.d+10
  real, parameter :: pi = 3.141592653589793
  real, parameter :: a = pi*pi*f0*f0



  ! Variables for elastic simulation
  ! global displacement, velocity and acceleration vectors
  real(kind=4), dimension(NDIM,NGLOB) :: displ,veloc,accel
  ! global diagonal mass matrix
  real(kind=4), dimension(NGLOB) :: rmass_inverse_elastic
  ! for stiffness matrix
  real(kind=4), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempz1,tempz2
  ! elastic stress
  real(kind=4) sigma_xx,sigma_zz,sigma_xz,sigma_zx
  real(kind=4) lambda,mu,lambdaplus2mu,lambdaplusmu


  !variables for acoustic simulation
  ! logical :: is_acoustic=.true.
  real(kind=4), dimension(NGLOB) :: potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic
  ! global diagonal mass matrix
  real(kind=4), dimension(NGLOB) :: rmass_inverse_acoustic
  real(kind=4), dimension(2,NGLOB) :: acoustic_displ


  ! for attenuation
  integer,  parameter :: N_SLS=3

  real(kind=4), dimension(:,:,:,:), allocatable :: e1,e11,e13
  real(kind=4), dimension(:,:,:,:), allocatable :: inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2
  real(kind=4), dimension(:,:,:) , allocatable :: Mu_nu1,Mu_nu2
  real(kind=4), dimension(:), allocatable :: tau_epsilon_nu1,tau_epsilon_nu2, &
       inv_tau_sigma_nu1_sent,inv_tau_sigma_nu2_sent, &
       phi_nu1_sent,phi_nu2_sent
  real(kind=4) :: Mu_nu1_sent,Mu_nu2_sent




  ! Variables for forcing
  logical, dimension(:), allocatable :: forced
  double precision :: xforced = 0d0
  integer :: nb_forced


  real(kind=4), dimension(NGLLX,NGLLZ,NSPEC) :: xix,xiz,gammax,gammaz,jacobian
  ! array with derivatives of Lagrange polynomials and precalculated products
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLZ) :: zigll,wzgll
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx,hprime_zz,hprimewgll_xx,hprimewgll_zz
  double precision :: xi,gamma,x,z
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl

  ! for gradient
  real(kind=4) dux_dxi,duz_dxi,dux_dgamma,duz_dgamma
  real(kind=4) dux_dxl,dux_dzl,duz_dxl,duz_dzl



  integer :: ispec,iglob,i,j,k,ix,iz,ipoin
  integer :: nglob_to_compute




  real(kind=4) Usolidnorm,current_value,time
  ! real(kind=4), parameter :: VERYSMALLVAL = 1.e-24
  ! displacement threshold above which we consider that the code became unstable
  real(kind=4), parameter :: STABILITY_THRESHOLD = 1.e+25


  !display and info
  integer, parameter :: NTSTEP_BETWEEN_OUTPUT_INFO = 200
  integer, parameter :: IIN = 40
  ! timer to count elapsed time
  character(len=8) :: datein
  character(len=10):: timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  integer ihours,iminutes,iseconds,int_tCPU
  double precision :: time_start,time_end,tCPU




  ! variables for locate receivers
  integer :: nrec = 7
  integer :: irec
  double precision :: hlagrange,dx,dz
  integer,dimension(:), allocatable :: ispec_selected_rec
  double precision, dimension(:), allocatable :: xi_receiver,gamma_receiver
  double precision, dimension(:), allocatable :: st_xval,st_zval
  double precision, dimension(:), allocatable  :: x_final_receiver, z_final_receiver
  ! record a seismogram to check that the simulation went well
  real(kind=4), dimension(:,:), allocatable :: seismogramX,seismogramZ
  character(len=100),dimension(:), allocatable :: seisNameX,seisNameZ
  ! Lagrange interpolators at receivers
  double precision, dimension(:), allocatable :: hxir,hgammar,hpxir,hpgammar
  double precision, dimension(:,:), allocatable :: hxir_store,hgammar_store




end module small_specfem_par
