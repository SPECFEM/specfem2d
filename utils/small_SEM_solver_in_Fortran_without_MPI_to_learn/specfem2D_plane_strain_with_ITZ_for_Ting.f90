!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

  program serial_specfem2D

  implicit none

!!!!!!!!!!
!!!!!!!!!! NGLLX and NGLLZ are set equal to 5,
!!!!!!!!!! therefore each element contains NGLLX * NGLLZ = 25 points.
!!!!!!!!!!

!!!!!!!!!!
!!!!!!!!!! All the calculations can be done in single precision.
!!!!!!!!!! We do not really need double precision in SPECFEM2D.
!!!!!!!!!! If you thus want to use single precision, just change the value in the include file from 8 to 4
!!!!!!!!!!
  include "precision.h"

! density, P wave velocity and S wave velocity of the geophysical medium under study
  real(kind=CUSTOM_REAL), parameter :: rho = 2700.
  real(kind=CUSTOM_REAL), parameter :: cp = 3000.
  real(kind=CUSTOM_REAL), parameter :: cs = cp / 1.732

! to create the mesh
! geometry of the model (origin lower-left corner = 0,0) and mesh description
  double precision, parameter :: xmin = 0.d0           ! abscissa of left side of the model
  double precision, parameter :: xmax = 4000.d0        ! abscissa of right side of the model
  double precision, parameter :: zmin = 0.d0           ! abscissa of bottom side of the model
  double precision, parameter :: zmax = 3000.d0        ! abscissa of top side of the model
  integer, parameter :: nelem_x = 80             ! number of spectral elements along X
  integer, parameter :: nelem_z = 60             ! number of spectral elements along Z

!! DK DK added this for example of contour integral for Ting
  integer, parameter :: EXTERNAL_SIZE_X = 6
  integer, parameter :: EXTERNAL_SIZE_Z = 6

! number of GLL integration points in each direction of an element (degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLZ = NGLLX

! use 4-node elements to describe the geometry of the mesh
  integer, parameter :: ngnod = 4

! number of spectral elements and of unique grid points of the mesh
  integer, parameter :: NSPEC = nelem_x * nelem_z
  integer, parameter :: NGLOB = (nelem_x*(NGLLX-1) + 1) * (nelem_z*(NGLLZ-1) + 1) &
               + nelem_x*(NGLLX-1) + 1  !! we now add this second line because the points on the ITZ interface are duplicated

! constant value of the time step in the main time loop, and total number of time steps to simulate
  real(kind=CUSTOM_REAL), parameter :: deltat = 1.1e-3_CUSTOM_REAL
  integer, parameter :: NSTEP = 1600

! we locate the source and the receiver in arbitrary elements here for this demo code
! the source is placed exactly in the middle of the grid here
  integer, parameter :: NSPEC_SOURCE = NSPEC/4 - nelem_x/2, IGLL_SOURCE = NGLLX, JGLL_SOURCE = NGLLZ

  integer, parameter :: NSPEC_RECEIVER = 2*NSPEC/3 - nelem_x/4, IGLL_RECEIVER = 1, JGLL_RECEIVER = 1

! for the source time function
  real(kind=CUSTOM_REAL), parameter :: f0 = 10.
  real(kind=CUSTOM_REAL), parameter :: t0 = 1.2 / f0
  real(kind=CUSTOM_REAL), parameter :: factor_amplitude = 1.e+10
  real(kind=CUSTOM_REAL), parameter :: pi = 3.141592653589793
  real(kind=CUSTOM_REAL), parameter :: a = pi*pi*f0*f0

  integer, parameter :: NTSTEP_BETWEEN_OUTPUT_INFO = 100

  integer, parameter :: IIN = 40

! 2-D simulation
  integer, parameter :: NDIM = 2

  real(kind=CUSTOM_REAL), parameter :: deltatover2 = 0.5*deltat, deltatsqover2 = 0.5*deltat*deltat

! real(kind=CUSTOM_REAL), parameter :: VERYSMALLVAL = 1.e-24

! displacement threshold above which we consider that the code became unstable
  real(kind=CUSTOM_REAL), parameter :: STABILITY_THRESHOLD = 1.e+25

! global displacement, velocity and acceleration vectors
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: displ,veloc,accel

! global diagonal mass matrix
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: rmass_inverse

! record a seismogram to check that the simulation went well
  real(kind=CUSTOM_REAL), dimension(NSTEP) :: seismogram

! time step
  integer it

! arrays with mesh topology
  integer, dimension(NGLLX,NGLLZ,NSPEC) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,NSPEC) :: xix,xiz,gammax,gammaz,jacobian

! arrays with the mesh in double precision
  double precision, dimension(NDIM,NGLOB) :: coord

!! DK DK added this for example of contour integral for Ting
  logical, dimension(NSPEC) :: is_on_Gamma_plus
  logical, dimension(NSPEC) :: is_on_Gamma_minus
  real(kind=CUSTOM_REAL) :: nx,nz,xxi,zxi,weight,jacobian1D
  real(kind=CUSTOM_REAL) :: my_function_to_integrate,test_integral_x,test_integral_z,exact

  double precision :: x_source,z_source
  double precision :: x_receiver,z_receiver

! array with derivatives of Lagrange polynomials and precalculated products
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLZ) :: zigll,wzgll
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx,hprime_zz,hprimewgll_xx,hprimewgll_zz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempz1,tempz2

  integer :: ispec,iglob,i,j,k,ix,iz,ipoin
  integer :: nglob_to_compute

  double precision :: xi,gamma,x,z
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) dux_dxi,duz_dxi,dux_dgamma,duz_dgamma
  real(kind=CUSTOM_REAL) dux_dxl,dux_dzl,duz_dxl,duz_dzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_zz,sigma_xz,sigma_zx
  real(kind=CUSTOM_REAL) lambda,mu,lambdaplus2mu

  real(kind=CUSTOM_REAL) Usolidnorm,current_value,time

! parameters and arrays needed for the simple mesh creation routine

! total number of geometrical points that describe the geometry
  integer, parameter :: npgeo = (nelem_x+1)*(nelem_z+1) &
               + nelem_x+1  !! we now add this second line because the points on the ITZ interface are duplicated

! numbering of the four corners of each mesh element
  integer, dimension(ngnod,NSPEC) :: knods

! coordinates of all the corners of the mesh elements in a first format
  double precision, dimension(0:nelem_x,0:nelem_z) :: xgrid,zgrid

! coordinates of all the corners of the mesh elements in another format
  double precision, dimension(NDIM,npgeo) :: coorg

! function that numbers the mesh points with a unique number
  integer, external :: num

!! height of the ITZ interface to create
  double precision :: h

  integer :: value_to_add_to_ipoin,value_to_add_to_knods
  double precision :: value_to_add_to_zgrid

! total number of sources and of receivers
  integer, parameter :: NSOURCES = 1,nrec = 1

! to create JPEG color snapshots of the results (optional)
  integer  :: NX_IMAGE_color,NZ_IMAGE_color,isnapshot_number = 0
  double precision  :: xmin_color_image,xmax_color_image,zmin_color_image,zmax_color_image
  integer, dimension(:,:), allocatable :: iglob_image_color
  double precision, dimension(:,:), allocatable :: image_color_data
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: vector_field_display
  integer, dimension(NSOURCES) :: ix_image_color_source,iy_image_color_source
  integer, dimension(nrec) :: ix_image_color_receiver,iy_image_color_receiver

! timer to count elapsed time
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  integer ihours,iminutes,iseconds,int_tCPU
  double precision :: time_start,time_end,tCPU

  print *
  print *,'NSPEC = ',NSPEC
  print *,'NGLOB = ',NGLOB
  print *

  print *,'NSTEP = ',NSTEP
  print *,'deltat = ',deltat
  print *

! set up Gauss-Lobatto-Legendre points, weights and also derivation matrices
  call define_derivation_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,hprimewgll_xx,hprimewgll_zz,NGLLX,NGLLZ)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! DK DK added this for example of contour integral for Ting
  is_on_Gamma_plus(:) = .false.
  is_on_Gamma_minus(:) = .false.

!!!! create the first half of the mesh (the bottom part)

  knods(:,:) = 0

  xgrid(:,:) = 0
  zgrid(:,:) = 0

  value_to_add_to_ipoin = 0

!--- definition of the mesh
  do iz = 0,nelem_z/2
    do ix = 0,nelem_x
        ! coordinates of the grid points (evenly spaced points along X and Z)
        xgrid(ix,iz) = xmin + (xmax - xmin) * dble(ix) / dble(nelem_x)
        zgrid(ix,iz) = zmin + (zmax - zmin) * dble(iz) / dble(nelem_z)
     enddo
  enddo

! create the coorg array
  do j = 0,nelem_z/2
    do i = 0,nelem_x
      ipoin = num(i,j,nelem_x)
      value_to_add_to_ipoin = max(value_to_add_to_ipoin,ipoin)
      coorg(1,ipoin) = xgrid(i,j)
      coorg(2,ipoin) = zgrid(i,j)
    enddo
  enddo

! create the knods array
  ispec = 0
  do j=0,nelem_z/2-1
    do i=0,nelem_x-1
      ispec = ispec + 1
      knods(1,ispec) = num(i,j,nelem_x)
      knods(2,ispec) = num(i+1,j,nelem_x)
      knods(3,ispec) = num(i+1,j+1,nelem_x)
      knods(4,ispec) = num(i,j+1,nelem_x)

!! DK DK create a flag for all the elements that are in contact with the Gamma_minus contour
      if (j == nelem_z/2-1) is_on_Gamma_minus(ispec) = .true.

    enddo
  enddo

!!!! create the second half of the mesh (the top part)

  value_to_add_to_knods = maxval(knods)

! total height of the first half of the mesh
  value_to_add_to_zgrid = maxval(zgrid)

!! height of the ITZ interface to create (set to 20% of the vertical size of an element here)
  h = 0.20d0 * (zmax - zmin) / dble(nelem_z)

!--- definition of the mesh
  do iz = 0,nelem_z/2
    do ix = 0,nelem_x
        ! coordinates of the grid points (evenly spaced points along X and Z)
        xgrid(ix,iz) = xmin + (xmax - xmin) * dble(ix) / dble(nelem_x)
        !! here is where we add the thickness of the ITZ interface
        zgrid(ix,iz) = zmin + (zmax - zmin) * dble(iz) / dble(nelem_z) + value_to_add_to_zgrid + h
     enddo
  enddo

! create the coorg array
  do j = 0,nelem_z/2
    do i = 0,nelem_x
      ipoin = num(i,j,nelem_x) + value_to_add_to_ipoin
      coorg(1,ipoin) = xgrid(i,j)
      coorg(2,ipoin) = zgrid(i,j)
    enddo
  enddo

! create the knods array
  do j=0,nelem_z/2-1
    do i=0,nelem_x-1
      ispec = ispec + 1
      knods(1,ispec) = num(i,j,nelem_x) + value_to_add_to_knods
      knods(2,ispec) = num(i+1,j,nelem_x) + value_to_add_to_knods
      knods(3,ispec) = num(i+1,j+1,nelem_x) + value_to_add_to_knods
      knods(4,ispec) = num(i,j+1,nelem_x) + value_to_add_to_knods

!! DK DK create a flag for all the elements that are in contact with the Gamma_plus contour
      if (j == 0) is_on_Gamma_plus(ispec) = .true.

    enddo
  enddo

  if (ispec /= NSPEC) stop 'error in the total number of spectral elements created in the mesh'

!
!---- generate the global numbering
!

  call createnum_slow(knods,ibool,nglob_to_compute,nspec,NGLLX,NGLLZ,ngnod) ! Create ibool and recompute nglob for checking
  if (nglob_to_compute /= NGLOB) stop 'error: incorrect total number of unique grid points found'
  if (minval(ibool) /= 1) stop 'error: incorrect minimum value of ibool'
  if (maxval(ibool) /= NGLOB) stop 'error: incorrect maximum value of ibool'

!
!----  set the coordinates of the points of the global grid
!
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          xi = xigll(i)
          gamma = zigll(j)

          call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl, &
                          jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo,NDIM)
          if (jacobianl <= 0.d0) stop 'error: negative Jacobian found'

          coord(1,ibool(i,j,ispec)) = x
          coord(2,ibool(i,j,ispec)) = z

          xix(i,j,ispec) = xixl
          xiz(i,j,ispec) = xizl
          gammax(i,j,ispec) = gammaxl
          gammaz(i,j,ispec) = gammazl
          jacobian(i,j,ispec) = jacobianl

        enddo
      enddo
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! build the mass matrix
  rmass_inverse(:) = 0.
  do ispec = 1,NSPEC
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        rmass_inverse(iglob) = rmass_inverse(iglob) + wxgll(i)*wzgll(j)*rho*jacobian(i,j,ispec)
      enddo
    enddo
  enddo

! we have built the real exactly diagonal mass matrix (not its inverse)
! therefore invert it here once and for all
  do i = 1,NGLOB
    rmass_inverse(i) = 1. / rmass_inverse(i)
  enddo

! compute the position of the source and of the receiver
  x_source = coord(1,ibool(IGLL_SOURCE,JGLL_SOURCE,NSPEC_SOURCE))
  z_source = coord(2,ibool(IGLL_SOURCE,JGLL_SOURCE,NSPEC_SOURCE))

  x_receiver = coord(1,ibool(IGLL_RECEIVER,JGLL_RECEIVER,NSPEC_RECEIVER))
  z_receiver = coord(2,ibool(IGLL_RECEIVER,JGLL_RECEIVER,NSPEC_RECEIVER))

  print *
  print *,'x_source = ',x_source
  print *,'z_source = ',z_source
  print *
  print *,'x_receiver = ',x_receiver
  print *,'z_receiver = ',z_receiver
  print *

! to prepare for the creation of color JPEG snapshots of the results (optional)
  call prepare_color_image_init(NDIM,NGLOB,NGLLX,coord,npgeo,NX_IMAGE_color,NZ_IMAGE_color, &
                          xmin_color_image,xmax_color_image,zmin_color_image,zmax_color_image)

  allocate(iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color))
  allocate(image_color_data(NX_IMAGE_color,NZ_IMAGE_color))

  call prepare_color_image_pixels(ngnod,npgeo,nspec,NDIM,NGLOB,NGLLX,NGLLZ,NX_IMAGE_color,NZ_IMAGE_color,NSOURCES,nrec, &
        xmin_color_image,xmax_color_image,zmin_color_image,zmax_color_image,x_source,z_source,x_receiver,z_receiver, &
        coord,coorg,knods,ibool,iglob_image_color,ix_image_color_source,iy_image_color_source, &
        ix_image_color_receiver,iy_image_color_receiver)
!
!---- compute a test 1D integral along Gamma_minus to create an example for Ting
!

!! DK DK added this for example of Gamma_minus integral for Ting

  test_integral_x = 0.
  test_integral_z = 0.

  do ispec = 1,nspec

  if (is_on_Gamma_minus(ispec)) then
! Gamma_minus is composed of the top edges (j == NGLLZ) of the spectral elements in contact with the mesh interface of height h
! (represented in yellow in the PostScript visualization files produced by the code)
    j = NGLLZ
    do i = 1,NGLLX
        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
!! BEWARE here: the normal is (0, +1), not (0, -1), because we are at the top of a spectral element
!! Ting, you may need to invert the sign of the normal here, depending if you want it to go down or up along Gamma_minus
        weight = jacobian1D * wxgll(i)

        x = coord(1,ibool(i,j,ispec))
        z = coord(2,ibool(i,j,ispec))
        !!!!!!!!!!!! print *,'z (should be 1500) = ',z
        my_function_to_integrate = x
        test_integral_x = test_integral_x + weight*nx*my_function_to_integrate
        test_integral_z = test_integral_z + weight*nz*my_function_to_integrate
    enddo
  endif

  enddo

  print *
! this one is exactly zero because the interface is horizontal and thus the normal vector (nx, nz) is always (0, 1)
  exact = 0.d0
  print *,'Value of test_integral_x  numerical = ',test_integral_x,'  exact = ',exact,'  difference = ',test_integral_x - exact
! the sign is positive here because the normal is (0, +1), not (0, -1), because we are at the top of a spectral element
  exact = (xmax**2 - xmin**2)/2
  print *,'Value of test_integral_z  numerical = ',test_integral_z,'  exact = ',exact,'  difference = ',test_integral_z - exact
  print *

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!---- compute a test 1D integral along Gamma_plus to create an example for Ting
!

!! DK DK added this for example of Gamma_plus integral for Ting

  test_integral_x = 0.
  test_integral_z = 0.

  do ispec = 1,nspec

  if (is_on_Gamma_plus(ispec)) then
! Gamma_plus is composed of the bottom edges (j == 1) of the spectral elements in contact with the mesh interface of height h
! (represented in green in the PostScript visualization files produced by the code)
    j = 1
    do i = 1,NGLLX
        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
!! BEWARE here: the normal is (0, -1), not (0, +1), because we are at the bottom of a spectral element
!! Ting, you may need to invert the sign of the normal here, depending if you want it to go up or down along Gamma_plus
        weight = jacobian1D * wxgll(i)

        x = coord(1,ibool(i,j,ispec))
        z = coord(2,ibool(i,j,ispec))
        !!!!!!!!!!!! print *,'z (should be ',1500 + h,') = ',z
        my_function_to_integrate = x
        test_integral_x = test_integral_x + weight*nx*my_function_to_integrate
        test_integral_z = test_integral_z + weight*nz*my_function_to_integrate
    enddo
  endif

  enddo

  print *
! this one is exactly zero because the interface is horizontal and thus the normal vector (nx, nz) is always (0, -1)
  exact = 0.d0
  print *,'Value of test_integral_x  numerical = ',test_integral_x,'  exact = ',exact,'  difference = ',test_integral_x - exact
! the sign is negative here because the normal is (0, -1), not (0, +1), because we are at the bottom of a spectral element
  exact = - (xmax**2 - xmin**2)/2
  print *,'Value of test_integral_z  numerical = ',test_integral_z,'  exact = ',exact,'  difference = ',test_integral_z - exact
  print *

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! clear initial vectors before starting the time loop
  displ(:,:) = 0. !!!!!!!!!! VERYSMALLVAL
  veloc(:,:) = 0.
  accel(:,:) = 0.

! count elapsed wall-clock time
  call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! this fails if we cross the end of the month
  time_start = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
               60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0

! start of the time loop (which must remain serial obviously)
  do it = 1,NSTEP

! compute maximum of norm of displacement from time to time and display it
! in order to monitor the simulation
    if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
      Usolidnorm = -1.
      do iglob = 1,NGLOB
        current_value = sqrt(displ(1,iglob)**2 + displ(2,iglob)**2)
        if (current_value > Usolidnorm) Usolidnorm = current_value
      enddo
      write(*,*) 'Time step # ',it,' out of ',NSTEP
! compute current time
      time = (it-1)*deltat
      write(*,*) 'Max norm displacement vector U in the solid (m) = ',Usolidnorm
! check stability of the code, exit if unstable
      if (Usolidnorm > STABILITY_THRESHOLD .or. Usolidnorm < 0) stop 'code became unstable and blew up'

! count elapsed wall-clock time
  call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! this fails if we cross the end of the month
  time_end = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
             60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0

! draw the displacement vector field in a PostScript file (optional)
      call plot_post(displ,coord,ibool,NGLOB,NSPEC,x_source,z_source,x_receiver,z_receiver,it,deltat,NGLLX,NGLLZ,NDIM, &
                          is_on_Gamma_minus,is_on_Gamma_plus)

! draw a color JPEG snapshot of the results (optional)
      call write_color_image_snaphot(it,NX_IMAGE_color,NZ_IMAGE_color,NDIM,NGLOB,displ,veloc,accel, &
                  vector_field_display,image_color_data,iglob_image_color,ix_image_color_source,iy_image_color_source, &
                  ix_image_color_receiver,iy_image_color_receiver,isnapshot_number,NSOURCES,nrec,cp)

! elapsed time since beginning of the simulation
  tCPU = time_end - time_start
  int_tCPU = int(tCPU)
  ihours = int_tCPU / 3600
  iminutes = (int_tCPU - 3600*ihours) / 60
  iseconds = int_tCPU - 3600*ihours - 60*iminutes
  write(*,*) 'Elapsed time in seconds = ',tCPU
  write(*,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
  write(*,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
  write(*,*)

    endif

! big loop over all the global points (not elements) in the mesh to update
! the displacement and velocity vectors and clear the acceleration vector
  displ(:,:) = displ(:,:) + deltat*veloc(:,:) + deltatsqover2*accel(:,:)
  veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)
  accel(:,:) = 0.

! big loop over all the elements in the mesh to localize data
! from the global vectors to the local mesh
! using indirect addressing (contained in array ibool)
! and then to compute the elemental contribution
! to the acceleration vector of each element of the finite-element mesh
  do ispec = 1,NSPEC

    tempx1(:,:) = 0.
    tempz1(:,:) = 0.
    tempx2(:,:) = 0.
    tempz2(:,:) = 0.

    !--- elastic spectral element

      ! get elastic parameters of current spectral element
      mu = rho*cs*cs
      lambda = rho*cp*cp - 2.*mu
      lambdaplus2mu = lambda + 2.*mu

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX

          ! derivative along x and along z
          dux_dxi = 0.
          duz_dxi = 0.
          dux_dgamma = 0.
          duz_dgamma = 0.

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displ(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ(2,ibool(i,k,ispec))*hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

          ! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          ! no attenuation
          sigma_xx = lambdaplus2mu*dux_dxl + lambda*duz_dzl
          sigma_xz = mu*(duz_dxl + dux_dzl)
          sigma_zz = lambdaplus2mu*duz_dzl + lambda*dux_dxl
          sigma_zx = sigma_xz

          ! weak formulation term based on stress tensor (non-symmetric form)
          ! also add GLL integration weights
          jacobianl = jacobian(i,j,ispec)

          tempx1(i,j) = wzgll(j)*jacobianl*(sigma_xx*xixl+sigma_zx*xizl) ! this goes to accel_x
          tempz1(i,j) = wzgll(j)*jacobianl*(sigma_xz*xixl+sigma_zz*xizl) ! this goes to accel_z

          tempx2(i,j) = wxgll(i)*jacobianl*(sigma_xx*gammaxl+sigma_zx*gammazl) ! this goes to accel_x
          tempz2(i,j) = wxgll(i)*jacobianl*(sigma_xz*gammaxl+sigma_zz*gammazl) ! this goes to accel_z

        enddo
      enddo  ! end of the loops on the collocation points i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !
      ! second double-loop over GLL to compute all the terms
      !
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          ! along x direction and z direction
          ! and assemble the contributions
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            accel(1,iglob) = accel(1,iglob) - (tempx1(k,j)*hprimewgll_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
            accel(2,iglob) = accel(2,iglob) - (tempz1(k,j)*hprimewgll_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))
          enddo
        enddo
      enddo ! second loop over the GLL points

  enddo   ! end of main loop on all the elements

! add the force source at a given grid point
    iglob = ibool(IGLL_SOURCE,JGLL_SOURCE,NSPEC_SOURCE)
! compute current time
    time = (it-1)*deltat
    accel(2,iglob) = accel(2,iglob) - factor_amplitude * (1.-2.*a*(time-t0)**2) * exp(-a*(time-t0)**2)

! big loop over all the global points (not elements) in the mesh to update
! the acceleration and velocity vectors.
! To compute acceleration from the elastic forces we need to divide them by the mass matrix, i.e. multiply by its inverse
    accel(1,:) = accel(1,:)*rmass_inverse(:)
    accel(2,:) = accel(2,:)*rmass_inverse(:)

    veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)

! record a seismogram to check that the simulation went well
    seismogram(it) = displ(2,ibool(IGLL_RECEIVER,JGLL_RECEIVER,NSPEC_RECEIVER))

  enddo ! end of the serial time loop

! save the seismogram at the end of the run
  open(unit=IIN,file='seismogram.txt',status='unknown')
  do it = 1,NSTEP
    write(IIN,*) (it-1)*deltat,seismogram(it)
  enddo
  close(IIN)

  end program serial_specfem2D

! ******************
! meshing subroutine
! ******************

!--- global node number

  integer function num(i,j,nx)

  implicit none

  integer i,j,nx

    num = j*(nx+1) + i + 1

  end function num

