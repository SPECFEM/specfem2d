!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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
!!!!!!!!!! All the calculations are done in single precision.
!!!!!!!!!! We do not need double precision in SPECFEM2D.
!!!!!!!!!!

! Use an axisymetric formulation or a line source (2D planar)
  logical :: AXISYM = .true.

! density, P wave velocity and S wave velocity of the geophysical medium under study
  real(kind=4), parameter :: rho = 2700.
  real(kind=4), parameter :: cp = 3000.
  real(kind=4), parameter :: cs = cp / 1.732

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

! number of GLJ integration points on the first element (for axisym)
  integer, parameter :: NGLJ = 5

! use 4-node elements to describe the geometry of the mesh
  integer, parameter :: ngnod = 4

! number of spectral elements and of unique grid points of the mesh
  integer, parameter :: NSPEC = nx * nz
  integer, parameter :: NGLOB = (nx*(NGLLX-1) + 1) * (nz*(NGLLZ-1) + 1)

! constant value of the time step in the main time loop, and total number of time steps to simulate
  real, parameter :: deltat = 1.1e-3
  integer, parameter :: NSTEP = 1600

! we locate the source and the receiver in arbitrary elements here for this demo code
  integer :: NSPEC_SOURCE = NSPEC/2 - nx/2, IGLL_SOURCE = 2, JGLL_SOURCE = 2
  integer, parameter :: NSPEC_RECEIVER = 2*NSPEC/3 - nx/4, IGLL_RECEIVER = 2, JGLL_RECEIVER = 2

! for the source time function
  real, parameter :: f0 = 10.
  real, parameter :: t0 = 1.2 / f0
  real, parameter :: factor_amplitude = 1.e+10
  real, parameter :: pi = 3.141592653589793
  real, parameter :: a = pi*pi*f0*f0

  integer, parameter :: NTSTEP_BETWEEN_OUTPUT_INFO = 100

  integer, parameter :: IIN = 40

! 2-D simulation
  integer, parameter :: NDIM = 2

  real(kind=4), parameter :: deltatover2 = 0.5*deltat, deltatsqover2 = 0.5*deltat*deltat

! real(kind=4), parameter :: VERYSMALLVAL = 1.e-24
  real(kind=4), parameter :: ONE = 1.0d0, TINYVAL = 1.0d-24

! displacement threshold above which we consider that the code became unstable
  real(kind=4), parameter :: STABILITY_THRESHOLD = 1.e+25

! global displacement, velocity and acceleration vectors
  real(kind=4), dimension(NDIM,NGLOB) :: displ,veloc,accel

! global diagonal mass matrix
  real(kind=4), dimension(NGLOB) :: rmass_inverse

! record a seismogram to check that the simulation went well
  real(kind=4), dimension(NSTEP) :: seismogram

! time step
  integer it

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLZ,NSPEC) :: ibool
  real(kind=4), dimension(NGLLX,NGLLZ,NSPEC) :: xix,xiz,gammax,gammaz,jacobian

! arrays with the mesh in double precision
  double precision, dimension(NDIM,NGLOB) :: coord

  double precision :: x_source,z_source
  double precision :: x_receiver,z_receiver

! array with derivatives of Lagrange polynomials and precalculated products
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLZ) :: zigll,wzgll
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx,hprime_zz,hprimewgll_xx,hprimewgll_zz

! array with derivatives of Lagrange polynomials for GLJ quadrature (axisym)
  double precision, dimension(NGLLZ) :: xiglj,wxglj
  double precision, dimension(NGLLX,NGLLX) :: hprimeBar_xx,hprimeBarwglj_xx

  real(kind=4), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempx3,tempz1,tempz2

  integer :: ispec,iglob,i,j,k,ix,iz,ipoin
  integer :: nglob_to_compute

  double precision :: xi,gamma,x,z
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=4) dux_dxi,duz_dxi,dux_dgamma,duz_dgamma
  real(kind=4) dux_dxl,dux_dzl,duz_dxl,duz_dzl

  real(kind=4) sigma_xx,sigma_zz,sigma_xz,sigma_zx
  real(kind=4) lambda,mu,lambdaplus2mu

  real(kind=4) Usolidnorm,current_value,time

! parameters and arrays needed for the simple mesh creation routine
  integer, parameter :: npgeo = (nx+1)*(nz+1) ! total number of geometrical points that describe the geometry
  integer, dimension(ngnod,NSPEC) :: knods ! numbering of the four corners of each mesh element
  double precision, dimension(0:nx,0:nz) :: xgrid,zgrid ! coordinates of all the corners of the mesh elements in a first format
  double precision, dimension(NDIM,npgeo) :: coorg ! coordinates of all the corners of the mesh elements in another format
  integer, external :: num ! function that numbers the mesh points with a unique number

! Flag used to determine if an element is on the axis or not (axisym)
  logical, dimension(NSPEC) :: is_on_the_axis

! Array that will contain axial elements indices
  integer, dimension(:), allocatable  ::ispec_of_axial_elements

! Axisym variables
  double precision, dimension(NGLJ,NGLLZ) :: r_xiplus1
  real(kind=4) :: xxi
  real(kind=4), dimension(NGLJ,NGLLZ) :: sigma_thetatheta

! timer to count elapsed time
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  integer ihours,iminutes,iseconds,int_tCPU
  double precision :: time_start,time_end,tCPU

  if (AXISYM) then
    print *,'AXISYM simulation ...'
  endif

! estimate of total memory size used
  print *
  print *,'NSPEC = ',NSPEC
  print *,'NGLOB = ',NGLOB
  print *

  print *,'NSTEP = ',NSTEP
  print *,'deltat = ',deltat
  print *

! set up Gauss-Lobatto-Legendre points, weights and also derivation matrices
  call define_derivation_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,hprimewgll_xx,hprimewgll_zz,NGLLX,NGLLZ)

! If AXISYM also define GLJ points and derivation matrices
  if (AXISYM) then
    call define_GLJ_derivation_matrix(xiglj,wxglj,hprimeBar_xx,hprimeBarwglj_xx,NGLJ)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--- definition of the mesh
  do iz = 0,nz
    do ix = 0,nx
        ! coordinates of the grid points (evenly spaced points along X and Z)
        xgrid(ix,iz) = xmin + (xmax - xmin) * dble(ix) / dble(nx)
        zgrid(ix,iz) = zmin + (zmax - zmin) * dble(iz) / dble(nz)
     enddo
  enddo

! create the coorg array
  do j = 0,nz
    do i = 0,nx
      ipoin = num(i,j,nx)
      coorg(1,ipoin) = xgrid(i,j)
      coorg(2,ipoin) = zgrid(i,j)
    enddo
  enddo

! create the knods array
  k = 0
  do j=0,nz-1
    do i=0,nx-1
      k = k + 1
      knods(1,k) = num(i,j,nx)
      knods(2,k) = num(i+1,j,nx)
      knods(3,k) = num(i+1,j+1,nx)
      knods(4,k) = num(i,j+1,nx)
    enddo
  enddo

!
!---- Build the symmetry axis (ispec_of_axial_elements, is_on_the_axis)
!

  if(AXISYM) then

    if(xmin * xmax < 0) stop 'in axisymmetric mode xmin and xmax must have the same sign, they cannot cross the symmetry axis'
    if(xmin < 0) stop 'in axisymmetric mode, case of symmetry axis on the right edge instead of left not supported yet'

    ! test if the left edge is on the symmetry axis
    if(abs(xmin) < TINYVAL) then

      !all the elements on the left edge are axial because that edge is vertical and located in x = 0

      allocate(ispec_of_axial_elements(nz))
      i = 1

      is_on_the_axis(:) = .false.
      do j = 1,nz
        ispec_of_axial_elements(j) = (j-1)*nx + (i-1) + 1
        is_on_the_axis(ispec_of_axial_elements(j)) = .true.
      enddo

      NSPEC_SOURCE = ispec_of_axial_elements(nz/2)
      IGLL_SOURCE = 1
      JGLL_SOURCE = 2

    else ! no elements on the symmetry axis
      allocate(ispec_of_axial_elements(1))
    endif

  else ! of AXISYM

    allocate(ispec_of_axial_elements(1))

  endif

!
!---- generate the global numbering
!
  call createnum_slow(knods,ibool,nglob_to_compute,nspec,NGLLX,NGLLZ,ngnod) ! Create ibool and recompute nglob for checking
  if(nglob_to_compute /= NGLOB) stop 'error: incorrect total number of unique grid points found'
  if(minval(ibool) /= 1) stop 'error: incorrect minimum value of ibool'
  if(maxval(ibool) /= NGLOB) stop 'error: incorrect maximum value of ibool'

!
!----  set the coordinates of the points of the global grid
!

  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX

        if(AXISYM) then
          if (is_on_the_axis(ispec)) then
            xi = xiglj(i)
          else
            xi = xigll(i)
          endif
        else
          xi = xigll(i)
        endif

        gamma = zigll(j)

        call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl, &
                        jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo,NDIM)
        if(jacobianl <= 0.d0) stop 'error: negative Jacobian found'

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

!
!---- build the mass matrix
!
  rmass_inverse(:) = 0.
  do ispec = 1,NSPEC
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            if (is_on_the_axis(ispec) .and. i == 1) then ! First GLJ point
              xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
              rmass_inverse(iglob) = rmass_inverse(iglob) + xxi*wxglj(i)*wzgll(j)*rho*jacobian(i,j,ispec)
            else
              rmass_inverse(iglob) = rmass_inverse(iglob) + coord(1,iglob)/(xiglj(i)+ONE)*wxglj(i)*wzgll(j)*rho*jacobian(i,j,ispec)
            endif
          else ! not on the axis
            rmass_inverse(iglob) = rmass_inverse(iglob) + coord(1,iglob)*wxgll(i)*wzgll(j)*rho*jacobian(i,j,ispec)
          endif
        else ! not axisym
          rmass_inverse(iglob) = rmass_inverse(iglob) + wxgll(i)*wzgll(j)*rho*jacobian(i,j,ispec)
        endif
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
    if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
      Usolidnorm = -1.
      do iglob = 1,NGLOB
        current_value = sqrt(displ(1,iglob)**2 + displ(2,iglob)**2)
        if(current_value > Usolidnorm) Usolidnorm = current_value
      enddo
      write(*,*) 'Time step # ',it,' out of ',NSTEP
! compute current time
      time = (it-1)*deltat
      write(*,*) 'Max norm displacement vector U in the solid (m) = ',Usolidnorm
! check stability of the code, exit if unstable
      if(Usolidnorm > STABILITY_THRESHOLD .or. Usolidnorm < 0) stop 'code became unstable and blew up'

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

! draw the displacement vector field in a PostScript file
      call plot_post(displ,coord,ibool,NGLOB,NSPEC,x_source,z_source,x_receiver,z_receiver,it,deltat,NGLLX,NGLLZ,NDIM)

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

          if( AXISYM ) then
            if (is_on_the_axis(ispec) ) then
              do k = 1,NGLJ
                dux_dxi = dux_dxi + displ(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                duz_dxi = duz_dxi + displ(2,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                dux_dgamma = dux_dgamma + displ(1,ibool(i,k,ispec))*hprime_zz(j,k)
                duz_dgamma = duz_dgamma + displ(2,ibool(i,k,ispec))*hprime_zz(j,k)
              enddo
            else
              do k = 1,NGLJ
                dux_dxi = dux_dxi + displ(1,ibool(k,j,ispec))*hprime_xx(i,k)
                duz_dxi = duz_dxi + displ(2,ibool(k,j,ispec))*hprime_xx(i,k)
                dux_dgamma = dux_dgamma + displ(1,ibool(i,k,ispec))*hprime_zz(j,k)
                duz_dgamma = duz_dgamma + displ(2,ibool(i,k,ispec))*hprime_zz(j,k)
              enddo
            endif
          else
            do k = 1,NGLLX
              dux_dxi = dux_dxi + displ(1,ibool(k,j,ispec))*hprime_xx(i,k)
              duz_dxi = duz_dxi + displ(2,ibool(k,j,ispec))*hprime_xx(i,k)
              dux_dgamma = dux_dgamma + displ(1,ibool(i,k,ispec))*hprime_zz(j,k)
              duz_dgamma = duz_dgamma + displ(2,ibool(i,k,ispec))*hprime_zz(j,k)
            enddo
          endif

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

          ! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          if (AXISYM) then
            if(is_on_the_axis(ispec) .and. i == 1) then ! d_uz/dr=0 on the axis
              duz_dxl = 0.d0
            endif
          endif

          ! no attenuation

          if( AXISYM ) then
            if (is_on_the_axis(ispec) ) then
              if (is_on_the_axis(ispec) .and. i == 1) then ! First GLJ point
                sigma_xx = 0.0d0
                sigma_zz = 0.0d0
                sigma_thetatheta(i,j) = 0.0d0
                xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                r_xiplus1(i,j) = xxi
                do k = 1,NGLJ
                  sigma_xx = sigma_xx + displ(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  sigma_zz = sigma_zz + displ(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  sigma_thetatheta(i,j) = sigma_thetatheta(i,j) + displ(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                enddo
                sigma_xx = lambdaplus2mu*dux_dxl + lambda*duz_dzl &
                           + lambda*sigma_xx/xxi
                sigma_zz = lambdaplus2mu*duz_dzl + lambda*dux_dxl &
                           + lambda*sigma_zz/xxi
                sigma_xz = mu*(duz_dxl + dux_dzl)
                sigma_zx = sigma_xz
                sigma_thetatheta(i,j) = lambda*duz_dzl + lambda*dux_dxl &
                                   + lambdaplus2mu*sigma_thetatheta(i,j)/xxi
              else ! Not first GLJ point
                sigma_xx = lambdaplus2mu*dux_dxl + lambda*duz_dzl &
                           + lambda*displ(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_zz = lambdaplus2mu*duz_dzl + lambda*dux_dxl &
                           + lambda*displ(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_xz = mu*(duz_dxl + dux_dzl)
                sigma_zx = sigma_xz
                sigma_thetatheta(i,j) = lambda*duz_dzl + lambda*dux_dxl &
                                        + lambdaplus2mu &
                                        * displ(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
              endif
            else ! Not on the axis
              sigma_xx = lambdaplus2mu*dux_dxl + lambda*duz_dzl &
                         + lambda*displ(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              sigma_zz = lambdaplus2mu*duz_dzl + lambda*dux_dxl &
                         + lambda*displ(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              sigma_xz = mu*(duz_dxl + dux_dzl)
              sigma_zx = sigma_xz
              sigma_thetatheta(i,j) = lambda*duz_dzl + lambda*dux_dxl &
                                      + lambdaplus2mu &
                                      * displ(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
            endif
          else ! Not axisym
            sigma_xx = lambdaplus2mu*dux_dxl + lambda*duz_dzl
            sigma_xz = mu*(duz_dxl + dux_dzl)
            sigma_zz = lambdaplus2mu*duz_dzl + lambda*dux_dxl
            sigma_zx = sigma_xz
          endif

          ! weak formulation term based on stress tensor (non-symmetric form)
          ! also add GLL integration weights
          jacobianl = jacobian(i,j,ispec)

          if (AXISYM) then
            if (is_on_the_axis(ispec) ) then
              ! This is normal, we always add a contribution depending on the value on the axis
              ! i.e. we purposely sum something at point (i,j) with something at point (1,j)
              tempx3(i,j) = wzgll(j)*jacobian(1,j,ispec)*sigma_thetatheta(1,j)*hprimeBarwglj_xx(1,i)

              ! not first GLJ point
              if ( abs(coord(1,ibool(i,j,ispec))) > TINYVAL ) then
                if ( i == 1 ) then
                  write(*,*) "Element number:",ispec
                  stop 'error: an axial element is rotated. The code should have been stopped before. Check that your &
                   &coordinates are >> TINYVAL. Maybe you should also have a look to &
                   &doc/problematic_case_that_we_exclude_for_axisymmetric.pdf '
                endif
                tempx3(i,j) = tempx3(i,j) + wzgll(j)*wxglj(i)*jacobian(i,j,ispec) &
                              * sigma_thetatheta(i,j)/(xiglj(i)+ONE) ! this goes to accel_x
              endif

              tempx2(i,j) = r_xiplus1(i,j)*wxglj(i)*jacobianl &
                            * (sigma_xx*gammaxl+sigma_zx*gammazl) ! this goes to accel_x
              tempz2(i,j) = r_xiplus1(i,j)*wxglj(i)*jacobianl &
                            * (sigma_xz*gammaxl+sigma_zz*gammazl) ! this goes to accel_z
              tempx1(i,j) = r_xiplus1(i,j)*wzgll(j)*jacobianl &
                            * (sigma_xx*xixl+sigma_zx*xizl) ! this goes to accel_x
              tempz1(i,j) = r_xiplus1(i,j)*wzgll(j)*jacobianl &
                            * (sigma_xz*xixl+sigma_zz*xizl) ! this goes to accel_z
            else ! axisym but not on the axis
              tempx2(i,j) = coord(1,ibool(i,j,ispec))*wxgll(i)*jacobianl &
                            *(sigma_xx*gammaxl+sigma_zx*gammazl) ! this goes to accel_x
              tempz2(i,j) = coord(1,ibool(i,j,ispec))*wxgll(i)*jacobianl &
                            *(sigma_xz*gammaxl+sigma_zz*gammazl) ! this goes to accel_z
              tempx1(i,j) = coord(1,ibool(i,j,ispec))*wzgll(j)*jacobianl &
                            *(sigma_xx*xixl+sigma_zx*xizl) ! this goes to accel_x
              tempz1(i,j) = coord(1,ibool(i,j,ispec))*wzgll(j)*jacobianl &
                            *(sigma_xz*xixl+sigma_zz*xizl) ! this goes to accel_z
              tempx3(i,j) = wxgll(i)*wzgll(j)*jacobianl*sigma_thetatheta(i,j) ! this goes to accel_x
            endif
          else

            tempx1(i,j) = wzgll(j)*jacobianl*(sigma_xx*xixl+sigma_zx*xizl) ! this goes to accel_x
            tempz1(i,j) = wzgll(j)*jacobianl*(sigma_xz*xixl+sigma_zz*xizl) ! this goes to accel_z

            tempx2(i,j) = wxgll(i)*jacobianl*(sigma_xx*gammaxl+sigma_zx*gammazl) ! this goes to accel_x
            tempz2(i,j) = wxgll(i)*jacobianl*(sigma_xz*gammaxl+sigma_zz*gammazl) ! this goes to accel_z

          endif

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
          if (AXISYM ) then
            if (is_on_the_axis(ispec) ) then
              do k = 1,NGLJ
                accel(1,iglob) = accel(1,iglob) &
                                         - (tempx1(k,j)*hprimeBarwglj_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
                accel(2,iglob) = accel(2,iglob) &
                                         - (tempz1(k,j)*hprimeBarwglj_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))
              enddo
              accel(1,iglob) = accel(1,iglob) - tempx3(i,j)
            else ! Axisym but not on the axis
              do k = 1,NGLLX
                accel(1,iglob) = accel(1,iglob) &
                                         - (tempx1(k,j)*hprimewgll_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
                accel(2,iglob) = accel(2,iglob) &
                                         - (tempz1(k,j)*hprimewgll_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))
              enddo
              accel(1,iglob) = accel(1,iglob) - tempx3(i,j)
            endif
          else !if AXISYM == false
            do k = 1,NGLLX
              accel(1,iglob) = accel(1,iglob) - (tempx1(k,j)*hprimewgll_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
              accel(2,iglob) = accel(2,iglob) - (tempz1(k,j)*hprimewgll_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))
            enddo
          endif
        enddo
      enddo ! second loop over the GLL points

  enddo   ! end of main loop on all the elements

! add the force source at a given grid point
    iglob = ibool(IGLL_SOURCE,JGLL_SOURCE,NSPEC_SOURCE)
! compute current time
    time = (it-1)*deltat
    accel(2,iglob) = accel(2,iglob) - factor_amplitude * (1.-2.*a*(time-t0)**2) * exp(-a*(time-t0)**2)

! big loop over all the global points (not elements) in the mesh to update
! the acceleration and velocity vectors
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

