
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) January 2005
!
!========================================================================

!====================================================================================
!
! An explicit 2D spectral element solver for the anelastic anisotropic wave equation
!
!====================================================================================

!
! version 5.1, January 2005 :
!               - more general mesher with any number of curved layers
!               - Dirac and Gaussian time sources and corresponding convolution routine
!               - option for acoustic medium instead of elastic
!               - receivers at any location, not only grid points
!               - moment-tensor source at any location, not only a grid point
!               - color PNM snapshots
!               - more flexible DATA/Par_file with any number of comment lines
!               - Xsu scripts for seismograms
!               - subtract t0 from seismograms
!
! version 5.0, May 2004 :
!               - got rid of useless routines, suppressed commons etc.
!               - weak formulation based explicitly on stress tensor
!               - implementation of full anisotropy
!               - implementation of attenuation based on memory variables
!
! based on SPECFEM2D version 4.2, June 1998
! (c) by Dimitri Komatitsch, Harvard University, USA
! and Jean-Pierre Vilotte, Institut de Physique du Globe de Paris, France
!
! itself based on SPECFEM2D version 1.0, 1995
! (c) by Dimitri Komatitsch and Jean-Pierre Vilotte,
! Institut de Physique du Globe de Paris, France
!

  program specfem2D

  implicit none

  include "constants.h"

  character(len=80) datlin

  integer source_type,time_function_type
  double precision x_source,z_source,xi_source,gamma_source,Mxx,Mzz,Mxz,f0,t0,factor,angleforce
  double precision, dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  double precision, dimension(:,:), allocatable :: coorg
  double precision, dimension(:), allocatable :: coorgread

! receiver information
  integer, dimension(:), allocatable :: ispec_selected_rec
  double precision, dimension(:), allocatable :: xi_receiver,gamma_receiver,st_xval,st_zval

  double precision, dimension(:,:), allocatable :: sisux,sisuz

  logical anyabs

  integer i,j,it,irec,ipoin,ip,id
  integer nbpoin,inump,n,npoinext,ispec,npoin,npgeo,iglob

  double precision dxd,dzd,valux,valuz,hlagrange,xdummy,zdummy
  double precision cpl,csl,rhol
  double precision cosrot,sinrot
  double precision xi,gamma,x,z

! coefficients of the explicit Newmark time scheme
  integer NSTEP
  double precision deltatover2,deltatsquareover2,time,deltat

! Gauss-Lobatto-Legendre points and weights
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLZ) :: zigll,wzgll

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz

! space derivatives
  double precision tempx1l,tempx2l,tempz1l,tempz2l
  double precision fac1,fac2,hp1,hp2
  double precision duxdxl,duzdxl,duxdzl,duzdzl
  double precision sigma_xx,sigma_xz,sigma_zx,sigma_zz

  double precision, dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempz1,tempz2

! for anisotropy
  double precision duydyl,duydzl,duzdyl,duxdyl,duydxl
  double precision duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  double precision duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

! Jacobian matrix and determinant
  double precision xixl,xizl,gammaxl,gammazl,jacobianl

! material properties of the elastic medium
  double precision mul_relaxed,lambdal_relaxed,lambdalplus2mul_relaxed,cpsquare
  double precision mul_unrelaxed,lambdal_unrelaxed,lambdalplus2mul_unrelaxed

  double precision, dimension(:,:), allocatable :: coord,accel,veloc,displ, &
    flagrange,xinterp,zinterp,Uxinterp,Uzinterp,elastcoef,vector_field_postscript

  double precision, dimension(:), allocatable :: rmass,density,vpext,vsext,rhoext,displread,velocread,accelread

  double precision, dimension(:,:,:), allocatable :: shape2D,shape2D_display,xix,xiz,gammax,gammaz,jacobian

  double precision, dimension(:,:,:,:), allocatable :: dershape2D,dershape2D_display

  integer, dimension(:,:,:), allocatable :: ibool
  integer, dimension(:,:), allocatable  :: knods
  integer, dimension(:), allocatable :: kmato,numabs,numsurface

  integer ie,k

  integer ispec_selected_source,iglob_source,ix_source,iz_source
  double precision a,displnorm_all
  double precision, dimension(:), allocatable :: source_time_function

  double precision rsizemin,rsizemax,cpoverdxmin,cpoverdxmax, &
    lambdaSmin,lambdaSmax,lambdaPmin,lambdaPmax,vpmin,vpmax

  integer colors,numbers,subsamp,vecttype,itaff,nrec,sismostype
  integer numat,ngnod,nspec,iptsdisp,nelemabs,nelemsurface

  logical interpol,meshvect,modelvect,boundvect,read_external_model,initialfield,abshaut, &
    outputgrid,gnuplot,ELASTIC,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON,output_postscript_image,output_PNM_image

  double precision cutvect,anglerec,xirec,gammarec

! for absorbing and free surface conditions
  integer ispecabs,ispecsurface,inum,numabsread,numsurfaceread,i1abs,i2abs
  logical codeabsread(4)
  double precision nx,nz,vx,vz,vn,rho_vp,rho_vs,tx,tz,weight,xxi,zgamma,kappal

  logical, dimension(:,:), allocatable  :: codeabs

! for attenuation
  integer nspec_allocate

  double precision tau_epsilon_nu1_mech1,tau_sigma_nu1_mech1, &
    tau_epsilon_nu2_mech1,tau_sigma_nu2_mech1,tau_epsilon_nu1_mech2, &
    tau_sigma_nu1_mech2,tau_epsilon_nu2_mech2,tau_sigma_nu2_mech2

  double precision Un,Unp1,tauinv,Sn,Snp1,theta_n,theta_np1
  double precision phi_nu1_mech1,phi_nu2_mech1,phi_nu1_mech2,phi_nu2_mech2
  double precision deltatsquare,deltatcube,deltatfourth
  double precision twelvedeltat,fourdeltatsquare
  double precision tauinvsquare,tauinvcube,tauinvUn
  double precision inv_tau_sigma_nu1_mech1,inv_tau_sigma_nu2_mech1
  double precision inv_tau_sigma_nu1_mech2,inv_tau_sigma_nu2_mech2

  double precision Mu_nu1,Mu_nu2

  double precision, dimension(:,:,:), allocatable :: &
    e1_mech1,e11_mech1,e13_mech1,e1_mech2,e11_mech2,e13_mech2, &
    duxdxl_n,duzdzl_n,duzdxl_n,duxdzl_n,duxdxl_np1,duzdzl_np1,duzdxl_np1,duxdzl_np1

! for color PNM images
  integer :: NX_IMAGE_PNM,NZ_IMAGE_PNM,iplus1,jplus1,iminus1,jminus1,nx_sem_PNM
  double precision :: xmin_PNM_image,xmax_PNM_image,zmin_PNM_image,zmax_PNM_image,taille_pixel_horizontal,taille_pixel_vertical
  integer, dimension(:,:), allocatable :: iglob_image_PNM_2D,copy_iglob_image_PNM_2D
  double precision, dimension(:,:), allocatable :: donnees_image_PNM_2D

! timing information for the stations
  character(len=MAX_LENGTH_STATION_NAME), allocatable, dimension(:) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), allocatable, dimension(:) :: network_name

! title of the plot
  character(len=60) stitle

! Lagrange interpolators at receivers
  double precision, dimension(:), allocatable :: hxir,hgammar,hpxir,hpgammar
  double precision, dimension(:,:), allocatable :: hxir_store,hgammar_store

! for Lagrange interpolants
  double precision, external :: hgll

!***********************************************************************
!
!             i n i t i a l i z a t i o n    p h a s e
!
!***********************************************************************

  open (IIN,file='OUTPUT_FILES/Database')

! uncomment this to write to file instead of standard output
! open (IOUT,file='results_simulation.txt')

!
!---  read job title and skip remaining titles of the input file
!
  read(IIN,40) datlin
  read(IIN,40) datlin
  read(IIN,40) datlin
  read(IIN,40) datlin
  read(IIN,40) datlin
  read(IIN,45) stitle

!
!---- print the date, time and start-up banner
!
  call datim(stitle)

  write(*,*)
  write(*,*)
  write(*,*) '*********************'
  write(*,*) '****             ****'
  write(*,*) '****  SPECFEM2D  ****'
  write(*,*) '****             ****'
  write(*,*) '*********************'

!
!---- read parameters from input file
!

  read(IIN,40) datlin
  read(IIN,*) npgeo

  read(IIN,40) datlin
  read(IIN,*) gnuplot,interpol

  read(IIN,40) datlin
  read(IIN,*) itaff,output_postscript_image,output_PNM_image,colors,numbers

  read(IIN,40) datlin
  read(IIN,*) meshvect,modelvect,boundvect,cutvect,subsamp,nx_sem_PNM
  cutvect = cutvect / 100.d0

  read(IIN,40) datlin
  read(IIN,*) anglerec

  read(IIN,40) datlin
  read(IIN,*) initialfield

  read(IIN,40) datlin
  read(IIN,*) sismostype,vecttype

  read(IIN,40) datlin
  read(IIN,*) read_external_model,outputgrid,ELASTIC,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON

!---- check parameters read
  write(IOUT,200) npgeo,NDIM
  write(IOUT,600) itaff,colors,numbers
  write(IOUT,700) sismostype,anglerec
  write(IOUT,750) initialfield,read_external_model,ELASTIC,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON,outputgrid
  write(IOUT,800) vecttype,100.d0*cutvect,subsamp

!---- read time step
  read(IIN,40) datlin
  read(IIN,*) NSTEP,deltat
  write(IOUT,703) NSTEP,deltat,NSTEP*deltat

!
!----  read source information
!
  read(IIN,40) datlin
  read(IIN,*) source_type,time_function_type,x_source,z_source,f0,t0,factor,angleforce,Mxx,Mzz,Mxz

!
!-----  check the input
!
 if(.not. initialfield) then
   if (source_type == 1) then
     write(IOUT,212) x_source,z_source,f0,t0,factor,angleforce
   else if(source_type == 2) then
     write(IOUT,222) x_source,z_source,f0,t0,factor,Mxx,Mzz,Mxz
   else
     stop 'Unknown source type number !'
   endif
 endif

! for the source time function
  a = pi*pi*f0*f0

!-----  convert angle from degrees to radians
  angleforce = angleforce * pi / 180.d0

!
!---- read the spectral macrobloc nodal coordinates
!
  allocate(coorg(NDIM,npgeo))

  ipoin = 0
  read(IIN,40) datlin
  allocate(coorgread(NDIM))
  do ip = 1,npgeo
   read(IIN,*) ipoin,(coorgread(id),id =1,NDIM)
   if(ipoin<1 .or. ipoin>npgeo) stop 'Wrong control point number'
   coorg(:,ipoin) = coorgread
  enddo
  deallocate(coorgread)

!
!---- read the basic properties of the spectral elements
!
  read(IIN,40) datlin
  read(IIN,*) numat,ngnod,nspec,iptsdisp,nelemabs,nelemsurface

!
!---- allocate arrays
!
  allocate(shape2D(ngnod,NGLLX,NGLLZ))
  allocate(dershape2D(NDIM,ngnod,NGLLX,NGLLZ))
  allocate(shape2D_display(ngnod,iptsdisp,iptsdisp))
  allocate(dershape2D_display(NDIM,ngnod,iptsdisp,iptsdisp))
  allocate(xix(NGLLX,NGLLZ,nspec))
  allocate(xiz(NGLLX,NGLLZ,nspec))
  allocate(gammax(NGLLX,NGLLZ,nspec))
  allocate(gammaz(NGLLX,NGLLZ,nspec))
  allocate(jacobian(NGLLX,NGLLZ,nspec))
  allocate(flagrange(NGLLX,iptsdisp))
  allocate(xinterp(iptsdisp,iptsdisp))
  allocate(zinterp(iptsdisp,iptsdisp))
  allocate(Uxinterp(iptsdisp,iptsdisp))
  allocate(Uzinterp(iptsdisp,iptsdisp))
  allocate(density(numat))
  allocate(elastcoef(4,numat))
  allocate(kmato(nspec))
  allocate(knods(ngnod,nspec))
  allocate(ibool(NGLLX,NGLLZ,nspec))

  if(read_external_model .and. outputgrid) stop 'cannot output the grid and read external model at the same time'

! for acoustic
  if(TURN_ANISOTROPY_ON .and. .not. ELASTIC) stop 'currently cannot have anisotropy in acoustic simulation'

  if(TURN_ATTENUATION_ON .and. .not. ELASTIC) stop 'currently cannot have attenuation in acoustic simulation'

  if(source_type == 2 .and. .not. ELASTIC) stop 'currently cannot have moment tensor source in acoustic simulation'

! for attenuation
  if(TURN_ANISOTROPY_ON .and. TURN_ATTENUATION_ON) stop 'cannot have anisotropy and attenuation both turned on in current version'

  if(TURN_ATTENUATION_ON) then
    nspec_allocate = nspec
  else
    nspec_allocate = 1
  endif

  allocate(e1_mech1(NGLLX,NGLLZ,nspec_allocate))
  allocate(e11_mech1(NGLLX,NGLLZ,nspec_allocate))
  allocate(e13_mech1(NGLLX,NGLLZ,nspec_allocate))
  allocate(e1_mech2(NGLLX,NGLLZ,nspec_allocate))
  allocate(e11_mech2(NGLLX,NGLLZ,nspec_allocate))
  allocate(e13_mech2(NGLLX,NGLLZ,nspec_allocate))
  allocate(duxdxl_n(NGLLX,NGLLZ,nspec_allocate))
  allocate(duzdzl_n(NGLLX,NGLLZ,nspec_allocate))
  allocate(duzdxl_n(NGLLX,NGLLZ,nspec_allocate))
  allocate(duxdzl_n(NGLLX,NGLLZ,nspec_allocate))
  allocate(duxdxl_np1(NGLLX,NGLLZ,nspec_allocate))
  allocate(duzdzl_np1(NGLLX,NGLLZ,nspec_allocate))
  allocate(duzdxl_np1(NGLLX,NGLLZ,nspec_allocate))
  allocate(duxdzl_np1(NGLLX,NGLLZ,nspec_allocate))

! --- allocate arrays for absorbing boundary conditions
  if(nelemabs <= 0) then
    nelemabs = 1
    anyabs = .false.
  else
    anyabs = .true.
  endif
  allocate(numabs(nelemabs))
  allocate(codeabs(4,nelemabs))

! --- allocate array for free surface condition in acoustic medium
  if(nelemsurface <= 0) nelemsurface = 1
  allocate(numsurface(nelemsurface))

!
!---- print element group main parameters
!
  write(IOUT,107)
  write(IOUT,207) nspec,ngnod,NGLLX,NGLLZ,NGLLX*NGLLZ,iptsdisp,numat,nelemabs

! set up Gauss-Lobatto-Legendre derivation matrices
  call define_derivative_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz)

!
!---- read the material properties
!
  call gmat01(density,elastcoef,numat)

!
!----  read spectral macrobloc data
!
  n = 0
  read(IIN,40) datlin
  do ie = 1,nspec
    read(IIN,*) n,kmato(n),(knods(k,n), k=1,ngnod)
  enddo

!
!----  read absorbing boundary data
!
  if(anyabs) then
    read(IIN,40) datlin
    do n=1,nelemabs
      read(IIN,*) inum,numabsread,codeabsread(1),codeabsread(2),codeabsread(3),codeabsread(4)
      if(inum < 1 .or. inum > nelemabs) stop 'Wrong absorbing element number'
      numabs(inum) = numabsread
      codeabs(ITOP,inum) = codeabsread(1)
      codeabs(IBOTTOM,inum) = codeabsread(2)
      codeabs(ILEFT,inum) = codeabsread(3)
      codeabs(IRIGHT,inum) = codeabsread(4)
    enddo
    write(*,*)
    write(*,*) 'Number of absorbing elements: ',nelemabs
  endif

!
!----  read free surface data
!
  read(IIN,40) datlin
  read(IIN,*) abshaut
  do n=1,nelemsurface
    read(IIN,*) inum,numsurfaceread
    if(inum < 1 .or. inum > nelemsurface) stop 'Wrong free surface element number'
    numsurface(inum) = numsurfaceread
  enddo
  write(*,*)
  write(*,*) 'Number of free surface elements: ',nelemsurface

!
!---- close input file
!
  close(IIN)

!
!---- compute shape functions and their derivatives for SEM grid
!
  do j = 1,NGLLZ
    do i = 1,NGLLX
      call define_shape_functions(shape2D(:,i,j),dershape2D(:,:,i,j),xigll(i),zigll(j),ngnod)
    enddo
  enddo

!
!---- generate the global numbering
!

! version "propre mais lente" ou version "sale mais rapide"
  if(FAST_NUMBERING) then
    call createnum_fast(knods,ibool,shape2D,coorg,npoin,npgeo,nspec,ngnod)
  else
    call createnum_slow(knods,ibool,npoin,nspec,ngnod)
  endif

!---- compute shape functions and their derivatives for regular !interpolated display grid
  do j = 1,iptsdisp
    do i = 1,iptsdisp
      xirec  = 2.d0*dble(i-1)/dble(iptsdisp-1) - 1.d0
      gammarec  = 2.d0*dble(j-1)/dble(iptsdisp-1) - 1.d0
      call define_shape_functions(shape2D_display(:,i,j),dershape2D_display(:,:,i,j),xirec,gammarec,ngnod)
    enddo
  enddo

!---- compute Lagrange interpolants on a regular interpolated grid in (xi,gamma)
!---- for display (assumes NGLLX = NGLLZ)
  do j=1,NGLLX
    do i=1,iptsdisp
      xirec  = 2.d0*dble(i-1)/dble(iptsdisp-1) - 1.d0
      flagrange(j,i) = hgll(j-1,xirec,xigll,NGLLX)
    enddo
  enddo


! read total number of receivers
  open(unit=IIN,file='DATA/STATIONS',status='old')
  read(IIN,*) nrec
  close(IIN)

  write(IOUT,*)
  write(IOUT,*) 'Total number of receivers = ',nrec
  write(IOUT,*)

  if(nrec < 1) stop 'need at least one receiver'

! allocate seismogram arrays
  allocate(sisux(NSTEP,nrec))
  allocate(sisuz(NSTEP,nrec))

! receiver information
  allocate(ispec_selected_rec(nrec))
  allocate(st_xval(nrec))
  allocate(st_zval(nrec))
  allocate(xi_receiver(nrec))
  allocate(gamma_receiver(nrec))
  allocate(station_name(nrec))
  allocate(network_name(nrec))

! allocate 1-D Lagrange interpolators and derivatives
  allocate(hxir(NGLLX))
  allocate(hpxir(NGLLX))
  allocate(hgammar(NGLLZ))
  allocate(hpgammar(NGLLZ))

! allocate Lagrange interpolators for receivers
  allocate(hxir_store(nrec,NGLLX))
  allocate(hgammar_store(nrec,NGLLZ))

!
!----  allocation des autres tableaux pour la grille globale et les bords
!

  allocate(coord(NDIM,npoin))

  allocate(accel(NDIM,npoin))
  allocate(displ(NDIM,npoin))
  allocate(veloc(NDIM,npoin))

! for acoustic medium
  if(ELASTIC) then
    allocate(vector_field_postscript(NDIM,1))
  else
    allocate(vector_field_postscript(NDIM,npoin))
  endif

  allocate(rmass(npoin))

  if(read_external_model) then
    npoinext = npoin
  else
    npoinext = 1
  endif
  allocate(vpext(npoinext))
  allocate(vsext(npoinext))
  allocate(rhoext(npoinext))

!
!----  set the coordinates of the points of the global grid
!
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX

        xi = xigll(i)
        gamma = zigll(j)

        call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl,jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo)

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
!--- save the grid of points in a file
!
  if(outputgrid) then
    print *
    print *,'Saving the grid in a text file...'
    print *
    open(unit=55,file='OUTPUT_FILES/grid_points_and_model.txt',status='unknown')
    write(55,*) npoin
    do n = 1,npoin
      write(55,*) (coord(i,n), i=1,NDIM)
    enddo
    close(55)
  endif

!
!-----   plot the GLL mesh in a Gnuplot file
!
  if(gnuplot) call plotgll(knods,ibool,coorg,coord,npoin,npgeo,ngnod,nspec)

!
!----   define coefficients of the Newmark time scheme
!
  deltatover2 = HALF*deltat
  deltatsquareover2 = HALF*deltat*deltat

!---- definir la position reelle des points source et recepteurs
  if(source_type == 1) then
! collocated force source
    call locate_source_force(coord,ibool,npoin,nspec,x_source,z_source,source_type, &
      ix_source,iz_source,ispec_selected_source,iglob_source)

  else if(source_type == 2) then
! moment-tensor source
    call locate_source_moment_tensor(ibool,coord,nspec,npoin,xigll,zigll,x_source,z_source, &
               ispec_selected_source,xi_source,gamma_source,coorg,knods,ngnod,npgeo)

! compute source array for moment-tensor source
    call compute_arrays_source(ispec_selected_source,xi_source,gamma_source,sourcearray, &
               Mxx,Mzz,Mxz,xix,xiz,gammax,gammaz,xigll,zigll,nspec)

  else
    stop 'incorrect source type'
  endif


! locate receivers in the mesh
  call locate_receivers(ibool,coord,nspec,npoin,xigll,zigll,nrec,st_xval,st_zval,ispec_selected_rec, &
                 xi_receiver,gamma_receiver,station_name,network_name,x_source,z_source,coorg,knods,ngnod,npgeo)

! define and store Lagrange interpolators at all the receivers
  do irec = 1,nrec
    call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
    call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
    hxir_store(irec,:) = hxir(:)
    hgammar_store(irec,:) = hgammar(:)
  enddo

!
!----  eventuellement lecture d'un modele externe de vitesse et de densite
!
  if(read_external_model) then
    print *
    print *,'Reading velocity and density model from external file...'
    print *
    open(unit=55,file='OUTPUT_FILES/grid_points_and_model.txt',status='old')
    read(55,*) nbpoin
    if(nbpoin /= npoin) stop 'Wrong number of points in input file'
    do n = 1,npoin
      read(55,*) xdummy,zdummy,rhoext(n),vpext(n),vsext(n)
    enddo
    close(55)
  endif

!
!---- define all arrays
!
  call defarrays(vpext,vsext,rhoext,density,elastcoef, &
          ibool,kmato,coord,npoin,rsizemin,rsizemax, &
          cpoverdxmin,cpoverdxmax,lambdaSmin,lambdaSmax,lambdaPmin,lambdaPmax, &
          vpmin,vpmax,read_external_model,nspec,numat)

! build the global mass matrix once and for all
  rmass(:) = ZERO
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
! if external density model
        if(read_external_model) then
          rhol = rhoext(iglob)
          cpsquare = vpext(iglob)**2
        else
          rhol = density(kmato(ispec))
          lambdal_relaxed = elastcoef(1,kmato(ispec))
          mul_relaxed = elastcoef(2,kmato(ispec))
          cpsquare = (lambdal_relaxed + 2.d0*mul_relaxed) / rhol
        endif
! for acoustic medium
        if(ELASTIC) then
          rmass(iglob) = rmass(iglob) + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)
        else
          rmass(iglob) = rmass(iglob) + wxgll(i)*wzgll(j)*jacobian(i,j,ispec) / cpsquare
        endif
      enddo
    enddo
  enddo

! convertir angle recepteurs en radians
  anglerec = anglerec * pi / 180.d0

!---- verifier le maillage, la stabilite et le nb de points par lambda
!---- seulement si la source en temps n'est pas un Dirac (sinon spectre non defini)
  if(time_function_type /= 4) call checkgrid(deltat,f0,t0,initialfield, &
      rsizemin,rsizemax,cpoverdxmin,cpoverdxmax,lambdaSmin,lambdaSmax,lambdaPmin,lambdaPmax)

!
!---- for color PNM images
!

! taille horizontale de l'image
  xmin_PNM_image = minval(coord(1,:))
  xmax_PNM_image = maxval(coord(1,:))

! taille verticale de l'image, augmenter un peu pour depasser de la topographie
  zmin_PNM_image = minval(coord(2,:))
  zmax_PNM_image = maxval(coord(2,:))
  zmax_PNM_image = zmin_PNM_image + 1.05d0 * (zmax_PNM_image - zmin_PNM_image)

! calculer le nombre de pixels en horizontal en fonction du nombre d'elements spectraux
  NX_IMAGE_PNM = nx_sem_PNM * (NGLLX-1) + 1

! calculer le nombre de pixels en vertical en fonction du rapport des tailles
  NZ_IMAGE_PNM = nint(NX_IMAGE_PNM * (zmax_PNM_image - zmin_PNM_image) / (xmax_PNM_image - xmin_PNM_image))

! allouer un tableau pour les donnees de l'image
  allocate(donnees_image_PNM_2D(NX_IMAGE_PNM,NZ_IMAGE_PNM))

! allouer un tableau pour le point de grille contenant cette donnee
  allocate(iglob_image_PNM_2D(NX_IMAGE_PNM,NZ_IMAGE_PNM))
  allocate(copy_iglob_image_PNM_2D(NX_IMAGE_PNM,NZ_IMAGE_PNM))

! creer tous les pixels
  print *
  print *,'localisation de tous les pixels des images PNM'

  taille_pixel_horizontal = (xmax_PNM_image - xmin_PNM_image) / dble(NX_IMAGE_PNM-1)
  taille_pixel_vertical = (zmax_PNM_image - zmin_PNM_image) / dble(NZ_IMAGE_PNM-1)

  iglob_image_PNM_2D(:,:) = -1

! boucle sur tous les points de grille pour leur affecter un pixel de l'image
      do n=1,npoin

! calculer les coordonnees du pixel
      i = nint((coord(1,n) - xmin_PNM_image) / taille_pixel_horizontal + 1)
      j = nint((coord(2,n) - zmin_PNM_image) / taille_pixel_vertical + 1)

! eviter les effets de bord
      if(i < 1) i = 1
      if(i > NX_IMAGE_PNM) i = NX_IMAGE_PNM

      if(j < 1) j = 1
      if(j > NZ_IMAGE_PNM) j = NZ_IMAGE_PNM

! affecter ce point a ce pixel
      iglob_image_PNM_2D(i,j) = n

      enddo

! completer les pixels manquants en les localisant par la distance minimum
  copy_iglob_image_PNM_2D(:,:) = iglob_image_PNM_2D(:,:)

  do j = 1,NZ_IMAGE_PNM
    do i = 1,NX_IMAGE_PNM

      if(copy_iglob_image_PNM_2D(i,j) == -1) then

        iplus1 = i + 1
        iminus1 = i - 1

        jplus1 = j + 1
        jminus1 = j - 1

! eviter les effets de bord
        if(iminus1 < 1) iminus1 = 1
        if(iplus1 > NX_IMAGE_PNM) iplus1 = NX_IMAGE_PNM

        if(jminus1 < 1) jminus1 = 1
        if(jplus1 > NZ_IMAGE_PNM) jplus1 = NZ_IMAGE_PNM

! utiliser les pixels voisins pour remplir les trous

! horizontales
        if(copy_iglob_image_PNM_2D(iplus1,j) /= -1) then
          iglob_image_PNM_2D(i,j) = copy_iglob_image_PNM_2D(iplus1,j)

        else if(copy_iglob_image_PNM_2D(iminus1,j) /= -1) then
          iglob_image_PNM_2D(i,j) = copy_iglob_image_PNM_2D(iminus1,j)

! verticales
        else if(copy_iglob_image_PNM_2D(i,jplus1) /= -1) then
          iglob_image_PNM_2D(i,j) = copy_iglob_image_PNM_2D(i,jplus1)

        else if(copy_iglob_image_PNM_2D(i,jminus1) /= -1) then
          iglob_image_PNM_2D(i,j) = copy_iglob_image_PNM_2D(i,jminus1)

! diagonales
        else if(copy_iglob_image_PNM_2D(iminus1,jminus1) /= -1) then
          iglob_image_PNM_2D(i,j) = copy_iglob_image_PNM_2D(iminus1,jminus1)

        else if(copy_iglob_image_PNM_2D(iplus1,jminus1) /= -1) then
          iglob_image_PNM_2D(i,j) = copy_iglob_image_PNM_2D(iplus1,jminus1)

        else if(copy_iglob_image_PNM_2D(iminus1,jplus1) /= -1) then
          iglob_image_PNM_2D(i,j) = copy_iglob_image_PNM_2D(iminus1,jplus1)

        else if(copy_iglob_image_PNM_2D(iplus1,jplus1) /= -1) then
          iglob_image_PNM_2D(i,j) = copy_iglob_image_PNM_2D(iplus1,jplus1)

        endif

      endif

    enddo
  enddo

  deallocate(copy_iglob_image_PNM_2D)

  print *,'fin localisation de tous les pixels des images PNM'

!
!---- initialiser sismogrammes
!
  sisux = ZERO
  sisuz = ZERO

  cosrot = cos(anglerec)
  sinrot = sin(anglerec)

! initialiser les tableaux a zero
  accel = ZERO
  veloc = ZERO
  displ = ZERO

!
!----  eventuellement lecture des champs initiaux dans un fichier
!
  if(initialfield) then
    print *
    print *,'Reading initial fields from external file...'
    print *
    open(unit=55,file='OUTPUT_FILES/wavefields.txt',status='unknown')
    read(55,*) nbpoin
    if(nbpoin /= npoin) stop 'Wrong number of points in input file'
    allocate(displread(NDIM))
    allocate(velocread(NDIM))
    allocate(accelread(NDIM))
    do n = 1,npoin
      read(55,*) inump, (displread(i), i=1,NDIM), &
          (velocread(i), i=1,NDIM), (accelread(i), i=1,NDIM)
      if(inump<1 .or. inump>npoin) stop 'Wrong point number'
      displ(:,inump) = displread
      veloc(:,inump) = velocread
      accel(:,inump) = accelread
    enddo
    deallocate(displread)
    deallocate(velocread)
    deallocate(accelread)
    close(55)
    print *,'Max norm of initial displacement = ',maxval(sqrt(displ(1,:)**2 + displ(2,:)**2))
  endif

! attenuation constants from Carcione 1993 Geophysics volume 58 pages 111 and 112
! for two memory-variables mechanisms.
! beware: these values implement specific values of the quality factor Q,
! see Carcione 1993 for details
  tau_epsilon_nu1_mech1 = 0.0334d0
  tau_sigma_nu1_mech1   = 0.0303d0
  tau_epsilon_nu2_mech1 = 0.0352d0
  tau_sigma_nu2_mech1   = 0.0287d0
  tau_epsilon_nu1_mech2 = 0.0028d0
  tau_sigma_nu1_mech2   = 0.0025d0
  tau_epsilon_nu2_mech2 = 0.0029d0
  tau_sigma_nu2_mech2   = 0.0024d0

  inv_tau_sigma_nu1_mech1 = ONE / tau_sigma_nu1_mech1
  inv_tau_sigma_nu2_mech1 = ONE / tau_sigma_nu2_mech1
  inv_tau_sigma_nu1_mech2 = ONE / tau_sigma_nu1_mech2
  inv_tau_sigma_nu2_mech2 = ONE / tau_sigma_nu2_mech2

  phi_nu1_mech1 = (ONE - tau_epsilon_nu1_mech1/tau_sigma_nu1_mech1) / tau_sigma_nu1_mech1
  phi_nu2_mech1 = (ONE - tau_epsilon_nu2_mech1/tau_sigma_nu2_mech1) / tau_sigma_nu2_mech1
  phi_nu1_mech2 = (ONE - tau_epsilon_nu1_mech2/tau_sigma_nu1_mech2) / tau_sigma_nu1_mech2
  phi_nu2_mech2 = (ONE - tau_epsilon_nu2_mech2/tau_sigma_nu2_mech2) / tau_sigma_nu2_mech2

  Mu_nu1 = ONE - (ONE - tau_epsilon_nu1_mech1/tau_sigma_nu1_mech1) - (ONE - tau_epsilon_nu1_mech2/tau_sigma_nu1_mech2)
  Mu_nu2 = ONE - (ONE - tau_epsilon_nu2_mech1/tau_sigma_nu2_mech1) - (ONE - tau_epsilon_nu2_mech2/tau_sigma_nu2_mech2)

  deltatsquare = deltat * deltat
  deltatcube = deltatsquare * deltat
  deltatfourth = deltatsquare * deltatsquare

  twelvedeltat = 12.d0 * deltat
  fourdeltatsquare = 4.d0 * deltatsquare

! --- compute the source time function and store it in a text file
  if(.not. initialfield) then

    allocate(source_time_function(NSTEP))

    print *
    print *,'Saving the source time function in a text file...'
    print *
    open(unit=55,file='OUTPUT_FILES/source.txt',status='unknown')

! boucle principale d'evolution en temps
    do it = 1,NSTEP

! compute current time
      time = (it-1)*deltat

! Ricker (second derivative of a Gaussian) source time function
      if(time_function_type == 1) then
        source_time_function(it) = - factor * (ONE-TWO*a*(time-t0)**2) * exp(-a*(time-t0)**2)

! first derivative of a Gaussian source time function
      else if(time_function_type == 2) then
        source_time_function(it) = - factor * TWO*a*(time-t0) * exp(-a*(time-t0)**2)

! Gaussian or Dirac (we use a very thin Gaussian instead) source time function
      else if(time_function_type == 3 .or. time_function_type == 4) then
        source_time_function(it) = factor * exp(-a*(time-t0)**2)

      else
        stop 'unknown source time function'
      endif

      write(55,*) sngl(time-t0),sngl(source_time_function(it))

    enddo

    close(55)

  else

    allocate(source_time_function(1))

  endif

!
!----          s t a r t   t i m e   i t e r a t i o n s
!

  write(IOUT,400)

! boucle principale d'evolution en temps
  do it = 1,NSTEP

! compute current time
    time = (it-1)*deltat

    if(mod(it,itaff) == 0) then
      write(IOUT,*)
      if(time >= 1.d-3) then
        write(IOUT,100) it,time
      else
        write(IOUT,101) it,time
      endif
    endif

! compute Grad(displ) at time step n for attenuation
  if(TURN_ATTENUATION_ON) call compute_gradient_attenuation(displ,duxdxl_n,duzdxl_n, &
      duxdzl_n,duzdzl_n,xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,NSPEC,npoin)

! update displacement using finite-difference time scheme (Newmark)
    displ(:,:) = displ(:,:) + deltat*veloc(:,:) + deltatsquareover2*accel(:,:)
    veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)
    accel(:,:) = ZERO


!--- free surface for an acoustic medium

! if acoustic, the free surface condition is a Dirichlet condition for the potential,
! not Neumann, in order to impose zero pressure at the surface. Also check that
! top absorbing boundary is not set because cannot be both absorbing and free surface
  if(.not. ELASTIC .and. .not. abshaut) then

    do ispecsurface=1,nelemsurface

      ispec = numsurface(ispecsurface)

      j = NGLLZ
      do i=1,NGLLX
        iglob = ibool(i,j,ispec)
        displ(:,iglob) = ZERO
        veloc(:,iglob) = ZERO
        accel(:,iglob) = ZERO
      enddo

    enddo

  endif  ! end of free surface condition for acoustic medium


!   integration over spectral elements
    do ispec = 1,NSPEC

! get relaxed elastic parameters of current spectral element
      lambdal_relaxed = elastcoef(1,kmato(ispec))
      mul_relaxed = elastcoef(2,kmato(ispec))
      lambdalplus2mul_relaxed = lambdal_relaxed + TWO*mul_relaxed

! first double loop over GLL to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX

!--- if external medium, get elastic parameters of current grid point
          if(read_external_model) then
            iglob = ibool(i,j,ispec)
            cpl = vpext(iglob)
            csl = vsext(iglob)
            rhol = rhoext(iglob)
            mul_relaxed = rhol*csl*csl
            lambdal_relaxed = rhol*cpl*cpl - TWO*mul_relaxed
            lambdalplus2mul_relaxed = lambdal_relaxed + TWO*mul_relaxed
          endif

! compute unrelaxed elastic coefficients from formulas in Carcione 1993 page 111
      lambdal_unrelaxed = (lambdal_relaxed + mul_relaxed) * Mu_nu1 - mul_relaxed * Mu_nu2
      mul_unrelaxed = mul_relaxed * Mu_nu2
      lambdalplus2mul_unrelaxed = lambdal_unrelaxed + TWO*mul_unrelaxed

! derivative along x
          tempx1l = ZERO
          tempz1l = ZERO
          do k = 1,NGLLX
            hp1 = hprime_xx(k,i)
            iglob = ibool(k,j,ispec)
            tempx1l = tempx1l + displ(1,iglob)*hp1
            tempz1l = tempz1l + displ(2,iglob)*hp1
          enddo

! derivative along z
          tempx2l = ZERO
          tempz2l = ZERO
          do k = 1,NGLLZ
            hp2 = hprime_zz(k,j)
            iglob = ibool(i,k,ispec)
            tempx2l = tempx2l + displ(1,iglob)*hp2
            tempz2l = tempz2l + displ(2,iglob)*hp2
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

! derivatives of displacement
          duxdxl = tempx1l*xixl + tempx2l*gammaxl
          duxdzl = tempx1l*xizl + tempx2l*gammazl

          duzdxl = tempz1l*xixl + tempz2l*gammaxl
          duzdzl = tempz1l*xizl + tempz2l*gammazl

! compute stress tensor (include attenuation or anisotropy if needed)

  if(TURN_ATTENUATION_ON) then

! compute the stress using the unrelaxed Lame parameters (Carcione page 111)
    sigma_xx = lambdalplus2mul_unrelaxed*duxdxl + lambdal_unrelaxed*duzdzl
    sigma_xz = mul_unrelaxed*(duzdxl + duxdzl)
    sigma_zz = lambdalplus2mul_unrelaxed*duzdzl + lambdal_unrelaxed*duxdxl

! add the memory variables using the relaxed parameters (Carcione page 111)
! beware: there is a bug in Carcione's equation for sigma_zz
    sigma_xx = sigma_xx + (lambdal_relaxed + mul_relaxed)* &
      (e1_mech1(i,j,k) + e1_mech2(i,j,k)) + TWO * mul_relaxed * (e11_mech1(i,j,k) + e11_mech2(i,j,k))
    sigma_xz = sigma_xz + mul_relaxed * (e13_mech1(i,j,k) + e13_mech2(i,j,k))
    sigma_zz = sigma_zz + (lambdal_relaxed + mul_relaxed)* &
      (e1_mech1(i,j,k) + e1_mech2(i,j,k)) - TWO * mul_relaxed * (e11_mech1(i,j,k) + e11_mech2(i,j,k))

  else

! no attenuation
    sigma_xx = lambdalplus2mul_relaxed*duxdxl + lambdal_relaxed*duzdzl
    sigma_xz = mul_relaxed*(duzdxl + duxdzl)
    sigma_zz = lambdalplus2mul_relaxed*duzdzl + lambdal_relaxed*duxdxl

  endif

! full anisotropy
  if(TURN_ANISOTROPY_ON) then

! implement anisotropy in 2D
     duydyl = ZERO
     duydzl = ZERO
     duzdyl = ZERO
     duxdyl = ZERO
     duydxl = ZERO

! precompute some sums
     duxdxl_plus_duydyl = duxdxl + duydyl
     duxdxl_plus_duzdzl = duxdxl + duzdzl
     duydyl_plus_duzdzl = duydyl + duzdzl
     duxdyl_plus_duydxl = duxdyl + duydxl
     duzdxl_plus_duxdzl = duzdxl + duxdzl
     duzdyl_plus_duydzl = duzdyl + duydzl

     sigma_xx = c11val*duxdxl + c16val*duxdyl_plus_duydxl + c12val*duydyl + &
        c15val*duzdxl_plus_duxdzl + c14val*duzdyl_plus_duydzl + c13val*duzdzl

!     sigma_yy = c12val*duxdxl + c26val*duxdyl_plus_duydxl + c22val*duydyl + &
!        c25val*duzdxl_plus_duxdzl + c24val*duzdyl_plus_duydzl + c23val*duzdzl

     sigma_zz = c13val*duxdxl + c36val*duxdyl_plus_duydxl + c23val*duydyl + &
        c35val*duzdxl_plus_duxdzl + c34val*duzdyl_plus_duydzl + c33val*duzdzl

!     sigma_xy = c16val*duxdxl + c66val*duxdyl_plus_duydxl + c26val*duydyl + &
!        c56val*duzdxl_plus_duxdzl + c46val*duzdyl_plus_duydzl + c36val*duzdzl

     sigma_xz = c15val*duxdxl + c56val*duxdyl_plus_duydxl + c25val*duydyl + &
        c55val*duzdxl_plus_duxdzl + c45val*duzdyl_plus_duydzl + c35val*duzdzl

!     sigma_yz = c14val*duxdxl + c46val*duxdyl_plus_duydxl + c24val*duydyl + &
!        c45val*duzdxl_plus_duxdzl + c44val*duzdyl_plus_duydzl + c34val*duzdzl

  endif

! stress tensor is symmetric
          sigma_zx = sigma_xz

          jacobianl = jacobian(i,j,ispec)

! weak formulation term based on stress tensor (non-symmetric form)
          tempx1(i,j) = jacobianl*(sigma_xx*xixl+sigma_zx*xizl)
          tempz1(i,j) = jacobianl*(sigma_xz*xixl+sigma_zz*xizl)

          tempx2(i,j) = jacobianl*(sigma_xx*gammaxl+sigma_zx*gammazl)
          tempz2(i,j) = jacobianl*(sigma_xz*gammaxl+sigma_zz*gammazl)

! for acoustic medium
          if(.not. ELASTIC) then
            tempx1(i,j) = jacobianl*(xixl*dUxdxl + xizl*dUxdzl)
            tempx2(i,j) = jacobianl*(gammaxl*dUxdxl + gammazl*dUxdzl)
          endif

        enddo
      enddo

!
! second double-loop over GLL to compute all terms
!
      do j = 1,NGLLZ
        do i = 1,NGLLX

! along x direction
          tempx1l = ZERO
          tempz1l = ZERO
          do k = 1,NGLLX
            fac1 = wxgll(k)*hprime_xx(i,k)
            tempx1l = tempx1l + tempx1(k,j)*fac1
            if(ELASTIC) tempz1l = tempz1l + tempz1(k,j)*fac1
          enddo

! along z direction
          tempx2l = ZERO
          tempz2l = ZERO
          do k = 1,NGLLZ
            fac2 = wzgll(k)*hprime_zz(j,k)
            tempx2l = tempx2l + tempx2(i,k)*fac2
            if(ELASTIC) tempz2l = tempz2l + tempz2(i,k)*fac2
          enddo

! GLL integration weights
          fac1 = wzgll(j)
          fac2 = wxgll(i)

! for acoustic medium
          iglob = ibool(i,j,ispec)
          accel(1,iglob) = accel(1,iglob) - (fac1*tempx1l + fac2*tempx2l)
          if(ELASTIC) then
            accel(2,iglob) = accel(2,iglob) - (fac1*tempz1l + fac2*tempz2l)
          else
            accel(2,iglob) = zero
          endif

        enddo ! second loop over the GLL points
      enddo

    enddo ! end of loop over all spectral elements

!
!--- absorbing boundaries
!
  if(anyabs) then

    do ispecabs=1,nelemabs

      ispec = numabs(ispecabs)

! get elastic parameters of current spectral element
      lambdal_relaxed = elastcoef(1,kmato(ispec))
      mul_relaxed = elastcoef(2,kmato(ispec))
      rhol  = density(kmato(ispec))
      kappal  = lambdal_relaxed + TWO*mul_relaxed/3.d0
      cpl = sqrt((kappal + 4.d0*mul_relaxed/3.d0)/rhol)
      csl = sqrt(mul_relaxed/rhol)


!--- left absorbing boundary
      if(codeabs(ILEFT,ispecabs)) then

        i = 1

        do j=1,NGLLZ

          iglob = ibool(i,j,ispec)

          zgamma = xix(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(read_external_model) then
            cpl = vpext(iglob)
            csl = vsext(iglob)
            rhol = rhoext(iglob)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          nx = -ONE
          nz = ZERO

          vx = veloc(1,iglob)
          vz = veloc(2,iglob)

          vn = nx*vx+nz*vz

          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

          weight = zgamma*wzgll(j)

! Clayton-Engquist condition if elastic, Sommerfeld condition if acoustic
          if(ELASTIC) then
            accel(1,iglob) = accel(1,iglob) - tx*weight
            accel(2,iglob) = accel(2,iglob) - tz*weight
          else
            accel(1,iglob) = accel(1,iglob) - veloc(1,iglob)*weight/cpl
          endif

        enddo

      endif  !  end of left absorbing boundary

!--- right absorbing boundary
      if(codeabs(IRIGHT,ispecabs)) then

        i = NGLLX

        do j=1,NGLLZ

          iglob = ibool(i,j,ispec)

          zgamma = xix(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(read_external_model) then
            cpl = vpext(iglob)
            csl = vsext(iglob)
            rhol = rhoext(iglob)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          nx = ONE
          nz = ZERO

          vx = veloc(1,iglob)
          vz = veloc(2,iglob)

          vn = nx*vx+nz*vz

          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

          weight = zgamma*wzgll(j)

! Clayton-Engquist condition if elastic, Sommerfeld condition if acoustic
          if(ELASTIC) then
            accel(1,iglob) = accel(1,iglob) - tx*weight
            accel(2,iglob) = accel(2,iglob) - tz*weight
          else
            accel(1,iglob) = accel(1,iglob) - veloc(1,iglob)*weight/cpl
          endif

        enddo

      endif  !  end of right absorbing boundary

!--- bottom absorbing boundary
      if(codeabs(IBOTTOM,ispecabs)) then

        j = 1

! exclude corners to make sure there is no contradiction on the normal
        i1abs = 1
        i2abs = NGLLX
        if(codeabs(ILEFT,ispecabs)) i1abs = 2
        if(codeabs(IRIGHT,ispecabs)) i2abs = NGLLX-1

        do i=i1abs,i2abs

          iglob = ibool(i,j,ispec)

          xxi = gammaz(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(read_external_model) then
            cpl = vpext(iglob)
            csl = vsext(iglob)
            rhol = rhoext(iglob)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          nx = ZERO
          nz = -ONE

          vx = veloc(1,iglob)
          vz = veloc(2,iglob)

          vn = nx*vx+nz*vz

          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

          weight = xxi*wxgll(i)

! Clayton-Engquist condition if elastic, Sommerfeld condition if acoustic
          if(ELASTIC) then
            accel(1,iglob) = accel(1,iglob) - tx*weight
            accel(2,iglob) = accel(2,iglob) - tz*weight
          else
            accel(1,iglob) = accel(1,iglob) - veloc(1,iglob)*weight/cpl
          endif

        enddo

      endif  !  end of bottom absorbing boundary

!--- top absorbing boundary
      if(codeabs(ITOP,ispecabs)) then

        j = NGLLZ

! exclude corners to make sure there is no contradiction on the normal
        i1abs = 1
        i2abs = NGLLX
        if(codeabs(ILEFT,ispecabs)) i1abs = 2
        if(codeabs(IRIGHT,ispecabs)) i2abs = NGLLX-1

        do i=i1abs,i2abs

          iglob = ibool(i,j,ispec)

          xxi = gammaz(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(read_external_model) then
            cpl = vpext(iglob)
            csl = vsext(iglob)
            rhol = rhoext(iglob)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          nx = ZERO
          nz = ONE

          vx = veloc(1,iglob)
          vz = veloc(2,iglob)

          vn = nx*vx+nz*vz

          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

          weight = xxi*wxgll(i)

! Clayton-Engquist condition if elastic, Sommerfeld condition if acoustic
          if(ELASTIC) then
            accel(1,iglob) = accel(1,iglob) - tx*weight
            accel(2,iglob) = accel(2,iglob) - tz*weight
          else
            accel(1,iglob) = accel(1,iglob) - veloc(1,iglob)*weight/cpl
          endif

        enddo

      endif  !  end of top absorbing boundary

    enddo

  endif  ! end of absorbing boundaries


! --- add the source
  if(.not. initialfield) then

! collocated force
! beware, for acoustic medium, source is a potential, therefore source time function
! gives shape of velocity, not displacement
  if(source_type == 1) then
    if(ELASTIC) then
      accel(1,iglob_source) = accel(1,iglob_source) - sin(angleforce)*source_time_function(it)
      accel(2,iglob_source) = accel(2,iglob_source) + cos(angleforce)*source_time_function(it)
    else
      accel(1,iglob_source) = accel(1,iglob_source) + source_time_function(it)
    endif

! moment tensor
  else if(source_type == 2) then

! add source array
    do j=1,NGLLZ
      do i=1,NGLLX
        iglob = ibool(i,j,ispec_selected_source)
        accel(:,iglob) = accel(:,iglob) + sourcearray(:,i,j)*source_time_function(it)
      enddo
    enddo

  endif

  else
    stop 'wrong source type'
  endif

! divide by the mass matrix
  accel(1,:) = accel(1,:) / rmass(:)
  accel(2,:) = accel(2,:) / rmass(:)

! update velocity
  veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)


!--- free surface for an acoustic medium

! if acoustic, the free surface condition is a Dirichlet condition for the potential,
! not Neumann, in order to impose zero pressure at the surface. Also check that
! top absorbing boundary is not set because cannot be both absorbing and free surface
  if(.not. ELASTIC .and. .not. abshaut) then

    do ispecsurface=1,nelemsurface

      ispec = numsurface(ispecsurface)

      j = NGLLZ
      do i=1,NGLLX
        iglob = ibool(i,j,ispec)
        displ(:,iglob) = ZERO
        veloc(:,iglob) = ZERO
        accel(:,iglob) = ZERO
      enddo

    enddo

  endif  ! end of free surface condition for acoustic medium


! implement attenuation
  if(TURN_ATTENUATION_ON) then

! compute Grad(displ) at time step n+1 for attenuation
    call compute_gradient_attenuation(displ,duxdxl_np1,duzdxl_np1, &
      duxdzl_np1,duzdzl_np1,xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,NSPEC,npoin)

! update memory variables with fourth-order Runge-Kutta time scheme for attenuation
  do k=1,nspec
  do j=1,NGLLZ
  do i=1,NGLLX

  theta_n   = duxdxl_n(i,j,k) + duzdzl_n(i,j,k)
  theta_np1 = duxdxl_np1(i,j,k) + duzdzl_np1(i,j,k)

! evolution e1_mech1
  Un = e1_mech1(i,j,k)
  tauinv = - inv_tau_sigma_nu1_mech1
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = theta_n * phi_nu1_mech1
  Snp1 = theta_np1 * phi_nu1_mech1
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e1_mech1(i,j,k) = Unp1

! evolution e1_mech2
  Un = e1_mech2(i,j,k)
  tauinv = - inv_tau_sigma_nu1_mech2
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = theta_n * phi_nu1_mech2
  Snp1 = theta_np1 * phi_nu1_mech2
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e1_mech2(i,j,k) = Unp1

! evolution e11_mech1
  Un = e11_mech1(i,j,k)
  tauinv = - inv_tau_sigma_nu2_mech1
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = (duxdxl_n(i,j,k) - theta_n/TWO) * phi_nu2_mech1
  Snp1 = (duxdxl_np1(i,j,k) - theta_np1/TWO) * phi_nu2_mech1
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e11_mech1(i,j,k) = Unp1

! evolution e11_mech2
  Un = e11_mech2(i,j,k)
  tauinv = - inv_tau_sigma_nu2_mech2
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = (duxdxl_n(i,j,k) - theta_n/TWO) * phi_nu2_mech2
  Snp1 = (duxdxl_np1(i,j,k) - theta_np1/TWO) * phi_nu2_mech2
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e11_mech2(i,j,k) = Unp1

! evolution e13_mech1
  Un = e13_mech1(i,j,k)
  tauinv = - inv_tau_sigma_nu2_mech1
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = (duxdzl_n(i,j,k) + duzdxl_n(i,j,k)) * phi_nu2_mech1
  Snp1 = (duxdzl_np1(i,j,k) + duzdxl_np1(i,j,k)) * phi_nu2_mech1
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e13_mech1(i,j,k) = Unp1

! evolution e13_mech2
  Un = e13_mech2(i,j,k)
  tauinv = - inv_tau_sigma_nu2_mech2
  tauinvsquare = tauinv * tauinv
  tauinvcube = tauinvsquare * tauinv
  tauinvUn = tauinv * Un
  Sn   = (duxdzl_n(i,j,k) + duzdxl_n(i,j,k)) * phi_nu2_mech2
  Snp1 = (duxdzl_np1(i,j,k) + duzdxl_np1(i,j,k)) * phi_nu2_mech2
  Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
      twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
      fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
      deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
  e13_mech2(i,j,k) = Unp1

  enddo
  enddo
  enddo

  endif ! end of test on attenuation

!----  display max of norm of displacement
  if(mod(it,itaff) == 0) then
    displnorm_all = maxval(sqrt(displ(1,:)**2 + displ(2,:)**2))
    print *,'Max norm of field = ',displnorm_all
! check stability of the code, exit if unstable
    if(displnorm_all > STABILITY_THRESHOLD) stop 'code became unstable and blew up'
  endif

! store the seismograms
  if(sismostype < 1 .or. sismostype > 3) stop 'Wrong field code for seismogram output'

  if(.not. ELASTIC) then
    if(sismostype == 1) then
      stop 'cannot store displacement field in acoustic medium because of potential formulation'
    else if(sismostype == 2) then
! for acoustic medium, compute gradient for display, displ represents the potential
      call compute_gradient_fluid(displ,vector_field_postscript, &
            xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,NSPEC,npoin)
    else
! for acoustic medium, compute gradient for display, veloc represents the first derivative of the potential
      call compute_gradient_fluid(veloc,vector_field_postscript, &
            xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,NSPEC,npoin)
    endif
  endif

  do irec=1,nrec

! perform the general interpolation using Lagrange polynomials
    valux = ZERO
    valuz = ZERO

    do j = 1,NGLLZ
      do i = 1,NGLLX

        iglob = ibool(i,j,ispec_selected_rec(irec))

        hlagrange = hxir_store(irec,i)*hgammar_store(irec,j)

        if(ELASTIC) then

          if(sismostype == 1) then
            dxd = displ(1,iglob)
            dzd = displ(2,iglob)
          else if(sismostype == 2) then
            dxd = veloc(1,iglob)
            dzd = veloc(2,iglob)
          else
            dxd = accel(1,iglob)
            dzd = accel(2,iglob)
          endif

        else

! for acoustic medium
          dxd = vector_field_postscript(1,iglob)
          dzd = vector_field_postscript(2,iglob)

        endif

! compute interpolated field
        valux = valux + dxd*hlagrange
        valuz = valuz + dzd*hlagrange

      enddo
    enddo

! rotate seismogram components if needed
    sisux(it,irec) =   cosrot*valux + sinrot*valuz
    sisuz(it,irec) = - sinrot*valux + cosrot*valuz

  enddo

!
!----  affichage des resultats a certains pas de temps
!
  if(mod(it,itaff) == 0 .or. it == 5 .or. it == NSTEP) then

  write(IOUT,*)
  if(time >= 1.d-3) then
    write(IOUT,110) time
  else
    write(IOUT,111) time
  endif
  write(IOUT,*)

!
!----  affichage postscript
!
  if(output_postscript_image) then

  write(IOUT,*) 'Dump PostScript'

! for elastic medium
  if(ELASTIC .and. vecttype == 1) then
    write(IOUT,*) 'drawing displacement field...'
    call plotpost(displ,coord,vpext,x_source,z_source,st_xval,st_zval, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,stitle,npoin,npgeo,vpmin,vpmax,nrec, &
          colors,numbers,subsamp,vecttype,interpol,meshvect,modelvect, &
          boundvect,read_external_model,cutvect,nelemabs,numat,iptsdisp,nspec,ngnod,ELASTIC)

  else if(ELASTIC .and. vecttype == 2) then
    write(IOUT,*) 'drawing velocity field...'
    call plotpost(veloc,coord,vpext,x_source,z_source,st_xval,st_zval, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,stitle,npoin,npgeo,vpmin,vpmax,nrec, &
          colors,numbers,subsamp,vecttype,interpol,meshvect,modelvect, &
          boundvect,read_external_model,cutvect,nelemabs,numat,iptsdisp,nspec,ngnod,ELASTIC)

  else if(ELASTIC .and. vecttype == 3) then
    write(IOUT,*) 'drawing acceleration field...'
    call plotpost(accel,coord,vpext,x_source,z_source,st_xval,st_zval, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,stitle,npoin,npgeo,vpmin,vpmax,nrec, &
          colors,numbers,subsamp,vecttype,interpol,meshvect,modelvect, &
          boundvect,read_external_model,cutvect,nelemabs,numat,iptsdisp,nspec,ngnod,ELASTIC)

! for acoustic medium
  else if(.not. ELASTIC .and. vecttype == 1) then
    stop 'cannot display displacement field in acoustic medium because of potential formulation'

  else if(.not. ELASTIC .and. vecttype == 2) then
    write(IOUT,*) 'drawing acoustic velocity field from velocity potential...'
! for acoustic medium, compute gradient for display, displ represents the potential
    call compute_gradient_fluid(displ,vector_field_postscript, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,NSPEC,npoin)
    call plotpost(vector_field_postscript,coord,vpext,x_source,z_source,st_xval,st_zval, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,stitle,npoin,npgeo,vpmin,vpmax,nrec, &
          colors,numbers,subsamp,vecttype,interpol,meshvect,modelvect, &
          boundvect,read_external_model,cutvect,nelemabs,numat,iptsdisp,nspec,ngnod,ELASTIC)

  else if(.not. ELASTIC .and. vecttype == 3) then
    write(IOUT,*) 'drawing acoustic acceleration field from velocity potential...'
! for acoustic medium, compute gradient for display, veloc represents the first derivative of the potential
    call compute_gradient_fluid(veloc,vector_field_postscript, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,NSPEC,npoin)
    call plotpost(vector_field_postscript,coord,vpext,x_source,z_source,st_xval,st_zval, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,stitle,npoin,npgeo,vpmin,vpmax,nrec, &
          colors,numbers,subsamp,vecttype,interpol,meshvect,modelvect, &
          boundvect,read_external_model,cutvect,nelemabs,numat,iptsdisp,nspec,ngnod,ELASTIC)

  else
    stop 'wrong field code for PostScript display'
  endif
  write(IOUT,*) 'Fin dump PostScript'

  endif

!
!----  affichage image PNM
!
  if(output_PNM_image) then

  write(IOUT,*) 'Creation image PNM de taille ',NX_IMAGE_PNM,' x ',NZ_IMAGE_PNM

  donnees_image_PNM_2D(:,:) = 0.d0

  do j = 1,NZ_IMAGE_PNM
    do i = 1,NX_IMAGE_PNM
      if(iglob_image_PNM_2D(i,j) /= -1) then
! display vertical component of vector
        if(ELASTIC) then
          if(vecttype == 1) then
            donnees_image_PNM_2D(i,j) = displ(2,iglob_image_PNM_2D(i,j))
          else if(vecttype == 2) then
            donnees_image_PNM_2D(i,j) = veloc(2,iglob_image_PNM_2D(i,j))
          else
            donnees_image_PNM_2D(i,j) = accel(2,iglob_image_PNM_2D(i,j))
          endif
        else
! for acoustic medium
          donnees_image_PNM_2D(i,j) = vector_field_postscript(2,iglob_image_PNM_2D(i,j))
        endif
      endif
    enddo
  enddo

  call cree_image_PNM(donnees_image_PNM_2D,iglob_image_PNM_2D,NX_IMAGE_PNM,NZ_IMAGE_PNM,it,cutvect)

  write(IOUT,*) 'Fin creation image PNM'

  endif

!----  save temporary seismograms
  call write_seismograms(sisux,sisuz,station_name,network_name,NSTEP,nrec,deltat,sismostype,st_xval,it,t0)


  endif

  enddo ! end of the main time loop

!----  save final seismograms
  call write_seismograms(sisux,sisuz,station_name,network_name,NSTEP,nrec,deltat,sismostype,st_xval,it,t0)

! print exit banner
  call datim(stitle)

!
!----  close output file
!
  close(IOUT)

!
!----  formats
!
 40   format(a80)
 45   format(a50)
 100  format('Pas de temps numero ',i5,'   t = ',f7.4,' s')
 101  format('Pas de temps numero ',i5,'   t = ',1pe10.4,' s')
 110  format('Sauvegarde deplacement temps t = ',f7.4,' s')
 111  format('Sauvegarde deplacement temps t = ',1pe10.4,' s')
 400  format(/1x,41('=')/,' =  T i m e  e v o l u t i o n  l o o p  ='/1x,41('=')/)

  200   format(//1x,'C o n t r o l',/1x,13('='),//5x,&
  'Number of spectral element control nodes. . .(npgeo) =',i8/5x, &
  'Number of space dimensions. . . . . . . . . . (NDIM) =',i8)
  600   format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Display frequency  . . . . . . . . . . . . . (itaff) = ',i5/ 5x, &
  'Color display . . . . . . . . . . . . . . . (colors) = ',i5/ 5x, &
  '        ==  0     black and white display              ',  / 5x, &
  '        ==  1     color display                        ',  /5x, &
  'Numbered mesh . . . . . . . . . . . . . . .(numbers) = ',i5/ 5x, &
  '        ==  0     do not number the mesh               ',  /5x, &
  '        ==  1     number the mesh                      ')
  700   format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Seismograms recording type . . . . . . .(sismostype) = ',i6/5x, &
  'Angle for first line of receivers. . . . .(anglerec) = ',f6.2)
  750   format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Read external initial field. . . . . .(initialfield) = ',l6/5x, &
  'Read external velocity model. .(read_external_model) = ',l6/5x, &
  'Elastic simulation or acoustic. . . . . . .(ELASTIC) = ',l6/5x, &
  'Turn anisotropy on or off. . . .(TURN_ANISOTROPY_ON) = ',l6/5x, &
  'Turn attenuation on or off. . .(TURN_ATTENUATION_ON) = ',l6/5x, &
  'Save grid in external file or not. . . .(outputgrid) = ',l6)
  800   format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Vector display type. . . . . . . . . . . .(vecttype) = ',i6/5x, &
  'Percentage of cut for vector plots. . . . .(cutvect) = ',f6.2/5x, &
  'Subsampling for velocity model display. . .(subsamp) = ',i6)

  703   format(//' I t e r a t i o n s '/1x,19('='),//5x, &
      'Number of time iterations . . . . .(NSTEP) =',i8,/5x, &
      'Time step increment. . . . . . . .(deltat) =',1pe15.6,/5x, &
      'Total simulation duration . . . . . (ttot) =',1pe15.6)

  107   format(/5x,'--> Isoparametric Spectral Elements <--',//)
  207   format(5x, &
           'Number of spectral elements . . . . .  (nspec) =',i7,/5x, &
           'Number of control nodes per element .  (ngnod) =',i7,/5x, &
           'Number of points in X-direction . . .  (NGLLX) =',i7,/5x, &
           'Number of points in Y-direction . . .  (NGLLZ) =',i7,/5x, &
           'Number of points per element. . .(NGLLX*NGLLZ) =',i7,/5x, &
           'Number of points for display . . . .(iptsdisp) =',i7,/5x, &
           'Number of element material sets . . .  (numat) =',i7,/5x, &
           'Number of absorbing elements . . . .(nelemabs) =',i7)

  212   format(//,5x, &
     'Source Type. . . . . . . . . . . . . . = Collocated Force',/5x, &
     'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
     'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
     'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
     'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
     'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
     'Angle from vertical direction (deg). . =',1pe20.10,/5x)
  222   format(//,5x, &
     'Source Type. . . . . . . . . . . . . . = Moment-tensor',/5x, &
     'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
     'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
     'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
     'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
     'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
     'Mxx. . . . . . . . . . . . . . . . . . =',1pe20.10,/5x, &
     'Mzz. . . . . . . . . . . . . . . . . . =',1pe20.10,/5x, &
     'Mxz. . . . . . . . . . . . . . . . . . =',1pe20.10)

  end program specfem2D

