
!========================================================================
!
!                   M E S H F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                         (c) January 2005
!
!========================================================================

!========================================================================
!
!  Basic mesh generator for SPECFEM2D
!
!========================================================================

  program meshfem2D

  implicit none

! coordinates of the grid points of the mesh
  double precision, dimension(:,:), allocatable :: x,z

! to compute the coordinate transformation
  integer :: ioffset
  double precision :: gamma,absx,a00,a01,bot0,top0

! stockage du modele de vitesse et densite
  double precision, dimension(:), allocatable :: rho,cp,cs,aniso3,aniso4
  integer, dimension(:), allocatable :: icodemat
  integer, dimension(:,:), allocatable :: num_modele

! interface data
  integer interface_current,ipoint_current,number_of_interfaces,npoints_interface_bottom,npoints_interface_top
  integer ilayer,number_of_layers,max_npoints_interface
  double precision xinterface_dummy,zinterface_dummy,xinterface_dummy_previous
  integer, dimension(:), allocatable :: nz_layer
  double precision, dimension(:), allocatable :: &
         xinterface_bottom,zinterface_bottom,coefs_interface_bottom, &
         xinterface_top,zinterface_top,coefs_interface_top

! for the source
  integer source_type,time_function_type
  double precision xs,zs,f0,t0,angleforce,Mxx,Mzz,Mxz,factor

! arrays for the receivers
  double precision, dimension(:), allocatable :: xrec,zrec

  character(len=50) interfacesfile,title

  integer imatnum,inumabs,inumsurface,inumelem
  integer nelemabs,nelemsurface,npgeo,nspec
  integer k,icol,ili,istepx,istepz,ix,iz,irec,i,j
  integer ixdebzone,ixfinzone,izdebzone,izfinzone,imodnum
  integer izone,imodele,nbzone,nbmodeles
  integer itaff,pointsdisp,subsamp,nrec
  integer sismostype,vecttype
  integer ngnod,nt,nx,nz,nxread,nzread
  integer icodematread

  logical codehaut,codebas,codegauche,codedroite

  double precision tang1,tangN,vpzone,vszone,poisson_ratio
  double precision cutvect,xspacerec,zspacerec
  double precision anglerec,xfin,zfin,xdeb,zdeb,xmin,xmax,dt
  double precision rhoread,cpread,csread,aniso3read,aniso4read

  logical interpol,gnuplot,readmodel,outputgrid
  logical abshaut,absbas,absgauche,absdroite
  logical source_surf,enreg_surf,meshvect,initialfield,modelvect,boundvect
  logical ELASTIC,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON

  integer, external :: num
  double precision, external :: value_spline

! flag to indicate an anisotropic material
  integer, parameter :: ANISOTROPIC_MATERIAL = 1

! file number for input DATA/Par_file and interface file
  integer, parameter :: IIN_PAR = 10
  integer, parameter :: IIN_INTERFACES = 15

! ignore variable name field (junk) at the beginning of each input line
  logical, parameter :: IGNORE_JUNK = .true.,DONT_IGNORE_JUNK = .false.

! ***
! *** read the parameter file
! ***

  print *,'Reading the parameter file ... '
  print *

  open(unit=10,file='DATA/Par_file',status='old')

! read file names and path for output
  call read_value_string(IIN_PAR,IGNORE_JUNK,title)
  call read_value_string(IIN_PAR,IGNORE_JUNK,interfacesfile)

  write(*,*) 'Titre de la simulation'
  write(*,*) title
  print *

! read grid parameters
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,xmin)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,xmax)
  call read_value_integer(IIN_PAR,IGNORE_JUNK,nx)
  call read_value_integer(IIN_PAR,IGNORE_JUNK,ngnod)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,initialfield)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,readmodel)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,ELASTIC)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,TURN_ANISOTROPY_ON)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,TURN_ATTENUATION_ON)

! get interface data from external file to count the spectral elements along Z
  print *,'Reading interface data from file DATA/',interfacesfile(1:len_trim(interfacesfile)),' to count the spectral elements'
  open(unit=15,file='DATA/'//interfacesfile,status='old')

  max_npoints_interface = -1

! read number of interfaces
  call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,number_of_interfaces)
  if(number_of_interfaces < 2) stop 'not enough interfaces (minimum is 2)'

! loop on all the interfaces
  do interface_current = 1,number_of_interfaces

    call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface_bottom)
    if(npoints_interface_bottom < 2) stop 'not enough interface points (minimum is 2)'
    max_npoints_interface = max(npoints_interface_bottom,max_npoints_interface)
    print *,'Reading ',npoints_interface_bottom,' points for interface ',interface_current

! loop on all the points describing this interface
    do ipoint_current = 1,npoints_interface_bottom
      read(IIN_INTERFACES,*) xinterface_dummy,zinterface_dummy
      if(ipoint_current > 1 .and. xinterface_dummy <= xinterface_dummy_previous) &
        stop 'interface points must be sorted in increasing X'
      xinterface_dummy_previous = xinterface_dummy
    enddo

  enddo

! define number of layers
  number_of_layers = number_of_interfaces - 1

  allocate(nz_layer(number_of_layers))

! loop on all the layers
  do ilayer = 1,number_of_layers

! read number of spectral elements in vertical direction in this layer
    call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,nz_layer(ilayer))
    if(nz_layer(ilayer) < 1) stop 'not enough spectral elements along Z in layer (minimum is 1)'
    print *,'There are ',nz_layer(ilayer),' spectral elements along Z in layer ',ilayer

  enddo

  close(15)

! compute total number of spectral elements in vertical direction
  nz = sum(nz_layer)

  print *
  print *,'Total number of spectral elements along Z = ',nz
  print *

  nxread = nx
  nzread = nz

! multiplier par 2 si elements 9 noeuds
  if(ngnod == 9) then
    nx = nx * 2
    nz = nz * 2
  endif

! read absorbing boundaries parameters
  call read_value_logical(IIN_PAR,IGNORE_JUNK,abshaut)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,absbas)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,absgauche)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,absdroite)

! read time step parameters
  call read_value_integer(IIN_PAR,IGNORE_JUNK,nt)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,dt)

! read source parameters
  call read_value_logical(IIN_PAR,IGNORE_JUNK,source_surf)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,xs)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,zs)
  call read_value_integer(IIN_PAR,IGNORE_JUNK,source_type)
  call read_value_integer(IIN_PAR,IGNORE_JUNK,time_function_type)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,f0)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,angleforce)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,Mxx)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,Mzz)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,Mxz)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,factor)

! if Dirac source time function, use a very thin Gaussian instead
  if(time_function_type == 4) f0 = 1.d0 / (5.d0 * dt)

! time delay of the source in seconds, use a 20 % security margin
  t0 = 1.20d0 / f0

  print *
  print *,'Source:'
  print *,'Position xs, zs = ',xs,zs
  print *,'Frequency, delay = ',f0,t0
  print *,'Source type (1=force, 2=explosion): ',source_type
  print *,'Time function type (1=Ricker, 2=First derivative, 3=Gaussian, 4=Dirac): ',time_function_type
  print *,'Angle of the source if force = ',angleforce
  print *,'Mxx of the source if moment tensor = ',Mxx
  print *,'Mzz of the source if moment tensor = ',Mzz
  print *,'Mxz of the source if moment tensor = ',Mxz
  print *,'Multiplying factor = ',factor

! read receivers line parameters
  call read_value_logical(IIN_PAR,IGNORE_JUNK,enreg_surf)
  call read_value_integer(IIN_PAR,IGNORE_JUNK,sismostype)
  call read_value_integer(IIN_PAR,IGNORE_JUNK,nrec)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,xdeb)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,zdeb)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,xfin)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,zfin)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,anglerec)

  allocate(xrec(nrec))
  allocate(zrec(nrec))

  print *
  print *,'There are ',nrec,' receivers'
  xspacerec = (xfin-xdeb) / dble(nrec-1)
  zspacerec = (zfin-zdeb) / dble(nrec-1)
  do i=1,nrec
    xrec(i) = xdeb + dble(i-1)*xspacerec
    zrec(i) = zdeb + dble(i-1)*zspacerec
  enddo

! read display parameters
  call read_value_integer(IIN_PAR,IGNORE_JUNK,itaff)
  call read_value_integer(IIN_PAR,IGNORE_JUNK,vecttype)
  call read_value_double_precision(IIN_PAR,IGNORE_JUNK,cutvect)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,meshvect)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,modelvect)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,boundvect)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,interpol)
  call read_value_integer(IIN_PAR,IGNORE_JUNK,pointsdisp)
  call read_value_integer(IIN_PAR,IGNORE_JUNK,subsamp)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,gnuplot)
  call read_value_logical(IIN_PAR,IGNORE_JUNK,outputgrid)

! lecture des differents modeles de materiaux
  call read_value_integer(IIN_PAR,IGNORE_JUNK,nbmodeles)
  if(nbmodeles <= 0) stop 'Negative number of models not allowed !!'

  allocate(icodemat(nbmodeles))
  allocate(rho(nbmodeles))
  allocate(cp(nbmodeles))
  allocate(cs(nbmodeles))
  allocate(aniso3(nbmodeles))
  allocate(aniso4(nbmodeles))
  allocate(num_modele(nx,nz))

  icodemat(:) = 0
  rho(:) = 0.d0
  cp(:) = 0.d0
  cs(:) = 0.d0
  aniso3(:) = 0.d0
  aniso4(:) = 0.d0
  num_modele(:,:) = 0

  do imodele=1,nbmodeles
    read(IIN_PAR,*) i,icodematread,rhoread,cpread,csread,aniso3read,aniso4read
    if(i < 1 .or. i > nbmodeles) stop 'Wrong model number!!'
    icodemat(i) = icodematread
    rho(i) = rhoread
    cp(i) = cpread
    cs(i) = csread

    if(rho(i) <= 0.d0 .or. cp(i) <= 0.d0 .or. cs(i) < 0.d0) stop 'negative value of velocity or density'

! check that Cs = 0 if acoustic simulation
    if(.not. ELASTIC .and. cs(i) > 0.0001) stop 'must have Cs = 0 for acoustic model'

    aniso3(i) = aniso3read
    aniso4(i) = aniso4read
  enddo

  print *
  print *, 'Nb de modeles de roche = ',nbmodeles
  print *
  do i=1,nbmodeles
    if(icodemat(i) /= ANISOTROPIC_MATERIAL) then
      print *,'Modele #',i,' isotrope'
      print *,'rho,cp,cs = ',rho(i),cp(i),cs(i)
    else
      print *,'Modele #',i,' anisotrope'
      print *,'rho,c11,c13,c33,c44 = ',rho(i),cp(i),cs(i),aniso3(i),aniso4(i)
    endif
  enddo
  print *

! lecture des numeros de modele des differentes zones
  call read_value_integer(IIN_PAR,IGNORE_JUNK,nbzone)

  if(nbzone <= 0) stop 'Negative number of zones not allowed !!'

  print *
  print *, 'Nb de zones du modele = ',nbzone
  print *

  do izone = 1,nbzone

    read(IIN_PAR,*) ixdebzone,ixfinzone,izdebzone,izfinzone,imodnum

    if(imodnum < 1) stop 'Negative model number not allowed !!'
    if(ixdebzone < 1) stop 'Left coordinate of zone negative !!'
    if(ixfinzone > nxread) stop 'Right coordinate of zone too high !!'
    if(izdebzone < 1) stop 'Bottom coordinate of zone negative !!'
    if(izfinzone > nzread) stop 'Top coordinate of zone too high !!'

    print *,'Zone ',izone
    print *,'IX from ',ixdebzone,' to ',ixfinzone
    print *,'IZ from ',izdebzone,' to ',izfinzone

  if(icodemat(imodnum) /= ANISOTROPIC_MATERIAL) then
    vpzone = cp(imodnum)
    vszone = cs(imodnum)
    print *,'Model # ',imodnum,' isotrope'
    print *,'vp = ',vpzone
    print *,'vs = ',vszone
    print *,'rho = ',rho(imodnum)
    poisson_ratio = 0.5d0*(vpzone*vpzone-2.d0*vszone*vszone) / (vpzone*vpzone-vszone*vszone)
    print *,'Poisson''s ratio = ',poisson_ratio
    if(poisson_ratio <= -1.00001d0 .or. poisson_ratio >= 0.50001d0) stop 'incorrect value of Poisson''s ratio'
  else
    print *,'Model # ',imodnum,' anisotrope'
    print *,'c11 = ',cp(imodnum)
    print *,'c13 = ',cs(imodnum)
    print *,'c33 = ',aniso3(imodnum)
    print *,'c44 = ',aniso4(imodnum)
    print *,'rho = ',rho(imodnum)
  endif
  print *,' -----'

! stocker le modele de vitesse et densite
   do i = ixdebzone,ixfinzone
     do j = izdebzone,izfinzone
       if(ngnod == 4) then
         num_modele(i,j) = imodnum
       else
         num_modele(2*(i-1)+1,2*(j-1)+1) = imodnum
         num_modele(2*(i-1)+1,2*(j-1)+2) = imodnum
         num_modele(2*(i-1)+2,2*(j-1)+1) = imodnum
         num_modele(2*(i-1)+2,2*(j-1)+2) = imodnum
       endif
     enddo
   enddo

  enddo

  if(minval(num_modele) <= 0) stop 'Velocity model not entirely set...'

  close(10)

  print *
  print *,' Parameter file successfully read... '

! --------- fin lecture fichier parametres --------------

  if(ngnod /= 4 .and. ngnod /= 9) stop 'erreur ngnod different de 4 ou 9 !!'

  print *
  if(ngnod == 4) then
    print *,'Le maillage comporte ',nx,' x ',nz,' elements'
  else
    print *,'Le maillage comporte ',nx/2,' x ',nz/2,' elements'
  endif
  print *
  print *,'Les elements de controle sont des elements ',ngnod,' noeuds'
  print *

!---

! perform basic checks on parameters read

! for acoustic
  if(TURN_ANISOTROPY_ON .and. .not. ELASTIC) stop 'currently cannot have anisotropy in acoustic simulation'

  if(TURN_ATTENUATION_ON .and. .not. ELASTIC) stop 'currently cannot have attenuation in acoustic simulation'

  if(source_type == 2 .and. .not. ELASTIC) stop 'currently cannot have moment tensor source in acoustic simulation'

! for attenuation
  if(TURN_ANISOTROPY_ON .and. TURN_ATTENUATION_ON) stop 'cannot have anisotropy and attenuation both turned on in current version'

!---

! allocate arrays for the grid
  allocate(x(0:nx,0:nz))
  allocate(z(0:nx,0:nz))

  x(:,:) = 0.d0
  z(:,:) = 0.d0

! get interface data from external file
  print *,'Reading interface data from file DATA/',interfacesfile(1:len_trim(interfacesfile))
  open(unit=15,file='DATA/'//interfacesfile,status='old')

  allocate(xinterface_bottom(max_npoints_interface))
  allocate(zinterface_bottom(max_npoints_interface))
  allocate(coefs_interface_bottom(max_npoints_interface))

  allocate(xinterface_top(max_npoints_interface))
  allocate(zinterface_top(max_npoints_interface))
  allocate(coefs_interface_top(max_npoints_interface))

! read number of interfaces
  call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,number_of_interfaces)

! read bottom interface
  call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface_bottom)

! loop on all the points describing this interface
  do ipoint_current = 1,npoints_interface_bottom
    read(IIN_INTERFACES,*) xinterface_bottom(ipoint_current),zinterface_bottom(ipoint_current)
  enddo

! boucle sur toutes les couches
  do ilayer = 1,number_of_layers

! read top interface
  call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface_top)

! loop on all the points describing this interface
  do ipoint_current = 1,npoints_interface_top
    read(IIN_INTERFACES,*) xinterface_top(ipoint_current),zinterface_top(ipoint_current)
  enddo

! calculer le spline pour l'interface du bas, imposer la tangente aux deux bords
  tang1 = (zinterface_bottom(2)-zinterface_bottom(1)) / (xinterface_bottom(2)-xinterface_bottom(1))
  tangN = (zinterface_bottom(npoints_interface_bottom)-zinterface_bottom(npoints_interface_bottom-1)) / &
          (xinterface_bottom(npoints_interface_bottom)-xinterface_bottom(npoints_interface_bottom-1))
  call spline(xinterface_bottom,zinterface_bottom,npoints_interface_bottom,tang1,tangN,coefs_interface_bottom)

! calculer le spline pour l'interface du haut, imposer la tangente aux deux bords
  tang1 = (zinterface_top(2)-zinterface_top(1)) / (xinterface_top(2)-xinterface_top(1))
  tangN = (zinterface_top(npoints_interface_top)-zinterface_top(npoints_interface_top-1)) / &
          (xinterface_top(npoints_interface_top)-xinterface_top(npoints_interface_top-1))
  call spline(xinterface_top,zinterface_top,npoints_interface_top,tang1,tangN,coefs_interface_top)

! tester si on est sur la derniere couche, qui contient la topographie
  if(ilayer == number_of_layers) then

! modifier position de la source si source exactement en surface
    if(source_surf) zs = value_spline(xs,xinterface_top,zinterface_top,coefs_interface_top,npoints_interface_top)

! modifier position des recepteurs si enregistrement exactement en surface
    if(enreg_surf) then
      do irec=1,nrec
        zrec(irec) = value_spline(xrec(irec),xinterface_top,zinterface_top,coefs_interface_top,npoints_interface_top)
      enddo
    endif

  endif

! calcul de l'offset de cette couche en nombre d'elements spectraux suivant Z
  if(ilayer > 1) then
    ioffset = sum(nz_layer(1:ilayer-1))
  else
    ioffset = 0
  endif

!--- definition du maillage

    do ix = 0,nx

! points regulierement espaces suivant X
      absx = xmin + (xmax - xmin) * dble(ix) / dble(nx)

! value of the bottom and top splines
      bot0 = value_spline(absx,xinterface_bottom,zinterface_bottom,coefs_interface_bottom,npoints_interface_bottom)
      top0 = value_spline(absx,xinterface_top,zinterface_top,coefs_interface_top,npoints_interface_top)

      do iz = 0,nz_layer(ilayer)

! linear interpolation between bottom and top
        gamma = dble(iz) / dble(nz_layer(ilayer))
        a00 = 1.d0 - gamma
        a01 = gamma

! coordinates of the grid points
        x(ix,iz + ioffset) = absx
        z(ix,iz + ioffset) = a00*bot0 + a01*top0

      enddo

    enddo

! l'interface du haut devient celle du bas pour passer a la couche suivante
    npoints_interface_bottom = npoints_interface_top
    xinterface_bottom(:) = xinterface_top(:)
    zinterface_bottom(:) = zinterface_top(:)

  enddo

  close(15)

! calculer min et max de X et Z sur la grille
  print *
  print *,'Valeurs min et max de X sur le maillage = ',minval(x),maxval(x)
  print *,'Valeurs min et max de Z sur le maillage = ',minval(z),maxval(z)
  print *

! afficher position de la source et des recepteurs
  print *
  print *,'Position (x,z) de la source'
  print *
  print *,'Source = ',xs,zs
  print *
  print *,'Position (x,z) des ',nrec,' recepteurs'
  print *
  do irec=1,nrec
    print *,'Receiver ',irec,' = ',xrec(irec),zrec(irec)
  enddo

! ***
! *** generer un fichier Gnuplot pour le controle de la grille ***
! ***

  print *
  print *,'Ecriture de la grille format Gnuplot...'

  open(unit=20,file='OUTPUT_FILES/gridfile.gnu',status='unknown')

! dessin des lignes horizontales de la grille
  print *,'Ecriture lignes horizontales'
  istepx = 1
  if(ngnod == 4) then
    istepz = 1
  else
    istepz = 2
  endif
  do ili=0,nz,istepz
    do icol=0,nx-istepx,istepx
      write(20,15) sngl(x(icol,ili)),sngl(z(icol,ili))
      write(20,15) sngl(x(icol+istepx,ili)),sngl(z(icol+istepx,ili))
      write(20,10)
    enddo
  enddo

! dessin des lignes verticales de la grille
  print *,'Ecriture lignes verticales'
  if(ngnod == 4) then
    istepx = 1
  else
    istepx = 2
  endif
  istepz = 1
  do icol=0,nx,istepx
    do ili=0,nz-istepz,istepz
      write(20,15) sngl(x(icol,ili)),sngl(z(icol,ili))
      write(20,15) sngl(x(icol,ili+istepz)),sngl(z(icol,ili+istepz))
      write(20,10)
    enddo
  enddo

  close(20)

! cree le script de dessin pour gnuplot
  open(unit=20,file='OUTPUT_FILES/plotgnu',status='unknown')
  write(20,*) '#set term postscript landscape monochrome solid "Helvetica" 22'
  write(20,*) '#set output "grille.ps"'
  write(20,*) 'plot "gridfile.gnu" title "Macrobloc mesh" w l'
  write(20,*) 'pause -1 "appuyez sur une touche"'
  close(20)

  print *,'Fin ecriture de la grille format Gnuplot'
  print *

! *** generation de la base de donnees

  open(unit=15,file='OUTPUT_FILES/Database',status='unknown')

  write(15,*) '#'
  write(15,*) '# Database for SPECFEM2D'
  write(15,*) '# Dimitri Komatitsch, (c) University of Pau, France'
  write(15,*) '#'

  write(15,*) 'Titre simulation'
  write(15,40) title

  npgeo = (nx+1)*(nz+1)
  if(ngnod == 4) then
    nspec = nx*nz
  else
    nspec = nx*nz/4
  endif
  write(15,*) 'npgeo'
  write(15,*) npgeo

  write(15,*) 'gnuplot interpol'
  write(15,*) gnuplot,interpol

  write(15,*) 'itaff colors numbers'
  write(15,*) itaff,' 1 0'

  write(15,*) 'meshvect modelvect boundvect cutvect subsamp nx_sem_PNM'
  write(15,*) meshvect,modelvect,boundvect,cutvect,subsamp,nxread

  write(15,*) 'anglerec'
  write(15,*) anglerec

  write(15,*) 'initialfield'
  write(15,*) initialfield

  write(15,*) 'sismostype vecttype'
  write(15,*) sismostype,vecttype

  write(15,*) 'readmodel outputgrid ELASTIC TURN_ANISOTROPY_ON TURN_ATTENUATION_ON'
  write(15,*) readmodel,outputgrid,ELASTIC,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON

  write(15,*) 'ncycl dtinc'
  write(15,*) nt,dt

  write(15,*) 'source'
  write(15,*) source_type,time_function_type,xs-xmin,zs,f0,t0,factor,angleforce,Mxx,Mzz,Mxz

  write(15,*) 'Coordinates of macrobloc mesh (coorg):'
  do j=0,nz
    do i=0,nx
      write(15,*) num(i,j,nx),x(i,j)-xmin,z(i,j)
    enddo
  enddo

!
!--- introduction des bords absorbants
!
  nelemabs = 0
  if(absbas) nelemabs = nelemabs + nx
  if(abshaut) nelemabs = nelemabs + nx
  if(absgauche) nelemabs = nelemabs + nz
  if(absdroite) nelemabs = nelemabs + nz

! on a deux fois trop d'elements si elements 9 noeuds
  if(ngnod == 9) nelemabs = nelemabs / 2

! enlever aussi les coins qui ont ete comptes deux fois
  if(absbas .and. absgauche) nelemabs = nelemabs - 1
  if(absbas .and. absdroite) nelemabs = nelemabs - 1
  if(abshaut .and. absgauche) nelemabs = nelemabs - 1
  if(abshaut .and. absdroite) nelemabs = nelemabs - 1

!
!--- introduction de la surface libre si milieu acoustique
!
  nelemsurface = nx

! on a deux fois trop d'elements si elements 9 noeuds
  if(ngnod == 9) nelemsurface = nelemsurface / 2

  write(15,*) 'numat ngnod nspec pointsdisp nelemabs nelemsurface'
  write(15,*) nbmodeles,ngnod,nspec,pointsdisp,nelemabs,nelemsurface

  write(15,*) 'Material sets (num 0 rho vp vs 0 0) or (num 1 rho c11 c13 c33 c44)'
  do i=1,nbmodeles
    write(15,*) i,icodemat(i),rho(i),cp(i),cs(i),aniso3(i),aniso4(i)
  enddo


  write(15,*) 'Arrays kmato and knods for each bloc:'

  k=0
  if(ngnod == 4) then
    do j=0,nz-1
    do i=0,nx-1

    k = k + 1
    imatnum = num_modele(i+1,j+1)
    write(15,*) k,imatnum,num(i,j,nx),num(i+1,j,nx),num(i+1,j+1,nx),num(i,j+1,nx)
    enddo
    enddo
  else
    do j=0,nz-2,2
    do i=0,nx-2,2

    k = k + 1
    imatnum = num_modele(i+1,j+1)
    write(15,*) k,imatnum,num(i,j,nx),num(i+2,j,nx),num(i+2,j+2,nx), &
              num(i,j+2,nx),num(i+1,j,nx),num(i+2,j+1,nx), &
              num(i+1,j+2,nx),num(i,j+1,nx),num(i+1,j+1,nx)

    enddo
    enddo
  endif

!
!--- sauvegarde des bords absorbants
!
  print *
  print *,'Au total il y a ',nelemabs,' elements absorbants'
  print *
  print *,'Bords absorbants actifs :'
  print *
  print *,'Haut   = ',abshaut
  print *,'Bas    = ',absbas
  print *,'Gauche = ',absgauche
  print *,'Droite = ',absdroite
  print *

! generer la liste des elements absorbants
  if(nelemabs > 0) then
  write(15,*) 'Liste des elements absorbants (haut bas gauche droite) :'
  inumabs = 0
  do iz = 1,nzread
  do ix = 1,nxread
    codehaut = .false.
    codebas = .false.
    codegauche = .false.
    codedroite = .false.
    inumelem = (iz-1)*nxread + ix
    if(abshaut   .and. iz == nzread) codehaut = .true.
    if(absbas    .and. iz == 1) codebas = .true.
    if(absgauche .and. ix == 1) codegauche = .true.
    if(absdroite .and. ix == nxread) codedroite = .true.
    if(codehaut .or. codebas .or. codegauche .or. codedroite) then
      inumabs = inumabs + 1
      write(15,*) inumabs,inumelem,codehaut,codebas,codegauche,codedroite
    endif
  enddo
  enddo
  endif

!
!--- sauvegarde de la surface libre
!
  print *
  print *,'Au total il y a ',nelemsurface,' elements a la surface libre'
  print *

! generer la liste des elements a la surface libre
  if(nelemsurface > 0) then
  write(15,*) 'Liste des elements a la surface libre'
  write(15,*) abshaut
  inumsurface = 0
  do iz = 1,nzread
  do ix = 1,nxread
    inumelem = (iz-1)*nxread + ix
    if(iz == nzread) then
      inumsurface = inumsurface + 1
      write(15,*) inumsurface,inumelem
    endif
  enddo
  enddo
  endif

  close(15)


!
!--- write the STATIONS file
!
  print *
  print *,'writing the DATA/STATIONS file'
  print *

  open(unit=15,file='DATA/STATIONS',status='unknown')
  write(15,*) nrec
  do irec=1,nrec
    write(15,100) irec,xrec(irec)-xmin,zrec(irec)
  enddo
  close(15)

 10 format('')
 15 format(e10.5,1x,e10.5)
 40 format(a50)

 100 format('S',i3.3,'    AA ',f20.7,1x,f20.7,'       0.0         0.0')

  end program meshfem2D

! ********************
! routines de maillage
! ********************

!--- numero global du noeud

  integer function num(i,j,nx)

  implicit none

  integer i,j,nx

    num = j*(nx+1) + i + 1

  end function num

!--- representation des interfaces par un spline

  double precision function value_spline(x,xinterface,zinterface,coefs_interface,npoints_interface)

  implicit none

  integer npoints_interface
  double precision x,xp
  double precision, dimension(npoints_interface) :: xinterface,zinterface,coefs_interface

  value_spline = 0.d0

  xp = x

! si on sort du modele, prolonger par continuite
  if(xp < xinterface(1)) xp = xinterface(1)
  if(xp > xinterface(npoints_interface)) xp = xinterface(npoints_interface)

  call splint(xinterface,zinterface,coefs_interface,npoints_interface,xp,value_spline)

  end function value_spline

! --------------------------------------

! routine de calcul des coefs du spline (adapted from Numerical Recipes)

  subroutine spline(x,y,n,yp1,ypn,y2)

  implicit none

  integer n
  double precision, dimension(n) :: x,y,y2
  double precision, dimension(:), allocatable :: u
  double precision yp1,ypn

  integer i,k
  double precision sig,p,qn,un

  allocate(u(n))

  y2(1)=-0.5d0
  u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)

  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.d0
    y2(i)=(sig-1.d0)/p
    u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
               /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  enddo

  qn=0.5d0
  un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo

  deallocate(u)

  end subroutine spline

! --------------

! routine d'evaluation du spline (adapted from Numerical Recipes)

  subroutine splint(XA,YA,Y2A,N,X,Y)

  implicit none

  integer n
  double precision, dimension(n) :: XA,YA,Y2A
  double precision x,y

  integer k,klo,khi
  double precision h,a,b

  KLO = 1
  KHI = N

  do while (KHI-KLO > 1)
    K=(KHI+KLO)/2
    if(XA(K) > X) then
      KHI=K
    else
      KLO=K
    endif
  enddo

  H = XA(KHI)-XA(KLO)
  IF (H == 0.d0) stop 'bad input in spline evaluation'

  A = (XA(KHI)-X) / H
  B = (X-XA(KLO)) / H

  Y = A*YA(KLO) + B*YA(KHI) + ((A**3-A)*Y2A(KLO) + (B**3-B)*Y2A(KHI))*(H**2)/6.d0

  end subroutine splint

!--------------------

  subroutine read_value_integer(iin,ignore_junk,value_to_read)

  implicit none

  integer iin
  logical ignore_junk

  integer value_to_read

  integer ios

  character(len=100) string_read
  character(len=34) junk

  do
    read(unit=iin,fmt=200,iostat=ios) string_read
    if(ios /= 0) stop 'error while reading input file'

! suppress leading white spaces, if any
    string_read = adjustl(string_read)

! exit loop when we find the first line that is not a comment
    if(string_read(1:1) /= '#') exit

  enddo

  if(ignore_junk) then
    read(string_read,100) junk,value_to_read
  else
    read(string_read,*) value_to_read
  endif

! format
 100 format(a,i8)
 200 format(a100)

  end subroutine read_value_integer

!--------------------

  subroutine read_value_double_precision(iin,ignore_junk,value_to_read)

  implicit none

  integer iin
  logical ignore_junk

  double precision value_to_read

  integer ios

  character(len=100) string_read
  character(len=34) junk

  do
    read(unit=iin,fmt=200,iostat=ios) string_read
    if(ios /= 0) stop 'error while reading input file'

! suppress leading white spaces, if any
    string_read = adjustl(string_read)

! exit loop when we find the first line that is not a comment
    if(string_read(1:1) /= '#') exit

  enddo

  if(ignore_junk) then
    read(string_read,100) junk,value_to_read
  else
    read(string_read,*) value_to_read
  endif

! format
 100 format(a,f12.5)
 200 format(a100)

  end subroutine read_value_double_precision

!--------------------

  subroutine read_value_logical(iin,ignore_junk,value_to_read)

  implicit none

  integer iin
  logical ignore_junk

  logical value_to_read

  integer ios

  character(len=100) string_read
  character(len=34) junk

  do
    read(unit=iin,fmt=200,iostat=ios) string_read
    if(ios /= 0) stop 'error while reading input file'

! suppress leading white spaces, if any
    string_read = adjustl(string_read)

! exit loop when we find the first line that is not a comment
    if(string_read(1:1) /= '#') exit

  enddo

  if(ignore_junk) then
    read(string_read,100) junk,value_to_read
  else
    read(string_read,*) value_to_read
  endif

! format
 100 format(a,l8)
 200 format(a100)

  end subroutine read_value_logical

!--------------------

  subroutine read_value_string(iin,ignore_junk,value_to_read)

  implicit none

  integer iin
  logical ignore_junk

  character(len=*) value_to_read

  integer ios

  character(len=100) string_read
  character(len=34) junk

  do
    read(unit=iin,fmt=200,iostat=ios) string_read
    if(ios /= 0) stop 'error while reading input file'

! suppress leading white spaces, if any
    string_read = adjustl(string_read)

! exit loop when we find the first line that is not a comment
    if(string_read(1:1) /= '#') exit

  enddo

  if(ignore_junk) then
    read(string_read,100) junk,value_to_read
  else
    read(string_read,*) value_to_read
  endif

! format
 100 format(a34,a)
 200 format(a100)

  end subroutine read_value_string

