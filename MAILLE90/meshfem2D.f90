
!========================================================================
!
!                   M E S H F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                         (c) December 2004
!
!========================================================================

!========================================================================
!
!  Basic mesh generator for SPECFEM2D
!
!========================================================================

  program meshfem2D

  implicit none

! definir les tableaux pour allocation dynamique

  integer, parameter :: ANISOTROPIC_MATERIAL = 1

! coordinates of the grid points of the mesh
  double precision, dimension(:,:), allocatable :: x,z

! to compute the coordinate transformation
  integer :: ioffset
  double precision :: eta,absx,a00,a01,bot0,top0

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
  integer isource_type
  double precision xs,zs,f0,t0,angle,factor

! arrays for the receivers
  double precision, dimension(:), allocatable :: xrec,zrec

  character(len=50) interfacesfile,title
  character(len=34) junk

  integer imatnum,inumabs,inumelem
  integer nelemabs,npgeo,nspec
  integer k,icol,ili,istepx,istepz,ix,iz,irec,i,j
  integer ixdebzone,ixfinzone,izdebzone,izfinzone,imodnum
  integer izone,imodele,nbzone,nbmodeles
  integer itaff,iptsdisp,isubsamp,nrec
  integer isismostype,ivecttype
  integer ngnod,nt,nx,nz,nxread,nzread
  integer icodematread

  logical codehaut,codebas,codegauche,codedroite

  double precision tang1,tangN,vpzone,vszone
  double precision cutvect,xspacerec,zspacerec
  double precision anglerec,xfin,zfin,xdeb,zdeb,xmin,xmax,dt
  double precision rhoread,cpread,csread,aniso3read,aniso4read

  logical interpol,ignuplot,ireadmodel,ioutputgrid
  logical abshaut,absbas,absgauche,absdroite
  logical isource_surf,ienreg_surf,imeshvect,initialfield,imodelvect,iboundvect
  logical TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON

  integer, external :: num
  double precision, external :: value_spline

! ***
! *** read the parameter file
! ***

  print *,'Reading the parameter file ... '
  print *

  open(unit=10,file='Par',status='old')

! formats
 1 format(a,f12.5)
 2 format(a,i8)
 3 format(a,a)
 4 format(a,l8)

! read the header
  do i=1,10
    read(10,*)
  enddo

! read file names and path for output
  read(10,3) junk,title
  read(10,3) junk,interfacesfile

  write(*,*) 'Titre de la simulation'
  write(*,*) title
  print *

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read grid parameters
  read(10,1) junk,xmin
  read(10,1) junk,xmax
  read(10,2) junk,nx
  read(10,2) junk,ngnod
  read(10,4) junk,initialfield
  read(10,4) junk,ireadmodel
  read(10,4) junk,TURN_ANISOTROPY_ON
  read(10,4) junk,TURN_ATTENUATION_ON

! get interface data from external file to count the spectral elements along Z
  print *,'Reading interface data from file ',interfacesfile(1:len_trim(interfacesfile)),' to count the spectral elements'
  open(unit=15,file=interfacesfile,status='old')

  max_npoints_interface = -1

! skip header
  read(15,*)
  read(15,*)
  read(15,*)

! read number of interfaces
  read(15,*) number_of_interfaces
  if(number_of_interfaces < 2) stop 'not enough interfaces (minimum is 2)'

! skip header
  read(15,*)
  read(15,*)
  read(15,*)

! loop on all the interfaces
  do interface_current = 1,number_of_interfaces

! skip header
    read(15,*)
    read(15,*)
    read(15,*)

    read(15,*) npoints_interface_bottom
    if(npoints_interface_bottom < 2) stop 'not enough interface points (minimum is 2)'
    max_npoints_interface = max(npoints_interface_bottom,max_npoints_interface)
    print *,'Reading ',npoints_interface_bottom,' points for interface ',interface_current

! loop on all the points describing this interface
    do ipoint_current = 1,npoints_interface_bottom
      read(15,*) xinterface_dummy,zinterface_dummy
      if(ipoint_current > 1 .and. xinterface_dummy <= xinterface_dummy_previous) &
        stop 'interface points must be sorted in increasing X'
      xinterface_dummy_previous = xinterface_dummy
    enddo

  enddo

! define number of layers
  number_of_layers = number_of_interfaces - 1

  allocate(nz_layer(number_of_layers))

! skip header
    read(15,*)
    read(15,*)
    read(15,*)

! loop on all the layers
  do ilayer = 1,number_of_layers

! skip header
    read(15,*)
    read(15,*)
    read(15,*)

! read number of spectral elements in vertical direction in this layer
    read(15,*) nz_layer(ilayer)
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

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read absorbing boundaries parameters
  read(10,4) junk,abshaut
  read(10,4) junk,absbas
  read(10,4) junk,absgauche
  read(10,4) junk,absdroite

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read time step parameters
  read(10,2) junk,nt
  read(10,1) junk,dt

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read source parameters
  read(10,4) junk,isource_surf
  read(10,1) junk,xs
  read(10,1) junk,zs
  read(10,1) junk,f0
  read(10,1) junk,t0
  read(10,2) junk,isource_type
  read(10,1) junk,angle
  read(10,1) junk,factor

  print *
  print *,'Source:'
  print *,'Position xs, zs = ',xs,zs
  print *,'Frequency, delay = ',f0,t0
  print *,'Source type (1=force 2=explosion) : ',isource_type
  print *,'Angle of the source if force = ',angle
  print *,'Multiplying factor = ',factor

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read receivers line parameters
  read(10,4) junk,ienreg_surf
  read(10,2) junk,isismostype
  read(10,2) junk,nrec
  read(10,1) junk,xdeb
  read(10,1) junk,zdeb
  read(10,1) junk,xfin
  read(10,1) junk,zfin
  read(10,1) junk,anglerec

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

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read display parameters
  read(10,2) junk,itaff
  read(10,2) junk,ivecttype
  read(10,1) junk,cutvect
  read(10,4) junk,imeshvect
  read(10,4) junk,imodelvect
  read(10,4) junk,iboundvect
  read(10,4) junk,interpol
  read(10,2) junk,iptsdisp
  read(10,2) junk,isubsamp
  read(10,4) junk,ignuplot
  read(10,4) junk,ioutputgrid

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! lecture des differents modeles de materiaux
  read(10,2) junk,nbmodeles
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
    read(10,*) i,icodematread,rhoread,cpread,csread,aniso3read,aniso4read
    if(i < 1 .or. i > nbmodeles) stop 'Wrong material point number'
    icodemat(i) = icodematread
    rho(i) = rhoread
    cp(i) = cpread
    cs(i) = csread
    aniso3(i) = aniso3read
    aniso4(i) = aniso4read
    if(i <= 0) stop 'Negative model number not allowed !!'
    if(rho(i) <= 0.d0 .or. cp(i) <= 0.d0 .or. cs(i) < 0.d0) stop 'negative value of velocity or density'
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
  read(10,2) junk,nbzone

  if(nbzone <= 0) stop 'Negative number of zones not allowed !!'

  print *
  print *, 'Nb de zones du modele = ',nbzone
  print *

  do izone = 1,nbzone

    read(10,*) ixdebzone,ixfinzone,izdebzone,izfinzone,imodnum

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
    print *,'Poisson''s ratio = ', &
      0.5d0*(vpzone*vpzone-2.d0*vszone*vszone) / (vpzone*vpzone-vszone*vszone)
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

! allocate arrays for the grid
  allocate(x(0:nx,0:nz))
  allocate(z(0:nx,0:nz))

  x(:,:) = 0.d0
  z(:,:) = 0.d0

! get interface data from external file
  print *,'Reading interface data from file ',interfacesfile(1:len_trim(interfacesfile))
  open(unit=15,file=interfacesfile,status='old')

  allocate(xinterface_bottom(max_npoints_interface))
  allocate(zinterface_bottom(max_npoints_interface))
  allocate(coefs_interface_bottom(max_npoints_interface))

  allocate(xinterface_top(max_npoints_interface))
  allocate(zinterface_top(max_npoints_interface))
  allocate(coefs_interface_top(max_npoints_interface))

! skip header
  read(15,*)
  read(15,*)
  read(15,*)

! read number of interfaces
  read(15,*) number_of_interfaces

! skip header
  read(15,*)
  read(15,*)
  read(15,*)

! skip header
  read(15,*)
  read(15,*)
  read(15,*)

! read bottom interface
  read(15,*) npoints_interface_bottom

! loop on all the points describing this interface
  do ipoint_current = 1,npoints_interface_bottom
    read(15,*) xinterface_bottom(ipoint_current),zinterface_bottom(ipoint_current)
  enddo

! boucle sur toutes les couches
  do ilayer = 1,number_of_layers

! skip header
  read(15,*)
  read(15,*)
  read(15,*)

! read top interface
  read(15,*) npoints_interface_top

! loop on all the points describing this interface
  do ipoint_current = 1,npoints_interface_top
    read(15,*) xinterface_top(ipoint_current),zinterface_top(ipoint_current)
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
    if(isource_surf) zs = value_spline(xs,xinterface_top,zinterface_top,coefs_interface_top,npoints_interface_top)

! modifier position des recepteurs si enregistrement exactement en surface
    if(ienreg_surf) then
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
        eta = dble(iz) / dble(nz_layer(ilayer))
        a00 = 1.d0 - eta
        a01 = eta

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

  open(unit=20,file='gridfile.gnu',status='unknown')

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
  open(unit=20,file='plotgnu',status='unknown')
  write(20,*) '#set term postscript landscape monochrome solid "Helvetica" 22'
  write(20,*) '#set output "grille.ps"'
  write(20,*) 'plot "gridfile.gnu" title "Macrobloc mesh" w l'
  write(20,*) 'pause -1 "appuyez sur une touche"'
  close(20)

  print *,'Fin ecriture de la grille format Gnuplot'
  print *

! *** generation de la base de donnees

  open(unit=15,file='../SPECFEM90/DataBase',status='unknown')

  write(15,*) '#'
  write(15,*) '# DataBase for SPECFEM2D'
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

  write(15,*) 'ignuplot interpol'
  write(15,*) ignuplot,interpol

  write(15,*) 'itaff icolor inumber'
  write(15,*) itaff,1,0

  write(15,*) 'imeshvect imodelvect iboundvect cutvect isubsamp'
  write(15,*) imeshvect,imodelvect,iboundvect,cutvect,isubsamp

  write(15,*) 'nrec anglerec'
  write(15,*) nrec,anglerec

  write(15,*) 'initialfield'
  write(15,*) initialfield

  write(15,*) 'isismostype ivecttype'
  write(15,*) isismostype,ivecttype

  write(15,*) 'ireadmodel ioutputgrid TURN_ANISOTROPY_ON TURN_ATTENUATION_ON'
  write(15,*) ireadmodel,ioutputgrid,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON

  write(15,*) 'ncycl dtinc'
  write(15,*) nt,dt

  write(15,*) 'source'
  write(15,*) isource_type,xs-xmin,zs,f0,t0,factor,angle,' 0'

  write(15,*) 'Receiver positions:'
  do irec=1,nrec
    write(15,*) irec,xrec(irec)-xmin,zrec(irec)
  enddo

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

  write(15,*) 'numat ngnod nspec iptsdisp ielemabs'
  write(15,*) nbmodeles,ngnod,nspec,iptsdisp,nelemabs

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

  close(15)

 10 format('')
 15 format(e10.5,1x,e10.5)
 40 format(a50)

  end program meshfem2D

! ********************
! routines de maillage
! ********************

! --- numero global du noeud

  integer function num(i,j,nx)

  implicit none

  integer i,j,nx

    num = j*(nx+1) + i + 1

  end function num

! --- representation des interfaces par un spline

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

