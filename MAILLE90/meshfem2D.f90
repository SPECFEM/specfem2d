
!========================================================================
!
!                   M E S H F E M 2 D  Version 5.0
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) May 2004
!
!========================================================================

!========================================================================
!
!  Mesh generator for SPECFEM2D version 5.0
!
!========================================================================

  program meshfem2D

  implicit none

! definir les tableaux pour allocation dynamique

! coordinates of the grid points
  double precision, allocatable :: x(:,:),z(:,:)

! variables needed to compute the transformation
  double precision, allocatable :: psi(:),eta(:),absx(:), &
      a00(:),a01(:),valeta(:),bot0(:),top0(:)

! stockage du modele de vitesse et densite
  double precision, allocatable :: rho(:),cp(:),cs(:),aniso3(:),aniso4(:)
  integer, allocatable :: icodemat(:)
  integer, allocatable :: num_modele(:,:)

! the topography data
  double precision, allocatable :: xtopo(:),ztopo(:),coefs_topo(:), &
      xinterf(:),zinterf(:),coefs_interf(:)

! for the source
  double precision xs,zs,f0,t0,angle,factor
  integer isource_type

! arrays for the receivers
  double precision, allocatable :: xrec(:),zrec(:)

! nom du fichier GNUPLOT contenant la grille
  character(len=50) file1

  character(len=50) interffile,topofile,title
  character(len=15) junk

  integer imatnum,inumabs,inumelem,netyp
  integer nelemabs,npgeo,nspec,ninterf,ntopo
  integer k,icol,ili,istepx,istepz,ncut,ix,iz,irec,i,j
  integer ixdebzone,ixfinzone,izdebzone,izfinzone,imodnum
  integer izone,imodele,nbzone,nbmodeles
  integer itaff,iptsdisp,isubsamp,nrec
  integer isismostype,ivecttype
  integer ngnod,nt,nx,nz,nxread,nzread
  integer icodematread

  logical codehaut,codebas,codegauche,codedroite

  double precision ratio
  double precision tang1,tangN,vpzone,vszone
  double precision cutvect
  double precision xspacerec,zspacerec
  double precision anglerec,xfin,zfin,xdeb,zdeb,xmin,xmax
  double precision dt
  double precision rhoread,cpread,csread,aniso3read,aniso4read

  logical interpol,ignuplot,ireadmodel,ioutputgrid
  logical abshaut,absbas,absgauche,absdroite
  logical isource_surf,ienreg_surf
  logical imeshvect
  logical initialfield
  logical imodelvect,iboundvect

  integer, external :: num
  double precision, external :: bottom,spl,dens

! ***
! *** read the parameter file
! ***

  print *,' Reading the parameter file ... '
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
  read(10,3)junk,title
  read(10,3)junk,topofile
  read(10,3)junk,interffile

  write(*,*) 'Titre de la simulation'
  write(*,*) title
  print *

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read grid parameters
  read(10,1)junk,xmin
  read(10,1)junk,xmax
  read(10,2)junk,nx
  read(10,2)junk,nz
  read(10,2)junk,ngnod
  read(10,1)junk,ratio
  read(10,4)junk,initialfield
  read(10,4)junk,ireadmodel

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
  read(10,4)junk,abshaut
  read(10,4)junk,absbas
  read(10,4)junk,absgauche
  read(10,4)junk,absdroite

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read time step parameters
  read(10,2)junk,nt
  read(10,1)junk,dt

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read source parameters
  read(10,4)junk,isource_surf
  read(10,1)junk,xs
  read(10,1)junk,zs
  read(10,1)junk,f0
  read(10,1)junk,t0
  read(10,2)junk,isource_type
  read(10,1)junk,angle
  read(10,1)junk,factor

  print *
  print *,'Source:'
  print *,'Position xs, zs = ',xs,zs
  print *,'Frequency, delay = ',f0,t0
  print *,'Source type (1=force 2=explo) : ',isource_type
  print *,'Angle of the source if force = ',angle
  print *,'Multiplying factor = ',factor

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read receivers line parameters
  read(10,4)junk,ienreg_surf
  read(10,2)junk,isismostype
  read(10,2)junk,nrec
  read(10,1)junk,xdeb
  read(10,1)junk,zdeb
  read(10,1)junk,xfin
  read(10,1)junk,zfin
  read(10,1)junk,anglerec

  allocate(xrec(nrec))
  allocate(zrec(nrec))

  print *
  print *,'There are ',nrec,' receivers'
  xspacerec=(xfin-xdeb)/dble(nrec-1)
  zspacerec=(zfin-zdeb)/dble(nrec-1)
  do i=1,nrec
    xrec(i) = xdeb + dble(i-1)*xspacerec
    zrec(i) = zdeb + dble(i-1)*zspacerec
  enddo

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read display parameters
  read(10,2)junk,itaff
  read(10,2)junk,ivecttype
  read(10,1)junk,cutvect
  read(10,4)junk,imeshvect
  read(10,4)junk,imodelvect
  read(10,4)junk,iboundvect
  read(10,4)junk,interpol
  read(10,2)junk,iptsdisp
  read(10,2)junk,isubsamp
  read(10,4)junk,ignuplot
  read(10,4)junk,ioutputgrid

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! lecture des differents modeles de materiaux

  read(10,2)junk,nbmodeles
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
      if(i<1 .or. i>nbmodeles) stop 'Wrong material point number'
      icodemat(i) = icodematread
      rho(i) = rhoread
      cp(i) = cpread
      cs(i) = csread
      aniso3(i) = aniso3read
      aniso4(i) = aniso4read
      if(i <= 0) stop 'Negative model number not allowed !!'
      if (rho(i) < 0.d0 .or. cp(i) < 0.d0 .or. cs(i) < 0.d0) &
          stop 'Negative value of velocity or density'
  enddo

  print *
  print *, 'Nb de modeles de roche = ',nbmodeles
  print *
  do i=1,nbmodeles
    if(icodemat(i) /= 2) then
      print *,'Modele #',i,' isotrope'
      print *,'rho,cp,cs = ',rho(i),cp(i),cs(i)
    else
      print *,'Modele #',i,' anisotrope'
      print *,'rho,c11,c13,c33,c44 = ',rho(i),cp(i),cs(i),aniso3(i),aniso4(i)
    endif
  enddo
  print *

! *** lecture des numeros de modele des differentes zones
  read(10,2)junk,nbzone

  if(nbzone <= 0) stop 'Negative number of zones not allowed !!'

  print *
  print *, 'Nb de zones du modele = ',nbzone
  print *

  do izone=1,nbzone
      read(10,*) ixdebzone,ixfinzone,izdebzone,izfinzone,imodnum

  if (imodnum < 1) stop 'Negative model number not allowed !!'
  if (ixdebzone < 1) stop 'Left coordinate of zone negative !!'
  if (ixfinzone > nxread) stop 'Right coordinate of zone too high !!'
  if (izdebzone < 1) stop 'Bottom coordinate of zone negative !!'
  if (izfinzone > nzread) stop 'Top coordinate of zone too high !!'

      print *,'Zone ',izone
      print *,'IX from ',ixdebzone,' to ',ixfinzone
      print *,'IZ from ',izdebzone,' to ',izfinzone
  if(icodemat(imodnum) /= 2) then
      vpzone = cp(imodnum)
      vszone = cs(imodnum)
      print *,'Model # ',imodnum,' isotrope'
      print *,'vp = ',vpzone
      print *,'vs = ',vszone
      print *,'rho = ',rho(imodnum)
      print *,'Poisson''s ratio = ', &
        0.5d0*(vpzone*vpzone-2.d0*vszone*vszone)/(vpzone*vpzone-vszone*vszone)
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
   do i=ixdebzone,ixfinzone
      do j=izdebzone,izfinzone
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

  if (minval(num_modele) <= 0) stop 'Velocity model not entirely set...'

  close(10)

  print *
  print *,' Parameter file successfully read... '

! --------- fin lecture fichier parametres --------------

  allocate(psi(0:nx))
  allocate(eta(0:nz))
  allocate(absx(0:nx))
  allocate(a00(0:nz))
  allocate(a01(0:nz))
  allocate(valeta(0:nz))
  allocate(bot0(0:nx))
  allocate(top0(0:nx))

! calcul des points regulierement espaces
  do i=0,nx
    psi(i) = i/dble(nx)
  enddo

  do j=0,nz
    eta(j) = j/dble(nz)
  enddo

! quelques verifications de base a faire

  if(ngnod /= 4.and.ngnod /= 9) stop 'erreur ngnod different de 4 ou 9 !!'

  print *
  if(ngnod == 4) then
    print *,'Le maillage comporte ',nx,' x ',nz,' elements'
  else
    print *,'Le maillage comporte ',nx/2,' x ',nz/2,' elements'
  endif
  print *
  print *,'Les elements de controle sont des elements ',ngnod,' noeuds'
  print *

!------------------------------------------------------

  allocate(x(0:nx,0:nz))
  allocate(z(0:nx,0:nz))

  x(:,:)=0.d0
  z(:,:)=0.d0

! get topography data from external file
  print *,'Reading topography from file ',topofile
  open(unit=15,file=topofile,status='old')
  read(15,*) ntopo
  if (ntopo < 2) stop 'Not enough topography points (min 2)'
  print *,'Reading ',ntopo,' points from topography file'
  print *

  allocate(xtopo(ntopo))
  allocate(ztopo(ntopo))
  allocate(coefs_topo(ntopo))

  do i=1,ntopo
    read(15,*) xtopo(i),ztopo(i)
  enddo
  close(15)

! get interface data from external file, if any
  if(interffile /= 'none') then
  print *,'Reading interface from file ',interffile
  open(unit=15,file=interffile,status='old')
  read(15,*) ninterf
  if (ninterf < 2) stop 'Not enough interface points (min 2)'
  print *,'Reading ',ninterf,' points from interface file'

  allocate(xinterf(ninterf))
  allocate(zinterf(ninterf))
  allocate(coefs_interf(ninterf))

  do i=1,ninterf
    read(15,*) xinterf(i),zinterf(i)
  enddo

  close(15)
  else
    print *,'*** No interface file specified ***'
  endif

! check the values read
  print *
  print *, 'Topography data points (x,z)'
  print *, '----------------------------'
  print *
  print *, 'Topo 1 = (',xtopo(1),',',ztopo(1),')'
  print *, 'Topo ntopo = (',xtopo(ntopo),',',ztopo(ntopo),')'

!--- calculate the spline function for the topography
!--- imposer les tangentes aux deux bords
  tang1 = (ztopo(2)-ztopo(1))/(xtopo(2)-xtopo(1))
  tangN = (ztopo(ntopo)-ztopo(ntopo-1))/(xtopo(ntopo)-xtopo(ntopo-1))
  call spline(xtopo,ztopo,ntopo,tang1,tangN,coefs_topo)

!--- calculate the spline function for the interface
!--- imposer les tangentes aux deux bords
  if (interffile /= 'none') then
  tang1 = (zinterf(2)-zinterf(1))/(xinterf(2)-xinterf(1))
  tangN = (zinterf(ntopo)-zinterf(ntopo-1))/(xinterf(ntopo)-xinterf(ntopo-1))
  call spline(xinterf,zinterf,Ninterf,tang1,tangN,coefs_interf)
  endif

! *** afficher limites du modele lu
  print *
  print *, 'Limites absolues modele fichier topo :'
  print *
  print *, 'Xmin = ',minval(xtopo),'   Xmax = ',maxval(xtopo)
  print *, 'Zmin = ',minval(ztopo),'   Zmax = ',maxval(ztopo)
  print *

! *** eventuellement modifier sources si sources en surface
  print *
  print *, 'Position (x,z) de la source'
  print *
  if(isource_surf) zs = spl(xs,xtopo,ztopo,coefs_topo,ntopo)
  print *, 'Source = ',xs,zs

! *** eventuellement modifier recepteurs si enregistrement en surface
  print *
  print *, 'Position (x,z) des ',nrec,' receivers'
  print *
  do irec=1,nrec
   if(ienreg_surf) zrec(irec) = spl(xrec(irec),xtopo,ztopo,coefs_topo,ntopo)
   print *, 'Receiver ',irec,' = ',xrec(irec),zrec(irec)
  enddo

!--- definition du maillage suivant X
  do ix=0,nx
    absx(ix) = dens(ix,psi,xmin,xmax,nx)
  enddo

  if (interffile == 'none') then

! *** une seule zone si pas d'interface specifiee

  do iz=0,nz
    valeta(iz) = eta(iz)
    a00(iz) = 1-valeta(iz)
    a01(iz) = valeta(iz)
  enddo

  do ix=0,nx
    bot0(ix) = bottom(absx(ix))
    top0(ix) = spl(absx(ix),xtopo,ztopo,coefs_topo,ntopo)
  enddo

! valeurs de x et y pour display domaine physique
  do ix=0,nx
    do iz=0,nz
      x(ix,iz) = absx(ix)
      z(ix,iz) = a00(iz)*bot0(ix) + a01(iz)*top0(ix)
    enddo
  enddo

  else

! *** deux zones si topo

  ncut = nint(nz*ratio)

! *** ZONE DU BAS ***

  do iz=0,ncut
    valeta(iz) = dble(iz)/(nz*ratio)
    a00(iz) = 1-valeta(iz)
    a01(iz) = valeta(iz)
  enddo

  do ix=0,nx
    bot0(ix) = bottom(absx(ix))
    top0(ix) = spl(absx(ix),xinterf,zinterf,coefs_interf,ninterf)
  enddo

! valeurs de x et y pour display domaine physique
  do ix=0,nx
    do iz=0,ncut
      x(ix,iz) = absx(ix)
      z(ix,iz) = a00(iz)*bot0(ix) + a01(iz)*top0(ix)
    enddo
  enddo

! *** ZONE DU HAUT ***

  do iz=0,nz-ncut
    valeta(iz) = dble(iz)/(nz*(1.d0-ratio))
    a00(iz) = 1-valeta(iz)
    a01(iz) = valeta(iz)
  enddo

  do ix=0,nx
    bot0(ix) = spl(absx(ix),xinterf,zinterf,coefs_interf,ninterf)
    top0(ix) = spl(absx(ix),xtopo,ztopo,coefs_topo,ntopo)
  enddo

! valeurs de x et y pour display domaine physique
  do ix=0,nx
    do iz=0,nz-ncut
      x(ix,iz+ncut) = absx(ix)
      z(ix,iz+ncut) = a00(iz)*bot0(ix) + a01(iz)*top0(ix)
    enddo
  enddo

  endif

! calculer min et max de X et Z sur la grille
  print *
  print *, 'Valeurs min et max de X sur le maillage = ',minval(x),maxval(x)
  print *, 'Valeurs min et max de Z sur le maillage = ',minval(z),maxval(z)
  print *

! ***
! *** generer un fichier 'GNUPLOT' pour le controle de la grille ***
! ***

  print *
  print *,' Ecriture de la grille format GNUPLOT...'

  file1='gridfile.gnu'

 open(unit=20,file=file1,status='unknown')

! dessin de la topo de surface (splines)
  do i=0,nx-1
    write(20,15) sngl(absx(i)),sngl(top0(i))
    write(20,15) sngl(absx(i+1)),sngl(top0(i+1))
    write(20,10)
  enddo

! dessin de l'interface du milieu
  if (interffile /= 'none') then
  do i=0,nx-1
    write(20,15) sngl(absx(i)),sngl(bot0(i))
    write(20,15) sngl(absx(i+1)),sngl(bot0(i+1))
    write(20,10)
  enddo
  endif

! dessin des lignes horizontales de la grille
  print *, 'Ecriture lignes horizontales'
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
  print *, 'Ecriture lignes verticales'
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
  open(unit=20,file='plotgrid.gnu',status='unknown')
  write(20,*) 'set term postscript landscape monochrome solid "Helvetica" 22'
  write(20,*) 'set output "grille.ps"'
  write(20,*) 'plot "gridfile.gnu" title "Macroblocs mesh" w l'
  close(20)

  print *,' Fin ecriture de la grille format GNUPLOT'
  print *

! *** generation de la base de donnees

  open(unit=15,file='../SPECFEM90/DataBase',status='unknown')

  write(15,*) '#'
  write(15,*) '# DataBase for SPECFEM2D version 5.0'
  write(15,*) '# Dimitri Komatitsch, (c) University of Pau, France, May 2004'
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
  write(15,*) itaff,0,0

  write(15,*) 'imeshvect imodelvect iboundvect cutvect isubsamp'
  write(15,*) imeshvect,imodelvect,iboundvect,cutvect,isubsamp

  write(15,*) 'nrec anglerec'
  write(15,*) nrec,anglerec

  write(15,*) 'initialfield'
  write(15,*) initialfield

  write(15,*) 'isismostype ivecttype'
  write(15,*) isismostype,ivecttype

  write(15,*) 'ireadmodel ioutputgrid'
  write(15,*) ireadmodel,ioutputgrid

  write(15,*) 'ncycl dtinc'
  write(15,*) nt,dt

  write(15,*) 'source'
  write(15,*) '6 ',isource_type,xs-xmin,zs,f0,t0,factor,angle,' 0'

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

  netyp = 2

  write(15,*) 'netyp numat ngnod nspec iptsdisp ielemabs'
  write(15,*) netyp,nbmodeles,ngnod,nspec,iptsdisp,nelemabs

  write(15,*) 'Material sets (num 0 rho vp vs 0 0) or (num 2 rho c11 c13 c33 c44)'
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
  do iz=1,nzread
  do ix=1,nxread
    codehaut = .false.
    codebas = .false.
    codegauche = .false.
    codedroite = .false.
    inumelem = (iz-1)*nxread + ix
    if(abshaut   .and. iz==nzread) codehaut = .true.
    if(absbas    .and. iz== 1) codebas = .true.
    if(absgauche .and. ix== 1) codegauche = .true.
    if(absdroite .and. ix==nxread) codedroite = .true.
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

! *****************
! routines maillage
! *****************

! --- numero global du noeud

  integer function num(i,j,nx)
  implicit none
  integer i,j,nx

    num = j*(nx+1) + i + 1

  end function num

! ------- definition des fonctions representant les interfaces -------

!
! --- bas du modele
!

  double precision function bottom(x)
  implicit none
  double precision x

    bottom = 0.d0

  end function bottom

!
! --- representation interfaces par un spline
!

!--- spline

  double precision function spl(x,xtopo,ztopo,coefs,ntopo)

  implicit none
  integer ntopo
  double precision x,xp
  double precision xtopo(ntopo),ztopo(ntopo)
  double precision coefs(ntopo)

  spl = 0.
  xp = x
  if (xp < xtopo(1)) xp = xtopo(1)
  if (xp > xtopo(ntopo)) xp = xtopo(ntopo)
  call splint(xtopo,ztopo,coefs,ntopo,xp,spl)

  end function spl

! --- fonction de densification du maillage horizontal

  double precision function dens(ix,psi,xmin,xmax,nx)

  implicit none
  integer ix,nx
  double precision psi(0:nx)
  double precision xmin,xmax

  dens = xmin + dble(xmax-xmin)*psi(ix)

  end function dens

! --------------------------------------

! routine de calcul des coefs du spline (Numerical Recipes)

  subroutine spline(x,y,n,yp1,ypn,y2)

  implicit none

  integer n
  double precision x(n),y(n),y2(n)
  double precision, dimension(:), allocatable :: u
  double precision yp1,ypn

  integer i,k
  double precision sig,p,qn,un

  allocate(u(n))

  y2(1)=-0.5
  u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.
    y2(i)=(sig-1.)/p
    u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
               /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  enddo
  qn=0.5
  un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo

  deallocate(u)

  end subroutine spline

! --------------

! routine d'evaluation du spline (Numerical Recipes)

  SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)

  implicit none

  integer n
  double precision XA(N),YA(N),Y2A(N)
  double precision x,y

  integer k,klo,khi
  double precision h,a,b

  KLO=1
  KHI=N
  do while (KHI-KLO > 1)
    K=(KHI+KLO)/2
    IF(XA(K) > X)THEN
      KHI=K
    ELSE
      KLO=K
    ENDIF
  enddo
  H=XA(KHI)-XA(KLO)
  IF (H == 0.d0) stop 'Bad input in spline evaluation'
  A=(XA(KHI)-X)/H
  B=(X-XA(KLO))/H

  Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+ (B**3-B)*Y2A(KHI))*(H**2)/6.d0

  end subroutine SPLINT

