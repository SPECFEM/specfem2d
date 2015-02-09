!=====================================================================
!
!             P r e m a i l l e u r    F o r t r a n  9 0
!             -------------------------------------------
!
!                           Version 3.0
!                           -----------
!
!                         Dimitri Komatitsch
!    Department of Earth and Planetary Sciences - Harvard University
!
!                         (c) August 1998
!
!=====================================================================

!
! *** Version optimisee avec maillage non structure Jacques Muller - Elf ***
! *** Raffinement d'un facteur 3 en surface ***
!

  program maille_non_struct_3

  implicit none

! definir les tableaux pour allocation dynamique

! coordinates of the grid points
  double precision, allocatable :: x(:,:),z(:,:)

! variables needed to compute the transformation
  double precision, allocatable :: psi(:),eta(:),absx(:), &
      a00(:),a01(:),valeta(:),bot0(:),top0(:)

! stockage du modele de vitesse et densite
  double precision, allocatable :: rho(:),cp(:),cs(:)

! the topography data
  double precision, allocatable :: xtopo(:),ztopo(:),coefs_topo(:)

! arrays for the source
  double precision, allocatable :: xs(:),zs(:),f0(:),tshift_src(:),angle(:),factor(:)
  integer, allocatable :: isource_type(:),itimetype(:)

! arrays for the receivers
  double precision, allocatable :: xrec(:),zrec(:)

  character(len=50) interffile,topofile,title
  character(len=15) junk

  integer imatnum,inumabs,inumelem,nelemperio,nxgll,netyp
  integer icodehaut,icodebas,icodegauche,icodedroite
  integer nelemabs,npgeo,nspec,ntopo,nspecvolume,nspecWz
  integer k,ix,iz,irec,i,j,iadd
  integer imodele,nbmodeles,iaffinfo
  integer itaff,itfirstaff,pointsdisp,isubsamp,nrec,n1ana,n2ana
  integer irepr,nrec1,nrec2,isamp,nbsources,isismostype,ivecttype
  integer ngnod,nt,niter,idegpoly,nx,nz
  integer icodematread

  double precision valseuil,freqmaxrep,ratio
  double precision tang1,tangN
  double precision orig_x,orig_z,sizemax,cutvect,scalex,scalez
  double precision factorxsu,factorana,xspacerec,zspacerec
  double precision anglerec,anglerec2,xmin,xmax
  double precision xfin,zfin,xfin2,zfin2,xdeb,zdeb,xdeb2,zdeb2
  double precision alphanewm,betanewm,gammanewm,dt
  double precision rhoread,cpread,csread,aniso3read,aniso4read

  logical interpol,ignuplot,ireadmodel,iavs,ivisual3,ioutputgrid
  logical abshaut,absbas,absgauche,absdroite,absstacey
  logical periohaut,periogauche
  logical sismos,isources_surf,ienreg_surf,ienreg_surf2,display
  logical ivectplot,imeshvect
  logical topoplane,iexec,initialfield
  logical imodelvect,iboundvect,usletter,compenergy

  integer, external :: num
  double precision, external :: bottom,spl,dens

  double precision, parameter :: zero = 0.d0, one = 1.d0

! simulation a 2D
  integer, parameter :: ndime = 2
  integer, parameter :: ndofn = 2

! --- code des numeros d'aretes pour les bords absorbants
  integer, parameter :: iaretebas    = 1
  integer, parameter :: iaretedroite = 2
  integer, parameter :: iaretehaut   = 3
  integer, parameter :: iaretegauche = 4

! DK DK DK ajout Elf : extraction de la topo du fichier SEP
!!  call system('rm -f topo_from_SEP.dat topo_SEP_maille90.dat ; xextract_topo')

  print *
  print *,' *** Version optimisee avec maillage non structure ***'
  print *,' *** Raffinement d''un facteur 3 en surface ***'
  print *

! ***
! *** read the parameter file
! ***

  print *,' Reading the parameter file ... '
  print *

  open(unit=10,file='DATA/Par_file',status='old')

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
  read(10,2)junk,idegpoly
  read(10,2)junk,ngnod
  read(10,1)junk,ratio
  read(10,4)junk,topoplane
  read(10,4)junk,initialfield
  read(10,4)junk,ireadmodel
  read(10,4)junk,iexec

! DK DK forcer pour Elf
  ngnod = 9
  topoplane = .false.
  initialfield = .false.

! pour le non structure, verifier la coherence du maillage
  if(nx < 2) stop 'nx must be greater or equal to 2'
  if(nz < 2) stop 'nz must be greater or equal to 2'
  if(mod(nx,2) /= 0) stop 'nx must be even'

! multiplier par 3 pour implementer le deraffinement non conforme
  nx = nx * 3
  nz = nz * 3

! multiplier par 2 pour elements 9 noeuds
  nx = nx * 2
  nz = nz * 2

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read absorbing boundaries parameters
  read(10,4)junk,abshaut
  read(10,4)junk,absbas
  read(10,4)junk,absgauche
  read(10,4)junk,absdroite
  read(10,4)junk,absstacey
  read(10,4)junk,periohaut
  read(10,4)junk,periogauche

! DK DK forcer pour Elf
  abshaut = .false.
  periohaut = .false.
  periogauche = .false.

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read time step parameters
  read(10,2)junk,nt
  read(10,1)junk,dt
  read(10,2)junk,niter
  read(10,1)junk,alphanewm
  read(10,1)junk,betanewm
  read(10,1)junk,gammanewm

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read source parameters
  read(10,2)junk,nbsources
  read(10,4)junk,isources_surf
  read(10,1)junk,valseuil
  read(10,1)junk,freqmaxrep
  print *,'Nb de sources a lire : ',nbsources

  allocate(xs(nbsources))
  allocate(zs(nbsources))
  allocate(f0(nbsources))
  allocate(tshift_src(nbsources))
  allocate(isource_type(nbsources))
  allocate(itimetype(nbsources))
  allocate(angle(nbsources))
  allocate(factor(nbsources))

  do i=1,nbsources
      read(10,*)
      read(10,1)junk,xs(i)
      read(10,1)junk,zs(i)
      read(10,1)junk,f0(i)
      read(10,1)junk,tshift_src(i)
      read(10,2)junk,isource_type(i)
      read(10,2)junk,itimetype(i)
      read(10,1)junk,angle(i)
      read(10,1)junk,factor(i)

      print *
      print *,' Source #',i
      print *,'Position xs, zs = ',xs(i),zs(i)
      print *,'Frequency, delay = ',f0(i),tshift_src(i)
      print *,'Source type (1=force 2=explo) : ', &
                    isource_type(i)
      print *,'Angle of the source if force = ',angle(i)
      print *,'Multiplying factor = ',factor(i)
  enddo

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read receivers line parameters
  read(10,4)junk,sismos
  read(10,2)junk,isamp
  read(10,2)junk,isismostype
  read(10,2)junk,irepr
  read(10,*)
  read(10,2)junk,nrec1
  read(10,1)junk,xdeb
  read(10,1)junk,zdeb
  read(10,1)junk,xfin
  read(10,1)junk,zfin
  read(10,4)junk,ienreg_surf
  read(10,1)junk,anglerec
  read(10,*)
  read(10,2)junk,nrec2
  read(10,1)junk,xdeb2
  read(10,1)junk,zdeb2
  read(10,1)junk,xfin2
  read(10,1)junk,zfin2
  read(10,4)junk,ienreg_surf2
  read(10,1)junk,anglerec2
  read(10,*)
  read(10,1)junk,factorxsu
  read(10,2)junk,n1ana
  read(10,2)junk,n2ana
  read(10,1)junk,factorana

! determination et affichage position ligne de receivers
  if(nrec2 < 0) stop 'negative value of nrec2 !'

  if(nrec2 == 0) then
    nrec = nrec1
  else
    nrec = nrec1 + nrec2
  endif

! DK DK forcer pour Elf
  n1ana = 1
  n2ana = nrec

  allocate(xrec(nrec))
  allocate(zrec(nrec))

  if(nrec2 == 0) then
  print *
  print *,'There are ',nrec,' receivers on a single line'
  xspacerec=(xfin-xdeb)/dble(nrec-1)
  zspacerec=(zfin-zdeb)/dble(nrec-1)
  do i=1,nrec
     xrec(i) = xdeb + dble(i-1)*xspacerec
     zrec(i) = zdeb + dble(i-1)*zspacerec
  enddo
  else
  print *
  print *,'There are ',nrec,' receivers on two lines'
  print *,'First line contains ',nrec1,' receivers'
  print *,'Second line contains ',nrec2,' receivers'
  xspacerec=(xfin-xdeb)/dble(nrec1-1)
  zspacerec=(zfin-zdeb)/dble(nrec1-1)
  do i=1,nrec1
     xrec(i) = xdeb + dble(i-1)*xspacerec
     zrec(i) = zdeb + dble(i-1)*zspacerec
  enddo
  xspacerec=(xfin2-xdeb2)/dble(nrec2-1)
  zspacerec=(zfin2-zdeb2)/dble(nrec2-1)
  do i=1,nrec2
     xrec(i+nrec1) = xdeb2 + dble(i-1)*xspacerec
     zrec(i+nrec1) = zdeb2 + dble(i-1)*zspacerec
  enddo
  endif

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! read display parameters
  read(10,4)junk,display
  read(10,2)junk,itaff
  read(10,2)junk,itfirstaff
  read(10,2)junk,iaffinfo
  read(10,4)junk,ivectplot
  read(10,2)junk,ivecttype
  read(10,1)junk,cutvect
  read(10,4)junk,imeshvect
  read(10,4)junk,imodelvect
  read(10,4)junk,iboundvect
  read(10,4)junk,interpol
  read(10,2)junk,pointsdisp
  read(10,2)junk,isubsamp
  read(10,1)junk,scalex
  read(10,1)junk,scalez
  read(10,1)junk,sizemax
  read(10,4)junk,usletter
  read(10,1)junk,orig_x
  read(10,1)junk,orig_z
  read(10,4)junk,ignuplot
  read(10,4)junk,iavs
  read(10,4)junk,ivisual3
  read(10,4)junk,ioutputgrid
  read(10,4)junk,compenergy

! DK DK forcer pour Elf
  ignuplot = .false.
  iavs = .false.
  ivisual3 = .false.
  compenergy = .false.

! skip comment
  read(10,*)
  read(10,*)
  read(10,*)

! lecture des differents modeles de materiaux

  read(10,2)junk,nbmodeles
  if(nbmodeles <= 0) stop 'Negative number of models not allowed !!'

  allocate(rho(nbmodeles))
  allocate(cp(nbmodeles))
  allocate(cs(nbmodeles))

  rho(:) = 0.d0
  cp(:) = 0.d0
  cs(:) = 0.d0

  do imodele=1,nbmodeles
      read(10,*) i,icodematread,rhoread,cpread,csread,aniso3read,aniso4read
      if(i<1 .or. i>nbmodeles) stop 'Wrong material set number'
      rho(i) = rhoread
      cp(i) = cpread
      cs(i) = csread
      if (rho(i) < 0.d0 .or. cp(i) < 0.d0 .or. cs(i) < 0.d0) &
          stop 'Negative value of velocity or density'
  enddo

  print *
  print *, 'Nb de modeles de roche = ',nbmodeles
  print *
  do i=1,nbmodeles
      print *,'Modele #',i,' isotrope'
      print *,'rho,cp,cs = ',rho(i),cp(i),cs(i)
  enddo
  print *

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

  if(ngnod /= 9) stop 'erreur ngnod different de 9 !!'

! calcul du nombre total d'elements spectraux, absorbants et periodiques
  nspecvolume = (nx/2/3)*((nz-6)/2/3)
  nspecWz = 5*(nx/2/3)
  nspec = nspecvolume + nspecWz
  nelemperio = 0

  if(absgauche .or. absdroite .or. absbas) then
    nelemabs = 2 * (nz/6 - 2) + nx/6 + 3 + 3
  else
    nelemabs = 0
  endif

  print *
  print *,'Le maillage comporte ',nspec,' elements spectraux (nx = ',nx/6, &
     ' nz = ',nz/6,')'
  print *,'soit ',nspecvolume,' elements spectraux dans le volume'
  print *,'et ',nspecWz,' elements spectraux dans la couche Wz'
  print *,'Chaque element comporte ',idegpoly+1,' points dans chaque direction'
  print *,'Le nombre maximum de points theorique est ',nspec*(idegpoly+1)**ndime
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

! *** afficher limites du modele lu
  print *
  print *, 'Limites absolues modele fichier topo :'
  print *
  print *, 'Xmin = ',minval(xtopo),'   Xmax = ',maxval(xtopo)
  print *, 'Zmin = ',minval(ztopo),'   Zmax = ',maxval(ztopo)
  print *

! *** modifier sources pour position par rapport a la surface
  print *
  print *, 'Position (x,z) des ',nbsources,' sources'
  print *
  do i=1,nbsources

! DK DK DK Elf : position source donnee en profondeur par rapport a la topo
   zs(i) = spl(xs(i),xtopo,ztopo,coefs_topo,ntopo) - zs(i)

   if(isources_surf) zs(i) = spl(xs(i),xtopo,ztopo,coefs_topo,ntopo)
   print *, 'Source ',i,' = ',xs(i),zs(i)
  enddo

! *** modifier recepteurs pour enregistrement en surface
  print *
  print *, 'Position (x,z) des ',nrec,' receivers'
  print *
  do irec=1,nrec

! DK DK DK Elf : distinguer les deux lignes de recepteurs
  if(irec <= nrec1) then
   if(ienreg_surf) zrec(irec) = spl(xrec(irec),xtopo,ztopo,coefs_topo,ntopo)
  else
   if(ienreg_surf2) zrec(irec) = spl(xrec(irec),xtopo,ztopo,coefs_topo,ntopo)
  endif
   print *, 'Receiver ',irec,' = ',xrec(irec),zrec(irec)

  enddo

!--- definition du maillage suivant X
  do ix=0,nx
          absx(ix) = dens(ix,psi,xmin,xmax,nx)
  enddo

! *** une seule zone

  do iz=0,nz

! DK DK DK densification sinusoidale ici en vertical
  valeta(iz) = eta(iz) + ratio * sin(3.14159265 * eta(iz))
  if(valeta(iz) < zero) valeta(iz) = zero
  if(valeta(iz) > one ) valeta(iz) = one
! DK DK DK densification sinusoidale ici en vertical

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

! calculer min et max de X et Z sur la grille
  print *
  print *, 'Valeurs min et max de X sur le maillage = ',minval(x),maxval(x)
  print *, 'Valeurs min et max de Z sur le maillage = ',minval(z),maxval(z)
  print *

! *** generation de la base de donnees

  print *
  print *,' Creation de la base de donnees pour SPECFEM...'
  print *

  open(unit=15,file='../SPECFEM90/DataBase',status='unknown')

  write(15,*) '#'
  write(15,*) '# Base de Donnees pour Specfem - Premailleur Fortran 90'
  write(15,*) '# ',title
  write(15,*) '# Dimitri Komatitsch, (c) EPS - Harvard August 1998'
  write(15,*) '#'

  write(15,*) 'Titre simulation'
  write(15,40) title

  npgeo = (nx+1)*(nz+1)
  write(15,*) 'ndofn ndime npgeo'
  write(15,*) ndofn,ndime,npgeo

  write(15,*) 'display ignuplot interpol'
  write(15,*) display,ignuplot,interpol

  write(15,*) 'itaff itfirstaff icolor inumber'
  write(15,*) itaff,itfirstaff,0,0

  write(15,*) 'ivectplot imeshvect imodelvect iboundvect cutvect isubsamp'
  write(15,*) ivectplot,imeshvect,imodelvect,iboundvect,cutvect,isubsamp

  write(15,*) 'scalex scalez sizemax angle rapport USletter'
  write(15,*) scalex,scalez,sizemax,20.,0.40,usletter

  write(15,*) 'orig_x orig_z isymbols'
  write(15,*) orig_x,orig_z,' T'

  write(15,*) 'valseuil freqmaxrep'
  write(15,*) valseuil,freqmaxrep

  write(15,*) 'sismos nrec nrec1 nrec2 isamp'
  write(15,*) sismos,nrec,nrec1,nrec2,isamp

  write(15,*) 'irepr anglerec anglerec2'
  write(15,*) irepr,anglerec,anglerec2

  write(15,*) 'topoplane absstacey compenergy'
  write(15,*) topoplane,absstacey,compenergy

  write(15,*) 'initialfield factorana factorxsu n1ana n2ana'
  write(15,*) initialfield,factorana,factorxsu,n1ana,n2ana

  write(15,*) 'isismostype ivecttype iaffinfo'
  write(15,*) isismostype,ivecttype,iaffinfo

  write(15,*) 'ireadmodel ioutputgrid iavs ivisual3'
  write(15,*) ireadmodel,ioutputgrid,iavs,ivisual3

  write(15,*) 'iexec iecho'
  if(iexec) then
    write(15,*) '1       1'
  else
    write(15,*) '0       1'
  endif

  write(15,*) 'ncycl dtinc niter'
  write(15,*) nt,dt,niter

  write(15,*) 'alpha beta gamma (alpha not used for the moment)'
  write(15,*) alphanewm,betanewm,gammanewm

  write(15,*) 'nltfl (number of force or pressure sources)'
  write(15,*) nbsources

  write(15,*) 'Collocated forces and/or pressure sources:'
  do i=1,nbsources
      write(15,*) itimetype(i),isource_type(i), &
         xs(i)-xmin ,zs(i), &
        f0(i),tshift_src(i),factor(i),angle(i),0
  enddo

  write(15,*) 'Receivers positions:'
  do irec=1,nrec
      write(15,*) irec,xrec(irec)-xmin ,zrec(irec)
  enddo

  write(15,*) 'Coordinates of macroblocs mesh (coorg):'
  do j=0,nz
      do i=0,nx
      write(15,*) num(i,j,nx),x(i,j)-xmin,z(i,j)
      enddo
  enddo

  netyp = 2
  nxgll = idegpoly + 1

  write(15,*) 'netyp numat ngnod nxgll nygll nspec pointsdisp ielemabs ielemperio'
  write(15,*) netyp,nbmodeles,ngnod,nxgll,nxgll,nspec,pointsdisp, &
                nelemabs,nelemperio

  write(15,*) 'Material sets (num 0 rho vp vs 0 0)'
  do i=1,nbmodeles
       write(15,*) i,0,rho(i),cp(i),cs(i),0,0
  enddo


  write(15,*) 'Arrays kmato and knods for each bloc:'

  imatnum = 1
  k=0

! zone structuree dans le volume
  do j=0,nz-12,6
  do i=0,nx-6,6
      k = k + 1
      write(15,*) k,imatnum,num(i,j,nx),num(i+6,j,nx),num(i+6,j+6,nx), &
              num(i,j+6,nx),num(i+3,j,nx),num(i+6,j+3,nx), &
              num(i+3,j+6,nx),num(i,j+3,nx),num(i+3,j+3,nx)
  enddo
  enddo

  if(k /= nspecvolume) stop 'nombre d''elements incoherent dans le volume'

! zone non structuree dans la couche Wz
  j=nz-6
  do i=0,nx-12,12

! element 1 du raccord
      k = k + 1
      write(15,*) k,imatnum,num(i,j,nx),num(i+6,j,nx),num(i+4,j+2,nx), &
              num(i,j+2,nx),num(i+3,j,nx),num(i+5,j+1,nx), &
              num(i+2,j+2,nx),num(i,j+1,nx),num(i+3,j+1,nx)

! element 2 du raccord
      k = k + 1
      write(15,*) k,imatnum,num(i,j+2,nx),num(i+4,j+2,nx),num(i+2,j+4,nx), &
              num(i,j+4,nx),num(i+2,j+2,nx),num(i+3,j+3,nx), &
              num(i+1,j+4,nx),num(i,j+3,nx),num(i+1,j+3,nx)

! element 3 du raccord
      k = k + 1
      write(15,*) k,imatnum,num(i,j+4,nx),num(i+2,j+4,nx),num(i+2,j+6,nx), &
              num(i,j+6,nx),num(i+1,j+4,nx),num(i+2,j+5,nx), &
              num(i+1,j+6,nx),num(i,j+5,nx),num(i+1,j+5,nx)

! element 4 du raccord
      k = k + 1
      write(15,*) k,imatnum,num(i+2,j+4,nx),num(i+4,j+2,nx),num(i+4,j+6,nx), &
              num(i+2,j+6,nx),num(i+3,j+3,nx),num(i+4,j+4,nx), &
              num(i+3,j+6,nx),num(i+2,j+5,nx),num(i+3,j+5,nx)

! element 5 du raccord
      k = k + 1
      write(15,*) k,imatnum,num(i+4,j+2,nx),num(i+6,j,nx),num(i+6,j+6,nx), &
              num(i+4,j+6,nx),num(i+5,j+1,nx),num(i+6,j+3,nx), &
              num(i+5,j+6,nx),num(i+4,j+4,nx),num(i+5,j+3,nx)

! element 6 du raccord
      k = k + 1
      write(15,*) k,imatnum,num(i+6,j,nx),num(i+8,j+2,nx),num(i+8,j+6,nx), &
              num(i+6,j+6,nx),num(i+7,j+1,nx),num(i+8,j+4,nx), &
              num(i+7,j+6,nx),num(i+6,j+3,nx),num(i+7,j+3,nx)

! element 7 du raccord
      k = k + 1
      write(15,*) k,imatnum,num(i+8,j+2,nx),num(i+10,j+4,nx),num(i+10,j+6,nx), &
              num(i+8,j+6,nx),num(i+9,j+3,nx),num(i+10,j+5,nx), &
              num(i+9,j+6,nx),num(i+8,j+4,nx),num(i+9,j+4,nx)

! element 8 du raccord
      k = k + 1
      write(15,*) k,imatnum,num(i+6,j,nx),num(i+12,j,nx),num(i+12,j+2,nx), &
              num(i+8,j+2,nx),num(i+9,j,nx),num(i+12,j+1,nx), &
              num(i+10,j+2,nx),num(i+7,j+1,nx),num(i+10,j+1,nx)

! element 9 du raccord
      k = k + 1
      write(15,*) k,imatnum,num(i+8,j+2,nx),num(i+12,j+2,nx),num(i+12,j+4,nx), &
              num(i+10,j+4,nx),num(i+10,j+2,nx),num(i+12,j+3,nx), &
              num(i+11,j+4,nx),num(i+9,j+3,nx),num(i+11,j+3,nx)

! element 10 du raccord
      k = k + 1
      write(15,*) k,imatnum,num(i+10,j+4,nx),num(i+12,j+4,nx),num(i+12,j+6,nx),&
              num(i+10,j+6,nx),num(i+11,j+4,nx),num(i+12,j+5,nx), &
              num(i+11,j+6,nx),num(i+10,j+5,nx),num(i+11,j+5,nx)

  enddo

  if(k /= nspec) stop 'nombre d''elements incoherent dans la couche Wz'

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
  print *,'Stacey = ',absstacey
  print *

! generer la liste des elements absorbants
  if(nelemabs > 0) then
  write(15,*) 'Liste des elements absorbants (haut bas gauche droite) :'

! repasser aux vrais valeurs de nx et nz
  nx = nx / 6
  nz = nz / 6

  inumabs = 0

! bord absorbant du bas sans les coins
  iz = 1
  do ix = 2,nx-1
    inumabs = inumabs + 1
    inumelem = (iz-1)*nx + ix
    icodehaut = 0
    icodebas = iaretebas
    icodegauche = 0
    icodedroite = 0
    write(15,*) inumabs,inumelem,icodehaut,icodebas,icodegauche,icodedroite
  enddo

! coin en bas a gauche
    inumabs = inumabs + 1
    inumelem = 1
    icodehaut = 0
    icodebas = iaretebas
    icodegauche = iaretegauche
    icodedroite = 0
    write(15,*) inumabs,inumelem,icodehaut,icodebas,icodegauche,icodedroite

! coin en bas a droite
    inumabs = inumabs + 1
    inumelem = nx
    icodehaut = 0
    icodebas = iaretebas
    icodegauche = 0
    icodedroite = iaretedroite
    write(15,*) inumabs,inumelem,icodehaut,icodebas,icodegauche,icodedroite

! partie structuree du bord de gauche
  ix = 1
  do iz = 2,nz-1
    inumabs = inumabs + 1
    inumelem = (iz-1)*nx + ix
    icodehaut = 0
    icodebas = 0
    icodegauche = iaretegauche
    icodedroite = 0
    write(15,*) inumabs,inumelem,icodehaut,icodebas,icodegauche,icodedroite
  enddo

! partie structuree du bord de droite
  ix = nx
  do iz = 2,nz-1
    inumabs = inumabs + 1
    inumelem = (iz-1)*nx + ix
    icodehaut = 0
    icodebas = 0
    icodegauche = 0
    icodedroite = iaretedroite
    write(15,*) inumabs,inumelem,icodehaut,icodebas,icodegauche,icodedroite
  enddo

! partie non structuree du bord de gauche (trois elements)
  do iadd = 1,3
    inumabs = inumabs + 1
    inumelem = nx*(nz-1) + iadd
    icodehaut = 0
    icodebas = 0
    icodegauche = iaretegauche
    icodedroite = 0
    write(15,*) inumabs,inumelem,icodehaut,icodebas,icodegauche,icodedroite
  enddo

! partie non structuree du bord de droite (trois elements)
  do iadd = 1,3
    inumabs = inumabs + 1
    inumelem = nspec - iadd + 1
    icodehaut = 0
    icodebas = 0
    icodegauche = 0
    icodedroite = iaretedroite
    write(15,*) inumabs,inumelem,icodehaut,icodebas,icodegauche,icodedroite
  enddo

  if(inumabs /= nelemabs) stop 'nombre d''elements absorbants incoherent'

  endif

! fermer la base de donnees

  close(15)

 40 format(a50)

  end program maille_non_struct_3

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
      endif
  enddo
  H=XA(KHI)-XA(KLO)
  IF (H == 0.d0) stop 'Bad input in spline evaluation'
  A=(XA(KHI)-X)/H
  B=(X-XA(KLO))/H

  Y=A*YA(KLO)+B*YA(KHI)+((A**3-A)*Y2A(KLO)+ (B**3-B)*Y2A(KHI))*(H**2)/6.d0

  end subroutine SPLINT

