
!=====================================================================
!
!                  P r e m a i l l e u r - 2 D
!                  ---------------------------
!
!                         Version 2.1
!                         -----------
!
!                       Dimitri Komatitsch
!
!                    Departement de Sismologie
!              Institut de Physique du Globe de Paris
!
!       (c) Institut de Physique du Globe de Paris, Octobre 1996
!
!=====================================================================

! program to mesh a 2D canyon for a study with Paco Sanchez-Sesma and comparison with Kawase (1988)

! DK DK Mexico August 1999 : mise a jour format base de donnees

  program circular_canyon

  implicit none

! max size of the model in elements
  integer, parameter :: mnx=7,mnz=7

  double precision, parameter :: pi=3.141592653589793d0

! seuil pour considerer deux points comme confondus
  double precision, parameter :: rseuil=1.d-2

! declare variables
  integer imaxabs,n2ana,itimetype,isource_type,nump1,nump2,nump3,nump4
  integer ndofn,ndime,ngnod,nnode,nbcnd,n1ana
  integer nofst,npgeo,nspel,nbmodeles,nbsources,nrec,lquad,isamp,nrec1,nrec2
  integer irec,imatnum,netyp,nxgll,nelemperio,nelemabs,nx,nz,i,j
  integer irepr,nrecsur3,nt,niter,itaff,itfirstaff,numerocourant,pointsdisp,isubsamp

  double precision R,theta_i,theta_init,delta_theta,eta_j,valseuil,freqmaxrep
  double precision f0,t0,xs,zs,angle,factor,dist,xoffs,zoffs
  double precision xrec,zrec,rho,cp,cs,anglerec
  double precision anglerec2,dt,alphanewm,betanewm,gammanewm
  double precision cutvect,cutcolor,scalex,scalez,sizemax,orig_x,orig_z
  double precision factorana,factorxsu

! stockage de la grille curvi (x et z)
  integer, parameter :: npoinz1=(4*mnx+1)*(mnz+1), nelemz1=(4*mnx)*mnz
  double precision x1(0:4*mnx,0:mnz)
  double precision z1(0:4*mnx,0:mnz)

  integer, parameter :: npoinz3=(2*mnx+1)*(4*mnz+1), nelemz3=(2*mnx)*(4*mnz)
  double precision x3(0:2*mnx,0:4*mnz)
  double precision z3(0:2*mnx,0:4*mnz)

  integer, parameter :: npoinz4=(2*mnx+1)*(2*mnz+1), nelemz4=(2*mnx)*(2*mnz)
  double precision x4(0:2*mnx,0:2*mnz)
  double precision z4(0:2*mnx,0:2*mnz)

  integer, parameter :: npoinz1b=(2*mnx+1)*(mnz+1), nelemz1b=(2*mnx)*mnz
  double precision x1b(0:2*mnx,0:mnz)
  double precision z1b(0:2*mnx,0:mnz)

  integer, parameter :: npoinz2b=(mnx+1)*(2*mnz+1), nelemz2b=mnx*(2*mnz)
  double precision x2b(0:mnx,0:2*mnz)
  double precision z2b(0:mnx,0:2*mnz)

  integer, parameter :: npoinz3b=(4*mnx+1)*(4*mnz+1), nelemz3b=(4*mnx)*(4*mnz)
  double precision x3b(0:4*mnx,0:4*mnz)
  double precision z3b(0:4*mnx,0:4*mnz)

  integer, parameter :: npoinz4b=(2*mnx+1)*(2*mnz+1), nelemz4b=(2*mnx)*(2*mnz)
  double precision x4b(0:2*mnx,0:2*mnz)
  double precision z4b(0:2*mnx,0:2*mnz)

! nombre max de points de maillage, et nombre exact d'elements
  integer, parameter :: npoin = npoinz1+npoinz3+npoinz4+npoinz1b+npoinz2b+npoinz3b+npoinz4b
  integer, parameter :: nelem = nelemz1+nelemz3+nelemz4+nelemz1b+nelemz2b+nelemz3b+nelemz4b

! coordonnees geometriques des points
  double precision xpoint(npoin)
  double precision zpoint(npoin)

! coordonnees des sommets de chaque element
  double precision x1e(nelem)
  double precision z1e(nelem)
  double precision x2e(nelem)
  double precision z2e(nelem)
  double precision x3e(nelem)
  double precision z3e(nelem)
  double precision x4e(nelem)
  double precision z4e(nelem)

! numero des points des elements
  integer numpoin1(nelem)
  integer numpoin2(nelem)
  integer numpoin3(nelem)
  integer numpoin4(nelem)

  character(len=50) title

  logical iexternal, aleatoire, topoplane, simulate, absstacey
  logical absorbhaut, absorbbas, absorbgauche, sismos
  logical absorbdroite, absorbstacey, absorbmodar, ifullabs

  logical display, ignuplot, ivectplot, icolorplot, imeshvect
  logical imeshcolor, imodelvect, iboundvect, interpol, isymbols, initialfield
  logical usletter,compenergy

  print *,'Number of elements = ',nelem
  print *,'Max number of points = ',npoin

  nx = mnx
  nz = mnz

  R = 1.

! ***************************************
! *** ZONE DE DROITE
! ***************************************

! generer les points de base de l'interpolation lineaire (zone 1)
  theta_init = 3 * pi / 2.
  delta_theta = pi / 2.
  do i=0,4*nx

! --- point de depart
  if (i < 2*nx) then
      x1(i,0) = 2.*R * dble(i) / dble(2*nx)
      z1(i,0) = - 2.*R
  else
      x1(i,0) = 2.*R
      z1(i,0) = - 2.*R * (1. - dble(i - 2*nx) / dble(2*nx))
  endif

! --- point d'arrivee
      theta_i = theta_init + delta_theta * dble(i) / dble(4*nx)
      x1(i,nz) = cos(theta_i)
      z1(i,nz) = sin(theta_i)

! --- points intermediaires par interpolation lineaire
      do j=1,nz-1
        eta_j = dble(j) / dble(nz)
        x1(i,j) = (1.-eta_j)*x1(i,0) + eta_j*x1(i,nz)
        z1(i,j) = (1.-eta_j)*z1(i,0) + eta_j*z1(i,nz)
      enddo
  enddo

! generer zone de gauche (zone 3)
  do i=0,2*nx
    do j=0,4*nz
      x3(i,j) = 5. * dble(i) / dble(2*nx) + 2.
      if (j <= 2*nz) then
        z3(i,j) = 7. * dble(j) / dble(2*nz) - 9.
      else
        z3(i,j) = 2. * dble(j-2*nz) / dble(2*nz) - 2.
      endif
    enddo
  enddo

! generer zone du bas (zone 4)
  do i=0,2*nx
    do j=0,2*nz
      x4(i,j) = 2. * dble(i) / dble(2*nx)
      z4(i,j) = 7. * dble(j) / dble(2*nz) - 9.
    enddo
  enddo

! ***************************************
! *** ZONE DE GAUCHE
! ***************************************

! generer les points de base de l'interpolation lineaire (zone 1)
  theta_init = pi / 4.
  delta_theta = pi / 4.

  do i=0,2*nx
! --- point de depart
      x1b(i,0) = 2.*R * (dble(i) / dble(2*nx) - 1.)
      z1b(i,0) = - 2.*R

! --- point d'arrivee
      theta_i = theta_init + delta_theta * dble(i) / dble(2*nx)
      x1b(i,nz) = - cos(theta_i)
      z1b(i,nz) = - sin(theta_i)

! --- points intermediaires par interpolation lineaire
      do j=1,nz-1
        eta_j = dble(j) / dble(nz)
        x1b(i,j) = (1.-eta_j)*x1b(i,0) + eta_j*x1b(i,nz)
        z1b(i,j) = (1.-eta_j)*z1b(i,0) + eta_j*z1b(i,nz)
      enddo
  enddo

! generer les points de base de l'interpolation lineaire (zone 2)
  theta_init = pi / 4.

  do j=0,2*nz
! --- point de depart
      x2b(0,j) = - 2.*R
      z2b(0,j) = 2.*R * (dble(j) / dble(2*nz) - 1.)

! --- point d'arrivee
      theta_i = theta_init - delta_theta * dble(j) / dble(2*nz)
      x2b(nx,j) = - cos(theta_i)
      z2b(nx,j) = - sin(theta_i)

! --- points intermediaires par interpolation lineaire
      do i=1,nx-1
            eta_j = dble(i) / dble(nx)
            x2b(i,j) = (1.-eta_j)*x2b(0,j) + eta_j*x2b(nx,j)
            z2b(i,j) = (1.-eta_j)*z2b(0,j) + eta_j*z2b(nx,j)
      enddo

  enddo

! generer zone de gauche (zone 3)
  do i=0,4*nx
    do j=0,4*nz
      x3b(i,j) = 10. * dble(i) / dble(4*nx) - 12.
      if (j <= 2*nz) then
        z3b(i,j) = 7. * dble(j) / dble(2*nz) - 9.
      else
        z3b(i,j) = 2. * dble(j-2*nz) / dble(2*nz) - 2.
      endif
    enddo
  enddo

! generer zone du bas (zone 4)
  do i=0,2*nx
    do j=0,2*nz
      x4b(i,j) = 2. * dble(i) / dble(2*nx) - 2.
      z4b(i,j) = 7. * dble(j) / dble(2*nz) - 9.
    enddo
  enddo

! ***
! *** generer un fichier 'GNUPLOT' pour le controle de la grille ***
! ***

  write(*,*)
  write(*,*) 'Ecriture de la grille format GNUPLOT...'

  open(unit=20,file='grid.gnu',status='unknown')

! *** dessiner la zone 1
  do j=0,nz
    do i=0,4*nx-1
      write(20,*) sngl(x1(i,j)),sngl(z1(i,j))
      write(20,*) sngl(x1(i+1,j)),sngl(z1(i+1,j))
      write(20,100)
    enddo
  enddo

  do i=0,4*nx
    do j=0,nz-1
      write(20,*) sngl(x1(i,j)),sngl(z1(i,j))
      write(20,*) sngl(x1(i,j+1)),sngl(z1(i,j+1))
      write(20,100)
    enddo
  enddo

! *** dessiner la zone 3
  do j=0,4*nz
    do i=0,2*nx-1
      write(20,*) sngl(x3(i,j)),sngl(z3(i,j))
      write(20,*) sngl(x3(i+1,j)),sngl(z3(i+1,j))
      write(20,100)
    enddo
  enddo

  do i=0,2*nx
    do j=0,4*nz-1
      write(20,*) sngl(x3(i,j)),sngl(z3(i,j))
      write(20,*) sngl(x3(i,j+1)),sngl(z3(i,j+1))
      write(20,100)
    enddo
  enddo

! *** dessiner la zone 4
  do j=0,2*nz
    do i=0,2*nx-1
      write(20,*) sngl(x4(i,j)),sngl(z4(i,j))
      write(20,*) sngl(x4(i+1,j)),sngl(z4(i+1,j))
      write(20,100)
    enddo
  enddo

  do i=0,2*nx
    do j=0,2*nz-1
      write(20,*) sngl(x4(i,j)),sngl(z4(i,j))
      write(20,*) sngl(x4(i,j+1)),sngl(z4(i,j+1))
      write(20,100)
    enddo
  enddo

! *** dessiner la zone 1
  do j=0,nz
    do i=0,2*nx-1
      write(20,*) sngl(x1b(i,j)),sngl(z1b(i,j))
      write(20,*) sngl(x1b(i+1,j)),sngl(z1b(i+1,j))
      write(20,100)
    enddo
  enddo

  do i=0,2*nx
    do j=0,nz-1
      write(20,*) sngl(x1b(i,j)),sngl(z1b(i,j))
      write(20,*) sngl(x1b(i,j+1)),sngl(z1b(i,j+1))
      write(20,100)
    enddo
  enddo

! *** dessiner la zone 2
  do j=0,2*nz
    do i=0,nx-1
      write(20,*) sngl(x2b(i,j)),sngl(z2b(i,j))
      write(20,*) sngl(x2b(i+1,j)),sngl(z2b(i+1,j))
      write(20,100)
    enddo
  enddo

  do i=0,nx
    do j=0,2*nz-1
      write(20,*) sngl(x2b(i,j)),sngl(z2b(i,j))
      write(20,*) sngl(x2b(i,j+1)),sngl(z2b(i,j+1))
      write(20,100)
    enddo
  enddo

! *** dessiner la zone 3
  do j=0,4*nz
    do i=0,4*nx-1
      write(20,*) sngl(x3b(i,j)),sngl(z3b(i,j))
      write(20,*) sngl(x3b(i+1,j)),sngl(z3b(i+1,j))
      write(20,100)
    enddo
  enddo

  do i=0,4*nx
    do j=0,4*nz-1
      write(20,*) sngl(x3b(i,j)),sngl(z3b(i,j))
      write(20,*) sngl(x3b(i,j+1)),sngl(z3b(i,j+1))
      write(20,100)
    enddo
  enddo

! *** dessiner la zone 4
  do j=0,2*nz
    do i=0,2*nx-1
      write(20,*) sngl(x4b(i,j)),sngl(z4b(i,j))
      write(20,*) sngl(x4b(i+1,j)),sngl(z4b(i+1,j))
      write(20,100)
    enddo
  enddo

  do i=0,2*nx
    do j=0,2*nz-1
      write(20,*) sngl(x4b(i,j)),sngl(z4b(i,j))
      write(20,*) sngl(x4b(i,j+1)),sngl(z4b(i,j+1))
      write(20,100)
    enddo
  enddo

  close(20)

  write(*,*) 'Fin ecriture de la grille format GNUPLOT'
  write(*,*)

 100  format('')

! ***
! *** generer la liste des points geometriques
! ***

  numerocourant = 1

! *** zone 1
  do j=0,nz
    do i=0,4*nx
      xpoint(numerocourant) = x1(i,j)
      zpoint(numerocourant) = z1(i,j)
      numerocourant = numerocourant + 1
    enddo
  enddo

! *** zone 3
  do j=0,4*nz
    do i=0,2*nx
      xpoint(numerocourant) = x3(i,j)
      zpoint(numerocourant) = z3(i,j)
      numerocourant = numerocourant + 1
    enddo
  enddo

! *** zone 4
  do j=0,2*nz
    do i=0,2*nx
      xpoint(numerocourant) = x4(i,j)
      zpoint(numerocourant) = z4(i,j)
      numerocourant = numerocourant + 1
    enddo
  enddo

! *** zone 1
  do j=0,nz
    do i=0,2*nx
      xpoint(numerocourant) = x1b(i,j)
      zpoint(numerocourant) = z1b(i,j)
      numerocourant = numerocourant + 1
    enddo
  enddo

! *** zone 2
  do j=0,2*nz
    do i=0,nx
      xpoint(numerocourant) = x2b(i,j)
      zpoint(numerocourant) = z2b(i,j)
      numerocourant = numerocourant + 1
    enddo
  enddo

! *** zone 3
  do j=0,4*nz
    do i=0,4*nx
      xpoint(numerocourant) = x3b(i,j)
      zpoint(numerocourant) = z3b(i,j)
      numerocourant = numerocourant + 1
    enddo
  enddo

! *** zone 4
  do j=0,2*nz
    do i=0,2*nx
      xpoint(numerocourant) = x4b(i,j)
      zpoint(numerocourant) = z4b(i,j)
      numerocourant = numerocourant + 1
    enddo
  enddo

  print *,'nb de points stockes = ',numerocourant - 1

! ***
! *** generer la liste des elements
! ***

  numerocourant = 1
  imaxabs = 0

! *** zone 1
  do j=0,nz-1
    do i=0,4*nx-1
      x1e(numerocourant) = x1(i,j)
      z1e(numerocourant) = z1(i,j)
      x2e(numerocourant) = x1(i+1,j)
      z2e(numerocourant) = z1(i+1,j)
      x3e(numerocourant) = x1(i+1,j+1)
      z3e(numerocourant) = z1(i+1,j+1)
      x4e(numerocourant) = x1(i,j+1)
      z4e(numerocourant) = z1(i,j+1)
      numerocourant = numerocourant + 1
    enddo
  enddo

! *** zone 3
  do j=0,4*nz-1
    do i=0,2*nx-1
      x1e(numerocourant) = x3(i,j)
      z1e(numerocourant) = z3(i,j)
      x2e(numerocourant) = x3(i+1,j)
      z2e(numerocourant) = z3(i+1,j)
      x3e(numerocourant) = x3(i+1,j+1)
      z3e(numerocourant) = z3(i+1,j+1)
      x4e(numerocourant) = x3(i,j+1)
      z4e(numerocourant) = z3(i,j+1)
      numerocourant = numerocourant + 1
    enddo
  enddo

! *** zone 4
  do j=0,2*nz-1
    do i=0,2*nx-1
      x1e(numerocourant) = x4(i,j)
      z1e(numerocourant) = z4(i,j)
      x2e(numerocourant) = x4(i+1,j)
      z2e(numerocourant) = z4(i+1,j)
      x3e(numerocourant) = x4(i+1,j+1)
      z3e(numerocourant) = z4(i+1,j+1)
      x4e(numerocourant) = x4(i,j+1)
      z4e(numerocourant) = z4(i,j+1)
      numerocourant = numerocourant + 1
    enddo
  enddo

! *** zone 1
  do j=0,nz-1
    do i=0,2*nx-1
      x1e(numerocourant) = x1b(i,j)
      z1e(numerocourant) = z1b(i,j)
      x2e(numerocourant) = x1b(i+1,j)
      z2e(numerocourant) = z1b(i+1,j)
      x3e(numerocourant) = x1b(i+1,j+1)
      z3e(numerocourant) = z1b(i+1,j+1)
      x4e(numerocourant) = x1b(i,j+1)
      z4e(numerocourant) = z1b(i,j+1)
      numerocourant = numerocourant + 1
    enddo
  enddo

! *** zone 2
  do j=0,2*nz-1
    do i=0,nx-1
      x1e(numerocourant) = x2b(i,j)
      z1e(numerocourant) = z2b(i,j)
      x2e(numerocourant) = x2b(i+1,j)
      z2e(numerocourant) = z2b(i+1,j)
      x3e(numerocourant) = x2b(i+1,j+1)
      z3e(numerocourant) = z2b(i+1,j+1)
      x4e(numerocourant) = x2b(i,j+1)
      z4e(numerocourant) = z2b(i,j+1)
      numerocourant = numerocourant + 1
    enddo
  enddo

! *** zone 3
  do j=0,4*nz-1
    do i=0,4*nx-1
      x1e(numerocourant) = x3b(i,j)
      z1e(numerocourant) = z3b(i,j)
      x2e(numerocourant) = x3b(i+1,j)
      z2e(numerocourant) = z3b(i+1,j)
      x3e(numerocourant) = x3b(i+1,j+1)
      z3e(numerocourant) = z3b(i+1,j+1)
      x4e(numerocourant) = x3b(i,j+1)
      z4e(numerocourant) = z3b(i,j+1)
      numerocourant = numerocourant + 1
    enddo
  enddo

! *** zone 4
  do j=0,2*nz-1
    do i=0,2*nx-1
      x1e(numerocourant) = x4b(i,j)
      z1e(numerocourant) = z4b(i,j)
      x2e(numerocourant) = x4b(i+1,j)
      z2e(numerocourant) = z4b(i+1,j)
      x3e(numerocourant) = x4b(i+1,j+1)
      z3e(numerocourant) = z4b(i+1,j+1)
      x4e(numerocourant) = x4b(i,j+1)
      z4e(numerocourant) = z4b(i,j+1)
      numerocourant = numerocourant + 1
    enddo
  enddo

  print *,'Number of elements stored = ',numerocourant - 1

! ***
! *** creation des elements sous forme topologique
! ***

  write(*,*)
  write(*,*) 'Creation de la topologie des elements...'

  do i=1,nelem

! recherche point 1
      do j=1,npoin
        dist = sqrt((x1e(i)-xpoint(j))**2 + (z1e(i)-zpoint(j))**2)
        if (dist <= rseuil) then
            nump1 = j
            goto 401
        endif
      enddo
      stop 'point not found !'
 401 continue

! recherche point 2
      do j=1,npoin
        dist = sqrt((x2e(i)-xpoint(j))**2 + (z2e(i)-zpoint(j))**2)
        if (dist <= rseuil) then
            nump2 = j
            goto 402
        endif
      enddo
      stop 'point not found !'
 402 continue

! recherche point 3
      do j=1,npoin
        dist = sqrt((x3e(i)-xpoint(j))**2 + (z3e(i)-zpoint(j))**2)
        if (dist <= rseuil) then
            nump3 = j
            goto 403
        endif
      enddo
      stop 'point not found !'
 403 continue

! recherche point 4
      do j=1,npoin
        dist = sqrt((x4e(i)-xpoint(j))**2 + (z4e(i)-zpoint(j))**2)
        if (dist <= rseuil) then
            nump4 = j
            goto 404
        endif
      enddo
      stop 'point not found !'
 404 continue

      numpoin1(i) = nump1
      numpoin2(i) = nump2
      numpoin3(i) = nump3
      numpoin4(i) = nump4

  enddo

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! *** generation de la base de donnees

  write(*,*)
  write(*,*) 'Generation de la base de donnees...'

  open(unit=15,file='OUTPUT_FILES/Database',status='unknown')

  title = 'Modele Canyon Paco'
  write(15,*) '#'
  write(15,*) '# Base de Donnees pour Specfem - Premailleur Fortran 90'
  write(15,*) '# ',title
  write(15,*) '# Dimitri Komatitsch, (c) EPS - Harvard February 1998'
  write(15,*) '#'

  write(15,*) 'Titre simulation'
  write(15,40) title

  ndofn = 2
  ndime = 2
  ngnod = 4
  nnode = 4
  nbcnd = 0
  nofst = 0
  npgeo = npoin
  nspel = nelem
  nbmodeles = 1
  nbsources = 1
  nrec = 150
  lquad = 1
  iexternal = .false.
  aleatoire = .false.
  topoplane = .false.
  simulate = .false.

  absorbhaut = .false.
  absorbbas = .false.
  absorbgauche = .false.
  absorbdroite = .false.
  absorbstacey = .true.
  absorbmodar = .false.
  ifullabs = .false.

  sismos = .true.
  isamp = 20
  nrec1 = nrec
  nrec2 = 0
  anglerec = 0.
  anglerec2 = 0.
  irepr = 1
  nrecsur3 = nrec / 3

  nt = 20000
  dt = 0.625e-3
  niter = 1
  alphanewm = 0.
  betanewm = 0.
  gammanewm = 0.5
  display = .true.
  ignuplot = .false.
  ivectplot = .true.
  icolorplot = .false.
  imeshvect = .true.
  imeshcolor = .false.
  imodelvect = .false.
  iboundvect = .false.
  interpol = .true.
  isymbols = .true.

!! DK DK Mexico August 1999, temporarily suppress external field
  initialfield = .true.
  initialfield = .false.

  itaff = 2000
  itfirstaff = 5
  cutvect = 1.
  cutcolor = 2.2
  scalex = 1.
  scalez = 1.
  sizemax = 1.
  pointsdisp = 7
  isubsamp = 2
  orig_x = 2.3
  orig_z = 3.4
  factorana = 50000.
  factorxsu = 3.5
  n1ana = 1
  n2ana = nrec1

  write(15,*) 'ndofn ndime npgeo'
  write(15,*) ndofn,ndime,npgeo

  write(15,*) 'display ignuplot interpol'
  write(15,*) display,ignuplot,interpol

  write(15,*) 'itaff itfirstaff icolor inumber'
  write(15,*) itaff,itfirstaff,0,0

  write(15,*) 'ivectplot imeshvect imodelvect iboundvect cutvect isubsamp'
  write(15,*) ivectplot,imeshvect,imodelvect,iboundvect,cutvect,isubsamp

  usletter = .true.
  write(15,*) 'scalex scalez sizemax angle rapport USletter'
  write(15,*) scalex,scalez,sizemax,20.,0.40,usletter

  write(15,*) 'orig_x orig_z isymbols'
  write(15,*) orig_x,orig_z,isymbols

  valseuil = 5.00
  freqmaxrep = 100.
  write(15,*) 'valseuil freqmaxrep'
  write(15,*) valseuil,freqmaxrep

  write(15,*) 'sismos nrec nrec1 nrec2 isamp'
  write(15,*) sismos,nrec,nrec1,nrec2,isamp

  write(15,*) 'irepr anglerec anglerec2'
  write(15,*) irepr,anglerec,anglerec2

  compenergy = .false.
  absstacey = .true.
  write(15,*) 'topoplane absstacey compenergy'
  write(15,*) topoplane,absstacey,compenergy

  write(15,*) 'initialfield factorana factorxsu n1ana n2ana'
  write(15,*) initialfield,factorana,factorxsu,n1ana,n2ana

  write(15,*) 'isismostype ivecttype iaffinfo'
  write(15,*) '1,  1,  40'
  write(15,*) 'ireadmodel ioutputgrid iavs ivisual3'
  write(15,*) 'F,  F,  F,  F'

  write(15,*) 'iexec iecho'
  write(15,*) '1       1'

  write(15,*) 'ncycl dtinc niter'
  write(15,*) nt,dt,niter

  write(15,*) 'alpha beta gamma'
  write(15,*) alphanewm,betanewm,gammanewm

  nbsources = 1
  write(15,*) 'nltfl (number of force or pressure sources)'
  write(15,*) nbsources

  itimetype = 6
  isource_type = 2
  f0 = 2.
  t0 = 0.55
  xs = +1.
  zs = -2.
  angle = 0.
  factor = 1.
  xoffs = 12.
  zoffs = 9.
  write(15,*) 'Collocated forces and/or pressure sources:'
  write(15,*) itimetype,isource_type,xs+xoffs,zs+zoffs,f0,t0,factor,angle,0

  write(15,*) 'Receivers (number, angle, position in meters)'
  do irec=1,nrec
  if (irec <= nrecsur3) then
    xrec = 2.*dble(irec-1)/dble(nrecsur3-1) + 9.
    zrec = 9.
  else if (irec >= 2*nrecsur3) then
    xrec = 2.*dble(irec-2*nrecsur3)/dble(nrecsur3) + 13.
    zrec = 9.
  else
    angle = pi + pi*dble(irec-nrecsur3)/dble(nrecsur3)
    xrec = 12. + cos(angle)
    zrec = 9. + sin(angle)
  endif
  write(15,*) irec,xrec,zrec
  enddo

  write(15,*) 'Coordinates of spectral control points'
  do i=1,npoin
    write(15,*) i,xpoint(i)+xoffs,zpoint(i)+zoffs
  enddo

  netyp = 2
  nxgll = 6
  nelemperio = 0
  nelemabs = 0

  write(15,*) 'params spectraux'
  write(15,*) netyp,nbmodeles,ngnod,nxgll,nxgll,nspel,pointsdisp,nelemabs,nelemperio

  write(15,*) 'Material sets (num 0 rho vp vs 0 0)'
  rho = 1.
  cp = 2.
  cs = 1.
  write(15,*) nbmodeles,0,rho,cp,cs,0,0

  write(15,*) 'Spectral-element topology'

  imatnum = 1

  do i=1,nspel
    write(15,*) i,imatnum,numpoin1(i),numpoin2(i),numpoin3(i),numpoin4(i)
  enddo

  close(15)

 40 format(a50)

  end program circular_canyon

