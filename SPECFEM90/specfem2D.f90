
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.0
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
!  An explicit 2D spectral element solver for the elastic wave equation
!
!========================================================================

! cleaned version 5.0 is based on SPECFEM2D version 4.2
! (c) June 1998 by Dimitri Komatitsch, Harvard University, USA
! and Jean-Pierre Vilotte, Institut de Physique du Globe de Paris, France

! version 5.0 : got rid of useless routines, suppressed commons etc.

  program specfem2D

  implicit none

  include "constants.h"

  character(len=80) datlin

  double precision gltfu(20)

  double precision, dimension(:,:), allocatable :: coorg,posrec
  double precision, dimension(:), allocatable :: coorgread
  double precision, dimension(:), allocatable :: posrecread

  double precision, dimension(:,:), allocatable :: sisux,sisuz

  logical anyabs

  integer i,j,it,irec,iglobrec,ipoin,ip,id,ip1,ip2,in,nnum
  integer nbpoin,inump,n,npoinext,netyp,ispec,npoin,npgeo

  double precision valux,valuz,rhoextread,vpextread,vsextread
  double precision dcosrot,dsinrot,xcor,zcor

! coefficients of the explicit Newmark time scheme
  integer NSTEP
  double precision deltatover2,deltatsqover2,time,deltat

  double precision, dimension(:), allocatable :: xigll,yigll,wxgll,wygll,xirec,etarec

  double precision, dimension(:,:), allocatable :: coord,accel,veloc,displ, &
    hprime,hTprime,flagrange,xinterp,zinterp,Uxinterp,Uzinterp,elastcoef

  double precision, dimension(:), allocatable :: rmass, &
    fglobx,fglobz,density,vpext,vsext,rhoext,displread,velocread,accelread

  double precision, dimension(:,:,:), allocatable :: shapeint,shape,dvolu, &
    a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a13x,a13z,Uxnewloc,Uznewloc

  double precision, dimension(:,:), allocatable :: a11,a12

  double precision, dimension(:,:,:,:), allocatable :: dershape

  double precision, dimension(:,:,:,:,:), allocatable :: xjaci

  integer, dimension(:,:,:), allocatable :: ibool
  integer, dimension(:,:), allocatable  :: knods,codeabs
  integer, dimension(:), allocatable :: kmato,numabs,is_bordabs

  integer ie,k
  integer inum,itourne,ntourne,idummy,numabsread
  integer isource,iexplo
  integer codeabsread(4)

  double precision rsizemin,rsizemax,cpoverdxmin,cpoverdxmax, &
    rlamdaSmin,rlamdaSmax,rlamdaPmin,rlamdaPmax,vpmin,vpmax

  integer icolor,inumber,isubsamp,ivecttype,itaff,nrec,isismostype
  integer numat,ngnod,nspec,iptsdisp,nelemabs

  logical interpol,imeshvect,imodelvect,iboundvect,ireadmodel,initialfield, &
    ioutputgrid,ignuplot

  double precision cutvect,anglerec

! title of the plot
  character(len=60) stitle

!
!***********************************************************************
!
!             i n i t i a l i z a t i o n    p h a s e
!
!***********************************************************************

  open (IIN,file='DataBase')

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
  write(*,*) '*********************************'
  write(*,*) '****                         ****'
  write(*,*) '****  SPECFEM2D VERSION 5.0  ****'
  write(*,*) '****                         ****'
  write(*,*) '*********************************'

!
!---- read parameters from input file
!

  read(IIN,40) datlin
  read(IIN,*) npgeo

  read(IIN,40) datlin
  read(IIN,*) ignuplot,interpol

  read(IIN,40) datlin
  read(IIN,*) itaff,icolor,inumber

  read(IIN,40) datlin
  read(IIN,*) imeshvect,imodelvect,iboundvect,cutvect,isubsamp
  cutvect = cutvect / 100.d0

  read(IIN,40) datlin
  read(IIN,*) nrec,anglerec

  read(IIN,40) datlin
  read(IIN,*) initialfield

  read(IIN,40) datlin
  read(IIN,*) isismostype,ivecttype

  read(IIN,40) datlin
  read(IIN,*) ireadmodel,ioutputgrid

!---- check parameters read
  write(IOUT,200) npgeo,NDIME
  write(IOUT,600) itaff,icolor,inumber
  write(IOUT,700) nrec,isismostype,anglerec
  write(IOUT,750) initialfield,ireadmodel,ioutputgrid
  write(IOUT,800) ivecttype,100.d0*cutvect,isubsamp

!---- read time step
  read(IIN,40) datlin
  read(IIN,*) NSTEP,deltat
  write(IOUT,703) NSTEP,deltat,NSTEP*deltat

!
!---- allocate first arrays needed
!
  if(nrec < 1) stop 'need at least one receiver'
  allocate(sisux(NSTEP,nrec))
  allocate(sisuz(NSTEP,nrec))
  allocate(posrec(NDIME,nrec))
  allocate(coorg(NDIME,npgeo))

!
!----  read source information
!
  read(IIN,40) datlin
  read(IIN,*) (gltfu(k), k=1,9)

!
!-----  check the input
!
 if(.not. initialfield) then
   if(nint(gltfu(1)) /= 6) stop 'Wrong function number in getltf !'
   isource = nint(gltfu(1))
   iexplo = nint(gltfu(2))
   if (iexplo == 1) then
     write(IOUT,212) (gltfu(k), k=3,8)
   else if(iexplo == 2) then
     write(IOUT,222) (gltfu(k), k=3,7)
   else
     stop 'Unknown source type number !'
   endif
 endif

!
!-----  convert angle from degrees to radians
!
   isource = nint(gltfu(1))
   iexplo = nint(gltfu(2))
   if(isource >= 4.and.isource <= 6.and.iexplo == 1) gltfu(8) = gltfu(8) * pi / 180.d0

!
!---- read receiver locations
!
  irec = 0
  read(IIN,40) datlin
  allocate(posrecread(NDIME))
  do i=1,nrec
   read(IIN ,*) irec,(posrecread(j),j=1,NDIME)
   if(irec<1 .or. irec>nrec) stop 'Wrong receiver number'
   posrec(:,irec) = posrecread
  enddo
  deallocate(posrecread)

!
!---- read the spectral macrobloc nodal coordinates
!
  ipoin = 0
  read(IIN,40) datlin
  allocate(coorgread(NDIME))
  do ip = 1,npgeo
   read(IIN,*) ipoin,(coorgread(id),id =1,NDIME)
   if(ipoin<1 .or. ipoin>npgeo) stop 'Wrong control point number'
   coorg(:,ipoin) = coorgread
  enddo
  deallocate(coorgread)

!
!---- read the basic properties of the spectral elements
!
  read(IIN ,40) datlin
  read(IIN ,*) netyp,numat,ngnod,nspec,iptsdisp,nelemabs

!
!---- allocate arrays
!

allocate(shape(ngnod,NGLLX,NGLLY))
allocate(shapeint(ngnod,iptsdisp,iptsdisp))
allocate(dershape(NDIME,ngnod,NGLLX,NGLLY))
allocate(dvolu(nspec,NGLLX,NGLLY))
allocate(xjaci(nspec,NDIME,NDIME,NGLLX,NGLLY))
allocate(hprime(NGLLX,NGLLY))
allocate(hTprime(NGLLX,NGLLY))
allocate(a1(NGLLX,NGLLY,nspec))
allocate(a2(NGLLX,NGLLY,nspec))
allocate(a3(NGLLX,NGLLY,nspec))
allocate(a4(NGLLX,NGLLY,nspec))
allocate(a5(NGLLX,NGLLY,nspec))
allocate(a6(NGLLX,NGLLY,nspec))
allocate(a7(NGLLX,NGLLY,nspec))
allocate(a8(NGLLX,NGLLY,nspec))
allocate(a9(NGLLX,NGLLY,nspec))
allocate(a10(NGLLX,NGLLY,nspec))
allocate(a11(NGLLX,NGLLY))
allocate(a12(NGLLX,NGLLY))
allocate(xigll(NGLLX))
allocate(yigll(NGLLY))
allocate(wxgll(NGLLX))
allocate(wygll(NGLLY))
allocate(Uxnewloc(NGLLX,NGLLY,nspec))
allocate(Uznewloc(NGLLX,NGLLY,nspec))
allocate(xirec(iptsdisp))
allocate(etarec(iptsdisp))
allocate(flagrange(NGLLX,iptsdisp))
allocate(xinterp(iptsdisp,iptsdisp))
allocate(zinterp(iptsdisp,iptsdisp))
allocate(Uxinterp(iptsdisp,iptsdisp))
allocate(Uzinterp(iptsdisp,iptsdisp))
allocate(density(numat))
allocate(elastcoef(4,numat))
allocate(kmato(nspec))
allocate(knods(ngnod,nspec))
allocate(ibool(NGLLX,NGLLY,nspec))

! --- allocate arrays for absorbing boundary conditions
  if(nelemabs <= 0) then
    nelemabs = 1
    anyabs = .false.
  else
    anyabs = .true.
  endif
  allocate(is_bordabs(nspec))
  allocate(numabs(nelemabs))
  allocate(codeabs(4,nelemabs))

!
!---- print element group main parameters
!
  write(IOUT,107)
  write(IOUT,207) nspec,ngnod,NGLLX,NGLLY,NGLLX*NGLLY,iptsdisp,numat,nelemabs

!
!----    set up coordinates of the Gauss-Lobatto-Legendre points
!
 call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
 call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)

!
!---- if nb of points is odd, the middle abscissa is exactly zero
!
  if(mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = ZERO
  if(mod(NGLLY,2) /= 0) yigll((NGLLY-1)/2+1) = ZERO

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
  read(IIN ,40) datlin
  do n=1,nelemabs
    read(IIN ,*) inum,numabsread,codeabsread(1), &
           codeabsread(2),codeabsread(3),codeabsread(4)
    if(inum < 1 .or. inum > nelemabs) stop 'Wrong absorbing element number'
    numabs(inum) = numabsread
    codeabs(ihaut,inum)   = codeabsread(1)
    codeabs(ibas,inum)    = codeabsread(2)
    codeabs(igauche,inum) = codeabsread(3)
    codeabs(idroite,inum) = codeabsread(4)

!----  eventuellement tourner element counterclockwise si condition absorbante

       if(codeabs(ibas,inum) == iaretebas .or. &
          codeabs(ihaut,inum) == iaretehaut .or. &
          codeabs(igauche,inum) == iaretegauche .or. &
          codeabs(idroite,inum) == iaretedroite) then
            ntourne = 0

  else if(codeabs(ibas,inum) == iaretegauche .or. &
          codeabs(ihaut,inum) == iaretedroite .or. &
          codeabs(igauche,inum) == iaretehaut .or. &
          codeabs(idroite,inum) == iaretebas) then
            ntourne = 3

  else if(codeabs(ibas,inum) == iaretehaut .or. &
          codeabs(ihaut,inum)  == iaretebas .or. &
          codeabs(igauche,inum) == iaretedroite .or. &
          codeabs(idroite,inum) == iaretegauche) then
            ntourne = 2

  else if(codeabs(ibas,inum) == iaretedroite .or. &
          codeabs(ihaut,inum) == iaretegauche .or. &
          codeabs(igauche,inum) == iaretebas .or. &
          codeabs(idroite,inum) == iaretehaut) then
            ntourne = 1
  else
     stop 'Error in absorbing conditions numbering'
  endif

!----  rotate element counterclockwise
  if(ntourne /= 0) then

  do itourne = 1,ntourne

      idummy = knods(1,numabs(inum))
      knods(1,numabs(inum)) = knods(2,numabs(inum))
      knods(2,numabs(inum)) = knods(3,numabs(inum))
      knods(3,numabs(inum)) = knods(4,numabs(inum))
      knods(4,numabs(inum)) = idummy

    if(ngnod == 9) then
      idummy = knods(5,numabs(inum))
      knods(5,numabs(inum)) = knods(6,numabs(inum))
      knods(6,numabs(inum)) = knods(7,numabs(inum))
      knods(7,numabs(inum)) = knods(8,numabs(inum))
      knods(8,numabs(inum)) = idummy
    endif

  enddo

  endif

  enddo
  write(*,*)
  write(*,*) 'Number of absorbing elements: ',nelemabs
  endif


!
!---- compute the spectral element shape functions and their local derivatives
!
  call q49shape(shape,dershape,xigll,yigll,ngnod,NGLLX,NGLLY,NDIME)

!
!---- generate the global numbering
!

! version "propre mais lente" ou version "sale mais rapide"
  if(fast_numbering) then
    call createnum_fast(knods,ibool,shape,coorg,npoin,npgeo,nspec,ngnod)
  else
    call createnum_slow(knods,ibool,npoin,nspec,ngnod)
  endif

!
!---- compute the spectral element jacobian matrix
!

  call q49spec(shapeint,dershape,dvolu,xjaci,xigll,coorg,knods,ngnod, &
          NGLLX,NGLLY,NDIME,nspec,npgeo,xirec,etarec,flagrange,iptsdisp)

!
!---- close input file
!
  close(IIN)

!
!----  allocation des autres tableaux pour la grille globale et les bords
!

  allocate(coord(NDIME,npoin))
  allocate(accel(NDIME,npoin))
  allocate(displ(NDIME,npoin))
  allocate(veloc(NDIME,npoin))
  allocate(rmass(npoin))
  allocate(fglobx(npoin))
  allocate(fglobz(npoin))

  if(ireadmodel) then
    npoinext = npoin
  else
    npoinext = 1
  endif
  allocate(vpext(npoinext))
  allocate(vsext(npoinext))
  allocate(rhoext(npoinext))

  allocate(a13x(NGLLX,NGLLY,nelemabs))
  allocate(a13z(NGLLX,NGLLY,nelemabs))

!
!----  set the coordinates of the points of the global grid
!
  do ispec = 1,nspec
  do ip1 = 1,NGLLX
  do ip2 = 1,NGLLY

      xcor = zero
      zcor = zero
      do in = 1,ngnod
        nnum = knods(in,ispec)
        xcor = xcor + shape(in,ip1,ip2)*coorg(1,nnum)
        zcor = zcor + shape(in,ip1,ip2)*coorg(2,nnum)
      enddo

      coord(1,ibool(ip1,ip2,ispec)) = xcor
      coord(2,ibool(ip1,ip2,ispec)) = zcor

   enddo
   enddo
   enddo

!
!--- save the grid of points in a file
!
  if(ioutputgrid) then
    print *
    print *,'Saving the grid in a text file...'
    print *
    open(unit=55,file='gridpoints.txt',status='unknown')
    write(55,*) npoin
    do n = 1,npoin
      write(55,*) n,(coord(i,n), i=1,NDIME)
    enddo
    close(55)
  endif

!
!-----   plot the GLL mesh in a Gnuplot file
!
  if(ignuplot) call plotgll(knods,ibool,coorg,coord,npoin,npgeo,ngnod,nspec)

!
!----   define coefficients of the Newmark time scheme
!
  deltatover2 = 0.5d0*deltat
  deltatsqover2 = deltat*deltat/2.d0

!
!---- definir la position reelle des points source et recepteurs
!
  call positsource(coord,ibool,gltfu,npoin,nspec)
  call positrec(coord,posrec,npoin,nrec)

!
!----  eventuellement lecture d'un modele externe de vitesse et de densite
!
  if(ireadmodel) then
    print *
    print *,'Reading velocity and density model from external file...'
    print *
    open(unit=55,file='extmodel.txt',status='unknown')
    read(55,*) nbpoin
    if(nbpoin /= npoin) stop 'Wrong number of points in input file'
    do n = 1,npoin
      read(55,*) inump,rhoextread,vpextread,vsextread
      if(inump<1 .or. inump>npoin) stop 'Wrong point number'
      rhoext(inump) = rhoextread
      vpext(inump) = vpextread
      vsext(inump) = vsextread
    enddo
    close(55)
  endif

!
!----   build the mass matrix for spectral elements
!
  call qmasspec(rhoext,wxgll,wygll,ibool,dvolu,rmass,density,kmato,npoin,ireadmodel,nspec,numat)

!
!---- definir les tableaux a1 a a13
!
  call defarrays(vpext,vsext,rhoext,density,elastcoef, &
          xigll,yigll,wxgll,wygll,hprime,hTprime, &
          a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13x,a13z, &
          ibool,kmato,dvolu,xjaci,coord,gltfu, &
          numabs,codeabs,anyabs,npoin,rsizemin,rsizemax, &
          cpoverdxmin,cpoverdxmax,rlamdaSmin,rlamdaSmax, &
          rlamdaPmin,rlamdaPmax,vpmin,vpmax,ireadmodel,nelemabs,nspec,numat)

! initialiser les tableaux a zero
  accel = zero
  veloc = zero
  displ = zero

!
!--- precalculer l'inverse de la matrice de masse pour efficacite
!
  rmass(:) = one / rmass(:)

! calculer la numerotation inverse pour les bords absorbants
  is_bordabs(:) = 0
  if(anyabs) then
    do ispec = 1,nelemabs
      is_bordabs(numabs(ispec)) = ispec
    enddo
  endif

! convertir angle recepteurs en radians
  anglerec = anglerec * pi / 180.d0

!
!----  eventuellement lecture des champs initiaux dans un fichier
!
  if(initialfield) then
      print *
      print *,'Reading initial fields from external file...'
      print *
      open(unit=55,file='wavefields.txt',status='unknown')
      read(55,*) nbpoin
      if(nbpoin /= npoin) stop 'Wrong number of points in input file'
      allocate(displread(NDIME))
      allocate(velocread(NDIME))
      allocate(accelread(NDIME))
      do n = 1,npoin
        read(55,*) inump, (displread(i), i=1,NDIME), &
            (velocread(i), i=1,NDIME), (accelread(i), i=1,NDIME)
        if(inump<1 .or. inump>npoin) stop 'Wrong point number'
        displ(:,inump) = displread
        veloc(:,inump) = velocread
        accel(:,inump) = accelread
      enddo
      deallocate(displread)
      deallocate(velocread)
      deallocate(accelread)
      close(55)
  endif

!
!---- afficher le max du deplacement initial
!
  print *,'Max norme U initial = ',maxval(sqrt(displ(1,:)**2 + displ(2,:)**2))

!
!---- verifier le maillage, la stabilite et le nb de points par lambda
!
  call checkgrid(deltat,gltfu,initialfield,rsizemin,rsizemax, &
      cpoverdxmin,cpoverdxmax,rlamdaSmin,rlamdaSmax,rlamdaPmin,rlamdaPmax)

!
!---- initialiser sismogrammes
!
  sisux = zero
  sisuz = zero

  dcosrot = dcos(anglerec)
  dsinrot = dsin(anglerec)

!
!----          s t a r t   t i m e   i t e r a t i o n s
!

  write(IOUT,400)

! boucle principale d'evolution en temps
  do it=1,NSTEP

! compute current time
  time = (it-1)*deltat

  if(mod(it-1,itaff)  ==  0) then
    if(time  >=  1.d-3) then
      write(IOUT,100) it,time
    else
      write(IOUT,101) it,time
    endif
  endif

! calculer le predictor
  displ(:,:) = displ(:,:) + deltat*veloc(:,:) + deltatsqover2*accel(:,:)
  veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)
  accel(:,:) = zero

!
!----  calcul du residu d'acceleration pour le corrector
!----  retourne dans accel le terme Fext - M*A(i,n+1) - K*D(i,n+1)
!
  call qsumspec(hprime,hTprime,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13x,a13z, &
          ibool,displ,veloc,accel,Uxnewloc,Uznewloc,rmass,npoin, &
          nspec,gltfu,initialfield,is_bordabs,nelemabs,anyabs,time)

!
!----  mise a jour globale du deplacement par corrector
!
  veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)

!
!----  display max of norm of displacement
!
  if(mod(it-1,itaff)  ==  0) &
     print *,'Max norme U = ',maxval(sqrt(displ(1,:)**2 + displ(2,:)**2))

! store the seismograms
  do irec=1,nrec
    iglobrec = nint(posrec(1,irec))

  if(isismostype  ==  1) then
    valux = displ(1,iglobrec)
    valuz = displ(2,iglobrec)
  else if(isismostype  ==  2) then
    valux = veloc(1,iglobrec)
    valuz = veloc(2,iglobrec)
  else if(isismostype  ==  3) then
    valux = accel(1,iglobrec)
    valuz = accel(2,iglobrec)
  else
    stop 'Wrong field code for seismogram output'
  endif

! rotation eventuelle des composantes
  sisux(it,irec) =   dcosrot*valux + dsinrot*valuz
  sisuz(it,irec) = - dsinrot*valux + dcosrot*valuz

  enddo

!
!----  affichage des resultats a certains pas de temps
!
  if(mod(it,itaff) == 0 .or. it == 5 .or. it == NSTEP) then

  write(IOUT,*)
  if(time  >=  1.d-3) then
    write(IOUT,110) time
  else
    write(IOUT,111) time
  endif
  write(IOUT,*)

!
!----  affichage postscript
!
  write(IOUT,*) 'Dump PostScript'
  if(ivecttype  ==  1) then
    write(IOUT,*) 'drawing displacement field...'
    call plotpost(displ,coord,vpext,gltfu,posrec, &
          it,deltat,coorg,xinterp,zinterp,shapeint, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,stitle,npoin,npgeo,vpmin,vpmax,nrec, &
          icolor,inumber,isubsamp,ivecttype,interpol,imeshvect,imodelvect, &
          iboundvect,ireadmodel,cutvect,nelemabs,numat,iptsdisp,nspec,ngnod)
  else if(ivecttype  ==  2) then
    write(IOUT,*) 'drawing velocity field...'
    call plotpost(veloc,coord,vpext,gltfu,posrec, &
          it,deltat,coorg,xinterp,zinterp,shapeint, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,stitle,npoin,npgeo,vpmin,vpmax,nrec, &
          icolor,inumber,isubsamp,ivecttype,interpol,imeshvect,imodelvect, &
          iboundvect,ireadmodel,cutvect,nelemabs,numat,iptsdisp,nspec,ngnod)
  else if(ivecttype  ==  3) then
    write(IOUT,*) 'drawing acceleration field...'
    call plotpost(accel,coord,vpext,gltfu,posrec, &
          it,deltat,coorg,xinterp,zinterp,shapeint, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,stitle,npoin,npgeo,vpmin,vpmax,nrec, &
          icolor,inumber,isubsamp,ivecttype,interpol,imeshvect,imodelvect, &
          iboundvect,ireadmodel,cutvect,nelemabs,numat,iptsdisp,nspec,ngnod)
  else
    stop 'wrong field code for PostScript display'
  endif
  write(IOUT,*) 'Fin dump PostScript'

!----  save temporary seismograms
  call write_seismograms(sisux,sisuz,NSTEP,nrec,deltat)

  endif

  enddo ! end of the main time loop

!----  save final seismograms
  call write_seismograms(sisux,sisuz,NSTEP,nrec,deltat)

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

  200   format(//1x,'C o n t r o l',/1x,34('='),//5x,&
  'Number of spectral element control nodes. . (npgeo) =',i8/5x, &
  'Number of space dimensions . . . . . . . . . (NDIME) =',i8)
  600   format(//1x,'C o n t r o l',/1x,34('='),//5x, &
  'Display frequency  . . . . . . . . . . . . . (itaff) = ',i5/ 5x, &
  'Color display . . . . . . . . . . . . . . . (icolor) = ',i5/ 5x, &
  '        ==  0     black and white display              ',  / 5x, &
  '        ==  1     color display                        ',  /5x, &
  'Numbered mesh . . . . . . . . . . . . . . .(inumber) = ',i5/ 5x, &
  '        ==  0     do not number the mesh               ',  /5x, &
  '        ==  1     number the mesh                      ')
  700   format(//1x,'C o n t r o l',/1x,34('='),//5x, &
  'Total number of receivers. . . . . . . . . . .(nrec) = ',i6/5x, &
  'Seismograms recording type. . . . . . .(isismostype) = ',i6/5x, &
  'Angle for first line of receivers. . . . .(anglerec) = ',f6.2)
  750   format(//1x,'C o n t r o l',/1x,34('='),//5x, &
  'Read external initial field or not . .(initialfield) = ',l6/5x, &
  'Read external velocity model or not. . .(ireadmodel) = ',l6/5x, &
  'Save grid in external file or not . . .(ioutputgrid) = ',l6)
  800   format(//1x,'C o n t r o l',/1x,34('='),//5x, &
  'Vector display type . . . . . . . . . . .(ivecttype) = ',i6/5x, &
  'Percentage of cut for vector plots. . . . .(cutvect) = ',f6.2/5x, &
  'Subsampling for velocity model display . .(isubsamp) = ',i6)

  703   format(//' I t e r a t i o n s '/1x,29('='),//5x, &
      'Number of time iterations . . . . .(NSTEP) =',i8,/5x, &
      'Time step increment . . . . . . . .(deltat) =',1pe15.6,/5x, &
      'Total simulation duration . . . . . (ttot) =',1pe15.6)

  107   format(/5x,'--> Isoparametric Spectral Elements <--',//)
  207   format(5x, &
           'Number of spectral elements . . . . .  (nspec) =',i7,/5x, &
           'Number of control nodes per element .  (ngnod) =',i7,/5x, &
           'Number of points in X-direction . . .  (NGLLX) =',i7,/5x, &
           'Number of points in Y-direction . . .  (NGLLY) =',i7,/5x, &
           'Number of points per element. . .(NGLLX*NGLLY) =',i7,/5x, &
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
     'Source Type. . . . . . . . . . . . . . = Explosion',/5x, &
     'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
     'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
     'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
     'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
     'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x)

  end program specfem2D

