
!=====================================================================
!
!                          S p e c f e m
!                          -------------
!
!                           Version 4.2
!                           -----------
!
!                         Dimitri Komatitsch
!    Department of Earth and Planetary Sciences - Harvard University
!
!                         Jean-Pierre Vilotte
!                 Departement de Sismologie - IPGP - Paris
!
!                           (c) June 1998
!
!=====================================================================
!
!            An explicit spectral element solver for the
!
!                      elastic wave equation
!
!=======================================================================

  program main
!
!=======================================================================
!
!  "m a i n" :  Allocate memory, initialize arrays and iterate in time
!   -------
!
! ======================================================================
!
  use iounit
  use captio
  use infos
  use mesh01
  use constspec
  use timeparams
  use defpi
  use spela202
  use energie
  use arraydir
  use loadft

  implicit none

  double precision, parameter :: zero = 0.d0, one = 1.d0

  character(len=80) datlin

  double precision, dimension(:,:), allocatable :: gltfu,force,coorg,posrec

! simple precision pour stockage sismogrammes au format SEP
  real, dimension(:,:), allocatable :: sisux,sisuz

  logical anyabs,anyperio

  integer i,it,irec,iter,itsis,iglobrec,iglobsource
  integer nbpoin,inump,n,npoinext,nseis,netyp,ipoin,ispec

  double precision valux,valuz,rhoextread,vpextread,vsextread
  double precision dcosrot,dsinrot,dcosrot1,dsinrot1,dcosrot2,dsinrot2

! coefficients of the explicit Newmark time scheme
  double precision deltatover2,deltatsqover2

!
!---- tableaux pour allocation dynamique
!

  double precision, dimension(:), allocatable :: &
    xi,yi,wx,wy,xirec,etarec

  double precision, dimension(:,:), allocatable :: &
    hprime,hTprime,flagrange,xinterp,zinterp,Uxinterp,Uzinterp,elastcoef, &
    coord,accel,veloc,displ,vpred

  double precision, dimension(:), allocatable :: rmass, &
    fglobx,fglobz,density,vpext,vsext,rhoext,displread,velocread,accelread

  double precision, dimension(:,:,:), allocatable :: shapeint,shape,dvolu, &
    a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13x,a13z,Uxnewloc,Uznewloc

  double precision, dimension(:,:,:,:), allocatable :: dershape

  double precision, dimension(:,:,:,:,:), allocatable :: xjaci

  integer, dimension(:,:,:), allocatable :: ibool,iboolori
  integer, dimension(:,:), allocatable  :: knods,codeabs,codeperio
  integer, dimension(:), allocatable :: kmato,numabs,is_bordabs

!
!***********************************************************************
!
!             i n i t i a l i z a t i o n    p h a s e
!
!***********************************************************************
!
!---- Assign unit numbers
!
  iin  = 8
  open (iin,file='DataBase')

  iout = 6
!  ecriture dans un fichier texte et non ecran
!  iout = 16
!  open (iout,file='results.txt')

! fichier pour le stockage des courbes d'energie
  ienergy = 17

!
!---  read job title and skip remaining titles of the input file
!
  read(iin ,40) datlin
  read(iin ,40) datlin
  read(iin ,40) jtitle
  read(iin ,40) datlin
  read(iin ,40) datlin
  read(iin ,40) datlin
  read(iin ,45) stitle
!
!---- Print the date, time and start-up banner
!
  call datim('Program  S P E C F E M : start',stitle,iout)

  write(*,*)
  write(*,*)
  write(*,*) '******************************************'
  write(*,*) '****                                  ****'
  write(*,*) '****  SPECFEM VERSION 4.2 FORTRAN 90  ****'
  write(*,*) '****                                  ****'
  write(*,*) '******************************************'

!
!***********************************************************************
!
!                      i n p u t   p h a s e
!
!***********************************************************************
!

!
!---- read first control parameters
!
  call contol
!
!---- read iteration parameters
!
  call intseq
!
!---- allocate first arrays needed
!

! mettre a zero la structure de stockage des tableaux
  nbarrays = 0
  arraysizes(:) = 0
  arraynames(:) = '                    '

  if(sismos) then
    nseis = ncycl/isamp
  else
    nseis = 1
  endif

  allocate(sisux(nseis,nrec))
  allocate(sisuz(nseis,nrec))
  allocate(posrec(ndime,max(nrec,1)))
  allocate(coorg(ndime,npgeo))
  allocate(force(ndime,max(nltfl,1)))
  allocate(gltfu(20,max(nltfl,1)))

  call storearray('sisux',nseis*nrec,isngl)
  call storearray('sisuz',nseis*nrec,isngl)
  call storearray('posrec',ndime*max(nrec,1),idouble)
  call storearray('coorg',ndime*npgeo,idouble)
  call storearray('force',ndime*max(nltfl,1),idouble)
  call storearray('gltfu',20*max(nltfl,1),idouble)

!-----------------------------------------------------------------------

!
!---- read load time functions
!

!
!----    Collocated forces or pressure sources
!
  if(nltfl > 0) call getltf(gltfu,nltfl,initialfield)

!
!----  lecture position receivers
!
  if(nrec > 0) call getrecepts(posrec,ndime,nrec)

!
!---- read the spectral macroblocs nodal coordinates
!
  call getspec(coorg,npgeo,ndime)

!
!***********************************************************************
!
!       S p e c t r a l    E l e m e n t s    P a r a m e t e r s
!
!***********************************************************************
!

!
!---- read the basic properties of the spectral elements
!
  read(iin ,40) datlin
  read(iin ,*) netyp,numat,ngnod,nxgll,nygll,nspec,iptsdisp,nelemabs,nelemperio

!
!---- check that the mesh is conform
!
  if(nxgll /= nygll) stop 'Non conform mesh in input'
!
!***********************************************************************
!
!              A l l o c a t e   a r r a y s
!
!***********************************************************************
!

allocate(shape(ngnod,nxgll,nygll))
allocate(shapeint(ngnod,iptsdisp,iptsdisp))
allocate(dershape(ndime,ngnod,nxgll,nygll))
allocate(dvolu(nspec,nxgll,nygll))
allocate(xjaci(nspec,ndime,ndime,nxgll,nygll))
allocate(hprime(nxgll,nygll))
allocate(hTprime(nxgll,nygll))
allocate(a1(nxgll,nygll,nspec))
allocate(a2(nxgll,nygll,nspec))
allocate(a3(nxgll,nygll,nspec))
allocate(a4(nxgll,nygll,nspec))
allocate(a5(nxgll,nygll,nspec))
allocate(a6(nxgll,nygll,nspec))
allocate(a7(nxgll,nygll,nspec))
allocate(a8(nxgll,nygll,nspec))
allocate(a9(nxgll,nygll,nspec))
allocate(a10(nxgll,nygll,nspec))
allocate(a11(nxgll,nygll,max(nltfl,1)))
allocate(a12(nxgll,nygll,max(nltfl,1)))
allocate(xi(nxgll))
allocate(yi(nygll))
allocate(wx(nxgll))
allocate(wy(nygll))
allocate(Uxnewloc(nxgll,nygll,nspec))
allocate(Uznewloc(nxgll,nygll,nspec))
allocate(xirec(iptsdisp))
allocate(etarec(iptsdisp))
allocate(flagrange(nxgll,iptsdisp))
allocate(xinterp(iptsdisp,iptsdisp))
allocate(zinterp(iptsdisp,iptsdisp))
allocate(Uxinterp(iptsdisp,iptsdisp))
allocate(Uzinterp(iptsdisp,iptsdisp))
allocate(density(numat))
allocate(elastcoef(4,numat))

allocate(kmato(nspec))
allocate(knods(ngnod,nspec))
allocate(ibool(nxgll,nygll,nspec))

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  call storearray('shape',ngnod*nxgll*nygll,idouble)
  call storearray('shapeint',ngnod*iptsdisp*iptsdisp,idouble)
  call storearray('dershape',ndime*ngnod*nxgll*nygll,idouble)
  call storearray('dvolu',nspec*nxgll*nygll,idouble)
  call storearray('xjaci',nspec*ndime*ndime*nxgll*nygll,idouble)
  call storearray('hprime',nxgll*nygll,idouble)
  call storearray('hTprime',nxgll*nygll,idouble)
  call storearray('a1',nxgll*nygll*nspec,idouble)
  call storearray('a2',nxgll*nygll*nspec,idouble)
  call storearray('a3',nxgll*nygll*nspec,idouble)
  call storearray('a4',nxgll*nygll*nspec,idouble)
  call storearray('a5',nxgll*nygll*nspec,idouble)
  call storearray('a6',nxgll*nygll*nspec,idouble)
  call storearray('a7',nxgll*nygll*nspec,idouble)
  call storearray('a8',nxgll*nygll*nspec,idouble)
  call storearray('a9',nxgll*nygll*nspec,idouble)
  call storearray('a10',nxgll*nygll*nspec,idouble)
  call storearray('a11',nxgll*nygll*max(nltfl,1),idouble)
  call storearray('a12',nxgll*nygll*max(nltfl,1),idouble)
  call storearray('xi',nxgll,idouble)
  call storearray('yi',nygll,idouble)
  call storearray('wx',nxgll,idouble)
  call storearray('wy',nygll,idouble)
  call storearray('Uxnewloc',nxgll*nygll*nspec,idouble)
  call storearray('Uznewloc',nxgll*nygll*nspec,idouble)
  call storearray('xirec',iptsdisp,idouble)
  call storearray('etarec',iptsdisp,idouble)
  call storearray('flagrange',nxgll*iptsdisp,idouble)
  call storearray('xinterp',iptsdisp*iptsdisp,idouble)
  call storearray('zinterp',iptsdisp*iptsdisp,idouble)
  call storearray('Uxinterp',iptsdisp*iptsdisp,idouble)
  call storearray('Uzinterp',iptsdisp*iptsdisp,idouble)
  call storearray('density',numat,idouble)
  call storearray('elastcoef',4*numat,idouble)

  call storearray('kmato',nspec,iinteg)
  call storearray('knods',ngnod*nspec,iinteg)
  call storearray('ibool',nxgll*nygll*nspec,iinteg)

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! --- allocate arrays for absorbing and periodic boundary conditions

  if(nelemabs <= 0) then
    nelemabs = 1
    anyabs = .false.
  else
    anyabs = .true.
  endif
  allocate(is_bordabs(nspec))
  allocate(numabs(nelemabs))
  allocate(codeabs(4,nelemabs))
  call storearray('is_bordabs',nspec,iinteg)
  call storearray('numabs',nelemabs,iinteg)
  call storearray('codeabs',4*nelemabs,iinteg)

  if(nelemperio <= 0) then
    nelemperio = 1
    anyperio = .false.
!!!!!!!!!!!!    allocate(iboolori(1,1,1))
!!! DK DK fix bug with Linux pgf90 compiler
    allocate(iboolori(nxgll,nygll,nspec))
    call storearray('iboolori',1,iinteg)
  else
    anyperio = .true.
    allocate(iboolori(nxgll,nygll,nspec))
    call storearray('iboolori',nxgll*nygll*nspec,iinteg)
  endif
  allocate(codeperio(4,nelemperio))
  call storearray('codeperio',4*nelemperio,iinteg)

!
!----  input element data and compute total number of points
!
  call qinpspec(density,elastcoef,xi,yi,wx,wy,knods, &
    ibool,kmato,shape,shapeint,dershape,dvolu,xjaci,coorg, &
    xirec,etarec,flagrange,numabs,codeabs,codeperio,anyabs,anyperio)

!
!---- close input file
!
  close(iin)

!
!----  allocation des autres tableaux pour la grille globale et les bords
!

  allocate(coord(ndime,npoin))
  allocate(accel(ndime,npoin))
  allocate(displ(ndime,npoin))
  allocate(veloc(ndime,npoin))
  allocate(vpred(ndime,npoin))
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

  allocate(a13x(nxgll,nygll,nelemabs))
  allocate(a13z(nxgll,nygll,nelemabs))

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  call storearray('coord',ndime*npoin,idouble)
  call storearray('accel',ndime*npoin,idouble)
  call storearray('displ',ndime*npoin,idouble)
  call storearray('veloc',ndime*npoin,idouble)
  call storearray('vpred',ndime*npoin,idouble)
  call storearray('rmass',npoin,idouble)
  call storearray('fglobx',npoin,idouble)
  call storearray('fglobz',npoin,idouble)
  call storearray('vpext',npoinext,idouble)
  call storearray('vsext',npoinext,idouble)
  call storearray('rhoext',npoinext,idouble)

  call storearray('a13x',nxgll*nygll*nelemabs,idouble)
  call storearray('a13z',nxgll*nygll*nelemabs,idouble)

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!
!----  list a short directory after input phase
!
  if(iecho /= 0) call dircty

!
!----  set the coordinates of the points of the global grid
!
  call setcor(coord,npoin,ndime,knods,shape,ibool,coorg,nxgll,nygll, &
                         nspec,npgeo,ngnod,ioutputgrid)

!
!-----   plot the GLL mesh in a Gnuplot file
!
  if (ignuplot) call plotgll(knods,ibool,coorg,coord)

!
!----   define coefficients of the Newmark time scheme
!
  deltatover2 = 0.5d0*deltat
  deltatsqover2 = deltat*deltat/2.d0

!
!----   mettre en oeuvre les periodic boundary conditions
!
  if(anyperio) call modifperio(ibool,iboolori,codeperio)

!
!---- definir la position reelle des points source et recepteurs
!
  call positsource(coord,ibool,gltfu,ndime,npoin,nltfl,nxgll,nygll,nspec)
  call positrec(coord,posrec,ndime,npoin,nrec)

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
  call qmasspec(rhoext,wx,wy,ibool,dvolu,rmass,density,kmato,npoin)

!
!---- definir les tableaux a1 a a13
!
  call defarrays(vpext,vsext,rhoext,density,elastcoef, &
                    xi,yi,wx,wy,hprime,hTprime, &
                    a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13x,a13z, &
                    ibool,iboolori,kmato,dvolu,xjaci,coord,gltfu, &
                    numabs,codeabs,anyabs,anyperio)

! initialiser les tableaux a zero
  accel = zero
  veloc = zero
  displ = zero
  vpred = zero
  force = zero

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
  anglerec2 = anglerec2 * pi / 180.d0

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
      allocate(displread(ndime))
      allocate(velocread(ndime))
      allocate(accelread(ndime))
      do n = 1,npoin
        read(55,*) inump, (displread(i), i=1,ndime), &
            (velocread(i), i=1,ndime), (accelread(i), i=1,ndime)
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
!---- verification des fonctions temporelles des sources
!
  if(.not. initialfield) call checksource(gltfu,nltfl,deltat,ncycl)

!
!---- verifier le maillage, la stabilite et le nb de points par lambda
!
  call checkgrid(deltat,gltfu,nltfl,initialfield)

!
!---- if data check mode then stop
!
  if(iexec  ==  0) then
   print *,'**********************************'
   print *,'* Aborting, data check mode only *'
   print *,'**********************************'
   call datim('Program  S P E C F E M :  end data checking',stitle,iout)
   stop
  endif

!
!---- initialiser sismogrammes
!
  sisux = sngl(zero)
  sisuz = sngl(zero)

  dcosrot1 = dcos(anglerec)
  dsinrot1 = dsin(anglerec)
  dcosrot2 = dcos(anglerec2)
  dsinrot2 = dsin(anglerec2)

!
!---- ouvrir fichier pour courbe d'energie
!
  if(compenergy) open(unit=ienergy,file='energy.gnu',status='unknown')

!
!----          s t a r t   t i m e   i t e r a t i o n s
!

  write(iout,400)

! boucle principale d'evolution en temps
  do it=1,ncycl

  if(mod(it-1,iaffinfo)  ==  0) then
  time = (it-1)*deltat
  if(time  >=  1.d-3) then
    write(iout,100) it,time
  else
    write(iout,101) it,time
  endif
  endif

! calculer le predictor
  displ(:,:) = displ(:,:) + deltat*veloc(:,:) + deltatsqover2*accel(:,:)
  vpred(:,:) = veloc(:,:) + deltatover2*accel(:,:)

! initialisation pour les iterations
  veloc(:,:) = vpred(:,:)

! calculer le terme source
  call calcforce(force,ndime,gltfu,nltfl,it*deltat)

! iteration sur le residu d'acceleration
  do iter = 1,niter

  accel(:,:) = zero

!
!----  calcul du residu d'acceleration pour le multicorrector
!----  retourne dans accel le terme Fext - M*A(i,n+1) - K*D(i,n+1)
!
  time = it*deltat
  call qsumspec(hprime,hTprime, &
          a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13x,a13z,force, &
          ibool,displ,veloc,accel, &
          Uxnewloc,Uznewloc,rmass,nxgll,npoin,ndime, &
          nspec,gltfu,nltfl,initialfield, &
          numabs,is_bordabs,nelemabs,anyabs)

!
!----  mise a jour globale du deplacement par multicorrector
!

  veloc(:,:) = vpred(:,:) + deltatover2*accel(:,:)

  enddo

!
!-----   calcul de l'energie cinetique et potentielle
!
  if(compenergy) &
      call calc_energie(hprime,hTprime,ibool,displ,veloc, &
      Uxnewloc,Uznewloc,kmato,dvolu,xjaci,density,elastcoef, &
      wx,wy,nxgll,npoin,ndime,nspec,numat)

!
!----  afficher le max du deplacement a certains pas de temps
!
  if(mod(it-1,iaffinfo)  ==  0) &
     print *,'Max norme U = ',maxval(sqrt(displ(1,:)**2 + displ(2,:)**2))

!
!----  affichage des resultats a certains pas de temps
!
  if (display .and. it  >  1 .and. (mod(it-1,itaff)  ==  0 .or. &
               it  ==  itfirstaff .or. it  ==  ncycl)) then

  time = it*deltat
  write(iout,*)
  if(time  >=  1.d-3) then
  write(iout,110) time
  else
  write(iout,111) time
  endif
  write(iout,*)

!
!----  affichage postscript
!
  if (ivectplot) then
  write(iout,*) 'Dump PostScript'
  if(ivecttype  ==  1) then
    write(iout,*) 'Drawing displacement field...'
    call plotpost(displ,coord,vpext,gltfu,posrec, &
          nltfl,it,deltat,coorg,xinterp,zinterp,shapeint, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,codeperio,anyabs,anyperio)
  else if(ivecttype  ==  2) then
    write(iout,*) 'Drawing velocity field...'
    call plotpost(veloc,coord,vpext,gltfu,posrec, &
          nltfl,it,deltat,coorg,xinterp,zinterp,shapeint, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,codeperio,anyabs,anyperio)
  else if(ivecttype  ==  3) then
    write(iout,*) 'Drawing acceleration field...'
    call plotpost(accel,coord,vpext,gltfu,posrec, &
          nltfl,it,deltat,coorg,xinterp,zinterp,shapeint, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,codeperio,anyabs,anyperio)
  else
    stop 'Wrong field code for PostScript display'
  endif
  write(iout,*) 'Fin dump PostScript'
  endif

!
!----  generation fichier AVS
!
  if(iavs) then
    if(anyperio) then
      call plotavs(displ,coord,kmato,iboolori,it)
    else
      call plotavs(displ,coord,kmato,ibool,it)
    endif
  endif

  endif

! stockage des sismogrammes
  if(sismos .and. (mod(it-1,isamp)  ==  0 .or. it  ==  ncycl)) then

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

! distinguer les deux lignes de recepteurs
  if(irec  <=  nrec1) then
    dcosrot = dcosrot1
    dsinrot = dsinrot1
  else
    dcosrot = dcosrot2
    dsinrot = dsinrot2
  endif

! rotation eventuelle des composantes
  itsis = min(it/isamp + 1,nseis)
  sisux(itsis,irec) =   sngl(dcosrot*valux + dsinrot*valuz)
  sisuz(itsis,irec) = - sngl(dsinrot*valux + dcosrot*valuz)

  enddo

  endif

  enddo

!
!----  sauvegarder sismogrammes en fin de simulation
!
  if(sismos) call writeseis(sisux,sisuz,coord,posrec,ndime,npoin,nseis,nrec, &
     isamp,deltat,factorxsu,n1ana,n2ana,irepr,nrec1,nrec2,isismostype)

!
!----  fermer fichier pour courbe d'energie et creer un petit script gnuplot
!
  if(compenergy) then
    close(ienergy)
    open(unit=ienergy,file='plotenergy',status='unknown')
    write(ienergy,*) 'set term postscript landscape color solid "Helvetica" 22'
    write(ienergy,*) 'set output "energy.ps"'
    write(ienergy,*) 'set xlabel "Time (s)"'
    write(ienergy,*) 'set ylabel "Energy (J)"'
    write(ienergy,*) 'plot "energy.gnu" us 1:4 t ''Total Energy'' w l 1, "energy.gnu" us 1:3 t ''Potential Energy'' w l 2'
    close(ienergy)
  endif

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!
!----  desallouer tous les tableaux avant de terminer l'execution
!
  deallocate(sisux)
  deallocate(sisuz)
  deallocate(posrec)
  deallocate(coorg)
  deallocate(coord)
  deallocate(force)
  deallocate(gltfu)
  deallocate(accel)
  deallocate(displ)
  deallocate(veloc)
  deallocate(vpred)
  deallocate(rmass)
  deallocate(fglobx)
  deallocate(fglobz)
  deallocate(shape)
  deallocate(shapeint)
  deallocate(dershape)
  deallocate(dvolu)
  deallocate(xjaci)
  deallocate(hprime)
  deallocate(hTprime)
  deallocate(ibool)
  deallocate(a1)
  deallocate(a2)
  deallocate(a3)
  deallocate(a4)
  deallocate(a5)
  deallocate(a6)
  deallocate(a7)
  deallocate(a8)
  deallocate(a9)
  deallocate(a10)
  deallocate(a13x)
  deallocate(a13z)
  deallocate(a11)
  deallocate(a12)
  deallocate(xi)
  deallocate(yi)
  deallocate(wx)
  deallocate(wy)
  deallocate(Uxnewloc)
  deallocate(Uznewloc)
  deallocate(xirec)
  deallocate(etarec)
  deallocate(flagrange)
  deallocate(xinterp)
  deallocate(zinterp)
  deallocate(Uxinterp)
  deallocate(Uzinterp)
  deallocate(density)
  deallocate(elastcoef)
  deallocate(kmato)
  deallocate(knods)
  deallocate(numabs)
  deallocate(codeabs)
  deallocate(codeperio)

  call datim('Program  S P E C F E M :  end',stitle,iout)

!
!----  close output file
!
  close(iout)

  stop

!
!----  formats
!
 40   format(a80)
 45   format(a50)
 100  format('Pas de temps numero ',i5,'   t = ',f7.4,' s')
 101  format('Pas de temps numero ',i5,'   t = ',1pe10.4,' s')
 110  format('Sauvegarde deplacement temps t = ',f7.4,' s')
 111  format('Sauvegarde deplacement temps t = ',1pe10.4,' s')
 400  format(/1x,41('=')/,' =  T i m e  ', &
  'e v o l u t i o n  l o o p  ='/1x,41('=')/)

  end program main
