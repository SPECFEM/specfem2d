
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

! version 5.0 : - got rid of useless routines, suppressed commons etc.
!               - weak formulation based explicitly on stress tensor
!               - implementation of full anisotropy

! version 5.0 is based on SPECFEM2D version 4.2
! (c) June 1998 by Dimitri Komatitsch, Harvard University, USA
! and Jean-Pierre Vilotte, Institut de Physique du Globe de Paris, France

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
  integer nbpoin,inump,n,npoinext,netyp,ispec,npoin,npgeo,iglob

  double precision valux,valuz,rhoextread,vpextread,vsextread
  double precision cpl,csl,rhol
  double precision dcosrot,dsinrot,xcor,zcor

! coefficients of the explicit Newmark time scheme
  integer NSTEP
  double precision deltatover2,deltatsquareover2,time,deltat

! Gauss-Lobatto-Legendre points and weights
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLZ) :: yigll,wzgll

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
  double precision mul,lambdal,lambdalplus2mul

  double precision, dimension(:), allocatable :: xirec,etarec

  double precision, dimension(:,:), allocatable :: coord,accel,veloc,displ, &
    flagrange,xinterp,zinterp,Uxinterp,Uzinterp,elastcoef

  double precision, dimension(:), allocatable :: rmass, &
    fglobx,fglobz,density,vpext,vsext,rhoext,displread,velocread,accelread

  double precision, dimension(:,:,:), allocatable :: shapeint,shape, &
    xix,xiz,gammax,gammaz,jacobian,a13x,a13z

  double precision, dimension(:,:), allocatable :: a11,a12

  double precision, dimension(:,:,:,:), allocatable :: dershape

  integer, dimension(:,:,:), allocatable :: ibool
  integer, dimension(:,:), allocatable  :: knods
  integer, dimension(:), allocatable :: kmato,numabs

  integer ie,k,isource,iexplo

  integer ielems,iglobsource
  double precision f0,t0,factor,a,angleforce,ricker

  double precision rsizemin,rsizemax,cpoverdxmin,cpoverdxmax, &
    lambdalSmin,lambdalSmax,lambdalPmin,lambdalPmax,vpmin,vpmax

  integer icolor,inumber,isubsamp,ivecttype,itaff,nrec,isismostype
  integer numat,ngnod,nspec,iptsdisp,nelemabs

  logical interpol,imeshvect,imodelvect,iboundvect,ireadmodel,initialfield, &
    ioutputgrid,ignuplot

  double precision cutvect,anglerec

! for absorbing conditions
  integer ispecabs,inum,numabsread,i1abs,i2abs
  logical codeabsread(4)
  double precision nx,nz,vx,vz,vn,rho_vp,rho_vs,tx,tz,weight,xxi,zeta,rKvol

  logical, dimension(:,:), allocatable  :: codeabs

! title of the plot
  character(len=60) stitle

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
   if(isource >= 4 .and. isource <= 6 .and. iexplo == 1) gltfu(8) = gltfu(8) * pi / 180.d0

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
  allocate(shape(ngnod,NGLLX,NGLLZ))
  allocate(shapeint(ngnod,iptsdisp,iptsdisp))
  allocate(dershape(NDIME,ngnod,NGLLX,NGLLZ))
  allocate(xix(NGLLX,NGLLZ,nspec))
  allocate(xiz(NGLLX,NGLLZ,nspec))
  allocate(gammax(NGLLX,NGLLZ,nspec))
  allocate(gammaz(NGLLX,NGLLZ,nspec))
  allocate(jacobian(NGLLX,NGLLZ,nspec))
  allocate(a11(NGLLX,NGLLZ))
  allocate(a12(NGLLX,NGLLZ))
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
  allocate(ibool(NGLLX,NGLLZ,nspec))

! --- allocate arrays for absorbing boundary conditions
  if(nelemabs <= 0) then
    nelemabs = 1
    anyabs = .false.
  else
    anyabs = .true.
  endif
  allocate(numabs(nelemabs))
  allocate(codeabs(4,nelemabs))

!
!---- print element group main parameters
!
  write(IOUT,107)
  write(IOUT,207) nspec,ngnod,NGLLX,NGLLZ,NGLLX*NGLLZ,iptsdisp,numat,nelemabs

! set up Gauss-Lobatto-Legendre derivation matrices
  call define_derivative_matrices(xigll,yigll,wxgll,wzgll,hprime_xx,hprime_zz)

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
      read(IIN ,*) inum,numabsread,codeabsread(1),codeabsread(2),codeabsread(3),codeabsread(4)
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
!---- compute the spectral element shape functions and their local derivatives
!
  call q49shape(shape,dershape,xigll,yigll,ngnod)

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

  call q49spec(shapeint,dershape,xix,xiz,gammax,gammaz,jacobian,xigll, &
          coorg,knods,ngnod,nspec,npgeo,xirec,etarec,flagrange,iptsdisp)

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

  allocate(a13x(NGLLX,NGLLZ,nelemabs))
  allocate(a13z(NGLLX,NGLLZ,nelemabs))

!
!----  set the coordinates of the points of the global grid
!
  do ispec = 1,nspec
  do ip1 = 1,NGLLX
  do ip2 = 1,NGLLZ

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
  deltatover2 = HALF*deltat
  deltatsquareover2 = HALF*deltat*deltat

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
!---- define all arrays
!
  call defarrays(vpext,vsext,rhoext,density,elastcoef, &
          xigll,yigll,xix,xiz,gammax,gammaz,a11,a12, &
          ibool,kmato,coord,gltfu,npoin,rsizemin,rsizemax, &
          cpoverdxmin,cpoverdxmax,lambdalSmin,lambdalSmax, &
          lambdalPmin,lambdalPmax,vpmin,vpmax,ireadmodel,nspec,numat)

! build the global mass matrix once and for all
  rmass(:) = ZERO
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
!--- if external density model
        if(ireadmodel) then
          rhol = rhoext(iglob)
        else
          rhol = density(kmato(ispec))
        endif
        rmass(iglob) = rmass(iglob) + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)
      enddo
    enddo
  enddo

! convertir angle recepteurs en radians
  anglerec = anglerec * pi / 180.d0

!
!---- verifier le maillage, la stabilite et le nb de points par lambda
!
  call checkgrid(deltat,gltfu,initialfield,rsizemin,rsizemax, &
      cpoverdxmin,cpoverdxmax,lambdalSmin,lambdalSmax,lambdalPmin,lambdalPmax)

!
!---- initialiser sismogrammes
!
  sisux = ZERO
  sisuz = ZERO

  dcosrot = dcos(anglerec)
  dsinrot = dsin(anglerec)

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
    print *,'Max norm of initial displacement = ',maxval(sqrt(displ(1,:)**2 + displ(2,:)**2))
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

! `predictor' update displacement using finite-difference time scheme (Newmark)
    displ(:,:) = displ(:,:) + deltat*veloc(:,:) + deltatsquareover2*accel(:,:)
    veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)
    accel(:,:) = ZERO

!   integration over spectral elements
    do ispec = 1,NSPEC

! get elastic parameters of current spectral element
      lambdal = elastcoef(1,kmato(ispec))
      mul = elastcoef(2,kmato(ispec))
      lambdalplus2mul = lambdal + 2.d0*mul

! first double loop over GLL to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX

!--- if external medium, get elastic parameters of current grid point
          if(ireadmodel) then
            iglob = ibool(i,j,ispec)
            cpl = vpext(iglob)
            csl = vsext(iglob)
            rhol = rhoext(iglob)
            mul = rhol*csl*csl
            lambdal = rhol*cpl*cpl - 2.d0*mul
            lambdalplus2mul = lambdal + 2.d0*mul
          endif

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

! compute stress tensor
          sigma_xx = lambdalplus2mul*duxdxl + lambdal*duzdzl
          sigma_xz = mul*(duzdxl + duxdzl)
          sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl

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
            tempz1l = tempz1l + tempz1(k,j)*fac1
          enddo

! along z direction
          tempx2l = ZERO
          tempz2l = ZERO
          do k = 1,NGLLZ
            fac2 = wzgll(k)*hprime_zz(j,k)
            tempx2l = tempx2l + tempx2(i,k)*fac2
            tempz2l = tempz2l + tempz2(i,k)*fac2
          enddo

          fac1 = wzgll(j)
          fac2 = wxgll(i)

          iglob = ibool(i,j,ispec)
          accel(1,iglob) = accel(1,iglob) - (fac1*tempx1l + fac2*tempx2l)
          accel(2,iglob) = accel(2,iglob) - (fac1*tempz1l + fac2*tempz2l)

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
      lambdal = elastcoef(1,kmato(ispec))
      mul = elastcoef(2,kmato(ispec))
      rhol  = density(kmato(ispec))
      rKvol  = lambdal + 2.d0*mul/3.d0
      cpl = dsqrt((rKvol + 4.d0*mul/3.d0)/rhol)
      csl = dsqrt(mul/rhol)


!--- left absorbing boundary
      if(codeabs(ILEFT,ispecabs)) then

        i = 1

        do j=1,NGLLZ

          iglob = ibool(i,j,ispec)

          zeta = xix(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(ireadmodel) then
            cpl = vpext(iglob)
            csl = vsext(iglob)
            rhol = rhoext(iglob)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          nx = -1.d0
          nz = 0.d0

          vx = veloc(1,iglob)
          vz = veloc(2,iglob)

          vn = nx*vx+nz*vz

          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

          weight = zeta*wzgll(j)

          accel(1,iglob) = accel(1,iglob) - tx*weight
          accel(2,iglob) = accel(2,iglob) - tz*weight

        enddo

      endif  !  end of left absorbing boundary

!--- right absorbing boundary
      if(codeabs(IRIGHT,ispecabs)) then

        i = NGLLX

        do j=1,NGLLZ

          iglob = ibool(i,j,ispec)

          zeta = xix(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(ireadmodel) then
            cpl = vpext(iglob)
            csl = vsext(iglob)
            rhol = rhoext(iglob)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          nx = 1.d0
          nz = 0.d0

          vx = veloc(1,iglob)
          vz = veloc(2,iglob)

          vn = nx*vx+nz*vz

          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

          weight = zeta*wzgll(j)

          accel(1,iglob) = accel(1,iglob) - tx*weight
          accel(2,iglob) = accel(2,iglob) - tz*weight

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
          if(ireadmodel) then
            cpl = vpext(iglob)
            csl = vsext(iglob)
            rhol = rhoext(iglob)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          nx = 0.d0
          nz = -1.d0

          vx = veloc(1,iglob)
          vz = veloc(2,iglob)

          vn = nx*vx+nz*vz

          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

          weight = xxi*wxgll(i)

          accel(1,iglob) = accel(1,iglob) - tx*weight
          accel(2,iglob) = accel(2,iglob) - tz*weight

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
          if(ireadmodel) then
            cpl = vpext(iglob)
            csl = vsext(iglob)
            rhol = rhoext(iglob)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          nx = 0.d0
          nz = 1.d0

          vx = veloc(1,iglob)
          vz = veloc(2,iglob)

          vn = nx*vx+nz*vz

          tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
          tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

          weight = xxi*wxgll(i)

          accel(1,iglob) = accel(1,iglob) - tx*weight
          accel(2,iglob) = accel(2,iglob) - tz*weight

        enddo

      endif  !  end of top absorbing boundary

    enddo

  endif  ! end of absorbing boundaries


! --- add the source
  if(.not. initialfield) then

  f0 = gltfu(5)
  t0 = gltfu(6)
  factor = gltfu(7)
  angleforce = gltfu(8)

! Ricker wavelet for the source time function
  a = pi*pi*f0*f0
  ricker = - factor * (1.d0-2.d0*a*(time-t0)**2)*exp(-a*(time-t0)**2)

! --- collocated force
  if(nint(gltfu(2)) == 1) then
    iglobsource = nint(gltfu(9))
    accel(1,iglobsource) = accel(1,iglobsource) - dsin(angleforce)*ricker
    accel(2,iglobsource) = accel(2,iglobsource) + dcos(angleforce)*ricker

!---- explosion
  else if(nint(gltfu(2)) == 2) then
!   recuperer numero d'element de la source
    ielems = nint(gltfu(12))
    do i=1,NGLLX
      do j=1,NGLLX
        iglob = ibool(i,j,ielems)
        accel(1,iglob) = accel(1,iglob) + a11(i,j)*ricker
        accel(2,iglob) = accel(2,iglob) + a12(i,j)*ricker
      enddo
    enddo
  endif

  else
    stop 'wrong source type'
  endif

! divide by the mass matrix
  accel(1,:) = accel(1,:) / rmass(:)
  accel(2,:) = accel(2,:) / rmass(:)

! `corrector' update velocity
  veloc(:,:) = veloc(:,:) + deltatover2*accel(:,:)

!
!----  display max of norm of displacement
!
  if(mod(it,itaff) == 0) print *,'Max norm of displacement = ',maxval(sqrt(displ(1,:)**2 + displ(2,:)**2))

! store the seismograms
  do irec=1,nrec
    iglobrec = nint(posrec(1,irec))

  if(isismostype == 1) then
    valux = displ(1,iglobrec)
    valuz = displ(2,iglobrec)
  else if(isismostype == 2) then
    valux = veloc(1,iglobrec)
    valuz = veloc(2,iglobrec)
  else if(isismostype == 3) then
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
  if(time >= 1.d-3) then
    write(IOUT,110) time
  else
    write(IOUT,111) time
  endif
  write(IOUT,*)

!
!----  affichage postscript
!
  write(IOUT,*) 'Dump PostScript'
  if(ivecttype == 1) then
    write(IOUT,*) 'drawing displacement field...'
    call plotpost(displ,coord,vpext,gltfu,posrec, &
          it,deltat,coorg,xinterp,zinterp,shapeint, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,stitle,npoin,npgeo,vpmin,vpmax,nrec, &
          icolor,inumber,isubsamp,ivecttype,interpol,imeshvect,imodelvect, &
          iboundvect,ireadmodel,cutvect,nelemabs,numat,iptsdisp,nspec,ngnod)
  else if(ivecttype == 2) then
    write(IOUT,*) 'drawing velocity field...'
    call plotpost(veloc,coord,vpext,gltfu,posrec, &
          it,deltat,coorg,xinterp,zinterp,shapeint, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,stitle,npoin,npgeo,vpmin,vpmax,nrec, &
          icolor,inumber,isubsamp,ivecttype,interpol,imeshvect,imodelvect, &
          iboundvect,ireadmodel,cutvect,nelemabs,numat,iptsdisp,nspec,ngnod)
  else if(ivecttype == 3) then
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
     'Source Type. . . . . . . . . . . . . . = Explosion',/5x, &
     'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
     'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
     'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
     'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
     'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x)

  end program specfem2D

