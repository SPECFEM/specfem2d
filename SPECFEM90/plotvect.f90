
!=====================================================================
!
!                 S p e c f e m  V e r s i o n  4 . 2
!                 -----------------------------------
!
!                         Dimitri Komatitsch
!    Department of Earth and Planetary Sciences - Harvard University
!                         Jean-Pierre Vilotte
!                 Departement de Sismologie - IPGP - Paris
!                           (c) June 1998
!
!=====================================================================

  subroutine plotvect(knods,coorg,coord,displ,density,elastcoef, &
     kmato,flagrange,xinterp,zinterp,shapeint,Uxinterp,Uzinterp,ibool,vpext, &
    numabs,codeabs,codeperio,anyabs,anyperio)
!
!=======================================================================
!
!     "p l o t v e c t" : Print the displacement vector field
!                         in a PostScript file together with
!                         the spectral elements boundaries
!
!=======================================================================
!
  use palette
  use constspec
  use verifs
  use mesh01
  use spela202
  use defpi
  use codebord

  implicit none

  double precision coorg(ndime,npgeo)
  double precision displ(ndofn,npoin)
  double precision coord(ndime,npoin)
  double precision xinterp(iptsdisp,iptsdisp),zinterp(iptsdisp,iptsdisp)
  double precision shapeint(ngnod,iptsdisp,iptsdisp)
  double precision Uxinterp(iptsdisp,iptsdisp)
  double precision Uzinterp(iptsdisp,iptsdisp)
  double precision flagrange(0:nxgll-1,iptsdisp)
  double precision density(numat),elastcoef(4,numat)
  double precision vpext(npoin)

  integer knods(ngnod,nspec),kmato(nspec)
  integer ibool(0:nxgll-1,0:nygll-1,nspec)

  integer numabs(nelemabs),codeabs(4,nelemabs)
  integer codeperio(4,nelemperio)
  logical anyabs,anyperio

  character(len=100) name
  character ch1(100),ch2(100)
  equivalence (name,ch1)
  logical first

  double precision convert,x1,rlamda,rmu,denst,rKvol,cploc,xw,zw,xa,za,xb,zb
  double precision z1,x2,z2,d,d1,d2,dummy,theta,thetaup,thetadown

  integer i,k,j,ispec,material,ispel,is,ir,nbcols,imat,icol,l,longueur
  integer indice,ii,ipoin,in,nnum,ispelabs,ideb,ifin,ibord
  integer numelem,nedgeloc,num2,nedgeother,n

!
!-----------------------------------------------------------------------
!
!---- print the spectral elements mesh in PostScript
!

  print *,'Shape functions based on ',ngnod,' control nodes'

  convert = pi/180.d0

!
!----  draw the velocity model in background
!
  if(imodelvect) then

  do ispec=1,nspec
    do i=0,nxgll-1-isubsamp,isubsamp
          do j=0,nxgll-1-isubsamp,isubsamp

  if((vpmax-vpmin)/vpmin  >  0.02d0) then
  if(ireadmodel) then
    x1 = (vpext(ibool(i,j,ispec))-vpmin)/ (vpmax-vpmin)
  else
 material = kmato(ispec)
 rlamda = elastcoef(1,material)
 rmu    = elastcoef(2,material)
 denst  = density(material)
 rKvol  = rlamda + 2.d0*rmu/3.d0
 cploc = dsqrt((rKvol + 4.d0*rmu/3.d0)/denst)
 x1 = (cploc-vpmin)/(vpmax-vpmin)
  endif
  else
    x1 = 0.5d0
  endif

! rescaler pour eviter gris trop sombre
  x1 = x1*0.7 + 0.2
  if (x1  >  1.d0) x1=1.d0

! inverser echelle : blanc = vpmin, gris = vpmax
  x1 = 1.d0 - x1

  xw = coord(1,ibool(i,j,ispec))
  zw = coord(2,ibool(i,j,ispec))
  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
    write(24,500) xw,zw
  xw = coord(1,ibool(i+isubsamp,j,ispec))
  zw = coord(2,ibool(i+isubsamp,j,ispec))
  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
    write(24,499) xw,zw
  xw = coord(1,ibool(i+isubsamp,j+isubsamp,ispec))
  zw = coord(2,ibool(i+isubsamp,j+isubsamp,ispec))
  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
    write(24,499) xw,zw
  xw = coord(1,ibool(i,j+isubsamp,ispec))
  zw = coord(2,ibool(i,j+isubsamp,ispec))
  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
    write(24,499) xw,zw
    write(24,604) x1

          enddo
    enddo
  enddo

  endif

!
!---- draw spectral element mesh
!

  if (imeshvect) then

  write(24,*) '%'
  write(24,*) '% spectral element mesh'
  write(24,*) '%'

  do ispel=1,nspec

  write(24,*) '% elem ',ispel

  do i=1,iptsdisp
  do j=1,iptsdisp
  xinterp(i,j) = 0.d0
  zinterp(i,j) = 0.d0
  do in = 1,ngnod
    nnum = knods(in,ispel)
      xinterp(i,j) = xinterp(i,j) + shapeint(in,i,j)*coorg(1,nnum)
      zinterp(i,j) = zinterp(i,j) + shapeint(in,i,j)*coorg(2,nnum)
  enddo
  enddo
  enddo

  is = 1
  ir = 1
  x1 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z1 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  write(24,*) 'MK'
  write(24,601) x1,z1

  if (ngnod  ==  4) then

! tracer des droites si elements 4 noeuds

  ir=iptsdisp
  x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,601) x2,z2

  ir=iptsdisp
  is=iptsdisp
  x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,601) x2,z2

  is=iptsdisp
  ir=1
  x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,601) x2,z2

  ir=1
  is=2
  x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,601) x2,z2

  else

! tracer des courbes si elements 9 noeuds
  do ir=2,iptsdisp
  x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,601) x2,z2
  enddo

  ir=iptsdisp
  do is=2,iptsdisp
  x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,601) x2,z2
  enddo

  is=iptsdisp
  do ir=iptsdisp-1,1,-1
  x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,601) x2,z2
  enddo

  ir=1
  do is=iptsdisp-1,2,-1
  x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,601) x2,z2
  enddo

  endif

  write(24,*) 'CO'

  if (icolor  ==  1) then

! For the moment 20 different colors max
  nbcols = 20

! Use a different color for each material set
  imat = kmato(ispel)
  icol = mod(imat - 1,nbcols) + 1

  write(24,600) red(icol),green(icol),blue(icol)

  endif

  if(imodelvect) then
  write(24,*) 'GC'
  else
  write(24,*) 'GG'
  endif

! write the element number, the group number and the
! material number inside the element
  if (inumber  ==  1) then

  xw = (coorg(1,knods(1,ispel)) + coorg(1,knods(2,ispel)) + &
          coorg(1,knods(3,ispel)) + coorg(1,knods(4,ispel))) / 4.d0
  zw = (coorg(2,knods(1,ispel)) + coorg(2,knods(2,ispel)) + &
          coorg(2,knods(3,ispel)) + coorg(2,knods(4,ispel))) / 4.d0
  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
  if (icolor  ==  1) write(24,*) '1 setgray'

  write(24,500) xw,zw

!--- ecriture numero de l'element
  write(24,502) ispel

  endif

  enddo

  endif

!
!----  draw the boundary conditions
!

  if((anyabs .or. anyperio) .and. iboundvect) then

  write(24,*) '%'
  write(24,*) '% boundary conditions on the mesh'
  write(24,*) '%'

  write(24,*) '0.05 CM SLW'

!--- bords absorbants

  if(anyabs) then

  do ispelabs = 1,nelemabs
  ispel = numabs(ispelabs)

!--- une couleur pour chaque condition absorbante
!--- bord absorbant de type "haut"   : orange
!--- bord absorbant de type "bas"    : vert clair
!--- bord absorbant de type "gauche" : rose clair
!--- bord absorbant de type "droite" : turquoise

  do ibord = 1,4

  if(codeabs(ibord,ispelabs)  /=  0) then

  if(ibord  ==  ihaut) then
    write(24,*) '1. .85 0. RG'
    ideb = 3
    ifin = 4
  else if(ibord  ==  ibas) then
    write(24,*) '.4 1. .4 RG'
    ideb = 1
    ifin = 2
  else if(ibord  ==  igauche) then
    write(24,*) '1. .43 1. RG'
    ideb = 4
    ifin = 1
  else if(ibord  ==  idroite) then
    write(24,*) '.25 1. 1. RG'
    ideb = 2
    ifin = 3
  else
    stop 'Wrong absorbing boundary code'
  endif

  x1 = (coorg(1,knods(ideb,ispel))-xmin)*rapp_page + orig_x
  z1 = (coorg(2,knods(ideb,ispel))-zmin)*rapp_page + orig_z
  x2 = (coorg(1,knods(ifin,ispel))-xmin)*rapp_page + orig_x
  z2 = (coorg(2,knods(ifin,ispel))-zmin)*rapp_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,602) x1,z1,x2,z2

  endif
  enddo

  enddo

  endif

!--- bords periodiques dessines en rouge

  if(anyperio) then

  write(24,*) '1. .15 0.25 RG'

  do n=1,nelemperio
    numelem    = codeperio(1,n)
    nedgeloc   = codeperio(2,n)
    num2       = codeperio(3,n)
    nedgeother = codeperio(4,n)

! dessin premiere arete
  if(nedgeloc  ==  iaretehaut) then
    ideb = 3
    ifin = 4
  else if(nedgeloc  ==  iaretebas) then
    ideb = 1
    ifin = 2
  else if(nedgeloc  ==  iaretegauche) then
    ideb = 4
    ifin = 1
  else if(nedgeloc  ==  iaretedroite) then
    ideb = 2
    ifin = 3
  endif

  x1 = (coorg(1,knods(ideb,numelem))-xmin)*rapp_page + orig_x
  z1 = (coorg(2,knods(ideb,numelem))-zmin)*rapp_page + orig_z
  x2 = (coorg(1,knods(ifin,numelem))-xmin)*rapp_page + orig_x
  z2 = (coorg(2,knods(ifin,numelem))-zmin)*rapp_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,602) x1,z1,x2,z2

! dessin arete correspondante
  if(nedgeother  ==  iaretehaut) then
    ideb = 3
    ifin = 4
  else if(nedgeother  ==  iaretebas) then
    ideb = 1
    ifin = 2
  else if(nedgeother  ==  iaretegauche) then
    ideb = 4
    ifin = 1
  else if(nedgeother  ==  iaretedroite) then
    ideb = 2
    ifin = 3
  endif

  x1 = (coorg(1,knods(ideb,num2))-xmin)*rapp_page + orig_x
  z1 = (coorg(2,knods(ideb,num2))-zmin)*rapp_page + orig_z
  x2 = (coorg(1,knods(ifin,num2))-xmin)*rapp_page + orig_x
  z2 = (coorg(2,knods(ifin,num2))-zmin)*rapp_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,602) x1,z1,x2,z2

  enddo

  endif

  write(24,*) '0 setgray'
  write(24,*) '0.01 CM SLW'

  endif

!
!----  draw the normalized displacement field
!

! return if the maximum displacement equals zero (no source)
  if (dispmax  ==  0.d0) then
    print *,' null displacement : returning !'
    return
  endif

  write(24,*) '%'
  write(24,*) '% vector field'
  write(24,*) '%'

! fleches en couleur si modele de vitesse en background
  if(imodelvect) then
        write(24,*) 'Colvects'
  else
        write(24,*) '0 setgray'
  endif

  if (interpol) then

  print *,'Interpolating the vector field...'

  do ispel=1,nspec

! interpolation sur grille reguliere
  if(mod(ispel,100)  ==  0) &
       write(*,*) 'Interpolation uniform grid element ',ispel

  do i=1,iptsdisp
  do j=1,iptsdisp

  xinterp(i,j) = 0.d0
  zinterp(i,j) = 0.d0
  do in = 1,ngnod
    nnum = knods(in,ispel)
      xinterp(i,j) = xinterp(i,j) + shapeint(in,i,j)*coorg(1,nnum)
      zinterp(i,j) = zinterp(i,j) + shapeint(in,i,j)*coorg(2,nnum)
  enddo

  Uxinterp(i,j) = 0.d0
  Uzinterp(i,j) = 0.d0

  do k=0,nxgll-1
  do l=0,nxgll-1

  Uxinterp(i,j) = Uxinterp(i,j) + &
                displ(1,ibool(k,l,ispel))*flagrange(k,i)*flagrange(l,j)
  Uzinterp(i,j) = Uzinterp(i,j) + &
                displ(2,ibool(k,l,ispel))*flagrange(k,i)*flagrange(l,j)

  enddo
  enddo

  x1 =(xinterp(i,j)-xmin)*rapp_page
  z1 =(zinterp(i,j)-zmin)*rapp_page

  x2 = Uxinterp(i,j)*sizemax/dispmax
  z2 = Uzinterp(i,j)*sizemax/dispmax

  d = dsqrt(x2**2 + z2**2)

! ignorer si vecteur trop petit
  if (d  >  cutvect*sizemax) then

  d1 = d * rapport
  d2 = d1 * dcos(angle*convert)

  dummy = x2/d
  if (dummy  >  0.9999d0) dummy = 0.9999d0
  if (dummy  <  -0.9999d0) dummy = -0.9999d0
  theta = dacos(dummy)

  if(z2  <  0.d0) theta = 360.d0*convert - theta
  thetaup = theta - angle*convert
  thetadown = theta + angle*convert

! tracer le vecteur proprement dit
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  x2 = x2 * centim
  z2 = z2 * centim
  xa = -d2*dcos(thetaup)
  za = -d2*dsin(thetaup)
  xa = xa * centim
  za = za * centim
  xb = -d2*dcos(thetadown)
  zb = -d2*dsin(thetadown)
  xb = xb * centim
  zb = zb * centim
  write(name,700) xb,zb,xa,za,x2,z2,x1,z1

! filtrer les blancs inutiles pour diminuer taille fichier PostScript
  longueur = 49
  indice = 1
  first = .false.
  do ii=1,longueur-1
    if(ch1(ii)  /=  ' '.or.first) then
    if(ch1(ii)  /=  ' '.or.ch1(ii+1)  /=  ' ') then
          ch2(indice) = ch1(ii)
          indice = indice + 1
          first = .true.
    endif
    endif
  enddo
  ch2(indice) = ch1(longueur)
  write(24,200) (ch2(ii),ii=1,indice)

  endif

  enddo
  enddo
  enddo

  else
! tracer les vecteurs deplacement aux noeuds du maillage

  do ipoin=1,npoin

  x1 =(coord(1,ipoin)-xmin)*rapp_page
  z1 =(coord(2,ipoin)-zmin)*rapp_page

  x2 = displ(1,ipoin)*sizemax/dispmax
  z2 = displ(2,ipoin)*sizemax/dispmax

  d = dsqrt(x2**2 + z2**2)

! ignorer si vecteur trop petit
  if (d  >  cutvect*sizemax) then

  d1 = d * rapport
  d2 = d1 * dcos(angle*convert)

  dummy = x2/d
  if (dummy  >  0.9999d0) dummy = 0.9999d0
  if (dummy  <  -0.9999d0) dummy = -0.9999d0
  theta = dacos(dummy)

  if(z2  <  0.d0) theta = 360.d0*convert - theta
  thetaup = theta - angle*convert
  thetadown = theta + angle*convert

! tracer le vecteur proprement dit
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  x2 = x2 * centim
  z2 = z2 * centim
  xa = -d2*dcos(thetaup)
  za = -d2*dsin(thetaup)
  xa = xa * centim
  za = za * centim
  xb = -d2*dcos(thetadown)
  zb = -d2*dsin(thetadown)
  xb = xb * centim
  zb = zb * centim
  write(name,700) xb,zb,xa,za,x2,z2,x1,z1

! filtrer les blancs inutiles pour diminuer taille fichier PostScript
  longueur = 49
  indice = 1
  first = .false.
  do ii=1,longueur-1
    if(ch1(ii)  /=  ' '.or.first) then
    if(ch1(ii)  /=  ' '.or.ch1(ii+1)  /=  ' ') then
          ch2(indice) = ch1(ii)
          indice = indice + 1
          first = .true.
    endif
    endif
  enddo
  ch2(indice) = ch1(longueur)
  write(24,200) (ch2(ii),ii=1,indice)

  endif

  enddo

  endif

  write(24,*) '0 setgray'

 200  format(80(a1))
 499  format(f5.1,1x,f5.1,' L')
 500  format(f5.1,1x,f5.1,' M')
 502  format('fN (',i4,') Cshow')
 600  format(f4.2,1x,f4.2,1x,f4.2,' RG GF')
 601  format(f6.2,1x,f6.2)
 602  format(f6.2,1x,f6.2,' M ',f6.2,1x,f6.2,' L ST')
 604  format('CP ',f4.2,' BK')
 700  format(8(f5.1,1x),'F')

  end subroutine plotvect
