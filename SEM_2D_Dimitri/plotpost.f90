
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

  subroutine plotpost(displ,coord,vpext,x_source,z_source,st_xval,st_zval,it,dt,coorg, &
          xinterp,zinterp,shapeint,Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs,stitle,npoin,npgeo,vpmin,vpmax,nrec, &
          colors,numbers,subsamp,vecttype,interpol,meshvect,modelvect, &
          boundvect,read_external_model,cutvect,nelemabs,numat,iptsdisp,nspec,ngnod,ELASTIC)

!
! routine affichage postscript
!

  implicit none

  include "constants.h"

! color palette
  integer, parameter :: MAXCOLORS = 100
  double precision, dimension(MAXCOLORS) :: red,green,blue

  integer it,nrec,nelemabs,numat,iptsdisp,nspec
  integer i,npoin,npgeo,ngnod

  integer kmato(nspec),knods(ngnod,nspec)
  integer ibool(NGLLX,NGLLZ,nspec)

  double precision xinterp(iptsdisp,iptsdisp),zinterp(iptsdisp,iptsdisp)
  double precision shapeint(ngnod,iptsdisp,iptsdisp)
  double precision Uxinterp(iptsdisp,iptsdisp)
  double precision Uzinterp(iptsdisp,iptsdisp)
  double precision flagrange(NGLLX,iptsdisp)
  double precision density(numat),elastcoef(4,numat)

  double precision dt,timeval,x_source,z_source
  double precision displ(NDIM,npoin),coord(NDIM,npoin)
  double precision vpext(npoin)

  double precision coorg(NDIM,npgeo)
  double precision, dimension(nrec) :: st_xval,st_zval

  integer numabs(nelemabs),codeabs(4,nelemabs)
  logical anyabs,ELASTIC

  double precision xmax,zmax,height,xw,zw,usoffset,sizex,sizez,vpmin,vpmax

  character(len=100) name
  character ch1(100),ch2(100)
  equivalence (name,ch1)
  logical first

  double precision convert,x1,rlamda,rmu,denst,rKvol,cploc,xa,za,xb,zb
  double precision z1,x2,z2,d,d1,d2,dummy,theta,thetaup,thetadown

  integer k,j,ispec,material,is,ir,nbcols,imat,icol,l,longueur
  integer indice,ii,ipoin,in,nnum,ispecabs,ideb,ifin,ibord

  integer colors,numbers,subsamp,vecttype
  logical interpol,meshvect,modelvect,boundvect,read_external_model
  double precision cutvect

  double precision rapp_page,dispmax,xmin,zmin

! title of the plot
  character(len=60) stitle

! papier A4 ou US letter
  if(US_LETTER) then
    usoffset = 1.75d0
    sizex = 27.94d0
    sizez = 21.59d0
  else
    usoffset = 0.d0
    sizex = 29.7d0
    sizez = 21.d0
  endif

! definition de la palette de couleur

! red
  red(1) = 1.d0
  green(1) = 0.d0
  blue(1) = 0.d0
! blue
  red(2) = 0.d0
  green(2) = 0.d0
  blue(2) = 1.d0
! violet
  red(3) = .93d0
  green(3) = .51d0
  blue(3) = .93d0
! medium orchid
  red(4) = .73d0
  green(4) = .33d0
  blue(4) = .83d0
! dark orchid
  red(5) = .6d0
  green(5) = .2d0
  blue(5) = .8d0
! blue violet
  red(6) = .54d0
  green(6) = .17d0
  blue(6) = .89d0
! slate blue
  red(7) = .42d0
  green(7) = .35d0
  blue(7) = .80d0
! deep pink
  red(8) = 1.d0
  green(8) = .08d0
  blue(8) = .58d0
! dodger blue
  red(9) = .12d0
  green(9) = .56d0
  blue(9) = 1.d0
! dark turquoise
  red(10) = 0.d0
  green(10) = .81d0
  blue(10) = .82d0
! turquoise
  red(11) = .25d0
  green(11) = .88d0
  blue(11) = .82d0
! lime green
  red(12) = .2d0
  green(12) = .8d0
  blue(12) = .2d0
! spring green
  red(13) = 0.d0
  green(13) = 1.d0
  blue(13) = .5d0
! chartreuse
  red(14) = .5d0
  green(14) = 1.d0
  blue(14) = 0.d0
! green yellow
  red(15) = .68d0
  green(15) = 1.d0
  blue(15) = .18d0
! yellow
  red(16) = 1.d0
  green(16) = 1.d0
  blue(16) = 0.d0
! lemon chiffon
  red(17) = 1.d0
  green(17) = .98d0
  blue(17) = .8d0
! gold
  red(18) = 1.d0
  green(18) = .84d0
  blue(18) = 0.d0
! mocassin
  red(19) = 1.d0
  green(19) = .89d0
  blue(19) = .71d0
! peach puff
  red(20) = 1.d0
  green(20) = .85d0
  blue(20) = .73d0

! recherche des positions maximales des points de la grille
  xmax=maxval(coord(1,:))
  zmax=maxval(coord(2,:))
  write(*,*) 'Max X = ',xmax
  write(*,*) 'Max Z = ',zmax

! limite du repere physique
  xmin=0.d0
  zmin=0.d0

! rapport taille page/taille domaine physique
  rapp_page = min(rpercentz*sizez/(zmax-zmin),rpercentx*sizex/(xmax-xmin)) / 100.d0

! recherche de la valeur maximum de la norme du deplacement
  dispmax = maxval(sqrt(displ(1,:)**2 + displ(2,:)**2))
  write(*,*) 'Max norme = ',dispmax

! hauteur des numeros de domaine en CM
  height = 0.25d0

!
!---- ouverture du fichier PostScript
!
  write(name,222) it
  open(unit=24,file=name,status='unknown')
  222 format('OUTPUT_FILES/vect',i5.5,'.ps')

!
!---- ecriture de l'entete du fichier PostScript
!
  write(24,10) stitle
  write(24,*) '/CM {28.5 mul} def'
  write(24,*) '/LR {rlineto} def'
  write(24,*) '/LT {lineto} def'
  write(24,*) '/L {lineto} def'
  write(24,*) '/MR {rmoveto} def'
  write(24,*) '/MV {moveto} def'
  write(24,*) '/M {moveto} def'
  write(24,*) '/MK {mark} def'
  write(24,*) '/ST {stroke} def'
  write(24,*) '/CP {closepath} def'
  write(24,*) '/RG {setrgbcolor} def'
  write(24,*) '/GF {gsave fill grestore} def'
  write(24,*) '/GG {0 setgray ST} def'
  write(24,*) '/GC {Colmesh ST} def'
  write(24,*) '/RF {setrgbcolor fill} def'
  write(24,*) '/SF {setgray fill} def'
  write(24,*) '/GS {gsave} def'
  write(24,*) '/GR {grestore} def'
  write(24,*) '/SLW {setlinewidth} def'
  write(24,*) '/SCSF {scalefont setfont} def'
  write(24,*) '% differents symboles utiles'
  write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
  write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/Cross {GS 0.05 CM SLW'
  write(24,*) 'GS 3 3 MR -6. -6. LR ST GR'
  write(24,*) 'GS 3 -3 MR -6. 6. LR ST GR'
  write(24,*) '0.01 CM SLW} def'
  write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
  write(24,*) '/Losange {GS 0.05 CM SLW 0 4.2 MR'
  write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
  write(24,*) 'GR 0.01 CM SLW} def'
  write(24,*) '%'
  write(24,*) '% niveaux de gris pour le modele de vitesse'
  write(24,*) '/BK {setgray fill} def'
  write(24,*) '% version noir et blanc'
  write(24,*) '%/BK {pop 1 setgray fill} def'
  write(24,*) '%'
  write(24,*) '% magenta pour les vecteurs deplacement'
  write(24,*) '/Colvects {0.01 CM SLW 1. 0. 1. RG} def'
  write(24,*) '% version noir et blanc'
  write(24,*) '%/Colvects {0.01 CM SLW 0. setgray} def'
  write(24,*) '%'
  write(24,*) '% chartreuse pour le maillage des macroblocs'
  write(24,*) '/Colmesh {0.02 CM SLW 0.5 1. 0. RG} def'
  write(24,*) '% version noir et blanc'
  write(24,*) '%/Colmesh {0.02 CM SLW 0. setgray} def'
  write(24,*) '%'
  write(24,*) '% cyan pour les sources et recepteurs'
  write(24,*) '/Colreceiv {0. 1. 1. RG} def'
  write(24,*) '% version noir et blanc'
  write(24,*) '%/Colreceiv {0. setgray} def'
  write(24,*) '%'
  write(24,*) '% macro dessin fleche'
  write(24,*) '/F {MV LR gsave LR ST grestore LR ST} def'
  write(24,*) '% macro dessin contour elements'
  write(24,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
  write(24,*) '%'
  write(24,*) '.01 CM SLW'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.35 CM SCSF'
  write(24,*) '%'
  write(24,*) '/vshift ',-height/2,' CM def'
  write(24,*) '/Rshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop neg vshift MR show } def'
  write(24,*) '/Cshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
  write(24,*) '/fN {/Helvetica-Bold findfont ',height,' CM SCSF} def'
  write(24,*) '%'
  write(24,*) 'gsave newpath 90 rotate'
  write(24,*) '0 ',-sizez,' CM translate 1. 1. scale'
  write(24,*) '%'

!
!--- ecriture des legendes du fichier PostScript
!
  write(24,*) '0 setgray'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.5 CM SCSF'

  if(legendes) then
  write(24,*) '24. CM 1.2 CM MV'
  write(24,610) usoffset,it
  write(24,*) '%'

  write(24,*) '24. CM 1.95 CM MV'
  timeval = it*dt
  if(timeval >= 1.d-3) then
    write(24,600) usoffset,timeval
  else
    write(24,601) usoffset,timeval
  endif
  write(24,*) '%'
  write(24,*) '24. CM 2.7 CM MV'
  write(24,640) usoffset,dispmax
  write(24,*) '%'
  write(24,*) '24. CM 3.45 CM MV'
  write(24,620) usoffset,cutvect*100.d0

  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.6 CM SCSF'
  if(colors == 1) write(24,*) '.4 .9 .9 setrgbcolor'
  write(24,*) '11 CM 1.1 CM MV'
  write(24,*) '(X axis) show'
  write(24,*) '%'
  write(24,*) '1.4 CM 9.5 CM MV'
  write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
  write(24,*) '(Y axis) show'
  write(24,*) 'grestore'
  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.7 CM SCSF'
  if(colors == 1) write(24,*) '.8 0 .8 setrgbcolor'
  write(24,*) '24.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  if(vecttype == 1) then
    write(24,*) '(Displacement vector field) show'
  else if(vecttype == 2) then
    write(24,*) '(Velocity vector field) show'
  else if(vecttype == 3) then
    write(24,*) '(Acceleration vector field) show'
  else
    stop 'Bad field code in PostScript display'
  endif
  write(24,*) 'grestore'
  write(24,*) '25.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(',stitle,') show'
  write(24,*) 'grestore'
  write(24,*) '26.45 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'

  if(ELASTIC) then
    write(24,*) '(Elastic Wave 2D - Spectral Element Method) show'
  else
    write(24,*) '(Acoustic Wave 2D - Spectral Element Method) show'
  endif

  write(24,*) 'grestore'

  endif

  write(24,*) '%'
  write(24,*) scalex,' ',scalez,' scale'
  write(24,*) '%'

!
!---- print the spectral elements mesh in PostScript
!

  print *,'Shape functions based on ',ngnod,' control nodes'

  convert = pi/180.d0

!
!----  draw the velocity model in background
!
  if(modelvect) then

  do ispec=1,nspec
    do i=1,NGLLX-subsamp,subsamp
          do j=1,NGLLX-subsamp,subsamp

  if((vpmax-vpmin)/vpmin > 0.02d0) then
  if(read_external_model) then
    x1 = (vpext(ibool(i,j,ispec))-vpmin)/ (vpmax-vpmin)
  else
    material = kmato(ispec)
    rlamda = elastcoef(1,material)
    rmu    = elastcoef(2,material)
    denst  = density(material)
    rKvol  = rlamda + 2.d0*rmu/3.d0
    cploc = sqrt((rKvol + 4.d0*rmu/3.d0)/denst)
    x1 = (cploc-vpmin)/(vpmax-vpmin)
  endif
  else
    x1 = 0.5d0
  endif

! rescaler pour eviter gris trop sombre
  x1 = x1*0.7 + 0.2
  if(x1 > 1.d0) x1=1.d0

! inverser echelle : blanc = vpmin, gris = vpmax
  x1 = 1.d0 - x1

  xw = coord(1,ibool(i,j,ispec))
  zw = coord(2,ibool(i,j,ispec))
  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
  write(24,500) xw,zw

  xw = coord(1,ibool(i+subsamp,j,ispec))
  zw = coord(2,ibool(i+subsamp,j,ispec))
  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
  write(24,499) xw,zw

  xw = coord(1,ibool(i+subsamp,j+subsamp,ispec))
  zw = coord(2,ibool(i+subsamp,j+subsamp,ispec))
  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
  write(24,499) xw,zw

  xw = coord(1,ibool(i,j+subsamp,ispec))
  zw = coord(2,ibool(i,j+subsamp,ispec))
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

  if(meshvect) then

  write(24,*) '%'
  write(24,*) '% spectral element mesh'
  write(24,*) '%'

  do ispec=1,nspec

  write(24,*) '% elem ',ispec

  do i=1,iptsdisp
  do j=1,iptsdisp
  xinterp(i,j) = 0.d0
  zinterp(i,j) = 0.d0
  do in = 1,ngnod
    nnum = knods(in,ispec)
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
  write(24,681) x1,z1

  if(ngnod == 4) then

! tracer des droites si elements 4 noeuds

  ir=iptsdisp
  x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  ir=iptsdisp
  is=iptsdisp
  x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  is=iptsdisp
  ir=1
  x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  ir=1
  is=2
  x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  else

! tracer des courbes si elements 9 noeuds
  do ir=2,iptsdisp
    x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
    z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim
    write(24,681) x2,z2
  enddo

  ir=iptsdisp
  do is=2,iptsdisp
    x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
    z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim
    write(24,681) x2,z2
  enddo

  is=iptsdisp
  do ir=iptsdisp-1,1,-1
    x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
    z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim
    write(24,681) x2,z2
  enddo

  ir=1
  do is=iptsdisp-1,2,-1
    x2 = (xinterp(ir,is)-xmin)*rapp_page + orig_x
    z2 = (zinterp(ir,is)-zmin)*rapp_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim
    write(24,681) x2,z2
  enddo

  endif

  write(24,*) 'CO'

  if(colors == 1) then

! For the moment 20 different colors max
  nbcols = 20

! Use a different color for each material set
  imat = kmato(ispec)
  icol = mod(imat - 1,nbcols) + 1

  write(24,680) red(icol),green(icol),blue(icol)

  endif

  if(modelvect) then
  write(24,*) 'GC'
  else
  write(24,*) 'GG'
  endif

! write the element number, the group number and the material number inside the element
  if(numbers == 1) then

  xw = (coorg(1,knods(1,ispec)) + coorg(1,knods(2,ispec)) + &
          coorg(1,knods(3,ispec)) + coorg(1,knods(4,ispec))) / 4.d0
  zw = (coorg(2,knods(1,ispec)) + coorg(2,knods(2,ispec)) + &
          coorg(2,knods(3,ispec)) + coorg(2,knods(4,ispec))) / 4.d0
  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
  if(colors == 1) write(24,*) '1 setgray'

  write(24,500) xw,zw

!--- ecriture numero de l'element
  write(24,502) ispec

  endif

  enddo

  endif

!
!----  draw the boundary conditions
!

  if(anyabs .and. boundvect) then

  write(24,*) '%'
  write(24,*) '% boundary conditions on the mesh'
  write(24,*) '%'

  write(24,*) '0.05 CM SLW'

!--- bords absorbants

  if(anyabs) then

  do ispecabs = 1,nelemabs
  ispec = numabs(ispecabs)

!--- une couleur pour chaque condition absorbante
!--- bord absorbant de type "haut"   : orange
!--- bord absorbant de type "bas"    : vert clair
!--- bord absorbant de type "gauche" : rose clair
!--- bord absorbant de type "droite" : turquoise

  do ibord = 1,4

  if(codeabs(ibord,ispecabs) /= 0) then

  if(ibord == ITOP) then
    write(24,*) '1. .85 0. RG'
    ideb = 3
    ifin = 4
  else if(ibord == IBOTTOM) then
    write(24,*) '.4 1. .4 RG'
    ideb = 1
    ifin = 2
  else if(ibord == ILEFT) then
    write(24,*) '1. .43 1. RG'
    ideb = 4
    ifin = 1
  else if(ibord == IRIGHT) then
    write(24,*) '.25 1. 1. RG'
    ideb = 2
    ifin = 3
  else
    stop 'Wrong absorbing boundary code'
  endif

  x1 = (coorg(1,knods(ideb,ispec))-xmin)*rapp_page + orig_x
  z1 = (coorg(2,knods(ideb,ispec))-zmin)*rapp_page + orig_z
  x2 = (coorg(1,knods(ifin,ispec))-xmin)*rapp_page + orig_x
  z2 = (coorg(2,knods(ifin,ispec))-zmin)*rapp_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,602) x1,z1,x2,z2

  endif
  enddo

  enddo

  endif

  write(24,*) '0 setgray'
  write(24,*) '0.01 CM SLW'

  endif

!
!----  draw the normalized displacement field
!

! return if the maximum displacement equals zero (no source)
  if(dispmax == 0.d0) then
    print *,' null displacement : returning !'
    return
  endif

  write(24,*) '%'
  write(24,*) '% vector field'
  write(24,*) '%'

! fleches en couleur si modele de vitesse en background
  if(modelvect) then
        write(24,*) 'Colvects'
  else
        write(24,*) '0 setgray'
  endif

  if(interpol) then

  print *,'Interpolating the vector field...'

  do ispec=1,nspec

! interpolation sur grille reguliere
  if(mod(ispec,1000) == 0) write(*,*) 'Interpolation uniform grid element ',ispec

  do i=1,iptsdisp
  do j=1,iptsdisp

  xinterp(i,j) = 0.d0
  zinterp(i,j) = 0.d0
  do in = 1,ngnod
    nnum = knods(in,ispec)
      xinterp(i,j) = xinterp(i,j) + shapeint(in,i,j)*coorg(1,nnum)
      zinterp(i,j) = zinterp(i,j) + shapeint(in,i,j)*coorg(2,nnum)
  enddo

  Uxinterp(i,j) = 0.d0
  Uzinterp(i,j) = 0.d0

  do k=1,NGLLX
  do l=1,NGLLX

  Uxinterp(i,j) = Uxinterp(i,j) + &
                displ(1,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
  Uzinterp(i,j) = Uzinterp(i,j) + &
                displ(2,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)

  enddo
  enddo

  x1 =(xinterp(i,j)-xmin)*rapp_page
  z1 =(zinterp(i,j)-zmin)*rapp_page

  x2 = Uxinterp(i,j)*sizemax/dispmax
  z2 = Uzinterp(i,j)*sizemax/dispmax

  d = sqrt(x2**2 + z2**2)

! ignorer si vecteur trop petit
  if(d > cutvect*sizemax) then

  d1 = d * rapport
  d2 = d1 * cos(angle*convert)

  dummy = x2/d
  if(dummy > 0.9999d0) dummy = 0.9999d0
  if(dummy < -0.9999d0) dummy = -0.9999d0
  theta = acos(dummy)

  if(z2 < 0.d0) theta = 360.d0*convert - theta
  thetaup = theta - angle*convert
  thetadown = theta + angle*convert

! tracer le vecteur proprement dit
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  x2 = x2 * centim
  z2 = z2 * centim
  xa = -d2*cos(thetaup)
  za = -d2*sin(thetaup)
  xa = xa * centim
  za = za * centim
  xb = -d2*cos(thetadown)
  zb = -d2*sin(thetadown)
  xb = xb * centim
  zb = zb * centim
  write(name,700) xb,zb,xa,za,x2,z2,x1,z1

! filtrer les blancs inutiles pour diminuer taille fichier PostScript
  longueur = 49
  indice = 1
  first = .false.
  do ii=1,longueur-1
    if(ch1(ii) /= ' ' .or. first) then
      if(ch1(ii) /= ' ' .or. ch1(ii+1) /= ' ') then
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

  d = sqrt(x2**2 + z2**2)

! ignorer si vecteur trop petit
  if(d > cutvect*sizemax) then

  d1 = d * rapport
  d2 = d1 * cos(angle*convert)

  dummy = x2/d
  if(dummy > 0.9999d0) dummy = 0.9999d0
  if(dummy < -0.9999d0) dummy = -0.9999d0
  theta = acos(dummy)

  if(z2 < 0.d0) theta = 360.d0*convert - theta
  thetaup = theta - angle*convert
  thetadown = theta + angle*convert

! tracer le vecteur proprement dit
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  x2 = x2 * centim
  z2 = z2 * centim
  xa = -d2*cos(thetaup)
  za = -d2*sin(thetaup)
  xa = xa * centim
  za = za * centim
  xb = -d2*cos(thetadown)
  zb = -d2*sin(thetadown)
  xb = xb * centim
  zb = zb * centim
  write(name,700) xb,zb,xa,za,x2,z2,x1,z1

! filtrer les blancs inutiles pour diminuer taille fichier PostScript
  longueur = 49
  indice = 1
  first = .false.
  do ii=1,longueur-1
    if(ch1(ii) /= ' ' .or. first) then
      if(ch1(ii) /= ' ' .or. ch1(ii+1) /= ' ') then
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

! sources et recepteurs en couleur si modele de vitesse
  if(modelvect) then
    write(24,*) 'Colreceiv'
  else
    write(24,*) '0 setgray'
  endif

!
!----  write position of the source
!
  xw = x_source
  zw = z_source
  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
  write(24,510) xw,zw
  if(isymbols) then
    write(24,*) 'Cross'
  else
    write(24,*) '(S) show'
  endif

!
!----  write position of the receivers
!
  do i=1,nrec
  if(i == 1) write(24,*) '% debut ligne recepteurs'
  if(i == nrec) write(24,*) '% fin ligne recepteurs'

  xw = st_xval(i)
  zw = st_zval(i)

  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
  write(24,510) xw,zw
  if(isymbols) then
    if(nrec > ndots .and. i /= 1 .and. i /= nrec) then
      write(24,*) 'VDot'
    else
      write(24,*) 'Losange'
    endif
  else
  write(24,*) '(R',i,') show'
  endif
  enddo

  write(24,*) '%'
  write(24,*) 'grestore'
  write(24,*) 'showpage'

  close(24)

 10   format('%!PS-Adobe-2.0',/,'%%',/,'%% Title: ',a50,/, &
          '%% Created by: Specfem Version 5.0',/, &
          '%% Author: Dimitri Komatitsch',/,'%%')
 510  format(f5.1,1x,f5.1,' M')
 600  format(f6.3,' neg CM 0 MR (Time =',f6.3,' s) show')
 601  format(f6.3,' neg CM 0 MR (Time =',1pe10.3,' s) show')
 610  format(f6.3,' neg CM 0 MR (Time step = ',i5,') show')
 620  format(f6.3,' neg CM 0 MR (Cut =',f5.2,' \%) show')
 640  format(f6.3,' neg CM 0 MR (Max norm =',1pe10.3,') show')

 200  format(80(a1))
 499  format(f5.1,1x,f5.1,' L')
 500  format(f5.1,1x,f5.1,' M')
 502  format('fN (',i4,') Cshow')
 680  format(f12.6,1x,f12.6,1x,f12.6,' RG GF')
 681  format(f6.2,1x,f6.2)
 602  format(f6.2,1x,f6.2,' M ',f6.2,1x,f6.2,' L ST')
 604  format('CP ',f12.6,' BK')
 700  format(8(f5.1,1x),'F')

  end subroutine plotpost

