
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

  subroutine plotpost(U,coord,vpext,gltfl,posrec,nltfl,it,dt,coorg, &
          xinterp,zinterp,shapeint, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,codeperio,anyabs,anyperio)

!
! routine affichage postscript
!

  use palette
  use captio
  use timeparams
  use constspec
  use mesh01
  use spela202

  implicit none

!--- ecrire legendes ou non
  logical, parameter :: legendes=.true.

  integer kmato(nspec),knods(ngnod,nspec)
  integer ibool(nxgll,nygll,nspec)

  double precision xinterp(iptsdisp,iptsdisp),zinterp(iptsdisp,iptsdisp)
  double precision shapeint(ngnod,iptsdisp,iptsdisp)
  double precision Uxinterp(iptsdisp,iptsdisp)
  double precision Uzinterp(iptsdisp,iptsdisp)
  double precision flagrange(0:nxgll-1,iptsdisp)
  double precision density(numat),elastcoef(4,numat)

  integer nltfl,it
  double precision dt,timeval
  double precision U(ndime,npoin),coord(ndime,npoin)
  double precision vpext(npoin)

  double precision coorg(ndime,npgeo)
  double precision gltfl(20,nltfl)
  double precision posrec(ndime,nrec)

  integer numabs(nelemabs),codeabs(4,nelemabs)
  integer codeperio(4,nelemperio)
  logical anyabs,anyperio

! limite pour afficher des points a la place des recepteurs
  integer, parameter :: ndots = 10

! taille de la fenetre de display Postscript en pourcentage de la feuille
  double precision, parameter :: rpercentx = 70.0d0, rpercentz = 77.0d0

  double precision xmax,zmax,height,xw,zw,usoffset
  integer i,iglobrec,iglobsource
  character(len=40) name

! papier A4 ou US letter
  if(usletter) then
  usoffset = 1.75d0
  sizex = 27.94d0
  sizez = 21.59d0
  else
  usoffset = 0.
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
  rapp_page = min(rpercentz*sizez/(zmax-zmin),rpercentx*sizex/(xmax-xmin)) &
                        / 100.d0

! recherche de la valeur maximum de la norme du deplacement
  dispmax = maxval(sqrt(U(1,:)**2 + U(2,:)**2))
  write(*,*) 'Max norme = ',dispmax

! hauteur des numeros de domaine en CM
  height = 0.25d0

!
!---- ouverture du fichier PostScript
!
  write(name,222) it
  open(unit=24,file=name,status='unknown')
  222 format('vect',i5.5,'.ps')

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

  if (legendes) then
  write(24,*) '24. CM 1.2 CM MV'
  write(24,610) usoffset,it
  write(24,*) '%'

  write(24,*) '24. CM 1.95 CM MV'
  timeval = it*dt
  if(timeval  >=  1.d-3) then
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
  write(24,*) '24. CM 4.2 CM MV'
  write(24,630) usoffset,niter

  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.6 CM SCSF'
  if (icolor  ==  1) write(24,*) '.4 .9 .9 setrgbcolor'
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
  if (icolor  ==  1) write(24,*) '.8 0 .8 setrgbcolor'
  write(24,*) '24.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  if (ivecttype  ==  1) then
      write(24,*) '(Displacement vector field) show'
  else if (ivecttype  ==  2) then
      write(24,*) '(Velocity vector field) show'
  else if (ivecttype  ==  3) then
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
  write(24,*) '(Elastic Wave 2D - Spectral Element Method) show'
  write(24,*) 'grestore'

  endif

  write(24,*) '%'
  write(24,*) scalex,' ',scalez,' scale'
  write(24,*) '%'

!
!----   plot mesh and displacement vector field in a PostScript file
!
  call plotvect(knods,coorg,coord,U, &
          density,elastcoef,kmato,flagrange,xinterp,zinterp,shapeint, &
          Uxinterp,Uzinterp,ibool,vpext, &
          numabs,codeabs,codeperio,anyabs,anyperio)

! sources et recepteurs en couleur si modele de vitesse
  if(imodelvect) then
    write(24,*) 'Colreceiv'
  else
    write(24,*) '0 setgray'
  endif

!
!----  write position of the sources
!
  do i=1,nltfl

  iglobsource = nint(gltfl(9,i))

  xw = coord(1,iglobsource)
  zw = coord(2,iglobsource)
  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
  write(24,510) xw,zw
  if (isymbols) then
  write(24,*) 'Cross'
  else
  write(24,*) '(S',i,') show'
  endif
  enddo

!
!----  write position of the receivers
!
  do i=1,nrec
  if(i  ==  n1ana .or. i  ==  n2ana) write(24,*) '% solution analytique trace ',i
  if(i  ==  1) write(24,*) '% debut ligne recepteurs'
  if(i  ==  nrec) write(24,*) '% fin ligne recepteurs'

  iglobrec = nint(posrec(1,i))
  xw = coord(1,iglobrec)
  zw = coord(2,iglobrec)

  xw = (xw-xmin)*rapp_page + orig_x
  zw = (zw-zmin)*rapp_page + orig_z
  xw = xw * centim
  zw = zw * centim
  write(24,510) xw,zw
  if (isymbols) then
  if(nrec  >  ndots.and.i  /=  1.and.i  /=  nrec &
              .and.i  /=  n1ana.and.i  /=  n2ana) then
        if(i  >  nrec1) then
   write(24,*) 'HDot'
        else
   write(24,*) 'VDot'
        endif
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
          '%% Created by: Specfem Version 4.2',/, &
          '%% Author: Dimitri Komatitsch',/,'%%')
 510  format(f5.1,1x,f5.1,' M')
 600  format(f6.3,' neg CM 0 MR (Time =',f6.3,' s) show')
 601  format(f6.3,' neg CM 0 MR (Time =',1pe10.3,' s) show')
 610  format(f6.3,' neg CM 0 MR (Time step = ',i5,') show')
 620  format(f6.3,' neg CM 0 MR (Cut =',f5.2,' \%) show')
 630  format(f6.3,' neg CM 0 MR (Niter =',i2,') show')
 640  format(f6.3,' neg CM 0 MR (Max norm =',1pe10.3,') show')

  end subroutine plotpost
