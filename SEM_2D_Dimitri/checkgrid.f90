
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!                         Dimitri Komatitsch
!                     University of Pau, France
!
!                          (c) April 2007
!
!========================================================================

  subroutine checkgrid(vpext,vsext,rhoext,density,elastcoef,ibool,kmato,coord,npoin,vpmin,vpmax, &
                 assign_external_model,nspec,numat,deltat,f0,t0,initialfield,time_function_type, &
                 coorg,xinterp,zinterp,shapeint,knods,simulation_title,npgeo,pointsdisp,ngnod,any_elastic)

! check the mesh, stability and number of points per wavelength

  implicit none

  include "constants.h"

  integer i,j,ispec,material,npoin,nspec,numat,time_function_type

  integer, dimension(nspec) :: kmato
  integer, dimension(NGLLX,NGLLX,nspec) :: ibool

  double precision, dimension(numat) :: density
  double precision, dimension(4,numat) :: elastcoef
  double precision, dimension(NGLLX,NGLLX,nspec) :: vpext,vsext,rhoext

  double precision coord(NDIM,npoin)

  double precision vpmin,vpmax,vsmin,vsmax,densmin,densmax,vpmax_local,vpmin_local,vsmin_local
  double precision lambdaplus2mu,mu,denst,cploc,csloc
  double precision distance_min,distance_max,distance_min_local,distance_max_local
  double precision courant_stability_number_max,lambdaPmin,lambdaPmax,lambdaSmin,lambdaSmax
  double precision f0,t0,deltat,distance_1,distance_2,distance_3,distance_4

  logical assign_external_model,initialfield,any_elastic

! for the stability condition
! maximum polynomial degree for which we can compute the stability condition
  integer, parameter :: NGLLX_MAX_STABILITY = 15
  double precision :: percent_GLL(NGLLX_MAX_STABILITY)

  integer pointsdisp,npgeo,ngnod,is,ir,in,nnum

  double precision :: xmax,zmax,height,usoffset,sizex,sizez,courant_stability_number
  double precision :: x1,z1,x2,z2,ratio_page,xmin,zmin,lambdaS_local,lambdaP_local

  integer knods(ngnod,nspec)

  double precision xinterp(pointsdisp,pointsdisp),zinterp(pointsdisp,pointsdisp)
  double precision shapeint(ngnod,pointsdisp,pointsdisp)

  double precision coorg(NDIM,npgeo)

! title of the plot
  character(len=60) simulation_title

! define percentage of smallest distance between GLL points for NGLLX points
! percentages were computed by calling the GLL points routine for each degree
  percent_GLL(2) = 100.d0
  percent_GLL(3) = 50.d0
  percent_GLL(4) = 27.639320225002102d0
  percent_GLL(5) = 17.267316464601141d0
  percent_GLL(6) = 11.747233803526763d0
  percent_GLL(7) = 8.4888051860716516d0
  percent_GLL(8) = 6.4129925745196719d0
  percent_GLL(9) = 5.0121002294269914d0
  percent_GLL(10) = 4.0233045916770571d0
  percent_GLL(11) = 3.2999284795970416d0
  percent_GLL(12) = 2.7550363888558858d0
  percent_GLL(13) = 2.3345076678918053d0
  percent_GLL(14) = 2.0032477366369594d0
  percent_GLL(15) = 1.7377036748080721d0

! convert to real percentage
  percent_GLL(:) = percent_GLL(:) / 100.d0

  if(NGLLX > NGLLX_MAX_STABILITY) stop 'cannot estimate the stability condition for that degree'

!---- compute parameters for the spectral elements

  vpmin = HUGEVAL
  vsmin = HUGEVAL
  vpmax = -HUGEVAL
  vsmax = -HUGEVAL
  densmin = HUGEVAL
  densmax = -HUGEVAL

  distance_min = HUGEVAL
  distance_max = -HUGEVAL

  courant_stability_number_max = -HUGEVAL

  lambdaPmin = HUGEVAL
  lambdaSmin = HUGEVAL
  lambdaPmax = -HUGEVAL
  lambdaSmax = -HUGEVAL

  do ispec=1,nspec

    material = kmato(ispec)

    mu = elastcoef(2,material)
    lambdaplus2mu  = elastcoef(3,material)
    denst = density(material)

    cploc = sqrt(lambdaplus2mu/denst)
    csloc = sqrt(mu/denst)

  vpmax_local = -HUGEVAL
  vpmin_local = HUGEVAL
  vsmin_local = HUGEVAL

  distance_min_local = HUGEVAL
  distance_max_local = -HUGEVAL

  do j=1,NGLLZ
    do i=1,NGLLX

!--- if heterogeneous formulation with external velocity model
    if(assign_external_model) then
      cploc = vpext(i,j,ispec)
      csloc = vsext(i,j,ispec)
      denst = rhoext(i,j,ispec)
    endif

!--- compute min and max of velocity and density models
    vpmin = min(vpmin,cploc)
    vpmax = max(vpmax,cploc)

! ignore fluid regions with Vs = 0
    if(csloc > 0.0001d0) vsmin = min(vsmin,csloc)
    vsmax = max(vsmax,csloc)

    densmin = min(densmin,denst)
    densmax = max(densmax,denst)

    vpmax_local = max(vpmax_local,cploc)
    vpmin_local = min(vpmin_local,cploc)
    vsmin_local = min(vsmin_local,csloc)

    enddo
  enddo

! compute minimum and maximum size of edges of this grid cell
  distance_1 = sqrt((coord(1,ibool(1,1,ispec)) - coord(1,ibool(NGLLX,1,ispec)))**2 + &
               (coord(2,ibool(1,1,ispec)) - coord(2,ibool(NGLLX,1,ispec)))**2)

  distance_2 = sqrt((coord(1,ibool(NGLLX,1,ispec)) - coord(1,ibool(NGLLX,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,1,ispec)) - coord(2,ibool(NGLLX,NGLLZ,ispec)))**2)

  distance_3 = sqrt((coord(1,ibool(NGLLX,NGLLZ,ispec)) - coord(1,ibool(1,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,NGLLZ,ispec)) - coord(2,ibool(1,NGLLZ,ispec)))**2)

  distance_4 = sqrt((coord(1,ibool(1,NGLLZ,ispec)) - coord(1,ibool(1,1,ispec)))**2 + &
               (coord(2,ibool(1,NGLLZ,ispec)) - coord(2,ibool(1,1,ispec)))**2)

  distance_min_local = min(distance_1,distance_2,distance_3,distance_4)
  distance_max_local = max(distance_1,distance_2,distance_3,distance_4)

  distance_min = min(distance_min,distance_min_local)
  distance_max = max(distance_max,distance_max_local)

  courant_stability_number_max = max(courant_stability_number_max,vpmax_local * deltat / (distance_min_local * percent_GLL(NGLLX)))

! ignore fluid regions with Vs = 0
  if(csloc > 0.0001d0) then
    lambdaSmin = min(lambdaSmin,vsmin_local / (distance_max_local / (NGLLX - 1)))
    lambdaSmax = max(lambdaSmax,vsmin_local / (distance_max_local / (NGLLX - 1)))
  endif

  lambdaPmin = min(lambdaPmin,vpmin_local / (distance_max_local / (NGLLX - 1)))
  lambdaPmax = max(lambdaPmax,vpmin_local / (distance_max_local / (NGLLX - 1)))

  enddo

  write(IOUT,*)
  write(IOUT,*) '********'
  write(IOUT,*) 'Model: P velocity min,max = ',vpmin,vpmax
  write(IOUT,*) 'Model: S velocity min,max = ',vsmin,vsmax
  write(IOUT,*) 'Model: density min,max = ',densmin,densmax
  write(IOUT,*) '********'
  write(IOUT,*)

  write(IOUT,*)
  write(IOUT,*) '*********************************************'
  write(IOUT,*) '*** Verification of simulation parameters ***'
  write(IOUT,*) '*********************************************'
  write(IOUT,*)
  write(IOUT,*) '*** Max grid size = ',distance_max
  write(IOUT,*) '*** Min grid size = ',distance_min
  write(IOUT,*) '*** Max/min ratio = ',distance_max/distance_min
  write(IOUT,*)
  write(IOUT,*) '*** Max stability for P wave velocity = ',courant_stability_number_max
  write(IOUT,*)

! only if time source is not a Dirac or Heaviside (otherwise maximum frequency of spectrum undefined)
! and if source is not an initial field, for the same reason
  if(.not. initialfield .and. time_function_type /= 4 .and. time_function_type /= 5) then

    write(IOUT,*) ' Onset time = ',t0
    write(IOUT,*) ' Fundamental period = ',1.d0/f0
    write(IOUT,*) ' Fundamental frequency = ',f0
    if(t0 <= 1.d0/f0) then
      stop 'Onset time too small'
    else
      write(IOUT,*) ' --> onset time ok'
    endif
    write(IOUT,*) '----'
    write(IOUT,*) ' Nb pts / lambdaPmin_fmax max = ',lambdaPmax/(2.5d0*f0)
    write(IOUT,*) ' Nb pts / lambdaPmin_fmax min = ',lambdaPmin/(2.5d0*f0)
    write(IOUT,*) '----'
    write(IOUT,*) ' Nb pts / lambdaSmin_fmax max = ',lambdaSmax/(2.5d0*f0)
    write(IOUT,*) ' Nb pts / lambdaSmin_fmax min = ',lambdaSmin/(2.5d0*f0)
    write(IOUT,*) '----'

  endif

!
!--------------------------------------------------------------------------------
!

! A4 or US letter paper
  if(US_LETTER) then
    usoffset = 1.75d0
    sizex = 27.94d0
    sizez = 21.59d0
  else
    usoffset = 0.d0
    sizex = 29.7d0
    sizez = 21.d0
  endif

! height of domain numbers in centimeters
  height = 0.25d0

! get minimum and maximum values of mesh coordinates
  xmin = minval(coord(1,:))
  zmin = minval(coord(2,:))
  xmax = maxval(coord(1,:))
  zmax = maxval(coord(2,:))

! ratio of physical page size/size of the domain meshed
  ratio_page = min(rpercentz*sizez/(zmax-zmin),rpercentx*sizex/(xmax-xmin)) / 100.d0

  print *
  print *,'Creating PostScript file with stability condition'

!
!---- open PostScript file
!
  open(unit=24,file='OUTPUT_FILES/mesh_stability.ps',status='unknown')

!
!---- write PostScript header
!
  write(24,10) simulation_title
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
  write(24,*) '% different useful symbols'
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
  write(24,*) '/Diamond {GS 0.05 CM SLW 0 4.2 MR'
  write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
  write(24,*) 'GR 0.01 CM SLW} def'
  write(24,*) '%'
  write(24,*) '% macro to draw the contour of the elements'
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
!--- write captions of PostScript figure
!
  write(24,*) '0 setgray'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.5 CM SCSF'

  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.6 CM SCSF'
  write(24,*) '.4 .9 .9 setrgbcolor'
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
  write(24,*) '.8 0 .8 setrgbcolor'
  write(24,*) '24.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(Mesh stability condition \(red = bad\)) show'
  write(24,*) 'grestore'
  write(24,*) '25.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(',simulation_title,') show'
  write(24,*) 'grestore'
  write(24,*) '26.45 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(2D Spectral Element Method) show'
  write(24,*) 'grestore'

  write(24,*) '%'
  write(24,*) '1 1 scale'
  write(24,*) '%'

!
!---- draw the spectral element mesh
!
  write(24,*) '%'
  write(24,*) '% spectral element mesh'
  write(24,*) '%'
  write(24,*) '0 setgray'

  do ispec = 1, nspec

  write(24,*) '% elem ',ispec

  do i=1,pointsdisp
  do j=1,pointsdisp
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
  x1 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z1 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  write(24,*) 'MK'
  write(24,681) x1,z1

! draw straight lines if elements have 4 nodes

  ir=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  ir=pointsdisp
  is=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  is=pointsdisp
  ir=1
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  ir=1
  is=2
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  write(24,*) 'CO'

    material = kmato(ispec)

    mu = elastcoef(2,material)
    lambdaplus2mu  = elastcoef(3,material)
    denst = density(material)

    cploc = sqrt(lambdaplus2mu/denst)
    csloc = sqrt(mu/denst)

  vpmax_local = -HUGEVAL
  vpmin_local = HUGEVAL
  vsmin_local = HUGEVAL

  distance_min_local = HUGEVAL
  distance_max_local = -HUGEVAL

  do j=1,NGLLZ
    do i=1,NGLLX

!--- if heterogeneous formulation with external velocity model
    if(assign_external_model) then
      cploc = vpext(i,j,ispec)
      csloc = vsext(i,j,ispec)
      denst = rhoext(i,j,ispec)
    endif

    vpmax_local = max(vpmax_local,cploc)
    vpmin_local = min(vpmin_local,cploc)
    vsmin_local = min(vsmin_local,csloc)

    enddo
  enddo

! compute minimum and maximum size of edges of this grid cell
  distance_1 = sqrt((coord(1,ibool(1,1,ispec)) - coord(1,ibool(NGLLX,1,ispec)))**2 + &
               (coord(2,ibool(1,1,ispec)) - coord(2,ibool(NGLLX,1,ispec)))**2)

  distance_2 = sqrt((coord(1,ibool(NGLLX,1,ispec)) - coord(1,ibool(NGLLX,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,1,ispec)) - coord(2,ibool(NGLLX,NGLLZ,ispec)))**2)

  distance_3 = sqrt((coord(1,ibool(NGLLX,NGLLZ,ispec)) - coord(1,ibool(1,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,NGLLZ,ispec)) - coord(2,ibool(1,NGLLZ,ispec)))**2)

  distance_4 = sqrt((coord(1,ibool(1,NGLLZ,ispec)) - coord(1,ibool(1,1,ispec)))**2 + &
               (coord(2,ibool(1,NGLLZ,ispec)) - coord(2,ibool(1,1,ispec)))**2)

  distance_min_local = min(distance_1,distance_2,distance_3,distance_4)
  distance_max_local = max(distance_1,distance_2,distance_3,distance_4)

  distance_min = min(distance_min,distance_min_local)
  distance_max = max(distance_max,distance_max_local)

  courant_stability_number = vpmax_local * deltat / (distance_min_local * percent_GLL(NGLLX))

! display bad elements that are above 80% of the threshold
  if(courant_stability_number >= 0.80 * courant_stability_number_max) then
    write(24,*) '1 0 0 RG GF GG'
  else
! do not color the elements if below the threshold
    write(24,*) 'ST'
  endif

  enddo ! end of loop on all the spectral elements

  write(24,*) '%'
  write(24,*) 'grestore'
  write(24,*) 'showpage'

  close(24)

  print *,'End of creation of PostScript file with stability condition'

!
!--------------------------------------------------------------------------------
!

  print *
  print *,'Creating PostScript file with mesh dispersion'

!
!---- open PostScript file
!
  if(any_elastic) then
    open(unit=24,file='OUTPUT_FILES/mesh_S_wave_dispersion.ps',status='unknown')
  else
    open(unit=24,file='OUTPUT_FILES/mesh_P_wave_dispersion.ps',status='unknown')
  endif

!
!---- write PostScript header
!
  write(24,10) simulation_title
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
  write(24,*) '% different useful symbols'
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
  write(24,*) '/Diamond {GS 0.05 CM SLW 0 4.2 MR'
  write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
  write(24,*) 'GR 0.01 CM SLW} def'
  write(24,*) '%'
  write(24,*) '% macro to draw the contour of the elements'
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
!--- write captions of PostScript figure
!
  write(24,*) '0 setgray'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.5 CM SCSF'

  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.6 CM SCSF'
  write(24,*) '.4 .9 .9 setrgbcolor'
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
  write(24,*) '.8 0 .8 setrgbcolor'
  write(24,*) '24.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  if(any_elastic) then
    write(24,*) '(Mesh elastic S-wave dispersion \(red = good, blue = bad\)) show'
  else
    write(24,*) '(Mesh acoustic P-wave dispersion \(red = good, blue = bad\)) show'
  endif
  write(24,*) 'grestore'
  write(24,*) '25.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(',simulation_title,') show'
  write(24,*) 'grestore'
  write(24,*) '26.45 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(2D Spectral Element Method) show'
  write(24,*) 'grestore'

  write(24,*) '%'
  write(24,*) '1 1 scale'
  write(24,*) '%'

!
!---- draw the spectral element mesh
!
  write(24,*) '%'
  write(24,*) '% spectral element mesh'
  write(24,*) '%'
  write(24,*) '0 setgray'

  do ispec = 1, nspec

  write(24,*) '% elem ',ispec

  do i=1,pointsdisp
  do j=1,pointsdisp
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
  x1 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z1 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  write(24,*) 'MK'
  write(24,681) x1,z1

! draw straight lines if elements have 4 nodes

  ir=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  ir=pointsdisp
  is=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  is=pointsdisp
  ir=1
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  ir=1
  is=2
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  write(24,*) 'CO'

    material = kmato(ispec)

    mu = elastcoef(2,material)
    lambdaplus2mu  = elastcoef(3,material)
    denst = density(material)

    cploc = sqrt(lambdaplus2mu/denst)
    csloc = sqrt(mu/denst)

  vpmax_local = -HUGEVAL
  vpmin_local = HUGEVAL
  vsmin_local = HUGEVAL

  distance_min_local = HUGEVAL
  distance_max_local = -HUGEVAL

  do j=1,NGLLZ
    do i=1,NGLLX

!--- if heterogeneous formulation with external velocity model
    if(assign_external_model) then
      cploc = vpext(i,j,ispec)
      csloc = vsext(i,j,ispec)
      denst = rhoext(i,j,ispec)
    endif

    vpmax_local = max(vpmax_local,cploc)
    vpmin_local = min(vpmin_local,cploc)
    vsmin_local = min(vsmin_local,csloc)

    enddo
  enddo

! compute minimum and maximum size of edges of this grid cell
  distance_1 = sqrt((coord(1,ibool(1,1,ispec)) - coord(1,ibool(NGLLX,1,ispec)))**2 + &
               (coord(2,ibool(1,1,ispec)) - coord(2,ibool(NGLLX,1,ispec)))**2)

  distance_2 = sqrt((coord(1,ibool(NGLLX,1,ispec)) - coord(1,ibool(NGLLX,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,1,ispec)) - coord(2,ibool(NGLLX,NGLLZ,ispec)))**2)

  distance_3 = sqrt((coord(1,ibool(NGLLX,NGLLZ,ispec)) - coord(1,ibool(1,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,NGLLZ,ispec)) - coord(2,ibool(1,NGLLZ,ispec)))**2)

  distance_4 = sqrt((coord(1,ibool(1,NGLLZ,ispec)) - coord(1,ibool(1,1,ispec)))**2 + &
               (coord(2,ibool(1,NGLLZ,ispec)) - coord(2,ibool(1,1,ispec)))**2)

  distance_min_local = min(distance_1,distance_2,distance_3,distance_4)
  distance_max_local = max(distance_1,distance_2,distance_3,distance_4)

  distance_min = min(distance_min,distance_min_local)
  distance_max = max(distance_max,distance_max_local)

! display mesh dispersion for S waves if there is at least one elastic element in the mesh
  if(any_elastic) then

! ignore fluid regions with Vs = 0
  if(csloc > 0.0001d0) then

    lambdaS_local = vsmin_local / (distance_max_local / (NGLLX - 1))

! display very good elements that are above 80% of the threshold in red
    if(lambdaS_local >= 0.80 * lambdaSmax) then
      write(24,*) '1 0 0 RG GF GG'

! display bad elements that are below 120% of the threshold in blue
    else if(lambdaS_local <= 1.20 * lambdaSmin) then
      write(24,*) '0 0 1 RG GF GG'

    else
! do not color the elements if not close to the threshold
      write(24,*) 'ST'
    endif

  else
! do not color the elements if S-wave velocity undefined
    write(24,*) 'ST'
  endif

! display mesh dispersion for P waves if there is no elastic element in the mesh
  else

    lambdaP_local = vpmin_local / (distance_max_local / (NGLLX - 1))

! display very good elements that are above 80% of the threshold in red
    if(lambdaP_local >= 0.80 * lambdaPmax) then
      write(24,*) '1 0 0 RG GF GG'

! display bad elements that are below 120% of the threshold in blue
    else if(lambdaP_local <= 1.20 * lambdaPmin) then
      write(24,*) '0 0 1 RG GF GG'

    else
! do not color the elements if not close to the threshold
      write(24,*) 'ST'
    endif

  endif

  enddo ! end of loop on all the spectral elements

  write(24,*) '%'
  write(24,*) 'grestore'
  write(24,*) 'showpage'

  close(24)

  print *,'End of creation of PostScript file with mesh dispersion'

!
!--------------------------------------------------------------------------------
!

  print *
  print *,'Creating PostScript file with velocity model'

!
!---- open PostScript file
!
  open(unit=24,file='OUTPUT_FILES/P_velocity_model.ps',status='unknown')

!
!---- write PostScript header
!
  write(24,10) simulation_title
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
  write(24,*) '% different useful symbols'
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
  write(24,*) '/Diamond {GS 0.05 CM SLW 0 4.2 MR'
  write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
  write(24,*) 'GR 0.01 CM SLW} def'
  write(24,*) '%'
  write(24,*) '% macro to draw the contour of the elements'
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
!--- write captions of PostScript figure
!
  write(24,*) '0 setgray'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.5 CM SCSF'

  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.6 CM SCSF'
  write(24,*) '.4 .9 .9 setrgbcolor'
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
  write(24,*) '.8 0 .8 setrgbcolor'
  write(24,*) '24.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(P-velocity model \(dark = fast, light = slow\)) show'
  write(24,*) 'grestore'
  write(24,*) '25.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(',simulation_title,') show'
  write(24,*) 'grestore'
  write(24,*) '26.45 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(2D Spectral Element Method) show'
  write(24,*) 'grestore'

  write(24,*) '%'
  write(24,*) '1 1 scale'
  write(24,*) '%'

!
!---- draw the spectral element mesh
!
  write(24,*) '%'
  write(24,*) '% spectral element mesh'
  write(24,*) '%'
  write(24,*) '0 setgray'

  do ispec = 1, nspec

  write(24,*) '% elem ',ispec

  do i=1,pointsdisp
  do j=1,pointsdisp
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
  x1 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z1 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  write(24,*) 'MK'
  write(24,681) x1,z1

! draw straight lines if elements have 4 nodes

  ir=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  ir=pointsdisp
  is=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  is=pointsdisp
  ir=1
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  ir=1
  is=2
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  write(24,681) x2,z2

  write(24,*) 'CO'

  if((vpmax-vpmin)/vpmin > 0.02d0) then
  if(assign_external_model) then
! use lower-left corner
    x1 = (vpext(1,1,ispec)-vpmin) / (vpmax-vpmin)
  else
    material = kmato(ispec)
    mu = elastcoef(2,material)
    lambdaplus2mu  = elastcoef(3,material)
    denst = density(material)
    cploc = sqrt(lambdaplus2mu/denst)
    x1 = (cploc-vpmin)/(vpmax-vpmin)
  endif
  else
    x1 = 0.5d0
  endif

! rescale to avoid very dark gray levels
  x1 = x1*0.7 + 0.2
  if(x1 > 1.d0) x1=1.d0

! invert scale: white = vpmin, dark gray = vpmax
  x1 = 1.d0 - x1

! display P-velocity model using gray levels
      write(24,*) sngl(x1),' setgray GF GG'

  enddo ! end of loop on all the spectral elements

  write(24,*) '%'
  write(24,*) 'grestore'
  write(24,*) 'showpage'

  close(24)

  print *,'End of creation of PostScript file with velocity model'

 10  format('%!PS-Adobe-2.0',/,'%%',/,'%% Title: ',a50,/,'%% Created by: Specfem2D',/,'%% Author: Dimitri Komatitsch',/,'%%')

 681 format(f6.2,1x,f6.2)

  end subroutine checkgrid

