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

  subroutine contol
!
!=======================================================================
!
!     "C o n t o l" :  Reads main control parameters
!      -----------
!
!=======================================================================
!

  use iounit
  use infos
  use mesh01
  use constspec
  use energie
  use verifs

  implicit none

  character(len=80) datlin

!
!-----------------------------------------------------------------------
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) ndofn,ndime,npgeo
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) display,ignuplot,interpol
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) itaff, itfirstaff, icolor, inumber
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) ivectplot,imeshvect,imodelvect,iboundvect,cutvect,isubsamp
  cutvect = cutvect / 100.d0
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) scalex,scalez,sizemax,angle,rapport,usletter
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) orig_x,orig_z,isymbols
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) valseuil,freqmaxrep
  valseuil = valseuil / 100.d0
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) sismos,nrec,nrec1,nrec2,isamp
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) irepr,anglerec,anglerec2
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) compenergy
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) initialfield,factorana,factorxsu,n1ana,n2ana
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) isismostype,ivecttype,iaffinfo
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) ireadmodel,ioutputgrid,iavs
!
!----
!
  read(iin , 40) datlin
  read(iin , * ) iexec,iecho
!
!---- output control parameters
!
  if(iecho  /=  0)  then
    write(iout,200) npgeo,ndofn,ndime
    write(iout,500) iexec,iecho
    write(iout,600) itaff,itfirstaff,iaffinfo,icolor,inumber
    write(iout,700) sismos,nrec,isamp,isismostype,nrec1,nrec2,anglerec, &
        anglerec2,compenergy,100.d0*valseuil,freqmaxrep
    write(iout,750) initialfield,ireadmodel,ioutputgrid,iavs
    write(iout,800) ivecttype,100.d0*cutvect,isubsamp,scalex,scalez, &
        sizemax,angle,rapport,orig_x,orig_z,usletter
  endif

  return

  40    format(a80)
  200   format(//1x,'C o n t r o l   c a r d   n o .  1',/1x,34('='),//5x,&
  'Number of spectral elements control nodes. . (npgeo) =',i8/5x, &
  'Number of d.o.f per node . . . . . . . . . . (ndofn) =',i8/5x, &
  'Number of space dimensions . . . . . . . . . (ndime) =',i8)
  500   format(//1x,'C o n t r o l   c a r d   n o .  2',/1x,34('='),//5x,&
  'Execution mode . . . . . . . . . . . . . . . (iexec) =',i5/ 5x, &
  '        ==  0     data check only                     ',  / 5x, &
  '        ==  1     resolution                          ',  / 5x, &
  'Data echoing . . . . . . . . . . . . . . . . (iecho) =',i5/ 5x, &
  '        ==  0     do not echo input data              ',  / 5x, &
  '        ==  1     echo input data - short listing     ',  / 5x, &
  '        ==  2     echo input data - full listing      ')
  600   format(//1x,'C o n t r o l   c a r d   n o .  3',/1x,34('='),//5x, &
  'Display frequency  . . . . . . . . . . . . . (itaff) = ',i5/ 5x, &
  'First display . . . . . . . . . . . . . (itfirstaff) = ',i5/ 5x, &
  'Basic info output frequency . . . . . . . (iaffinfo) = ',i5/ 5x, &
  'Color display . . . . . . . . . . . . . . . (icolor) = ',i5/ 5x, &
  '        ==  0     black and white display              ',  / 5x, &
  '        ==  1     color display                        ',  /5x, &
  'Numbered mesh . . . . . . . . . . . . . . .(inumber) = ',i5/ 5x, &
  '        ==  0     do not number the mesh               ',  /5x, &
  '        ==  1     number the mesh                      ')
  700   format(//1x,'C o n t r o l   c a r d   n o .  4',/1x,34('='),//5x, &
  'Record seismograms or not. . . . . . . . . .(sismos) = ',l6/5x, &
  'Total number of receivers. . . . . . . . . . .(nrec) = ',i6/5x, &
  'Subsampling for seismograms recording . . . .(isamp) = ',i6/5x, &
  'Seismograms recording type. . . . . . .(isismostype) = ',i6/5x, &
  'Number of receivers on first line . . . . . .(nrec1) = ',i6/5x, &
  'Number of receivers on second line. . . . . .(nrec2) = ',i6/5x, &
  'Angle for first line of receivers. . . . .(anglerec) = ',f6.2/5x, &
  'Angle for second line of receivers. . . .(anglerec2) = ',f6.2/5x, &
  'Compute total and potential energy . . .(compenergy) = ',l6/5x, &
  'Threshold for maximum frequency in % . . .(valseuil) = ',f6.2/5x, &
  'Maximal frequency plotted in spectrum. .(freqmaxrep) = ',1pe8.2)
  750   format(//1x,'C o n t r o l   c a r d   n o .  5',/1x,34('='),//5x, &
  'Read external initial field or not . .(initialfield) = ',l6/5x, &
  'Read external velocity model or not. . .(ireadmodel) = ',l6/5x, &
  'Save grid in external file or not . . .(ioutputgrid) = ',l6/5x, &
  'Save results in AVS file or not. . . . . . . .(iavs) = ',l6)
  800   format(//1x,'C o n t r o l   c a r d   n o .  6',/1x,34('='),//5x, &
  'Vector display type . . . . . . . . . . .(ivecttype) = ',i6/5x, &
  'Percentage of cut for vector plots. . . . .(cutvect) = ',f6.2/5x, &
  'Subsampling for velocity model display . .(isubsamp) = ',i6/5x, &
  'X-Scaling of plot for PostScript . . . . . .(scalex) = ',f6.2/5x, &
  'Z-Scaling of plot for PostScript . . . . . .(scalez) = ',f6.2/5x, &
  'Max size of arrows for PostScript . . . . .(sizemax) = ',f6.2/5x, &
  'Angle of vector arrows. . . . . . . . . . . .(angle) = ',f6.2/5x, &
  'Head to body ratio for arrows . . . . . . .(rapport) = ',f6.2/5x, &
  'X origin for Postscript display. . . . . . .(orig_x) = ',f6.2/5x, &
  'Z origin for Postscript display. . . . . . .(orig_z) = ',f6.2/5x, &
  'US letter format or French A4. . . . . . .(usletter) = ',l6)
!
!----
!
  end subroutine contol
