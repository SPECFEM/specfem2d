
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

  subroutine writeseis(sisux,sisuz,coord,posrec,ndime, &
         npoin,nseis,nrec,isamp,deltat,factorxsu, &
          n1ana,n2ana,irepr,nrec1,nrec2,isismostype)

!
!----  sauvegarde des sismogrammes en fin de simulation
!

  implicit none

  integer ndime,npoin,nseis
  integer nrec,isamp,n1ana,n2ana,irepr,nrec1,nrec2,isismostype
  double precision deltat,factorxsu

! simple precision pour le stockage au format SEP
  real sisux(nseis,nrec)
  real sisuz(nseis,nrec)

  double precision coord(ndime,npoin)
  double precision posrec(ndime,nrec)

  logical invert
  integer nt,irec,i,iana,it
  double precision xval,xvaladd

  write(*,*) 'Sauvegarde sismos sur disk ...'

  nt = nseis

  write(*,*)
  write(*,*) ' valeur de isamp = ',isamp
  write(*,*) ' nb d''echantillons stockes en temps = ',nt
  write(*,*) ' nb de recepteurs = ',nrec
  write(*,*)

  write(*,*) 'Sauvegarde sismos sur disk ...'

!----

  write(*,*)
  write(*,*) 'Sauvegarde traces format SEP...'
  write(*,*) 'DK DK using ASCII instead of SEP'

  goto 333

! ecriture au format binaire deplacement horizontal
  open(unit=11,file='Ux_file',status='unknown', &
                        access='direct',recl=nt*nrec*4)
  write(11,rec=1) (sisux(i,1),i=1,nt*nrec)
  close(11)

! ecriture au format binaire deplacement vertical
  open(unit=11,file='Uz_file',status='unknown', &
                        access='direct',recl=nt*nrec*4)
  write(11,rec=1) (sisuz(i,1),i=1,nt*nrec)
  close(11)

 333  continue

! ecriture au format ASCII
  open(unit=11,file='Ux_file.dat',status='unknown')
!!!!!!! DK DK UUUUUUU only one receiver for tests  do irec=1,nrec
  do irec=1,1
    do it=1,nt
      write(11,*) sngl(dble(it-1)*dble(isamp)*deltat),sisux(it,irec)
    enddo
  enddo
  close(11)
  open(unit=11,file='Uz_file.dat',status='unknown')
!!!!!!! DK DK UUUUUUU only one receiver for tests  do irec=1,nrec
  do irec=1,1
    do it=1,nt
      write(11,*) sngl(dble(it-1)*dble(isamp)*deltat),sisuz(it,irec)
    enddo
  enddo
  close(11)

!----

  write(*,*)
  write(*,*) 'Sauvegarde headers pour visu...'

!----
!---- ligne de recepteurs pour Xwindow
!----

  open(unit=12,file='xline',status='unknown')

  write(12,100) factorxsu,nseis,deltat*isamp,nrec
! inverser representation si recepteurs orientes negativement
  invert = .false.
  if(irepr  ==  1.and.coord(1,nint(posrec(1,nrec)))  <  &
          coord(1,nint(posrec(1,1)))) then
    invert = .true.
  endif
  if(irepr  ==  2.and.coord(2,nint(posrec(1,nrec)))  <  &
          coord(2,nint(posrec(1,1)))) then
    invert = .true.
  endif

!--- premiere partie de la ligne de recepteurs
  do irec=1,nrec1
! recepteurs en distance
    if(irepr  ==  3.or.nrec2  >  0) then
          xval = dsqrt((coord(1,nint(posrec(1,irec))) - &
                coord(1,nint(posrec(1,1))))**2 + &
                (coord(2,nint(posrec(1,irec))) - &
                coord(2,nint(posrec(1,1))))**2)
! recepteurs suivant coordonnee X
    else if(irepr  ==  1) then
      if(invert) then
          xval = coord(1,nint(posrec(1,1))) - coord(1,nint(posrec(1,irec)))
      else
          xval = coord(1,nint(posrec(1,irec)))
      endif
! recepteurs suivant coordonnee Z
    else if(irepr  ==  2) then
      if(invert) then
          xval = coord(2,nint(posrec(1,1))) - coord(2,nint(posrec(1,irec)))
      else
          xval = coord(2,nint(posrec(1,irec)))
      endif
    else
          stop 'wrong value of irepr !'
    endif

    write(12,140) xval

    if (irec  <  nrec1) write(12,*) ','
  enddo

!--- deuxieme partie de la ligne de recepteurs
  if(nrec2  >  0) then
  write(12,*) ','
  xvaladd = xval
  do irec=nrec1+1,nrec
          xval = &
    dsqrt((coord(1,nint(posrec(1,irec))) - coord(1,nint(posrec(1,nrec1))))**2 + &
        (coord(2,nint(posrec(1,irec))) - coord(2,nint(posrec(1,nrec1))))**2)
    write(12,140) xval + xvaladd
    if (irec  <  nrec) write(12,*) ','
  enddo
  endif

  if(isismostype  ==  1) then
      write(12,*) '@title="Ux@displacement@component"@<@Ux_file'
  else if(isismostype  ==  2) then
      write(12,*) '@title="Ux@velocity@component"@<@Ux_file'
  else
      write(12,*) '@title="Ux@acceleration@component"@<@Ux_file'
  endif

  close(12)

!----
!---- script de visualisation
!----

  open(unit=12,file='showseis',status='unknown')
  write(12,110)
  write(12,*)
  write(12,*) '/bin/rm -f tempfile psline'
  write(12,*) '# concatener toutes les lignes'
  write(12,*) 'tr -d ''\012'' <xline >tempfile'
  write(12,*) '# remettre fin de ligne'
  write(12,*) 'echo " " >> tempfile'
  write(12,*) '# supprimer espaces, changer arobas, dupliquer'
  write(12,137)
  write(12,*) '/bin/rm -f tempfile'
  write(12,*) '# copier fichier pour sortie postscript'
  write(12,130)
  write(12,*) '/bin/rm -f tempfile'
  write(12,*) 'echo ''rm -f uxpoly.ps uzpoly.ps'' > tempfile'
  write(12,*) 'cat tempfile psline > tempfile2'
  write(12,*) '/bin/mv -f tempfile2 psline'
  write(12,*) '/bin/rm -f tempfile'
  write(12,*) '# executer commande xsu'
  write(12,*) 'sh xline'
  write(12,*) '/bin/rm -f tempfile tempfile2'
  close(12)

!----
!---- une trace pour Xwindow
!----

  open(unit=12,file='xtrace',status='unknown')
  write(12,110)
  write(12,*)
  write(12,*) 'set nt = ',nseis
  write(12,*) 'set dt = ',sngl(deltat*isamp)
  write(12,*) '@ trace=10000'
  write(12,*) 'while ($trace > -1)'
  write(12,*) 'echo Donnez le numero de trace a visualiser :'
  write(12,*) 'set rep=$<'
  write(12,*) '@ trace = $rep'
  write(12,*) 'echo Trace demandee : $trace'
  write(12,*) '# traces commencent a zero dans format SEP'
  write(12,*) '@ septrace = $trace - 1'
  if(isismostype  ==  1) then
      write(12,120)
      write(12,125)
  else if(isismostype  ==  2) then
      write(12,121)
      write(12,126)
  else
      write(12,122)
      write(12,127)
  endif
  write(12,*) 'end'
  close(12)

!----
!---- une trace pour postscript
!----

  open(unit=12,file='pstrace',status='unknown')
  write(12,110)
  write(12,*)
  write(12,*) 'set nt = ',nseis
  write(12,*) 'set dt = ',sngl(deltat*isamp)
  write(12,*) '@ trace=10000'
  write(12,*) 'while ($trace > -1)'
  write(12,*) 'echo Donnez le numero de trace a tranformer en postscript :'
  write(12,*) 'set rep=$<'
  write(12,*) '@ trace = $rep'
  write(12,*) 'echo Trace demandee : $trace'
  write(12,*) '# traces commencent a zero dans format SEP'
  write(12,*) '@ septrace = $trace - 1'
  write(12,*) 'rm -f uxtrace{$trace}.ps uztrace{$trace}.ps'
  if(isismostype  ==  1) then
      write(12,220)
      write(12,225)
  else if(isismostype  ==  2) then
      write(12,221)
      write(12,226)
  else
      write(12,222)
      write(12,227)
  endif
  write(12,*) 'end'
  close(12)

!----
!---- une trace avec comparaison analytique pour Xwindow
!----

  open(unit=12,file='xcomptrace',status='unknown')
  write(12,110)
  write(12,*)
  write(12,*) 'set nt = ',nseis
  write(12,*) 'set dt = ',sngl(deltat*isamp)
  write(12,*) 'set traceana1 = ',n1ana
  write(12,*) 'set traceana2 = ',n2ana
  write(12,*) '# traces commencent a zero dans format SEP'
  write(12,*) '@ septraceana1 = $traceana1 - 1'
  write(12,*) '@ septraceana2 = $traceana2 - 1'
  write(12,*) '# premiere trace analytique'
  write(12,*) '@ septraceref = 0'
  write(12,*) '/bin/rm -f tutuan tutucomp'
  write(12,320) 1,'x','x'
  write(12,330) 'x',1
  write(12,*) '/bin/rm -f tutuan tutucomp'
  write(12,320) 1,'z','z'
  write(12,330) 'z',1
  write(12,*) '# deuxieme trace analytique'
  write(12,*) '@ septraceref = 1'
  write(12,*) '/bin/rm -f tutuan tutucomp'
  write(12,320) 2,'x','x'
  write(12,330) 'x',2
  write(12,*) '/bin/rm -f tutuan tutucomp'
  write(12,320) 2,'z','z'
  write(12,330) 'z',2
  write(12,*) '/bin/rm -f tutuan tutucomp'
  close(12)

!----
!---- une trace avec comparaison analytique pour postscript
!----

  open(unit=12,file='pscomptrace',status='unknown')
  write(12,110)
  write(12,*)
  write(12,*) 'set nt = ',nseis
  write(12,*) 'set dt = ',sngl(deltat*isamp)
  write(12,*) 'set traceana1 = ',n1ana
  write(12,*) 'set traceana2 = ',n2ana
  write(12,*) '# traces commencent a zero dans format SEP'
  write(12,*) '@ septraceana1 = $traceana1 - 1'
  write(12,*) '@ septraceana2 = $traceana2 - 1'
  write(12,*) 'echo Generating PostScript files...'
  write(12,*) '/bin/rm -f uxtracecompana1.ps uztracecompana1.ps'
  write(12,*) '/bin/rm -f uxtracecompana2.ps uztracecompana2.ps'
  write(12,*) '# premiere trace analytique'
  write(12,*) '@ septraceref = 0'
  write(12,*) '/bin/rm -f tutuan tutucomp'
  write(12,320) 1,'x','x'
  write(12,340) 'x',1,'x',1
  write(12,*) '/bin/rm -f tutuan tutucomp'
  write(12,320) 1,'z','z'
  write(12,340) 'z',1,'z',1
  write(12,*) '# deuxieme trace analytique'
  write(12,*) '@ septraceref = 1'
  write(12,*) '/bin/rm -f tutuan tutucomp'
  write(12,320) 2,'x','x'
  write(12,340) 'x',2,'x',2
  write(12,*) '/bin/rm -f tutuan tutucomp'
  write(12,320) 2,'z','z'
  write(12,340) 'z',2,'z',2
  write(12,*) '/bin/rm -f tutuan tutucomp'
  close(12)

!----
!---- residus trace analytique pour Xwindow
!----

  open(unit=12,file='xresid',status='unknown')
  write(12,110)
  write(12,*)
  write(12,*) 'set nt = ',nseis
  write(12,*) 'set dt = ',sngl(deltat*isamp)
  write(12,*) '@ trace=',n1ana
  write(12,*) '@ septrace = $trace - 1'
  iana = 0
  write(12,170)
  write(12,150) iana,iana
  write(12,160)
  write(12,*) '@ trace=',n2ana
  write(12,*) '@ septrace = $trace - 1'
  iana = 1
  write(12,170)
  write(12,150) iana,iana
  write(12,160)
  write(12,170)
  close(12)

!----
!---- residus trace analytique pour PostScript (utilise Gnuplot)
!----

! facteur d'amplification des residus
  open(unit=12,file='psresid',status='unknown')
  write(12,110)
  write(12,*)
  write(12,*) 'set ampli = 5'
  write(12,*) 'set nt = ',nseis
  write(12,*) 'set dt = ',sngl(deltat*isamp)
  write(12,200)
  write(12,*) '@ trace=',n1ana
  write(12,*) '@ septrace = $trace - 1'
  iana = 0
  write(12,170)
  write(12,171)
  write(12,151) iana,iana
  write(12,152)
  write(12,154)
  write(12,155)
  write(12,*) '@ trace=',n2ana
  write(12,*) '@ septrace = $trace - 1'
  iana = 1
  write(12,170)
  write(12,171)
  write(12,151) iana,iana
  write(12,152)
  write(12,154)
  write(12,155)
  write(12,170)
  write(12,171)
  close(12)

  100  format('xwigb@xcur=',f8.2,'@n1=',i5, &
     '@d1=',f15.8,'@label1="Time@(s)"@label2="x@(m)"@n2=', &
     i5,'@x2=')
  110  format('#!/bin/csh -f')
  120  format('subset < Ux_file n1=$nt', &
     ' if2s=$septrace n2s=1 | xgraph -geometry 1085x272', &
     ' label1="Time (s)" label2="Amplitude (m)" ', &
     ' n=$nt style=normal d1=$dt', &
     ' title="Ux displacement component trace "$trace &')
  121  format('subset < Ux_file n1=$nt', &
     ' if2s=$septrace n2s=1 | xgraph -geometry 1085x272', &
     ' label1="Time (s)" label2="Amplitude (m)" ', &
     ' n=$nt style=normal d1=$dt', &
     ' title="Ux velocity component trace "$trace &')
  122  format('subset < Ux_file n1=$nt', &
     ' if2s=$septrace n2s=1 | xgraph -geometry 1085x272', &
     ' label1="Time (s)" label2="Amplitude (m)" ', &
     ' n=$nt style=normal d1=$dt', &
     ' title="Ux acceleration component trace "$trace &')
  125  format('subset < Uz_file n1=$nt', &
     ' if2s=$septrace n2s=1 | xgraph -geometry 1085x272', &
     ' label1="Time (s)" label2="Amplitude (m)" ', &
     ' n=$nt style=normal d1=$dt', &
     ' title="Uz displacement component trace "$trace &')
  126  format('subset < Uz_file n1=$nt', &
     ' if2s=$septrace n2s=1 | xgraph -geometry 1085x272', &
     ' label1="Time (s)" label2="Amplitude (m)" ', &
     ' n=$nt style=normal d1=$dt', &
     ' title="Uz velocity component trace "$trace &')
  127  format('subset < Uz_file n1=$nt', &
     ' if2s=$septrace n2s=1 | xgraph -geometry 1085x272', &
     ' label1="Time (s)" label2="Amplitude (m)" ', &
     ' n=$nt style=normal d1=$dt', &
     ' title="Uz acceleration component trace "$trace &')
  220  format('subset < Ux_file n1=$nt', &
     ' if2s=$septrace n2s=1 | psgraph ', &
     ' label1="Time (s)" label2="Amplitude (m)" ', &
     ' n=$nt style=normal d1=$dt hbox=4.0 wbox=6.0', &
     ' title="Ux displacement component trace "$trace > uxtrace{$trace}.ps')
  221  format('subset < Ux_file n1=$nt', &
     ' if2s=$septrace n2s=1 | psgraph ', &
     ' label1="Time (s)" label2="Amplitude (m)" ', &
     ' n=$nt style=normal d1=$dt hbox=4.0 wbox=6.0', &
     ' title="Ux velocity component trace "$trace > uxtrace{$trace}.ps')
  222  format('subset < Ux_file n1=$nt', &
     ' if2s=$septrace n2s=1 | psgraph ', &
     ' label1="Time (s)" label2="Amplitude (m)" ', &
     ' n=$nt style=normal d1=$dt hbox=4.0 wbox=6.0', &
     ' title="Ux acceleration component trace "$trace > uxtrace{$trace}.ps')
  225  format('subset < Uz_file n1=$nt', &
     ' if2s=$septrace n2s=1 | psgraph ', &
     ' label1="Time (s)" label2="Amplitude (m)" ', &
     ' n=$nt style=normal d1=$dt hbox=4.0 wbox=6.0', &
     ' title="Uz displacement component trace "$trace > uztrace{$trace}.ps')
  226  format('subset < Uz_file n1=$nt', &
     ' if2s=$septrace n2s=1 | psgraph ', &
     ' label1="Time (s)" label2="Amplitude (m)" ', &
     ' n=$nt style=normal d1=$dt hbox=4.0 wbox=6.0', &
     ' title="Uz velocity component trace "$trace > uztrace{$trace}.ps')
  227  format('subset < Uz_file n1=$nt', &
     ' if2s=$septrace n2s=1 | psgraph ', &
     ' label1="Time (s)" label2="Amplitude (m)" ', &
     ' n=$nt style=normal d1=$dt hbox=4.0 wbox=6.0', &
     ' title="Uz acceleration component trace "$trace > uztrace{$trace}.ps')
  130  format('sed -e ''1,$s/xwigb/pswigp/g'' ', &
       '-e ''1,$s/Ux_file/Ux_file > uxpoly.ps/g'' ', &
       '-e ''1,$s/Uz_file/Uz_file > uzpoly.ps/g'' ', &
       'xline > psline')
  137  format('sed -e ''1,$s/ //g'' -e ''1,$s/@/ /g'' ', &
      '-e ''1,1p'' -e ''$,$s/Ux/Uz/g'' <tempfile > xline')
  140  format(f9.2)
  150  format('echo Extracting trace $trace...'/, &
    'subset < Ux_file_an n1=$nt if2s=',i1,' n2s=1 > Ux_num1'/, &
    'subset < Uz_file_an n1=$nt if2s=',i1,' n2s=1 > Uz_num1'/, &
    'subset < Ux_file n1=$nt if2s=$septrace n2s=1 > Ux_num2'/, &
    'subset < Uz_file n1=$nt if2s=$septrace n2s=1 > Uz_num2'/, &
    'cat Ux_num1 Ux_num2 > Ux_num'/, &
    'cat Uz_num1 Uz_num2 > Uz_num'/, &
    'suaddhead ns=$nt ftn=0 < Ux_num1 > Ux_num1_segy'/, &
    'suaddhead ns=$nt ftn=0 < Uz_num1 > Uz_num1_segy'/, &
    'suaddhead ns=$nt ftn=0 < Ux_num2 > Ux_num2_segy'/, &
    'suaddhead ns=$nt ftn=0 < Uz_num2 > Uz_num2_segy'/, &
    'echo Computing residuals...'/, &
    'suop2 Ux_num2_segy Ux_num1_segy op=diff > Ux_num_segy'/, &
    'suop2 Uz_num2_segy Uz_num1_segy op=diff > Uz_num_segy'/, &
    'sustrip head=/dev/null outpar=/dev/null ftn=0 ', &
    '<Ux_num_segy > Ux_num_resid'/, &
    'sustrip head=/dev/null outpar=/dev/null ftn=0 ', &
    '<Uz_num_segy > Uz_num_resid'/, &
    'cat Ux_num_resid >> Ux_num'/, &
    'cat Uz_num_resid >> Uz_num')
  151  format('echo Extracting trace $trace...'/, &
    'subset < Ux_file_an n1=$nt if2s=',i1,' n2s=1 > Ux_num1'/, &
    'subset < Uz_file_an n1=$nt if2s=',i1,' n2s=1 > Uz_num1'/, &
    'subset < Ux_file n1=$nt if2s=$septrace n2s=1 > Ux_num'/, &
    'subset < Uz_file n1=$nt if2s=$septrace n2s=1 > Uz_num'/, &
    'suaddhead ns=$nt ftn=0 < Ux_num1 > Ux_num1_segy'/, &
    'suaddhead ns=$nt ftn=0 < Uz_num1 > Uz_num1_segy'/, &
    'suaddhead ns=$nt ftn=0 < Ux_num > Ux_num2_segy'/, &
    'suaddhead ns=$nt ftn=0 < Uz_num > Uz_num2_segy'/, &
    'echo Computing residuals...'/, &
    'suop2 Ux_num2_segy Ux_num1_segy op=diff > Ux_num_segy'/, &
    'suop2 Uz_num2_segy Uz_num1_segy op=diff > Uz_num_segy'/, &
    'sustrip head=/dev/null outpar=/dev/null ftn=0 ', &
    '<Ux_num_segy > Ux_num_resid'/, &
    'sustrip head=/dev/null outpar=/dev/null ftn=0 ', &
    '<Uz_num_segy > Uz_num_resid')
  152  format('echo Multiplying residuals by $ampli ...',/, &
    '/bin/rm -f prog_awk',/, &
    'echo \{print NR\*$dt ,  \$1\*=$ampli\} > prog_awk',/, &
    'b2a n1=1 outpar=/dev/null < Ux_num_resid | awk -f prog_awk', &
    ' > Ux_num_resid_asc_mul',/, &
    'b2a n1=1 outpar=/dev/null < Uz_num_resid | awk -f prog_awk', &
    ' > Uz_num_resid_asc_mul',/, &
    '/bin/rm -f prog_awk',/, &
    'echo \{print NR\*$dt ,  \$1\} > prog_awk',/, &
    'b2a n1=1 outpar=/dev/null < Ux_num | awk -f prog_awk', &
    ' > Ux_num_asc_txt',/, &
    'b2a n1=1 outpar=/dev/null < Uz_num | awk -f prog_awk', &
    ' > Uz_num_asc_txt')
  154  format('echo Generating PostScript files...',/, &
    'gnuplot << EOF',/, &
    'set output "uxresid$trace.ps"',/, &
    'set term postscript landscape color solid "Helvetica" 22',/, &
    'set xrange [0:$tottime]',/, &
    'set title "Ux residuals trace $trace"',/, &
    'set xlabel "Time (s)"',/, &
    'set ylabel "Amplitude (m)"',/, &
    'set nozeroaxis',/, &
    'set data style lines',/, &
    'plot "Ux_num_asc_txt" us 1:2 title ', &
    '"Numerical results" w l 1,', &
    ' "Ux_num_resid_asc_mul" us 1:2 title ', &
    '"Residuals * $ampli" w l 2',/, &
    'EOF')
  155  format('gnuplot << EOF',/, &
    'set output "uzresid$trace.ps"',/, &
    'set term postscript landscape color solid "Helvetica" 22',/, &
    'set xrange [0:$tottime]',/, &
    'set title "Uz residuals trace $trace"',/, &
    'set xlabel "Time (s)"',/, &
    'set ylabel "Amplitude (m)"',/, &
    'set nozeroaxis',/, &
    'set data style lines',/, &
    'plot "Uz_num_asc_txt" us 1:2 title ', &
    '"Numerical results" w l 1,', &
    ' "Uz_num_resid_asc_mul" us 1:2 title ', &
    '"Residuals * $ampli" w l 2',/, &
    'EOF')
  160  format('xgraph -geometry 1085x272 label1="Time (s)" ', &
    'label2="Amplitude (m)"  title="Ux residuals trace $trace ', &
    '(blue=Numerical red=Analytical green=Residuals)" ', &
    'n=$nt style=normal d1=$dt nplot=3 linecolor=2,4,3 < Ux_num &',/, &
    'xgraph -geometry 1085x272 label1="Time (s)" ', &
    'label2="Amplitude (m)"  title="Uz residuals trace $trace ', &
    '(blue=Numerical red=Analytical green=Residuals)" ', &
    'n=$nt style=normal d1=$dt nplot=3 linecolor=2,4,3 < Uz_num &')
  170 format('/bin/rm -f Ux_num1 Ux_num2 Ux_num Ux_num1_segy Ux_num2_segy ', &
    'Ux_num_segy Ux_num_resid Ux_num_segy Ux_num_resid_asc_mul ', &
    'Ux_num_asc_txt Uz_num1 Uz_num2 Uz_num Uz_num1_segy Uz_num2_segy ', &
    'Uz_num_segy Uz_num_resid Uz_num_segy Uz_num_resid_asc_mul Uz_num_asc_txt')
  171 format('/bin/rm -f prog_awk')
  200 format('set tottime = `echo $dt | awk ''{ print $1*i }'' i=$nt `', &
    /,'echo Total time $tottime seconds...')
  320 format('subset n1=$nt if2s=$septraceana',i1,' n2s=1 < U',a1,'_file ', &
    '> tutucomp ; subset n1=$nt if2s=$septraceref n2s=1 < U',a1,'_file_an ', &
    '> tutuan')
  330 format('cat tutuan tutucomp | xgraph -geometry 1085x272 ', &
    'linecolor=2,4 label1="Time (s)" label2="Amplitude (m)" nplot=2 ', &
    'n=$nt,$nt style=normal d1=$dt title="U',a1,' component numerical ', &
    '(blue) and analytical (red) trace "$traceana',i1,' &')
  340 format('cat tutuan tutucomp | psgraph hbox=4.0 wbox=6.0 ', &
    'linecolor=red,blue label1="Time (s)" label2="Amplitude (m)" nplot=2 ', &
    'n=$nt,$nt style=normal d1=$dt title="U',a1,' numerical (blue) ', &
    'and analytical (red) trace "$traceana',i1,' > u',a1,'tracecompana',i1,'.ps')

  return
  end subroutine writeseis
