#!/usr/bin/perl -w

#==========================================================
#
# plot_wavefield.pl
# Carl Tape, 11-April-2011
#
# This perl script writes a shell script that calls GMT to plot kernels or snapshots of the wavefield.
# Most of the input commands are associated with plotting.
# Anyone is welcome/encouraged to improve this script.
# It has been tested for the Tromp2005 and Tape2007 examples.
#
# NOTES:
#   1. must have output_wavefield_snapshot = .true. in Par_file to generate wavefield snapshots
#   2. boundaries between material regions are NOT plotted here
#
# WAVEFIELD EXAMPLES:
#    plot_wavefield.pl 100/1600/400 200/800/200  0/4/0/4      1/0.1/1/0.1  -3/-3/-3 6/6/6  0.0/0.001 1/0/1/1/1 1.7/1/1 1/0/1/120 M2_UPPA PSV
#    plot_wavefield.pl 400/2800/400 400/2000/400 0/200/0/80    50/10/40/10 -3/-3/-3 4/4/4  -8.0/0.02 1/0/1/1/1 3.0/1/0 1/0/1/200 Tromp2005 PSV_homo
#    plot_wavefield.pl 400/2800/400 400/2000/400 0/200/0/80    50/10/40/10 -2/-2/-2 1/1/1  -8.0/0.02 0/1/0/0/1 3.0/1/0 1/0/1/200 Tromp2005 SH_homo
#    plot_wavefield.pl 400/4800/400 400/2800/800 0/480/0/480 120/20/120/20 -3/-3/-3 6/6/6 -48.0/0.06 0/1/0/0/1 1.7/1/0 1/0/1/120 Tape2007 onerec_homo
#    plot_wavefield.pl 400/4800/400 800/2000/400 0/480/0/480 120/20/120/20 -3/-3/-3 6/6/6 -48.0/0.06 0/1/0/0/1 1.7/1/0 1/0/1/120 Tape2007 132rec_checker
#
# KERNEL EXAMPLES:
#    plot_wavefield.pl 400/3000/400 400/2800/800 0/480/0/480 120/20/120/20 -8/-8/-8 1/1/1 -48.0/0.06 0/1/0/0/0 2.0/1/0 1/0/1/120 Tape2007_kernel onerec_homo
#    plot_wavefield.pl 400/2800/400 400/2000/400 0/200/0/80    50/10/40/10 -8/-8/-8 1/1/1  -8.0/0.02 1/0/1/1/0 4.0/1/0 1/0/1/200 Tromp2005_kernel PSV_homo
#
#
#----------------------------------------------
# USER INPUT

if (@ARGV < 1) {die("Usage: plot_wavefield.pl xxx\n");}
($range,$ftemp,$R0,$ticks,$ptemp,$ctemp,$tinc,$ipar1,$ipar2,$ipar3,$tlab,$ttag) = @ARGV;

($rfirst,$rend,$rint) = split("/",$range);                 # all available frames
($pfirst,$pend,$pint) = split("/",$ftemp);                 # frames to plot
($xmin,$xmax,$zmin,$zmax) = split("/",$R0);                # region, km
($btick1x,$btick2x,$btick1z,$btick2z) = split("/",$ticks); # ticks, km
@pwr     = split("/",$ptemp);     # PWR  : increase (larger negative power) for more contrast
@cmax    = split("/",$ctemp);     # CMAX : decrease for more contrast
($t0,$dt) = split("/",$tinc);
($ipx,$ipy,$ipz,$ipsv,$itype) = split("/",$ipar1);   # itype =0 (kernel), =1 (forward wavefield), (=2) adjoint wavefield
($xwid,$iportrait,$irecsurf) = split("/",$ipar2);
($igrd,$imask,$msize,$xpix) = split("/",$ipar3);

#@frames  = split("/",$ftemp);
#$numf = @frames;
$numf = ($pend - $pfirst)/$pint + 1;

###############################################################
##
## setup
##
###############################################################
# directory with data files (USER CHANGE THIS)
$outdir = "OUTPUT_FILES";

# kernel type
$iktype = 2;   # =1 rho-kappa-mu, =2 rhop-Vp-Vs

# plot the color frames or not
$icolor = 1;    # ccc
$ibound = 0;    # outlines of model units (yahtse, yakutat)
$iseis = 0;     # seismicity (yakutat)
$iplotrec = 1;
$iplotreclab = 1;
$iheader = 1;
#$igrd = 0;
#$imask = 1;

# region -- in km (but the files are in meters)
#$xmin = 0; $xmax = 160; $zmin = -60; $zmax = 0;
#$btick1x = 40; $btick2x = 10;
#$btick1z = 20; $btick2z = 10;

$xgap = 0.25;
$ygap = 1.0;
$originY = 1.25;
#$msize = 0.25;    # KEY: marker size for plotting points or pixels

# gmt version (=1 new, =0 older gmt4 version)
$gmt5 = 1;
###############################################################

# you can specify you own path by changing $bdir and $idir1
#$bdir = "/data2/SVN/seismo/2D/SPECFEM2D_20120420";
#$idir1 = "$bdir/OUTPUT_FILES";                 # if running from the default directory
#$idir1 = "$bdir/EXAMPLES/$tlab/OUTPUT_FILES";   # if running from an examples directory
$bdir = $ENV{'PWD'};
$idir1 = "$bdir/$outdir";
if (not -e $idir1) {die("check if directory $idir1 exist or not\n");}


$R = "-R$R0";
$xran = $xmax - $xmin;
$zran = $zmax - $zmin;
$xzratio = $xran/$zran;

# if ($xzratio > 4) {
#   $orient = " ";
#   $xwid = 4.5;
# } else {
#   $orient = "-P";
#   if ($xzratio > 2) {$xwid = 3.0;} else {$xwid = 2.0;}
# }
#$xwid = 1.5;

print "\n xzratio = $xzratio ; xwid = $xwid --\n";

#$orient = " "; $xwid = 9.0;
#if($itype==0) {$orient = "-P"; $xwid = $xwid*1.5; $sc = 3;}

# resolution of color plots -- OR you can plot the pixels as circles
$ypix = int($xpix/$xzratio);
$interp = "-I${xpix}+/${ypix}+";
$grdfile = "temp.grd";
$minfo = "-Sc${msize}p";

if($igrd==1) {print "interpolation for grd file is $interp\n";}

#----------------------------------------------

if($iportrait==1) {
  $orient = "-P";
} else {
  $orient = "";
}

$ncol = $ipx + $ipy + $ipz;
if($ncol==0) {
  print "Error $ncol (ncol) to display: $ipx (X), $ipy (Y), $ipz (Z)\n";
  die("Must specify at least one component to display");
}

if($ncol==1 && $ipx==1) {@comps = (1);}
if($ncol==1 && $ipy==1) {@comps = (2);}
if($ncol==1 && $ipz==1) {@comps = (3);}
if($ncol==2 && $ipx==0) {@comps = (2,3);}
if($ncol==2 && $ipy==0) {@comps = (1,3);}
if($ncol==2 && $ipz==0) {@comps = (1,2);}
if($ncol==3) {@comps = (1,2,3);}

print "$ncol (ncol) to display: $ipx (X), $ipy (Y), $ipz (Z)\n";
print "comps: @comps\n";
#die("TESTING");

if($ncol==3) {
  $originX = 0.75;
} else {
  $originX = 1.5;
}
$origin = "-X$originX -Y$originY";

# yakutat
#if($ncol==1) {$xwid = 9; $orient = " ";}
#if($ncol==1 && $itype==0) {$xwid = 6.5; $orient = "-P";}

$ywid = $xwid * $zran/$xran;
$dX0 = $xwid + $xgap;
$dY0 = $ywid + $ygap;
$dX2 = 2*$dX0;
$dX     = "-X$dX0";
$mdX    = "-X-$dX0";
$mdX2   = "-X-$dX2";
$dY     = "-Y$dY0";

# plotting specifications
$fsize0 = "18";
$fsize1 = "10";
$fsize2 = "10";
$fontno = "1";
$tick   = "0.1c";
$cgray = 200;

if($imask==1) {
  $BG = " ";
} else {
  if($gmt5==1){
    $BG = "+g$cgray";
  }else{
    # older gmt
    $BG = "-G$cgray";
  }
};
$stitype = sprintf("%02d",$itype);

# plot symbols for sources, receivers
$rfill = "-G255";
$rfill = "";
$sfill = "-G255";
$sfill = "";
#$src = "-W0.5p -Sa0.2 $sfill";
$src = "-W1.0p $sfill -Sa12p";
$rdx = 0; $rdy = 0; $tdx = 0; $tdy = 10;
if($irecsurf == 1) {$rdy = 6; $tdy = 15;}

# -N or not
$rec = "-W1p $rfill -Si10p -D${rdx}p/${rdy}p";
#$rec2 = "-W1p,0/255/255 $rfill -Si10p -D${rdx}p/${rdy}p";
$textrec = "-D${tdx}p/${tdy}p -W255 -C1p -N";

# source and receivers
$srcfile = "$idir1/SOURCE";
$recfile = "$idir1/STATIONS";
#$recfile2 = "$idir1/STATIONS_target";
if (not -e $srcfile) {die("check if sfile $srcfile exist or not\n");}
if (not -e $recfile) {die("check if rfile $recfile exist or not\n");}
$srcx = `grep "xs                              =" $srcfile | awk '{print \$3}'`; chomp($srcx);
$srcz = `grep "zs                              =" $srcfile | awk '{print \$3}'`; chomp($srcz);
$srcx = $srcx/1000;
$srcz = $srcz/1000;

print "\nsource at ($srcx, $srcz)\n";

$R = "-R$xmin/$xmax/$zmin/$zmax";
print "\nregion is    : $R\n";

# projection
#$ywid = $xwid*($zran/$xran);
$J = "-JX${xwid}i/${ywid}i";
print "projection is: $J \n\n";

#=================================================================
#
# SCRIPT
#
#=================================================================

$shfile = "plot_wavefield.sh";
open(CSH,">$shfile");
print "\nWriting CSH file: $shfile\n";
print CSH "#!/bin/bash\n\n";
print CSH "echo \"##################################\"\n";
print CSH "echo \"run script $shfile\"\n";
print CSH "echo \"##################################\"\n";
print CSH "echo\n";

if($gmt5==1){
  print CSH "gmtset MAP_FRAME_TYPE plain PS_MEDIA letter PROJ_LENGTH_UNIT inch FORMAT_GEO_MAP D MAP_TICK_LENGTH_PRIMARY $tick  PS_CHAR_ENCODING Standard+ COLOR_NAN $cgray\n";
}else{
  # older gmt
  print CSH "gmtset BASEMAP_TYPE plain PAPER_MEDIA letter MEASURE_UNIT inch PLOT_DEGREE_FORMAT D TICK_LENGTH $tick LABEL_FONT_SIZE $fsize2 ANOT_FONT_SIZE $fsize2  HEADER_FONT $fontno ANOT_FONT $fontno LABEL_FONT $fontno HEADER_FONT_SIZE $fsize1 CHAR_ENCODING Standard+ COLOR_NAN $cgray\n";
}


# KEY: scaling for color
$scale_color = 101; # 21;
$colorbar = "seis";
@norm = ("1e$pwr[0]","1e$pwr[1]","1e$pwr[2]");
#$fac = 5;      # KEY: enhance the interaction field (increase for more contrast)
print "norm = @norm \n";

$numw = @cmax;
for ($k = 0; $k < $numw; $k++ ){

  $ss[$k] = $cmax[$k];  # to ensure that the max is displayed in the scalebar
  $ds[$k] = 2*$ss[$k]/$scale_color;
  #$bs[$k] = sprintf("%3.3e",0.9*$ss[$k]);  # colorbar
  $bs[$k] = sprintf("%3.3e",$ss[$k]);  # colorbar
  $Ts[$k] = sprintf("-T%3.3f/%3.3f/%3.3f",-$ss[$k]*1.05,$ss[$k]*1.05,$ds[$k]);
  print "color range Ts = $Ts[$k]\n";

  print CSH "makecpt -C$colorbar $Ts[$k] -D > color_${k}.cpt\n";
}

# colorbar
$Dx = ($xwid - $originX)/2;
if($gmt5==1){
  $Dscale = "-Dx$Dx/-0.5+w1.5/0.15+h+e10p";
}else{
  # older gmt
  $Dscale = "-D$Dx/-0.3/1.5/0.15h -E10p";
}
print "color scale position D = $Dscale\n";

#=================================================================

if ($itype != 0) {

#=================================================================
#
# WAVEFIELDS
#
#=================================================================

  print "\nwavefields\n\n";

  # labels
  @wavefield = ("wavefield","wavefield","wavefield");
  @titles    = ("Ux","Uy","Uz");

  $name = "wavefield_${tlab}_${ttag}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  print CSH "echo\n";
  print CSH "echo \"postscript file = $psfile\";echo\n";
  print CSH "echo\n";
  # clean
  print CSH "rm -f $psfile\n\n";

  #$numf = 1;

  $imin = 0; $imax = $numf-1;  # default
  #$imin = 4; $imax = $imin;

  #for ($j = $rfirst; $j <= $rend; $j = $j + $rint) {
  for ($i = $imin; $i <= $imax; $i++) {

   #$j1 = $frames[$i];           # forward frame
   $j1 = $pfirst + $i*$pint;
   $j2 = $rend + $rfirst - $j1;   # corresponding adjoint frame
   $snap1 = sprintf("%07d",$j1);
   $snap2 = sprintf("%07d",$j2);
   #$time = sprintf("%04d",$j1*$dt);
   $time = sprintf("%.3f",$t0 + $j1*$dt);
   #$time = sprintf("%.0f",$j1);     # snapshot
   #print "\n--time,t0,j1,i -- $time, $t0, $j1, $i--\n";

   $snapshot_f = "${idir1}/$wavefield[0]${snap1}_${stitype}_000.txt";   # 000 assumes one processor only
   #if($igrd==1) {$snapshot_f = "${idir1}/$wavefield[0]${snap1}_${stitype}_pixel.txt";}

   #$snapshot_a = "${idir1}/$wavefield[1]${snap2}.txt";
   #$snapshot_k = "${idir1}/$wavefield[2]${snap2}.txt";

   if (not -f $snapshot_f) {die("check if snapshot_f $snapshot_f exist or not\n");}
   #if (not -f $snapshot_a) {die("check if snapshot_a $snapshot_a exist or not\n");}
   #if (not -f $snapshot_k) {die("check if snapshot_k $snapshot_k exist or not\n");}

   print CSH "echo $psfile\n";
   print CSH "echo $snapshot_f\n";

   if($gmt5==1){
    $B0 = sprintf("-Bxa${btick1x}f${btick2x} -Bya${btick1z}f${btick2z}+l\"t = $time s\"");
    $B      = "$B0 -BWsne";
    $B_row1 = "$B0 -BWSne";
   }else{
    # older gmt
    $B0 = sprintf("-Ba${btick1x}f${btick2x}/a${btick1z}f${btick2z}:\"t = $time s\"::.\"  \"");
    $B      = "$B0:Wsne";
    $B_row1 = "$B0:WSne";
   }
   if ($i == $imin) { $B = $B_row1;}

   #-------------------------
   # loop over components

   #$kmin = 0; $kmax = 1;
   #$kmin = 1; $kmax = $kmin;  # zcomp
   #$kmin = 0; $kmax = $kmin;  # xcomp

   for ($k = 1; $k <= $ncol; $k++) {

     $nm = $norm[$k-1];
     $comp = $comps[$k-1];
     $cfile = sprintf("color_%i.cpt",$k-1);
     if ($k==1) {
       if($ncol==1) {$shift = "$dY";}
       if($ncol==2) {$shift = "$dY $mdX";}
       if($ncol==3) {$shift = "$dY $mdX2";}
     } else {
       $shift = "$dX";
     }

     if($gmt5==1){
      $BscaleS = sprintf("-B%2.2f+l\"%s (x, z, t), 10\@+%2.2i\@+  m\"",$bs[$comp-1],$titles[$comp-1],$pwr[$comp-1]);
     }else{
      $BscaleS = sprintf("-B%2.2e:\"%s (x, z, t), 10\@+%2.2i\@+  m\":",$bs[$comp-1],$titles[$comp-1],$pwr[$comp-1]);
     }

     if($gmt5==1){
       $B = "$B0 -BWesn";
       if($k > 1) {$B = "$B0 -Bwesn";}
       if($i==$imin && $k==1) {$B = "$B0 -BWeSn";}
       if($i==$imin && $k > 1) {$B = "$B0 -BweSn";}
     }else{
       $B = "$B0:Wesn";
       if($k > 1) {$B = "$B0:wesn";}
       if($i==$imin && $k==1) {$B = "$B0:WeSn";}
       if($i==$imin && $k > 1) {$B = "$B0:weSn";}
     }
     #if ($k==0) {
     #  $comp=1; $shift = "$dY $mdX"; $BscaleS = $BscaleS1;
     #} else {
     #  $comp=3; $shift = $dX; $B = "$B0:wesn"; $BscaleS = $BscaleS2;
     #}

     #if ($kmin==$kmax) {
     #  $shift = $dY; $B = "$B0:WeSn";
     #}       # one column only
     ##if ($i == $imin) {$B = "$B0:weSn";} else {$B = "$B0:wesn";}    # temp, zcomp
     ##if ($i == $imin) {$B = $B_row1;} else {$B = "$B0:weSn";}    # temp
     #if ($i == $imin && $k == 1) {
     #  $B = "$B0:weSn";
     #}       # temp

     if ($i == $imin && $k==1) {
       print CSH "psbasemap $J $R $B $BG -K -V $orient $origin > $psfile\n";
     }        # START
     else {
       print CSH "psbasemap $J $R $B $BG -K -O -V $shift >> $psfile\n";
     }

     # PLOT THE FORWARD WAVEFIELD
     if ($icolor==1) {
       if ($igrd==1) {
          #print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2) / $nm}' $snapshot_f > dfile\n";
          #print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2) / $nm}' $snapshot_f | nearneighbor -G$grdfile $R $interp\n";
          print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2) / $nm}' $snapshot_f | xyz2grd -G$grdfile $R $interp\n";
          print CSH "grdimage $grdfile -C$cfile $J -K -O -V -Q >> $psfile\n";
       } else {
          print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2) / $nm}' $snapshot_f | psxy $R $J $minfo -C$cfile -K -O -V >> $psfile\n";
       }

       # mask points above topography
       if ($imask==1) {
          #print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2)}' $snapshot_f > maskpts\n";
          printf CSH "grep NaN $snapshot_f > maskpts\n";
          print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2)}' maskpts | psxy $R $J $minfo -C$cfile -K -O -V >> $psfile\n";
       }
     }

     # plot the boundaries of the geometrical units in the model
     if ($ibound==1) {
      $bfile = "/home/carltape/PROJECTS/yakutat/profiles/yakutat_p1_m3.xy";
      #$bfile = "/home/carltape/PROJECTS/yahtse/data/icybay_chris/profiles_specfem2d/profile_p1.xy";
      #$bfile = "/home/carltape/PROJECTS/yahtse/data/icybay_chris/profiles_specfem2d/profile_p1_top.xy";
      if (not -f $bfile) {
       die("check if bfile $bfile exist or not\n");
      }
      #print CSH "awk '{print \$1/1000,\$2/1000}' $bfile | psxy $J $R -W1.5p -m -K -O -V >> $psfile\n";
      print CSH "awk '{print \$1,\$2}' $bfile | psxy $J $R -W1p -m -K -O -V >> $psfile\n";
     }

     # plot seismicity
     if ($iseis==1 && $i==$imin) {
        $bfile = "/home/carltape/PROJECTS/yakutat/data/yakutat_p1_seis_km_xyz.dat";
        if (not -f $bfile) {
          die("check if bfile $bfile exist or not\n");
        }
        print CSH "awk '{print \$1,\$2}' $bfile | psxy $J $R -Sc3p -G0 -K -O -V >> $psfile\n";
     }

     if ($i == $imin) {print CSH "psscale -C$cfile $Dscale $BscaleS -K -O -V >> $psfile \n";}
     print CSH "psbasemap $J $R $B -K -O -V >> $psfile\n";

     # label the figures
     if (0==1) {
       if ($i==$imin) {
          if($gmt5==1){
          }else{
            print CSH "pstext $J $R -K -O -V >>$psfile<<EOF\n -2 -0.5 14 0 $fontno CM rock\nEOF\n";
            print CSH "pstext $J $R -K -O -V >>$psfile<<EOF\n -1 0.75 14 0 $fontno CM ice\nEOF\n";
            print CSH "pstext $J $R -K -O -V >>$psfile<<EOF\n 1.0 0.75 14 0 $fontno CM water\nEOF\n";
          }
          print CSH "psxy -W1p $J $R -K -O -V >>$psfile<<EOF\n1.0 0.6\n 0.7 -0.1\nEOF\n";
          print CSH "psxy -W1p $J $R -K -O -V >>$psfile<<EOF\n-1.0 0.6\n -0.7 -0.05\nEOF\n";
       }
     }

     # plot stations and labels
     if($iplotrec==1) {
        #print CSH "awk '{print \$3/1000,\$4/1000}' $recfile2 | psxy $J $R -K -O -V $rec2 >> $psfile\n";
        print CSH "awk '{print \$3/1000,\$4/1000}' $recfile | psxy $J $R -K -O -V $rec >> $psfile\n";
        #if($i==$imin && $iplotreclab==1) {print CSH "awk '{print \$3/1000,\$4/1000,8,0,$fontno,\"CM\",\"S\"\$7}' $recfile | pstext $textrec $J $R -K -O -V >> $psfile\n";}
        if($i==$imin && $iplotreclab==1 && $k==1) {print CSH "awk '{print \$3/1000,\$4/1000,8,0,$fontno,\"CM\",\$1}' $recfile | pstext $textrec $J $R -K -O -V >> $psfile\n";}
     }
     print CSH "psxy $J $R -K -O -V $src >> $psfile<<EOF\n $srcx $srcz\nEOF\n";

     # plot the time of the snapshot (for some reason, it won't work inside the B command)
     #$xtext = $xmin-0.3*$xran;
     #$ztext = $zmin+0.5*$zran;
     #$tstr = "t = $time s";
     #print CSH "pstext -N $J $R -K -O -V >>$psfile<<EOF\n $xtext $ztext $fsize1 90 $fontno CM $tstr\nEOF\n";

     #$xtx = $xmin+0.5*$xran; $ztx = $zmin+1.1*$zran;
     #if ($i == $imax) {print CSH "pstext -N $J $R -K -O -V >>$psfile<<EOF\n $xtx $ztx $fsize1 0 $fontno CM $titles[$k]\nEOF\n";}

     #-------------------------

#    $comp = 3;
#    $k = 1;

#    $B      = "$B0:wesn";
#    $B_row1 = "$B0:weSn";
#    if ($i == $imin) { $B = $B_row1;}

#    # PLOT THE FORWARD WAVEFIELD -- Uz
#    print CSH "psbasemap $J $R $B $BG -K -O -V $dX >> $psfile\n";
#    if($icolor==1) {
#       if($igrd==1) {
#         #print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2) / $nm}' $snapshot_f > dfile\n";
#         #print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2) / $nm}' $snapshot_f | nearneighbor -G$grdfile $R $interp\n";
#         print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2) / $nm}' $snapshot_f | xyz2grd -G$grdfile $R $interp\n";
#         print CSH "grdimage $grdfile -C$cfile $J -K -O -V -Q >> $psfile\n";
#       } else {
#         print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2) / $nm}' $snapshot_f | psxy $R $J $minfo -C$cfile -K -O -V >> $psfile\n";
#       }

#       # mask points above topography
#       if($imask==1) {
#       #print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2)}' $snapshot_f > maskpts\n";
#       printf CSH "grep NaN $snapshot_f > maskpts\n";
#       print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2)}' maskpts | psxy $R $J $minfo -C$cfile -K -O -V >> $psfile\n";
#       }
#    }
#    #print CSH "pscoast $J $R $B -W1p -Na/1p -Dh -K -O -V >> $psfile\n";
#    #print CSH "awk '{print \$2,\$1}' INPUT/oms_shelf |psxy $J $R $Wshelf -K -O -V >> $psfile\n";
#    print CSH "psbasemap $J $R $B -K -O -V >> $psfile\n";
#    if ($i == $imin) {print CSH "psscale -C$cfile $Dscale $BscaleS2 -K -O -V >> $psfile \n";}

#    # plot simple numerical index for the station label (otherwise, use column 1)
#    print CSH "awk '{print \$3/1000,\$4/1000}' $recfile | psxy $J $R -K -O -V $rec >> $psfile\n";
#    if($i == $imin) {print CSH "awk '{print \$3/1000,\$4/1000,8,0,$fontno,\"CM\",\"S\"\$7}' $recfile | pstext $textrec $J $R -K -O -V >> $psfile\n";}
#    print CSH "psxy $J $R -K -O -V $src >> $psfile<<EOF\n $srcx $srcz\nEOF\n";

#    if ($i == $imax) {print CSH "pstext -N $J $R -K -O -V >>$psfile<<EOF\n $xtx $ztx $fsize1 0 $fontno CM $titles[1]\nEOF\n";}

    }  # for $m

  }

  # plot title and GMT header
  if($iheader==1) {
    $plabel = "plot_wavefield.pl";
    $ux = -$xwid;
    $ux = 0;
    $uy = $ywid + 0.3;
    $Utag = "-U/$ux/$uy/$plabel"; # GMT header
    $shift = "-X0i -Y0.7i";
    if($ncol==2) {$shift = "-X-${dX0}i -Y0.7i";}
    if($ipsv==1) {$tlab0 = "PSV"} else {$tlab0 = "SH"}
    $title = "$tlab0 wavefield -- $tlab";
    if($gmt5==1){
      print CSH "pstext -N $J $R $Utag -K -O -V $shift -F+f+a+j >>$psfile<<EOF\n $xmin $zmax $fsize0,$fontno 0 LM $title\nEOF\n";
    }else{
      print CSH "pstext -N $J $R $Utag -K -O -V $shift >>$psfile<<EOF\n $xmin $zmax $fsize0 0 $fontno LM $title\nEOF\n";  # LM or CM
    }
  }


} else {

#=================================================================
#
# KERNELS
#
#=================================================================

  print "\nkernels\n\n";

  # labels for kernels
  $ik = 1;     # index for kernel, e.g., the time window for a single seismogram (Yahtse)
  $sik = $ik;
  $sik = "";

  if ($iktype==1) {
    @titles    = ("K$sik-\@~r\@","K$sik-\@~k\@","K$sik-\@~m\@");
    @ytitles    = ("Krho","Kkappa","Kmu");
    $kfile1 = "${idir1}/proc000000_rho_kappa_mu_kernel.dat";
  } else {
    @titles    = ("K$sik-\@~r\@","K$sik-Vp","K$sik-Vs");
    @ytitles    = ("Krho","Kalpha","Kbeta");
    $kfile1 = "${idir1}/proc000000_rhop_alpha_beta_kernel.dat";
  }

  # kernel file for plotting
  if (not -f $kfile1) {die("check if kernel file $kfile1 exist or not\n");}
  $kfile = $kfile1;

  $name = "kernel_${tlab}_${ttag}";
  $psfile  = "$name.ps";
  $jpgfile = "$name.jpg";

  print CSH "echo\n";
  print CSH "echo \"kernel file     = $kfile\"\n";
  print CSH "echo \"postscript file = $psfile\";echo\n";
  print CSH "echo\n";
  # clean
  print CSH "rm -f $psfile\n\n";

  if($gmt5==1){
    $BscaleS = sprintf("-B%2.2e+l\"%s (x, z, t), 10\@+%2.2i\@+  m\@+-2\@+\" -BS",$bs[0],"K",$pwr[0]);
  }else{
    # older gmt
    $BscaleS = sprintf("-B%2.2e:\"%s (x, z, t), 10\@+%2.2i\@+  m\@+-2\@+\":",$bs[0],"K",$pwr[0]);
    #$BscaleS = sprintf("-B%2.2e:\" \":",$bs[0],"K",$pwr[0]);
  }
  print CSH "echo \"scale B option \"$BscaleS\" \"\n";

  if($gmt5==1){
    $B1 = sprintf("-Bxa${btick1x}f${btick2x} -Bya${btick1z}f${btick2z}+l\"$ytitles[0]\" -BWSne$BG");
    $B2 = sprintf("-Bxa${btick1x}f${btick2x} -Bya${btick1z}f${btick2z}+l\"$ytitles[1]\" -BWsne$BG");
    $B3 = sprintf("-Bxa${btick1x}f${btick2x} -Bya${btick1z}f${btick2z}+l\"$ytitles[2]\" -BWsne$BG");
  }else{
    # older gmt
    #$B0 = sprintf("-Ba${btick1x}f${btick2x}/a${btick1z}f${btick2z}:\"t = $time s\"::.\"  \"");
    $B1 = sprintf("-Ba${btick1x}f${btick2x}/a${btick1z}f${btick2z}:\"$ytitles[0]\"::.\" \":WSne $BG");
    $B2 = sprintf("-Ba${btick1x}f${btick2x}/a${btick1z}f${btick2z}:\"$ytitles[1]\"::.\" \":Wsne $BG");
    $B3 = sprintf("-Ba${btick1x}f${btick2x}/a${btick1z}f${btick2z}:\"$ytitles[2]\"::.\" \":Wsne $BG");
  }

  $i = 0;
  # forward frame
  $j1 = $pfirst + $i*$pint;
  # corresponding adjoint frame
  $j2 = $rend + $rfirst - $j1;
  $snap1 = sprintf("%07d",$j1);
  $snap2 = sprintf("%07d",$j2);
  $time = sprintf("%.3f",$j1*$dt);

  $kmin = 0; $kmax = 2;
  for ($k = $kmin; $k <= $kmax; $k++) {

    $comp = $k+1;
    if($k==0) {$B=$B1}; if($k==1) {$B=$B2}; if($k==2) {$B=$B3};

    print CSH "echo;echo \"B option $B\";echo;\n";

    # START
    if ($k == $kmin) {
      print CSH "psbasemap $J $R $B -K -V $orient $origin > $psfile\n";
    } else {
      print CSH "psbasemap $J $R $B -K -O -V $dY >> $psfile\n";
    }

    print CSH "echo;echo \"norm = $norm[$k]\"\n";

    if ($icolor==1) {
      print CSH "echo;echo \"using color\";\n";
      if ($igrd==1) {
        print CSH "echo;echo \"plotting by grdimage\";echo;\n";
        #print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2) / $norm[$k]}' $kfile > dfile\n";
        #print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2) / $norm[$k]}' $kfile | nearneighbor -G$grdfile $R $interp\n";
        print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2) / $norm[$k]}' $kfile | xyz2grd -G$grdfile $R $interp\n";
        print CSH "grdinfo $grdfile\n";
        print CSH "grdimage $grdfile -Ccolor_${k}.cpt $J -K -O -V -Q >> $psfile\n";
      } else {
        print CSH "echo;echo \"plotting by psxy\";\n";
        print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2) / $norm[$k]}' $kfile | psxy $R $J $minfo -Ccolor_${k}.cpt -K -O -V >> $psfile\n";
      }

      # mask points above topography
      if($imask==1) {
        print CSH "echo;echo \"using mask\";\n";
        #print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2)}' $kfile > maskpts\n";
        printf CSH "grep NaN $kfile > maskpts\n";
        print CSH "awk '{print \$1/1000,\$2/1000,\$($comp+2)}' maskpts | psxy $R $J $minfo -Ccolor_${k}.cpt -K -O -V >> $psfile\n";
      }
    }

    # plot the boundaries of the geometrical units in the model
    if ($ibound==1) {
      print CSH "echo;echo \"using boundaries\";\n";
      $bfile = "/home/carltape/PROJECTS/yakutat/profiles/yakutat_p1_m3.xy";
      #$bfile = "/home/carltape/PROJECTS/yahtse/data/icybay_chris/profiles_specfem2d/profile_p1.xy";
      #$bfile = "/home/carltape/PROJECTS/yahtse/data/icybay_chris/profiles_specfem2d/profile_p1_top.xy";
      if (not -f $bfile) {
        die("check if bfile $bfile exist or not\n");
      }
      #print CSH "awk '{print \$1/1000,\$2/1000}' $bfile | psxy $J $R -W1.5p -m -K -O -V >> $psfile\n";
      print CSH "awk '{print \$1,\$2}' $bfile | psxy $J $R -W1p -m -K -O -V >> $psfile\n";
    }

    #print CSH "pscoast $J $R $B -W1p -Na/1p -Dh -K -O -V >> $psfile\n";
    #print CSH "awk '{print \$2,\$1}' INPUT/oms_shelf |psxy $J $R $Wshelf -K -O -V >> $psfile\n";
    if ($k == $kmin) {
      print CSH "psscale -Ccolor_${k}.cpt $Dscale $BscaleS -K -O -V >> $psfile \n";
    }
    #print CSH "awk '{print \$3/1000,\$4/1000}' $recfile | psxy $J $R -K -O -V $rec >> $psfile\n";
    #print CSH "psxy -N $J $R -K -O -V $src >> $psfile<<EOF\n $srcx $srcz\nEOF\n";
    #print CSH "psbasemap $J $R $B -K -O -V >> $psfile\n";

    # label the figures
    if(0==1) {
      print CSH "echo;echo \"using labels\";\n";
      if($gmt5==1){
      }else{
        if($ik==1 && $k==1) {print CSH "pstext $J $R -K -O -V >>$psfile<<EOF\n -4.2 0 14 0 $fontno CM P\nEOF\n";}
        if($ik==2 && $k==2) {print CSH "pstext $J $R -K -O -V >>$psfile<<EOF\n -4.2 0 14 0 $fontno CM S\nEOF\n";}
        if($ik==3 && $k==1) {print CSH "pstext $J $R -K -O -V >>$psfile<<EOF\n -2.8 -0.3 14 0 $fontno CM P\nEOF\n";}
        if($ik==3 && $k==2) {print CSH "pstext $J $R -K -O -V >>$psfile<<EOF\n -2.4 0.4 14 -5 $fontno CM Rayleigh\nEOF\n";}
        if($ik==4 && $k==2) {
          print CSH "pstext $J $R -K -O -V >>$psfile<<EOF\n -4.2 0 14 0 $fontno CM S\nEOF\n";
          print CSH "pstext $J $R -K -O -V >>$psfile<<EOF\n -2.4 0.4 14 -5 $fontno CM Rayleigh (x2)\nEOF\n";
        }
      }
    }

    # plot title
    $xtx = $xmin+0.98*$xran; $ztx = $zmin+0.93*$zran;
    $textinfo = "-G255"; #"-G0 -W255 -C2p -N";
    print "\ntitle plotting at $xtx, $ztx info: $textinfo $titles[$k]\n";
    if($gmt5==1){
      print CSH "pstext $textinfo $J $R -K -O -V -F+f+a+j >>$psfile<<EOF\n $xtx $ztx 16,$fontno 0 RT $titles[$k]\nEOF\n";
    }else{
      # older gmt
      print CSH "pstext $textinfo $J $R -K -O -V >>$psfile<<EOF\n $xtx $ztx 16 0 $fontno RT $titles[$k]\nEOF\n";
    }

    print CSH "echo \"iteration $k done\";echo\n\n\n";
  }

  if($iheader==1) {
    # plot title and GMT header
    print CSH "echo;echo \"plotting header\";echo;\n";
    $plabel = "plot_wavefield.pl";
    $ux = 0;
    $uy = $ywid + 0.3;
    $Utag = "-U/$ux/$uy/$plabel"; # GMT header
    $shift = "-X0i -Y0.7i";
    if ($ipsv==1) {$tlab0 = "PSV";} else {$tlab0 = "SH";}
    $title = "$tlab0 kernel -- $tlab";
    if($gmt5==1){
      print CSH "pstext -N $J $R $Utag -K -O -V $shift -F+f+a+j >>$psfile<<EOF\n $xmin $zmax $fsize0,$fontno 0 LM $title\nEOF\n";
    }else{
      # older gmt
      print CSH "pstext -N $J $R $Utag -K -O -V $shift >>$psfile<<EOF\n $xmin $zmax $fsize0 0 $fontno LM $title\nEOF\n";
    }
  }
}

#=================================================================

#-------------------------
# FINISH
#print CSH "pstext $J -R0/1/0/1 -O -V >>$psfile<<EOF\n 10 10 $fsize0 0 $fontno CM junk \nEOF\n";
print CSH "psxy $J -R0/1/0/1 -O -V >>$psfile<<EOF\nEOF\n\n";
# cleanup
print CSH "rm -f temp.grd color_0.cpt color_1.cpt color_2.cpt\n";
# user output
print CSH "\necho\necho \"see output psfile: $psfile\"\necho\n";

# converts (using ImageMagick command)
#print CSH "convert $psfile $jpgfile\n";
close (CSH);

# executes script
system("sh -f $shfile");

# opens output ps-file
#system("gv $psfile &");

#=================================================================
