## sets terminal
set term png size 1000,1000
# w/ font arial
#set term png size 1000,1000 enhanced font "arial,10"

## sets color range

#set palette defined ( -8e-15 "red", 0 "yellow", 8e-15 "blue")
#set cbrange [-8e-15:8e-15]

set palette defined ( -1e-4 "red", 0 "yellow", 1e-4 "blue")
set cbrange [-1e-4:1e-4]

set pm3d map
unset xtics
unset ytics
unset key
unset grid
set samples 2
set isosamples 2

# point size
set pointsize 0.2

set output "OUTPUT_FILES/rho_and_kappa_kernels.png"

set multiplot

set size 0.5,0.5
set origin 0,0
set title "rho kernel acoustic only"
splot "OUTPUT_FILES/proc000000_rho_kappa_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? 1/0 : $3) w points palette pt 5

set size 0.5,0.5
set origin 0.47,0
set title "kappa kernel acoustic only"
splot "OUTPUT_FILES/proc000000_rho_kappa_kernel.dat" using 1:2:(($1>=2000.0 &&$1<=3000.0)? 1/0 : $4)   w points palette pt 5

set size 0.5,0.5
set origin 0,0.5
set title "rho kernel acoustic and elastic"
splot "OUTPUT_FILES/proc000000_rho_kappa_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? 1/0 : $3) w points palette pt 5,"OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? $3 : 1/0) w points palette pt 5

set size 0.5,0.5
set origin 0.47,0.5
set title "kappa kernel acoustic and elastic"
splot "OUTPUT_FILES/proc000000_rho_kappa_kernel.dat" using 1:2:(($1>=2000.0 &&$1<=3000.0)? 1/0 : $4)   w points palette pt 5,"OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? $4 : 1/0) w points palette pt 5

unset multiplot

set print "-"
print ""
print "gnuplot plotted: OUTPUT_FILES/rho_and_kappa_kernels.png"


############## add the drawing of the first two pictures as separate picture and without any caption
############## in order to use them as a comparison and validation image to check adjoint runs and kernel calculations in BuildBot

unset xtics
unset ytics
unset key
unset grid
unset cbtics
unset title
unset colorbox
unset border
unset tics

# for output compatibility linux/mac/.. to compare images
# sets margins <left>, <right>, <bottom>, <top> to place two figures in the middle of the canvas (or screen)
# gnuplot versions >= 5
#set margins -1, -1, screen 0.25, screen 0.75
# gnuplot versions < 5
set bmargin screen 0.25
set tmargin screen 0.75

set output "OUTPUT_FILES/rho_and_kappa_kernels_acoustic_and_elastic_only_no_caption.png"

set multiplot

set lmargin screen 0.0
set rmargin screen 0.5

set size 0.5,0.5
set origin 0,0.5
splot "OUTPUT_FILES/proc000000_rho_kappa_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? 1/0 : $3) w points palette pt 5,"OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? $3 : 1/0) w points palette pt 5

set lmargin screen 0.5
set rmargin screen 1.0

set size 0.5,0.5
set origin 0.5,0.5
splot "OUTPUT_FILES/proc000000_rho_kappa_kernel.dat" using 1:2:(($1>=2000.0 &&$1<=3000.0)? 1/0 : $4)   w points palette pt 5,"OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? $4 : 1/0) w points palette pt 5

unset multiplot

#debug view
#show variables all

set print "-"
print "gnuplot plotted: OUTPUT_FILES/rho_and_kappa_kernels_acoustic_and_elastic_only_no_caption.png"
print ""

