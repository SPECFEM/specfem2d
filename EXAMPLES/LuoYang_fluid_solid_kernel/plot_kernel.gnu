set term png

set output "rho_and_kappa_kernels.png"

set palette defined ( -8e-15 "red", 0 "yellow", 8e-15 "blue")
set pm3d map
unset xtics
unset ytics
unset key
unset grid
set samples 2
set cbrange [-8e-15:8e-15]
set isosamples 2

set multiplot

set size 0.5,0.5
set origin 0,0
set title "rho kernel acoustic only"
splot "OUTPUT_FILES/proc000000_rho_kappa_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? 1/0 : $3) w points palette ps 0.02 pt 5

set size 0.5,0.5
set origin 0.47,0
set title "kappa kernel acoustic only"
splot "OUTPUT_FILES/proc000000_rho_kappa_kernel.dat" using 1:2:(($1>=2000.0 &&$1<=3000.0)? 1/0 : $4)   w points palette ps 0.02 pt 5

set size 0.5,0.5
set origin 0,0.5
set title "rho kernel acoustic and elastic"
splot "OUTPUT_FILES/proc000000_rho_kappa_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? 1/0 : $3) w points palette ps 0.02 pt 5,"OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? $3 : 1/0) w points palette ps 0.02 pt 5

set size 0.5,0.5
set origin 0.47,0.5
set title "kappa kernel acoustic and elastic"
splot "OUTPUT_FILES/proc000000_rho_kappa_kernel.dat" using 1:2:(($1>=2000.0 &&$1<=3000.0)? 1/0 : $4)   w points palette ps 0.02 pt 5,"OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? $4 : 1/0) w points palette ps 0.02 pt 5

unset multiplot

############## add the drawing of the first two pictures as separate picture and without any caption
############## in order to use them as a comparison and validation image to check adjoint runs and kernel calculations in BuildBot

unset xtics
unset ytics
unset key
unset grid
unset cbtics
unset title

set output "rho_and_kappa_kernels_acoustic_and_elastic_only_no_caption.png"

set multiplot

set size 0.5,0.5
set origin 0,0.5
splot "OUTPUT_FILES/proc000000_rho_kappa_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? 1/0 : $3) w points palette ps 0.02 pt 5,"OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? $3 : 1/0) w points palette ps 0.02 pt 5

set size 0.5,0.5
set origin 0.5,0.5
splot "OUTPUT_FILES/proc000000_rho_kappa_kernel.dat" using 1:2:(($1>=2000.0 &&$1<=3000.0)? 1/0 : $4)   w points palette ps 0.02 pt 5,"OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat" using 1:2:(($1>=2000.0 && $1<=3000.0)? $4 : 1/0) w points palette ps 0.02 pt 5

unset multiplot

