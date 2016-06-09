set term png size 1200,400

set output "image_rho_kappa_mu_kernels.png"

set palette defined ( -8e-15 "red", 0 "yellow", 8e-15 "blue")
set pm3d map
unset xtics
unset ytics
unset key
unset grid
set samples 2
set cbrange [-1.e-8:1.e-8]
set isosamples 2

set multiplot

set size 0.333,1
set origin 0,0
set title "rho kernel"
splot "OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat" using 1:2:3 w points palette pt 7 ps 0.5

set size 0.333,1
set origin 0.33,0
set title "kappa kernel"
splot "OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat" using 1:2:4 w points palette pt 7 ps 0.5

set size 0.333,1
set origin 0.64,0
set title "mu kernel"
splot "OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat" using 1:2:5 w points palette pt 7 ps 0.5

unset multiplot


