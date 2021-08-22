
set term x11

set xrange [0:0.600]


plot "attenuation_viscoacoustic_NSLS_3/pressure_time_analytical_solution_viscoacoustic_800.dat" title 'Quasi-analytical pressure Carcione NSLS = 3' w l lc 3, "OUTPUT_FILES/AA.S0001.PRE.semp" title 'SPECFEM2D pressure with viscoacoustic' w l lc 2, "../viscoacoustic_attenuation_off_versus_analytical/OUTPUT_FILES/AA.S0001.PRE.semp" title 'SPECFEM2D pressure without viscoacoustic' w l lc 4
pause -1 "Hit any key..."

plot "attenuation_viscoacoustic_NSLS_3/pressure_time_analytical_solution_viscoacoustic_800.dat" title 'Quasi-analytical pressure Carcione NSLS = 3' w l lc 3, "OUTPUT_FILES/AA.S0001.PRE.semp" title 'SPECFEM2D pressure with viscoacoustic' w l lc 2
pause -1 "Hit any key..."

