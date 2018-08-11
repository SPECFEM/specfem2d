
set term x11
set term qt

set xlabel "Time (s)"
set ylabel "Pressure"

f0 = 18.
factor = - 1. / (2 * 3.141592653589793 * 3.141592653589793 * f0 * f0)

set xrange [0:0.8]

plot "pressure_time_analytical_solution_acoustic_time_domain.dat" w l lc 4, "pressure_time_analytical_solution_acoustic_500.dat" w l lc 3, "OUTPUT_FILES/AA.S0001.PRE.semp" using 1:(factor*$2) w l lc 1
pause -1 "Hit any key..."

