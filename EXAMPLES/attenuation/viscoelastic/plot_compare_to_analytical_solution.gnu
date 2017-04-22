
set term x11

set xlabel "Time (s)"
set ylabel "Amplitude of displacement component (m)"

set xrange [0:1.4]

plot "OUTPUT_FILES/AA.S0001.BXX.semd" t 'Numerical Ux' w l lc 1, "Ux_time_analytical_solution_viscoelastic.dat" t 'Quasi-analytical Ux' w l lc 3
pause -1 "Hit any key..."

plot "OUTPUT_FILES/AA.S0001.BXZ.semd" t 'Numerical Uz' w l lc 1, "Uz_time_analytical_solution_viscoelastic.dat" t 'Quasi-analytical Uz' w l lc 3
pause -1 "Hit any key..."

