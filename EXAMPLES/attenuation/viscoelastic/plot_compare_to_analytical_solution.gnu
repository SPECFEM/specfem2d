
set term x11

set xlabel "Time (s)"
set ylabel "Amplitude of displacement component (m)"

set xrange [0:0.8]

plot "REF_SEIS/elastic_AA.S0001.BXX.semd" t 'Numerical Ux elastic' w l lc 4, "Ux_time_analytical_solution_elastic.dat" t 'Quasi-analytical Ux elastic' w l lc 5, "OUTPUT_FILES/AA.S0001.BXX.semd" t 'Numerical Ux viscoelastic' w l lc 1, "Ux_time_analytical_solution_viscoelastic.dat" t 'Quasi-analytical Ux viscoelastic' w l lc 3
pause -1 "Hit any key..."

plot "REF_SEIS/elastic_AA.S0001.BXZ.semd" t 'Numerical Uz elastic' w l lc 4, "Uz_time_analytical_solution_elastic.dat" t 'Quasi-analytical Uz elastic' w l lc 5, "OUTPUT_FILES/AA.S0001.BXZ.semd" t 'Numerical Uz viscoelastic' w l lc 1, "Uz_time_analytical_solution_viscoelastic.dat" t 'Quasi-analytical Uz viscoelastic' w l lc 3
pause -1 "Hit any key..."

