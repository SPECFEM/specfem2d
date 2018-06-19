
set term x11

set xlabel "Time (s)"
set ylabel "Amplitude of displacement component (m)"

set xrange [0:0.8]

# the analytical code from Denmark uses depth as the Z axis, i.e. its X axis points right but its Z axis points down rather than up.
# thus, to compare to the standard reference, one needs to invert the sign of Ux (but leave the sign of Uz unchanged)
factor_to_invert_Ux_Denmark = -1

plot "OUTPUT_FILES/AA.S0001.BXX.semd" t 'Numerical Ux elastic' w l lc 5, "Ux_time_analytical_solution_elastic_time_domain.dat" t 'Ux my analytical code in the time domain' w l lc 4, "Ux_file_analytical_from_Denmark.dat" using 1:(factor_to_invert_Ux_Denmark*$2) t 'Ux other analytical code from Denmark' w l lc 1, "Ux_time_analytical_solution_elastic.dat" t 'Quasi-analytical Ux elastic Carcione in frequency' w l lc 3
pause -1 "Hit any key..."

plot "OUTPUT_FILES/AA.S0001.BXZ.semd" t 'Numerical Uz elastic' w l lc 5, "Uz_time_analytical_solution_elastic_time_domain.dat" t 'Uz my analytical code in the time domain' w l lc 4, "Uz_file_analytical_from_Denmark.dat" t 'Uz other analytical code from Denmark' w l lc 1, "Uz_time_analytical_solution_elastic.dat" t 'Quasi-analytical Uz elastic Carcione in frequency' w l lc 3
pause -1 "Hit any key..."

