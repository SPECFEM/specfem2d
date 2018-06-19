
set term x11
set term qt

set xlabel "Time (s)"
set ylabel "Amplitude of displacement component (m)"

set xrange [-0.1:0.8]

plot "pressure_time_analytical_solution_acoustic_time_domain_0point1m.dat" w l lc 4, "pressure_time_analytical_solution_acoustic_0point1m.dat" w l lc 3
pause -1 "Hit any key..."

plot "pressure_time_analytical_solution_acoustic_time_domain_1m.dat" w l lc 4, "pressure_time_analytical_solution_acoustic_1m.dat" w l lc 3
pause -1 "Hit any key..."

set xrange [0:0.8]

plot "pressure_time_analytical_solution_acoustic_time_domain_500m.dat" w l lc 4, "pressure_time_analytical_solution_acoustic_500m.dat" w l lc 3
pause -1 "Hit any key..."




####################### plot "pressure_time_analytical_solution_acoustic_time_domain_0point1m.dat" using ($1+.21441454147432084501):(factor_to_invert_Ux_Denmark*$2/fact) w l lc 1, "pressure_time_analytical_solution_acoustic_time_domain_1m.dat" using ($1+.21441454147432084501):(factor_to_invert_Ux_Denmark*$2/fact) w l lc 2, "pressure_time_analytical_solution_acoustic_time_domain_500m.dat" using 1:(factor_to_invert_Ux_Denmark*$2) w l lc 4, "pressure_time_analytical_solution_acoustic_500_500m.dat" using 1:(factor_to_invert_Ux_Denmark*$2) w l lc 3
####################### pause -1 "Hit any key..."

####################### plot "OUTPUT_FILES_inside/AA.S0001.PRE.semp" t 'Numerical Ux acoustic inside' w l lc 5, "OUTPUT_FILES_corner/AA.S0001.PRE.semp" t 'Numerical Ux acoustic corner' w l lc 3, "pressure_time_analytical_solution_acoustic_500.dat" using 1:(factor_to_invert_Ux_Denmark*$2) w l lc 1
####################### pause -1 "Hit any key..."

####################### plot "OUTPUT_FILES_inside/AA.S0001.PRE.semp" t 'Numerical Ux acoustic inside' w l lc 5, "OUTPUT_FILES_corner/AA.S0001.PRE.semp" t 'Numerical Ux acoustic corner' w l lc 3
####################### pause -1 "Hit any key..."

