
set term x11
#set term qt

set xlabel "Time (s)"
set ylabel "Pressure"

set xrange [0:0.8]

plot "pressure_time_analytical_solution_acoustic_time_domain.dat" w l lc 4, "pressure_time_analytical_solution_acoustic_500.dat" w l lc 3, "OUTPUT_FILES/AA.S0001.PRE.semp" w l lc 1
pause -1 "Hit any key..."

