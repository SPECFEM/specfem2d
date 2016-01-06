
set term x11

set xlabel "Time (s)"
set ylabel "Amplitude of displacement component (m)"

set xrange [-0.15:0.7]

plot "OUTPUT_FILES/AA.S0001.BXX.semd" t 'Numerical Ux' w l lc 1
pause -1 "Hit any key..."

plot "Ux_file.dat" t 'Analytic Ux' w l lc 3
pause -1 "Hit any key..."

plot "OUTPUT_FILES/AA.S0001.BXX.semd" t 'Numerical Ux' w l lc 1, "Ux_file.dat" t 'Analytic Ux' w l lc 3
pause -1 "Hit any key..."

plot "OUTPUT_FILES/AA.S0001.BXZ.semd" t 'Numerical Uz' w l lc 1
pause -1 "Hit any key..."

plot "Uz_file.dat" t 'Analytic Uz' w l lc 3
pause -1 "Hit any key..."

plot "OUTPUT_FILES/AA.S0001.BXZ.semd" t 'Numerical Uz' w l lc 1, "Uz_file.dat" t 'Analytic Uz' w l lc 3
pause -1 "Hit any key..."


#plot "Ux_file.dat" t 'Analytic Ux' w l lc 3
#pause -1 "Hit any key..."

#plot "OUTPUT_FILES/AA.S0001.BXX.semd" t 'Numerical Ux' w l lc 1
#pause -1 "Hit any key..."

#plot "OUTPUT_FILES/AA.S0001.BXZ.semd" t 'Numerical Uz' w l lc 1
#pause -1 "Hit any key..."

#plot "OUTPUT_FILES/plot_source_time_function.txt" t 'Numerical Uz' w l lc 1
#pause -1 "Hit any key..."

#plot "OUTPUT_FILES/AA.S0001.BXX.semd" t 'Numerical Ux' w l lc 1, "AA.S0001.BXX.rk.semd" t 'RKUx' w l lc 3
#pause -1 "Hit any key..."
#
#plot "OUTPUT_FILES/AA.S0001.BXZ.semd" t 'Numerical Uz' w l lc 1, "AA.S0001.BXZ.rk.semd" t 'RKUw' w l lc 3
#pause -1 "Hit any key..."
