
set term x11

set xlabel "Time (s)"
set ylabel "Amplitude of displacement component (m)"

set xrange [-0.15:0.7]

plot "OUTPUT_FILES/AA.S0001.BXX.semd" t 'Numerical Ux' w l lc 1, "x_displacement.dat" t 'Analytic Ux' w l lc 2
pause -1 "Hit any key..."

plot "OUTPUT_FILES/AA.S0001.BXZ.semd" t 'Numerical Uz' w l lc 1, "z_displacement.dat" t 'Analytic Uz' w l lc 2
pause -1 "Hit any key..."

