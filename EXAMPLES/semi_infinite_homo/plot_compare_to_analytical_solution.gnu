
set term x11

set xlabel "Time (s)"
set ylabel "Amplitude of displacement component (m)"

set xrange [-0.15:0.7]

plot "OUTPUT_FILES/S0001.AA.BXX.semd" t 'Numerical Ux' w l lc 1
pause -1 "Hit any key..."

plot "Ux_file.dat" t 'Analytic Ux' w l lc 3
pause -1 "Hit any key..."

plot "OUTPUT_FILES/S0001.AA.BXX.semd" t 'Numerical Ux' w l lc 1, "Ux_file.dat" t 'Analytic Ux' w l lc 3
pause -1 "Hit any key..."

plot "OUTPUT_FILES/S0001.AA.BXZ.semd" t 'Numerical Uz' w l lc 1
pause -1 "Hit any key..."

plot "Uz_file.dat" t 'Analytic Uz' w l lc 3
pause -1 "Hit any key..."

plot "OUTPUT_FILES/S0001.AA.BXZ.semd" t 'Numerical Uz' w l lc 1, "Uz_file.dat" t 'Analytic Uz' w l lc 3
pause -1 "Hit any key..."


#plot "Ux_file.dat" t 'Analytic Ux' w l lc 3
#pause -1 "Hit any key..."

#plot "OUTPUT_FILES/S0001.AA.BXX.semd" t 'Numerical Ux' w l lc 1
#pause -1 "Hit any key..."

#plot "OUTPUT_FILES/S0001.AA.BXZ.semd" t 'Numerical Uz' w l lc 1
#pause -1 "Hit any key..."

#plot "source.txt" t 'Numerical Uz' w l lc 1 
#pause -1 "Hit any key..."

#plot "OUTPUT_FILES/S0001.AA.BXX.semd" t 'Numerical Ux' w l lc 1, "S0001.AA.BXX.rk.semd" t 'RKUx' w l lc 3
#pause -1 "Hit any key..."
#
#plot "OUTPUT_FILES/S0001.AA.BXZ.semd" t 'Numerical Uz' w l lc 1, "S0001.AA.BXZ.rk.semd" t 'RKUw' w l lc 3
#pause -1 "Hit any key..."
