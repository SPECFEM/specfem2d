
set term x11

set xrange [0:0.600]


# this second factor fixes a -1.e+3 unit difference between the two analytical solutions (from Carcione and from Gar6more2D) for some reason
factor = - 1.e+3

# onset time of the source (1.2 / f0, with f0 = 50 Hz)
t0 = .024000

plot "no_attenuation_analytical/pressure_time_analytical_solution_acoustic_800.dat" title 'Quasi-analytical pressure Carcione' w l lc 3, "Gar6more2D_without_attenuation/P_Garcimore_at_800m.dat" using ($1 - t0):(factor*$2) title 'Quasi-analytical pressure Gar6more2D' w l lc 4, "OUTPUT_FILES/AA.S0001.PRE.semp" title 'SPECFEM2D pressure' w l lc 1

pause -1 "Hit any key..."


