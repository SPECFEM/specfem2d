
#set term x11
set term qt

set term pdf color solid
set output "comparison_with_analytical_solution.pdf"

set xrange [1.3:2.1]

# in the original paper from the year 2000 https://specfem.github.io/komatitsch.free.fr/published_papers/geophysics_fluid_solid_Dirichlet_added.pdf
# we did not subtract the t0 of the source from the seismograms, and t0 is 0.104 for this test, thus we add it here
# in the Gnuplot script because the current version of SPECFEM2D subtracts it from seismograms automatically

# the sign of the horizontal component is also inverted compared to the analytical solution, which was for a receiver
# located on the other side compared to the horizontal position of the source, thus we multiply the horizontal seismogram by -1 here

# We also used a different absolute amplitude for the source, but this does not matter because the wave equation is linear
# and thus absolute amplitude does not really matter in this test; we thus multiply by this constant factor here in the script

amplitude_of_source = 44705882000.

plot "OUTPUT_FILES/AA.S0040.BXX.semv" using ($1 + 0.104):(- amplitude_of_source * $2) title 'Ux SEM' w l lc 1, "trace2an.dat" using 1:2 title 'Ux analytical' w l lc 3
pause -1 "Hit any key..."

plot "OUTPUT_FILES/AA.S0040.BXZ.semv" using ($1 + 0.104):(+ amplitude_of_source * $2) title 'Uz SEM' w l lc 1, "trace2an.dat" using 1:3 title 'Uz analytical' w l lc 3
pause -1 "Hit any key..."

