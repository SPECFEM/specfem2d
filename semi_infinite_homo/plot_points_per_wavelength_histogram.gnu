 set term x11
 #set term gif
 #set output "points_per_wavelength_histogram_P_in_fluid.gif"
 
 set boxwidth   1.5000000E-02
 set xlabel "Range of min number of points per P wavelength in fluid"
 set ylabel "Percentage of elements (%)"
 plot "points_per_wavelength_histogram_P_in_fluid.txt" with boxes
 pause -1 "hit any key..."
