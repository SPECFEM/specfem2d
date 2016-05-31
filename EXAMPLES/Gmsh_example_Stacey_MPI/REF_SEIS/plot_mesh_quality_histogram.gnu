 set term wxt
 #set term gif
 #set output "mesh_quality_histogram.gif"

 set xrange [0:1]
 set xtics 0,0.1,1
 set boxwidth    5.00000007E-02
 set xlabel "Skewness range"
 set ylabel "Percentage of elements (%)"
 set loadpath "./OUTPUT_FILES"
 plot "mesh_quality_histogram.txt" with boxes
 pause -1 "hit any key..."
