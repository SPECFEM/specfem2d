#!/opt/local/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.0 patchlevel 3    last modified 2016-02-21 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2016
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
## set terminal x11  nopersist enhanced font "Arial,10" size 800,400 
# set output
#############################################

set terminal png 
set output 'image.png'

#############################################

set xlabel 'Time (s)'
set ylabel 'vertical velocity'

set yrange [ -0.165 : 0.165 ] 
set xtics 0.04

plot [-0.04:0.2] 'OUTPUT_FILES/AA.S0001.BXZ.semv' w l
