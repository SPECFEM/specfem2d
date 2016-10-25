#!/bin/bash
#
######################################################
## INPUT parameters

# type (=0 kernel or =1 wavefield)
type=0

# title PSV (=1 P-SV or =0 for SH)
psv=1

# normalization factor (divides kernel values by 10**norm)
norm=-8

#######################################################

# plotting script
rm -f plot_wavefield.pl
ln -s ../../utils/Visualization/plot_wavefield.pl


# plots kernels
#
# note: script by default plots kernels K_(rho/vp/vs)
./plot_wavefield.pl 400/2800/400 400/2000/400 0/200/0/80 50/10/40/10 $norm/$norm/$norm 1/1/1  -8.0/0.02 1/0/1/$psv/$type 4.0/1/0 1/0/1/200 Tromp2005_kernel PSV_homo
