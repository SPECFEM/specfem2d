#!/bin/bash
#
# This script make the FFTs in order to benchmark specfem2d with COMSOL and ROTVARS.
# It compute transmission losses for specfem seismograms and compare them to the ones computed with COMSOL
# This script has to be run from the directory specfem2d/EXAMPLES/paper_axisymmetry_example
#

factor=0.00000000000055 # TODO change that to obtain a good fit, it is just a multiplying factor

# Run specfem2d:
#./run_this_example.sh

# Compute Transmission losses using the signals at the receivers:
./computeTL.py OUTPUT_FILES/ --freq 5.0 --NdownSampled 4 --outputPath OUTPUT_FILES/lossesGOOD --noplot
#/home/bottero/bin/mpiexec -n 12 python ./computeTLmpi.py OUTPUT_FILES/ --freq 5.0 --NdownSampled 4 --outputPath OUTPUT_FILES/lossesGOOD --noplot

# Convert kilometers to meters:
awk '{ print $1/1000," ",$2}' TLcomsol.txt > TLcomsolKM.txt

# Scale amplitudes:
./processAmplitudes.py OUTPUT_FILES/lossesGOODampl --ref $factor --noplot

# Shift range (rotvars starts at range >0): 
awk '{ print $1-0.03386," ",-$2}' rotvars_Jensen_2007_Figure7a.txt > rotvars_Jensen_2007_Figure7a2.txt

# Plot comparison
../../utils/Visualization/plot.py rotvars_Jensen_2007_Figure7a2.txt TLcomsolKM.txt OUTPUT_FILES/lossesGOODamplLosses --hold -w 3 --fontsize 16 --xlabel "Range (km)" --ylabel "Transmission losses (dB)" --invert_yaxis -l 'ROTVARS','COMSOL','Spectral elements' -c 0,0,0.8+0,0.8,0+0.8,0,0
