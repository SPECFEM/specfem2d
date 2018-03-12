#!/bin/bash
#

./run_this_example.sh # Run specfem2d
eog OUTPUT_FILES/image0001800.jpg &
./plot.py OUTPUT_FILES/AA.S0001.PRE.semp OUTPUT_FILES/AA.S0003.PRE.semp --hold -n --shift -0.05 --xlabel "Shifted and normalized signals at source and first receiver" &
./plot.py OUTPUT_FILES/AA.S0003.PRE.semp OUTPUT_FILES/AA.S0002.PRE.semp --hold -n --shift -0.05 --xlabel "Shifted and normalized signals at first and second receiver" &

