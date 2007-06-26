#! /bin/sh


# This is a script for submitting serial jobs on the UPPA cluster.

# BEWARE : by default, old Databases along with new Databases in ./OUTPUT_FILES directory are removed;
# comment "rm ./OUTPUT_FILES/Database*" in this file and in qsub_UPPA.sh if you wish to change this behavior.

# note : xmeshfem2D DOES NOT run on a node, its runs on a frontal and outside of SGE. It was done so as to provide a way to
# submit jobs when the Databases were previously created.


rm ./OUTPUT_FILES/Database*

./xmeshfem2D

qsub ./qsub_UPPA_serial.sh

