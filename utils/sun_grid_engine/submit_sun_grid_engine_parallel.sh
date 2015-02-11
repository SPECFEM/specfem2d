#! /bin/sh


# This is a script for submitting parallel jobs on the UPPA cluster.
# DO NOT FORGET to check file qsub_openmpi_UPPA.sh to change the number of
# proccesses according to the value specified in Par_file.
# The path to OpenMPI libs should or should not be specified, depending
# on the user environement.

# The allocation rule for parallel jobs is round robin when using "openmpi" PE;
# specify "openmpi_fillup" to have the fill up allocation rule.
# "openmpi_fillup" should reduce the communication overhead (more intra-nodal communications) with most MPI libraries.
# "openmpi" should leave a few CPUs per node to manage the system without hampering the run. You should try both, and decide depending
# on your code and setup.

# BEWARE : by default, old Databases along with new Databases in ./OUTPUT_FILES directory are removed;
# comment "rm ./OUTPUT_FILES/Database*" in this file and in qsub_openmpi_UPPA.sh if you wish to change this behavior.

# note : xmeshfem2D DOES NOT run on a node, its runs on a frontal and outside of SGE. It was done so as to provide a way to
# submit jobs when the Databases were previously created.


rm ./OUTPUT_FILES/Database*

LD_LIBRARY_PATH=/opt/openmpi-1.2.2/pgi64/lib ./xmeshfem2D

qsub -pe openmpi_fillup 8 ./qsub_UPPA_parallel.sh

