#!/bin/bash

# define your Fortran compiler here
FC="gfortran"
FLAGS=-Wall

## forward simulation
./change_simulation_type.pl -F

./run_a_single_run_of_this_two_step_example.sh
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# adjoint source
echo
echo "creating adjoint source..."
echo
rm -f xadj_source
$FC $FLAGS adj_source.f90 -o xadj_source
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# clean adjoint source
mkdir -p SEM
rm -f SEM/*
echo
echo "***************************************"
echo
./xadj_source
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
echo
echo "***************************************"
echo

rm -f xadj_source

## adjoint kernel simulation
./change_simulation_type.pl -b

./run_a_single_run_of_this_two_step_example.sh do_not_erase_disk_files
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

## kernel images
gnuplot plot_kernel.gnu

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi
