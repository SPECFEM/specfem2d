#!/bin/bash

# define your Fortran compiler here
FC="gfortran"

# forward simulation
./change_simulation_type.pl -F
./run_a_single_run_of_this_two_step_example.sh

# adjoint source
rm -f xadj_source
$FC adj_source.f90 -o xadj_source
./xadj_source
rm -f xadj_source

# adjoint kernel simulation
./change_simulation_type.pl -b
./run_a_single_run_of_this_two_step_example.sh do_not_erase_disk_files

# kernel images
gnuplot plot_kernel.gnu

