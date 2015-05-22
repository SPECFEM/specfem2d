#!/bin/bash

# define your Fortran compiler here
FC="gfortran"

cp -f Par_file_forward Par_file
./run_a_single_run_of_this_two_step_example.sh

rm -f xadj_source
$FC adj_source.f90 -o xadj_source
./xadj_source
rm -f xadj_source

cp -f Par_file_adjoint Par_file
./run_a_single_run_of_this_two_step_example.sh do_not_erase_disk_files

gnuplot plot_kernel.gnu

