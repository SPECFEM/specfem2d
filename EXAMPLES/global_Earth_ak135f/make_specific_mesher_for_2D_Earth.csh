#!/bin/csh

rm -f xcreate_mesh_files

#ifort -o xcreate_mesh_files -O0 -vec-report0 -std03 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check all -align sequence -assume byterecl -ftz -traceback -ftrapuv create_mesh_AK135F_2D_with_central_cube_no_PML.F90
#ifort -o xcreate_mesh_files -O3 -xHost -vec-report0 -std03 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds -mcmodel=medium -shared-intel create_mesh_AK135F_2D_with_central_cube_no_PML.F90
gfortran -o xcreate_mesh_files -O3 create_mesh_AK135F_2D_with_central_cube_no_PML.F90

./xcreate_mesh_files

rm -f xcreate_mesh_files

