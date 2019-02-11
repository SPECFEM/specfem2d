#!/bin/csh

rm -f bin/xcreate_mesh_files

#ifort -o xcreate_mesh_files -O0 -vec-report0 -std03 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check all -align sequence -assume byterecl -ftz -traceback -ftrapuv create_mesh_AK135F_2D_with_central_cube_no_PML.F90
#ifort -o xcreate_mesh_files -O3 -xHost -vec-report0 -std03 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds create_mesh_AK135F_2D_with_central_cube_no_PML.F90
gfortran -o bin/xcreate_mesh_files -O3 create_mesh_AK135F_2D_with_central_cube_no_PML.F90

./bin/xcreate_mesh_files
# checks exit code
if ( $? != 0 ) then
  exit 1
endif

rm -f bin/xcreate_mesh_files

