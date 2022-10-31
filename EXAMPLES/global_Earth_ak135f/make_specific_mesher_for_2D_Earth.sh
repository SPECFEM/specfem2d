#!/bin/bash

# gets the default compiler
f90=`grep "^FC " ../../Makefile | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`


## intel debug flags
#flags="-O0 -vec-report0 -std03 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check all -align sequence -assume byterecl -ftz -traceback -ftrapuv"
## intel run flags
#flags="-O3 -xHost -vec-report0 -std03 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds"
## default flags
flags="-O3"

# clean
rm -f bin/xcreate_mesh_files

# compile
$f90 $flags -o bin/xcreate_mesh_files create_mesh_AK135F_2D_with_central_cube_no_PML.F90

./bin/xcreate_mesh_files
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

#rm -f bin/xcreate_mesh_files

echo
echo "done w/ script: make_specific_mesher_for_2D_Earth.sh"
echo

