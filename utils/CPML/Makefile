
# Fortran90 compiler
F90 = gfortran
#F90 = ifort

FFLAGS = -O3

# for Intel ifort
#FFLAGS = -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -std08 -diag-disable 6477 -implicitnone -gen-interfaces -warn all -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv

.NOTPARALLEL:

default: xconvert_external_layers_of_a_given_mesh_to_CPML_layers xadd_CPML_layers_to_an_existing_mesh cleanobj

all: clean default

xconvert_external_layers_of_a_given_mesh_to_CPML_layers: convert_external_layers_of_a_given_mesh_to_CPML_layers.o
	 ${F90} $(FFLAGS) -o xconvert_external_layers_of_a_given_mesh_to_CPML_layers convert_external_layers_of_a_given_mesh_to_CPML_layers.o

xadd_CPML_layers_to_an_existing_mesh: add_CPML_layers_to_an_existing_mesh.o
	 ${F90} $(FFLAGS) -o xadd_CPML_layers_to_an_existing_mesh add_CPML_layers_to_an_existing_mesh.o

clean:
	rm -rf *.o *__genmod.*90 *.mod xconvert_external_layers_of_a_given_mesh_to_CPML_layers xadd_CPML_layers_to_an_existing_mesh

cleanobj:
	rm -rf *.o *__genmod.*90 *.mod

# OBJECTS:

convert_external_layers_of_a_given_mesh_to_CPML_layers.o: convert_external_layers_of_a_given_mesh_to_CPML_layers.f90
	${F90} $(FFLAGS) -c -o convert_external_layers_of_a_given_mesh_to_CPML_layers.o convert_external_layers_of_a_given_mesh_to_CPML_layers.f90

add_CPML_layers_to_an_existing_mesh.o: add_CPML_layers_to_an_existing_mesh.f90
	${F90} $(FFLAGS) -c -o add_CPML_layers_to_an_existing_mesh.o add_CPML_layers_to_an_existing_mesh.f90

