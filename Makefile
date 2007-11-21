#
# Makefile for SPECFEM2D version 5.2
#
# Dimitri Komatitsch, University of Pau, April 2007
#
SHELL=/bin/sh

O = obj

# Portland
#F90 = /opt/openmpi-1.2.2/pgi64/bin/mpif90 -DUSE_MPI -DUSE_METIS -DUSE_SCOTCH
F90 = pgf90
CC = pgcc
FLAGS_NOCHECK=-fast -Mnobounds -Minline -Mneginfo -Mdclchk -Knoieee -Minform=warn -fastsse -tp amd64e -Msmart
FLAGS_CHECK=-fast -Mbounds -Mneginfo -Mdclchk -Minform=warn

# Intel
#F90 = ifort
#CC = gcc
#FLAGS_NOCHECK=-O3 -implicitnone -warn stderrors -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -assume byterecl -check nobounds
#FLAGS_CHECK = $(FLAGS_NOCHECK) -check bounds

# GNU gfortran
#F90 = /opt/openmpi-1.2.1/gfortran64/bin/mpif90 -DUSE_MPI -DUSE_METIS -DUSE_SCOTCH
#F90 = gfortran
#CC = gcc
#FLAGS_NOCHECK = -O3 -march=opteron -m64 -mfpmath=sse,387
#FLAGS_NOCHECK = -std=gnu -fimplicit-none -frange-check -O2 -Wunused-labels -Waliasing -Wampersand -Wsurprising -Wline-truncation -Wunderflow
#FLAGS_CHECK = $(FLAGS_NOCHECK) -fbounds-check

LINK = $(F90)

#LIB = /opt/metis-4.0/gcc64/lib/libmetis.a /opt/scotch-4.0/gcc64/lib/libscotch.a  /opt/scotch-4.0/gcc64/lib/libscotcherr.a
LIB = 

OBJS_MESHFEM2D = $O/part_unstruct.o $O/meshfem2D.o $O/read_value_parameters.o $O/numerical_recipes.o

OBJS_SPECFEM2D = $O/checkgrid.o $O/datim.o $O/enforce_acoustic_free_surface.o\
        $O/compute_forces_acoustic.o $O/compute_forces_elastic.o\
        $O/lagrange_poly.o $O/gmat01.o $O/gll_library.o $O/plotgll.o $O/define_derivation_matrices.o\
        $O/plotpost.o $O/locate_receivers.o $O/locate_source_force.o $O/compute_gradient_attenuation.o\
        $O/specfem2D.o $O/write_seismograms.o $O/define_external_model.o $O/createnum_fast.o $O/createnum_slow.o\
        $O/define_shape_functions.o $O/attenuation_model.o $O/create_color_image.o $O/compute_vector_field.o $O/compute_pressure.o\
        $O/recompute_jacobian.o $O/compute_arrays_source.o $O/locate_source_moment_tensor.o $O/netlib_specfun_erf.o\
        $O/construct_acoustic_surface.o $O/assemble_MPI.o $O/compute_energy.o\
        $O/attenuation_compute_param.o

default: clean meshfem2D specfem2D convolve_source_timefunction

all: default

clean:
	/bin/rm -r -f xmeshfem2D xmeshfem2D.trace xspecfem2D xspecfem2D.trace $O/*.o *.o $O/*.il *.mod core xconvolve_source_timefunction *.oo *.ipo

meshfem2D: $(OBJS_MESHFEM2D)
	$(LINK) $(FLAGS_CHECK) -o xmeshfem2D $(OBJS_MESHFEM2D) $(LIB)

### use optimized compilation option for solver only
specfem2D: $(OBJS_SPECFEM2D)
	$(LINK) $(FLAGS_NOCHECK) -o xspecfem2D $(OBJS_SPECFEM2D)

convolve_source_timefunction: $O/convolve_source_timefunction.o
	${F90} $(FLAGS_CHECK) -o xconvolve_source_timefunction $O/convolve_source_timefunction.o

$O/checkgrid.o: checkgrid.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/checkgrid.o checkgrid.F90
    
$O/meshfem2D.o: meshfem2D.F90
	${F90} $(FLAGS_CHECK) -c -o $O/meshfem2D.o meshfem2D.F90

$O/createnum_fast.o: createnum_fast.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/createnum_fast.o createnum_fast.f90
    
$O/createnum_slow.o: createnum_slow.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/createnum_slow.o createnum_slow.f90
    
$O/convolve_source_timefunction.o: convolve_source_timefunction.f90
	${F90} $(FLAGS_CHECK) -c -o $O/convolve_source_timefunction.o convolve_source_timefunction.f90

$O/read_value_parameters.o: read_value_parameters.f90
	${F90} $(FLAGS_CHECK) -c -o $O/read_value_parameters.o read_value_parameters.f90

$O/datim.o: datim.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/datim.o datim.f90
    
$O/lagrange_poly.o: lagrange_poly.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/lagrange_poly.o lagrange_poly.f90
    
$O/gmat01.o: gmat01.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/gmat01.o gmat01.f90
    
$O/gll_library.o: gll_library.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/gll_library.o gll_library.f90
    
$O/define_derivation_matrices.o: define_derivation_matrices.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/define_derivation_matrices.o define_derivation_matrices.f90
    
$O/plotgll.o: plotgll.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/plotgll.o plotgll.f90
    
$O/plotpost.o: plotpost.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/plotpost.o plotpost.F90
    
$O/locate_receivers.o: locate_receivers.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/locate_receivers.o locate_receivers.F90
    
$O/recompute_jacobian.o: recompute_jacobian.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/recompute_jacobian.o recompute_jacobian.f90
    
$O/locate_source_force.o: locate_source_force.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/locate_source_force.o locate_source_force.F90
    
$O/locate_source_moment_tensor.o: locate_source_moment_tensor.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/locate_source_moment_tensor.o locate_source_moment_tensor.F90
    
$O/define_shape_functions.o: define_shape_functions.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/define_shape_functions.o define_shape_functions.f90
    
$O/attenuation_model.o: attenuation_model.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/attenuation_model.o attenuation_model.f90
    
### use optimized compilation option for solver only
$O/specfem2D.o: specfem2D.F90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/specfem2D.o specfem2D.F90
    
### use optimized compilation option for solver only
$O/enforce_acoustic_free_surface.o: enforce_acoustic_free_surface.f90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/enforce_acoustic_free_surface.o enforce_acoustic_free_surface.f90
    
### use optimized compilation option for solver only
$O/compute_forces_acoustic.o: compute_forces_acoustic.f90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/compute_forces_acoustic.o compute_forces_acoustic.f90
    
### use optimized compilation option for solver only
$O/compute_forces_elastic.o: compute_forces_elastic.f90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/compute_forces_elastic.o compute_forces_elastic.f90
    
### use optimized compilation option for solver only
$O/compute_gradient_attenuation.o: compute_gradient_attenuation.f90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/compute_gradient_attenuation.o compute_gradient_attenuation.f90
    
$O/compute_energy.o: compute_energy.f90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/compute_energy.o compute_energy.f90
    
$O/compute_vector_field.o: compute_vector_field.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/compute_vector_field.o compute_vector_field.f90
    
$O/compute_pressure.o: compute_pressure.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/compute_pressure.o compute_pressure.f90
    
$O/compute_arrays_source.o: compute_arrays_source.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/compute_arrays_source.o compute_arrays_source.f90
    
$O/create_color_image.o: create_color_image.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/create_color_image.o create_color_image.f90
    
$O/numerical_recipes.o: numerical_recipes.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/numerical_recipes.o numerical_recipes.f90
    
$O/netlib_specfun_erf.o: netlib_specfun_erf.f90
	${F90} $(FLAGS_CHECK) -c -o $O/netlib_specfun_erf.o netlib_specfun_erf.f90
    
$O/define_external_model.o: define_external_model.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/define_external_model.o define_external_model.f90
    
$O/write_seismograms.o: write_seismograms.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/write_seismograms.o write_seismograms.F90
    
$O/part_unstruct.o: part_unstruct.F90 constants_unstruct.h 
	${F90} $(FLAGS_CHECK) -c -o $O/part_unstruct.o part_unstruct.F90

$O/construct_acoustic_surface.o: construct_acoustic_surface.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/construct_acoustic_surface.o construct_acoustic_surface.f90

$O/assemble_MPI.o: assemble_MPI.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/assemble_MPI.o assemble_MPI.F90

$O/attenuation_compute_param.o: attenuation_compute_param.c
	${CC} -c -o $O/attenuation_compute_param.o attenuation_compute_param.c
