#
# Makefile for SPECFEM2D version 5.1
#
# Dimitri Komatitsch, Universite de Pau et des Pays de l'Adour, December 2004
# 
SHELL=/bin/sh

O = obj

# Portland Linux
#F90 = pgf90
#FLAGS_CHECK=-O0 -Mbounds -Mneginfo -Mdclchk
#FLAGS_NOCHECK=-fast -Mnobounds -Minline -Mneginfo -Mdclchk

# Intel Linux
F90 = ifort
FLAGS_CHECK=-O0 -implicitnone -warn stderrors -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -assume byterecl -check bounds
FLAGS_NOCHECK=-O3 -implicitnone -warn stderrors -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -assume byterecl -check nobounds
#FLAGS_NOCHECK = $(FLAGS_CHECK)

#
# g95 (free f95 compiler from http://www.g95.org, still under development, but works)
#
#F90 = g95
#FLAGS_CHECK = -O
#FLAGS_NOCHECK = $(FLAGS_CHECK)

# Dec Alpha
#F90 = f90
#FLAGS_CHECK=-O0 -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nounderflow -check bounds -C
#FLAGS_NOCHECK=-fast -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nounderflow -check nobounds

LINK = $(F90) 

OBJS_MESHFEM2D = $O/meshfem2D.o $O/read_value_parameters.o

OBJS_SPECFEM2D = $O/checkgrid.o $O/datim.o $O/defarrays.o\
        $O/lagrange_poly.o $O/gmat01.o $O/gll_library.o $O/plotgll.o $O/define_derivative_matrices.o\
        $O/plotpost.o $O/locate_receivers.o $O/locate_source_force.o $O/compute_gradient_attenuation.o\
        $O/specfem2D.o $O/write_seismograms.o $O/createnum_fast.o $O/createnum_slow.o\
        $O/define_shape_functions.o $O/cree_image_PNM.o $O/compute_gradient_fluid.o\
        $O/recompute_jacobian.o $O/compute_arrays_source.o $O/locate_source_moment_tensor.o $O/numerical_recipes.o

default: meshfem2D specfem2D convolve_source_timefunction

all: default

clean:
	/bin/rm -r -f xmeshfem2D xmeshfem2D.trace xspecfem2D xspecfem2D.trace $O/*.o *.o $O/*.il *.mod core *.gnu *.ps Ux*.bin Uz*.bin image*.pnm xconvolve_source_timefunction *receiver_line_* plotgnu source.txt *.sem* xcreate_earth_model

meshfem2D: $(OBJS_MESHFEM2D)
	$(LINK) $(FLAGS_CHECK) -o xmeshfem2D $(OBJS_MESHFEM2D)

### use optimized compilation option for solver only
specfem2D: $(OBJS_SPECFEM2D)
	$(LINK) $(FLAGS_NOCHECK) -o xspecfem2D $(OBJS_SPECFEM2D)

convolve_source_timefunction: $O/convolve_source_timefunction.o
	${F90} $(FLAGS_CHECK) -o xconvolve_source_timefunction $O/convolve_source_timefunction.o

create_earth_model: $O/create_earth_model.o
	${F90} $(FLAGS_CHECK) -o xcreate_earth_model $O/create_earth_model.o

$O/checkgrid.o: checkgrid.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/checkgrid.o checkgrid.f90
    
$O/meshfem2D.o: meshfem2D.f90
	${F90} $(FLAGS_CHECK) -c -o $O/meshfem2D.o meshfem2D.f90

$O/createnum_fast.o: createnum_fast.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/createnum_fast.o createnum_fast.f90
    
$O/createnum_slow.o: createnum_slow.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/createnum_slow.o createnum_slow.f90
    
$O/convolve_source_timefunction.o: convolve_source_timefunction.f90
	${F90} $(FLAGS_CHECK) -c -o $O/convolve_source_timefunction.o convolve_source_timefunction.f90

$O/create_earth_model.o: create_earth_model.f90
	${F90} $(FLAGS_CHECK) -c -o $O/create_earth_model.o create_earth_model.f90

$O/read_value_parameters.o: read_value_parameters.f90
	${F90} $(FLAGS_CHECK) -c -o $O/read_value_parameters.o read_value_parameters.f90

$O/datim.o: datim.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/datim.o datim.f90
    
$O/defarrays.o: defarrays.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/defarrays.o defarrays.f90
    
$O/lagrange_poly.o: lagrange_poly.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/lagrange_poly.o lagrange_poly.f90
    
$O/gmat01.o: gmat01.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/gmat01.o gmat01.f90
    
$O/gll_library.o: gll_library.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/gll_library.o gll_library.f90
    
$O/define_derivative_matrices.o: define_derivative_matrices.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/define_derivative_matrices.o define_derivative_matrices.f90
    
$O/plotgll.o: plotgll.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/plotgll.o plotgll.f90
    
$O/plotpost.o: plotpost.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/plotpost.o plotpost.f90
    
$O/locate_receivers.o: locate_receivers.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/locate_receivers.o locate_receivers.f90
    
$O/recompute_jacobian.o: recompute_jacobian.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/recompute_jacobian.o recompute_jacobian.f90
    
$O/locate_source_force.o: locate_source_force.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/locate_source_force.o locate_source_force.f90
    
$O/locate_source_moment_tensor.o: locate_source_moment_tensor.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/locate_source_moment_tensor.o locate_source_moment_tensor.f90
    
$O/define_shape_functions.o: define_shape_functions.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/define_shape_functions.o define_shape_functions.f90
    
### use optimized compilation option for solver only
$O/specfem2D.o: specfem2D.f90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/specfem2D.o specfem2D.f90
    
$O/compute_gradient_attenuation.o: compute_gradient_attenuation.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/compute_gradient_attenuation.o compute_gradient_attenuation.f90
    
$O/compute_gradient_fluid.o: compute_gradient_fluid.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/compute_gradient_fluid.o compute_gradient_fluid.f90
    
$O/compute_arrays_source.o: compute_arrays_source.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/compute_arrays_source.o compute_arrays_source.f90
    
$O/cree_image_PNM.o: cree_image_PNM.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/cree_image_PNM.o cree_image_PNM.f90
    
$O/numerical_recipes.o: numerical_recipes.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/numerical_recipes.o numerical_recipes.f90
    
$O/write_seismograms.o: write_seismograms.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/write_seismograms.o write_seismograms.f90
    
