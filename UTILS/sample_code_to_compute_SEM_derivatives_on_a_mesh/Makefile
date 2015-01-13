#
# Makefile
#
# Dimitri Komatitsch, CNRS, Marseille, France, August 2014
# 
SHELL=/bin/sh

B = .
O = obj
SRC = .

F90 = gfortran
FLAGS=-std=f2003 -fimplicit-none -frange-check -O2 -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow

#F90 = ifort
#FLAGS=-check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -vec-report0 -std03 -implicitnone -gen-interfaces -warn all

LINK = $(F90) 
EXEC = xcompute_gradient_Earth
EXEC2 = xcompute_gradient_rectangle
OBJS = $O/define_derivation_matrices.o $O/compute_gradient_of_a_field_with_the_SEM_Earth.o $O/define_shape_functions.o\
	$O/gll_library.o $O/lagrange_poly.o $O/recompute_jacobian.o
OBJS2 = $O/define_derivation_matrices.o $O/compute_gradient_of_a_field_with_the_SEM_rectangle.o $O/define_shape_functions.o\
	$O/gll_library.o $O/lagrange_poly.o $O/recompute_jacobian.o
DIRS = .

.f90.o:
	$(F90) $(FLAGS) $*.f90

default : $(OBJS) $(OBJS2) make_the_code make_the_code2 cleanmods

make_the_code : $(DIRS) $(OBJS)
	$(LINK) $(OBJS) -o $(EXEC)

make_the_code2 : $(DIRS) $(OBJS2)
	$(LINK) $(OBJS2) -o $(EXEC2)

all : $(OBJS) $(OBJS2) make_the_code make_the_code2 cleanmods

clean :
	/bin/rm -f $(EXEC) $(EXEC).trace $(EXEC2) $(EXEC2).trace $O/*.o *.o *.mod core *__genmod.*;

cleanmods : make_the_code make_the_code2
	/bin/rm -f *.mod core *__genmod.*;

$O/compute_gradient_of_a_field_with_the_SEM_rectangle.o: $(SRC)/compute_gradient_of_a_field_with_the_SEM_rectangle.f90
	${F90} $(FLAGS) -c -o $O/compute_gradient_of_a_field_with_the_SEM_rectangle.o $(SRC)/compute_gradient_of_a_field_with_the_SEM_rectangle.f90
    
$O/compute_gradient_of_a_field_with_the_SEM_Earth.o: $(SRC)/compute_gradient_of_a_field_with_the_SEM_Earth.f90
	${F90} $(FLAGS) -c -o $O/compute_gradient_of_a_field_with_the_SEM_Earth.o $(SRC)/compute_gradient_of_a_field_with_the_SEM_Earth.f90
    
$O/define_derivation_matrices.o: $(SRC)/define_derivation_matrices.f90
	${F90} $(FLAGS) -c -o $O/define_derivation_matrices.o $(SRC)/define_derivation_matrices.f90
    
$O/define_shape_functions.o: $(SRC)/define_shape_functions.f90
	${F90} $(FLAGS) -c -o $O/define_shape_functions.o $(SRC)/define_shape_functions.f90
    
$O/gll_library.o: $(SRC)/gll_library.f90
	${F90} $(FLAGS) -c -o $O/gll_library.o $(SRC)/gll_library.f90
    
$O/lagrange_poly.o: $(SRC)/lagrange_poly.f90
	${F90} $(FLAGS) -c -o $O/lagrange_poly.o $(SRC)/lagrange_poly.f90
    
$O/recompute_jacobian.o: $(SRC)/recompute_jacobian.f90
	${F90} $(FLAGS) -c -o $O/recompute_jacobian.o $(SRC)/recompute_jacobian.f90
    
