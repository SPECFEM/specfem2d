#
# Makefile for SPECFEM2D version 5.0
#
# Dimitri Komatitsch, Universite de Pau et des Pays de l'Adour, May 2004
# 
SHELL=/bin/sh

O = obj

# Portland Linux
#F90 = pgf90
#FLAGS=-fast -Mnobounds -Minline -Mneginfo -Mdclchk

# Intel Linux
F90 = ifort
#FLAGS=-O0 -implicitnone -warn stderrors -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check bounds -C
FLAGS=-O2 -implicitnone -warn stderrors -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nobounds

#
# g95 (free f95 compiler from http://www.g95.org, still under development, but works)
#
#F90 = g95
#FLAGS = -O

# Dec Alpha
#F90 = f90
#FLAGS=-O0 -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nounderflow -check bounds -C
##FLAGS=-fast -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nounderflow -check nobounds

LINK = $(F90) 
EXEC = xspecfem2D
OBJS = $O/checkgrid.o $O/datim.o $O/defarrays.o\
        $O/lagrange_poly.o $O/gmat01.o $O/gll_library.o $O/plotgll.o $O/define_derivative_matrices.o\
        $O/plotpost.o $O/positrec.o $O/positsource.o $O/q49spec.o $O/compute_gradient_attenuation.o\
        $O/specfem2D.o $O/write_seismograms.o $O/createnum_fast.o $O/createnum_slow.o $O/q49shape.o

default: all

clean:
	/bin/rm -f $(EXEC) $(EXEC).trace $O/*.o *.o $O/*.il *.mod core *.gnu *.ps Ux* Uz* sources;

all: $(OBJS)
	$(LINK) $(FLAGS) -o $(EXEC) $(OBJS)

$O/checkgrid.o: checkgrid.f90 constants.h
	${F90} $(FLAGS) -c -o $O/checkgrid.o checkgrid.f90
    
$O/createnum_fast.o: createnum_fast.f90 constants.h
	${F90} $(FLAGS) -c -o $O/createnum_fast.o createnum_fast.f90
    
$O/createnum_slow.o: createnum_slow.f90 constants.h
	${F90} $(FLAGS) -c -o $O/createnum_slow.o createnum_slow.f90
    
$O/datim.o: datim.f90 constants.h
	${F90} $(FLAGS) -c -o $O/datim.o datim.f90
    
$O/defarrays.o: defarrays.f90 constants.h
	${F90} $(FLAGS) -c -o $O/defarrays.o defarrays.f90
    
$O/lagrange_poly.o: lagrange_poly.f90 constants.h
	${F90} $(FLAGS) -c -o $O/lagrange_poly.o lagrange_poly.f90
    
$O/gmat01.o: gmat01.f90 constants.h
	${F90} $(FLAGS) -c -o $O/gmat01.o gmat01.f90
    
$O/gll_library.o: gll_library.f90 constants.h
	${F90} $(FLAGS) -c -o $O/gll_library.o gll_library.f90
    
$O/define_derivative_matrices.o: define_derivative_matrices.f90 constants.h
	${F90} $(FLAGS) -c -o $O/define_derivative_matrices.o define_derivative_matrices.f90
    
$O/plotgll.o: plotgll.f90 constants.h
	${F90} $(FLAGS) -c -o $O/plotgll.o plotgll.f90
    
$O/plotpost.o: plotpost.f90 constants.h
	${F90} $(FLAGS) -c -o $O/plotpost.o plotpost.f90
    
$O/positrec.o: positrec.f90 constants.h
	${F90} $(FLAGS) -c -o $O/positrec.o positrec.f90
    
$O/positsource.o: positsource.f90 constants.h
	${F90} $(FLAGS) -c -o $O/positsource.o positsource.f90
    
$O/q49shape.o: q49shape.f90 constants.h
	${F90} $(FLAGS) -c -o $O/q49shape.o q49shape.f90
    
$O/q49spec.o: q49spec.f90 constants.h
	${F90} $(FLAGS) -c -o $O/q49spec.o q49spec.f90
    
$O/specfem2D.o: specfem2D.f90 constants.h
	${F90} $(FLAGS) -c -o $O/specfem2D.o specfem2D.f90
    
$O/compute_gradient_attenuation.o: compute_gradient_attenuation.f90 constants.h
	${F90} $(FLAGS) -c -o $O/compute_gradient_attenuation.o compute_gradient_attenuation.f90
    
$O/write_seismograms.o: write_seismograms.f90 constants.h
	${F90} $(FLAGS) -c -o $O/write_seismograms.o write_seismograms.f90
    
