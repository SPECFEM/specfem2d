#
# Makefile for MESHFEM2D version 5.1
#
# Dimitri Komatitsch, Universite de Pau et des Pays de l'Adour, December 2004
# 
SHELL=/bin/sh

O = obj

# Portland Linux
#F90 = pgf90
#FLAGS=-fast -Mnobounds -Minline -Mneginfo -Mdclchk

# Intel Linux
F90 = ifort
FLAGS=-O0 -implicitnone -warn stderrors -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check bounds
#FLAGS=-O2 -implicitnone -warn stderrors -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nobounds

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
EXEC = xmeshfem2D
OBJS = $O/meshfem2D.o

default: all

clean:
	/bin/rm -f $(EXEC) $(EXEC).trace $O/*.o *.o $O/*.il *.mod core *.gnu *.ps plotgnu

all: $(OBJS)
	$(LINK) $(FLAGS) -o $(EXEC) $(OBJS)

$O/meshfem2D.o: meshfem2D.f90
	${F90} $(FLAGS) -c -o $O/meshfem2D.o meshfem2D.f90
    
