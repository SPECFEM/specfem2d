#
# Makefile simple pour f90
#
# Dimitri Komatitsch, Harvard University, May 1998
# 
SHELL=/bin/sh

B = .
O = obj
SRC = .
TRNDIR = ./bak

# Portland
#F90 = pgf90
#FLAGS=-c -fast -Mnobounds -Minline -Mneginfo -Mdclchk

# Dec Alpha
#F90 = f90
#FLAGS=-c -fast -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nounderflow -check bounds

# Intel Linux
F90 = ifort
FLAGS=-c -O3 -e95 -implicitnone

LINK = $(F90) 
EXEC = xspecfem
OBJS = $O/modules.o $O/calc_energie.o $O/calcdivcurl.o $O/calcforce.o $O/checkgrid.o\
	$O/checksource.o $O/contol.o $O/datim.o\
	$O/defarrays.o $O/dirac.o $O/dircty.o $O/endw1.o $O/endw2.o $O/gammaf.o\
	$O/getelspec.o $O/getltf.o $O/getrecepts.o $O/getspec.o $O/gmat01.o $O/hdgll.o $O/hgll.o\
	$O/intseq.o $O/jacg.o $O/jacobf.o $O/modifperio.o $O/plotavs.o $O/plotgll.o\
	$O/plotpost.o $O/plotvect.o $O/pndleg.o $O/pnleg.o\
	$O/pnormj.o $O/positrec.o $O/positsource.o $O/q49spec.o $O/qinpspec.o $O/qmasspec.o\
	$O/qsumspec.o $O/ricker.o $O/setcor.o $O/specfem.o $O/storearray.o \
	$O/writeseis.o $O/zwgjd.o $O/zwgljd.o $O/createnum_fast.o $O/createnum_slow.o\
	$O/q49shape.o
DIRS = .

.f90.o:
	$(F90) $(FLAGS) $*.f90

default : $(OBJS) mxspecfem

TIMESTAMP: $(SRC)/*.f90 $(SRC)/Makefile
	cp $? ${TRNDIR}
	touch TIMESTAMP

mxspecfem : $(DIRS)
	$(LINK) $(OBJS) -o $(EXEC)

all : clean $(OBJS) mxspecfem

clean :
	/bin/rm -f $(EXEC) $(EXEC).trace $O/*.o *.o *.mod core *.gnu *.ps Ux* Uz* sources;

$O/calc_energie.o: $(SRC)/calc_energie.f90
	${F90} $(FLAGS) -c -o $O/calc_energie.o $(SRC)/calc_energie.f90
    
$O/calcdivcurl.o: $(SRC)/calcdivcurl.f90
	${F90} $(FLAGS) -c -o $O/calcdivcurl.o $(SRC)/calcdivcurl.f90
    
$O/calcforce.o: $(SRC)/calcforce.f90
	${F90} $(FLAGS) -c -o $O/calcforce.o $(SRC)/calcforce.f90
    
$O/checkgrid.o: $(SRC)/checkgrid.f90
	${F90} $(FLAGS) -c -o $O/checkgrid.o $(SRC)/checkgrid.f90
    
$O/checksource.o: $(SRC)/checksource.f90
	${F90} $(FLAGS) -c -o $O/checksource.o $(SRC)/checksource.f90
    
$O/contol.o: $(SRC)/contol.f90
	${F90} $(FLAGS) -c -o $O/contol.o $(SRC)/contol.f90
    
$O/createnum_fast.o: $(SRC)/createnum_fast.f90
	${F90} $(FLAGS) -c -o $O/createnum_fast.o $(SRC)/createnum_fast.f90
    
$O/createnum_slow.o: $(SRC)/createnum_slow.f90
	${F90} $(FLAGS) -c -o $O/createnum_slow.o $(SRC)/createnum_slow.f90
    
$O/datim.o: $(SRC)/datim.f90
	${F90} $(FLAGS) -c -o $O/datim.o $(SRC)/datim.f90
    
$O/defarrays.o: $(SRC)/defarrays.f90
	${F90} $(FLAGS) -c -o $O/defarrays.o $(SRC)/defarrays.f90
    
$O/dirac.o: $(SRC)/dirac.f90
	${F90} $(FLAGS) -c -o $O/dirac.o $(SRC)/dirac.f90
    
$O/dircty.o: $(SRC)/dircty.f90
	${F90} $(FLAGS) -c -o $O/dircty.o $(SRC)/dircty.f90
    
$O/endw1.o: $(SRC)/endw1.f90
	${F90} $(FLAGS) -c -o $O/endw1.o $(SRC)/endw1.f90
    
$O/endw2.o: $(SRC)/endw2.f90
	${F90} $(FLAGS) -c -o $O/endw2.o $(SRC)/endw2.f90
    
$O/ezfftf.o: $(SRC)/ezfftf.f90
	${F90} $(FLAGS) -c -o $O/ezfftf.o $(SRC)/ezfftf.f90
    
$O/gammaf.o: $(SRC)/gammaf.f90
	${F90} $(FLAGS) -c -o $O/gammaf.o $(SRC)/gammaf.f90
    
$O/getelspec.o: $(SRC)/getelspec.f90
	${F90} $(FLAGS) -c -o $O/getelspec.o $(SRC)/getelspec.f90
    
$O/getltf.o: $(SRC)/getltf.f90
	${F90} $(FLAGS) -c -o $O/getltf.o $(SRC)/getltf.f90
    
$O/getrecepts.o: $(SRC)/getrecepts.f90
	${F90} $(FLAGS) -c -o $O/getrecepts.o $(SRC)/getrecepts.f90
    
$O/getspec.o: $(SRC)/getspec.f90
	${F90} $(FLAGS) -c -o $O/getspec.o $(SRC)/getspec.f90
    
$O/gmat01.o: $(SRC)/gmat01.f90
	${F90} $(FLAGS) -c -o $O/gmat01.o $(SRC)/gmat01.f90
    
$O/hdgll.o: $(SRC)/hdgll.f90
	${F90} $(FLAGS) -c -o $O/hdgll.o $(SRC)/hdgll.f90
    
$O/hgll.o: $(SRC)/hgll.f90
	${F90} $(FLAGS) -c -o $O/hgll.o $(SRC)/hgll.f90
    
$O/intseq.o: $(SRC)/intseq.f90
	${F90} $(FLAGS) -c -o $O/intseq.o $(SRC)/intseq.f90
    
$O/jacg.o: $(SRC)/jacg.f90
	${F90} $(FLAGS) -c -o $O/jacg.o $(SRC)/jacg.f90
    
$O/jacobf.o: $(SRC)/jacobf.f90
	${F90} $(FLAGS) -c -o $O/jacobf.o $(SRC)/jacobf.f90
    
$O/modifperio.o: $(SRC)/modifperio.f90
	${F90} $(FLAGS) -c -o $O/modifperio.o $(SRC)/modifperio.f90
    
$O/modules.o: $(SRC)/modules.f90
	${F90} $(FLAGS) -c -o $O/modules.o $(SRC)/modules.f90
    
$O/plotavs.o: $(SRC)/plotavs.f90
	${F90} $(FLAGS) -c -o $O/plotavs.o $(SRC)/plotavs.f90
    
$O/plotgll.o: $(SRC)/plotgll.f90
	${F90} $(FLAGS) -c -o $O/plotgll.o $(SRC)/plotgll.f90
    
$O/plotpost.o: $(SRC)/plotpost.f90
	${F90} $(FLAGS) -c -o $O/plotpost.o $(SRC)/plotpost.f90
    
$O/plotvect.o: $(SRC)/plotvect.f90
	${F90} $(FLAGS) -c -o $O/plotvect.o $(SRC)/plotvect.f90
    
$O/pndleg.o: $(SRC)/pndleg.f90
	${F90} $(FLAGS) -c -o $O/pndleg.o $(SRC)/pndleg.f90
    
$O/pnleg.o: $(SRC)/pnleg.f90
	${F90} $(FLAGS) -c -o $O/pnleg.o $(SRC)/pnleg.f90
    
$O/pnormj.o: $(SRC)/pnormj.f90
	${F90} $(FLAGS) -c -o $O/pnormj.o $(SRC)/pnormj.f90
    
$O/positrec.o: $(SRC)/positrec.f90
	${F90} $(FLAGS) -c -o $O/positrec.o $(SRC)/positrec.f90
    
$O/positsource.o: $(SRC)/positsource.f90
	${F90} $(FLAGS) -c -o $O/positsource.o $(SRC)/positsource.f90
    
$O/q49shape.o: $(SRC)/q49shape.f90
	${F90} $(FLAGS) -c -o $O/q49shape.o $(SRC)/q49shape.f90
    
$O/q49spec.o: $(SRC)/q49spec.f90
	${F90} $(FLAGS) -c -o $O/q49spec.o $(SRC)/q49spec.f90
    
$O/qinpspec.o: $(SRC)/qinpspec.f90
	${F90} $(FLAGS) -c -o $O/qinpspec.o $(SRC)/qinpspec.f90
    
$O/qmasspec.o: $(SRC)/qmasspec.f90
	${F90} $(FLAGS) -c -o $O/qmasspec.o $(SRC)/qmasspec.f90
    
$O/qsumspec.o: $(SRC)/qsumspec.f90
	${F90} $(FLAGS) -c -o $O/qsumspec.o $(SRC)/qsumspec.f90
    
$O/ricker.o: $(SRC)/ricker.f90
	${F90} $(FLAGS) -c -o $O/ricker.o $(SRC)/ricker.f90
    
$O/setcor.o: $(SRC)/setcor.f90
	${F90} $(FLAGS) -c -o $O/setcor.o $(SRC)/setcor.f90
    
$O/specfem.o: $(SRC)/specfem.f90
	${F90} $(FLAGS) -c -o $O/specfem.o $(SRC)/specfem.f90
    
$O/storearray.o: $(SRC)/storearray.f90
	${F90} $(FLAGS) -c -o $O/storearray.o $(SRC)/storearray.f90
    
$O/writeseis.o: $(SRC)/writeseis.f90
	${F90} $(FLAGS) -c -o $O/writeseis.o $(SRC)/writeseis.f90
    
$O/zwgjd.o: $(SRC)/zwgjd.f90
	${F90} $(FLAGS) -c -o $O/zwgjd.o $(SRC)/zwgjd.f90
    
$O/zwgljd.o: $(SRC)/zwgljd.f90
	${F90} $(FLAGS) -c -o $O/zwgljd.o $(SRC)/zwgljd.f90
    
