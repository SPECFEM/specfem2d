%.o: %.f90
	$(F90) $(F90FLAGS1) -c $< -o $@ 
%.o: %.F90
	$(F90) $(F90FLAGS1) -c $< -o $@
%.out: %.f90
	$(F90) $(F90FLAGS2) -o $@ $^  $(LIBS)
%.out: %.F90
	$(F90) $(F90FLAGS2) -o $@ $^  $(LIBS)
#LIBS=/usr/lib64/liblapack.so.3.2.1 /usr/lib64/libblas.so.3.2.1
LIBS=/home/komatits/lapack-3.7.0/liblapack.a /home/komatits/lapack-3.7.0/librefblas.a

F90 = gfortran
F90FLAGS1 =  -O3 -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8  -Iincludefic/ -Iincludemod/
F90FLAGS2 =  -O3 -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8  -Iincludefic/ -Iincludemod/ 
RM = \rm -f
ARFLAGS = rvU

#############
# VARIABLES #
#############
# 1. LIBRAIRIES INDEPENDANTES
# general
CHEGEN = lib/libgeneral
LSTGEN = $(foreach DIR,$(CHEGEN),$(wildcard $(DIR)/*.F90))
OBJGEN = $(LSTGEN:.F90=.o)
LIBGEN = lib/libgeneral/libgeneral.a
# acousacous
CHEAA = lib/libacousacous
LSTAA = $(foreach DIR,$(CHEAA),$(wildcard $(DIR)/*.F90))
OBJAA = $(LSTAA:.F90=.o)
LIBAA = lib/libacousacous/libacousacous.a
# acouselasto
CHEAE = lib/libacouselasto
LSTAE = $(foreach DIR,$(CHEAE),$(wildcard $(DIR)/*.F90))
OBJAE = $(LSTAE:.F90=.o)
LIBAE = lib/libacouselasto/libacouselasto.a
# acousporo
CHEAP = lib/libacousporo
LSTAP = $(foreach DIR,$(CHEAP),$(wildcard $(DIR)/*.F90))
OBJAP = $(LSTAP:.F90=.o)
LIBAP = lib/libacousporo/libacousporo.a
# poroporo
CHEPP = lib/libporoporo
LSTPP = $(foreach DIR,$(CHEPP),$(wildcard $(DIR)/*.F90))
OBJPP = $(LSTPP:.F90=.o)
LIBPP = lib/libporoporo/libporoporo.a
# 2. LIBRAIRIE DES MODULES DE DECLARATION DES VAR GLOBALES
# Variables necessaires a la creation de la librairie des mod
CHEMOD = mod
LSTMOD = $(foreach DIR,$(CHEMOD),$(wildcard $(DIR)/*.F90)) 
OBJMOD = $(LSTMOD:.F90=.o)
LIBMOD = mod/libmod.a
# 3. LIBRAIRIES DEPENDANT DES MODULES ET DES LIBRAIRIES 1.
# Variables necessaires a la creation de l'executable
CHEBIN = lib/bin
LSTBIN = $(foreach DIR,$(CHEBIN),$(wildcard $(DIR)/*.[Ff]90))
BINBIN = $(addsuffix .out,$(basename $(LSTBIN)))
BIN = $(subst $(CHEBIN),lib,$(BINBIN))

###########
# CHEMINS #
###########
vpath %.F90 $(CHEBIN)
vpath %.m includemod/
vpath %.a lib/libgeneral lib/libacousacous  lib/libacouselasto lib/libacousporo lib/libporoporo mod/   

###############
# COMPILATION #
###############
# 1. LIBRAIRIES INDEPENDANTES
srcgeneral: $(LIBGEN)
$(LIBGEN): $(LIBGEN)($(OBJGEN))
srcacousacous: $(LIBAA)
$(LIBAA): $(LIBAA)($(OBJAA))
srcacouselasto: $(LIBAE)
$(LIBAE): $(LIBAE)($(OBJAE))
srcacousporo: $(LIBAP)
$(LIBAP): $(LIBAP)($(OBJAP))
srcporoporo: $(LIBPP)
$(LIBPP): $(LIBPP)($(OBJPP))
# Creation des sources 1
src1:  srcgeneral srcacousacous  srcacouselasto srcacousporo srcporoporo
# 2. LIBRAIRIE DES MODULES DE DECLARATION DES VAR GLOBALES
# Creation de la librairie des mod
srcmod: $(LIBMOD)
$(LIBMOD): $(LIBMOD)($(OBJMOD))
	mv *.mod includemod/.
# 3. LIBRAIRIES DEPENDANT DES MODULES ET DES LIBRAIRIES 1.
# 5. PROGRAMME PRINCIPAL
# Creation de l'executable du programme principal
all:  $(LIBMOD) $(LIBGEN) $(LIBAA) $(LIBAE) $(LIBAP) $(LIBPP)   $(BIN)
$(BIN): $(BINBIN)
	mv $(CHEBIN)/$(@F) .
$(BINBIN): -lmod  -lacousacous -lacouselasto -lacousporo -lporoporo -lgeneral



###############
# SUPPRESSION #
###############
# 1. LIBRAIRIES INDEPENDANTES
# Supression de la librairie des Tab
clean-srcgen:	
	$(RM) lib/libgeneral/*.a
clean-srcaa:	
	$(RM) lib/libacousacous/*.a
clean-srcae:	
	$(RM) lib/libacouselasto/*.a
clean-srcap:	
	$(RM) lib/libacousporo/*.a
clean-srcpp:	
	$(RM) lib/libporoporo/*.a
# Suppression des librairies 1
clean-src1: clean-srcgen  clean-srcaa clean-srcae clean-srcap clean-srcpp
# 2. LIBRAIRIE DES MODULES DE DECLARATION DES VAR GLOBALES
# Suppression de la librairie des mod
clean-srcmod:	
	$(RM) mod/*.a
	$(RM) mod/*.o
	$(RM) includemod/*.mod
# 3. LIBRAIRIES DEPENDANT DES MODULES ET DES LIBRAIRIES 1.
# Suppression de l'executable
clean-out:
	$(RM) *.out

# Suppression de tout
clean: clean-src1 clean-srcmod  clean-out
	$(RM) lib/*/*.o





