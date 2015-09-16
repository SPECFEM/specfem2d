#========================================================================
#
#                   S P E C F E M 2 D  Version 7 . 0
#                   --------------------------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                        Princeton University, USA
#                and CNRS / University of Marseille, France
#                 (there are currently many more authors!)
# (c) Princeton University and CNRS / University of Marseille, April 2014
#
# This software is a computer program whose purpose is to solve
# the two-dimensional viscoelastic anisotropic or poroelastic wave equation
# using a spectral-element method (SEM).
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and Inria at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The full text of the license is available in file "LICENSE".
#
#========================================================================

## compilation directories
S := ${S_TOP}/src/tomography/postprocess_sensitivity_kernels
$(postprocess_OBJECTS): S := ${S_TOP}/src/tomography/postprocess_sensitivity_kernels

#######################################

postprocess_TARGETS = \
	$E/xcombine_sem \
	$E/xsum_kernels_ascii \
	$E/xsmooth_sem \
	$(EMPTY_MACRO)

postprocess_OBJECTS = \
	$(xcombine_sem_OBJECTS) \
	$(xsum_kernels_ascii_OBJECTS) \
	$(xsmooth_sem_OBJECTS) \
	$(EMPTY_MACRO)

postprocess_MODULES = \
	$(FC_MODDIR)/postprocess_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

postprocess_SHARED_OBJECTS = \
	$(EMPTY_MACRO)


#######################################

####
#### rules for executables
####

postprocess: $(postprocess_TARGETS)

combine_sem: xcombine_sem
xcombine_sem: $E/xcombine_sem

sum_kernels_ascii: xsum_kernels_ascii
xsum_kernels_ascii: $E/xsum_kernels_ascii

smooth_sem: xsmooth_sem
xsmooth_sem: $E/xsmooth_sem

#######################################

##
## combine_sem
##

xcombine_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/combine_sem.postprocess.o \
	$O/parse_kernel_names.postprocess.o \
	$(EMPTY_MACRO)

xcombine_sem_SHARED_OBJECTS = \
	$(EMPTY_MACRO)

${E}/xcombine_sem: $(xcombine_sem_OBJECTS) $(xcombine_sem_SHARED_OBJECTS)
	${LINK} $(DEF_FFLAGS) -o $@ $+


#######################################

##
## sum_kernels_ascii
##

xsum_kernels_ascii_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/sum_kernels_ascii.postprocess.o \
	$(EMPTY_MACRO)

xsum_kernels_ascii_SHARED_OBJECTS = \
	$(EMPTY_MACRO)

${E}/xsum_kernels_ascii: $(xsum_kernels_ascii_OBJECTS) $(xsum_kernels_ascii_SHARED_OBJECTS)
	${LINK} $(DEF_FFLAGS) -o $@ $+

#######################################

##
## smooth_sum
##

xsmooth_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/smooth_sem.postprocess.o \
	$O/parse_kernel_names.postprocess.o \
	$(EMPTY_MACRO)

${E}/xsmooth_sem: $(xsmooth_sem_OBJECTS) 
	${LINK} $(DEF_FFLAGS) -o $@ $+

#######################################

####
#### rule for each .o file below
####

##
## postprocess
##

$O/%.postprocess_module.o: $S/%.f90 ${SETUP}/constants.h $O/specfem2D_par.spec.o
	${F90} ${DEF_FFLAGS} -c -o $@ $<

$O/%.postprocess.o: $S/%.f90
	${F90} ${DEF_FFLAGS} -c -o $@ $<

