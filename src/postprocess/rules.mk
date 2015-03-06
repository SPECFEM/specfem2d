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
S := ${S_TOP}/src/postprocess
$(postprocess_OBJECTS): S := ${S_TOP}/src/postprocess

#######################################

postprocess_TARGETS = \
	$E/xsum_kernels_ascii \
	$(EMPTY_MACRO)

postprocess_OBJECTS = \
	$(EMPTY_MACRO)

postprocess_MODULES = \
	$(EMPTY_MACRO)

# These files come from the shared directory
postprocess_SHARED_OBJECTS = \
	$(EMPTY_MACRO)

#######################################

##
## sum_kernels_ascii
##
sum_kernels_ascii_postprocess_OBJECTS = \
	$O/sum_kernels_ascii.aux.o \
	$(EMPTY_MACRO)

sum_kernels_ascii_postprocess_SHARED_OBJECTS = \
	$(EMPTY_MACRO)

postprocess_OBJECTS += $(sum_kernels_ascii_postprocess_OBJECTS)
postprocess_SHARED_OBJECTS += $(sum_kernels_ascii_postprocess_SHARED_OBJECTS)


#######################################

####
#### rules for executables
####

aux: $(postprocess_TARGETS)

sum_kernels_ascii: xsum_kernels_ascii
xsum_kernels_ascii: $E/xsum_kernels_ascii

$E/xsum_kernels_ascii: $(sum_kernels_ascii_postprocess_OBJECTS) $(sum_kernels_ascii_postprocess_SHARED_OBJECTS)
	${LINK} $(DEF_FFLAGS) -o $@ $+

#######################################

####
#### rule for each .o file below
####

##
## postprocess
##

$O/%.aux.o: $S/%.f90
	${F90} ${DEF_FFLAGS} -c -o $@ $<

