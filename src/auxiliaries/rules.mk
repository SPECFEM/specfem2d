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
S := ${S_TOP}/src/auxiliaries
$(auxiliaries_OBJECTS): S := ${S_TOP}/src/auxiliaries

#######################################

auxiliaries_TARGETS = \
	$E/xadj_seismogram \
	$E/xcheck_quality_external_mesh \
	$E/xconvolve_source_timefunction \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS = \
	$(EMPTY_MACRO)

auxiliaries_MODULES = \
	$(EMPTY_MACRO)

# These files come from the shared directory
auxiliaries_SHARED_OBJECTS = \
	$(EMPTY_MACRO)

#######################################

##
## adj_seismogram
##
adj_seismogram_auxiliaries_OBJECTS = \
	$O/adj_seismogram.aux.o \
	$(EMPTY_MACRO)

adj_seismogram_auxiliaries_SHARED_OBJECTS = \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS += $(adj_seismogram_auxiliaries_OBJECTS)
auxiliaries_SHARED_OBJECTS += $(adj_seismogram_auxiliaries_SHARED_OBJECTS)

##
## check_quality_external_mesh
##
check_quality_external_mesh_auxiliaries_OBJECTS = \
	$O/check_quality_external_mesh.aux.o \
	$(EMPTY_MACRO)

check_quality_external_mesh_auxiliaries_SHARED_OBJECTS = \
	$O/define_shape_functions.shared.o \
	$O/param_reader.cc.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS += $(check_quality_external_mesh_auxiliaries_OBJECTS)
auxiliaries_SHARED_OBJECTS += $(check_quality_external_mesh_auxiliaries_SHARED_OBJECTS)

##
## convolve_source_timefunction
##
convolve_source_timefunction_auxiliaries_OBJECTS = \
	$O/convolve_source_timefunction.aux.o \
	$(EMPTY_MACRO)

convolve_source_timefunction_auxiliaries_SHARED_OBJECTS = \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS += $(convolve_source_timefunction_auxiliaries_OBJECTS)
auxiliaries_SHARED_OBJECTS += $(convolve_source_timefunction_auxiliaries_SHARED_OBJECTS)

#######################################

####
#### rules for executables
####

aux: $(auxiliaries_TARGETS)

adj_seismogram: xadj_seismogram
xadj_seismogram: $E/xadj_seismogram

check_quality_external_mesh: xcheck_quality_external_mesh
xcheck_quality_external_mesh: $E/xcheck_quality_external_mesh

convolve_source_timefunction: xconvolve_source_timefunction
xconvolve_source_timefunction: $E/xconvolve_source_timefunction


$E/xadj_seismogram: $(adj_seismogram_auxiliaries_OBJECTS) $(adj_seismogram_auxiliaries_SHARED_OBJECTS)
	${FCLINK} -o $@ $+

$E/xcheck_quality_external_mesh: $(check_quality_external_mesh_auxiliaries_OBJECTS) $(check_quality_external_mesh_auxiliaries_SHARED_OBJECTS)
	${FCLINK} -o $@ $+

$E/xconvolve_source_timefunction: $(convolve_source_timefunction_auxiliaries_OBJECTS) $(convolve_source_timefunction_auxiliaries_SHARED_OBJECTS)
	${FCLINK} -o $@ $+

#######################################

####
#### rule for each .o file below
####

##
## auxiliaries
##

$O/%.aux.o: $S/%.f90
	${F90} ${FCFLAGS_f90} -c -o $@ $<

