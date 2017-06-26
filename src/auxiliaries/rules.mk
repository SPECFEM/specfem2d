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
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
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
	$O/shared_par.shared_module.o \
	$O/define_shape_functions.shared.o \
	$O/exit_mpi.shared.o \
	$O/parallel.shared.o \
	$O/read_parameter_file.mesh.o \
	$O/read_value_parameters.shared.o \
	$O/read_material_table.mesh.o \
	$O/read_interfaces_file.mesh.o \
	$O/read_regions.mesh.o \
	$O/param_reader.cc.o \
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
	@echo ""
	@echo "building xadj_seismogram"
	@echo ""
	${FCLINK} -o $@ $+
	@echo ""

$E/xcheck_quality_external_mesh: $(check_quality_external_mesh_auxiliaries_OBJECTS) $(check_quality_external_mesh_auxiliaries_SHARED_OBJECTS)
	@echo ""
	@echo "building xcheck_quality_external_mesh"
	@echo ""
	${FCLINK} -o $@ $+
	@echo ""

$E/xconvolve_source_timefunction: $(convolve_source_timefunction_auxiliaries_OBJECTS) $(convolve_source_timefunction_auxiliaries_SHARED_OBJECTS)
	@echo ""
	@echo "building xconvolve_source_timefunction"
	@echo ""
	${FCLINK} -o $@ $+
	@echo ""

#######################################

####
#### rule for each .o file below
####

##
## auxiliaries
##

$O/%.aux.o: $S/%.f90
	${F90} ${FCFLAGS_f90} -c -o $@ $<

