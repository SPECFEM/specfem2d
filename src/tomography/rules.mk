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
S := ${S_TOP}/src/tomography
$(tomography_OBJECTS): S := ${S_TOP}/src/tomography

#######################################

tomography_TARGETS = \
	$E/xsum_kernels \
	$(EMPTY_MACRO)

tomography_OBJECTS = \
	$(xsum_kernels_OBJECTS) \
	$(EMPTY_MACRO)

# These files come from the shared directory
tomography_SHARED_OBJECTS = \
	$(xsum_kernels_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

tomography_MODULES = \
	$(FC_MODDIR)/tomography_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

####
#### rules for executables
####

.PHONY: all_tomo tomo tomography

all_tomo: $(tomography_TARGETS)

tomo: $(tomography_TARGETS)

tomography: $(tomography_TARGETS)


### single targets
sum_kernels: xsum_kernels
xsum_kernels: $E/xsum_kernels


#######################################

####
#### rules for each program follow
####

#######################################


##
## xsum_kernels
##
xsum_kernels_OBJECTS = \
	$O/tomography_par.tomo_module.o \
	$O/sum_kernels.tomo.o \
	$O/read_model.tomo.o \
	$(EMPTY_MACRO)

xsum_kernels_SHARED_OBJECTS = \
	$O/specfem2D_par.spec_module.o \
	$O/shared_par.shared_module.o \
	$O/exit_mpi.shared.o \
	$O/parallel.shared.o \
	$O/read_parameter_file.mesh.o \
	$O/read_value_parameters.shared.o \
	$O/read_material_table.mesh.o \
	$O/read_interfaces_file.mesh.o \
	$O/read_regions.mesh.o \
	$O/param_reader.cc.o \
	$(EMPTY_MACRO)

${E}/xsum_kernels: $(xsum_kernels_OBJECTS) $(xsum_kernels_SHARED_OBJECTS) $(COND_MPI_OBJECTS)
	@echo ""
	@echo "building xsum_kernels"
	@echo ""
	${FCLINK} -o $@ $+ $(MPILIBS)
	@echo ""



#######################################

###
### Module dependencies
###

####
#### rule for each .o file below
####

$O/%.tomo_module.o: $S/%.f90 ${SETUP}/constants_tomography.h $O/specfem2D_par.spec_module.o $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.tomo.o: $S/%.f90 ${SETUP}/constants_tomography.h $O/tomography_par.tomo_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

