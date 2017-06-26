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
S := ${S_TOP}/src/shared
$(shared_OBJECTS): S = ${S_TOP}/src/shared

#######################################

shared_TARGETS = \
	$(shared_OBJECTS) \
	$(EMPTY_MACRO)


shared_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/define_shape_functions.shared.o \
	$O/exit_mpi.shared.o \
	$O/force_ftz.cc.o \
	$O/gll_library.shared.o \
	$O/lagrange_poly.shared.o \
	$O/parallel.shared.o \
	$O/param_reader.cc.o \
	$O/read_value_parameters.shared.o \
	$O/set_color_palette.shared.o \
	$(EMPTY_MACRO)


shared_MODULES = \
	$(FC_MODDIR)/constants.$(FC_MODEXT) \
	$(FC_MODDIR)/shared_input_parameters.$(FC_MODEXT) \
	$(FC_MODDIR)/shared_parameters.$(FC_MODEXT) \
	$(FC_MODDIR)/my_mpi_communicator.$(FC_MODEXT) \
	$(EMPTY_MACRO)


#######################################

####
#### rule for each .o file below
####


##
## shared
##

$O/%.shared_module.o: $S/%.f90 ${SETUP}/constants.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared_module.o: $S/%.F90 ${SETUP}/constants.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared.o: $S/%.f90 $O/shared_par.shared_module.o
	${F90} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared.o: $S/%.F90 $O/shared_par.shared_module.o
	${F90} ${FCFLAGS_f90} -c -o $@ $<


##
## C compilation
##

$O/%.cc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) -o $@ $<
