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
S := ${S_TOP}/src/meshfem2D
$(meshfem2D_OBJECTS): S := ${S_TOP}/src/meshfem2D

#######################################

####
#### targets
####

meshfem2D_TARGETS = \
	$E/xmeshfem2D \
	$(EMPTY_MACRO)

meshfem2D_OBJECTS = \
	$O/get_node_number.mesh.o \
	$O/part_unstruct.mesh.o \
	$O/read_interfaces_file.mesh.o \
	$O/read_materials.mesh.o \
	$O/read_parameter_file.mesh.o \
	$O/read_regions.mesh.o \
	$O/read_source_file.mesh.o \
	$O/save_databases.mesh.o \
	$O/save_gnuplot_file.mesh.o \
	$O/save_stations_file.mesh.o \
	$O/spline_routines.mesh.o \
	$O/meshfem2D.mesh.o \
	$(EMPTY_MACRO)

meshfem2D_MODULES = \
	$(FC_MODDIR)/interfaces_file.$(FC_MODEXT) \
	$(FC_MODDIR)/parameter_file.$(FC_MODEXT) \
	$(FC_MODDIR)/part_unstruct.$(FC_MODEXT) \
	$(FC_MODDIR)/source_file.$(FC_MODEXT) \
	$(EMPTY_MACRO)

meshfem2D_SHARED_OBJECTS = \
	$O/read_value_parameters.shared.o \
	$O/param_reader.cc.o \
	$(EMPTY_MACRO)

$(SCOTCH_INCDIR)/scotchf.h: scotch_library
scotch_library:
ifeq ($(USE_BUNDLED_SCOTCH),1)
	@echo "Using bundled Scotch"
	$(MAKE) -C "$(SCOTCH_DIR)/src"
else
	@echo "Not using bundled Scotch"
endif


#######################################

####
#### rules for executables
####

mesh: $(meshfem2D_TARGETS)

meshfem2D: xmeshfem2D
xmeshfem2D: $E/xmeshfem2D

$E/xmeshfem2D: $(meshfem2D_OBJECTS) $(meshfem2D_SHARED_OBJECTS)
	$(FCLINK) -o ${E}/xmeshfem2D $(meshfem2D_OBJECTS) $(meshfem2D_SHARED_OBJECTS) $(MPILIBS)


#######################################

###
### Module dependencies
###

$O/meshfem2D.mesh.o: $O/part_unstruct.mesh.o $O/read_interfaces_file.mesh.o $O/read_parameter_file.mesh.o $O/read_source_file.mesh.o

ifdef SCOTCH_INCDIR
$O/part_unstruct.mesh.o: $(SCOTCH_INCDIR)/scotchf.h
endif

$O/save_databases.mesh.o: $O/part_unstruct.mesh.o $O/read_parameter_file.mesh.o $O/read_source_file.mesh.o

####
#### rule to build each .o file below
####

$O/%.mesh.o: $S/%.f90 ${SETUP}/constants.h
	${F90} ${FCFLAGS_f90} -c -o $@ $<

$O/%.mesh.o: $S/%.F90 ${SETUP}/constants.h
	${F90} ${FCFLAGS_f90} -c -o $@ $<
