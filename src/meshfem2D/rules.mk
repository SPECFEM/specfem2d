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
	$O/meshfem2D_par.mesh_module.o \
	$O/compute_elements_load_par.mesh.o \
	$O/decompose_mesh.mesh.o \
	$O/determine_abs_surface.mesh.o \
	$O/determine_acoustic_surface.mesh.o \
	$O/get_node_number.mesh.o \
	$O/metis_partitioning.mesh.o \
	$O/part_unstruct.mesh.o \
	$O/read_external_mesh_files.mesh.o \
	$O/read_interfaces_file.mesh.o \
	$O/read_material_table.mesh.o \
	$O/read_parameter_file.mesh.o \
	$O/read_regions.mesh.o \
	$O/read_source_file.mesh.o \
	$O/read_mesh_files.mesh.o \
	$O/repartition_coupling.mesh.o \
	$O/rotate_mesh.mesh.o \
	$O/save_databases.mesh.o \
	$O/save_gnuplot_file.mesh.o \
	$O/save_stations_file.mesh.o \
	$O/scotch_partitioning.mesh.o \
	$O/spline_routines.mesh.o \
	$O/meshfem2D.mesh.o \
	$(EMPTY_MACRO)

meshfem2D_MODULES = \
	$(FC_MODDIR)/decompose_par.$(FC_MODEXT) \
	$(FC_MODDIR)/part_unstruct_par.$(FC_MODEXT) \
	$(FC_MODDIR)/compute_elements_load_par.$(FC_MODEXT) \
	$(FC_MODDIR)/source_file_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

meshfem2D_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/read_value_parameters.shared.o \
	$O/exit_mpi.shared.o \
	$O/parallel.shared.o \
	$O/param_reader.cc.o \
	$(EMPTY_MACRO)

# default mesher flags
FCFLAGS_f90_MESH = $(FCFLAGS_f90)
MPILIBS_MESH = $(MPILIBS)

# only mesher needs scotch for compilation of parallel version
ifeq ($(SCOTCH),yes)
  FCFLAGS_f90_MESH += $(SCOTCH_FLAGS)
  ## scotch libraries
  MPILIBS_MESH = $(MPILIBS) $(SCOTCH_LIBS)
endif

#######################################

####
#### rules for executables
####

mesh: $(meshfem2D_TARGETS)

meshfem2D: xmeshfem2D
xmeshfem2D: $E/xmeshfem2D

$E/xmeshfem2D: $(meshfem2D_OBJECTS) $(meshfem2D_SHARED_OBJECTS)
	@echo ""
	@echo "building xmeshfem2D"
	@echo ""
	$(FCLINK) -o ${E}/xmeshfem2D $(meshfem2D_OBJECTS) $(meshfem2D_SHARED_OBJECTS) $(MPILIBS_MESH)
	@echo ""

# target for SCOTCH
$(SCOTCH_INCDIR)/scotchf.h: scotch_library
scotch_library:
ifeq ($(USE_BUNDLED_SCOTCH),1)
	@echo "Using bundled Scotch"
	$(MAKE) -C "$(SCOTCH_DIR)/src"
else
	@echo "Not using bundled Scotch"
endif


#######################################

###
### Module dependencies
###

$O/decompose_mesh.mesh.o: $O/compute_elements_load_par.mesh.o
$O/meshfem2D.mesh.o: $O/compute_elements_load_par.mesh.o
$O/metis_partitioning.mesh.o: $O/compute_elements_load_par.mesh.o
$O/read_external_mesh_files.mesh.o: $O/compute_elements_load_par.mesh.o
$O/scotch_partitioning.mesh.o: $O/compute_elements_load_par.mesh.o

ifdef SCOTCH_INCDIR
$O/scotch_partitioning.mesh.o: $(SCOTCH_INCDIR)/scotchf.h
endif


# Version file
$O/meshfem2D.mesh.o: ${SETUP}/version.fh


####
#### rule to build each .o file below
####

$O/%.mesh_module.o: $S/%.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o
	${F90} ${FCFLAGS_f90_MESH} -c -o $@ $<

$O/%.mesh.o: $S/%.f90 ${SETUP}/constants.h $O/meshfem2D_par.mesh_module.o
	${F90} ${FCFLAGS_f90_MESH} -c -o $@ $<

$O/%.mesh.o: $S/%.F90 ${SETUP}/constants.h $O/meshfem2D_par.mesh_module.o
	${F90} ${FCFLAGS_f90_MESH} -c -o $@ $<
