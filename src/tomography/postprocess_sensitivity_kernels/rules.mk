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
S := ${S_TOP}/src/tomography/postprocess_sensitivity_kernels
$(tomography/postprocess_sensitivity_kernels_OBJECTS): S := ${S_TOP}/src/tomography/postprocess_sensitivity_kernels

#######################################

tomography/postprocess_sensitivity_kernels_TARGETS = \
	$E/xcombine_sem \
	$E/xsmooth_sem \
	$(EMPTY_MACRO)

tomography/postprocess_sensitivity_kernels_OBJECTS = \
	$(xcombine_sem_OBJECTS) \
	$(xsmooth_sem_OBJECTS) \
	$(EMPTY_MACRO)

tomography/postprocess_sensitivity_kernels_MODULES = \
	$(FC_MODDIR)/postprocess_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

tomography/postprocess_sensitivity_kernels_SHARED_OBJECTS = \
	$(xcombine_sem_SHARED_OBJECTS) \
	$(xsmooth_sem_SHARED_OBJECTS) \
	$(EMPTY_MACRO)


#######################################

####
#### rules for executables
####

postprocess: $(tomography/postprocess_sensitivity_kernels_TARGETS)

combine_sem: xcombine_sem
xcombine_sem: $E/xcombine_sem

smooth_sem: xsmooth_sem
xsmooth_sem: $E/xsmooth_sem

#######################################

####
#### rules for each program follow
####

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

${E}/xcombine_sem: $(xcombine_sem_OBJECTS) $(xcombine_sem_SHARED_OBJECTS)
	@echo ""
	@echo "building xcombine_sem"
	@echo ""
	${FCLINK} -o $@ $+
	@echo ""

#######################################

##
## smooth_sum
##

xsmooth_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/smooth_sem.postprocess.o \
	$O/parse_kernel_names.postprocess.o \
	$(EMPTY_MACRO)

xsmooth_sem_SHARED_OBJECTS = \
	$O/specfem2D_par.spec_module.o \
	$O/shared_par.shared_module.o \
	$O/exit_mpi.shared.o \
	$O/gll_library.shared.o \
	$O/parallel.shared.o \
	$O/read_parameter_file.mesh.o \
	$O/read_value_parameters.shared.o \
	$O/read_material_table.mesh.o \
	$O/read_interfaces_file.mesh.o \
	$O/read_regions.mesh.o \
	$O/param_reader.cc.o \
	$(EMPTY_MACRO)

cuda_smooth_sem_STUBS = \
	$O/smooth_sem_cuda_stubs.postprocess.o \
	$(EMPTY_MACRO)

cuda_smooth_sem_OBJECTS = \
	$O/smooth_cuda.postprocess.cuda.o \
	$O/check_fields_cuda.cuda.o \
	$O/initialize_cuda.cuda.o \
	$(EMPTY_MACRO)

cuda_smooth_sem_DEVICE_OBJ = \
	$O/cuda_device_smooth_obj.o \
	$(EMPTY_MACRO)

ifeq ($(CUDA),yes)
## cuda version
xsmooth_sem_OBJECTS += $(cuda_smooth_sem_OBJECTS)
ifeq ($(CUDA_PLUS),yes)
xsmooth_sem_OBJECTS += $(cuda_smooth_sem_DEVICE_OBJ)
endif
## libs
xsmooth_sem_LIBS = $(MPILIBS) $(CUDA_LINK)
INFO_CUDA_SEM="building xsmooth_sem with CUDA support"
else
## non-cuda version
xsmooth_sem_OBJECTS += $(cuda_smooth_sem_STUBS)
## libs
xsmooth_sem_LIBS = $(MPILIBS)
INFO_CUDA_SEM="building xsmooth_sem without CUDA support"
endif


${E}/xsmooth_sem: $(xsmooth_sem_OBJECTS) $(xsmooth_sem_SHARED_OBJECTS)
	@echo ""
	@echo $(INFO_CUDA_SEM)
	@echo ""
	$(FCLINK) -o $@ $+ $(xsmooth_sem_LIBS)
	@echo ""

#######################################

###
### Module dependencies
###


####
#### rule for each .o file below
####

##
## postprocess
##

$O/%.postprocess_module.o: $S/%.f90 $O/specfem2D_par.spec_module.o $O/shared_par.shared_module.o
	${F90} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.f90 $O/postprocess_par.postprocess_module.o
	${F90} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.F90 $O/postprocess_par.postprocess_module.o
	${F90} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<


###
### CUDA
###
$O/%.postprocess.cuda.o: $S/%.cu ${SETUP}/config.h $S/smooth_cuda.h
	${NVCC} -c $< -o $@ $(NVCC_FLAGS)

$(cuda_smooth_sem_DEVICE_OBJ): $(cuda_smooth_sem_OBJECTS)
	${NVCCLINK} -o $(cuda_smooth_sem_DEVICE_OBJ) $(cuda_smooth_sem_OBJECTS)

