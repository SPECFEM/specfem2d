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
S := ${S_TOP}/src/cuda
$(cuda_OBJECTS): S = ${S_TOP}/src/cuda

#######################################

cuda_TARGETS = \
	$(cuda_OBJECTS) \
	$(EMPTY_MACRO)

cuda_OBJECTS = \
	$O/assemble_MPI_scalar_cuda.cuda.o \
	$O/assemble_MPI_vector_cuda.cuda.o \
	$O/check_fields_cuda.cuda.o \
	$O/compute_add_sources_viscoacoustic_cuda.cuda.o \
	$O/compute_add_sources_viscoelastic_cuda.cuda.o \
	$O/compute_coupling_cuda.cuda.o \
	$O/compute_forces_acoustic_cuda.cuda.o \
	$O/compute_forces_viscoelastic_cuda.cuda.o \
	$O/compute_kernels_cuda.cuda.o \
	$O/compute_stacey_acoustic_cuda.cuda.o \
	$O/compute_stacey_viscoelastic_cuda.cuda.o \
	$O/initialize_cuda.cuda.o \
	$O/prepare_mesh_constants_cuda.cuda.o \
	$O/transfer_fields_cuda.cuda.o \
	$O/update_displacement_cuda.cuda.o \
	$O/write_seismograms_cuda.cuda.o \
	$(EMPTY_MACRO)


cuda_STUBS = \
	$O/specfem2D_gpu_cuda_method_stubs.cudacc.o \
	$(EMPTY_MACRO)


cuda_DEVICE_OBJ = \
	$O/cuda_device_obj.o \
	$(EMPTY_MACRO)


#######################################

####
#### rule for each .o file below
####

###
### CUDA compilation
###

$O/%.cuda.o: $S/%.cu ${SETUP}/config.h $S/mesh_constants_cuda.h $S/prepare_constants_cuda.h
	${NVCC} -c $< -o $@ $(NVCC_FLAGS)

$O/%.cudacc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<

