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
	$(EMPTY_MACRO)

${E}/xcombine_sem: $(xcombine_sem_OBJECTS) $(xcombine_sem_SHARED_OBJECTS)
	${FCLINK} -o $@ $+

#######################################

##
## smooth_sum
##

xsmooth_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/smooth_sem.postprocess.o \
	$O/parse_kernel_names.postprocess.o \
	$O/gll_library.spec.o \
	$O/exit_mpi.spec.o \
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
INFO_CUDA="building xsmooth_sem with CUDA support"
else
## non-cuda version
xsmooth_sem_OBJECTS += $(cuda_smooth_sem_STUBS)
## libs
xsmooth_sem_LIBS = $(MPILIBS)
INFO_CUDA="building xsmooth_sem without CUDA support"
endif


${E}/xsmooth_sem: $(xsmooth_sem_OBJECTS)
	@echo ""
	@echo $(INFO_CUDA)
	@echo ""
	$(FCLINK) -o ${E}/xsmooth_sem $(xsmooth_sem_OBJECTS) $(xsmooth_sem_LIBS)
	@echo ""

#######################################


####
#### rule for each .o file below
####

##
## postprocess
##

$O/%.postprocess_module.o: $S/%.f90 ${SETUP}/constants.h $O/specfem2D_par.spec.o
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

