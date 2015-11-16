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
S := ${S_TOP}/src/specfem2D
$(specfem2D_OBJECTS): S := ${S_TOP}/src/specfem2D

# libjpeg file directory
LIBJPEG = ${S_TOP}/src/specfem2D/libjpeg
$(JPEGLIB_OBJECTS): S := ${S_TOP}/src/specfem2D/libjpeg


####
#### targets
####

# default targets
specfem2D_TARGETS = \
	$E/xspecfem2D \
	$(EMPTY_MACRO)


specfem2D_OBJECTS = \
	$O/acoustic_forcing_boundary.spec.o \
	$O/assemble_MPI.spec.o \
	$O/attenuation_model.spec.o \
	$O/axisymmetric_routines.spec.o \
	$O/calendar.spec.o \
	$O/check_grid.spec.o \
	$O/check_stability.spec.o \
	$O/compute_Bielak_conditions.spec.o \
	$O/compute_arrays_source.spec.o \
	$O/compute_coupling_acoustic_el.spec.o \
	$O/compute_coupling_acoustic_po.spec.o \
	$O/compute_curl_one_element.spec.o \
	$O/compute_energy.spec.o \
	$O/compute_forces_acoustic.spec.o \
	$O/compute_forces_acoustic_backward.spec.o \
	$O/compute_add_sources_acoustic.spec.o \
	$O/compute_forces_gravitoacoustic.spec.o \
	$O/compute_forces_poro_fluid.spec.o \
	$O/compute_attenuation_poro_fluid_part.spec.o \
	$O/compute_forces_poro_solid.spec.o \
	$O/compute_coupling_poro_ac.spec.o \
	$O/compute_coupling_poro_viscoelastic.spec.o \
	$O/compute_add_sources_poro.spec.o \
	$O/compute_forces_viscoelastic.spec.o \
	$O/compute_forces_viscoelastic_backward.spec.o \
	$O/compute_coupling_viscoelastic_ac.spec.o \
	$O/compute_coupling_viscoelastic_po.spec.o \
	$O/compute_add_sources_viscoelastic.spec.o \
	$O/compute_gradient_attenuation.spec.o \
	$O/compute_normal_vector.spec.o \
	$O/compute_pressure.spec.o \
	$O/compute_vector_field.spec.o \
	$O/construct_acoustic_surface.spec.o \
	$O/convert_time.spec.o \
	$O/create_color_image.spec.o \
	$O/createnum_fast.spec.o \
	$O/createnum_slow.spec.o \
	$O/datim.spec.o \
	$O/define_derivation_matrices.spec.o \
	$O/define_external_model.spec.o \
	$O/define_external_model_from_tomo_file.spec.o \
	$O/enforce_acoustic_free_surface.spec.o \
	$O/exit_mpi.spec.o \
	$O/finalize_simulation.spec.o \
	$O/force_ftz.cc.o \
	$O/get_MPI.spec.o \
	$O/get_global.spec.o \
	$O/get_poroelastic_velocities.spec.o \
	$O/gll_library.spec.o \
	$O/gmat01.spec.o \
	$O/initialize_simulation.spec.o \
	$O/invert_mass_matrix.spec.o \
	$O/is_in_convex_quadrilateral.spec.o \
	$O/iterate_time.spec.o \
	$O/iterate_time_undoatt.spec.o \
	$O/lagrange_poly.spec.o \
	$O/locate_receivers.spec.o \
	$O/locate_source_force.spec.o \
	$O/locate_source_moment_tensor.spec.o \
	$O/netlib_specfun_erf.spec.o \
	$O/noise_tomography.spec.o \
	$O/paco_beyond_critical.spec.o \
	$O/paco_convolve_fft.spec.o \
	$O/plot_gll.spec.o \
	$O/plot_post.spec.o \
	$O/pml_init.spec.o \
	$O/pml_compute.spec.o \
	$O/prepare_absorb.spec.o \
	$O/prepare_assemble_MPI.spec.o \
	$O/prepare_color_image.spec.o \
	$O/prepare_initial_field.spec.o \
	$O/prepare_source_time_function.spec.o \
	$O/prepare_timerun.spec.o \
	$O/prepare_timerun_body.spec.o \
	$O/read_databases.spec.o \
	$O/read_external_model.spec.o \
	$O/recompute_jacobian.spec.o \
	$O/save_adjoint_kernels.spec.o \
	$O/save_openDX_jacobian.spec.o \
	$O/set_sources.spec.o \
	$O/setup_sources_receivers.spec.o \
	$O/sort_array_coordinates.spec.o \
	$O/specfem2D.spec.o \
	$O/specfem2D_par.spec.o \
	$O/update_displacement_scheme.spec.o \
	$O/compute_kernels.spec.o \
	$O/write_jpeg_image.cc.o \
	$O/write_output_SU.spec.o \
	$O/write_seismograms.spec.o \
	$O/write_postscript_snapshot.spec.o \
	$O/write_color_image_snaphot.spec.o \
	$O/write_wavefield_dumps.spec.o \
	$O/save_read_array_for_undoatt.spec.o \
	$(EMPTY_MACRO)

specfem2D_MODULES = \
	$(FC_MODDIR)/constants.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par.$(FC_MODEXT) \
	$(FC_MODDIR)/interpolation.$(FC_MODEXT) \
	$(FC_MODDIR)/model_tomography_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

specfem2D_SHARED_OBJECTS = \
	$O/define_shape_functions.shared.o \
	$(EMPTY_MACRO)

JPEGLIB_OBJECTS = \
	$O/jaricom.cc.o \
	$O/jcapimin.cc.o \
	$O/jcapistd.cc.o \
	$O/jcarith.cc.o \
	$O/jccoefct.cc.o \
	$O/jccolor.cc.o \
	$O/jcdctmgr.cc.o \
	$O/jchuff.cc.o \
	$O/jcinit.cc.o \
	$O/jcmainct.cc.o \
	$O/jcmarker.cc.o \
	$O/jcmaster.cc.o \
	$O/jcomapi.cc.o \
	$O/jcparam.cc.o \
	$O/jcprepct.cc.o \
	$O/jcsample.cc.o \
	$O/jctrans.cc.o \
	$O/jdapimin.cc.o \
	$O/jdapistd.cc.o \
	$O/jdarith.cc.o \
	$O/jdatadst.cc.o \
	$O/jdatasrc.cc.o \
	$O/jdcoefct.cc.o \
	$O/jdcolor.cc.o \
	$O/jddctmgr.cc.o \
	$O/jdhuff.cc.o \
	$O/jdinput.cc.o \
	$O/jdmainct.cc.o \
	$O/jdmarker.cc.o \
	$O/jdmaster.cc.o \
	$O/jdmerge.cc.o \
	$O/jdpostct.cc.o \
	$O/jdsample.cc.o \
	$O/jdtrans.cc.o \
	$O/jerror.cc.o \
	$O/jfdctflt.cc.o \
	$O/jfdctfst.cc.o \
	$O/jfdctint.cc.o \
	$O/jidctflt.cc.o \
	$O/jidctfst.cc.o \
	$O/jidctint.cc.o \
	$O/jmemmgr.cc.o \
	$O/jmemnobs.cc.o \
	$O/jquant1.cc.o \
	$O/jquant2.cc.o \
	$O/jutils.cc.o \
	$(EMPTY_MACRO)

specfem2D_OBJECTS += $(JPEGLIB_OBJECTS)


###
### CUDA
###

cuda_specfem2D_OBJECTS = \
	$O/assemble_MPI_scalar_cuda.cuda.o \
	$O/assemble_MPI_vector_cuda.cuda.o \
	$O/check_fields_cuda.cuda.o \
	$O/compute_add_sources_acoustic_cuda.cuda.o \
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
	$O/acoustic_cuda.spec.o \
	$O/elastic_cuda.spec.o \
	$O/init_host_to_dev_variable.spec.o \
	$O/prepare_timerun_gpu.spec.o \
	$(EMPTY_MACRO)


cuda_specfem2D_STUBS = \
	$O/specfem2D_wrapper_cuda_method_stubs.cudaf90.o \
	$O/specfem2D_gpu_cuda_method_stubs.cudacc.o \
	$(EMPTY_MACRO)

cuda_specfem2D_DEVICE_OBJ = \
	$O/cuda_device_obj.o \
	$(EMPTY_MACRO)

ifeq ($(CUDA),yes)
specfem2D_OBJECTS += $(cuda_specfem2D_OBJECTS)
ifeq ($(CUDA_PLUS),yes)
specfem2D_OBJECTS += $(cuda_specfem2D_DEVICE_OBJ)
endif
else
specfem2D_OBJECTS += $(cuda_specfem2D_STUBS)
endif

#######################################


####
#### rules for executables
####

spec: $(specfem2D_TARGETS)

specfem2D: xspecfem2D
xspecfem2D: $E/xspecfem2D


ifeq ($(CUDA),yes)
## cuda version
ifeq ($(CUDA_PLUS),yes)
## cuda 5x & 6x version
INFO_CUDA="building xspecfem2D with CUDA support"
else
## cuda 4 version
INFO_CUDA="building xspecfem2D with CUDA 4 support"
endif

${E}/xspecfem2D: $(specfem2D_OBJECTS) $(specfem2D_SHARED_OBJECTS)
	@echo ""
	@echo $(INFO_CUDA)
	@echo ""
	$(FCLINK) -o ${E}/xspecfem2D $(specfem2D_OBJECTS) $(specfem2D_SHARED_OBJECTS) $(MPILIBS) $(CUDA_LINK)
	@echo ""

else

## non-cuda version
${E}/xspecfem2D: $(specfem2D_OBJECTS) $(specfem2D_SHARED_OBJECTS)
	@echo ""
	@echo "building xspecfem2D without CUDA support"
	@echo ""
	$(FCLINK) -o ${E}/xspecfem2D $(specfem2D_OBJECTS) $(specfem2D_SHARED_OBJECTS) $(MPILIBS)
	@echo ""

endif




#######################################

###
### Module dependencies
###

$O/acoustic_forcing_boundary.spec.o: $O/specfem2D_par.spec.o
$O/acoutic_cuda.spec.o: $O/specfem2D_par.spec.o
$O/assemble_MPI.spec.o: $O/specfem2D_par.spec.o
$O/attenuation_model.spec.o: $O/specfem2D_par.spec.o
$O/axisymmetric_routines.spec.o: $O/specfem2D_par.spec.o
$O/check_grid.spec.o: $O/specfem2D_par.spec.o
$O/check_stability.spec.o: $O/specfem2D_par.spec.o
$O/compute_add_sources_acoustic.spec.o: $O/specfem2D_par.spec.o
$O/compute_add_sources_poro.spec.o: $O/specfem2D_par.spec.o
$O/compute_add_sources_viscoelastic.spec.o: $O/specfem2D_par.spec.o
$O/compute_arrays_source.spec.o: $O/specfem2D_par.spec.o
$O/compute_attenuation_poro_fluid_part.spec.o: $O/specfem2D_par.spec.o
$O/compute_coupling_acoustic_el.spec.o: $O/specfem2D_par.spec.o
$O/compute_coupling_acoustic_po.spec.o: $O/specfem2D_par.spec.o
$O/compute_coupling_poro_ac.spec.o: $O/specfem2D_par.spec.o
$O/compute_coupling_poro_viscoelastic.spec.o: $O/specfem2D_par.spec.o
$O/compute_coupling_viscoelastic_ac.spec.o: $O/specfem2D_par.spec.o
$O/compute_coupling_viscoelastic_po.spec.o: $O/specfem2D_par.spec.o
$O/compute_curl_one_element.spec.o: $O/specfem2D_par.spec.o
$O/compute_energy.spec.o: $O/specfem2D_par.spec.o
$O/compute_forces_acoustic.spec.o: $O/specfem2D_par.spec.o
$O/compute_forces_acoustic_backward.spec.o: $O/specfem2D_par.spec.o
$O/compute_forces_gravitoacoustic.spec.o: $O/specfem2D_par.spec.o
$O/compute_forces_poro_fluid.spec.o: $O/specfem2D_par.spec.o
$O/compute_forces_poro_solid.spec.o: $O/specfem2D_par.spec.o
$O/compute_forces_viscoelastic.spec.o: $O/specfem2D_par.spec.o
$O/compute_forces_viscoelastic_backward.spec.o: $O/specfem2D_par.spec.o
$O/compute_gradient_attenuation.spec.o: $O/specfem2D_par.spec.o
$O/compute_kernels.spec.o: $O/specfem2D_par.spec.o
$O/compute_pressure.spec.o: $O/specfem2D_par.spec.o
$O/compute_vector_field.spec.o: $O/specfem2D_par.spec.o
$O/construct_acoustic_surface.spec.o: $O/specfem2D_par.spec.o
$O/create_color_image.spec.o: $O/specfem2D_par.spec.o
$O/createnum_fast.spec.o: $O/specfem2D_par.spec.o
$O/createnum_slow.spec.o: $O/specfem2D_par.spec.o
$O/define_derivation_matrices.spec.o: $O/specfem2D_par.spec.o
$O/define_external_model.spec.o: $O/specfem2D_par.spec.o
$O/define_external_model_from_tomo_file.spec.o: $O/specfem2D_par.spec.o
$O/elastic_cuda.spec.o: $O/specfem2D_par.spec.o
$O/enforce_acoustic_free_surface.spec.o: $O/specfem2D_par.spec.o
$O/finalize_simulation.spec.o: $O/specfem2D_par.spec.o
$O/get_MPI.spec.o: $O/specfem2D_par.spec.o
$O/get_global.spec.o: $O/specfem2D_par.spec.o
$O/gmat01.spec.o: $O/specfem2D_par.spec.o
$O/init_host_to_dev_variable.spec.o: $O/specfem2D_par.spec.o
$O/initialize_simulation.spec.o: $O/specfem2D_par.spec.o
$O/invert_mass_matrix.spec.o: $O/specfem2D_par.spec.o
$O/iterate_time.spec.o: $O/specfem2D_par.spec.o
$O/iterate_time_undoatt.spec.o: $O/specfem2D_par.spec.o
$O/locate_receivers.spec.o: $O/specfem2D_par.spec.o
$O/noise_tomography.spec.o: $O/specfem2D_par.spec.o
$O/paco_beyond_critical.spec.o: $O/specfem2D_par.spec.o
$O/plot_gll.spec.o: $O/specfem2D_par.spec.o
$O/plot_post.spec.o: $O/specfem2D_par.spec.o
$O/pml_compute.spec.o: $O/specfem2D_par.spec.o
$O/pml_init.spec.o: $O/specfem2D_par.spec.o
$O/prepare_absorb.spec.o: $O/specfem2D_par.spec.o
$O/prepare_assemble_MPI.spec.o: $O/specfem2D_par.spec.o
$O/prepare_color_image.spec.o: $O/specfem2D_par.spec.o
$O/prepare_initial_field.spec.o: $O/specfem2D_par.spec.o
$O/prepare_source_time_function.spec.o: $O/specfem2D_par.spec.o
$O/prepare_timerun.spec.o: $O/specfem2D_par.spec.o
$O/prepare_timerun_body.spec.o: $O/specfem2D_par.spec.o
$O/prepare_timerun_gpu.spec.o: $O/specfem2D_par.spec.o
$O/read_databases.spec.o: $O/specfem2D_par.spec.o
$O/read_external_model.spec.o: $O/specfem2D_par.spec.o
$O/recompute_jacobian.spec.o: $O/specfem2D_par.spec.o
$O/save_adjoint_kernels.spec.o: $O/specfem2D_par.spec.o
$O/save_read_array_for_undoatt.spec.o: $O/specfem2D_par.spec.o
$O/set_sources.spec.o: $O/specfem2D_par.spec.o
$O/setup_sources_receivers.spec.o: $O/specfem2D_par.spec.o
$O/specfem2D.spec.o: $O/specfem2D_par.spec.o
$O/update_displacement_scheme.spec.o: $O/specfem2D_par.spec.o
$O/write_color_image_snaphot.spec.o: $O/specfem2D_par.spec.o
$O/write_output_SU.spec.o: $O/specfem2D_par.spec.o
$O/write_postscript_snapshot.spec.o: $O/specfem2D_par.spec.o
$O/write_seismograms.spec.o: $O/specfem2D_par.spec.o
$O/write_wavefield_dumps.spec.o: $O/specfem2D_par.spec.o


##
## object files
##

####
#### rule to build each .o file below
####

$O/%.spec.o: $S/%.f90 ${SETUP}/constants.h
	${F90} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec.o: $S/%.F90 ${SETUP}/constants.h
	${F90} ${FCFLAGS_f90} -c -o $@ $<

$O/%.cc.o: $S/%.c ${SETUP}/config.h
	${CC} ${CFLAGS} -c -o $@ $<

###
### CUDA 5 only
###

$(cuda_specfem2D_DEVICE_OBJ): $(cuda_OBJECTS)
	${NVCCLINK} -o $(cuda_specfem2D_DEVICE_OBJ) $(cuda_OBJECTS)


##
## shared
##

$O/%.cc.o: $S/libjpeg/%.c
	${CC} -c $(CFLAGS) -I${LIBJPEG} -o $@ $<
