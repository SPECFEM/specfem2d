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
	$O/specfem2D_par.spec_module.o \
	$O/acoustic_forcing_boundary.spec.o \
	$O/add_acoustic_forcing.spec.o \
	$O/add_manual_crack.spec.o \
	$O/assemble_MPI.spec.o \
	$O/attenuation_model.spec.o \
	$O/axisymmetric_routines.spec.o \
	$O/calendar.spec.o \
	$O/check_grid.spec.o \
	$O/check_stability.spec.o \
	$O/comp_source_time_function.spec.o \
	$O/compute_add_sources_acoustic.spec.o \
	$O/compute_add_sources_poro.spec.o \
	$O/compute_add_sources_viscoelastic.spec.o \
	$O/compute_arrays_source.spec.o \
	$O/compute_attenuation_poro_fluid_part.spec.o \
	$O/compute_attenuation_viscoelastic.spec.o \
	$O/compute_Bielak_conditions.spec.o \
	$O/compute_coupling_acoustic_el.spec.o \
	$O/compute_coupling_acoustic_po.spec.o \
	$O/compute_coupling_poro_ac.spec.o \
	$O/compute_coupling_poro_viscoelastic.spec.o \
	$O/compute_coupling_viscoelastic_ac.spec.o \
	$O/compute_coupling_viscoelastic_po.spec.o \
	$O/compute_curl_one_element.spec.o \
	$O/compute_energy.spec.o \
	$O/compute_forces_viscoacoustic.spec.o \
	$O/compute_attenuation_viscoacoustic.spec.o \
	$O/compute_forces_viscoacoustic_calling_routine.spec.o \
	$O/compute_forces_poroelastic_calling_routine.spec.o \
	$O/compute_forces_poro_fluid.spec.o \
	$O/compute_forces_poro_solid.spec.o \
	$O/compute_forces_poro_viscous_damping.spec.o \
	$O/compute_forces_viscoelastic.spec.o \
	$O/compute_forces_viscoelastic_calling_routine.spec.o \
	$O/compute_gpu_acoustic.spec.o \
	$O/compute_gpu_elastic.spec.o \
	$O/compute_gradient_attenuation.spec.o \
	$O/compute_interpolated_dva.spec.o \
	$O/compute_normal_vector.spec.o \
	$O/compute_pressure.spec.o \
	$O/compute_stacey_acoustic.spec.o \
	$O/compute_stacey_elastic.spec.o \
	$O/compute_stacey_poroelastic.spec.o \
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
	$O/define_external_model_from_marmousi.spec.o \
	$O/enforce_acoustic_free_surface.spec.o \
	$O/enforce_fields.spec.o \
	$O/finalize_simulation.spec.o \
	$O/get_coupling_edges.spec.o \
	$O/get_MPI.spec.o \
	$O/get_global.spec.o \
	$O/get_poroelastic_velocities.spec.o \
	$O/get_simulation_domains.spec.o \
	$O/initialize_simulation.spec.o \
	$O/invert_mass_matrix.spec.o \
	$O/is_in_convex_quadrilateral.spec.o \
	$O/iterate_time.spec.o \
	$O/iterate_time_undoatt.spec.o \
	$O/locate_receivers.spec.o \
	$O/locate_source.spec.o \
	$O/netlib_specfun_erf.spec.o \
	$O/noise_tomography.spec.o \
	$O/paco_beyond_critical.spec.o \
	$O/paco_convolve_fft.spec.o \
	$O/plot_gll.spec.o \
	$O/plot_post.spec.o \
	$O/pml_init.spec.o \
	$O/pml_compute.spec.o \
	$O/pml_compute_accel_contribution.spec.o \
	$O/pml_compute_memory_variables.spec.o \
	$O/prepare_absorb.spec.o \
	$O/prepare_assemble_MPI.spec.o \
	$O/prepare_color_image.spec.o \
	$O/prepare_gpu.spec.o \
	$O/prepare_initial_field.spec.o \
	$O/prepare_pml.spec.o \
	$O/prepare_source_time_function.spec.o \
	$O/prepare_timerun.spec.o \
	$O/prepare_wavefields.spec.o \
	$O/read_materials.spec.o \
	$O/read_mesh_databases.spec.o \
	$O/read_external_model.spec.o \
	$O/read_forward_arrays.spec.o \
	$O/read_save_binary_database.spec.o \
	$O/recompute_jacobian.spec.o \
	$O/save_adjoint_kernels.spec.o \
	$O/save_openDX_jacobian.spec.o \
	$O/set_source_parameters.spec.o \
	$O/setup_GLL_points.spec.o \
	$O/setup_mesh.spec.o \
	$O/setup_sources_receivers.spec.o \
	$O/sort_array_coordinates.spec.o \
	$O/specfem2D.spec.o \
	$O/update_displacement_LDDRK.spec.o \
	$O/update_displacement_Newmark.spec.o \
	$O/update_displacement_RK.spec.o \
	$O/compute_kernels.spec.o \
	$O/write_jpeg_image.cc.o \
	$O/attenuation_compute_param.cc.o \
	$O/write_movie_output.spec.o \
	$O/write_output_SU.spec.o \
	$O/write_seismograms.spec.o \
	$O/write_postscript_snapshot.spec.o \
	$O/write_color_image_snaphot.spec.o \
	$O/write_wavefield_dumps.spec.o \
	$O/save_read_array_for_undoatt.spec.o \
	$(EMPTY_MACRO)

specfem2D_MODULES = \
	$(FC_MODDIR)/constants.$(FC_MODEXT) \
	$(FC_MODDIR)/enforce_par.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_gpu.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_movie.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_noise.$(FC_MODEXT) \
	$(FC_MODDIR)/interpolation.$(FC_MODEXT) \
	$(FC_MODDIR)/model_tomography_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

specfem2D_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/define_shape_functions.shared.o \
	$O/exit_mpi.shared.o \
	$O/force_ftz.cc.o \
	$O/gll_library.shared.o \
	$O/lagrange_poly.shared.o \
	$O/parallel.shared.o \
	$O/read_parameter_file.mesh.o \
	$O/read_value_parameters.shared.o \
	$O/read_material_table.mesh.o \
	$O/read_interfaces_file.mesh.o \
	$O/read_regions.mesh.o \
	$O/param_reader.cc.o \
	$O/set_color_palette.shared.o \
	$(EMPTY_MACRO)

JPEGLIB_OBJECTS = \
	$O/jaricom.cc_jpeg.o \
	$O/jcapimin.cc_jpeg.o \
	$O/jcapistd.cc_jpeg.o \
	$O/jcarith.cc_jpeg.o \
	$O/jccoefct.cc_jpeg.o \
	$O/jccolor.cc_jpeg.o \
	$O/jcdctmgr.cc_jpeg.o \
	$O/jchuff.cc_jpeg.o \
	$O/jcinit.cc_jpeg.o \
	$O/jcmainct.cc_jpeg.o \
	$O/jcmarker.cc_jpeg.o \
	$O/jcmaster.cc_jpeg.o \
	$O/jcomapi.cc_jpeg.o \
	$O/jcparam.cc_jpeg.o \
	$O/jcprepct.cc_jpeg.o \
	$O/jcsample.cc_jpeg.o \
	$O/jctrans.cc_jpeg.o \
	$O/jdapimin.cc_jpeg.o \
	$O/jdapistd.cc_jpeg.o \
	$O/jdarith.cc_jpeg.o \
	$O/jdatadst.cc_jpeg.o \
	$O/jdatasrc.cc_jpeg.o \
	$O/jdcoefct.cc_jpeg.o \
	$O/jdcolor.cc_jpeg.o \
	$O/jddctmgr.cc_jpeg.o \
	$O/jdhuff.cc_jpeg.o \
	$O/jdinput.cc_jpeg.o \
	$O/jdmainct.cc_jpeg.o \
	$O/jdmarker.cc_jpeg.o \
	$O/jdmaster.cc_jpeg.o \
	$O/jdmerge.cc_jpeg.o \
	$O/jdpostct.cc_jpeg.o \
	$O/jdsample.cc_jpeg.o \
	$O/jdtrans.cc_jpeg.o \
	$O/jerror.cc_jpeg.o \
	$O/jfdctflt.cc_jpeg.o \
	$O/jfdctfst.cc_jpeg.o \
	$O/jfdctint.cc_jpeg.o \
	$O/jidctflt.cc_jpeg.o \
	$O/jidctfst.cc_jpeg.o \
	$O/jidctint.cc_jpeg.o \
	$O/jmemmgr.cc_jpeg.o \
	$O/jmemnobs.cc_jpeg.o \
	$O/jquant1.cc_jpeg.o \
	$O/jquant2.cc_jpeg.o \
	$O/jutils.cc_jpeg.o \
	$(EMPTY_MACRO)

specfem2D_OBJECTS += $(JPEGLIB_OBJECTS)


###
### CUDA
###

cuda_specfem2D_OBJECTS = \
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


cuda_specfem2D_STUBS = \
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

# mostly as example how to specify special dependencies
# the general dependency on the specfem module is handled by the rules below

$O/specfem2D.spec.o: $O/specfem2D_par.spec_module.o

# Version file
$O/initialize_simulation.spec.o: ${SETUP}/version.fh

##
## object files
##

####
#### rule to build each .o file below
####

$O/%.spec_module.o: $S/%.f90 ${SETUP}/constants.h $O/shared_par.shared_module.o
	${F90} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec.o: $S/%.f90 ${SETUP}/constants.h $O/specfem2D_par.spec_module.o
	${F90} ${FCFLAGS_f90} -c -o $@ $<

$O/%.spec.o: $S/%.F90 ${SETUP}/constants.h $O/specfem2D_par.spec_module.o
	${F90} ${FCFLAGS_f90} -c -o $@ $<

$O/%.cc.o: $S/%.c ${SETUP}/config.h
	${CC} ${CFLAGS} -c -o $@ $<

###
### CUDA 5 only
###

$(cuda_specfem2D_DEVICE_OBJ): $(cuda_OBJECTS)
	${NVCCLINK} -o $(cuda_specfem2D_DEVICE_OBJ) $(cuda_OBJECTS)


##
## JPEG library files
##

$O/%.cc_jpeg.o: $S/libjpeg/%.c
	${CC} -c $(CFLAGS) -I${LIBJPEG} -o $@ $<
