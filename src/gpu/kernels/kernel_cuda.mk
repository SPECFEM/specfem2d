cuda_kernels_OBJS := \
	$O/add_sources_ac_SIM_TYPE_2_OR_3_kernel.cuda-kernel.o \
	$O/add_sources_el_SIM_TYPE_2_OR_3_kernel.cuda-kernel.o \
	$O/assemble_boundary_accel_on_device_kernel.cuda-kernel.o \
	$O/assemble_boundary_potential_on_device_kernel.cuda-kernel.o \
	$O/compute_acoustic_seismogram_kernel.cuda-kernel.o \
	$O/compute_add_moving_sources_acoustic_kernel.cuda-kernel.o \
	$O/compute_add_sources_acoustic_kernel.cuda-kernel.o \
	$O/compute_add_sources_kernel.cuda-kernel.o \
	$O/compute_coupling_acoustic_el_kernel.cuda-kernel.o \
	$O/compute_coupling_elastic_ac_kernel.cuda-kernel.o \
	$O/compute_elastic_seismogram_kernel.cuda-kernel.o \
	$O/compute_gradient_kernel.cuda-kernel.o \
	$O/compute_kernels_acoustic_kernel.cuda-kernel.o \
	$O/compute_kernels_cuda_kernel.cuda-kernel.o \
	$O/compute_kernels_hess_ac_cuda_kernel.cuda-kernel.o \
	$O/compute_stacey_acoustic_kernel.cuda-kernel.o \
	$O/compute_stacey_elastic_kernel.cuda-kernel.o \
	$O/enforce_free_surface_cuda_kernel.cuda-kernel.o \
	$O/get_maximum_kernel.cuda-kernel.o \
	$O/get_maximum_vector_kernel.cuda-kernel.o \
	$O/Kernel_2_acoustic_impl.cuda-kernel.o \
	$O/Kernel_2_viscoelastic_impl.cuda-kernel.o \
	$O/kernel_3_accel_cuda_device.cuda-kernel.o \
	$O/kernel_3_acoustic_cuda_device.cuda-kernel.o \
	$O/kernel_3_cuda_device.cuda-kernel.o \
	$O/kernel_3_veloc_cuda_device.cuda-kernel.o \
	$O/pml_boundary_acoustic_cuda_kernel.cuda-kernel.o \
	$O/prepare_boundary_accel_on_device_kernel.cuda-kernel.o \
	$O/prepare_boundary_potential_on_device_kernel.cuda-kernel.o \
	$O/process_smooth.cuda-kernel.o \
	$O/UpdateDispVeloc_kernel.cuda-kernel.o \
	$O/UpdatePotential_kernel.cuda-kernel.o \
	$(EMPTY_MACRO)

