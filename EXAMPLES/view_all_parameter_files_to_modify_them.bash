#!/bin/bash

# I get the list of files using " ls */Par_fi* */*/Par_fi* */*/*/Par_fi* */*/*/*/Par_fi* | grep -v "~" | tr '\n' ' ' "
# and I add ../DATA/Par_file manually

#
# Note that if you want to change the name of a parameter or anything else in all the Par_files below in an automatic way
# you can also use a "find" command combined with a "sed" command, here is an example:
#
#     find . -type f -name "*Par_file*" -exec sed -i 's/ATTENUATION_VISCOELASTIC_SOLID/ATTENUATION_VISCOELASTIC/g' {} \;
#
# However this will *not* make the change in ../DATA/Par_file, which you will need to do separately.
#

if [ -z "$EDITOR" ]
then
  EDITOR=vi
fi

$EDITOR ../DATA/Par_file anisotropic_zinc_crystal/DATA/Par_file attenuation/Par_file_attenuation_2D axisymmetric_case_AXISYM_option/DATA/Par_file canyon/Par_file_canyon fluid_solid/fluid_solid_external_mesh/Par_file_fluid_solid fluid_solid/from_2000_Geophysics_paper_flat_ocean_bottom/Par_file_fluid_solid fluid_solid/from_2000_Geophysics_paper_sinusoidal_ocean_bottom/Par_file_fluid_solid global_Earth_ak135f/DATA/Par_file Gmsh_example_Stacey_MPI/DATA/Par_file Gmsh_example_Stacey_MPI/DATA/Par_file.serial infinite_homogeneous_moment_tensor_vertical_dip_slip/Par_file_elastic_2D infinite_homogeneous_plane_wave/Par_file_elastic_2D initial_plane_wave_with_free_surface/Par_file_Slave initial_plane_wave_with_free_surface/Par_file_Slave_for initial_plane_wave_with_free_surface/Par_file_Slave_kernel LuoYang_fluid_solid_kernel/DATA/Par_file noise_uniform/DATA/Par_file_noise_1 noise_uniform/DATA/Par_file_noise_2 noise_uniform/DATA/Par_file_noise_3 noise_uniform/REF_SEIS/Par_file poroelastic_acoustic/DATA/Par_file poroelastic_semi_infinite_homogeneous/DATA/Par_file Rayleigh_wave_no_crack/Par_file_Rayleigh_2D Rayleigh_wave_with_crack/Par_file_Rayleigh_2D salt_dome_CUBIT_mesh/CPML_homogeneous/Par_file salt_dome_CUBIT_mesh/CPML_normal_fluid_solid/Par_file salt_dome_CUBIT_mesh/CPML_normal_solid_only/Par_file salt_dome_CUBIT_mesh/Stacey_homogeneous/Par_file salt_dome_CUBIT_mesh/Stacey_normal_fluid_solid/Par_file semi_infinite_homogeneous/DATA/Par_file simple_topography_and_also_a_simple_fluid_layer/DATA/Par_file simple_topography_and_also_a_simple_fluid_layer/DATA/Par_file_with_PML simple_topography_and_also_a_simple_fluid_layer/REF_SEIS/Par_file simple_topography_and_also_a_simple_fluid_layer/REF_SEIS_with_PML/Par_file Tape2007_kernel/Par_file_Tape2007_onerec Tape2007/Par_file_Tape2007_132rec_checker Tape2007/Par_file_Tape2007_onerec thermocline/Par_file_Abel_Balanche_bathy_source_solid Tromp2005/DATA/Par_file Tromp2005/DATA/Par_file_Tromp2005_s100 Tromp2005_kernel/DATA/Par_file ZZZ_currently_broken_or_obsolete_examples_but_do_not_remove/INDUSTRIAL_FORMAT/Par_file ZZZ_currently_broken_or_obsolete_examples_but_do_not_remove/older_attenuation_with_Carcione_1988_1_over_N_problem/Par_file_attenuation_2D
