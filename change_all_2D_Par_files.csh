#!/bin/csh 

# I get the list of files using 
# alias list_all_2D_Par_files_to_change " ls DATA/Par_file EXAMPLES/*/Par_fi* EXAMPLES/*/*/Par_fi* "
#
# and then I manually remove EXAMPLES/unused_older_examples_DATA_to_sort  , which is by definition in an obsolete/incompatible format
#

foreach file ( DATA/Par_file EXAMPLES/Abel_Brest/Par_file_Abel_Balanche_bathy_source_solid EXAMPLES/acoustic_poroelastic/Par_file_acoustic_poroelastic EXAMPLES/attenuation/Par_file_attenuation_2D EXAMPLES/canyon/Par_file_canyon EXAMPLES/fluid_solid/fluid_solid_external_mesh/Par_file_fluid_solid EXAMPLES/fluid_solid/from_2000_Geophysics_paper_flat_ocean_bottom/Par_file_fluid_solid EXAMPLES/fluid_solid/from_2000_Geophysics_paper_sinusoidal_ocean_bottom/Par_file_fluid_solid EXAMPLES/Gmsh_example_MPI/Par_file_Gmsh_SqrCircles.in EXAMPLES/Gmsh_example_serial/Par_file_Gmsh_SqrCircles.in EXAMPLES/INDUSTRIAL_FORMAT/Par_file EXAMPLES/init_plane/Par_file_Slave EXAMPLES/init_plane/Par_file_Slave_for EXAMPLES/init_plane/Par_file_Slave_kernel EXAMPLES/M2_UPPA/Par_file_M2_UPPA EXAMPLES/noise_layered/model_0/Par_file_fair EXAMPLES/noise_layered/model_0/Par_file_good EXAMPLES/noise_layered/model_1/Par_file_best EXAMPLES/noise_layered/model_1/Par_file_fair EXAMPLES/noise_layered/model_1/Par_file_good EXAMPLES/noise_layered/model_2/Par_file_fair EXAMPLES/noise_layered/model_2/Par_file_good EXAMPLES/noise_uniform/Par_file_noise_1 EXAMPLES/noise_uniform/Par_file_noise_2 EXAMPLES/noise_uniform/Par_file_noise_3 EXAMPLES/Rayleigh_wave_no_crack/Par_file_Rayleigh_2D EXAMPLES/Rayleigh_wave_with_crack/Par_file_Rayleigh_2D EXAMPLES/semi_infinite_homo/Par_file_elastic_2D EXAMPLES/Tape2007_kernel/Par_file_Tape2007_onerec EXAMPLES/Tape2007/Par_file_Tape2007_132rec_checker EXAMPLES/Tape2007/Par_file_Tape2007_onerec EXAMPLES/Tromp2005_kernel/Par_file_Tromp2005 EXAMPLES/Tromp2005/Par_file_Tromp2005 EXAMPLES/Tromp2005/Par_file_Tromp2005_s100 EXAMPLES/unused_older_examples_DATA_to_sort/Par_file_no_canyon )

set newfile = "../"$file

# DK DK leave the white spaces, in order to have the right alignment with other variables
###sed -e '1,$s/PERFORM_CUTHILL_MCKEE           = .true.         # perform inverse Cuthill-McKee (1969) optimization\/permutation for mesh numbering/PERFORM_CUTHILL_MCKEE           = .false.        # perform inverse Cuthill-McKee (1969) optimization\/permutation for mesh numbering (can be very costly and not very useful)/g' < $newfile > __________temp_27_zzzyyyyxxxx__________
sed -e '1,$s/PERFORM_CUTHILL_MCKEE           = .true.         # perform inverse Cuthill-McKee (1969) optimization\/permutation for mesh numbering/PERFORM_CUTHILL_MCKEE           = .false.        # perform inverse Cuthill-McKee (1969) optimization\/permutation for mesh numbering (can be very costly and not very useful)/g' < $newfile > __________temp_27_zzzyyyyxxxx__________

clean_a_Par_file.py __________temp_27_zzzyyyyxxxx__________ $file

#mv __________temp_27_zzzyyyyxxxx__________ $newfile

end

