#!/bin/bash
#
# Run Cubit/Trelis (http://www.csimsoft.com/trelis) and cubit2specfem2d to create the mesh in a new directory called MESH by default
# This is how it looks:
#  ________________________________v__________________________________
# |   ^         ___________________|______________________________|   | ^
# |   | d1   __/                   ^                              | C | |
# |___v_____/                      d2                             | P | | lz
# |<------>--->                                                   | M | |
# |   l1    l2                                                    | L | |
# |_______________________________________________________________|   | v
# |_____________________________C_P_M_L_______________________________|
#  <------------------------------------------------------------->
#                                   lx
# !! BEFORE RUNNING THAT SCRIPT !!
# !! BEWARE !! 1. This script has to be run from the directory specfem2d/EXAMPLES/paper_axisymmetry_example
#              2. Update variable PATH_TO_CUBIT below
#              3. Read the comments in file JensenMesh.py, you will have to modify one line in that file
#              4. If you rebuild the mesh you will have to change the Par_file. Lines:
#                   mesh_file                       = MESH_SMALL/mesh_file         # file containing the mesh
#                   nodes_coords_file               = MESH_SMALL/nodes_coords_file # file containing the nodes coordinates
#                   materials_file                  = MESH_SMALL/materials_file    # file containing the material number for each element
#                   free_surface_file               = MESH_SMALL/free_surface_file # file containing the free surface
#                   axial_elements_file             = MESH_SMALL/elements_axis          # file containing the axial elements if AXISYM is true
#                   absorbing_surface_file          = MESH_SMALL/absorbing_surface_file # file containing the absorbing surface
#                   acoustic_forcing_surface_file   = ./DATA_PML/MSH/Surf_acforcing_Bottom_enforcing_mesh   # file containing the acoustic forcing surface
#                   CPML_element_file               = MESH_SMALL/elements_cpml_list     # file containing the CPML element numbers
#                   tangential_detection_curve_file = ./DATA/courbe_eros_nodes # file containing the curve delimiting the velocity model
#
#                 Has to be changed to:
#                   mesh_file                       = MESH/mesh_file         # file containing the mesh
#                   nodes_coords_file               = MESH/nodes_coords_file # file containing the nodes coordinates
#                   materials_file                  = MESH/materials_file    # file containing the material number for each element
#                   free_surface_file               = MESH/free_surface_file # file containing the free surface
#                   axial_elements_file             = MESH/elements_axis          # file containing the axial elements if AXISYM is true
#                   absorbing_surface_file          = MESH/absorbing_surface_file # file containing the absorbing surface
#                   acoustic_forcing_surface_file   = ./DATA_PML/MSH/Surf_acforcing_Bottom_enforcing_mesh   # file containing the acoustic forcing surface
#                   CPML_element_file               = MESH/elements_cpml_list     # file containing the CPML element numbers
#                   tangential_detection_curve_file = ./DATA/courbe_eros_nodes # file containing the curve delimiting the velocity model
#

PATH_TO_CUBIT=/opt/Trelis-15.0/bin/trelis # TODO update that line: path to cubit executable

while true; do
    read -p "Have you read all the comments in files make_mesh.sh and JensenMesh.py before running that script? [y|n] " yn
    case $yn in
        [Yy]* ) echo "Great!"; break;;
        [Nn]* ) echo "You should!"; exit;;
        * ) echo "Please answer yes (Y or y) or no (N or n)";;
    esac
done

mkdir -p MESH
rm -rf MESH/*

echo "Making mesh ..."

$PATH_TO_CUBIT -nographics -batch ./JensenMesh.py # >> /dev/null 2>&1

echo "Done! Mesh files has been written in directory MESH"

rm -rf history* free_surface_file materials_file mesh_file nodes_coords_file elements_cpml_list absorbing_surface_file

