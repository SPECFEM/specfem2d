#!/bin/bash
#
# create the absorbing and free surface files from the Gmsh file
#

python ../../UTILS/Gmsh/LibGmsh2Specfem_convert_Gmsh_to_Specfem2D_official.py SqrCirc.msh -t F -l A -b A -r A

