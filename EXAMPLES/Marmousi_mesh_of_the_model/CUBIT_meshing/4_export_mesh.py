#!/usr/bin/env python-cubit
#
# creates block entities to export mesh with the cubit2specfem2d.py script
#
from __future__ import print_function
import os
import sys

# to run this script from command line, python must have its PATH environment set such that it
# includes the path to CUBIT/Trelis (cubit.py).
#
# you can also explicitly set it here e.g. like:
#sys.path.append('/opt/Trelis-15.0/bin/')

try:
    import cubit
except ImportError:
    print("Error: Importing cubit as python module failed")
    print("could not import cubit, please check your PYTHONPATH settings...")
    print("")
    print("current path: ")
    print(sys.path)
    print("")
    print("try to include path to directory which includes file cubit.py, e.g. /opt/Trelis-15.0/bin/")
    print("")
    sys.exit("Import cubit failed")

#debug
#print(sys.path)

# initializes cubit
cubit.init(["-noecho","-nojournal"])

# gets version string
cubit_version = cubit.get_version()
print("# version: ",cubit_version)

# extracts major number
v = cubit_version.split('.')
cubit_version_major = int(v[0])
cubit_version_minor = int(v[1])
print("# major version number: ",cubit_version_major)

# current work directory
cubit.cmd('pwd')

# Creating the volumes
cubit.cmd('reset')

#checks if new file available
if not os.path.exists("mesh_3_marmousi_smoothed.cub"):
    print("#")
    print("# file mesh_3_marmousi_smoothed.cub not found, please check ...")
    print("#")
    cubit.cmd('pause')
cubit.cmd('open "mesh_3_marmousi_smoothed.cub"')

# infos
# get 2D plane of coordinates
# determined based on bounding box dimensions
list_vol = cubit.parse_cubit_list("volume","all")
print("# volume list length = ",len(list_vol))

xmin_box = cubit.get_total_bounding_box("volume",list_vol)[0]
xmax_box = cubit.get_total_bounding_box("volume",list_vol)[1]

ymin_box = cubit.get_total_bounding_box("volume",list_vol)[3]
ymax_box = cubit.get_total_bounding_box("volume",list_vol)[4]

zmin_box = cubit.get_total_bounding_box("volume",list_vol)[6] #it is the z_min of the box ... box= xmin,xmax,d,ymin,ymax,d,zmin...
zmax_box = cubit.get_total_bounding_box("volume",list_vol)[7]

print("# bounding box xmin/xmax = ",xmin_box,xmax_box)
print("# bounding box ymin/ymax = ",ymin_box,ymax_box)
print("# bounding box zmin/zmax = ",zmin_box,zmax_box)
print("#")

# plane identifier: 1 == XZ-plane, 2 == XY-plane, 3 == YZ-plane
if abs(ymax_box-ymin_box) < 0.001:
    print("# model in XZ plane")
    plane_id = 1
elif abs(zmax_box-zmin_box) < 0.001:
    print("# model in XY plane")
    plane_id = 2
elif abs(xmax_box-xmin_box) < 0.001:
    print("# model in YZ plane")
    plane_id = 3
else:
    print("# WARNING: model is not 2D, will ignored Y-dimension")
    plane_id = 1
print("#   plane identifier: ",plane_id)
print("#")

# moves 2D model to XZ-plane
# (not necessary anymore as cubit2specfem2d.py script will recognize in which plane the model lies)
if 1 == 0:
    # Get a list of all vertices
    node_list = cubit.parse_cubit_list('node','all') # Import all the nodes of the model
    print("# number of nodes: ",len(node_list))

    # Iterate over all vertices and swap Y and Z coordinates
    cubit.cmd('set echo off')
    cubit.cmd('set node constraint off')
    for vertex in node_list:
        # get its coordinates (3 coordinates even for a 2D model in cubit)
        x, y, z = cubit.get_nodal_coordinates(vertex)
        if plane_id == 1:
            # YZ -> XZ
            tmp = x
            x = y
            y = tmp
        elif plane_id == 2:
            # XZ all ok
            pass

        elif plane_id == 3:
            # XY -> XZ
            tmp = z
            z = y
            y = tmp
        else:
            print("# Error: plane id not recognized :",plane_id)
            cubit.cmd('pause')
        cubit.cmd('node '+str(vertex)+'  move X '+str(x)+' Y '+str(y)+' Z '+str(z))

    # check with new bounding box
    xmin_box = cubit.get_total_bounding_box("volume",list_vol)[0]
    xmax_box = cubit.get_total_bounding_box("volume",list_vol)[1]

    ymin_box = cubit.get_total_bounding_box("volume",list_vol)[3]
    ymax_box = cubit.get_total_bounding_box("volume",list_vol)[4]

    zmin_box = cubit.get_total_bounding_box("volume",list_vol)[6] #it is the z_min of the box ... box= xmin,xmax,d,ymin,ymax,d,zmin...
    zmax_box = cubit.get_total_bounding_box("volume",list_vol)[7]

    print("# updated bounding box xmin/xmax = ",xmin_box,xmax_box)
    print("# updated bounding box ymin/ymax = ",ymin_box,ymax_box)
    print("# updated bounding box zmin/zmax = ",zmin_box,zmax_box)


# surface count
nsurf = cubit.get_surface_count()
print("#")
print("# number of surfaces = ",nsurf)
if nsurf != 435:
    print("Error: wrong number of surfaces, must be 435 for this Marmousi model")
    sys.exit(1)
else:
    print("# surface number okay")
print("#")

# gets list of surface ids
surfaces = cubit.parse_cubit_list("surface", "all")
if len(surfaces) != nsurf:
    print("Error: wrong number of surface ids received ",len(surfaces))
    sys.exit(1)

# total number of quads in mesh
num_elems = cubit.get_quad_count()
print("#")
print("# number of quad elements: ",num_elems)
print("#")

# assign material properties as block attributes
print("#")
print("# material properties: assign block attributes")
print("#")
print("#### DEFINE MATERIAL PROPERTIES #######################")
print("#")

# single block for all surfaces
if 1 == 0:
    print("# single block material definition")
    print("#")

    # sets the id of the volume block
    id_block = 1

    print("# cubit block:")
    print("#  volume block id = " + str(id_block))
    print("#")

    # Define material properties
    cubit.cmd('block '+str(id_block)+' face in volume all')
    cubit.cmd('block '+str(id_block)+' element type QUAD4') # convert if elements types are SHELL4 to QUAD4

    # elastic material
    if 1 == 1:
        print("# using elastic material definition")
        cubit.cmd('block '+str(id_block)+' name "elastic 1" ')          # elastic material region
        cubit.cmd('block '+str(id_block)+' attribute count 6')
        cubit.cmd('block '+str(id_block)+' attribute index 1 1')        # flag for material: 1 for 1. material
        cubit.cmd('block '+str(id_block)+' attribute index 2 2800')     # vp
        cubit.cmd('block '+str(id_block)+' attribute index 3 1500')     # vs
        cubit.cmd('block '+str(id_block)+' attribute index 4 2300')     # rho
        cubit.cmd('block '+str(id_block)+' attribute index 5 9000.0')   # Qmu
        cubit.cmd('block '+str(id_block)+' attribute index 6 0')        # anisotropy_flag

    # acoustic material
    if 1 == 0:
        print("# using acoustic material definition")
        cubit.cmd('block '+str(id_block)+' name "acoustic 1" ')         # acoustic material region
        cubit.cmd('block '+str(id_block)+' attribute count 4')
        cubit.cmd('block '+str(id_block)+' attribute index 1 1  ')      # material 1
        cubit.cmd('block '+str(id_block)+' attribute index 2 1480 ')    # vp
        cubit.cmd('block '+str(id_block)+' attribute index 3 0 ')       # vs
        cubit.cmd('block '+str(id_block)+' attribute index 4 1028 ')    # rho (ocean salt water density:

# creates for each surface its own block
if 1 == 1:
    # note: having a block definition per surface allows to create a MESH/materials_file with different material ids per surface.
    #       this will help when creating a new nummaterial_velocity_file for Marmousi2 with different properties for each surface.
    print("# each surface with its own block material definition")
    print("#")
    print("# using elastic material definition")    
    print("#")
    
    for i in surfaces:
        isurf = i  # surface ids start at 1
        print("# surface: ",isurf)

        # sets the id of the block
        id_block = isurf   # using same block number as surface number

        # Define material properties
        cubit.cmd('block '+str(id_block)+' face in surface '+str(isurf))
        cubit.cmd('block '+str(id_block)+' element type QUAD4') # convert if elements types are SHELL4 to QUAD4

        # checks with number of quad faces in this block
        # no empty blocks or surfaces should appear in mesh
        quads = cubit.get_block_faces(id_block)
        print("#   number of faces: ",len(quads))
        if len(quads) == 0:
            print("Error: surface {} has no mesh quads/faces {}".format(isurf,quads))
            sys.exit(1)

        # at this point, it doesn't matter to define the same properties for all blocks.
        # we will create a new nummaterial_velocity_file using the script ./convert_surface_rock_to_velocities.py
        # with specific velocity properties for each material/surface
        cubit.cmd('block '+str(id_block)+' name "elastic '+str(isurf)+'" ')          # elastic material region
        cubit.cmd('block '+str(id_block)+' attribute count 6')
        cubit.cmd('block '+str(id_block)+' attribute index 1 '+str(isurf))  # flag for material: 1 for 1. material, 2 for 2. material, ..
        cubit.cmd('block '+str(id_block)+' attribute index 2 2800')     # vp
        cubit.cmd('block '+str(id_block)+' attribute index 3 1500')     # vs
        cubit.cmd('block '+str(id_block)+' attribute index 4 2300')     # rho
        cubit.cmd('block '+str(id_block)+' attribute index 5 9999.0')   # Qmu
        cubit.cmd('block '+str(id_block)+' attribute index 6 0')        # anisotropy_flag

### Remove duplicates and fix problems
# not necessary...
#print("")
#print("cleaning mesh:")
#print("")
#cubit.cmd('set dev on') # Set development commands on (to delete)
#cubit.cmd('delete face all except face in block all') # This is necessary to remove duplicated faces
#cubit.cmd('compress')   # Remove gap in numbering

# backup cubit
cubit.cmd('save as "MESH/final_mesh.cub" overwrite')

print("#")
print("# meshing done: see final cubit file MESH/final_mesh.cub")
print("#")

# export to specfem
#
# adds path to scripts (if not setup yet)
sys.path.append('../../../utils/cubit2specfem2d')

## creates blocks with edges lying on boundaries
## for determining absorbing boundaries and free surface
try:
    import boundary_definition
except:
    print("Error: Importing script boundary_definition failed")
    print("could not import from `../../../utils/cubit2specfem2d` , please check your working dir and PYTHONPATH settings...")
    print("")
    print("current path: ")
    print(sys.path)
    print("")
    print("try to include path to directory which includes file boundary_definition.py and cubit2specfem2d.py, e.g. SPECFEM2D/utils/cubit2specfem2d/")
    print("")
    sys.exit("Import boundary_definition failed")

print("#")
print("#### DEFINE BOUNDARIES #######################")
print("#")
print("# creating boundary blocks")
print("#")

# avoids assigning empty blocks
cubit.cmd('set duplicate block elements on')

print("# boundary definition:")
# bounding edges
boundary_definition.define_bc_edges()

# creates MESH/ directory for file output
os.system('mkdir -p MESH')

# conversion of mesh to specfem-format
try:
    import cubit2specfem2d
except:
    print("Error: Importing script cubit2specfem2d failed")
    print("could not import from `../../../utils/cubit2specfem2d` , please check your working dir and PYTHONPATH settings...")
    print("")
    print("current path: ")
    print(sys.path)
    print("")
    print("try to include path to directory which includes file boundary_definition.py and cubit2specfem2d.py, e.g. SPECFEM2D/utils/cubit2specfem2d/")
    print("")
    sys.exit("Import cubit2specfem2d failed")

print("#")
print("#### MESH EXPORT #######################")
print("#")
print("# exporting to SPECFEM2D-format:")
print("#")
# Export to SPECFEM2D format
#from cubit2specfem2d import mtools,block_tools,mesh_tools,mesh
#profile = mesh() # Store the mesh from Cubit
#profile.write("./MESH") # Write it into files (in specfem2d format)
#or
cubit2specfem2d.export2SPECFEM2D('MESH/')

# screen shot
# (could crash version < 16.4)
if cubit_version_major >= 16 and cubit_version_minor >= 4:
    if cubit.graphics_enabled():
        print("#")
        print("# screenshot")
        print("#")
        print("# mesh visible: ",cubit.is_mesh_visibility_on())
        if plane_id == 1:
            # XZ-plane
            cubit.cmd('view top')
        elif plane_id == 2:
            # XY-plane
            cubit.cmd('view front')
        elif plane_id == 3:
            # YZ-plane
            cubit.cmd('view left')
        cubit.cmd('hardcopy "MESH/image-final-mesh.jpg" jpg')

# all files needed by mesher & solver are now in directory MESH/

print("#")
print("# all mesh files in folder: MESH/")
print("#")

