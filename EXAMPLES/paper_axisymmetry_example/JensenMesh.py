#!/usr/bin/env python
#######################################################################################
#### Script to create the axisymmetric 2.5d mesh necessary to reproduce
#### Jensen et al. 2007.
#### This script can be played as a journal file in
#### Trelis/Cubit (http://www.csimsoft.com/trelis) directly or through the script:
#### make_mesh.sh.
#### It first run cubit commands and then export the mesh under specfem2d format
#### using utils/cubit2specfem2d/cubit2specfem2d.py
####
####  ________________________________v__________________________________
#### |   ^         ___________________|______________________________|   | ^
#### |   | d1   __/                   ^                              | C | |
#### |___v_____/                      d2                             | P | | lz
#### |<------>--->                                                   | M | |
#### |   l1    l2                                                    | L | |
#### |_______________________________________________________________|   | v
#### |_____________________________C_P_M_L_______________________________|
####  <------------------------------------------------------------->
####                                   lx
####
#### To understand that script:
####   _make sure the script tab (python command line) is open in Trelis/Cubit :
####      Trelis 15.0: Tools -> Options -> Layout --> Show Script Tab
####      Cubit 15.1: Tools -> Options -> Layout --> Cubit Layout --> Show Script Tab
####
####   _copy paste one line of this script after another in this script tab and
####    observe what happens.
####
#### !!! Beware !!! There are some (annoying) differences between cubit built-in Python
#### and the actual Python langage:
####
#### 1. "aString" + 'anotherString' can cause problems even after being stored:
####    a = "aString"
####    b = a + 'anotherString'
####    Example which is not working:
####        pathToMeshDir = pathToSpecfem + 'EXAMPLES/paper_axisymmetry_example/MESH'
####        cubit.cmd('cd \"'+pathToMeshDir+'\"')
####
#### 2. No comments after double dots:
####    Example which is not working:
####        if True: # Just a dummy comment
####            print "Ok!"
####
####          -> This example works without the comment
####
#### 3. os.makedirs("~/aDirectory/") does not work. It creates a directory named "~"
####    !!!!! To remove that do: rm -R ./~ AND NEVER rm -rf ~ !!!!!
####
#### 4. sys.argv can not be used
####
#### 5. No comments """ """ at the beginning of a script
####
#### Author alexis dot bottero at gmail dot com
#######################################################################################
from __future__ import print_function
############################ IMPORT MODULES ############################
# Module to interact with cubit:
import cubit
# System modules:
import sys,os
# Locate specfem2d:
pathToSpecfem = '/home/bottero/specfem2d' # TODO update that line: do not use "~"
# If the path given don't finish by slash we add one:
if pathToSpecfem[-1] != '/':
    pathToSpecfem = pathToSpecfem+'/'
# Module to convert the mesh to specfem2d format:
sys.path.append(pathToSpecfem+'utils/cubit2specfem2d')
from cubit2specfem2d import mtools,block_tools,mesh_tools,mesh
########################################################################

############################# SET GEOMETRY #############################
# All lengths are set in meters
# See the schema above

lx = 15000.0
lz = 3000.0
d1 = 600.0
d2 = 100.0
l1 = 2000.0
l2 = 6000.0

############################### SET MESH ###############################
nPml = 3 # PMLs Thickness (number of elements in pmls, 3 is the best option)
elementSize = 70 # Mean element size in meters
# Path where the files have to be written (TODO update this line if wanted)
#pathToMeshDir = pathToSpecfem + 'EXAMPLES/paper_axisymmetry_example/MESH'

########################################################################
#################### NOTHING HAS TO BE CHANGED BELOW ###################
########################################################################

### Create pathToMeshDir if needed
#if not os.path.exists(pathToMeshDir):
#    os.makedirs(pathToMeshDir)
#    print 'Directory '+pathToMeshDir+' created!'

### Change Cubit working directory :
#cubit.cmd('cd \"'+pathToMeshDir+'\"')

### Reset
cubit.cmd('reset')

### Create the rectangle :
cubit.cmd('create surface rectangle width '+str(lx)+' height '+str(lz)+' yplane') # Create a lx*lz rectangle
cubit.cmd('volume 1 move x '+str(lx/2.0)+' z '+str(-lz/2.0)) # To get correct coordinates

### Cut the rectangle along sea bottom :
cubit.cmd('create vertex 0 0 '+str(-d1)+' on curve 2')
cubit.cmd('create vertex '+str(l1)+' 0 '+str(-d1)+' on surface 1')
cubit.cmd('create vertex '+str(l2)+' 0 '+str(-d2)+' on surface 1')
cubit.cmd('create vertex '+str(lx)+' 0 '+str(-d2)+' on curve 4')
cubit.cmd('create curve polyline vertex 5 6 7 8')
#cubit.cmd('partition create surface 1 curve 5 6 7') # Too bad... that create virtual surfaces
cubit.cmd('create surface vertex 3 5 6 10 8 4')
cubit.cmd('subtract volume 2 from volume 1 imprint  keep')
cubit.cmd('delete body 1')

### Clean
cubit.cmd('compress') # Fill the gaps in the numbering
cubit.cmd('color body 1 deepskyblue')
cubit.cmd('color body 2 orange')

# Imprinting is the process of projecting curves from one surface onto an overlapping surface. Merging is
# the process of taking two overlapping surfaces and merging them into one surface shared by two
# volumes, creating non-manifold geometry. Both imprinting and merging are necessary to make adjacent
# volumes have identical meshes at their intersection.

cubit.cmd('merge vol all')
cubit.cmd('imprint vol all')

### Meshing the surfaces :
cubit.cmd('surface 1 2 scheme pave')
cubit.cmd('surface 1 2 size '+str(elementSize)) # Define the size of an element in surfaces
cubit.cmd('mesh surface 1 2')

### End of meshing

### Define material properties
cubit.cmd('####################### DEFINE MATERIAL PROPERTIES #######################')
cubit.cmd('block 1 face in surface 1')
cubit.cmd('block 1 name "Water acoustic"')    # acoustic region
cubit.cmd('block 1 attribute count 1')        # number of attributes
cubit.cmd('block 1 attribute index 1 2')      # material index
cubit.cmd('block 1 element type QUAD4')       # element type. Must be QUAD4 for cubit2specfem2d

cubit.cmd('block 2 face in surface 2')
cubit.cmd('block 2 name "Elastic bottom"')    # elastic region
cubit.cmd('block 2 attribute count 1')        # number of attributes
cubit.cmd('block 2 attribute index 1 1')      # material index
cubit.cmd('block 2 element type QUAD4')       # element type. Must be QUAD4 for cubit2specfem2d

### Create PMLs elements (it is far better to create them directly before meshing but I did that long time ago!):
# z pmls :
cubit.cmd('create element extrude edge in curve 12 direction 0 0 -1 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 100 (face all except face in block 1 2)')
# x pml acoustic :
cubit.cmd('create element extrude edge in curve 8 with x_coord > '+str(lx-elementSize/100.0)+' direction 1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 101 (face all except face in block 1 2 100)')
# x pml elastic :
cubit.cmd('create element extrude edge in face all in block 100 2 with x_coord > '+str(lx-elementSize/100.0)+' direction 1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 102 (face all except face in block 1 2 100 101)')
# xz pmls :
cubit.cmd('set duplicate block elements on')
cubit.cmd('block 103 face in block 102 with z_coord < '+str(-lz))

### Create geometry from free mesh (necessary to merge) :
cubit.cmd('create mesh geometry Face in block 100 feature_angle 135.0 keep')
cubit.cmd('create mesh geometry Face in block 101 feature_angle 135.0 keep')
cubit.cmd('create mesh geometry Face in block 102 feature_angle 135.0 keep')
cubit.cmd('create mesh geometry Face in block 103 feature_angle 135.0 keep')

### Merging
cubit.cmd('merge all')
cubit.cmd('delete block 100 101 102 103')

########################################################################
####### FROM HERE WE CREATE BLOCKS GATHERING ELEMENTS, THEY HAVE #######
##### TABULATED NAME ALLOWING THEM TO BE READ BY cubit2specfem2d.py ####
########################################################################

### Create PMLs blocks
cubit.cmd('block 1000 face in surf 3')
cubit.cmd('block 1000 name "pml_z_elast"')
cubit.cmd('block 1000 attribute count 1')      # number of attributes
cubit.cmd('block 1000 attribute index 1 1')    # material index
cubit.cmd('block 1000 element type QUAD4')

cubit.cmd('block 1001 face in surf 4')
cubit.cmd('block 1001 name "pml_x_acoust"')
cubit.cmd('block 1001 attribute count 1')      # number of attributes
cubit.cmd('block 1001 attribute index 1 2')    # material index
cubit.cmd('block 1001 element type QUAD4')

cubit.cmd('block 1002 face in surf 6')
cubit.cmd('block 1002 name "pml_x_elast"')
cubit.cmd('block 1002 attribute count 1')      # number of attributes
cubit.cmd('block 1002 attribute index 1 1')    # material index
cubit.cmd('block 1002 element type QUAD4')

cubit.cmd('block 1003 face in surf 5')
cubit.cmd('block 1003 name "pml_xz_elast"')
cubit.cmd('block 1003 attribute count 1')      # number of attributes
cubit.cmd('block 1003 attribute index 1 1')    # material index
cubit.cmd('block 1003 element type QUAD4')

### Remove duplicates and fix problems
cubit.cmd('set dev on') # Set development commands on (to delete)
cubit.cmd('delete face all except face in block all') # This is necessary to remove duplicated faces

### Remove gap in numbering :
cubit.cmd('compress')

### Creating absorbing surfaces block (after the merging!)
cubit.cmd('block 3 edge in curve 15 18 19')
cubit.cmd('block 3 name "abs_right"')
cubit.cmd('block 3 element type BAR2')
cubit.cmd('block 4 edge in curve 12 17')
cubit.cmd('block 4 name "abs_bottom"')
cubit.cmd('block 4 element type BAR2')

### Creating topo surface block (after the merging!)
cubit.cmd('block 5 edge in curve 6 14')
cubit.cmd('block 5 name "topo"')
cubit.cmd('block 5 element type BAR2')

### Create Axis block (must include pmls)
cubit.cmd('block 6 edge in curve 4 9 11')
cubit.cmd('block 6 name "elements_axis"')
cubit.cmd('block 6 element type BAR2')

########################################################################
### Export the mesh under specfem2d format:
profile=mesh() # Store the mesh from Cubit
profile.write("./MESH") # Write it into files (in specfem2d format)

