#!/usr/bin/env python

###########################################################################
#### Simple axisymmetric 2.5d mesh with 3 pmls
#### Play this script as a python journal file in Cubit13+/Trelis
#### Then you can run cubit2specfem2d.py
#### Author alexis dot bottero at gmail dot com
###########################################################################

import cubit

lx=50000 #5000
lz=50000  #500
nPml=3 # PMLs Thickness (number of elements in pmls)

######################

elementSize = 1000 # Mean element size

#cubit.cmd('cd "/home1/bottero/Desktop/mesh2d"') # Change Cubit working directory
cubit.cmd('reset')

### Create the rectangle :
cubit.cmd('create surface rectangle width '+str(lx)+' height '+str(lz)+' yplane') # Create a lx*lz rectangle
cubit.cmd('volume 1 move x '+str(lx/2.0)+' z '+str(-lz/2.0)) # To get correct coordinates

### Meshing the surface :
cubit.cmd('surface 1 size '+str(elementSize)) # Define the size of an element in surface 1
cubit.cmd('mesh surface 1') # Mesh the surface

### Define blocks and material properties
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')
cubit.cmd(' block 1 face in surface 1')
cubit.cmd('block 1 name "Elastic region"')    # elastic region
cubit.cmd('block 1 attribute count 1')        # number of attributes
cubit.cmd('block 1 attribute index 1 1')      # material index
cubit.cmd('block 1 element type QUAD4')       # element type. Must be QUAD4 for cubit2specfem2d

### Create PMLs quads :
# z pmls :
cubit.cmd('create element extrude edge in curve 1 direction 0 0 -1 distance '+str(3*elementSize)+' layers '+str(nPml))
cubit.cmd('create element extrude edge in curve 3 direction 0 0 1 distance '+str(3*elementSize)+' layers '+str(nPml))
cubit.cmd('block 100 (face all except face in block 1)')
# x pmls :
cubit.cmd('create element extrude edge all with x_coord > '+str(lx-elementSize/100.0)+' direction 1 0 0 distance '+str(3*elementSize)+' layers '+str(nPml))
cubit.cmd('block 101 (face all except face in block 1 100)')
# xz pmls :
cubit.cmd('set duplicate block elements on')
cubit.cmd('block 200 face in block 101 with z_coord > 0')
cubit.cmd('block 201 face in block 101 with z_coord < '+str(-lz))
cubit.cmd('block 102 face in block 200 201')
#cubit.cmd('block 102 group (face in block 101 with z_coord > 0) and (face in block 101 with z_coord < '+str(-lz)) # Does not work??!!

### Create geometry from free mesh (necessary to merge) :
cubit.cmd('create mesh geometry face in block 100 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 101 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 102 feature_angle 135.0')
cubit.cmd('merge all')
cubit.cmd('delete block 100 101 102 200 201')

### Create PMLs blocks
# z pmls :
cubit.cmd('block 1000 face in surf 2 3)')
cubit.cmd('block 1000 name "pml_z_elast"')
cubit.cmd('block 1000 attribute count 1')        # number of attributes
cubit.cmd('block 1000 attribute index 1 1')      # material index
cubit.cmd('block 1000 element type QUAD4')
# x pmls :
cubit.cmd('block 1001 face in surf 4)')
cubit.cmd('block 1001 name "pml_x_elast"')
cubit.cmd('block 1001 attribute count 1')        # number of attributes
cubit.cmd('block 1001 attribute index 1 1')      # material index
cubit.cmd('block 1001 element type QUAD4')
# xz pmls :
cubit.cmd('block 1002 face in surf 5 6)')
cubit.cmd('block 1002 name "pml_xz_elast"')
cubit.cmd('block 1002 attribute count 1')        # number of attributes
cubit.cmd('block 1002 attribute index 1 1')      # material index
cubit.cmd('block 1002 element type QUAD4')

### Creating absorbing surfaces
cubit.cmd('block 2 edge in curve 22 23 24')
cubit.cmd('block 2 name "abs_right"')
cubit.cmd('block 2 element type BAR2')
cubit.cmd('block 3 edge in curve 8 18')
cubit.cmd('block 3 name "abs_bottom"')
cubit.cmd('block 3 element type BAR2')
cubit.cmd('block 4 edge in curve 12 21')
cubit.cmd('block 4 name "abs_top"')
cubit.cmd('block 4 element type BAR2')

### Create Axis (must include pmls)

cubit.cmd('block 5 edge with x_coord < '+str(elementSize/100.0))
cubit.cmd('block 5 name "elements_axis"')
cubit.cmd('block 5 element type BAR2')

### Remove duplicates and fix problems
cubit.cmd('set dev on') # Set development commands on (to be able to delete)
cubit.cmd('delete face all except face in block all') # Delete all duplicate faces
