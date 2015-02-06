#!/usr/bin/env python

###########################################################################
#### Simple axisymmetric 2.5d mesh with pmls
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
cubit.cmd('mesh surface 1')

### End of meshing

### Define material properties
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')
cubit.cmd(' block 1 face in surface 1')
cubit.cmd('block 1 name "Material"')         # elastic region
cubit.cmd('block 1 attribute count 6')        # number of attributes
cubit.cmd('block 1 attribute index 1 1')    # material index
cubit.cmd('block 1 attribute index 2 2500 ')  # vp
cubit.cmd('block 1 attribute index 3 1000 ')  # vs
cubit.cmd('block 1 attribute index 4 2200 ')  # rho
cubit.cmd('block 1 attribute index 5 0 ')     # Q_flag
cubit.cmd('block 1 attribute index 6 0 ')     # anisotropy_flag
cubit.cmd('block 1 element type QUAD4')       # element type. Must be QUAD4 for cubit2specfem2d

### Creating absorbing surfaces
cubit.cmd('block 2 edge in curve 4') 
cubit.cmd('block 2 name "abs_right"')
cubit.cmd('block 2 element type BAR2')
cubit.cmd('block 3 edge in curve 1') 
cubit.cmd('block 3 name "abs_bottom"')
cubit.cmd('block 3 element type BAR2')
cubit.cmd('block 4 edge in curve 3') 
cubit.cmd('block 4 name "abs_top"')
cubit.cmd('block 4 element type BAR2')

### Create PMLs :
# z pmls :
cubit.cmd('create element extrude edge in block 3 direction 0 0 -1 distance '+str(3*elementSize)+' layers '+str(nPml))
cubit.cmd('create element extrude edge in block 4 direction 0 0 1 distance '+str(3*elementSize)+' layers '+str(nPml))
cubit.cmd('block 100 (face all except face in block 1)')
cubit.cmd('block 100 name "pml_z"')
cubit.cmd('block 100 attribute count 6')        # number of attributes
cubit.cmd('block 100 attribute index 1 1')    # material index
cubit.cmd('block 100 attribute index 2 2500 ')  # vp
cubit.cmd('block 100 attribute index 3 1000 ')     # vs
cubit.cmd('block 100 attribute index 4 2200 ')  # rho
cubit.cmd('block 100 attribute index 5 0 ')     # Q_flag
cubit.cmd('block 100 attribute index 6 0 ')     # anisotropy_flag
cubit.cmd('block 100 element type QUAD4')

# x pmls :
cubit.cmd('create element extrude edge all with x_coord > '+str(lx-elementSize/100.0)+' direction 1 0 0 distance '+str(3*elementSize)+' layers '+str(nPml))
cubit.cmd('block 101 (face all except face in block 1 100)')
cubit.cmd('block 101 name "pml_x"')
cubit.cmd('block 101 attribute count 6')        # number of attributes
cubit.cmd('block 101 attribute index 1 1')    # material index
cubit.cmd('block 101 attribute index 2 2500 ')  # vp
cubit.cmd('block 101 attribute index 3 1000 ')     # vs
cubit.cmd('block 101 attribute index 4 2200 ')  # rho
cubit.cmd('block 101 attribute index 5 0 ')     # Q_flag
cubit.cmd('block 101 attribute index 6 0 ')     # anisotropy_flag
cubit.cmd('block 101 element type QUAD4')

# xz pmls :
cubit.cmd('set duplicate block elements on')
cubit.cmd('block 300 face in block 101 with z_coord > 0')
cubit.cmd('block 301 face in block 101 with z_coord < '+str(-lz)) 
cubit.cmd('block 102 face in block 300 301') 
cubit.cmd('delete block 300 301')
#cubit.cmd('block 102 group (face in block 101 with z_coord > 0) and (face in block 101 with z_coord < '+str(-lz)) # Does not work??!!
cubit.cmd('block 102 name "pml_xz"')
cubit.cmd('block 102 attribute count 6')        # number of attributes
cubit.cmd('block 102 attribute index 1 1')    # material index
cubit.cmd('block 102 attribute index 2 2500 ')  # vp
cubit.cmd('block 102 attribute index 3 1000 ')     # vs
cubit.cmd('block 102 attribute index 4 2200 ')  # rho
cubit.cmd('block 102 attribute index 5 0 ')     # Q_flag
cubit.cmd('block 102 attribute index 6 0 ')     # anisotropy_flag
cubit.cmd('block 102 element type QUAD4')

### Create geometry from free mesh (necessary to merge) :
cubit.cmd('create mesh geometry Face in block 100 feature_angle 135.0')
cubit.cmd('create mesh geometry Face in block 101 feature_angle 135.0')
cubit.cmd('create mesh geometry Face in block 102 feature_angle 135.0')
cubit.cmd('merge all')

### Create Axis (must include pmls)

cubit.cmd('block 6 edge with x_coord < '+str(elementSizeSolid/100.0))
cubit.cmd('block 6 name "elements_axis"')

### Remove duplicates and fix problems
cubit.cmd('block 101 remove face in block 102')
cubit.cmd('set dev on') # Set development commands on (to delete)
cubit.cmd('delete face all except face in block all') # Sometimes this could be usefull!
#cubit.cmd('block 3 element type BAR2') # with some examples this is necessary
