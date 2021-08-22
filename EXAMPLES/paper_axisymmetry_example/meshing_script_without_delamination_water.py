#!/usr/bin/env python
###########################################################################
#### This example create the model, mesh it and add PMLs
#### It should be run as a Python journal file in cubit
####
####      ________________________________________________________________
####     |    _________________________________________________________   |
####     |   |  ^                                                     |   |
####     | C |  | lwUp                Water                           |   |
####     | P |__v_____________________________________________________|   |
####     | M |  ^                                                     |   |
####     | L |  | lzUp                                                |   |
####     |   |  v                       lGlue                         | C |
####     |   |----------------------------^---------------------------| P |
####     |   X (0,0)                      |                           | M |
####     |   |----------------------------v---------------------------| L |
####     | C |  ^                                                     |   |
####     | P |  | lzDown                                              |   |
####     | M |__v_____________________________________________________|   |
####     | L |  ^                                                     |   |
####     |   |  |                     Water                           |   |
####     |   |__v_____________________________________________________|   |
####     |________________________________________________________________|
####         <-------------------------------------------------------->
####                                   lx
####
####      Point (0,0) is in the middle of the glue at the left
####
###########################################################################
from __future__ import print_function
import cubit
import sys

sys.path.append('/home/bottero/specfem2d/utils/cubit2specfem2d')
from cubit2specfem2d import mtools,block_tools,mesh_tools,mesh

# All the distances are given in meters. See the header for a draft of the model
#f0 = 250e3
#fmax = f0 * 2.5 # maximal frequency used in specfem to calculte the resolution (true for a Ricker source) fmax = f0*2,5
#lambda_min = 1036/fmax # smaller wavelength due to the small value of shear wave speed in epoxy ... thus normally the model is enough sampled
e#lementSize = lambda_min # Average spectral elements size
elementSize = 1.326e-3
lx = 370.0e-3 		# TODO For the moment, I use the same dimension as the transmision through the homogeneous plate
lzUp = 3.0e-3 		# the value lx is used in runSpecFem.py, more precisely in the script CreateReceivers to generate a receiver line over the whole domain
lzGlue = 0.25e-3 	# it will be a good thing to give these geometric parameters as argument into the make_mesh.sh script !!!!
lzDown = 3.0e-3
lwUp  =  110.0e-3   # same dimension as the transmision through the homogeneous plate
lwDown = 280.0e-3   # same dimension as the transmision through the homogeneous plate

# I think that 3 is the best option below:
nPml = 3  # Width of PML layer in term of elements

#########################################################
############ NOTHING HAS TO BE CHANGED BELOW ############
#########################################################

# To debug
if False:
    elementSize = 0.5e-3
    lx = 10.0e-3
    lzUp = 2.0e-3
    lzGlue = 0.25e-3
    lzDown = 2.0e-3
    lwUp  = 3.0e-3
    lwDown = 2.0e-3
    nPml = 3

### Create the geometry ###

cubit.cmd('merge tolerance '+str(elementSize/1000.0)) # Distance below which we consider that two points are the same

# Create rectangles
cubit.cmd('create surface rectangle width '+str(lx)+' height '+str(lzGlue)+' yplane')
cubit.cmd('create surface rectangle width '+str(lx)+' height '+str(lzUp)+' yplane')
cubit.cmd('volume 2 move x 0.0 y 0.0 z '+str(lzUp/2.0+lzGlue/2.0))
cubit.cmd('create surface rectangle width '+str(lx)+' height '+str(lzDown)+' yplane')
cubit.cmd('volume 3 move x 0.0 y 0.0 z -'+str(lzDown/2.0+lzGlue/2.0))
cubit.cmd('volume all move x '+str(lx/2.0)+' y 0 z 0') # To get correct x coordinates

# Merge what must be merged
cubit.cmd('merge all')

# Compress numbering
cubit.cmd('compress')

cubit.cmd('create surface rectangle width '+str(lx)+' height '+str(lzUp+lzGlue+lzDown+lwUp+lwDown)+' yplane')
cubit.cmd('move vertex 11 location 0 0 '+str(lzUp+lwUp+lzGlue/2.0)+' include_merged')
cubit.cmd('subtract body 1 2 3 from body 4 keep ')
cubit.cmd('delete body 4')

# Merge what must be merged
cubit.cmd('merge all')

# Compress numbering
cubit.cmd('compress')

# Cosmetic stuff
cubit.cmd('color body 1 grey')
cubit.cmd('color body 2 khaki')
cubit.cmd('color body 3 khaki')
cubit.cmd('color body 4 lightskyblue')

### Meshing ###

cubit.cmd('surface all size '+str(elementSize))
cubit.cmd('mesh surf all')

### Define material properties ###

cubit.cmd('####################### DEFINE MATERIAL PROPERTIES #######################')
cubit.cmd('block 1 face in surf 2')           # Create block number 1 containing all faces in surface 1
cubit.cmd('block 1 name "FiberUp"')           # The name of the block is...
cubit.cmd('block 1 attribute count 1')        # It will contain ... attributes
cubit.cmd('block 1 attribute index 1 1')      # First attribute is material index
cubit.cmd('block 1 element type QUAD4')       # Element type. Must be QUAD4 for cubit2specfem2d

cubit.cmd('block 2 face in surf 3')           # Create block number 1 containing all faces in surface 2
cubit.cmd('block 2 name "FiberDown"')         # The name of the block is...
cubit.cmd('block 2 attribute count 1')        # It will contain ... attributes
cubit.cmd('block 2 attribute index 1 2')      # First attribute is material index
cubit.cmd('block 2 element type QUAD4')       # Element type. Must be QUAD4 for cubit2specfem2d

cubit.cmd('block 3 face in surf 1')           # Create block number 1 containing all faces in surface 3
cubit.cmd('block 3 name "Resin"')             # The name of the block is...
cubit.cmd('block 3 attribute count 1')        # It will contain ... attributes
cubit.cmd('block 3 attribute index 1 3')      # First attribute is material index
cubit.cmd('block 3 element type QUAD4')       # Element type. Must be QUAD4 for cubit2specfem2d

cubit.cmd('block 4 face in surf 4 5')         # Create block number 1 containing all faces in surface 4 5
cubit.cmd('block 4 name "Water"')             # The name of the block is...
cubit.cmd('block 4 attribute count 1')        # It will contain ... attributes
cubit.cmd('block 4 attribute index 1 4')      # First attribute is material index
cubit.cmd('block 4 element type QUAD4')       # Element type. Must be QUAD4 for cubit2specfem2d

cubit.cmd('set duplicate block elements on')
cubit.cmd('block 999 (face all)')

### Create PMLs layers ###

cubit.cmd('create element extrude edge in curve 7 direction 1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 100 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in curve 4 direction 1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 101 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in curve 10 direction 1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 102 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in curve 5 direction -1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 103 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in curve 2 direction -1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 104 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in curve 9 direction -1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 105 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in curve 11 direction 1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 106 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in curve 16 direction 1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 107 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in curve 13 direction -1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 108 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in curve 14 direction -1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 109 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in curve 15 direction 0 0 1 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 110 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in curve 12 direction 0 0 -1 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 111 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in face all in block 107 with z_coord > '+str(lzGlue/2.0+lzUp+lwUp-elementSize/100.0)+' direction 0 0 1 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 112 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in face all in block 109 with z_coord > '+str(lzGlue/2.0+lzUp+lwUp-elementSize/100.0)+' direction 0 0 1 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 113 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in face all in block 106 with z_coord < '+str(-lzGlue/2.0-lzDown-lwDown+elementSize/100.0)+' direction 0 0 -1 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 114 (face all except face in block 999)')
cubit.cmd('block 999 (face all)')

cubit.cmd('create element extrude edge in face all in block 108 with z_coord < '+str(-lzGlue/2.0-lzDown-lwDown+elementSize/100.0)+' direction 0 0 -1 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 115 (face all except face in block 999)')

cubit.cmd('create mesh geometry face in block 100 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 101 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 102 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 103 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 104 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 105 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 106 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 107 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 108 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 109 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 110 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 111 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 112 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 113 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 114 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 115 feature_angle 135.0')

cubit.cmd('merge all')

cubit.cmd('compress')

# Delete temporary blocks
cubit.cmd('delete block 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 999')

# Create definitive pml blocks (that will be read by cubit2specfem2d.py)

cubit.cmd('block 1000 face in surf 6 9')
cubit.cmd('block 1000 name "pml_x_elast_1"')
cubit.cmd('block 1000 attribute count 1')
cubit.cmd('block 1000 attribute index 1 1')
cubit.cmd('block 1000 element type QUAD4')

cubit.cmd('block 1001 face in surf 8 11')
cubit.cmd('block 1001 name "pml_x_elast_2"')
cubit.cmd('block 1001 attribute count 1')
cubit.cmd('block 1001 attribute index 1 2')
cubit.cmd('block 1001 element type QUAD4')

cubit.cmd('block 1002 face in surf 7 10')
cubit.cmd('block 1002 name "pml_x_elast_3"')
cubit.cmd('block 1002 attribute count 1')
cubit.cmd('block 1002 attribute index 1 3')
cubit.cmd('block 1002 element type QUAD4')

cubit.cmd('block 1003 face in surf 12 13 14 15')
cubit.cmd('block 1003 name "pml_x_acoust"')
cubit.cmd('block 1003 attribute count 1')
cubit.cmd('block 1003 attribute index 1 4')
cubit.cmd('block 1003 element type QUAD4')

cubit.cmd('block 1004 face in surf 16 17')
cubit.cmd('block 1004 name "pml_z_acoust"')
cubit.cmd('block 1004 attribute count 1')
cubit.cmd('block 1004 attribute index 1 4')
cubit.cmd('block 1004 element type QUAD4')

cubit.cmd('block 1005 face in surf 18 19 20 21')
cubit.cmd('block 1005 name "pml_xz_acoust"')
cubit.cmd('block 1005 attribute count 1')
cubit.cmd('block 1005 attribute index 1 4')
cubit.cmd('block 1005 element type QUAD4')

cubit.cmd('block 10 edge in curve 45 34 19 21 23 32 49')
cubit.cmd('block 10 name "abs_right"')
cubit.cmd('block 10 element type BAR2')

cubit.cmd('block 11 edge in curve 52 44 50')
cubit.cmd('block 11 name "abs_bottom"')
cubit.cmd('block 11 element type BAR2')

cubit.cmd('block 12 edge in curve 47 38 26 28 30 36 51')
cubit.cmd('block 12 name "abs_left"')
cubit.cmd('block 12 element type BAR2')

cubit.cmd('block 13 edge in curve 48 41 46')
cubit.cmd('block 13 name "abs_top"')
cubit.cmd('block 13 element type BAR2')

cubit.cmd('set dev on') # This is needed for next command:
cubit.cmd('delete face all except face in block all')

profile=mesh() # Store the mesh from Cubit
profile.write("./MESH/") # Write it into files (in specfem2d format). profile.write(/path/to/directory)


