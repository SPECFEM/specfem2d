#!/usr/bin/env python
###########################################################################
#### This example create a square, mesh it and add PMLs
#### It should be run as a Python journal file in cubit
#### Then just play cubit2specfem2d.py as a journal file and it will create
#### specfem files
###########################################################################

# to run this script from command line, python must have its PATH environment set such that it
# includes the path to CUBIT/Trelis (cubit.py).
#
# you can also explicitly set it here e.g. like:
#sys.path.append('/opt/Trelis-15.0/bin/')
import sys
try:
    import cubit
except ImportError:
    print "Error: Importing cubit as python module failed"
    print "could not import cubit, please check your PYTHONPATH settings..."
    print ""
    print "current path: "
    print sys.path
    print ""
    print "try to include path to directory which includes file cubit.py, e.g. /opt/Trelis-15.0/bin/"
    print ""
    sys.exit("Import cubit failed")

#print sys.path

cubit.init([""])

elementSize=0.5
nPml=3

# Create the geometry

cubit.cmd('create surface rectangle width 10 height 10 yplane')
cubit.cmd('volume 1 move x 5 z -5') # To get correct coordinates
cubit.cmd('color volume 1 grey')

# Imprinting # This is not necessary here as we have just one volume but it has to be done for more than one!

cubit.cmd('imprint all')

# Meshing #

cubit.cmd('surface 1 size '+str(elementSize))
cubit.cmd('mesh surf all')

# Define material properties #
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')
cubit.cmd('block 1 face in surf 1')           # Create block number 1 containing all faces in surface 1
cubit.cmd('block 1 name "Granit"')            # The name of the block is...
cubit.cmd('block 1 attribute count 1')        # It will contain ... attributes
cubit.cmd('block 1 attribute index 1 1')      # First attribute is material index
cubit.cmd('block 1 element type QUAD4')       # Element type. Must be QUAD4 for cubit2specfem2d

# Create PMLs layers #
cubit.cmd('create element extrude edge in curve 1 direction 0 0 -1 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('create element extrude edge in curve 3 direction 0 0 1 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 100 (face all except face in block 1)')

cubit.cmd('create element extrude edge all with x_coord > '+str(10-elementSize/100)+' direction 1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('create element extrude edge all with x_coord < '+str(elementSize/100)+' direction -1 0 0 distance '+str(nPml*elementSize)+' layers '+str(nPml))
cubit.cmd('block 101 (face all except face in block 1 100)')

cubit.cmd('set duplicate block elements on') # To be able to create blocks containing elements already contained in others
cubit.cmd('block 200 face in block 101 with z_coord > 0')
cubit.cmd('block 201 face in block 101 with z_coord < '+str(-10))
cubit.cmd('block 102 face in block 200 201')

# Create geometry (surfaces, curves... from free elements just created) #
cubit.cmd('create mesh geometry face in block 100 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 101 feature_angle 135.0')
cubit.cmd('create mesh geometry face in block 102 feature_angle 135.0')

cubit.cmd('merge all')

# Delete temporqry blocks #
cubit.cmd('delete block 100 101 102 200 201')

# Create definitive blocks (that will be read by cubit2specfem2d.py) #
cubit.cmd('block 1000 face in surf 2 3')
cubit.cmd('block 1000 name "pml_z_elast"')
cubit.cmd('block 1000 attribute count 1')
cubit.cmd('block 1000 attribute index 1 1') # 1 if you speak about elastic PML, 2 if it is acoustic
cubit.cmd('block 1000 element type QUAD4')

cubit.cmd('block 1001 face in surf 4 7')
cubit.cmd('block 1001 name "pml_x_elast"')
cubit.cmd('block 1001 attribute count 1')
cubit.cmd('block 1001 attribute index 1 1')
cubit.cmd('block 1001 element type QUAD4')

cubit.cmd('block 1002 face in surf 5 6 8 9')
cubit.cmd('block 1002 name "pml_xz_elast"')
cubit.cmd('block 1002 attribute count 1')
cubit.cmd('block 1002 attribute index 1 1')
cubit.cmd('block 1002 element type QUAD4')

cubit.cmd('block 3 edge in curve 30 8 18')
cubit.cmd('block 3 name "abs_bottom"')
cubit.cmd('block 3 element type BAR2')

cubit.cmd('block 4 edge in curve 34 35 36')
cubit.cmd('block 4 name "abs_left"')
cubit.cmd('block 4 element type BAR2')

cubit.cmd('block 5 edge in curve 12 21 33')
cubit.cmd('block 5 name "abs_top"')
cubit.cmd('block 5 element type BAR2')

cubit.cmd('block 6 edge in curve 22 23 24')
cubit.cmd('block 6 name "abs_right"')
cubit.cmd('block 6 element type BAR2')

cubit.cmd('set dev on') # This is needed for next command:
cubit.cmd('delete face all except face in block all')


