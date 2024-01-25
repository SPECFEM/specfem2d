#!/usr/bin/env python
import os
import sys

## CUBIT/Trelis
import cubit
try:
    #cubit.init([""])
    cubit.init(["-noecho","-nojournal"])
except:
    pass

version = cubit.get_version()
version_major = int(version.split(".")[0])
version_minor = int(version.split(".")[1])
print "cubit version: ",version

os.system('mkdir -p MESH')

cubit.cmd('set echo on')

## This is boundary_definition.py 
# extracts the bounding edges and defines them into blocks
import boundary_definition
reload(boundary_definition)

boundary_definition.define_bc_edges()

## Define material properties for model
cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')
#cubit.cmd('set duplicate block elements on') # allow for multiple occurrences of entities in different blocks 
cubit.cmd('block 1 face in volume all')
cubit.cmd('block 1 element type QUAD4')

cubit.cmd('#### DEFINE MATERIAL PROPERTIES #######################')
cubit.cmd('block 1 name "elastic 1" ')        # elastic material region
cubit.cmd('block 1 attribute count 6')
cubit.cmd('block 1 attribute index 1 1  ')      # volume 1
cubit.cmd('block 1 attribute index 2 2800 ')   # vp
cubit.cmd('block 1 attribute index 3 1500 ')   # vs
cubit.cmd('block 1 attribute index 4 2300 ')   # rho
cubit.cmd('block 1 attribute index 5 9000.0 ')       # Q_mu
cubit.cmd('block 1 attribute index 6 0 ')     # anisotropy_flag

cubit.cmd('save as "top.cub" overwrite')

## export mesh
import cubit2specfem2d
reload(cubit2specfem2d)
cubit2specfem2d.export2SPECFEM2D('MESH/')

cubit.cmd('set echo off')

