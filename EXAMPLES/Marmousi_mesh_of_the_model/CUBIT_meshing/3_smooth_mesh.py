#!/usr/bin/env python-cubit
#
# performs additional mesh smoothing to improve mesh quality
# this can take a while...
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
print("# major version number: ",cubit_version_major)

# current work directory
cubit.cmd('pwd')
cubit.cmd('reset')

#checks if new file available
if not os.path.exists("mesh_2_marmousi.cub"):
    print("#")
    print("# file mesh_2_marmousi.cub not found, please check ...")
    print("#")
    cubit.cmd('pause')
cubit.cmd('open "mesh_2_marmousi.cub"')

# surface count
nsurf = cubit.get_surface_count()
print("#")
print("# number of surfaces = ",nsurf)
print("#")

# smoothing
# laplacian
#cubit.cmd('surface all smooth scheme untangle smart laplacian cpu 2')
#cubit.cmd('smooth surface all')
#cubit.cmd('quality surface all condition no. global')
print("#")

# improves mesh by smoothing
cubit.cmd('surface all smooth scheme untangle condition number')
cubit.cmd('smooth surface all')
#cubit.cmd('surface all smooth scheme condition number beta 2.0 cpu 2')
#cubit.cmd('smooth surface all global')
#cubit.cmd('surface all smooth scheme untangle mean ratio cpu 2')
#cubit.cmd('smooth surface all global')

#cubit.cmd('surface all smooth scheme equipotential')
#cubit.cmd('smooth surface all')


# mesh quality summary
cubit.cmd('quality surface all condition no. global')
print("#")

# save
cubit.cmd('save as "mesh_3_marmousi_smoothed.cub" overwrite')

print('#')
print('# done: see file mesh_3_marmousi_smoothed.cub')
print('#')


