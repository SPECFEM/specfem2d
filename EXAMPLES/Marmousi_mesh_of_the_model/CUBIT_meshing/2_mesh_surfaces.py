#!/usr/bin/env python-cubit
#
# creates a quad mesh for all surfaces
#
from __future__ import print_function
import os
import sys

# to run this script from command line, python must have its PATH environment set such that it
# includes the path to CUBIT/Trelis (cubit.py).
#
# you can also explicitly set it here e.g. like:
#sys.path.append('/Applications/Coreform-Cubit-2023.11.app/Contents/lib')

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

#######################################################################
# USER PARAMETERS

# small curve length threshold
THRESHOLD_CURVE_LENGTH = 1.0

# mesh minimum element size
# hi-res
#elem_size_min = 1.0
# low-res
elem_size_min = 10.0

# mesh maximum element size
# hi-res
#elem_size_max = 5.0
# low-res
elem_size_max = 25.0   # must be <= 25

#######################################################################

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
if not os.path.exists("mesh_1_surfaces.cub"):
    print("#")
    print("# file mesh_1_surfaces.cub not found, please check ...")
    print("#")
    cubit.cmd('pause')
cubit.cmd('open "mesh_1_surfaces.cub"')


# surface count
nsurf = cubit.get_surface_count()
print("#")
print("# number of surfaces = ",nsurf)
if nsurf != 435:
    print("Error: wrong number of surfaces, must be 435 for this Marmousi model")
    cubit.cmd('pause')
else:
    print("# surface number okay")
print("#")

# note: the imprint and merge commands will create a very small curve (0.15m length)
#       due to an offset at one of the faults.
#
#       this small curve is a challenge for meshing, as it will create distorted elements
#
# small curves
surface_list = []
small_curve_list = []
for i in range(1,nsurf+1):
    #print("surface: ",i)
    curves = cubit.get_relatives("surface",i,"curve")
    has_small_curve = 0
    min_length = 1.e30
    min_curve = 0
    for c in curves:
        length = cubit.get_curve_length(c)
        # stores curve if below threshold
        if length < THRESHOLD_CURVE_LENGTH:
            has_small_curve = 1
            # total minimum length
            if length < min_length:
                min_length = length
                min_curve = c
            if c not in small_curve_list:
                small_curve_list.append(c)
            #print("curve: ",c," is smaller than threshold and has length = ",length)
    # stores surface if it contains a small curve
    if has_small_curve == 1:
        #print("# surface: ",i)
        #print("#  has small curves of minimum length ",min_length)
        #print("#")
        surface_list.append([i,min_curve,min_length])

# output
if len(surface_list) > 0:
    print("# small feature surfaces:")
    for s in surface_list:
        print("#  surface ",s[0]," has curve ",s[1]," with minimum length ",s[2])
    print("#")

# collapses minimum curves
# will not work for current topology...
if len(small_curve_list) > 0:
    print("# small curves:")
    for c in small_curve_list:
        vert = cubit.get_relatives("curve",c,"vertex")
        print("# curve ",c," has vertices ",vert)
    # not working...
    #cubit.cmd('collapse curve '+str(c)+' vertex '+str(vert[0]))
    print("#")

# removing small curves
# will not work...
#cubit.cmd('tweak remove_topology curve all small_curve_size 0.5 backoff_distance 1')

# validation (closed surface?)
cubit.cmd('validate surface all')

# curve count
#n = cubit.get_curve_count()
#print("")
#print("# number of curves = ",n)
#print("")

# skew control will fail...
#cubit.cmd('control skew surface all')

########################################################################
# meshing

# paving
cubit.cmd('surface all scheme pave')

#original
cubit.cmd('surface 186  size '+str(elem_size_min)) # otherwise a meshing error occurs..
cubit.cmd('surface all except 186 size '+str(elem_size_max))


# geometry adaptive
#cubit.silent_cmd('surface all sizing function type skeleton' + \
#                 ' min_size '+ str(elem_size_min) + \
#                 ' max_size ' + str(elem_size_max) + \
#                 ' max_gradient 1.0 min_num_layers_2d 1 min_num_layers_1d 1')


# problematic surface
if 0 == 1:
    # large curve
    cubit.cmd('curve 227 to 229 size ' + str(elem_size))
    cubit.cmd('curve 227 to 229 scheme equal')
    cubit.cmd('mesh curve 227 to 229')
    # progressive sizes
    cubit.cmd('curve 230 scheme bias fine size 0.2 coarse size ' + str(elem_size) + ' start vertex 290')
    cubit.cmd('mesh curve 230')
    # small curve
    cubit.cmd('curve 231 size 1')
    cubit.cmd('curve 231 scheme equal')
    cubit.cmd('mesh curve 231')
    # progressive sizes
    cubit.cmd('curve 581 scheme bias fine size 0.2 coarse size '  + str(elem_size) + ' start vertex 291')
    cubit.cmd('mesh curve 581')

    # problematic surfaces
    cubit.cmd('mesh surface 149')
    cubit.cmd('quality surface 149 condition no. global')
#else:
    #cubit.cmd('mesh surface 149')
    #cubit.cmd('quality surface 149 condition no. global')

# meshes small curves
if 0 == 1:
    for c in small_curve_list:
        length = cubit.get_curve_length(c)
        print("# meshing small curve: ",c,"length = ",length," - mesh size min/max = ",elem_size_min,elem_size_max)
        if length < elem_size_min:
            # equal sizes
            cubit.cmd('curve ' + str(c) + ' size ' + str(elem_size_min))
            cubit.cmd('curve ' + str(c) + ' scheme equal')
        else:
            # equal sizes
            cubit.cmd('curve ' + str(c) + ' size ' + str(elem_size_max))
            cubit.cmd('curve ' + str(c) + ' scheme equal')

            # progressive sizes
            #vert = cubit.get_relatives("curve",c,"vertex")
            #cubit.cmd('curve ' + str(c) + ' scheme bias' + \
            #          ' fine size ' + str(elem_size_min) + \
            #          ' coarse size ' + str(elem_size_max) + ' start vertex ' + str(vert[0]))
        cubit.cmd('mesh curve ' + str(c))

def mesh_surface(i,elem_size_min,elem_size_max):
    """
    meshes a surface
    """
    print("# meshing surface: ",i)
    #print("# surface area = ",cubit.get_surface_area(i))
    if cubit.is_meshed("surface",i):
        print("# already done")
        print("#")
        return

    elem_size1 = elem_size_min
    elem_size2 = elem_size_max
    cubit.cmd('mesh surface '+str(i))
    # retry if there was an error
    itry = 0
    while cubit.get_error_count() > 0 and itry < 8:
        # retry with smaller size
        #cubit.cmd('delete mesh surface '+str(i)+' propagate')
        cubit.cmd('delete mesh surface '+str(i))
        cubit.cmd('reset errors')
        itry += 1
        elem_size2 = elem_size2 / 2.0
        print("# ",itry,"retrying with maximum element size = ",elem_size2)
        cubit.cmd('surface '+str(i)+' size '+str(elem_size2))

        # geometry adaptive
        #if elem_size1 > elem_size2:
        #    elem_size1 = elem_size2
        #    cubit.cmd('surface '+str(i)+' size '+str(elem_size1))
        #else:
        #    cubit.cmd('surface '+str(i)+' sizing function type skeleton min_size '+str(elem_size1) + ' max_size '+str(elem_size2) + ' max_gradient 2.0 min_num_layers_2d 1 min_num_layers_1d 1')
        cubit.cmd('mesh surface '+str(i))

    # check mesh
    if 0 == 1:
        #cubit.cmd('group "myfaces" add face in surface ' + str(i))
        #group1 = cubit.get_id_from_name("myfaces")
        #quads = cubit.get_group_quads(group1)
        #print("# mesh quads ",quads)
        #for j in quads:
        #    print("# mesh quad ",j)
        #not working...
        #print("# mesh condition number quality = ",cubit.get_quality_value("quad",group1,"condition"))
        #print("# mesh skew quality             = ",cubit.get_quality_value("quad",group1,"skew"))
        #print("# mesh shape quality            = ",cubit.get_quality_value("quad",group1,"shape"))
        #print("")
        #cubit.cmd('del group myfaces')
        # counts quads with low shape metric (between 0 and 0.1)
        shape_min = 0.0
        shape_max = 0.01
        #
        cubit.silent_cmd('group "poor" add quality surface ' + str(i) + \
                         ' shape low ' + str(shape_min) + ' high ' + str(shape_max))
        group2 = cubit.get_id_from_name("poor")
        quads = cubit.get_group_quads(group2)
        cubit.delete_group(group2)
        print("# number of quads with poor shape = ",len(quads))
        print("#")
        # retry
        itry = 0
        while len(quads) > 0 and itry < 8:
            # retry with smaller size
            #cubit.cmd('delete mesh surface '+str(i)+' propagate')
            cubit.cmd('delete mesh surface '+str(i))
            itry += 1
            elem_size2 = elem_size2 / 2.0
            print("# ",itry,"retrying with maximum element size = ",elem_size2)
            if elem_size1 > elem_size2:
                elem_size1 = elem_size2
                cubit.cmd('surface '+str(i)+' size '+str(elem_size1))
            else:
                cubit.cmd('surface '+str(i)+' sizing function type skeleton min_size '+str(elem_size1) + ' max_size '+str(elem_size2) + ' max_gradient 2.0 min_num_layers_2d 1 min_num_layers_1d 1')
            cubit.cmd('mesh surface '+str(i))
            # check mesh
            cubit.silent_cmd('group "poor" add quality surface ' + str(i) + \
                             ' shape low ' + str(shape_min) + ' high ' + str(shape_max))
            group2 = cubit.get_id_from_name("poor")
            quads = cubit.get_group_quads(group2)
            cubit.delete_group(group2)
            print("# number of quads with poor shape = ",len(quads))
            print("")

    # safety check
    if cubit.get_error_count() > 0:
        print("# mesh: generated error")
        print("#")
        cubit.cmd('stop')
        sys.exit(1)


cubit.cmd('reset errors')

# small curve surface meshing
#for s in surface_list:
#    i = s[0]
#    print('# surface with small curve ',i)
#    mesh_surface(i,elem_size_min,elem_size_max)

# gets list
surface_area = []
for i in range(1,nsurf+1):
    area = cubit.get_surface_area(i)
    surface_area.append(area)

# gets index of sorted areas
surface_sorted = [i[0] for i in sorted(enumerate(surface_area), key=lambda x:x[1])]
icount = 0
for i in surface_sorted:
    icount += 1
    # surface id = list id + 1
    id = i + 1
    print("# ",icount," surface ",id,"area = ",surface_area[i])
    #mesh_surface(id,elem_size_min,elem_size_max)

cubit.cmd('mesh surface all')

# checks if all surfaces are meshed
for i in range(1,nsurf+1):
    if not cubit.is_meshed("surface",i):
        print("# mesh failed for surface ",i)
        print("#")
        sys.exit(1)
    else:
        print("# meshed surface ",i)

# error check
if cubit.get_error_count() > 0:
    print("#")
    print("# mesh: generated errors, aborting... ")
    print("#")
    cubit.cmd('stop')
    sys.exit(1)
else:
    print("#")
    print("# mesh: all done")
    print("#")


# warning info - not supported anymore...
#if cubit.get_warning_count() > 0:
#    print("")
#    print("meshing generated warnings: ",n)
#    print("")

# validate
cubit.cmd('validate surface all mesh')

#c_min = cubit.get_quality_value("QUAD",149,"condition no")
#print("")
#print("# mesh: minimum surface mesh quality value = ",c_min)
#print("")

print("#")
print("# mesh: total number of QUAD elements = ",cubit.get_quad_count())
print("#")

# mesh quality summary
cubit.cmd('quality surface all condition no. global')
print('#')

# save
cubit.cmd('save as "mesh_2_marmousi.cub" overwrite')

print("#")
print("# done: see file mesh_2_marmousi.cub")
print("#")

