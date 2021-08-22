# boundary_definition2D.py
#
# this file is based on the 3D-version boundary_defnition.py,
# part of GEOCUBIT (see SPECFEM3D/CUBIT_GEOCUBIT/ folder)
from __future__ import print_function
try:
    import start as start
    cubit                   = start.start_cubit()
except:
    try:
        import cubit
    except:
        print('error importing cubit, check if cubit is installed')
        pass


##-----------------------------------------------------------------------------
##
## defines blocks for edges along model boundaries
##
##-----------------------------------------------------------------------------



def define_bc_edges():
    """
    define the absorbing edges for a layered topological box where boundary are edges parallel to the axis.
    it returns absorbing_surf,absorbing_surf_xmin,absorbing_surf_xmax,absorbing_surf_bottom,topo_surf
    where
    absorbing_surf is the list of all the absorbing boundary edges
    absorbing_surf_xmin is the list of the absorbing boundary edges that correnspond to x=xmin
    ...
    absorbing_surf_bottom is the list of the absorbing boundary edges that correspond to z=zmin
    """
    #from utilities import get_v_h_list
    #
    list_vol = cubit.parse_cubit_list("volume","all",1)
    print("#define_edge: volume list length = ",len(list_vol))

    xmin_box = cubit.get_total_bounding_box("volume",list_vol)[0]
    xmax_box = cubit.get_total_bounding_box("volume",list_vol)[1]

    ymin_box = cubit.get_total_bounding_box("volume",list_vol)[3]
    ymax_box = cubit.get_total_bounding_box("volume",list_vol)[4]

    zmin_box = cubit.get_total_bounding_box("volume",list_vol)[6] #it is the z_min of the box ... box= xmin,xmax,d,ymin,ymax,d,zmin...
    zmax_box = cubit.get_total_bounding_box("volume",list_vol)[7]

    print("#define_edge: bounding box xmin/xmax = ",xmin_box,xmax_box)
    print("#define_edge: bounding box ymin/ymax = ",ymin_box,ymax_box)
    print("#define_edge: bounding box zmin/zmax = ",zmin_box,zmax_box)

    # cubit2specfem2d.py needs one block per border (abs_bottom, abs_right, abs_left, abs_top, topo, axis)
    #
    # 2D model should be defined in with top at maximum, bottom at minimum (e.g. z pointing up in positive direction, not depth)
    if abs(xmax_box-xmin_box) < 0.001:
        print("")
        print("# model in YZ plane")
        print("")
        border_t = zmax_box - 0.1 # top minus margin
        border_b = zmin_box + 0.1 # bottom
        border_l = ymin_box + 0.1 # left
        border_r = ymax_box - 0.1 # right

        cubit.cmd('block 1001 edge in surf all with z_coord > ' + str(border_t))
        cubit.cmd('block 1001 name "topo"')

        cubit.cmd('block 1002 edge in surf all with z_coord < ' + str(border_b))
        cubit.cmd('block 1002 name "abs_bottom"')

        cubit.cmd('block 1003 edge in surf all with y_coord < ' + str(border_l))
        cubit.cmd('block 1003 name "abs_left"')

        cubit.cmd('block 1004 edge in surf all with y_coord > ' + str(border_r))
        cubit.cmd('block 1004 name "abs_right"')

    if abs(ymax_box-ymin_box) < 0.001:
        print("")
        print("# model in XZ plane")
        print("")
        border_t = zmax_box - 0.1 # top minus margin
        border_b = zmin_box + 0.1 # bottom
        border_l = xmin_box + 0.1 # left
        border_r = xmax_box - 0.1 # right

        cubit.cmd('block 1001 edge in surf all with z_coord > ' + str(border_t))
        cubit.cmd('block 1001 name "topo"')

        cubit.cmd('block 1002 edge in surf all with z_coord < ' + str(border_b))
        cubit.cmd('block 1002 name "abs_bottom"')

        cubit.cmd('block 1003 edge in surf all with x_coord < ' + str(border_l))
        cubit.cmd('block 1003 name "abs_left"')

        cubit.cmd('block 1004 edge in surf all with x_coord > ' + str(border_r))
        cubit.cmd('block 1004 name "abs_right"')

    if abs(zmax_box-zmin_box) < 0.001:
        print("")
        print("# model in XY plane")
        print("")
        border_t = ymax_box - 0.1 # top minus margin
        border_b = ymin_box + 0.1 # bottom
        border_l = xmin_box + 0.1 # left
        border_r = xmax_box - 0.1 # right

        cubit.cmd('block 1001 edge in surf all with y_coord > ' + str(border_t))
        cubit.cmd('block 1001 name "topo"')

        cubit.cmd('block 1002 edge in surf all with y_coord < ' + str(border_b))
        cubit.cmd('block 1002 name "abs_bottom"')

        cubit.cmd('block 1003 edge in surf all with x_coord < ' + str(border_l))
        cubit.cmd('block 1003 name "abs_left"')

        cubit.cmd('block 1004 edge in surf all with x_coord > ' + str(border_r))
        cubit.cmd('block 1004 name "abs_right"')

    # outputs block infos
    list_ids = cubit.get_block_id_list()
    print("#blocks: ",list_ids)
    for id in list_ids:
        print("#block id: ",id)
        name = cubit.get_exodus_entity_name("block", id)
        print("#block name: ",name)
        num_a = cubit.get_block_attribute_count(id)
        print("#block number of attributes: ",num_a)
        #attr = cubit.get_block_attribute_value(id,1)
        edges = cubit.get_block_edges(id)
        print("#    number of edges = ",len(edges))
        print("")

    return
