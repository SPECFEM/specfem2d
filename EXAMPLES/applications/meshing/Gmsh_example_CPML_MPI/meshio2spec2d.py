# convert mesh to specfem format
import numpy as np
import numba
import time
import os

@numba.jit(nopython=True)
def get_cpml_cells_except_damping(
                PML_X_elms_orig,
                PML_Y_elms_orig,
                PML_XY_elms_orig,
                cells_quad_total,
                cells_line_total,
                cells_line_inner_boundary,
                elm_id_offset):
    PML_X_elms = []
    PML_Y_elms = []
    PML_XY_elms = []

    # create a node list of the inner boundary
    inner_boundary_nodes = []
    for _ielm in cells_line_inner_boundary:
        inner_boundary_nodes.extend(cells_line_total[_ielm])

    for _ielm in PML_X_elms_orig:
        # offcet
        ielm = int(_ielm - elm_id_offset)
        # check if the element is not in the inner boundary
        this_nodes_elm = cells_quad_total[ielm]

        # check if the nodes are not in the inner boundary
        if_this_element_touches_inner_boundary = False
        for this_node in this_nodes_elm:
            if this_node in inner_boundary_nodes:
                if_this_element_touches_inner_boundary = True
                break

        if not if_this_element_touches_inner_boundary:
            PML_X_elms.append(_ielm)

    for _ielm in PML_Y_elms_orig:
        # offcet
        ielm = int(_ielm - elm_id_offset)
        # check if the element is not in the inner boundary
        this_nodes_elm = cells_quad_total[ielm]

        # check if the nodes are not in the inner boundary
        if_this_element_touches_inner_boundary = False
        for this_node in this_nodes_elm:
            if this_node in inner_boundary_nodes:
                if_this_element_touches_inner_boundary = True
                break

        if not if_this_element_touches_inner_boundary:
            PML_Y_elms.append(_ielm)

    for _ielm in PML_XY_elms_orig:
        # offcet
        ielm = int(_ielm - elm_id_offset)
        # check if the element is not in the inner boundary
        this_nodes_elm = cells_quad_total[ielm]

        # check if the nodes are not in the inner boundary
        if_this_element_touches_inner_boundary = False
        for this_node in this_nodes_elm:
            if this_node in inner_boundary_nodes:
                if_this_element_touches_inner_boundary = True
                break

        if not if_this_element_touches_inner_boundary:
            PML_XY_elms.append(_ielm)

    # check the changes in the number of elements
    diff_x = len(PML_X_elms_orig) - len(PML_X_elms)
    diff_y = len(PML_Y_elms_orig) - len(PML_Y_elms)
    diff_xy = len(PML_XY_elms_orig) - len(PML_XY_elms)

    print("Number of PML_X elements: ", len(PML_X_elms_orig), " -> ", len(PML_X_elms), " (", diff_x, " elements removed)")
    print("Number of PML_Y elements: ", len(PML_Y_elms_orig), " -> ", len(PML_Y_elms), " (", diff_y, " elements removed)")
    print("Number of PML_XY elements: ", len(PML_XY_elms_orig), " -> ", len(PML_XY_elms), " (", diff_xy, " elements removed)")

    return PML_X_elms, PML_Y_elms, PML_XY_elms


@numba.jit(nopython=True)
def write_str_surface(cells_line, cells_quad, bound_edges, _node_ids, flag_abs):
    """ Search element id from edge id and write to string list
        This function may be a heavy bottleneck for large meshes, as
        it is a double loop over the number of edges and elements.
        Thus it is implemented in numba to speed up the process.

        cells_line: list of all line cells
        cells_quad: list of all quad cells (elements)
        bound_edges: list of target line cells on one boundary (e.g. on top bound)
        bound_nodes: list of node ids in one single element (node ordering rule)
    """

    str_lines = []

    for line_elm in bound_edges:
        # get node ids from edge id
        nodes_ids_edge = cells_line[line_elm]
        node_0 = nodes_ids_edge[0]
        node_1 = nodes_ids_edge[1]

        # get element id and node ids from edge id
        for ielm, _node_ids_elm in enumerate(cells_quad):
            if node_0 in _node_ids_elm and node_1 in _node_ids_elm:
                break
        else:
            # If no break occurred, continue to the next iteration
            continue

        # specfem requires only 2 two edge nodes
        out_str = f"{ielm} 2 {node_0} {node_1}"

        if flag_abs is not None:
            out_str += f" {flag_abs}"

        str_lines.append(out_str)


    return str_lines



@numba.jit(nopython=True)
def get_cpml_data(mesh_sets_x, mesh_sets_y, mesh_sets_xy, cell_id_offset):

    n_elm_pml = 0
    str_lines = []

    n_elm_pml += len(mesh_sets_x)
    if (len(mesh_sets_x) > 0):
        for elm in mesh_sets_x:
            str_lines.append(str(int(elm-cell_id_offset+1)) + " " + str(1))

    n_elm_pml += len(mesh_sets_y)
    if (len(mesh_sets_y) > 0):
        for elm in mesh_sets_y:
            str_lines.append(str(int(elm-cell_id_offset+1)) + " " + str(2))

    n_elm_pml += len(mesh_sets_xy)
    if (len(mesh_sets_xy) > 0):
        for elm in mesh_sets_xy:
            str_lines.append(str(int(elm-cell_id_offset+1)) + " " + str(3))

    # add number of pml elements to the first line
    str_lines.insert(0, str(n_elm_pml))

    return str_lines


class Meshio2Specfem2D:

    outdir = "./MESH"

    # node ordering in meshio is the same as vtk
    # https://raw.githubusercontent.com/Kitware/vtk-examples/gh-pages/src/Testing/Baseline/Cxx/GeometricObjects/TestIsoparametricCellsDemo.png

    fhead_Nodes = "Nodes"
    fhead_Mesh = "Mesh"
    fhead_Material = "Material"
    fhead_Surf_abs = "Surf_abs"
    fhead_Surf_free = "Surf_free"
    fhead_CPML = "EltPML"
    fname_out = "TEST"

    fname_Nodes = ""
    fname_Mesh = ""
    fname_Material = ""
    fname_Surf_abs = ""
    fname_Surf_free = ""
    fname_CPML = ""

    n_nodes = 0
    n_cells = 0
    n_edges = 0
    if_second_order = False
    key_line = "line" # line3 for second order
    key_quad = "quad" # quad9 for second order

    # stacy boundary flags
    top_abs = False
    bot_abs = True
    left_abs = True
    right_abs = True

    use_cpml = False
    cell_id_offset = 0

    # pml transition layer
    pml_transition_layer = True


    def __init__(self, mesh, top_abs=False, bot_abs=True, left_abs=True, right_abs=True):
        self.top_abs = top_abs
        self.bot_abs = bot_abs
        self.left_abs = left_abs
        self.right_abs = right_abs
        self.mesh = mesh

        # create output directory
        import os
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # check if PML_X is included in the physical groups
        if "PML_X" in self.mesh.cell_sets_dict:
            self.use_cpml = True
            # When CPML is used, all outer PML boundaries need absorbing edges
            # for Dirichlet boundary conditions (see pml_init.F90:
            # determine_boundary_abs_points_PML). Auto-set all abs flags to True.
            self.top_abs = True
            self.bot_abs = True
            self.left_abs = True
            self.right_abs = True

        # check if second order elements are included in the mesh
        if "quad9" in self.mesh.cells_dict:
            self.if_second_order = True

        if self.if_second_order:
            self.key_line = "line3"
            self.key_quad = "quad9"
        else:
            self.key_line = "line"
            self.key_quad = "quad"


    def write_nodes(self):
        # number of nodes
        # coord x y
        nodes = self.mesh.points
        self.n_nodes = len(nodes)
        with open(self.fname_Nodes, "w") as f:
            f.write(f"{self.n_nodes}\n")
            np.savetxt(f, nodes[:,0:2], fmt="%f %f")


    def write_mesh(self):
        # number of elements
        # node id1 id2 id3 id4 (id5 id6 id7 id8 id9 for second order elements)

        with open(self.fname_Mesh, "w") as f:
            cell_data = self.mesh.cells_dict[self.key_quad] # id starts from 0
            self.n_cells = len(cell_data)
            f.write(str(self.n_cells) + "\n")
            if self.if_second_order:
                np.savetxt(f, cell_data[:,0:9], fmt="%d %d %d %d %d %d %d %d %d")
            else:
                np.savetxt(f, cell_data[:,0:4], fmt="%d %d %d %d")


    def write_material(self):

        # find all keys starting with "M" from mesh.cell_sets
        M_keys = [key for key in self.mesh.cell_sets_dict if key[0] == "M"]
        print("material keys: ", M_keys)

        # reconstruct material flag array
        arr_mflag = np.ones(self.n_cells, dtype=int) * -1

        # id offset for quad (subtract the number of lines)
        self.cell_id_offset = int(np.min(self.mesh.cell_sets_dict["M1"][self.key_quad]))

        print("cell_id_offset: ", self.cell_id_offset)

        for key in M_keys:
            for cell in self.mesh.cell_sets_dict[key][self.key_quad]:
                try:
                    icell = int(cell)-self.cell_id_offset
                    # skip the first character "M"
                    arr_mflag[icell] = int(key.lstrip("M"))
                except:
                    print(
                        "icell, cell, key, key.lstrip('M'), self.n_cells, len(arr_mflag)")
                    print(icell, cell, key, key.lstrip(
                        "M"), self.n_cells, len(arr_mflag))
                    raise Exception("Error: cell id out of range")

        # check if all cells are assigned to a material
        if -1 in arr_mflag:
            raise Exception("Error: not all cells are assigned to a material")

        # write material file
        with open(self.fname_Material, "w") as f:
            np.savetxt(f, arr_mflag, fmt="%d")


    def _write_str_surface(self, bound_edges, _node_ids, flag_abs=None):
        # pipline function to avoid passing dict to numba (unsupported)

        cells_line = self.mesh.cells_dict[self.key_line]
        cells_quad = self.mesh.cells_dict[self.key_quad]

        return write_str_surface(cells_line, cells_quad, bound_edges, _node_ids, flag_abs)


    def write_surf_free_and_abs(self):
        # free surface
        # number of free surface edges
        # elemnt id, num nodes, node id1, node id2 (, node id3 for second order elements)

        # abs boundary
        # number of abs boundary edges
        # elemnt id, num nodes, node id1, node id2 (, node id3 for second order elements), boundary flag (1=bottom, 2=right, 3=top, 4=left)

        # node id order in meshio is the same as vtk
        #
        #
        # 3      6      2
        #
        #
        # 7      8      5
        #
        #
        # 0      4      1

        bot_nodes = np.array([0, 1, 4])
        right_nodes = np.array([1, 2, 5])
        top_nodes = np.array([2, 3, 6])
        left_nodes = np.array([3, 0, 7])

        # write number of free surface edges
        edges_top = self.mesh.cell_sets_dict["Top"][self.key_line]
        edges_bot = self.mesh.cell_sets_dict["Bottom"][self.key_line]
        edges_left = self.mesh.cell_sets_dict["Left"][self.key_line]
        edges_right = self.mesh.cell_sets_dict["Right"][self.key_line]

        n_edges_free = 0
        n_edges_abs = 0

        str_lines = []

        # write number of free surface edges
        if self.bot_abs != True:
            str_lines.extend(self._write_str_surface(edges_bot, bot_nodes))
        if self.right_abs != True:
            str_lines.extend(self._write_str_surface(edges_right, right_nodes))
        if self.top_abs != True:
            str_lines.extend(self._write_str_surface(edges_top, top_nodes))
        if self.left_abs != True:
            str_lines.extend(self._write_str_surface(edges_left, left_nodes))

        n_edges_free = len(str_lines)

        # add number of free surface edges to the first line
        str_lines.insert(0, str(n_edges_free))

        np.savetxt(self.fname_Surf_free, str_lines, fmt="%s")

        # write number of abs boundary edges
        str_lines = []
        if self.bot_abs:
            str_lines.extend(self._write_str_surface(edges_bot, bot_nodes, 1))
        if self.right_abs:
            str_lines.extend(self._write_str_surface(edges_right, right_nodes, 2))
        if self.top_abs:
            str_lines.extend(self._write_str_surface(edges_top, top_nodes, 3))
        if self.left_abs:
            str_lines.extend(self._write_str_surface(edges_left, left_nodes, 4))
        n_edges_abs = len(str_lines)

        # add number of abs boundary edges to the first line
        str_lines.insert(0, str(n_edges_abs))

        np.savetxt(self.fname_Surf_abs, str_lines, fmt="%s")


    def write_cpml(self):
        """
        Write cpml data to file.
        In default, this function creates a layer of damping elements between
        the PML and the main domain.
        """

        if self.pml_transition_layer == False:
            cells_PML_X  = self.mesh.cell_sets_dict["PML_X"][self.key_quad]
            cells_PML_Y  = self.mesh.cell_sets_dict["PML_Y"][self.key_quad]
            cells_PML_XY = self.mesh.cell_sets_dict["PML_XY"][self.key_quad]
        else:
            # prepare list of cells for cpml excluding the damping layers
            _cells_inner_boundary = []

            # find keys which are _Top, _Bottom, _Left, _Right
            for key in self.mesh.cell_sets_dict:
                if key in ["_Top", "_Bottom", "_Left", "_Right"]:
                    _cells_inner_boundary.append(self.mesh.cell_sets_dict[key][self.key_line])

            # concatenate a list of numpy arrays to one numpy array
            cells_line_inner_boundary = np.concatenate(_cells_inner_boundary)

            cells_quad_total = self.mesh.cells_dict[self.key_quad]
            cells_line_total = self.mesh.cells_dict[self.key_line]

            cells_PML_X, cells_PML_Y, cells_PML_XY = get_cpml_cells_except_damping(
                self.mesh.cell_sets_dict["PML_X"][self.key_quad],
                self.mesh.cell_sets_dict["PML_Y"][self.key_quad],
                self.mesh.cell_sets_dict["PML_XY"][self.key_quad],
                cells_quad_total,
                cells_line_total,
                cells_line_inner_boundary,
                self.cell_id_offset
            )

        # n_elme pml
        # elm_id cpml_flag
        str_lines = get_cpml_data(cells_PML_X
                                , cells_PML_Y
                                , cells_PML_XY
                                , self.cell_id_offset)

        np.savetxt(self.fname_CPML, str_lines, fmt="%s")


    def write(self, filename_out="TEST", pml_transition_layer=False):

        # measure time
        start_time = time.time()

        # set output file name
        self.filename_out = filename_out

        # set pml transition layer
        self.pml_transition_layer = pml_transition_layer

        # construct file names
        self.fname_Nodes = os.path.join(self.outdir, self.fhead_Nodes + "_" + self.filename_out)
        self.fname_Mesh = os.path.join(self.outdir, self.fhead_Mesh + "_" + self.filename_out)
        self.fname_Material = os.path.join(self.outdir, self.fhead_Material + "_" + self.filename_out)
        self.fname_Surf_abs = os.path.join(self.outdir, self.fhead_Surf_abs + "_" + self.filename_out)
        self.fname_Surf_free = os.path.join(self.outdir, self.fhead_Surf_free + "_" + self.filename_out)
        self.fname_CPML = os.path.join(self.outdir, self.fhead_CPML + "_" + self.filename_out)

        # write mesh files in specfem format from meshio object
        t_start_node = time.time()
        self.write_nodes()
        t_end_node = time.time()
        print("Time elapsed for nodes: ", t_end_node - t_start_node, " seconds")

        # write mesh
        t_start_mesh = time.time()
        self.write_mesh()
        t_end_mesh = time.time()
        print("Time elapsed for mesh: ", t_end_mesh - t_start_mesh, " seconds")

        # write material
        t_start_material = time.time()
        self.write_material()
        t_end_material = time.time()
        print("Time elapsed for material: ", t_end_material - t_start_material, " seconds")

        # write free surface
        t_start_surf = time.time()
        self.write_surf_free_and_abs()
        t_end_surf = time.time()
        print("Time elapsed for surf: ", t_end_surf - t_start_surf, " seconds")

        # write cpml
        if self.use_cpml:
            t_start_cpml = time.time()
            self.write_cpml()
            t_end_cpml = time.time()
            print("Time elapsed for cpml: ", t_end_cpml - t_start_cpml, " seconds")

        print("Time elapsed: ", time.time() - start_time, " seconds")

