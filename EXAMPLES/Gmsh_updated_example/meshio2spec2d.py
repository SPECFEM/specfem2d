# convert mesh to specfem format
import numpy as np
import numba
import time
import os


@numba.jit(nopython=True)
def write_str_surface(str_lines, cells_line, cells_quad, bound_edges, bound_nodes, flag_abs):

    for line_elm in bound_edges:
        # get node ids from edge id
        nodes_ids_edge = cells_line[line_elm]
        node_0 = nodes_ids_edge[bound_nodes[0]]
        node_1 = nodes_ids_edge[bound_nodes[1]]

        # get element id and node ids from edge id
        for ielm, _node_ids_elm in enumerate(cells_quad):
            if node_0 in _node_ids_elm and node_1 in _node_ids_elm:
                nodes = _node_ids_elm
                break
        else:
            # If no break occurred, continue to the next iteration
            continue

        # specfem requires only 2 two edge nodes
        out_str = f"{ielm} 2 {nodes[0]} {nodes[1]}"

        if flag_abs is not None:
            out_str += f" {flag_abs}"

        str_lines.append(out_str)

    return str_lines



@numba.jit(nopython=True)
def get_cpml_data(mesh_sets_x, mesh_sets_y, mesh_sets_xy, cell_id_offset):

    n_elm_pml = 0
    #pml_x, pml_y, pml_xy = [], [], []

    str_lines = []

    n_elm_pml += len(mesh_sets_x)
    #pml_x = [int(elm - cell_id_offset) for elm in mesh_sets["PML_X"][1]]
    if (len(mesh_sets_x) > 0):
        for elm in mesh_sets_x:
            str_lines.append(str(int(elm-cell_id_offset)) + " " + str(1) + "\n")

    n_elm_pml += len(mesh_sets_y)
    #pml_y = [int(elm - cell_id_offset) for elm in mesh_sets["PML_Y"][1]]
    if (len(mesh_sets_y) > 0):
        for elm in mesh_sets_y:
            str_lines.append(str(int(elm-cell_id_offset)) + " " + str(2) + "\n")

    n_elm_pml += len(mesh_sets_xy)
    #pml_xy = [int(elm - cell_id_offset) for elm in mesh_sets["PML_XY"][1]]
    if (len(mesh_sets_xy) > 0):
        for elm in mesh_sets_xy:
            str_lines.append(str(int(elm-cell_id_offset)) + " " + str(3) + "\n")

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
    key_line = "line"
    key_quad = "quad"

    # stacy boundary flags
    top_abs = False
    bot_abs = True
    left_abs = True
    right_abs = True

    use_cpml = False
    #cpml_flags = {"PML_X": 1, "PML_Y": 2, "PML_XY": 3} # specfem2d/setup/constants.h
    cell_id_offset = 0


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
        if "PML_X" in self.mesh.cell_sets:
            self.use_cpml = True

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
            cell_data = self.mesh.cells_dict[self.key_quad]
            self.n_cells = len(cell_data)
            f.write(str(self.n_cells) + "\n")
            if self.if_second_order:
                np.savetxt(f, cell_data[:,0:9], fmt="%d %d %d %d %d %d %d %d %d")
            else:
                np.savetxt(f, cell_data[:,0:4], fmt="%d %d %d %d")


    def write_material(self):

        # find all keys starting with "M" from mesh.cell_sets
        M_keys = [key for key in self.mesh.cell_sets if key[0] == "M"]
        print("material keys: ", M_keys)

        # reconstruct material flag array
        arr_mflag = np.ones(self.n_cells, dtype=int) * -1

        # id offset for quad (subtract the number of lines)
        self.cell_id_offset = int(np.min(self.mesh.cell_sets["M1"][1]))-1

        print("cell_id_offset: ", self.cell_id_offset)

        for key in M_keys:
            for cell in self.mesh.cell_sets[key][1]:
                try:
                    icell = int(cell)-1-self.cell_id_offset
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


    def _write_str_surface(self, str_lines, bound_edges, bound_nodes, flag_abs=None):
        # pipline function to avoid passing dict to numba (unsupported)

        cells_line = self.mesh.cells_dict[self.key_line]
        cells_quad = self.mesh.cells_dict[self.key_quad]

        return write_str_surface(str_lines, cells_line, cells_quad, bound_edges[0], bound_nodes, flag_abs)


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
        elems_top = self.mesh.cell_sets["Top"]
        elems_bot = self.mesh.cell_sets["Bottom"]
        elems_left = self.mesh.cell_sets["Left"]
        elems_right = self.mesh.cell_sets["Right"]

        n_edges_free = 0
        n_edges_abs = 0

        if self.bot_abs != False:
            n_edges_free += len(elems_bot)
        else:
            n_edges_abs += len(elems_bot)
        if self.right_abs != False:
            n_edges_free += len(elems_right)
        else:
            n_edges_abs += len(elems_right)
        if self.top_abs != False:
            n_edges_free += len(elems_top)
        else:
            n_edges_abs += len(elems_top)
        if self.left_abs != False:
            n_edges_free += len(elems_left)
        else:
            n_edges_abs += len(elems_left)

        str_lines = []

        # write number of free surface edges
        str_lines.append(str(n_edges_free))
        if self.bot_abs != False:
            str_lines = self._write_str_surface(str_lines, elems_bot, bot_nodes)
        if self.right_abs != False:
            str_lines = self._write_str_surface(str_lines, elems_right, right_nodes)
        if self.top_abs != False:
            str_lines = self._write_str_surface(str_lines, elems_top, top_nodes)
        if self.left_abs != False:
            str_lines = self._write_str_surface(str_lines, elems_left, left_nodes)

        np.savetxt(self.fname_Surf_free, str_lines, fmt="%s")

        # write number of abs boundary edges
        str_lines = []
        str_lines.append(str(n_edges_abs))
        if self.bot_abs:
            str_lines = self._write_str_surface(str_lines, elems_bot, bot_nodes, 1)
        if self.right_abs:
            str_lines = self._write_str_surface(str_lines, elems_right, right_nodes, 2)
        if self.top_abs:
            str_lines = self._write_str_surface(str_lines, elems_top, top_nodes, 3)
        if self.left_abs:
            str_lines = self._write_str_surface(str_lines, elems_left, left_nodes, 4)

        np.savetxt(self.fname_Surf_abs, str_lines, fmt="%s")


    def write_cpml(self):
        # n_elme pml
        # elm_id cpml_flag
        str_lines = get_cpml_data(self.mesh.cell_sets["PML_X"][1]
                                , self.mesh.cell_sets["PML_Y"][1]
                                , self.mesh.cell_sets["PML_XY"][1]
                                , self.cell_id_offset)

        np.savetxt(self.fname_CPML, str_lines, fmt="%s")


    def write(self, filename_out="TEST"):

        # measure time
        start_time = time.time()

        self.filename_out = filename_out

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

