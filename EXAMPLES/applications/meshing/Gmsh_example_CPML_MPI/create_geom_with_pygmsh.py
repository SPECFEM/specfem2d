import pygmsh
import numpy as np
import sys

class line_store:
    def __init__(self, rect_id, line_id, point_ids):
        self.rect_id = rect_id
        self.point_ids = point_ids # [point_id1, point_id2]

class rectangles:


    list_rects = [] # list of one_rect objects
    dict_points_created = {} # key: (x, y, z), value: [id_rect, id_point]
    dict_lines_created = {} # key: (xy1, xy2), value: [id_rect, id_line]

    n_rects = 0

    def __init__(self, geom):
        # reset list and dict
        self.list_rects = []
        self.dict_points_created = {}
        self.dict_lines_created = {}

    def add_one_rect(self, geom, xy1, xy2, xy3, xy4, lc, nelm_h=None, nelm_v=None, topo=None,
                     transfinite=False, mat_tag=None):
        # add z component to the coordinates
        xy1 = (xy1[0], xy1[1], 0)
        xy2 = (xy2[0], xy2[1], 0)
        xy3 = (xy3[0], xy3[1], 0)
        xy4 = (xy4[0], xy4[1], 0)

        transfinite_1d = False
        if nelm_h is not None or nelm_v is not None:
            transfinite_1d = True

        self.list_rects.append(one_rect(geom, self.n_rects, xy1, xy2, xy3, xy4, lc,
                                        transfinite=transfinite, transfinite_1d=transfinite_1d,
                                        mat_tag=mat_tag, nelm_h=nelm_h, nelm_v=nelm_v,
                                        topo=topo))
        self.n_rects += 1

    def add_one_rect_1d(self, geom, xy1, xy2, xy3, xy4, lc, nelm_h, nelm_v, transfinite=True,
                        mat_tag=None, pml_tag=None, bound_tag=None):
        # add z component to the coordinates
        xy1 = (xy1[0], xy1[1], 0)
        xy2 = (xy2[0], xy2[1], 0)
        xy3 = (xy3[0], xy3[1], 0)
        xy4 = (xy4[0], xy4[1], 0)

        self.list_rects.append(one_rect(geom, self.n_rects, xy1, xy2, xy3, xy4, lc,
                                        nelm_h=nelm_h, nelm_v=nelm_v, transfinite=transfinite, transfinite_1d=transfinite,
                                        mat_tag=mat_tag, pml_tag=pml_tag, bound_tag=bound_tag))
        self.n_rects += 1


    # add pml layers
    # Currently only works for one column of rectangles
    # but it is easy to update for multiple columns, which only requires
    # additional lines to detect the leftmost and rightmost rectangles
    # + multiple top/bottom rectangles
    def add_pml_layers(self, geom, w_pml, lc_pml, top_pml=False):
        # find the bottom and top rectangle
        ave_y_min = 1e10
        ave_y_max = -1e10
        irect_bottom = -1
        irect_top = -1
        nelm_pml = int(w_pml / lc_pml)+1

        for reg_rect in self.list_rects:
            # calculate the average y coordinate
            ave_y = (reg_rect.list_xy[0][1] + reg_rect.list_xy[1][1] + reg_rect.list_xy[2][1] + reg_rect.list_xy[3][1]) / 4
            if ave_y < ave_y_min:
                ave_y_min = ave_y
                irect_bottom = reg_rect.id
            if ave_y > ave_y_max:
                ave_y_max = ave_y
                irect_top = reg_rect.id

        # error if top or bottom rectangle is not found
        if irect_bottom == -1:
            print("Error: bottom rectangle is not found")
            sys.exit()
        if irect_top == -1:
            print("Error: top rectangle is not found")
            sys.exit()

        print("bottom rectangle id", irect_bottom)
        print("top rectangle id", irect_top)

        # add top tag if top_pml is False
        if not top_pml:
            self.list_rects[irect_top].bound_tag = "TopFree"

        # side of the main rectangles
        tmp_list_rects = self.list_rects.copy()
        for rect in tmp_list_rects:
            nelm_v = rect.nelm_v
            # left side
            self.add_one_rect_1d(geom,(rect.list_xy[0][0]-w_pml, rect.list_xy[0][1], 0),
                                      (rect.list_xy[0][0],       rect.list_xy[0][1], 0),
                                      (rect.list_xy[3][0],       rect.list_xy[3][1], 0),
                                      (rect.list_xy[3][0]-w_pml, rect.list_xy[3][1], 0),
                                       rect.list_lc, nelm_pml, nelm_v,
                                       mat_tag=rect.mat_tag, pml_tag="PML_X", bound_tag="Left")

            # right side
            self.add_one_rect_1d(geom, (rect.list_xy[1][0],       rect.list_xy[1][1], 0),
                                       (rect.list_xy[1][0]+w_pml, rect.list_xy[1][1], 0),
                                       (rect.list_xy[2][0]+w_pml, rect.list_xy[2][1], 0),
                                       (rect.list_xy[2][0],       rect.list_xy[2][1], 0),
                                        rect.list_lc, nelm_pml, nelm_v,
                                        mat_tag=rect.mat_tag, pml_tag="PML_X", bound_tag="Right")

        # add pml layers
        # bottom left
        xy_base = self.list_rects[irect_bottom].list_xy[0]
        list_lc_bot = self.list_rects[irect_bottom].list_lc
        lc_corner = list_lc_bot[0]
        self.add_one_rect_1d(geom, (xy_base[0]-w_pml, xy_base[1]-w_pml, 0),
                                  (xy_base[0],       xy_base[1]-w_pml, 0),
                                  (xy_base[0],       xy_base[1],       0),
                                  (xy_base[0]-w_pml, xy_base[1],       0),
                                  lc_corner, nelm_pml, nelm_pml-2,
                                  mat_tag=self.list_rects[irect_bottom].mat_tag, pml_tag="PML_XY",
                                  bound_tag="Bottom,Left")

        # bottom
        xy_base_left  = self.list_rects[irect_bottom].list_xy[0]
        xy_base_right = self.list_rects[irect_bottom].list_xy[1]
        nelm_h = self.list_rects[irect_bottom].nelm_h
        self.add_one_rect_1d(geom, (xy_base_left[0],  xy_base_left[1] -w_pml, 0),
                                   (xy_base_right[0], xy_base_right[1]-w_pml, 0),
                                   (xy_base_right[0], xy_base_right[1]      , 0),
                                   (xy_base_left[0],  xy_base_left[1]       , 0),
                                   list_lc_bot, nelm_h, nelm_pml-2,
                                   mat_tag=self.list_rects[irect_bottom].mat_tag, pml_tag="PML_Y",
                                   bound_tag="Bottom")

        # bottom right
        xy_base = self.list_rects[irect_bottom].list_xy[1]
        lc_corner = list_lc_bot[1]
        self.add_one_rect_1d(geom, (xy_base[0],       xy_base[1]-w_pml, 0),
                                  (xy_base[0]+w_pml, xy_base[1]-w_pml, 0),
                                  (xy_base[0]+w_pml, xy_base[1],       0),
                                  (xy_base[0],       xy_base[1],       0),
                                  lc_corner, nelm_pml, nelm_pml-2,
                                  mat_tag=self.list_rects[irect_bottom].mat_tag, pml_tag="PML_XY",
                                  bound_tag="Bottom,Right")

        if top_pml:
            list_lc_top = self.list_rects[irect_top].list_lc
            lc_corner = list_lc_top[3]
            # top left
            xy_base = self.list_rects[irect_top].list_xy[3]
            self.add_one_rect_1d(geom, (xy_base[0]-w_pml, xy_base[1],       0),
                                      (xy_base[0],       xy_base[1],       0),
                                      (xy_base[0],       xy_base[1]+w_pml, 0),
                                      (xy_base[0]-w_pml, xy_base[1]+w_pml, 0),
                                      lc_corner, nelm_pml, nelm_pml,
                                      mat_tag=self.list_rects[irect_top].mat_tag, pml_tag="PML_XY", bound_tag="Top,Left")
            # top
            xy_base_left = self.list_rects[irect_top].list_xy[3]
            xy_base_right = self.list_rects[irect_top].list_xy[2]
            nelm_h = self.list_rects[irect_top].nelm_h
            self.add_one_rect_1d(geom, (xy_base_left[0],  xy_base_left[1],        0),
                                      (xy_base_right[0], xy_base_right[1],       0),
                                      (xy_base_right[0], xy_base_right[1]+w_pml, 0),
                                      (xy_base_left[0],  xy_base_left[1] +w_pml, 0),
                                      list_lc_top, nelm_h, nelm_pml,
                                      mat_tag=self.list_rects[irect_top].mat_tag, pml_tag="PML_Y", bound_tag="Top")
            # top right
            xy_base = self.list_rects[irect_top].list_xy[2]
            self.add_one_rect_1d(geom, (xy_base[0],       xy_base[1],       0),
                                      (xy_base[0]+w_pml, xy_base[1],       0),
                                      (xy_base[0]+w_pml, xy_base[1]+w_pml, 0),
                                      (xy_base[0],       xy_base[1]+w_pml, 0),
                                      list_lc_top, nelm_pml, nelm_pml,
                                      mat_tag=self.list_rects[irect_top].mat_tag, pml_tag="PML_XY", bound_tag="Top,Right")


    def build_points_edges(self, geom, debug=False):

        # loop each one_rect object and creat the point line elements
        # the register the created objects to a list the use it to check if
        # each point/line should be newly created or passed from the adjacent rectangle
        for rect in self.list_rects:
            for ixy, xy in enumerate(rect.list_xy):
                if xy not in self.dict_points_created:
                    # create a point
                    rect.list_points[ixy] = geom.add_point(xy, rect.list_lc[ixy])
                    self.dict_points_created[xy] = [rect.id, ixy]
                    #print ("created point", xy, "id", rect.list_points[ixy] )
                    #print ("stored as ", self.dict_points_created[xy], " acccessed as ", self.list_rects[rect.id].list_points[ixy])
                else:
                    # use precreated point
                    id_rect, id_point = self.dict_points_created[xy]
                    rect.list_points[ixy] = self.list_rects[id_rect].list_points[id_point]
                    #print ("used point", xy, "id", rect.list_points[ixy] )

        # debug
        if debug:
            print("dict_points_created", self.dict_points_created)
            for rect in self.list_rects:
                print("rect.id", rect.id)
                print("rect.list_points", rect.list_points)

        # create lines
        for rect in self.list_rects:
            for iline, (xy1, xy2) in enumerate([(rect.list_xy[0], rect.list_xy[1]),
                                                (rect.list_xy[1], rect.list_xy[2]),
                                                (rect.list_xy[2], rect.list_xy[3]),
                                                (rect.list_xy[3], rect.list_xy[0])]):
                if debug:
                    print ("rect.id, iline", rect.id, iline)

                if (xy1, xy2) not in self.dict_lines_created and (xy2, xy1) not in self.dict_lines_created:
                    # create a line
                    ipt = iline
                    ipt2 = (iline + 1) % 4 # point id is 0, 1, 2, 3

                    if rect.topo is None or iline != 2: # topo only top line
                        rect.list_lines[iline] = geom.add_line(rect.list_points[ipt], rect.list_points[ipt2])
                    elif rect.topo is not None and iline == 2:
                        try:
                            # create a line with topology
                            x_arr = rect.topo["x"]
                            y_arr = rect.topo["z"]
                        except:
                            print("Error: topo is not set. Set topo as a dictionary with x and z keys")
                            sys.exit()

                        #try:
                        p_list = [geom.add_point((x, z, 0), rect.list_lc[2]) for x, z in zip(x_arr, y_arr)]
                        p_list = p_list[::-1]
                        # add p3 and p4 to the list
                        p_list.insert(0, rect.list_points[2])
                        p_list.append(rect.list_points[3])
                        #p_list[0] = rect.list_points[2]
                        #p_list[-1] = rect.list_points[3]
                        rect.list_lines[iline] = geom.add_spline(p_list)

                        #except:
                        #    print("Error: failed to create a line with topo")
                        #    sys.exit()

                    self.dict_lines_created[(xy1, xy2)] = [rect.id, iline]
                else:
                    if (xy1, xy2) in self.dict_lines_created:
                        id_rect, id_line = self.dict_lines_created[(xy1, xy2)]
                        rect.list_lines[iline] = self.list_rects[id_rect].list_lines[id_line]
                    else:
                        id_rect, id_line = self.dict_lines_created[(xy2, xy1)]
                        rect.list_lines[iline] = -self.list_rects[id_rect].list_lines[id_line]
                    if debug:
                        print("taken from id_rect, id_line", id_rect, id_line)

        # create line loops and surfaces
        for rect in self.list_rects:
            print("creating surface for rect.id", rect.id)
            rect.create_surface(geom)

        # do recombile
        tmp_ps = []
        for rect in self.list_rects:
            tmp_ps.append(rect.ps)
        geom.set_recombined_surfaces(tmp_ps)

        # create physical groups
        # set material groups
        list_mat_tags = []
        for rect in self.list_rects:
            if rect.mat_tag not in list_mat_tags:
                list_mat_tags.append(rect.mat_tag)

        for mat_tag in list_mat_tags:
            list_ps = []
            for rect in self.list_rects:
                if rect.mat_tag == mat_tag:
                    list_ps.append(rect.ps)
            geom.add_physical(list_ps, label=mat_tag)

        # set pml groups (PML_X, PML_Y, PML_XY)
        list_pml_tags = ["PML_X", "PML_Y", "PML_XY"]
        for pml_tag in list_pml_tags:
            list_ps = []
            for rect in self.list_rects:
                if rect.pml_tag == pml_tag:
                    list_ps.append(rect.ps)
            print("pml_tag", pml_tag, "list_ps", list_ps)
            if (len(list_ps) > 0):
                geom.add_physical(list_ps, label=pml_tag)
            else:
                # warning if no pml is found
                print("Warning: no PML is found for ", pml_tag)

        # set boundary groups
        list_l_left = []
        list_l_right = []
        list_l_top = []
        list_l_bottom = []
        list_l_left_inner = []
        list_l_right_inner = []
        list_l_top_inner = []
        list_l_bottom_inner = []
        for rect in self.list_rects:
            if rect.bound_tag == "Left":
                list_l_left.append(rect.list_lines[3])
                list_l_left_inner.append(rect.list_lines[1])
            elif rect.bound_tag == "Right":
                list_l_right.append(rect.list_lines[1])
                list_l_right_inner.append(rect.list_lines[3])
            elif rect.bound_tag == "Top":
                list_l_top.append(rect.list_lines[2])
                list_l_top_inner.append(rect.list_lines[0])
            elif rect.bound_tag == "TopFree":
                # Top boundary of the domain when top_pml is False.
                # Only add outer edge to Top group; do NOT create _Top inner
                # boundary since there is no top PML layer.
                list_l_top.append(rect.list_lines[2])
            elif rect.bound_tag == "Bottom":
                list_l_bottom.append(rect.list_lines[0])
                list_l_bottom_inner.append(rect.list_lines[2])
            elif rect.bound_tag == "Bottom,Left":
                list_l_bottom.append(rect.list_lines[0])
                list_l_left.append(rect.list_lines[3])
            elif rect.bound_tag == "Bottom,Right":
                list_l_bottom.append(rect.list_lines[0])
                list_l_right.append(rect.list_lines[1])
            elif rect.bound_tag == "Top,Left":
                list_l_top.append(rect.list_lines[2])
                list_l_left.append(rect.list_lines[3])
            elif rect.bound_tag == "Top,Right":
                list_l_top.append(rect.list_lines[2])
                list_l_right.append(rect.list_lines[1])

        if len(list_l_left) > 0:
            geom.add_physical(list_l_left, label="Left")
        if len(list_l_right) > 0:
            geom.add_physical(list_l_right, label="Right")
        if len(list_l_bottom) > 0:
            geom.add_physical(list_l_bottom, label="Bottom")
        if len(list_l_top) > 0:
            geom.add_physical(list_l_top, label="Top")

        if len(list_l_left_inner) > 0:
            geom.add_physical(list_l_left_inner, label="_Left")
        if len(list_l_right_inner) > 0:
            geom.add_physical(list_l_right_inner, label="_Right")
        if len(list_l_bottom_inner) > 0:
            geom.add_physical(list_l_bottom_inner, label="_Bottom")
        if len(list_l_top_inner) > 0:
            geom.add_physical(list_l_top_inner, label="_Top")


class one_rect:


    def __init__(self, geom, id_rect, xy1, xy2, xy3, xy4, lc,
                 transfinite=False, transfinite_1d=False,
                 nelm_h=None, nelm_v=None,
                 mat_tag=None, pml_tag=None, bound_tag=None, topo=None):
        # xy1, xy2, xy3, xy4: four corners of the rectangle
        # e.g. xy1 = (0,0,0) # dummy z
        # points should be placed as
        """
        xy4 -l3 - xy3
        |         |
        l4       l2
        |         |
        xy1 -l1- xy2
        """
        self.id = id_rect
        self.list_xy = [xy1, xy2, xy3, xy4]
        self.list_points = [None, None, None, None]
        self.list_lines = [None, None, None, None] # l1, l2, l3, l4
        self.list_line_direction = [1, 1, 1, 1]
        self.ll = None # line loop
        self.ps = None # plane surface
        self.transfinite = transfinite
        self.transfinite_1d = transfinite_1d

        # if nelm_h or nelm_v is scalar, set those values for top and bottom lines
        if isinstance(nelm_h, (int, float)):
            nelm_h = [nelm_h, nelm_h]
        if isinstance(nelm_v, (int, float)):
            nelm_v = [nelm_v, nelm_v]

        self.nelm_h = nelm_h
        self.nelm_v = nelm_v
        # print
        print(f"rectangle id {self.id} is meshed with nelm_h (bottom,top) = {nelm_h} \
                and nelm_v (left, right) = {nelm_v}")

        # material tag
        self.mat_tag = mat_tag
        # pml tag
        self.pml_tag = pml_tag
        # bound tag
        self.bound_tag = bound_tag

        # topo: topology (top)
        self.topo = topo

        # lc: mesh size
        # if lc is a scalar, all points and lines are meshed with the same size
        if isinstance(lc, (int, float)):
            lc1 = lc
            lc2 = lc
            lc3 = lc
            lc4 = lc
        # if lc is a list, each point and line is meshed with different size
        elif isinstance(lc, list):
            # check the length of the list
            if len(lc) != 4:
                print ("Error: length of lc is not 4")
                sys.exit()
            lc1 = lc[0]
            lc2 = lc[1]
            lc3 = lc[2]
            lc4 = lc[3]

        self.list_lc = [lc1, lc2, lc3, lc4]


    def create_surface(self, geom):

        # create line loop
        try:
            self.ll = geom.add_curve_loop(self.list_lines)
        except:
            print("Error: line loop is not created")
            for line in self.list_lines:
                print("line", line)

            sys.exit()
        # create surface
        self.ps = geom.add_plane_surface(self.ll)

        # set recombine
        #geom.set_recombined_surfaces([self.ps])

        if self.transfinite_1d:
            if self.nelm_h is None and self.nelm_v is None:
                print("Error: nelm_h or nelm_v is not set")
                sys.exit()
            if self.nelm_h is not None:
                geom.set_transfinite_curve(self.list_lines[0], self.nelm_h[0], "Progression", 1.0)
                geom.set_transfinite_curve(self.list_lines[2], self.nelm_h[1], "Progression", 1.0)
            if self.nelm_v is not None:
                geom.set_transfinite_curve(self.list_lines[1], self.nelm_v[1], "Progression", 1.0)
                geom.set_transfinite_curve(self.list_lines[3], self.nelm_v[0], "Progression", 1.0)

        # transfinite surface
        if self.transfinite:# or self.transfinite_1d:
            geom.set_transfinite_surface(self.ps, "Left", corner_pts=self.list_points)

