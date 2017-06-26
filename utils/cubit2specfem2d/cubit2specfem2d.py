#!/usr/bin/env python
#
# Script to export a Cubit13+/Trelis 2D mesh in specfem2d format
# Author unknown, comments and modifications by Alexis Bottero (alexis dot bottero At gmail dot com)
#
# Create your mesh in Cubit (or build the simpleAxisym2dMesh.py example) and play this script within Cubit as a Python journal file.
# Instructions for mesh creation :
# _The mesh must be in XZ plane!
# _One block per material :
#      cubit.cmd('block 1 name "Acoustic channel" ') # acoustic material region
#      cubit.cmd('block 1 attribute count 6')        # number of attributes
#      cubit.cmd('block 1 attribute index 1 1')    # material index
#      cubit.cmd('block 1 attribute index 2 1500 ')  # vp
#      cubit.cmd('block 1 attribute index 3 0 ')     # vs
#      cubit.cmd('block 1 attribute index 4 1000 ')  # rho
#      cubit.cmd('block 1 attribute index 5 0 ')     # Q_flag
#      cubit.cmd('block 1 attribute index 6 0 ')     # anisotropy_flag
#      cubit.cmd('block 1 element type QUAD4')
#
# _One block per border (abs_bottom, abs_right, abs_left, abs_top, topo, axis). If axisymmetric simulation don't create a block
#  abs_left but a block axis.
#  Ex:
#      cubit.cmd('block 3 edge in surf all with z_coord > -0.1') # topo
#      cubit.cmd('block 3 name "topo"')
#
#_ One block per pml layer of a given type (acoustic or elastic) : pml_x_acoust,pml_z_acoust,pml_xz_acoust,pml_x_elast,pml_z_elast,pml_xz_elast
#      !! Warning !! pml blocks don't have faces in common
#      !! Warning !! you must create the corresponding absorbing surface blocks (abs_bottom, abs_right, abs_left, abs_top)!
#
# Ideas to improve that script (ctrl+f for TODO also): _Allow 2D models built in XY and ZY planes
#
# The names of the block and the entities types must match the ones given during the definition of the class mesh on this file :
# Below :
# class mesh(object,mesh_tools):
#     """ A class to store the mesh """
#     def __init__(self):
#
#!! Warning : a block in cubit != quad !! A block is a group of something (quads, edges, volumes, surfaces...)
# On this case the blocks are used to gather faces corresponding to different materials and edges corresponding to free surfaces,
# absorbing surfaces, topography or axis

import cubit

class mtools(object):
    """docstring for mtools"""
    def __init__(self,frequency,list_surf,list_vp):
        super(mtools, self).__init__()
        self.frequency = frequency
        self.list_surf = list_surf
        self.list_vp = list_vp
        self.ngll = 5
        self.percent_gll = 0.172
        self.point_wavelength = 5
    def __repr__(self):
        txt = 'Meshing for frequency up to '+str(self.frequency)+'Hz\n'
        for surf,vp in zip(self.list_surf,self.list_vp):
            txt = txt+'surface '+str(surf)+', vp ='+str(vp)+'  -> size '+str(self.freq2meshsize(vp)[0])+' -> dt '+str(self.freq2meshsize(vp)[0])+'\n'
        return txt
    def freq2meshsize(self,vp):
        velocity = vp*.5
        self.size = (1/2.5)*velocity/self.frequency*(self.ngll-1)/self.point_wavelength
        self.dt = .4*self.size/vp*self.percent_gll
        return self.size,self.dt
    def mesh_it(self):
        for surf,vp in zip(self.list_surf,self.list_vp):
            command = "surface "+str(surf)+" size "+str(self.freq2meshsize(vp)[0])
            cubit.cmd(command)
            command = "surface "+str(surf)+ 'scheme pave'
            cubit.cmd(command)
            command = "mesh surf "+str(surf)
            cubit.cmd(command)

class block_tools:
    def __int__(self):
        pass
    def create_blocks(self,mesh_entity,list_entity = None,):
        if mesh_entity =='surface':
            txt = ' face in surface '
        elif mesh_entity == 'curve':
            txt = ' edge in curve '
        elif mesh_entity == 'group':
            txt = ' face in group '
        if list_entity:
            if not isinstance(list_entity,list):
                list_entity = [list_entity]
        for entity in list_entity:
            iblock = cubit.get_next_block_id()
            command = "block "+str(iblock)+ txt +str(entity)
            cubit.cmd(command)
    def material_file(self,filename):
        matfile = open(filename,'w')
        material = []
        for record in matfile:
            mat_name,vp_str = record.split()
            vp = float(vp_str)
            material.append([mat_name,vp])
        self.material = dict(material)
    def assign_block_material(self,id_block,mat_name,vp = None):
        try:
            material = self.material
        except:
            material = None
        cubit.cmd('block '+str(id_block)+' attribute count 2')
        cubit.cmd('block '+str(id_block)+'  attribute index 1 '+str(id_block))
        if material:
            if material.has_key(mat_name):
                cubit.cmd('block '+str(id_block)+'  attribute index 2 '+str(material[mat_name]))
                print 'block '+str(id_block)+' - material '+mat_name+' - vp '+str(material[mat_name])+' from database'
        elif vp:
            cubit.cmd('block '+str(id_block)+'  attribute index 2 '+str(vp))
            print 'block '+str(id_block)+' - material '+mat_name+' - vp '+str(vp)
        else:
            print 'assignment impossible: check if '+mat_name+' is in the database or specify vp'

class mesh_tools(block_tools):
    """Tools for the mesh
    #########
    dt,edge_dt,freq,edge_freq = seismic_resolution(edges,velocity,bins_d = None,bins_u = None,sidelist = None,ngll = 5,np = 8)
        Given the velocity of a list of edges, seismic_resolution provides the minimum Dt required for the stability condition (and the corrisponding edge).
        Furthermore, given the number of gll point in the element (ngll) and the number of GLL point for wavelength, it provide the maximum resolved frequency.
    #########
    length = edge_length(edge)
        return the length of a edge
    #########
    edge_min,length = edge_min_length(surface)
        given the cubit id of a surface, it return the edge with minimun length
    #########
    """
    def __int__(self):
        pass
    def seismic_resolution(self,edges,velocity,bins_d = None,bins_u = None,sidelist = None):
        """
        dt,edge_dt,freq,edge_freq = seismic_resolution(edges,velocity,bins_d = None,bins_u = None,sidelist = None,ngll = 5,np = 8)
            Given the velocity of a list of edges, seismic_resolution provides the minimum Dt required for the stability condition (and the corrisponding edge).
            Furthermore, given the number of gll point in the element (ngll) and the number of GLL point for wavelength, it provide the maximum resolved frequency.
        """
        ratiostore = 1e10
        dtstore = 1e10
        edgedtstore = -1
        edgeratiostore = -1
        for edge in edges:
            d = self.edge_length(edge)
            ratio = (1/2.5)*velocity/d*(self.ngll-1)/self.point_wavelength
            dt = .4*d/velocity*self.percent_gll
            if dt<dtstore:
               dtstore = dt
               edgedtstore = edge
            if ratio < ratiostore:
                ratiostore = ratio
                edgeratiostore = edge
            try:
                for bin_d,bin_u,side in zip(bins_d,bins_u,sidelist):
                    if ratio >= bin_d and ratio < bin_u:
                        command = "sideset "+str(side)+" edge "+str(edge)
                        cubit.cmd(command)
                        #print command
                        break
            except:
                pass
        return dtstore,edgedtstore,ratiostore,edgeratiostore
    def edge_length(self,edge):
        """
        length = edge_length(edge)
            return the length of a edge
        """
        from math import sqrt
        nodes = cubit.get_connectivity('edge',edge)
        x0,y0,z0 = cubit.get_nodal_coordinates(nodes[0])
        x1,y1,z1 = cubit.get_nodal_coordinates(nodes[1])
        d = sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
        return d
    def edge_min_length(self,surface):
        """
        edge_min,length = edge_min_length(surface)
            given the cubit id of a surface, it return the edge with minimun length
        """
        from math import sqrt
        self.dmin = 99999
        edge_store = 0
        command = "group 'list_edge' add edge in surf "+str(surface)
        command = command.replace("["," ").replace("]"," ")
        #print command
        cubit.cmd(command)
        group = cubit.get_id_from_name("list_edge")
        edges = cubit.get_group_edges(group)
        command = "delete group "+ str(group)
        cubit.cmd(command)
        for edge in edges:
            d = self.edge_length(edge)
            if d<dmin:
                self.dmin = d
                edge_store = edge
        self.edgemin = edge_store
        return self.edgemin,self.dmin
    def jac_check(self,nodes):
        x0 = cubit.get_nodal_coordinates(nodes[0])
        x1 = cubit.get_nodal_coordinates(nodes[1])
        x2 = cubit.get_nodal_coordinates(nodes[2])
        xv1 = x1[0]-x0[0]
        xv2 = x2[0]-x1[0]
        zv1 = x1[2]-x0[2]
        zv2 = x2[2]-x1[2]
        jac = -xv2*zv1+xv1*zv2
        if  jac > 0:
            return nodes
        elif jac < 0:
            return nodes[0],nodes[3],nodes[2],nodes[1]
        else:
            print 'error, jacobian = 0', jac,nodes
    def mesh_analysis(self,frequency):
        cubit.cmd('set info off') # Turn off return messages from Cubit commands
        cubit.cmd('set echo off') # Turn off echo of Cubit commands
        cubit.cmd('set journal off') # Do not save journal file
        bins_d = [0.0001]+range(0,int(frequency)+1)+[1000]
        bins_u = bins_d[1:]
        dt = []
        ed_dt = []
        r = []
        ed_r = []
        nstart = cubit.get_next_sideset_id()
        command = "del sideset all"
        cubit.cmd(command)
        for bin_d,bin_u in zip(bins_d,bins_u):
            nsideset = cubit.get_next_sideset_id()
            command = 'create sideset '+str(nsideset)
            cubit.cmd(command)
            command  =  "sideset "+str(nsideset)+ " name "+ "'ratio-["+str(bin_d)+"_"+str(bin_u)+"['"
            cubit.cmd(command)
        nend = cubit.get_next_sideset_id()
        sidelist = range(nstart,nend)
        for block in self.block_mat:
            name = cubit.get_exodus_entity_name('block',block)
            velocity = self.material[name][1]
            if velocity > 0:
                faces = cubit.get_block_faces(block)
                edges = []
                for face in faces:
                    es = cubit.get_sub_elements("face", face, 1)
                    edges = edges+list(es)
                dtstore,edgedtstore,ratiostore,edgeratiostore = self.seismic_resolution(edges,velocity,bins_d,bins_u,sidelist)
                dt.append(dtstore)
                ed_dt.append(edgedtstore)
                r.append(ratiostore)
                ed_r.append(edgeratiostore)
        self.ddt = zip(ed_dt,dt)
        self.dr = zip(ed_r,r)
        def sorter(x, y):
            return cmp(x[1],y[1])
        self.ddt.sort(sorter)
        self.dr.sort(sorter)
        print self.ddt,self.dr
        print 'Deltat minimum => edge:'+str(self.ddt[0][0])+' dt: '+str(self.ddt[0][1])
        print 'Minimum frequency resolved => edge:'+str(self.dr[0][0])+' frequency: '+str(self.dr[0][1])
        return self.ddt[0],self.dr[0]

class mesh(object,mesh_tools):
    """ A class to store the mesh """
    def __init__(self):
        super(mesh, self).__init__()
        self.mesh_name = 'mesh_file'
        self.axisymmetric_mesh = False # Will be set to true if a group self.pml_boun_name is found
        self.topo_mesh = False # Will be set to true if a group self.topo is found
        self.forcing_mesh = False # Will be set to true if a group self.forcing_boun_name is found
        self.abs_mesh = False # Will be set to true if a group self.pml_boun_name or self.abs_boun_name is found
        self.pml_layers = False # Will be set to true if a group self.pml_boun_name is found
        self.write_nummaterial_velocity_file = False # Will be set to True if 2d blocks have 6 attributes
        self.nodecoord_name = 'nodes_coords_file' # Name of nodes coordinates file to create
        self.material_name = 'materials_file' # Name of material file to create
        self.nummaterial_name = 'nummaterial_velocity_file'
        self.absname = 'absorbing_surface_file' # Name of absorbing surface file to create
        self.forcname = 'forcing_surface_file' # Name of forcing surface file to create
        self.freename = 'free_surface_file' # Name of free surface file to create
        self.pmlname = 'elements_cpml_list' # Name of cpml file to create
        self.axisname = 'elements_axis' # Name of axial elements file to create and name of the block containing axial edges
        self.recname = 'STATIONS'
        self.face = 'QUAD4'  # Faces' type
        self.edge = 'BAR2'   # Edges' type
        self.topo = 'topo'   # Name of the block containing topography edges
        self.pml_boun_name = ['pml_x_acoust','pml_z_acoust','pml_xz_acoust','pml_x_elast','pml_z_elast','pml_xz_elast']  # Name of the block containing pml layers elements
        self.abs_boun_name = ['abs_bottom','abs_right','abs_top','abs_left'] # Name of the block containing absorbing layer edges
        self.forcing_boun_name = ['forcing_bottom','forcing_right','forcing_top','forcing_left'] # Name of the block containing forcing layer edges
        self.abs_boun = []  # block numbers for abs boundaries
        self.pml_boun = []  # block numbers for pml boundaries
        self.nforc = 4 # Maximum number of forcing surfaces (4)
        self.nabs = 4 # Maximum number of absorbing surfaces (4)
        self.rec = 'receivers'
        self.block_definition() # Import blocks features from Cubit
        self.ngll = 5
        self.percent_gll = 0.172
        self.point_wavelength = 5
        cubit.cmd('compress') # Fill the gaps in the numbering of the entities
    def __repr__(self):
        pass
    def block_definition(self):
        """ Import blocks features from Cubit """
        block_flag = [] # Will contain material id (1 if fluid 2 if solid)
        block_mat = [] # Will contain face block ids
        block_bc = [] # Will contain edge block ids
        block_bc_flag = [] # Will contain edge id -> 2
        abs_boun = [-1] * self.nabs # total 4 sides of absorbing boundaries (index 0 : bottom, index 1 : right, index 2 : top, index 3 : left)
        forcing_boun = [-1] * self.nforc # 4 possible forcing boundaries (index 0 : bottom, index 1 : right, index 2 : top, index 3 : left)
        #pml_boun = [-1] * 6 # To store pml layers id (for each pml layer : x_acoust, z_acoust, xz_acoust, x_elast, z_elast, xz_elast)
        pml_boun = [[] for _ in range(6)] # To store the block id corresponding to pml layers id (arbitrary number of blocks for each pml layer : x_acoust, z_acoust, xz_acoust, x_elast, z_elast, xz_elast)
        material = {} # Will contain each material name and their properties
        bc = {} # Will contains each boundary name and their connectivity -> 2
        blocks = cubit.get_block_id_list() # Load the blocks list
        for block in blocks: # Loop on the blocks
            name = cubit.get_exodus_entity_name('block',block) # Contains the name of the blocks
            ty = cubit.get_block_element_type(block) # Contains the block element type (QUAD4...)
            if ty == self.face: # If we are dealing with a block containing faces
                nAttributes = cubit.get_block_attribute_count(block)
                if (nAttributes != 1 and nAttributes != 6):
                    print 'Blocks not properly defined, 2d blocks must have one attribute (material id) or 6 attributes'
                    return None,None,None,None,None,None,None,None
                flag=int(cubit.get_block_attribute_value(block,0)) # Fetch the first attribute value (containing material id)
                print "nAttributes : ",nAttributes
                if nAttributes == 6:
                    self.write_nummaterial_velocity_file = True
                    velP = cubit.get_block_attribute_value(block,1)  # Fetch the first attribute value (containing P wave velocity)
                    velS = cubit.get_block_attribute_value(block,2)  # Fetch the second attribute value (containing S wave velocity)
                    rho = cubit.get_block_attribute_value(block,3)  # Fetch the third attribute value (containing material density)
                    qFlag = cubit.get_block_attribute_value(block,4)  # Fetch the first attribute value (containing Qflag)
                    anisotropy_flag = cubit.get_block_attribute_value(block,5)  # Fetch the first attribute value (containing anisotropy_flag)
                    # Store (material_id,rho,velP,velS,Qflag,anisotropy_flag) in par :
                    par = tuple([flag,rho,velP,velS,qFlag,anisotropy_flag])
                    material[name] = par # associate the name of the block to its id and properties
                block_flag.append(int(flag)) # Append material id to block_flag
                block_mat.append(block)  # Append block id to block_mat
                for pml_idx,pml_name in enumerate(self.pml_boun_name):
                    if pml_name in name:
                        pml_boun[pml_idx].append(block)
                        self.abs_mesh = True
                        self.pml_layers = True
                    # -> Put it at the correct position in pml_boun
                    # (index 0 : pml_x_acoust, index 1 : pml_z_acoust, index 2 : pml_xz_acoust,
                    #  index 3 : pml_x_elast, index 4 : pml_z_elast, index 5 : pml_xz_elast)
                #if name in self.pml_boun_name : # If the block considered refered to one of the pml layer
                #    self.abs_mesh = True
                #    self.pml_layers = True
                #    pml_boun[self.pml_boun_name.index(name)] = block
                #    # -> Put it at the correct position in pml_boun
                #    # (index 0 : pml_x_acoust, index 1 : pml_z_acoust, index 2 : pml_xz_acoust,
                #    #  index 3 : pml_x_elast, index 4 : pml_z_elast, index 5 : pml_xz_elast)
            elif ty == self.edge: # If we are dealing with a block containing edges
                block_bc_flag.append(2) # Append "2" to block_bc_flag
                block_bc.append(block) # Append block id to block_bc
                bc[name] = 2 # Associate the name of the block with its connectivity : an edge has connectivity = 2
                if name == self.topo:
                    self.topo_mesh = True
                    topography = block # If the block considered refered to topography store its id in "topography"
                if name in self.forcing_boun_name:
                    self.forcing_mesh = True
                    forcing_boun[self.forcing_boun_name.index(name)] = block
                    # -> Put it at the correct position in abs_boun (index 0 : bottom, index 1 : right, index 2 : top, index 3 : left)
                if name == self.axisname:
                    self.axisymmetric_mesh = True
                    axisId = block # AXISYM If the block considered refered to the axis store its id in "axisId"
                if name in self.abs_boun_name : # If the block considered refered to one of the boundaries
                    self.abs_mesh = True
                    abs_boun[self.abs_boun_name.index(name)] = block
                    # -> Put it at the correct position in abs_boun (index 0 : bottom, index 1 : right, index 2 : top, index 3 : left)
            else:
                print 'Blocks not properly defined', ty
                return None,None,None,None,None,None,None,None
        nsets = cubit.get_nodeset_id_list() # Get the list of all nodeset
        if len(nsets) == 0: self.receivers = None # If this list is empty : put None in self.receivers
        for nset in nsets:
            name = cubit.get_exodus_entity_name('nodeset',nset) # Contains the name of the nodeset
            if name == self.rec: # If the name considered match self.rec (receivers)
                self.receivers = nset # Store the id of the nodeset in self.receivers
            else:
                print 'nodeset '+name+' not defined'
                self.receivers = None
        # Store everything in the object :
        try:
            self.block_mat = block_mat
            self.block_flag = block_flag
            self.block_bc = block_bc
            self.block_bc_flag = block_bc_flag
            self.bc = bc
            if self.write_nummaterial_velocity_file:
                self.material = material
            if self.abs_mesh:
                self.abs_boun = abs_boun
            if self.topo_mesh:
                self.topography = topography
            if self.forcing_mesh:
                self.forcing_boun = forcing_boun
            if self.axisymmetric_mesh:
                self.axisId = axisId
            if self.pml_layers:
                self.pml_boun = pml_boun
        except:
            print 'Blocks not properly defined'
#    def tomo(self,flag,vel):
#        vp = vel/1000
#        rho = (1.6612*vp-0.472*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**4)*1000
#        txt = '%3i %1i %20f %20f %20f %1i %1i\n' % (flag,1,rho,vel,vel/(3**.5),0,0)
#        return txt
    def nummaterial_write(self,nummaterial_name):
        """ Write material features on file : nummaterial_name """
        print 'Writing '+nummaterial_name+'.....'
        nummaterial = open(nummaterial_name,'w')  # Create the file "nummaterial_name" and open it
        for block in self.block_mat: # For each 2D block
            name = cubit.get_exodus_entity_name('block',block) # Extract the name of the block
            lineToWrite = str(self.material[name][0])+" 1 "+str(self.material[name][1])+" " \
                         +str(self.material[name][2])+" "+str(self.material[name][3])+" "+str(self.material[name][4])+" " \
                         +str(self.material[name][5])+"\n" # flag rho vp vs rho Qflag anisotropy_flag
            nummaterial.write(lineToWrite)
            #nummaterial.write(self.tomo(self.material[name][0],self.material[name][2]))
        nummaterial.close()
        print 'Ok'
    def mesh_write(self,mesh_name):
        """ Write mesh (quads ids with their corresponding nodes ids) on file : mesh_name """
        meshfile = open(mesh_name,'w')
        print 'Writing '+mesh_name+'.....'
        num_elems = cubit.get_quad_count() # Store the number of elements
        toWritetoFile = [""]*(num_elems+1)
        toWritetoFile[0] = str(num_elems)+'\n'
        #meshfile.write(str(num_elems)+'\n') # Write it on first line
        num_write = 0
        for block,flag in zip(self.block_mat,self.block_flag): # for each 2D block
            quads = cubit.get_block_faces(block) # Import quads ids
            for inum,quad in enumerate(quads): # For each of these quads
                nodes = cubit.get_connectivity('face',quad) # Get the nodes
                nodes = self.jac_check(nodes) # Check the jacobian
                txt = ('%10i %10i %10i %10i\n')% nodes
                toWritetoFile[quad] = txt
                #meshfile.write(txt) # Write a line to mesh file
            num_write = num_write+inum+1
            print 'block', block, 'number of ',self.face,' : ', inum+1
        meshfile.writelines(toWritetoFile)
        meshfile.close()
        print 'Ok num elements/write =',str(num_elems), str(num_write)
    def material_write(self,mat_name):
        """ Write quads material on file : mat_name """
        mat = open(mat_name,'w')
        print 'Writing '+mat_name+'.....'
        num_elems = cubit.get_quad_count() # Store the number of elements
        toWritetoFile = [""]*num_elems
        print 'block_mat:',self.block_mat
        print 'block_flag:',self.block_flag
        for block,flag in zip(self.block_mat,self.block_flag): # for each 2D block
            print 'mat: ',block,' flag: ',flag
            quads = cubit.get_block_faces(block) # Import quads id
            for quad in quads: # For each quad
                toWritetoFile[quad-1] = ('%10i\n') % flag
                #mat.write(('%10i\n') % flag) # Write its id in the file
        mat.writelines(toWritetoFile)
        mat.close()
        print 'Ok'
    def pmls_write(self,pml_name):
        """ Write pml elements on file : mat_name """
        cubit.cmd('set info off') # Turn off return messages from Cubit commands
        cubit.cmd('set echo off') # Turn off echo of Cubit commands
        pml_file = open(pml_name,'w')
        print 'Writing '+pml_name+'.....'
        npml_elements = 0
        #id_element = 0 # Global id
        indexFile = 1
        faces_all = [[] for _ in range(6)]
        for block,flag in zip(self.block_mat,self.block_flag): # For each 2D block
            for ipml in range(6): # ipml = 0,1,2,3,4,5 : for each pml layer (x_acoust, z_acoust, xz_acoust,x_elast, z_elast, xz_elast)
               if block in self.pml_boun[ipml]: # If the block considered correspond to the pml
                    faces_all[ipml] = faces_all[ipml] + list(cubit.get_block_faces(block)) # Concatenation
        npml_elements = sum(map(len, faces_all))
        toWritetoFile = [""]*(npml_elements+1)
        toWritetoFile[0] = '%10i\n' % npml_elements # Print the number of faces on the pmls
        #pml_file.write('%10i\n' % npml_elements) # Print the number of faces on the pmls
        print 'Number of elements in all PMLs :',npml_elements
        for block,flag in zip(self.block_mat,self.block_flag): # For each 2D block
            quads = cubit.get_block_faces(block) # Import quads id
            for quad in quads: # For each quad
                #id_element = id_element+1 # global id of this quad
                for ipml in range(0, 6): # iabs = 0,1,2,3,4,5 : for each pml layer (x_acoust, z_acoust, xz_acoust,x_elast, z_elast, xz_elast)
                    if faces_all[ipml] != []: #type(faces_all[ipml]) is not int: # ~ if there are elements in that pml
                        if quad in faces_all[ipml]: # If this quad is belong to that pml
                          #  nodes = cubit.get_connectivity('face',quad) # Import the nodes describing the quad
                          #  nodes = self.jac_check(list(nodes)) # Check the jacobian of the quad
                            toWritetoFile[indexFile] = ('%10i %10i\n') % (quad,ipml%3+1)
                            indexFile = indexFile + 1
                            #pml_file.write(('%10i %10i\n') % (id_element,ipml%3+1)) # Write its id in the file next to its type
        # ipml%3+1 = 1 -> element belongs to a X CPML layer only (either in Xmin or in Xmax)
        # ipml%3+1 = 2 -> element belongs to a Z CPML layer only (either in Zmin or in Zmax)
        # ipml%3+1 = 3 -> element belongs to both a X and a Y CPML layer (i.e., to a CPML corner)
        pml_file.writelines(toWritetoFile)
        pml_file.close()
        print 'Ok'
        cubit.cmd('set info on') # Turn on return messages from Cubit commands
        cubit.cmd('set echo on') # Turn on echo of Cubit commands
    def nodescoord_write(self,nodecoord_name):
        """ Write nodes coordinates on file : nodecoord_name """
        nodecoord = open(nodecoord_name,'w')
        print 'Writing '+nodecoord_name+'.....'
        node_list = cubit.parse_cubit_list('node','all') # Import all the nodes of the model
        num_nodes = len(node_list) # Total number of nodes
        nodecoord.write('%10i\n' % num_nodes) # Write the number of nodes on the first line
        for node in node_list: # For all nodes
            x,y,z = cubit.get_nodal_coordinates(node) # Import its coordinates (3 coordinates even for a 2D model in cubit)
            txt = ('%20f %20f\n') % (x,z)
            nodecoord.write(txt) # Write x and z coordinates on the file -> Model must be in x,z coordinates. TODO
        nodecoord.close()
        print 'Ok'
    def free_write(self,freename): #freename = None):
        """ Write free surface on file : freename """
        cubit.cmd('set info off') # Turn off return messages from Cubit commands
        cubit.cmd('set echo off') # Turn off echo of Cubit commands
        cubit.cmd('set journal off') # Do not save journal file
        from sets import Set
        # if not freename: freename = self.freename
        freeedge = open(freename,'w')
        print 'Writing '+freename+'.....'
        if self.topo_mesh:
            for block,flag in zip(self.block_bc,self.block_bc_flag): # For each 1D block
                if block == self.topography: # If the block correspond to topography
                    edges_all = Set(cubit.get_block_edges(block)) # Import all topo edges id as a Set
            toWritetoFile = [] #[""]*(len(edges_all)+1)
            toWritetoFile.append('%10i\n' % len(edges_all)) # Print the number of edges on the free surface
            for block,flag in zip(self.block_mat,self.block_flag): # For each 2D block
                print block,flag
                quads = cubit.get_block_faces(block) # Import quads id
                for quad in quads: # For each quad
                    edges = Set(cubit.get_sub_elements("face", quad, 1)) # Get the lower dimension entities associated with a higher dimension entities.
                    # Here it gets the 1D edges associates with the face of id "quad". Store it as a Set
                    intersection = edges & edges_all # Contains the edges of the considered quad that is on the free surface
                    if len(intersection) != 0: # If this quad touch the free surface
                        #print "  ",quad," -> this quad touch the free surface!"
                        nodes = cubit.get_connectivity('face',quad) # Import the nodes describing the quad
                        #print "    it is described by nodes:",nodes," and edges :",edges
                        #print "      edges:",intersection," is/are on the free surface"
                        nodes = self.jac_check(list(nodes)) # Check the jacobian of the quad
                        for e in intersection: # For each edge on the free surface
                            node_edge = cubit.get_connectivity('edge',e) # Import the nodes describing the edge
                            #print "      edge",e,"is composed of nodes",node_edge
                            nodes_ok = []
                            for i in nodes: # Loop on the nodes of the quad
                                if i in node_edge: # If this node is belonging to the free surface
                                    nodes_ok.append(i) # Put it in nodes_ok
                            #print "    nodes:",nodes_ok,"belong to free surface"
                            txt = '%10i %10i %10i %10i\n' % (quad,2,nodes_ok[0],nodes_ok[1])
                            toWritetoFile.append(txt)
                            # Write the id of the quad, 2 (number of nodes describing a free surface elements), and the nodes
            freeedge.writelines(toWritetoFile)
        else:
            freeedge.write('0') # Even without any free surface specfem2d need a file with a 0 in first line
        freeedge.close()
        print 'Ok'
        cubit.cmd('set info on') # Turn on return messages from Cubit commands
        cubit.cmd('set echo on') # Turn on echo of Cubit commands
    def forcing_write(self,forcname):
        """ Write forcing surfaces on file : forcname """
        cubit.cmd('set info off') # Turn off return messages from Cubit commands
        cubit.cmd('set echo off') # Turn off echo of Cubit commands
        cubit.cmd('set journal off') # Do not save journal file
        from sets import Set
        forceedge = open(forcname,'w')
        print 'Writing '+forcname+'.....'
        edges_forc = [Set()]*self.nforc # edges_forc[0] will be a Set containing the nodes describing the forcing boundary
        # (index 0 : bottom, index 1 : right, index 2 : top, index 3 : left)
        nedges_all = 0 # To count the total number of forcing edges
        for block,flag in zip(self.block_bc,self.block_bc_flag): # For each 1D block
            for iforc in range(0, self.nforc): # iforc = 0,1,2,3 : for each forcing boundaries
                if block == self.forcing_boun[iforc]: # If the block considered correspond to the boundary
                    edges_forc[iforc] = Set(cubit.get_block_edges(block)) # Store each edge on edges_forc
                    nedges_all = nedges_all+len(edges_forc[iforc]) # add the number of edges to nedges_all
        toWritetoFile = [""]*(nedges_all+1)
        toWritetoFile[0] = '%10i\n' % nedges_all # Write the total number of forcing edges to the first line of file
        #forceedge.write('%10i\n' % nedges_all) # Write the total number of forcing edges to the first line of file
        print 'Number of edges', nedges_all
        #id_element = 0
        indexFile = 1
        for block,flag in zip(self.block_mat,self.block_flag): # For each 2D block
                quads = cubit.get_block_faces(block) # Import quads id
                for quad in quads: # For each quad
                    #id_element = id_element+1 # id of this quad
                    edges = Set(cubit.get_sub_elements("face", quad, 1)) # Get the lower dimension entities associated with a higher dimension entities.
                    # Here it gets the 1D edges associates with the face of id "quad". Store it as a Set
                    for iforc in range(0,self.nforc): # iforc = 0,1,2,3 : for each forcing boundaries
                        intersection = edges & edges_forc[iforc]  # Contains the edges of the considered quad that is on the forcing boundary considered
                        if len(intersection) != 0: # If this quad touch the forcing boundary considered
                            nodes = cubit.get_connectivity('face',quad) # Import the nodes describing the quad
                            nodes = self.jac_check(list(nodes)) # Check the jacobian of the quad
                            for e in intersection: # For each edge on the forcing boundary considered
                                node_edge = cubit.get_connectivity('edge',e) # Import the nodes describing the edge
                                nodes_ok = []
                                for i in nodes: # Loop on the nodes of the quad
                                    if i in node_edge: # If this node is belonging to forcing surface
                                        nodes_ok.append(i) # add it to nodes_ok
                                # forcname contains 1/ element number, 2/ number of nodes that form the acoustic forcing edge
                                # (which currently must always be equal to two, see comment below),
                                # 3/ first node on the acforcing surface, 4/ second node on the acforcing surface
                                # 5/ 1 = IBOTTOME, 2 = IRIGHT, 3 = ITOP, 4 = ILEFT
                                #txt = '%10i %10i %10i %10i %10i\n' % (id_element,2,nodes_ok[0],nodes_ok[1],iforc+1)
                                txt = '%10i %10i %10i %10i %10i\n' % (quad,2,nodes_ok[0],nodes_ok[1],iforc+1)
                                # Write the id of the quad, 2 (number of nodes describing a free surface elements), the nodes and the type of boundary
                                #print indexFile
                                toWritetoFile[indexFile] = txt
                                indexFile = indexFile + 1
                                #forceedge.write(txt)
        forceedge.writelines(toWritetoFile)
        forceedge.close()
        print 'Ok'
        cubit.cmd('set info on') # Turn on return messages from Cubit commands
        cubit.cmd('set echo on') # Turn on echo of Cubit commands
    def abs_write(self,absname):
        """ Write absorbing surfaces on file : absname """
        cubit.cmd('set info off') # Turn off return messages from Cubit commands
        cubit.cmd('set echo off') # Turn off echo of Cubit commands
        cubit.cmd('set journal off') # Do not save journal file.
        from sets import Set
        # if not absname: absname = self.absname
        absedge = open(absname,'w')
        print 'Writing '+absname+'.....'
        edges_abs = [Set()]*self.nabs # edges_abs[0] will be a Set containing the nodes describing bottom adsorbing boundary
        # (index 0 : bottom, index 1 : right, index 2 : top, index 3 : left)
        nedges_all = 0 # To count the total number of absorbing edges
        for block,flag in zip(self.block_bc,self.block_bc_flag): # For each 1D block
            for iabs in range(0, self.nabs): # iabs = 0,1,2,3 : for each absorbing boundaries
                if block == self.abs_boun[iabs]: # If the block considered correspond to the boundary
                    edges_abs[iabs] = Set(cubit.get_block_edges(block)) # Store each edge on edges_abs
                    nedges_all = nedges_all+len(edges_abs[iabs]); # add the number of edges to nedges_all
        toWritetoFile = [""]*(nedges_all+1)
        toWritetoFile[0] = '%10i\n' % nedges_all # Write the total number of absorbing edges to the first line of file
        #absedge.write('%10i\n' % nedges_all) # Write the total number of absorbing edges to the first line of file
        print 'Number of edges', nedges_all
        #id_element = 0
        indexFile = 1
        for block,flag in zip(self.block_mat,self.block_flag): # For each 2D block
                quads = cubit.get_block_faces(block) # Import quads id
                for quad in quads: # For each quad
                    #id_element = id_element+1 # id of this quad
                    edges = Set(cubit.get_sub_elements("face", quad, 1)) # Get the lower dimension entities associated with a higher dimension entities.
                    # Here it gets the 1D edges associates with the face of id "quad". Store it as a Set
                    for iabs in range(0,self.nabs): # iabs = 0,1,2,3 : for each absorbing boundaries
                        intersection = edges & edges_abs[iabs]  # Contains the edges of the considered quad that is on the absorbing boundary considered
                        if len(intersection) != 0: # If this quad touch the absorbing boundary considered
                            nodes = cubit.get_connectivity('face',quad) # Import the nodes describing the quad
                            nodes = self.jac_check(list(nodes)) # Check the jacobian of the quad
                            for e in intersection: # For each edge on the absorbing boundary considered
                                node_edge = cubit.get_connectivity('edge',e) # Import the nodes describing the edge
                                nodes_ok = []
                                for i in nodes: # Loop on the nodes of the quad
                                    if i in node_edge: # If this node is belonging to absorbing surface
                                        nodes_ok.append(i) # Add it to nodes_ok
                                #txt = '%10i %10i %10i %10i %10i\n' % (id_element,2,nodes_ok[0],nodes_ok[1],iabs+1)
                                txt = '%10i %10i %10i %10i %10i\n' % (quad,2,nodes_ok[0],nodes_ok[1],iabs+1)
                                # Write the id of the quad, 2 (number of nodes describing a free surface elements), the nodes and the type of boundary
                                toWritetoFile[indexFile] = txt
                                indexFile = indexFile + 1
                                #absedge.write(txt)
        absedge.writelines(toWritetoFile)
        absedge.close()
        print 'Ok'
        cubit.cmd('set info on') # Turn on return messages from Cubit commands
        cubit.cmd('set echo on') # Turn on echo of Cubit commands
    def axis_write(self,axis_name):
        """ Write axis on file """
        cubit.cmd('set info off') # Turn off return messages from Cubit commands
        cubit.cmd('set echo off') # Turn off echo of Cubit commands
        cubit.cmd('set journal off') # Do not save journal file
        from sets import Set
        axisedge = open(axis_name,'w')
        print 'Writing '+axis_name+'.....'
        for block,flag in zip(self.block_bc,self.block_bc_flag): # For each 1D block
            if block == self.axisId: # If the block correspond to the axis
                edges_all = Set(cubit.get_block_edges(block)) # Import all axis edges id as a Set
        toWritetoFile = [""]*(len(edges_all)+1)
        toWritetoFile[0] = '%10i\n' % len(edges_all) # Write the number of edges on the axis
        #axisedge.write('%10i\n' % len(edges_all)) # Write the number of edges on the axis
        print 'Number of edges on the axis :',len(edges_all)
        #id_element = 0
        indexFile = 1
        for block,flag in zip(self.block_mat,self.block_flag): # For each 2D block
                quads = cubit.get_block_faces(block) # Import quads id
                for quad in quads: # For each quad
                    #id_element = id_element+1 # id of this quad
                    edges = Set(cubit.get_sub_elements("face", quad, 1)) # Get the lower dimension entities associated with a higher dimension entities.
                    # Here it gets the 1D edges associates with the face of id "quad". Store it as a Set
                    intersection = edges & edges_all # Contains the edges of the considered quad that are on the axis
                    if len(intersection) != 0: # If this quad touch the axis
                        nodes = cubit.get_connectivity('face',quad) # Import the nodes describing the quad
                        nodes = self.jac_check(list(nodes)) # Check the jacobian of the quad
                        for e in intersection: # For each edge on the axis
                            node_edge = cubit.get_connectivity('edge',e) # Import the nodes describing the edge
                            nodes_ok = []
                            for i in nodes: # Loop on the nodes of the quad
                                if i in node_edge: # If this node is belonging to the axis
                                    nodes_ok.append(i) # Add it to nodes_ok
                            txt = '%10i %10i %10i %10i\n' % (quad,2,nodes_ok[0],nodes_ok[1])
                            #txt = '%10i %10i %10i %10i\n' % (id_element,2,nodes_ok[0],nodes_ok[1])
                            # Write the id of the quad, 2 (number of nodes describing a free surface elements), the nodes
                            toWritetoFile[indexFile] = txt
                            indexFile = indexFile + 1
                            #axisedge.write(txt)
        axisedge.writelines(toWritetoFile)
        axisedge.close()
        print 'Ok'
        cubit.cmd('set info on') # Turn on return messages from Cubit commands
        cubit.cmd('set echo on') # Turn on echo of Cubit commands
    def rec_write(self,recname):
        """ Write receivers coordinates on file recname """
        print 'Writing '+self.recname+'.....'
        recfile = open(self.recname,'w')
        nodes = cubit.get_nodeset_nodes(self.receivers) # Import nodes in nodeset containing receiver positions
        for i,n in enumerate(nodes): # For each receiver
            x,y,z = cubit.get_nodal_coordinates(n) # Import its coordinates (3 coordinates even for a 2D model in cubit)
            recfile.write('ST%i XX %20f %20f 0.0 0.0 \n' % (i,x,z))  # Write x and z coordinates on the file -> Model must be in x,z coordinates. TODO
        recfile.close()
        print 'Ok'
    def write(self,path = ''):
        """ Write mesh in specfem2d format """
        print 'Writing '+self.recname+'.....'
        import os
        cubit.cmd('set info off') # Turn off return messages from Cubit commands
        cubit.cmd('set echo off') # Turn off echo of Cubit commands
        cubit.cmd('set journal off') # Do not save journal file
        if len(path) != 0: # If a path is supplied add a / at the end if needed
            if path[-1] != '/': path = path+'/'
        else:
            path = os.getcwd()+'/'
        self.mesh_write(path+self.mesh_name) # Write mesh file
        self.material_write(path+self.material_name) # Write material file
        self.nodescoord_write(path+self.nodecoord_name) # Write nodes coord file
        self.free_write(path+self.freename) # Write free surface file (specfem2d needs it even if there is no free surface)
        if self.abs_mesh:
            self.abs_write(path+self.absname) # Write absorbing surface file
        if self.forcing_mesh:
            self.forcing_write(path+self.forcname) # Write forcing surface file
        if self.axisymmetric_mesh:
            self.axis_write(path+self.axisname) # Write axis on file
        if self.pml_layers:
            self.pmls_write(path+self.pmlname) # Write axis on file
        if self.write_nummaterial_velocity_file:
            self.nummaterial_write(path+self.nummaterial_name) # Write nummaterial file
        if self.receivers:
            self.rec_write(path+self.recname) # If receivers has been set (as nodeset) write receiver file as well
        print 'Mesh files has been writen in '+path
        cubit.cmd('set info on') # Turn on return messages from Cubit commands
        cubit.cmd('set echo on') # Turn on echo of Cubit commands

profile = mesh() # Store the mesh from Cubit
profile.write() # Write it into files (in specfem2d format). profile.write(/path/to/directory)
