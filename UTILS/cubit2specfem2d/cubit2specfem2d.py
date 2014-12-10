#!python
#!/usr/bin/env python

class mtools(object):
    """docstring for mtools"""
    def __init__(self,frequency,list_surf,list_vp):
        super(mtools, self).__init__()
        self.frequency = frequency
        self.list_surf = list_surf
        self.list_vp = list_vp
        self.ngll=5
        self.percent_gll=0.172
        self.point_wavelength=5
    def __repr__(self):
        txt='Meshing for frequency up to '+str(self.frequency)+'Hz\n'
        for surf,vp in zip(self.list_surf,self.list_vp):
            txt=txt+'surface '+str(surf)+', vp ='+str(vp)+'  -> size '+str(self.freq2meshsize(vp)[0])+' -> dt '+str(self.freq2meshsize(vp)[0])+'\n' 
        return txt
    def freq2meshsize(self,vp):
        velocity=vp*.5
        self.size=(1/2.5)*velocity/self.frequency*(self.ngll-1)/self.point_wavelength
        self.dt=.4*self.size/vp*self.percent_gll
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
    def create_blocks(self,mesh_entity,list_entity=None,):
        if mesh_entity =='surface':
            txt=' face in surface '
        elif mesh_entity == 'curve':
            txt=' edge in curve '
        elif mesh_entity == 'group':
            txt=' face in group '
        if list_entity:
            if not isinstance(list_entity,list):
                list_entity=[list_entity]
        for entity in list_entity:
            iblock=cubit.get_next_block_id()
            command = "block "+str(iblock)+ txt +str(entity)
            cubit.cmd(command)
    def material_file(self,filename):
        matfile=open(filename,'w')
        material=[]
        for record in matfile:
            mat_name,vp_str=record.split()
            vp=float(vp_str)
            material.append([mat_name,vp])
        self.material=dict(material)
    def assign_block_material(self,id_block,mat_name,vp=None):
        try:
            material=self.material
        except:
            material=None
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
    dt,edge_dt,freq,edge_freq=seismic_resolution(edges,velocity,bins_d=None,bins_u=None,sidelist=None,ngll=5,np=8)
        Given the velocity of a list of edges, seismic_resolution provides the minimum Dt required for the stability condition (and the corrisponding edge).
        Furthermore, given the number of gll point in the element (ngll) and the number of GLL point for wavelength, it provide the maximum resolved frequency.
    #########
    length=edge_length(edge)
        return the length of a edge
    #########
    edge_min,length=edge_min_length(surface)
        given the cubit id of a surface, it return the edge with minimun length 
    #########
    """
    def __int__(self):
        pass
    def seismic_resolution(self,edges,velocity,bins_d=None,bins_u=None,sidelist=None):
        """
        dt,edge_dt,freq,edge_freq=seismic_resolution(edges,velocity,bins_d=None,bins_u=None,sidelist=None,ngll=5,np=8)
            Given the velocity of a list of edges, seismic_resolution provides the minimum Dt required for the stability condition (and the corrisponding edge).
            Furthermore, given the number of gll point in the element (ngll) and the number of GLL point for wavelength, it provide the maximum resolved frequency.
        """
        ratiostore=1e10
        dtstore=1e10
        edgedtstore=-1
        edgeratiostore=-1
        for edge in edges:
            d=self.edge_length(edge)
            ratio=(1/2.5)*velocity/d*(self.ngll-1)/self.point_wavelength
            dt=.4*d/velocity*self.percent_gll
            if dt<dtstore:
               dtstore=dt
               edgedtstore=edge
            if ratio < ratiostore:
                ratiostore=ratio
                edgeratiostore=edge
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
        length=edge_length(edge)
            return the length of a edge
        """
        from math import sqrt
        nodes=cubit.get_connectivity('edge',edge)
        x0,y0,z0=cubit.get_nodal_coordinates(nodes[0])
        x1,y1,z1=cubit.get_nodal_coordinates(nodes[1])
        d=sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
        return d
    def edge_min_length(self,surface):
        """
        edge_min,length=edge_min_length(surface)
            given the cubit id of a surface, it return the edge with minimun length 
        """
        from math import sqrt
        self.dmin=99999
        edge_store=0
        command = "group 'list_edge' add edge in surf "+str(surface)
        command = command.replace("["," ").replace("]"," ")
        #print command
        cubit.cmd(command)
        group=cubit.get_id_from_name("list_edge")
        edges=cubit.get_group_edges(group)
        command = "delete group "+ str(group)
        cubit.cmd(command)
        for edge in edges:
            d=self.edge_length(edge)
            if d<dmin:
                self.dmin=d
                edge_store=edge
        self.edgemin=edge_store
        return self.edgemin,self.dmin
    def jac_check(self,nodes):
        x0=cubit.get_nodal_coordinates(nodes[0])
        x1=cubit.get_nodal_coordinates(nodes[1])
        x2=cubit.get_nodal_coordinates(nodes[2])
        xv1=x1[0]-x0[0]
        xv2=x2[0]-x1[0]
        zv1=x1[2]-x0[2]
        zv2=x2[2]-x1[2]
        jac=-xv2*zv1+xv1*zv2
        if  jac > 0:
            return nodes
        elif jac < 0:
            return nodes[0],nodes[3],nodes[2],nodes[1]
        else:
            print 'error, jacobian=0', jac,nodes 
    def mesh_analysis(self,frequency):
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        bins_d=[0.0001]+range(0,int(frequency)+1)+[1000]
        bins_u=bins_d[1:]
        dt=[]
        ed_dt=[]
        r=[]
        ed_r=[]
        nstart=cubit.get_next_sideset_id()
        command = "del sideset all"
        cubit.cmd(command)
        for bin_d,bin_u in zip(bins_d,bins_u):
            nsideset=cubit.get_next_sideset_id()
            command='create sideset '+str(nsideset)
            cubit.cmd(command)
            command = "sideset "+str(nsideset)+ " name "+ "'ratio-["+str(bin_d)+"_"+str(bin_u)+"['"
            cubit.cmd(command)
        nend=cubit.get_next_sideset_id()            
        sidelist=range(nstart,nend)
        for block in self.block_mat:
            name=cubit.get_exodus_entity_name('block',block)
            velocity=self.material[name][1]
            if velocity > 0:
                faces=cubit.get_block_faces(block)
                edges=[]
                for face in faces:
                    es=cubit.get_sub_elements("face", face, 1)
                    edges=edges+list(es)
                dtstore,edgedtstore,ratiostore,edgeratiostore=self.seismic_resolution(edges,velocity,bins_d,bins_u,sidelist)
                dt.append(dtstore)
                ed_dt.append(edgedtstore)
                r.append(ratiostore)
                ed_r.append(edgeratiostore)
        self.ddt=zip(ed_dt,dt)
        self.dr=zip(ed_r,r)
        def sorter(x, y):
            return cmp(x[1],y[1])
        self.ddt.sort(sorter)
        self.dr.sort(sorter)
        print self.ddt,self.dr
        print 'Deltat minimum => edge:'+str(self.ddt[0][0])+' dt: '+str(self.ddt[0][1])
        print 'minimum frequency resolved => edge:'+str(self.dr[0][0])+' frequency: '+str(self.dr[0][1])
        return self.ddt[0],self.dr[0]

class mesh(object,mesh_tools):
    def __init__(self):
        super(mesh, self).__init__()
        self.mesh_name='mesh_file'
        self.nodecoord_name='nodes_coords_file'
        self.material_name='materials_file'
        self.nummaterial_name='nummaterial_velocity_file'
        self.absname='absorbing_surface_file'
        self.freename='free_surface_file'
        self.recname='STATIONS'
        self.face='QUAD4'
        self.edge='BAR2'
        self.topo='topo'
        self.abs_boun_name=['abs_bottom','abs_right','abs_top','abs_left']
        self.abs_boun=[]  # block numbers for abs bouns
        self.nabs=4
        self.rec='receivers'
        self.block_definition() # export blocks from cubit
        self.ngll=5
        self.percent_gll=0.172
        self.point_wavelength=5
        cubit.cmd('compress')
    def __repr__(self):
        pass
    def block_definition(self):
        block_flag=[]
        block_mat=[]
        block_bc=[]
        block_bc_flag=[]
        abs_boun=[-1] * self.nabs # total 4 sides of absorbing bouns
        material={}
        bc={}
        blocks=cubit.get_block_id_list()
        for block in blocks:
            name=cubit.get_exodus_entity_name('block',block)
            ty=cubit.get_block_element_type(block)
            if ty == self.face:
                flag=int(cubit.get_block_attribute_value(block,0))
                vel=cubit.get_block_attribute_value(block,1)
                block_flag.append(int(flag))
                block_mat.append(block)
                par=tuple([flag,vel])
                material[name]=par
            elif ty == self.edge:
                block_bc_flag.append(2)
                block_bc.append(block)
                bc[name]=2 #edge has connectivity = 2
                if name == self.topo: topography=block
                if name in self.abs_boun_name : 
                    abs_boun[self.abs_boun_name.index(name)]=block
            else:
                print 'blocks not properly defined', ty
                return None, None,None,None,None,None,None,None
        nsets=cubit.get_nodeset_id_list()
        if len(nsets) == 0: self.receivers=None
        for nset in nsets:
            name=cubit.get_exodus_entity_name('nodeset',nset)
            if name == self.rec:
                self.receivers=nset
            else:
                print 'nodeset '+name+' not defined'
                self.receivers=None
        try:
            self.block_mat=block_mat
            self.block_flag=block_flag
            self.block_bc=block_bc
            self.block_bc_flag=block_bc_flag
            self.material=material
            self.bc=bc
            self.topography=topography
            self.abs_boun=abs_boun
        except:
            print 'blocks not properly defined'
        print 'abs_boun=',abs_boun
    def tomo(self,flag,vel):
        vp=vel/1000
        rho=(1.6612*vp-0.472*vp**2+0.0671*vp**3-0.0043*vp**4+0.000106*vp**4)*1000
        txt='%3i %1i %20f %20f %20f %1i %1i\n' % (flag,1,rho,vel,vel/(3**.5),0,0)
        return txt
    def nummaterial_write(self,nummaterial_name):
        print 'Writing '+nummaterial_name+'.....'
        nummaterial=open(nummaterial_name,'w')
        for block in self.block_mat:
            name=cubit.get_exodus_entity_name('block',block)
            nummaterial.write(self.tomo(self.material[name][0],self.material[name][1]))
        nummaterial.close()
        print 'Ok'
    def mesh_write(self,mesh_name):
        meshfile=open(mesh_name,'w')
        print 'Writing '+mesh_name+'.....'
        num_elems=cubit.get_quad_count()
        meshfile.write(str(num_elems)+'\n')
        num_write=0
        for block,flag in zip(self.block_mat,self.block_flag):
            quads=cubit.get_block_faces(block)
            for inum,quad in enumerate(quads):
                nodes=cubit.get_connectivity('face',quad)
                nodes=self.jac_check(nodes)
                txt=('%10i %10i %10i %10i\n')% nodes
                meshfile.write(txt)
            num_write=num_write+inum+1
            print 'block', block, 'number of quad4', inum
        meshfile.close()
        print 'Ok num elements/write=',str(num_elems), str(num_write)
    def material_write(self,mat_name):
        mat=open(mat_name,'w')
        print 'Writing '+mat_name+'.....'
        for block,flag in zip(self.block_mat,self.block_flag):
                quads=cubit.get_block_faces(block)
                for quad in quads:
                    mat.write(('%10i\n') % flag)
        mat.close()
        print 'Ok'
    def nodescoord_write(self,nodecoord_name):
        nodecoord=open(nodecoord_name,'w')
        print 'Writing '+nodecoord_name+'.....'
        node_list=cubit.parse_cubit_list('node','all')
        num_nodes=len(node_list)
        nodecoord.write('%10i\n' % num_nodes)
        #
        for node in node_list:
            x,y,z=cubit.get_nodal_coordinates(node)
            txt=('%20f %20f\n') % (x,z)
            nodecoord.write(txt)
        nodecoord.close()
        print 'Ok'
    def free_write(self,freename=None):
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        from sets import Set
        if not freename: freename=self.freename
        freeedge=open(freename,'w')
        print 'Writing '+freename+'.....'
        #
        #
        for block,flag in zip(self.block_bc,self.block_bc_flag):
            if block == self.topography:
                edges_all=Set(cubit.get_block_edges(block))
        freeedge.write('%10i\n' % len(edges_all))
        print len(edges_all)
        id_element=0
        for block,flag in zip(self.block_mat,self.block_flag):
                quads=cubit.get_block_faces(block)
                for quad in quads:
                    id_element=id_element+1
                    edges=Set(cubit.get_sub_elements("face", quad, 1))
                    intersection=edges & edges_all
                    if len(intersection) != 0:
                        for e in intersection:
                            node_edge=cubit.get_connectivity('edge',e)
                            nodes=cubit.get_connectivity('face',quad)
                            nodes=self.jac_check(list(nodes))
                            nodes_ok=[]
                            for i in nodes:
                                if i in node_edge:
                                    nodes_ok.append(i)
                            txt='%10i %10i %10i %10i\n' % (id_element,2,nodes_ok[0],nodes_ok[1])
                            freeedge.write(txt)
        freeedge.close()
#        print edges_all
        print 'Ok'
        cubit.cmd('set info on')
        cubit.cmd('set echo on')
    def abs_write(self,absname=None):
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        from sets import Set
        if not absname: absname=self.absname
        absedge=open(absname,'w')
        print 'Writing '+absname+'.....'
        #
        edges_abs=[Set()]*self.nabs
        nedges_all=0
        for block,flag in zip(self.block_bc,self.block_bc_flag):
            for iabs in range(0, self.nabs):
                if block == self.abs_boun[iabs]:
                    edges_abs[iabs]=Set(cubit.get_block_edges(block))
                    nedges_all=nedges_all+len(edges_abs[iabs]);
        absedge.write('%10i\n' % nedges_all)
        print 'number of edges', nedges_all
        id_element=0
        for block,flag in zip(self.block_mat,self.block_flag):
                quads=cubit.get_block_faces(block)
                for quad in quads:
                    id_element=id_element+1
                    edges=Set(cubit.get_sub_elements("face", quad, 1))
                    for iabs in range(0,self.nabs):
                        intersection=edges & edges_abs[iabs]
                        if len(intersection) != 0:
                            #print intersection, iabs
                            for e in intersection:
                                node_edge=cubit.get_connectivity('edge',e)
                                nodes=cubit.get_connectivity('face',quad)
                                nodes=self.jac_check(list(nodes))
                                nodes_ok=[]
                                for i in nodes:
                                    if i in node_edge:
                                        nodes_ok.append(i)
                                txt='%10i %10i %10i %10i %10i\n' % (id_element,2,nodes_ok[0],nodes_ok[1],iabs+1)
                                absedge.write(txt)
        absedge.close()
#        print edges_abs[0:self.nabs],'- 2 corners'
        print 'Ok'
        cubit.cmd('set info on')
        cubit.cmd('set echo on')
    def rec_write(self,recname):
        print 'Writing '+self.recname+'.....'
        recfile=open(self.recname,'w')
        nodes=cubit.get_nodeset_nodes(self.receivers)
        for i,n in enumerate(nodes):
            x,y,z=cubit.get_nodal_coordinates(n)
            recfile.write('ST%i XX %20f %20f 0.0 0.0 \n' % (i,x,z))
        recfile.close()
        print 'Ok'
    def write(self,path=''):
        cubit.cmd('set info off')
        cubit.cmd('set echo off')
        cubit.cmd('set journal off')
        if len(path) != 0:
            if path[-1] != '/': path=path+'/'
        self.mesh_write(path+self.mesh_name)
        self.material_write(path+self.material_name)
        self.nodescoord_write(path+self.nodecoord_name)
        self.free_write(path+self.freename)
        self.abs_write(path+self.absname)
        self.nummaterial_write(path+self.nummaterial_name)
        if self.receivers: self.rec_write(path+self.recname)
        cubit.cmd('set info on')
        cubit.cmd('set echo on')

profile=mesh()
profile.write()
