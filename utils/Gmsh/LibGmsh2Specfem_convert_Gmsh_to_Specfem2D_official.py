#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Python code to link gmsh with specfem 
#  New version with identification of the zone
#  and handling of PML 
#
#@author: Cristini Paul, 
#  Laboratoire de Mecanique et d'Acoustique, CNRS, Marseille, France
#
# September 2012, 6rd
#
import sys, string
from os.path import splitext

try:
    from numpy import *
except ImportError:
    print "error: package python-numpy is not installed"
#
def SauvFicSpecfem(Ng, Ct, Var, Fv):
    # Sauvegarde au format ascii
    # Ng est le nom generique
    # Ct le nombre de lignes a lire
    # Var est le nom de la variable contenant les informations a ecrire
    # Fv est le format d'ecriture '%f' pour les noeuds, '%i' pour les autres fichiers
    savetxt(Ng,(Ct,), fmt='%i') # on met le nombre de lignes en entete
    fd = open(Ng,'a')
    savetxt(fd, Var, fmt=Fv)
    fd.close()
    return
#
def OuvreGmsh(Dir,Nom,Bords):
    # Reading file .msh created with Gmsh
    if splitext(Nom)[-1]=='.msh':
       fic=Nom
       Nom=splitext(Nom)[0]
    elif splitext(Nom)[-1]=='':
       fic=Nom+'.msh'
    else:
        print 'File extension is not correct'
        print 'script aborted'
        sys.exit()
    print '#'*20+' STARTING PROCESSING '+'#'*20
    print 'File : ', fic
    print '-'*60
    #
    # Open the file and get the lines
    # 
    f = file(Dir+fic,'rU')
    lignes= f.readlines()
    f.close()
    # Looking for positions
    #
    for ii in range(len(lignes)):
        if lignes[ii]=='$Nodes\n': PosNodes=ii 
        if lignes[ii]=='$PhysicalNames\n': PosPhys=ii 
        if lignes[ii]=='$Elements\n':
            PosElem=ii
            break
    # Elements type 4 nodes or 9 nodes
    TypElem1D = int(string.split(lignes[PosElem+2])[1])
    if TypElem1D==1:
        Ngnod, LinElem, SurfElem = 4, 1, 3
        len1D, len2D = 2, 4
    elif TypElem1D==8:
        Ngnod, LinElem, SurfElem = 9, 8, 10
        len1D, len2D = 3, 9
    else:
        print 'The number of nodes of elemnts is not 4 nor 9'
    #-------------------------------------------------------------------------
    # Conditions on sides of the domain
    # Possible choices: Abso, Free
    Bord_abso, Bord_free = [], []  # Initialization
    PML = False
    print Bords
    #-------------------------------------------------------------------------
    # PHYSICAL NAMES
    NbPhysNames = int(string.split(lignes[PosPhys+1])[0])
    # Detection of the presence of PML
    for Ip in range(NbPhysNames):
        Nam = string.split(lignes[PosPhys+2+Ip])[2][1:-1]
        if Nam.endswith('PML'):
            PML = True
            break
    print 'Number of physical Names : ', NbPhysNames
    print 'PML : ', PML
    # Structure de la variable : 1 entier et une chaine de caractere de taille 16
    dt = dtype([('dimension',int), ('zone', int), ('name', str, 16)])
    PhysCar=zeros((NbPhysNames,), dtype=dt)
    Med=zeros(NbPhysNames+1)   # Les milieux: domaines dont le nom commence par M
    #--------------------------------------------------------------------------
    # Initialisation des zones (celles qui ne seront pas affectees auront 
    #  une valeur nulle)
    PML_right, PML_right_bottom, PML_right_top = 0, 0, 0
    PML_left, PML_left_bottom, PML_left_top = 0, 0, 0
    PML_top, PML_bottom = 0, 0
    
    for Ip in range(NbPhysNames):
        Dim = int(string.split(lignes[PosPhys+2+Ip])[0])
        Zon = int(string.split(lignes[PosPhys+2+Ip])[1])
        Nam = string.split(lignes[PosPhys+2+Ip])[2][1:-1]
        PhysCar[Ip] = (Dim, Zon, Nam)
        print PhysCar[Ip]
        if Nam.startswith('M'): Med[Zon]=int(Nam[1:])   # Milieux
        if Bords.has_key(Nam):
            if Bords[Nam] == 'Abso': Bord_abso.append(Zon)
            if Bords[Nam] == 'Free': Bord_free.append(Zon)
        if PML:
            #------------------------------------------------------------------
            # Les bords du domaine PML
            if Nam == 'Right_PML' : Bord_right=Zon
            if Nam == 'Left_PML'  : Bord_left=Zon
            if Nam == 'Top_PML'   : Bord_top=Zon
            if Nam == 'Bottom_PML': Bord_bottom=Zon
        else:
            #------------------------------------------------------------------
            # Les bords du domaine initial
            if Nam == 'Right' : Bord_right=Zon
            if Nam == 'Left'  : Bord_left=Zon
            if Nam == 'Top'   : Bord_top=Zon
            if Nam == 'Bottom': Bord_bottom=Zon
        if PML:
            #------------------------------------------------------------------
            #  Les differents domaines PML
            if Nam == 'R'     : PML_right, Med[Zon] = Zon, 1
            if Nam == 'RB'    : PML_right_bottom, Med[Zon] = Zon, 1    # Corner
            if Nam == 'L'     : PML_left, Med[Zon] = Zon, 1
            if Nam == 'LB'    : PML_left_bottom, Med[Zon] = Zon, 1     # Corner
            if Nam == 'T'     : PML_top, Med[Zon] = Zon, 1
            if Nam == 'RT'    : PML_right_top, Med[Zon] = Zon, 1       # Corner
            if Nam == 'LT'    : PML_left_top, Med[Zon] = Zon, 1        # Corner
            if Nam == 'B'     : PML_bottom, Med[Zon] = Zon, 1
    #--------------------------------------------------------------------------
    print '-'*60
    print 'Absorbing boundaries : ', Bord_abso
    print 'Free boundaries : ', Bord_free
    print '-'*60    
    print 'Right boundaries : ', Bord_right
    print 'Left boundaries : ', Bord_left
    print 'Top boundaries : ', Bord_top
    print 'Bottom boundaries : ', Bord_bottom
    print '-'*60
    print 'Med : ',Med
    #---------------------------------------------------------------------------
    # Infos sur le fichier Gmsh
    Ver=float(string.split(lignes[1])[0])
    File_Type=int(string.split(lignes[1])[1])
    Data_Size=int(string.split(lignes[1])[2])
    # Nodes
    NbNodes=int(string.split(lignes[PosNodes+1])[0])
    print 'Number of nodes: ',NbNodes
    Nodes=zeros((NbNodes,2),dtype=float)
    for Ninc in range(NbNodes):
        Nodes[Ninc]= [float(val) for val in \
                      (string.split(lignes[PosNodes+2+Ninc])[1:3])]
    #
    # Save to SPECFEM format
    SauvFicSpecfem('Nodes_'+Nom, NbNodes, Nodes, '%f')
    # Elements
    DecElem=12+NbNodes
    NbElements=int(string.split(lignes[PosElem+1])[0])
    print 'Number of elements: ', NbElements 
    #--------------------------------------------------------------------------
    #      Initialization
    #
    Elements           = empty((NbElements,len2D),dtype=int)
    Milieu             = empty((NbElements,1),dtype=int)
    Elements1D         = empty((NbElements,len1D),dtype=int)
    #--------------------------------------------------------------------------
    Elements1DBordAbso  = empty((NbElements,len1D),dtype=int)
    Elements1DBordFree  = empty((NbElements,len1D),dtype=int)
    Elements2DBordAbso  = zeros((NbElements,5),dtype=int) # only 2 nodes + zone
    Elements2DBordFree  = zeros((NbElements,4),dtype=int)
    #--------------------------------------------------------------------------
    Elements1DBordTop     = empty((NbElements,len1D),dtype=int)
    Elements1DBordBottom  = empty((NbElements,len1D),dtype=int)
    #--------------------------------------------------------------------------
    Elements1DBordRight  = empty((NbElements,len1D),dtype=int)
    Elements1DBordLeft   = empty((NbElements,len1D),dtype=int)
    #--------------------------------------------------------------------------
    if PML:
        ElementsPML      = empty((NbElements,2),dtype=int)  # elements + codage
    DecElem+=1  #
    #--------------------------------------------------------------------------
    # Initialization of counters
    N1D, N2D, N2DPML, N1DBord, N1DBordAbso, N1DBordFree = 0, 0, 0, 0, 0 ,0
    N1DBordLeft, N1DBordRight, N1DBordTop, N1DBordBottom = 0, 0, 0, 0
    #--------------------------------------------------------------------------
    # Loop over all elements
    # 
    for Ninc in range(NbElements):
        Pos     = PosElem+Ninc+2
        Ispec   = int(string.split(lignes[Pos])[0])
        TypElem = int(string.split(lignes[Pos])[1])
        ZonP    = int(string.split(lignes[Pos])[3])
        #
        if TypElem==LinElem: # Elements 1D (lines)
            Elements1D[N1D] = [int(val) for val in \
                               (string.split(lignes[Pos])[5:])] 
            #  Bottom
            if ZonP==Bord_bottom:
                Elements1DBordBottom[N1DBordBottom] = Elements1D[N1D]
                N1DBordBottom+=1
            # Top
            if ZonP==Bord_top:
                Elements1DBordTop[N1DBordTop] = Elements1D[N1D]
                N1DBordTop+=1
            # Left
            if ZonP==Bord_left:
                Elements1DBordLeft[N1DBordLeft] = Elements1D[N1D]
                N1DBordLeft+=1
            # Right
            if ZonP==Bord_right:
                Elements1DBordRight[N1DBordRight] = Elements1D[N1D]
                N1DBordRight+=1
            N1D+=1
        
        #----------------------------------------------------------------------

#
# flags for CPML absorbing boundaries:
#  CPML_left = 1
#  CPML_right = 2
#  CPML_bottom = 3
#  CPML_top = 4
#  CPML_top_left = 5
#  CPML_top_right = 6
#  CPML_bottom_left = 7
#  CPML_bottom_right = 8
#

        if TypElem==SurfElem:
            Elements[N2D]= [int(val) for val in \
                            (string.split(lignes[Pos])[5:])]
            Milieu[N2D]= Med[ZonP]
            if PML:
                # PML Bottom
                if ZonP==PML_bottom:
                    ElementsPML[N2DPML,:] = [N2D+1, 3]
                    N2DPML+=1
                # PML Right
                if ZonP==PML_right:
                    ElementsPML[N2DPML,:] = [N2D+1, 2]
                    N2DPML+=1
                # PML Top
                if ZonP==PML_top:
                    ElementsPML[N2DPML,:] = [N2D+1, 4]
                    N2DPML+=1
                # PML Left
                if ZonP==PML_left:
                    ElementsPML[N2DPML,:] = [N2D+1, 1]
                    N2DPML+=1
                #--------------------------------------------------------------
                # PML Corner Right Bottom
                if ZonP==PML_right_bottom:
                    ElementsPML[N2DPML,:] = [N2D+1, 8]
                    N2DPML+=1
                # PML Corner Right Top
                if ZonP==PML_right_top:
                    ElementsPML[N2DPML,:] = [N2D+1, 6]
                    N2DPML+=1
                # PML Corner Left Top
                if ZonP==PML_left_top:
                    ElementsPML[N2DPML,:] = [N2D+1, 5]
                    N2DPML+=1
                # PML Corner Left Bottom
                if ZonP==PML_left_bottom:
                    ElementsPML[N2DPML,:] = [N2D+1, 7]
                    N2DPML+=1
            N2D+=1
    #--------------------------------------------------------------------------
    print '-'*60
    print "Elements 1D, 2D : ", N1D, N2D
    if PML: print "Elements PML : ", N2DPML
    #--------------------------------------------------------------------------
    Elements      = Elements[:N2D,:]
    Milieu        = Milieu[:N2D,:]
    if PML: ElementsPML   = ElementsPML[:N2DPML,:]
    #  removal of unnecessary zeros created at initialization
    Elements1D              = Elements1D[:N1D,:]
    #
    Elements1DBordAbso      = Elements1DBordAbso[:N1DBordAbso,:]
    Elements1DBordFree      = Elements1DBordFree[:N1DBordFree,:]
    #
    Elements1DBordRight     = Elements1DBordRight[:N1DBordRight,:]
    Elements1DBordLeft      = Elements1DBordLeft[:N1DBordLeft,:]
    Elements1DBordTop       = Elements1DBordTop[:N1DBordTop,:]
    Elements1DBordBottom    = Elements1DBordBottom[:N1DBordBottom,:]

    # Creation of a vector containing only the nodes of one side
    # We take the first two nodes (end nodes) not intermediate (high order)
    #--------------------------------------------------------------------------
    NodesBordRight  = set(ravel(Elements1DBordRight[:,:2]))
    NodesBordLeft   = set(ravel(Elements1DBordLeft[:,:2]))
    NodesBordTop    = set(ravel(Elements1DBordTop[:,:2]))
    NodesBordBottom = set(ravel(Elements1DBordBottom[:,:2]))
    #--------------------------------------------------------------------------
    ctf, cta = 0, 0
    #
    print '#'*60
    print 'Identification of elements in contact with sides'
    print '# Corners  '
    #--------------------------------------------------------------------------
    #  Detection of elements in contact with sides
    
    for Ct2D in xrange(N2D):
        #
        jj=set(Elements[Ct2D,:])  # take the nodes of the element
        ct_corner = 0
        
        # Bottom
        if not set.isdisjoint(jj, NodesBordBottom):
            rr = set.intersection(jj, NodesBordBottom)
            lrr = len(rr)
            if lrr > 1:
                ela = concatenate(([Ct2D+1], [lrr], list(rr), [1]))
                elf = concatenate(([Ct2D+1], [lrr], list(rr)))
                if Bords['Bottom']=='Abso':
                    Elements2DBordAbso[cta,:] = ela
                    cta+=1
                if Bords['Bottom']=='Free':
                    Elements2DBordFree[ctf,:] = elf
                    ctf+=1
                ct_corner +=1
            
        # Right
        if not set.isdisjoint(jj, NodesBordRight):
            rr = set.intersection(jj, NodesBordRight)
            lrr = len(rr)
            if lrr > 1:
                ela = concatenate(([Ct2D+1], [lrr], list(rr), [2]))
                elf = concatenate(([Ct2D+1], [lrr], list(rr)))        
                if Bords['Right']=='Abso':
                    Elements2DBordAbso[cta,:] = ela
                    cta+=1
                if Bords['Right']=='Free':
                    Elements2DBordFree[ctf,:] = elf
                    ctf+=1
                ct_corner +=1
            
        # Top
        if not set.isdisjoint(jj, NodesBordTop):
            rr = set.intersection(jj, NodesBordTop)
            lrr = len(rr)
            if lrr > 1:
                ela = concatenate(([Ct2D+1], [lrr], list(rr), [3]))
                elf = concatenate(([Ct2D+1], [lrr], list(rr)))
                if Bords['Top']=='Abso':
                    Elements2DBordAbso[cta,:] = ela
                    cta+=1
                if Bords['Top']=='Free':
                    Elements2DBordFree[ctf,:] = elf
                    ctf+=1
                ct_corner +=1
            
        # Left
        if not set.isdisjoint(jj, NodesBordLeft):
            rr = set.intersection(jj, NodesBordLeft)
            lrr = len(rr)
            if lrr > 1:
                ela = concatenate(([Ct2D+1], [lrr], list(rr), [4]))
                elf = concatenate(([Ct2D+1], [lrr], list(rr)))
                if Bords['Left']=='Abso':
                    Elements2DBordAbso[cta,:] = ela
                    cta+=1
                if Bords['Left']=='Free':
                    Elements2DBordFree[ctf,:] = elf
                    ctf+=1
                ct_corner +=1
            
        if ct_corner > 1: print 'Corner', ct_corner, Ct2D+1
    
    #
    print '-'*60
    print 'Number of elements in contact with a free surface', ctf
    print 'Number of elements in contact with an absorbing surface', cta
    print '#'*20+' END OF PROCESSING '+'#'*20    
    # remove unnecessary zeros created at initialization
    #----------------------------------------------------------------------
    Elements2DBordAbso   = Elements2DBordAbso[:cta,:]
    Elements2DBordFree   = Elements2DBordFree[:ctf,:]
    #-----------------------------------------------------------------------
    # Save to SPECFEM format
    SauvFicSpecfem('Mesh_'+Nom, N2D, Elements, '%i')
    #
    SauvFicSpecfem('Surf_abs_'+Nom, cta, Elements2DBordAbso, '%i')
    #
    SauvFicSpecfem('Surf_free_'+Nom, ctf, Elements2DBordFree, '%i')
    #
    savetxt('Material_'+Nom,Milieu, fmt='%i')
    #
    if PML: SauvFicSpecfem('EltPML_'+Nom, N2DPML, ElementsPML, '%i')
    #
    return
if __name__=='__main__':
    set_printoptions(precision=6, threshold=None, edgeitems=None, \
    linewidth=200, suppress=None, nanstr=None, infstr=None)
    #
    # Reading input parameter if provided
    #
    Fic = sys.argv[1];                          del sys.argv[1]
    #
    # Default values
    Bords   = {'Top':'Abso', 'Bottom':'Abso', 'Left':'Abso' , 'Right':'Abso' }
    #--------------------------------------------------------------------------
    while len(sys.argv) > 1:
        opt=sys.argv[1];                           del sys.argv[1]
        if opt == '-t':
            if sys.argv[1]== 'F':
                Bords['Top']='Free'
            elif sys.argv[1]== 'A':
                Bords['Top']='Abso'
            else:
                print 'Wrong condition'
            del sys.argv[1]
        elif opt == '-b':
            if sys.argv[1]== 'F':
                Bords['Bottom']='Free'
            elif sys.argv[1]== 'A':
                Bords['Bottom']='Abso'
            else:
                print 'Wrong condition'
            del sys.argv[1]
        elif opt == '-l':
            if sys.argv[1]== 'F':
                Bords['Left']='Free'
            elif sys.argv[1]== 'A':
                Bords['Left']='Abso'
            else:
                print 'Wrong condition'
            del sys.argv[1]
        elif opt == '-r':
            if sys.argv[1]== 'F':
                Bords['Right']='Free'
            elif sys.argv[1]== 'A':
                Bords['Right']='Abso'
            else:
                print 'Wrong condition'
            del sys.argv[1]
        else:
            print sys.argv[0], ': invalid option', option
            sys.exit(1)
    #
    OuvreGmsh('',Fic, Bords)

