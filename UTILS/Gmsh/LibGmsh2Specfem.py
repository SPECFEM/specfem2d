#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Python code to link gmsh with specfem
#
#@author: Cristini Paul, 
#  Laboratoire de Mecanique et d'Acoustique, CNRS, Marseille, France
#
# Feb, 28, 2012
#
import sys, string, time
from os.path import splitext, isfile
try:
    from numpy import *
except ImportError:
    print "numpy is not installed"
#
def SauvFicSpecfem(Ng, Ct, Var, Fv):
    # Sauvegarde au format ascii
    # Ng est le nom générique
    # Ct le nombre de lignes à lire
    # Var est le nom de la variable contenant les informations à écrire
    # Fv est le format d'écriture '%f' pour les noeuds, '%i' pour les autres fichiers
    savetxt(Ng,(Ct,), fmt='%i')
    fd = open(Ng,'a')
    savetxt(fd, Var, fmt=Fv)
    fd.close()
    return
#
def OuvreGmsh(Dir,Nom,Bords):
    # Lecture de fichiers .msh genere avec Gmsh
    if splitext(Nom)[-1]=='.msh':
       fic=Nom
    elif splitext(Nom)[-1]=='':
       fic=Nom+'.msh'
    else:
        print 'File extension is not correct'
        print 'script aborted'
        sys.exit()
    #
    # Open the file and get the lines
    # 
    f = file(Dir+fic,'r')
    lignes= f.readlines()
    f.close()
    # Recherche des positions
    #MotsCles=['','']
    for ii in range(len(lignes)):
        if lignes[ii]=='$Nodes\n': PosNodes=ii 
        if lignes[ii]=='$PhysicalNames\n': PosPhys=ii 
        if lignes[ii]=='$Elements\n':
            PosElem=ii
            break
    # Type d'elements 4 noeuds ou 9 noeuds
    TypElem1D = int(string.split(lignes[PosElem+2])[1])
    if TypElem1D==1:
        Ngnod, LinElem, SurfElem = 4, 1, 3
        len1D, len2D = 2, 4
    elif TypElem1D==8:
        Ngnod, LinElem, SurfElem = 9, 8, 10
        len1D, len2D = 3, 9
    else:
        print 'Element type is not 4 or 9 nodes'
    #-------------------------------------------------------------------------
    # Conditions aux bords du domaine
    # Possible choices: Abso, Free or Perio
    Bord_abso, Bord_free = [], []  # Initialisation
    print Bords
    #-------------------------------------------------------------------------
    # PHYSICAL NAMES
    NbPhysNames = int(string.split(lignes[PosPhys+1])[0])
    print 'PhysNames', NbPhysNames
    # Structure de la variable : 1 entier et une chaine de caractere de taille 16
    dt = dtype([('dimension',int), ('zone', int), ('name', str, 16)])
    PhysCar=zeros((NbPhysNames,), dtype=dt)
    for Ip in range(NbPhysNames):
        Dim = int(string.split(lignes[PosPhys+2+Ip])[0])
        Zon = int(string.split(lignes[PosPhys+2+Ip])[1])
        Nam = string.split(lignes[PosPhys+2+Ip])[2][1:-1]
        PhysCar[Ip] = (Dim, Zon, Nam)
        if Bords.has_key(Nam):
            if Bords[Nam] == 'Abso': Bord_abso.append(Zon)
            if Bords[Nam] == 'Free': Bord_free.append(Zon)
            if Nam == 'Right':  Bord_right=Zon
            if Nam == 'Left':   Bord_left=Zon
            if Nam == 'Top':    Bord_top=Zon
            if Nam == 'Bottom': Bord_bottom=Zon
    #--------------------------------------
    print 'Physical Names', PhysCar
    print 'Bords absorbants', Bord_abso
    print 'Bords libres', Bord_free
    print 'Bords droit', Bord_right
    print 'Bords gauche', Bord_left
    print 'Bords haut', Bord_top
    print 'Bords bas', Bord_bottom
    #---------------------------------------------------------------------------
    # Infos sur le fichier Gmsh
    Ver=float(string.split(lignes[1])[0])
    File_Type=int(string.split(lignes[1])[1])
    Data_Size=int(string.split(lignes[1])[2])
    # Lecture de noeuds
    NbNodes=int(string.split(lignes[PosNodes+1])[0])
    print 'Number of nodes : ',NbNodes
    Nodes=zeros((NbNodes,2),dtype=float)
    for Ninc in range(NbNodes):
        Nodes[Ninc]= [float(val) for val in (string.split(lignes[PosNodes+2+Ninc])[1:3])]
    #
    # Sauvegarde au format SPECFEM
    SauvFicSpecfem('Nodes_'+Nom, NbNodes, Nodes, '%f')
    # Lecture des elements
    DecElem=12+NbNodes
    NbElements=int(string.split(lignes[PosElem+1])[0])
    print 'Number of elements : ', NbElements
    # depend de l'ordre
    Elements        = empty((NbElements,len2D),dtype=int)
    Milieu          = empty((NbElements,1),dtype=int)
    Elements2DBord  = empty((NbElements),dtype=int)
    Elements1D      = empty((NbElements,len1D),dtype=int)
    Elements1DBord  = empty((NbElements,len1D),dtype=int)
    #---------------------------------------------------------------------------
    Elements1DBordAbso  = empty((NbElements,len1D),dtype=int)
    Elements1DBordFree  = empty((NbElements,len1D),dtype=int)
    Elements2DBordAbso  = zeros((NbElements,4),dtype=int)
    Elements2DBordFree  = zeros((NbElements,4),dtype=int)
    #---------------------------------------------------------------------------
    Elements1DBordTop  = empty((NbElements,len1D),dtype=int)
    Elements1DBordBottom  = empty((NbElements,len1D),dtype=int)
    Elements2DBordTop     = zeros((NbElements,4),dtype=int)
    Elements2DBordBottom  = zeros((NbElements,4),dtype=int)
    #---------------------------------------------------------------------------
    Elements1DBordRight  = empty((NbElements,len1D),dtype=int)
    Elements1DBordLeft  = empty((NbElements,len1D),dtype=int)
    Elements2DBordRight = zeros((NbElements,4),dtype=int)
    Elements2DBordLeft  = zeros((NbElements,4),dtype=int)
    #---------------------------------------------------------------------------
    DecElem+=1
    Ninc1D, Ninc2D, Ninc1DBord, Ninc1DBordAbso, Ninc1DBordFree = 0, 0, 0, 0, 0
    Ninc1DBordLeft, Ninc1DBordRight, Ninc1DBordTop, Ninc1DBordBottom = 0, 0, 0, 0
    for Ninc in range(NbElements):
        Pos = PosElem+Ninc+2
        TypElem = int(string.split(lignes[Pos])[1])
        ZonP    = int(string.split(lignes[Pos])[3])
        Milieu[Ninc2D]= 1
        if TypElem==LinElem: 
            Elements1D[Ninc1D] = [int(val) for val in (string.split(lignes[Pos])[5:])]
            # Bord droit
            if ZonP==Bord_right:
                Elements1DBordRight[Ninc1DBordRight] = Elements1D[Ninc1D]
                Ninc1DBordRight+=1
                Elements1DBord[Ninc1DBord] = Elements1D[Ninc1D]
                Ninc1DBord+=1
                if Bords['Right']=='Abso':
                    Elements1DBordAbso[Ninc1DBordAbso] = Elements1D[Ninc1D]
                    Ninc1DBordAbso+=1
                else:
                    Elements1DBordFree[Ninc1DBordFree] = Elements1D[Ninc1D]
                    Ninc1DBordFree+=1
            # Bord gauche
            if ZonP==Bord_left:
                Elements1DBordLeft[Ninc1DBordLeft] = Elements1D[Ninc1D]
                Ninc1DBordLeft+=1
                Elements1DBord[Ninc1DBord] = Elements1D[Ninc1D]
                Ninc1DBord+=1
                if Bords['Left']=='Abso':
                    Elements1DBordAbso[Ninc1DBordAbso] = Elements1D[Ninc1D]
                    Ninc1DBordAbso+=1
                else:
                    Elements1DBordFree[Ninc1DBordFree] = Elements1D[Ninc1D]
                    Ninc1DBordFree+=1
            # Bord haut
            if ZonP==Bord_top:
                Elements1DBordTop[Ninc1DBordTop] = Elements1D[Ninc1D]
                Ninc1DBordTop+=1
                Elements1DBord[Ninc1DBord] = Elements1D[Ninc1D]
                Ninc1DBord+=1
                if Bords['Top']=='Abso':
                    Elements1DBordAbso[Ninc1DBordAbso] = Elements1D[Ninc1D]
                    Ninc1DBordAbso+=1
                else:
                    Elements1DBordFree[Ninc1DBordFree] = Elements1D[Ninc1D]
                    Ninc1DBordFree+=1
            # Bord bas
            if ZonP==Bord_bottom:
                Elements1DBordBottom[Ninc1DBordBottom] = Elements1D[Ninc1D]
                Ninc1DBordBottom+=1
                Elements1DBord[Ninc1DBord] = Elements1D[Ninc1D]
                Ninc1DBord+=1
                if Bords['Bottom']=='Abso':
                    Elements1DBordAbso[Ninc1DBordAbso] = Elements1D[Ninc1D]
                    Ninc1DBordAbso+=1
                else:
                    Elements1DBordFree[Ninc1DBordFree] = Elements1D[Ninc1D]
                    Ninc1DBordFree+=1          
            Ninc1D+=1
        if TypElem==SurfElem:
            Elements[Ninc2D]= [int(val) for val in (string.split(lignes[Pos])[5:])]
            Milieu[Ninc2D]= ZonP-4
            Ninc2D+=1
    Elements = Elements[:Ninc2D,:]
    Milieu   = Milieu[:Ninc2D,:]
    #
    Elements1D=Elements1D[:Ninc1D,:]
    Elements1DBord=Elements1DBord[:Ninc1DBord,:]
    Elements1DBordAbso=Elements1DBordAbso[:Ninc1DBordAbso,:]
    Elements1DBordFree=Elements1DBordFree[:Ninc1DBordFree,:]
    Elements1DBordRight=Elements1DBordRight[:Ninc1DBordRight,:]
    Elements1DBordLeft=Elements1DBordLeft[:Ninc1DBordLeft,:]
    Elements1DBordTop=Elements1DBordTop[:Ninc1DBordTop,:]
    Elements1DBordBottom=Elements1DBordBottom[:Ninc1DBordBottom,:]
    # Modification la matrice contenant les éléments 1D pour travailler sur les noeuds
    Elements1DBordFlat=ravel(Elements1DBord)
    # Noeuds appartenant aux bords du domaine
    NodesBordC=set(Elements1DBordFlat)       # permet d'enlever les elements dupliqués
    NodesBord2n=set(Elements1DBordFlat[::3]) # Noeuds aux bords sans intermediaire
    #-------------------------------------------------------
    NodesBordRight  = set(ravel(Elements1DBordRight))
    NodesBordLeft   = set(ravel(Elements1DBordLeft))
    NodesBordTop    = set(ravel(Elements1DBordTop))
    NodesBordBottom = set(ravel(Elements1DBordBottom))
    #-------------------------------------------------------
    Elt=ravel(Elements1DBordRight)
    NodesBordRight2n  = set(hstack((Elt,Elt[-2]))[::3])
    Elt=ravel(Elements1DBordLeft)
    NodesBordLeft2n   = set(hstack((Elt,Elt[-2]))[::3])
    Elt=ravel(Elements1DBordTop)
    NodesBordTop2n    = set(hstack((Elt,Elt[-2]))[::3])
    Elt=ravel(Elements1DBordBottom)
    NodesBordBottom2n = set(hstack((Elt,Elt[-2]))[::3])
    #
    ctBord=0
    # Detection des elements qui sont aux bords droit gauche haut bas
    t0=time.clock()
    ctf, cta, ctr, ctl, ctt, ctb = 0, 0, 0, 0, 0, 0
    for Ct2D in xrange(Ninc2D):
        # Attention aux éléments qui ont un point correspondant à un coin
        jj=set(Elements[Ct2D])
        if not set.isdisjoint(jj, NodesBordC):
            ctBord+=1
            # Bord Right
            if not set.isdisjoint(jj, NodesBordRight2n):
                rr = set.intersection(jj, NodesBordRight2n)
                if len(rr)==2:
                    el = concatenate(([Ct2D+1], [len(rr)], list(rr)))
                    Elements2DBordRight[ctr,:] = el
                    ctr+=1
                    if Bords['Right']=='Abso':
                       Elements2DBordAbso[cta,:] = el
                       cta+=1
                    if Bords['Right']=='Free':
                        Elements2DBordFree[ctf,:] = el
                        ctf+=1
                # Cas du bord double
                if len(rr)==3:
                    ctdb=0
                    # Recherche des deux éléments 1D
                    for Nn in xrange(Ninc1DBordRight):
                        kk  = set(Elements1DBordRight[Nn])
                        if not set.isdisjoint(rr, kk):
                            if kk.issubset(jj):
                                exec 'kk'+str(ctdb)+'=kk'
                                ctdb+=1
                    Elc = set.intersection(kk0, kk1)
                    rrm=rr.difference(Elc)
                    rr1=list(rrm)
                    El1=concatenate(([rr1[0]],list(Elc)))
                    El2=concatenate((list(Elc),[rr1[1]]))
                    el1 = concatenate(([Ct2D+1], [2], El1))
                    el2 = concatenate(([Ct2D+1], [2], El2))
                    Elements2DBordRight[ctr,:] = el1
                    ctr+=1
                    Elements2DBordRight[ctr,:] = el2
                    ctr+=1
                    if Bords['Right']=='Abso':
                        Elements2DBordAbso[cta,:] = el1
                        cta+=1
                        Elements2DBordAbso[cta,:] = el2
                        cta+=1
                    if Bords['Right']=='Free':
                        Elements2DBordFree[ctf,:] = el1
                        ctf+=1
                        Elements2DBordFree[ctf,:] = el2
                        ctf+=1
            # Bord Left
            if not set.isdisjoint(jj, NodesBordLeft2n):
                rr = set.intersection(jj, NodesBordLeft2n)
                if len(rr)==2:
                    el = concatenate(([Ct2D+1], [len(rr)], list(rr)))
                    Elements2DBordLeft[ctl,:] = el
                    ctl+=1
                    if Bords['Left']=='Abso':
                        Elements2DBordAbso[cta,:] = el
                        cta+=1
                    if Bords['Left']=='Free':
                        Elements2DBordFree[ctf,:] = el
                        ctf+=1
                # Cas du bord double
                if len(rr)==3:
                    ctdb=0
                    # Recherche des deux éléments 1D
                    for Nn in xrange(Ninc1DBordLeft):
                        kk  = set(Elements1DBordLeft[Nn])
                        if not set.isdisjoint(rr, kk):
                            if kk.issubset(jj):
                                exec 'kk'+str(ctdb)+'=kk'
                                ctdb+=1
                    Elc = set.intersection(kk0, kk1)
                    rrm=rr.difference(Elc)
                    rr1=list(rrm)
                    El1=concatenate(([rr1[0]],list(Elc)))
                    El2=concatenate((list(Elc),[rr1[1]]))
                    el1 = concatenate(([Ct2D+1], [2], El1))
                    el2 = concatenate(([Ct2D+1], [2], El2))
                    Elements2DBordLeft[ctl,:] = el1
                    ctl+=1
                    Elements2DBordLeft[ctl,:] = el2
                    ctl+=1
                    if Bords['Left']=='Abso':
                        Elements2DBordAbso[cta,:] = el1
                        cta+=1
                        Elements2DBordAbso[cta,:] = el2
                        cta+=1
                    if Bords['Left']=='Free':
                        Elements2DBordFree[ctf,:] = el1
                        ctf+=1
                        Elements2DBordFree[ctf,:] = el2
                        ctf+=1
            # Bord Top
            if not set.isdisjoint(jj, NodesBordTop2n):
                rr = set.intersection(jj, NodesBordTop2n)
                if len(rr)==2:
                    el = concatenate(([Ct2D+1], [len(rr)], list(rr)))
                    Elements2DBordTop[ctt,:len(el)] = el
                    ctt+=1
                    if Bords['Top']=='Abso':
                        Elements2DBordAbso[cta,:] = el
                        cta+=1
                    if Bords['Top']=='Free':
                        Elements2DBordFree[ctf,:] = el
                        ctf+=1
            # Bord Bottom
            if not set.isdisjoint(jj, NodesBordBottom2n):
                rr = set.intersection(jj, NodesBordBottom2n)
                if len(rr)==2:
                    el = concatenate(([Ct2D+1], [len(rr)], list(rr)))
                    Elements2DBordBottom[ctb,:len(el)] = el
                    ctb+=1
                    if Bords['Bottom']=='Abso':
                        Elements2DBordAbso[cta,:] = el
                        cta+=1
                    if Bords['Bottom']=='Free':
                        Elements2DBordFree[ctf,:] = el
                        ctf+=1
    Elements2DBord=Elements2DBord[:ctBord]
    #----------------------------------------------------------------------
    Elements2DBordAbso   = Elements2DBordAbso[:cta,:]
    Elements2DBordFree   = Elements2DBordFree[:ctf,:]
    Elements2DBordTop    = Elements2DBordTop[:ctt,:]
    Elements2DBordBottom = Elements2DBordBottom[:ctb,:]
    Elements2DBordRight  = Elements2DBordRight[:ctr,:]
    Elements2DBordLeft   = Elements2DBordLeft[:ctl,:]
    #-----------------------------------------------------------------------
    # Sauvegarde au format SPECFEM
    SauvFicSpecfem('Mesh_'+Nom, Ninc2D, Elements, '%i')
    #
    SauvFicSpecfem('Surf_abs_'+Nom, cta, Elements2DBordAbso, '%i')
    #
    SauvFicSpecfem('Surf_free_'+Nom, ctf, Elements2DBordFree, '%i')
    #
    savetxt('Material_'+Nom,Milieu, fmt='%i')
    #
    SauvFicSpecfem('Surf_top_'+Nom, ctt, Elements2DBordTop, '%i')
    #
    SauvFicSpecfem('Surf_bottom_'+Nom, ctb, Elements2DBordBottom, '%i')
    #
    SauvFicSpecfem('Surf_right_'+Nom, ctr, Elements2DBordRight, '%i')
    #
    SauvFicSpecfem('Surf_left_'+Nom, ctl, Elements2DBordLeft, '%i')
    return
if __name__=='__main__':
    set_printoptions(precision=6, threshold=None, edgeitems=None, linewidth=200, suppress=None, nanstr=None, infstr=None)
    #
    # Lecture des paramètres d'entrée éventuels
    #
    Fic = sys.argv[1];                          del sys.argv[1]
    #
    Bords={'Top':'Abso', 'Bottom':'Abso', 'Left':'Abso' , 'Right':'Abso' }
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
        elif opt == '-r':            # On affiche le résltat ou pas
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
