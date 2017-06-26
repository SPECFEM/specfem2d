#!python

import cubit

# pour creer un maillage 2D avec elements de 9 noeuds sans CPML dans Cubit
# Author: Ting Yu, CNRS and EDF, France, December 2015

# From Ting Yu: type d'elements : QUAD9 pour les elements surface et BAR3 pour les bords (bottom, right, top, left)
#  (si j'ai bien compris ton dernier mail, on a besoin seulement les deux bouts du edge pour
#  les conditions aux limites des bords (absorbants, libres). Je suis d'accord, mais pour que tous les elements
#  de surface dans le maillage ont 9 noeuds, je prefere changer le type des bords a BAR3, qui donne la
#  simplification de modifier le script.)

lx=0.08
lz=0.13
rayon=0.006

elementSize = 0.001  # matrice
elementSize2 = 0.001  # cercle

cubit.cmd('reset')
cubit.cmd('creat surf rectangle width '+str(lx)+' height '+str(lz)+' yplane')
cubit.cmd('creat surf circle radius '+str(rayon)+' yplane')
cubit.cmd('merge all')
cubit.cmd('imprint all')
cubit.cmd('curve 6 size  '+str(elementSize2))
cubit.cmd('curve 6 scheme equal')
cubit.cmd('curve 1 2 3 4 size  '+str(elementSize))
cubit.cmd('curve 1 2 3 4 scheme equal')
cubit.cmd('mesh surface 3 4 ')

cubit.cmd('block 1 face in surface 3')
cubit.cmd('block 1 name "Elastic region"')
cubit.cmd('block 1 attribute count 1')
cubit.cmd('block 1 attribute index 1 1')
cubit.cmd('block 1 element type QUAD9')

cubit.cmd('block 2 face in surface 4')
cubit.cmd('block 2 name "Inclusion"')
cubit.cmd('block 2 attribute count 1')
cubit.cmd('block 2 attribute index 1 2')
cubit.cmd('block 2 element type QUAD9')

cubit.cmd('block 3 edge in curve 1')
cubit.cmd('block 3 name "abs_bottom"')
cubit.cmd('block 3 element type BAR3')

cubit.cmd('block 4 edge in curve 2')
cubit.cmd('block 4 name "abs_left"')
cubit.cmd('block 4 element type BAR3')

cubit.cmd('block 5 edge in curve 3')
cubit.cmd('block 5 name "abs_top"')
cubit.cmd('block 5 element type BAR3')

cubit.cmd('block 6 edge in curve 4')
cubit.cmd('block 6 name "abs_right"')
cubit.cmd('block 6 element type BAR3')

cubit.cmd('set dev on')
cubit.cmd('delete face all except face in block all')

cubit.cmd('compress ids face edg node')












