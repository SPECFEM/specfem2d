# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 09:00:00 2011

Process Par_file to update them to the new format of release 19201

Usage : "python PathTo/SPECFEM2D/UTILS/ProcessParFileTor19201.py PathTo/filename"

@author: Cristini Paul, Laboratoire de Mecanique et d'Acoustique, CNRS, Marseille, France
"""

import sys
from shutil import move
from os.path import exists

def ProcessParfileTor19201(fic):
    # New additions to the Par_file
    a1='PERFORM_CUTHILL_MCKEE           = .true.         # perform inverse Cuthill-McKee (1969) optimization/permutation for mesh numbering\n'
    a2='USER_T0                         = 0.0d0          # use this t0 as earliest starting time rather than the automatically calculated one\n'
    a3='SU_FORMAT                       = .false.        # output seismograms in Seismic Unix format (adjoint traces will be read in the same format)\n'
    a4='factor_subsample_image          = 1              # factor to subsample color images output by the code (useful for very large models)\n'+ \
    'POWER_DISPLAY_COLOR             = 0.30d0         # non linear display to enhance small amplitudes in color images\n'+ \
    'DRAW_WATER_CONSTANT_BLUE_IN_JPG = .true.         # display acoustic layers as constant blue in JPEG images, because they likely correspond to water\n'
    a5='US_LETTER                       = .false.        # US letter paper or European A4\n'+ \
    'USE_SNAPSHOT_NUMBER_IN_FILENAME = .false.        # use snapshot number in the file name of JPEG color snapshots instead of the time step\n'
    a6='\n# for horizontal periodic conditions: detect common points between left and right edges\n'+ \
    'ADD_PERIODIC_CONDITIONS         = .false.\n\n'+ \
    '# horizontal periodicity distance for periodic conditions\n'+ \
    'PERIODIC_horiz_dist             = 0.3597d0\n\n'+ \
    '# grid point detection tolerance for periodic conditions\n'+ \
    "PERIODIC_DETECT_TOL             = 3.3334d-6\n"
    #
    f = open(fic,'r')
    ligs= f.readlines()
    f.close()
    #
    # Ajout des parametres supplementaires
    # On verifie si le fichier n'a pas deja ete traite
    if not (ligs[0].endswith('r19201\n')):
        ligs[0]=ligs[0][:-1]+' r19201\n'  # On indique que le fichier est traite pour cette release
        #
        Ct=0
        for ilg, lig in enumerate(ligs):
            if lig.startswith('partitioning'):
                ligs.insert(ilg+1,a1)
            if lig.startswith('deltat'):
                ligs.insert(ilg+1,a2)
            if lig.startswith('rec_normal'):
                ligs.insert(ilg+1,a3)
            if lig.startswith('subsamp'):
                ligs[ilg]=string.replace(ligs[ilg],'subsamp           ','subsamp_postscript',1)
                ligs.insert(ilg+1,a4)
            if lig.startswith('sizemax'):
                ligs.insert(ilg+1,a5)
            if lig.startswith('absorbing_conditions'):
                ligs.insert(ilg+1,a6)
#        #
        move(fic,fic+'.old')
        fm = open(fic,'w')
        fm.writelines(ligs)
        fm.close()
        print 'File : '+fic+' processed'
    else:
        print 'File : '+fic+' already processed'
    return
if __name__=='__main__':
    ProcessParfileTor19201(sys.argv[1])
