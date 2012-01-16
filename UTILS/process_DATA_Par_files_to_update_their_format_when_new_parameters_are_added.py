# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 09:00:00 2011
Updated on Fri Jan 13 2012

Processing of Par_file to update them to new format

Usage : "python PathTo/SPECFEM2D/UTILS/ProcessParFileParametersToNewRelease.py"
This will process all Par_file starting from current directory

@author: Cristini Paul, Laboratoire de Mecanique et d'Acoustique, CNRS, Marseille, France
"""
import os, string, sys
from os import listdir, walk
from shutil import copy, move
from os.path import exists
#------------------------------------------------------------------------------
def LoadLig(Fichier):
    f = open(Fichier,'r')
    ligs= f.readlines()
    f.close()
    return ligs
#------------------------------------------------------------------------------
def mylister(currdir):
    for file in os.listdir(currdir):
        path=os.path.join(currdir, file)
        if not os.path.isdir(path):
            #print path
            Fichiers.append(path)
        else:
            mylister(path)
#------------------------------------------------------------------------------
def ProcessParfile_r19201(fic):
    # Open the file and get all lines from Par_file
    ligs= LoadLig(fic)
    # Test pour voir si le traitement a déjà été fait
    for lig in ligs:
        if 'ADD_PERIODIC_CONDITIONS' in lig:
            print '----> '+fic+' already processed to r19201'            
            return
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
    'PERIODIC_DETECT_TOL             = 3.3334d-6\n'  
    #--------------------------------------------------------------------------
    # Ajout des parametres supplementaires
    # 
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
    #
    move(fic,fic+'.before_update_to_r19201')
    fm = open(fic,'w')
    fm.writelines(ligs)
    fm.close()
    #
    print 'xxxxx------> '+fic+' processed to r19201'
    return
#------------------------------------------------------------------------------
def ProcessParfile_r19340(fic):
    # Open the file and get all lines from Par_file
    ligs= LoadLig(fic)
    # Teste si le traitement a déjà été fait
    for lig in ligs:
        if 'nreceiversets' in lig:
            print '----> '+fic+' already processed to r19340'            
            return
    #
    # Ajout des parametres supplementaires
    # 
    for ilg, lig in enumerate(ligs):
        if lig.startswith('nreceiverlines'):
            ligs[ilg]=ligs[ilg].replace('lines','sets ')
    #
    move(fic,fic+'.before_update_to_r19340')
    #
    fm = open(fic,'w')
    fm.writelines(ligs)
    fm.close()
    #
    print 'xxxxx------> '+fic+' processed to r19340'
    return
#------------------------------------------------------------------------------
def ProcessParfile_r19346(fic):
    # Open the file and get all lines from Par_file
    ligs= LoadLig(fic)
    # Teste si le traitement a déjà été fait
    for lig in ligs:
        if 'ATTENUATION_PORO_FLUID_PART' in lig:
            print '----> '+fic+' already processed to r19346'            
            return
    #--------------------------------------------------------------------------
    # Ajout des parametres supplementaires
    # 
    for ilg, lig in enumerate(ligs):
        if lig.startswith('TURN_ATTENUATION_ON'):
            ligs[ilg]=ligs[ilg].replace('TURN_ATTENUATION_ON           ', \
                            'ATTENUATION_VISCOELASTIC_SOLID')
        if lig.startswith('TURN_VISCATTENUATION_ON'):
            ligs[ilg]=ligs[ilg].replace('TURN_VISCATTENUATION_ON    ', \
                            'ATTENUATION_PORO_FLUID_PART')
    #
    move(fic,fic+'.before_update_to_r19346')
    #
    fm = open(fic,'w')
    fm.writelines(ligs)
    fm.close()
    #
    print 'xxxxx------> '+fic+' processed to r19346'
    return
#------------------------------------------------------------------------------
def ProcessParfile_r19xxx(fic):
    # Open the file and get all lines from Par_file
    ligs= LoadLig(fic)
    # Teste si le traitement a déjà été fait
    for lig in ligs:
        if 'time_stepping_scheme' in lig:
            print '----> '+fic+' already processed to r19xxx'            
            return
    #
    a1='time_stepping_scheme            = 1   # 1 = Newmark (2nd order), \
    2 = LDDRK4-6 (4th-order 6-stage low storage Runge-Kutta), \
    3 = classical 4th-order 4-stage Runge-Kutta\n'

    #--------------------------------------------------------------------------
    # Ajout des parametres supplementaires
    # 
    for ilg, lig in enumerate(ligs):
        if lig.startswith('USER_T0'):
            ligs.insert(ilg+1,a1)
    #
    move(fic,fic+'.before_update_to_r19xxx')
    #
    fm = open(fic,'w')
    fm.writelines(ligs)
    fm.close()
    #
    print 'xxxxx------> '+fic+' processed to r19xxx'
    return 
#------------------------------------------------------------------------------
if __name__=='__main__':
    ## Liste de tous les fichiers à partir du répertoire courant
    Fichiers=[]
    mylister('.')
    #
    print '~'*80
    Ct_Par_file=0
    for fic in Fichiers:
        repert, ficname = os.path.split(fic)
        if not( ('.svn' in repert) or ('unused' in repert) or \
                '.before_update_to_' in ficname):
            if ficname.startswith('Par_file'):
                print 'Analysis of file : '+fic
                if not (ficname.endswith('~')):
                    Ct_Par_file+=1
                    ProcessParfile_r19201(fic)
                    ProcessParfile_r19340(fic)
                    ProcessParfile_r19346(fic)
                    ProcessParfile_r19xxx(fic)
                print '~'*80
    #                
    print 'Number of Par_file analysed : ', Ct_Par_file   
    print 'END OF Par_file PROCESSING'
    