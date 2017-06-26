# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 09:40:11 2011

Add images with time step in file names if number was originally used or add
images with number in file names if time step was originally used

Usage : change directory to the the directory where the images were created and
type "python PathTo/SPECFEM2D/UTILS/ProcessJpgFiles.py"

@author: Cristini Paul, Laboratoire de Mecanique et d'Acoustique, CNRS, Marseille, France
"""
from os import listdir
from shutil import copy
from os.path import exists

def LecDataBaseNstep():
    # Open the first Database
    f = file('Database00000','r')
    lignes= f.readlines()
    f.close()
    #
    nstep=int(lignes[13])
    return nstep

def ConvertImage(ListFich):
    ct=0
    for filename in ListFich:
        if filename.startswith('image'):
            ct+=1
            fname='img'+str(ct)+'.jpg'
            copy(filename,fname)
    return

def ConvertImage000(ListFich):
    nstep=LecDataBaseNstep()
    ct=0
    for filename in ListFich:
        if filename.startswith('img'):
            ct+=1
            if ct==1:
                tstep=5
            else:
                tstep=(ct-1)*nstep
            fname='image%07i.jpg'% tstep
            copy(filename,fname)
    return

if __name__=='__main__':
     ListFich=sorted(listdir('.'))
     if not (exists('img1.jpg') & exists('image0000005.jpg')):
         if exists('img1.jpg'): ConvertImage000(ListFich)
         if  exists('image0000005.jpg'): ConvertImage(ListFich)
     else:
        print 'Nothing to do !'
