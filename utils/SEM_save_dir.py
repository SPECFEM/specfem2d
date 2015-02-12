#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Python code to save data obtained after running SPECFEM with the configuration files.
This file must be made executable (chmod +x) and put in the root directory of the SPECFEM package.
The new directory where information will be stored is automatically created.
The OUTPUT_FILES directory should be empty before running the xmeshfem + xspecfem sequence.

syntax : ./SEM_save_dir.py dirname

Created on Sun Nov 7 2010

@author: Paul Cristini, Laboratoire de Mecanique et d'Acoustique, CNRS, Marseille, France
"""

from os import *
import sys, shutil, string
import os.path as op

def SemSave(rep):
    SEM=getcwd()
    nvdir=op.join(SEM,rep)
    if not path.isdir(rep):
        # Copy of the entire OUTPUT_FILES directory to the new directory
        shutil.copytree(op.join(SEM,'OUTPUT_FILES'),nvdir,symlinks=False)
        # Copy of Par_file file
        shutil.copyfile(op.join(SEM,'DATA','Par_file'),op.join(nvdir,'Par_file'))
        # Copy of SOURCE file
        shutil.copyfile(op.join(SEM,'DATA','SOURCE'),op.join(nvdir,'SOURCE'))
        # Par_file reading
        filename=SEM+'/DATA/Par_file'
        f = file(filename,'r')
        lignes= f.readlines()
        f.close()
        # Save stations if generated
        if GetValuePar('generate_STATIONS',lignes)=='.true.':
            shutil.copyfile(op.join(SEM,'DATA','STATIONS'),op.join(nvdir,'STATIONS'))
        # Save configuration files
        if GetValuePar('read_external_mesh',lignes)=='.true.':
            fic=GetValuePar('mesh_file',lignes)
            shutil.copyfile(fic,op.join(nvdir,op.split(fic)[1]))
            fic=GetValuePar('nodes_coords_file',lignes)
            shutil.copyfile(fic,op.join(nvdir,op.split(fic)[1]))
            fic=GetValuePar('materials_file',lignes)
            shutil.copyfile(fic,op.join(nvdir,op.split(fic)[1]))
            fic=GetValuePar('free_surface_file',lignes)
            shutil.copyfile(fic,op.join(nvdir,op.split(fic)[1]))
            fic=GetValuePar('absorbing_surface_file',lignes)
            shutil.copyfile(fic,op.join(nvdir,op.split(fic)[1]))
        else:
            fic=GetValuePar('interfacesfile',lignes)
            shutil.copyfile(op.join(SEM,'DATA',fic),op.join(nvdir,fic))
    else:
        print 'Unable to save, directory /'+rep+' already exists. Change name !'
        
def GetValuePar(VAR,lignes):
    """ Return the values of a parameter present in the lines of a file"""
    for ligne in lignes:
        lsplit=string.split(ligne)
        if lsplit!=[]:
            if lsplit[0]==VAR:
                val=lsplit[2]
                break
    return val
    
if __name__=='__main__':
    SemSave(sys.argv[1])
