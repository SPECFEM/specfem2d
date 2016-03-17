#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 17:36:51 2016

Draw integrated energy plot

@author: bottero
"""
from __future__ import (absolute_import, division, print_function)
import numpy as np # NumPy (multidimensional arrays, linear algebra, ...)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import os,sys,glob,shutil
import argparse # To deal with arguments :
# https://docs.python.org/2/library/argparse.html
import scipy.ndimage

####################### PARSE ARGUMENTS #######################
# Here we read the argument given and we check them

parser = argparse.ArgumentParser(
    description="This script plot the files integrated_energy_fieldXXXXX representing the energy that have crossed each point")
parser.add_argument("--input_directory","-d",type=str,default="./",
                    help="input_directory : directory where we can find the files integrated_energy_fieldXXXXX")
parser.add_argument("--name_of_files","-n",type=str,default="integrated_energy_field",
                    help="name_of_files : directory where we can find the files integrated_energy_fieldXXXXX")
parser.add_argument('-p','--profiles', action='store_true',
    help='profiles: calculate energy profiles')
parser.add_argument('-nc','--no_concatenate_files', action='store_true',
    help='no_concatenate_files: don t concatenate files at the beginning of the script')
parser.add_argument("--title","-t",type=str,default="",
                    help="title : title of the figures")
parser.add_argument('-v','--verbose', action='store_true',
    help='display more information')

args = parser.parse_args()
directory=args.input_directory

# Check
if not os.path.isdir(directory):
    print("Wrong directory!")
    sys.exit(0)

if directory[-1] != '/': #If the path given don't finish by slash...
    directory=directory+'/' #... we add one

if not glob.glob(directory+args.name_of_files+"*"): # If we don't find any matching energy file...
    print("No files "+directory+args.name_of_files+"* were found!")
    sys.exit(0)

# Concatenate all the files

if not args.no_concatenate_files:
    with open(args.input_directory+args.name_of_files+"All", 'w') as outfile:
        for infile in glob.glob(directory+args.name_of_files+"0*"):
            shutil.copyfileobj(open(infile), outfile)

plt.close('all')

# Load data
x,z,intEnergy=np.loadtxt(directory+args.name_of_files+"All").T

intEnergy=10*np.log(intEnergy)
mask0=~np.isinf(intEnergy)
intEnergy[~mask0]=min(intEnergy[mask0])

# Color map to use:
cmap = cm.Greys #cm.BuPu
climMin = -200
climMax = -90

# Display limits:
xmin=0
xmax=95000
zmin=-9800
zmax=-100

# Quick ways to plot the energy using the non homogeneous grid:
##plt.tricontourf(x,z,intEnergy,20,shading='gouraud',extend="both",cmap=cmap)
#plt.tripcolor(x,z,intEnergy,shading='gouraud',cmap=cmap)
#plt.axis([xmin, xmax, zmin, zmax])
#plt.colorbar()
#plt.clim(climMin,climMax)
##cmap.set_bad('w',1.)

#%%
# Interpolation on a regular grid

# Size of regular grid
nx, nz = 500, 100

# Generate a regular grid to interpolate the data.
xil = np.linspace(xmin, xmax, nx)
zil = np.linspace(zmin, zmax, nz)
xi, zi = np.meshgrid(xil, zil)

# Interpolate using delaunay triangularization:
intEnergyi = mlab.griddata(x,z,intEnergy,xi,zi)

#Plot:
plt.figure(1)
plt.pcolormesh(xi,zi,intEnergyi,shading='gouraud',cmap=cmap)
plt.colorbar()
plt.axis([xmin, xmax, zmin, zmax])
plt.clim(climMin,climMax)
plt.title(args.title)

#%%

if args.profiles:
    # Plot energy profiles:

    def find_index(x,z,xil,zil):
        """Return the indices of the closest point in 2D array"""
        idxX=np.searchsorted(xil,x)
        idxZ=np.searchsorted(zil,z)
        return idxX,idxZ

    num=1000
    cmap2 = plt.get_cmap('prism')

    plt.figure(2)
    zVector=np.arange(-9500,0,1000)
    colors = [cmap2(i) for i in np.linspace(0, 1, len(zVector))]

    for i,zz in enumerate(zVector):
        x0,z0=0,zz
        x1,z1=95000,zz

        idxX0,idxZ0 = find_index(x0,z0,xil,zil)
        idxX1,idxZ1 = find_index(x1,z1,xil,zil)

        plt.figure(1)
        plt.hold(True)
        xLine, zLine = np.linspace(idxX0,idxX1, num), np.linspace(idxZ0, idxZ1, num)
        # Extract the values along the line, using cubic interpolation

        zi = intEnergyi[zLine.astype(np.int),xLine.astype(np.int)]

        #zi = scipy.ndimage.map_coordinates(intEnergyi, np.vstack((xLine,zLine)))
        plt.plot([x0, x1], [z0, z1], 'o-',color=colors[i])

        plt.figure(2)
        plt.plot(zi,color=colors[i])
        plt.title(args.title)

    plt.figure(3)
    xVector=np.arange(10000,100000,10000)
    colors = [cmap2(i) for i in np.linspace(0, 1, len(xVector))]
    z0=-8000
    z1=-100
    zvect=np.linspace(z0,z1,num)
    depthIntegratedEnergy=np.zeros(len(xVector))

    for i,xx in enumerate(xVector):
        x0=xx
        x1=xx

        idxX0,idxZ0 = find_index(x0,z0,xil,zil)
        idxX1,idxZ1 = find_index(x1,z1,xil,zil)

        plt.figure(1)
        plt.hold(True)
        xLine, zLine = np.linspace(idxX0,idxX1, num), np.linspace(idxZ0, idxZ1, num)

        zi = intEnergyi[zLine.astype(np.int),xLine.astype(np.int)]

        # Extract the values along the line, using cubic interpolation
        zi2 = scipy.ndimage.map_coordinates(intEnergyi.filled(), np.vstack((zLine,xLine)))
        depthIntegratedEnergy[i]=zi2.sum()

        plt.plot([x0, x1], [z0, z1], 'o-',color=colors[i])

        plt.figure(3)
        #plt.plot(zi,zvect,color=colors[i]) # Without filtering
        plt.plot(zi2,zvect,color=colors[i])
        plt.ylim([z0,z1])
        plt.title(args.title)

    plt.figure(4)
    plt.plot(xVector,depthIntegratedEnergy,'o-')
    plt.title(args.title)

plt.show()