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
from scipy import interpolate

def representsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False
def representsFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class ParFile:
    """ This class is used to store the data contained on a specfem2d
    par_file"""

    def __init__(self,pathToParFile):
        """ Constructor (what happen when we do aParFile=ParFile(path)) """
        self.path=pathToParFile
        self.nt=''
        self.dt=''
        if os.path.exists(self.path):
            self.readParFile()
        else:
            raise IOError('Impossible to read '+pathToParFile)

    def readParFile(self):
        """ Open a the par_file and read some values """
        with open(self.path) as parFile:
            for line in parFile:
                if 'NSTEP=' in line.replace(" ", ""):
                    self.nt=line.split(" = ")[1].split("#")[0].strip()
                    if representsInt(self.nt): # Verify that the string extracted is a int
                        self.nt = int(self.nt)
                    else:
                        raise ValueError('Incorrect value of NSTEP read')
                if 'DT=' in line.replace(" ", ""):
                    self.dt=line.split(" = ")[1].split("#")[0].strip()
                    self.dt=self.dt.replace("d","e").replace("D","E") # Convert to real scientific fomat
                    if representsFloat(self.dt): # Verify that the string extracted is a int
                        self.dt = float(self.dt)
                    else:
                        raise ValueError('Incorrect value of DT read')

def find_index(x,z,xil,zil):
    """Return the indices of the closest point in 2D array"""
    idxX=np.searchsorted(xil,x)
    idxZ=np.searchsorted(zil,z)
    return idxX,idxZ

def interpolateValue(array,xil,zil,x,z):
    """Return the value of the 2D field described by array,xil and zil at (x,z)"""
    idxX,idxZ = find_index(x,z,xil,zil)
    xLine, zLine = np.array([idxX]), np.array([idxZ])
    # Extract the values along the line, using cubic interpolation
    if type(intEnergyi) == np.ndarray:
        zi = scipy.ndimage.map_coordinates(array, np.vstack((zLine,xLine)),order=1)
    else:
        zi = scipy.ndimage.map_coordinates(array.filled(), np.vstack((zLine,xLine)),order=1)

    return zi[0]

####################### PARSE ARGUMENTS #######################
# Here we read the argument given and we check them

parser = argparse.ArgumentParser(
    description="This script plot the files total_integrated_energy_fieldXXXXX representing the energy that have crossed each point")
parser.add_argument("--input_directory","-d",type=str,default="./",
    help="input_directory: directory where we can find the files total_integrated_energy_fieldXXXXX")
parser.add_argument("--par_file_directory",type=str,default="../",
    help="par_file_directory: directory where we can find the Par_file of the run. Default: input_directory/../")
parser.add_argument("--name_of_files","-n",type=str,default="total_integrated_energy_field",
    help="name_of_files: to plot something different than total_integrated_energy_fieldXXXXX")
parser.add_argument('-p','--profiles', action='store_true',
    help='profiles: calculate energy profiles')
parser.add_argument('-nc','--no_concatenate_files', action='store_true',
    help='no_concatenate_files: don t concatenate files at the beginning of the script')
parser.add_argument('-nl','--nolog', action='store_true',
    help='nolog: no apply log')
parser.add_argument("--title","-t",type=str,default="",
    help="title : title of the figures")
parser.add_argument("-nx",type=int,default=300,
    help="nx: number of x points for interpolated field")
parser.add_argument("-nz",type=int,default=300,
    help="nz: number of z points for interpolated field")
parser.add_argument('-w','--writeInterpolatedField', action='store_true',
    help='writeInterpolatedField: write Interpolated field on file and integrated energy a a function of range if needed')
parser.add_argument('-sp','--saveOneProfile',nargs=4,
    help='saveOneProfile: Save one X or Z profile to file (set it in source code). -sp x0 x1 z0 z1. Ex: -sp 0.0 10000.0 -300 -300')
parser.add_argument('-q','--quickDisplay', action='store_true',
    help='quickDisplay: display after a simple linear interpolation')
parser.add_argument('--displayPoints', action='store_true',
    help='displayPoints: plot the data points')
parser.add_argument('--clim',nargs=2,default=[],
    help='clim: set limits of the color bar')
parser.add_argument('--xlim',nargs=2,default=[],
    help='xlim: set the x limits of the plots')
parser.add_argument('--zlim',nargs=2,default=[],
    help='zlim: set the z limits of the plots')
parser.add_argument('--ref',nargs=2,default=[300,-2500],
    help='ref: Point that we set to 1')
parser.add_argument('--rgs', action='store_true',
    help='rgs: Compensate geometrical speading by multiplying by r')
parser.add_argument('--nxnzProfiles',nargs=2,default=[10,10],
    help='nxnz: set the number of x and z profiles to plot')
parser.add_argument('--noplot', action='store_true',
    help='noplot: do not plot anything')
parser.add_argument("--substract","-s",type=str,default="",
    help="substract : substract the field with the field given here")
parser.add_argument('-r','--reset', action='store_true',
    help='reset: delete all field previously built')
parser.add_argument('-v','--verbose', action='store_true',
    help='verbose: display more information')

args = parser.parse_args()
directory=args.input_directory

fontsize = 14
zminProfiles = -10000
zmaxProfiles = -100 #-200 # TODO (-650m for 0.5Hz, -300m for 2Hz...)
num = 2000 # Number of points to describe the profiles

# Check
if not os.path.isdir(directory):
    print("Wrong directory! "+directory)
    sys.exit(0)

if directory[-1] != '/': #If the path given don't finish by slash...
    directory=directory+'/' #... we add one

if args.par_file_directory:
    if args.par_file_directory[-1] != '/': #If the path given don't finish by slash...
        par_file_directory=args.par_file_directory+'/' #... we add one
    else:
        par_file_directory = args.par_file_directory
else:
    par_file_directory = directory+"../../"

if args.name_of_files[0] == '/': # If the full path has been given
    directory = ""

if not glob.glob(directory+args.name_of_files+"*"): # If we don't find any matching energy file...
    print("No files "+directory+args.name_of_files+"* were found!")
    sys.exit(0)

# Concatenate all the files (if specfem has been run on parallel each proc has created its own files)

if not os.path.isfile(directory+args.name_of_files+"All") or args.reset: # If the concatenation has not already been done
    if args.verbose:
        print("Concatenate files...")
    if not args.no_concatenate_files:
        with open(directory+args.name_of_files+"All", 'w') as outfile:
            for infile in glob.glob(directory+args.name_of_files+"0*"): # !!Warning!! The 0 is important!! Otherwise infinite loop
                shutil.copyfileobj(open(infile), outfile)
else:
    print(directory+args.name_of_files+"All has been found!")

#### OPTION SUBSTRACT ####
if args.substract:
    if not glob.glob(args.substract+"*"): # If we don't find any matching energy file...
        print("No files "+args.substract+"* were found!")
        sys.exit(0)
    # Concatenate all the files (if specfem has been run on parallel each proc has created its own files)

    if not os.path.isfile(args.substract+"All") or args.reset: # If the concatenation has not already been done
        if args.verbose:
            print("Concatenate files of substracted field...")
        if not args.no_concatenate_files:
            with open(args.substract+"All", 'w') as outfile:
                for infile in glob.glob(args.substract+"0*"): # !!Warning!! The 0 is important!! Otherwise infinite loop
                    shutil.copyfileobj(open(infile), outfile)
    else:
        print(args.substract+"All has been found!")

##########################

if args.verbose:
    print("Done")
plt.close('all')

# Load data
if args.verbose:
    print("Load data in "+directory+args.name_of_files+"All")
x,z,intEnergy = np.loadtxt(directory+args.name_of_files+"All").T

#### OPTION SUBSTRACT ####
if args.substract:
    if args.verbose:
        print("Load data in "+args.substract+"All")
    xSubstract,zSubstract,intEnergySubstract = np.loadtxt(args.substract+"All").T
##########################

#if args.verbose:
#    print("Load seismograms "+directory+"AA.S0001.BXX.semv AA.S0001.BXZ.semv at 300m from the source")
#t,vx0=np.loadtxt(directory+"AA.S0001.BXX.semv").T
#t,vz0=np.loadtxt(directory+"AA.S0001.BXZ.semv").T
if args.verbose:
    print("Done")

factorGs = 1.0 # Factor to compensate geometrical spreading if asked
factorGsSubstract = 1.0 # Factor to compensate geometrical spreading if asked
if args.rgs: # Remove geometrical spreading
    factorGs = x
    if args.substract:
        factorGsSubstract = xSubstract

if "integrated" in args.name_of_files: # We have to multiply by dt
    #scalingFactor = (vx0**2+vz0**2).sum()
    if args.verbose:
        print("Opening Par_file in ",par_file_directory,"...")
    par_file=ParFile(par_file_directory+'Par_file') # Open and read the Par_file
    intEnergy = intEnergy * par_file.dt * factorGs
    #intEnergy = intEnergy/scalingFactor
    if args.substract:
        intEnergySubstract = intEnergySubstract * par_file.dt * factorGsSubstract

if "max" in args.name_of_files:
    #scalingFactor = (vx0**2+vz0**2).max()
    intEnergy = intEnergy * factorGs
    #intEnergy = intEnergy/scalingFactor
    if args.substract:
        intEnergySubstract = intEnergySubstract * factorGsSubstract

mask0=~np.isinf(intEnergy)
intEnergy[~mask0]=min(intEnergy[mask0])
if args.substract:
    mask0substract=~np.isinf(intEnergySubstract)
    intEnergySubstract[~mask0substract]=min(intEnergySubstract[mask0substract])

nxProfiles = int(args.nxnzProfiles[0])
nzProfiles = int(args.nxnzProfiles[1])

# Color map to use:
cmap = cm.BuPu #cm.Greys #cm.BuPu

if args.clim:
    climMin = float(args.clim[0])
    climMax = float(args.clim[1])

# Display limits:
if args.xlim:
    xmin=float(args.xlim[0])
    xmax=float(args.xlim[1])
else:
    xmin=x.min()
    xmax=0.98*x.max()

if args.zlim:
    zmin=float(args.zlim[0])
    zmax=float(args.zlim[1])
else:
    zmin=z.min()+0.001*(z.max()-z.min())
    zmax=z.max()-0.001*(z.max()-z.min())
#print("zmin:",zmin,"zmax:",zmax)

if args.displayPoints:
    if not args.noplot:
        plt.plot(x,z,'o')
        plt.show()

if args.quickDisplay: # Quick ways to plot the energy using the non homogeneous grid:
    if not args.nolog:
         intEnergy = 10.0*np.log10(intEnergy)
    if not args.noplot:
        #plt.tricontourf(x,z,intEnergy,20,shading='gouraud',extend="both",cmap=cmap)
        plt.figure(figsize=(6,2))
        plt.tripcolor(x,z,intEnergy,shading='gouraud',cmap=cmap)
        plt.axis([xmin, xmax, zmin, zmax])
        plt.colorbar()
        if args.clim:
            plt.clim(climMin,climMax)
        plt.show()
    sys.exit()
    #cmap.set_bad('w',1.)

#%%
# Interpolation on a regular grid

# Size of regular grid
nx, nz = args.nx,args.nz

# Margins around the model
xmargin = (xmax - xmin) / 1000.0
zmargin = (zmax - zmin) / 1000.0

# Generate a regular grid to interpolate the data.
xil = np.linspace(xmin-xmargin, xmax+xmargin, nx)
zil = np.linspace(zmin-zmargin, zmax+zmargin, nz)
xi, zi = np.meshgrid(xil, zil)

#print("TODO max zil:",zil.max(),"  min zil:",zil.min())

if os.path.isfile(directory+args.name_of_files+"AllInterpolatedx"+str(nx)+"z"+str(nz)) and not args.reset: # If the interpolation has already been done and written
    if args.verbose:
        print("Interpolated field file has been found. Loading...")
    intEnergyi = np.load(directory+args.name_of_files+"AllInterpolatedx"+str(nx)+"z"+str(nz))
    if args.verbose:
        print("Done")
else:
    # Interpolate using delaunay triangularization:
    if args.verbose:
        print("Interpolation...")
    intEnergyi = mlab.griddata(x,z,intEnergy,xi,zi,interp="linear")
    if args.verbose:
        print("Done")
    if args.writeInterpolatedField:
        if args.verbose:
            print("Writing the interpolated field to file..."+directory+args.name_of_files+"AllInterpolatedx"+str(nx)+"z"+str(nz))
        intEnergyi.dump(directory+args.name_of_files+"AllInterpolatedx"+str(nx)+"z"+str(nz))
        if args.verbose:
            print("Done")

#### OPTION SUBSTRACT ####

if args.substract:
    if os.path.isfile(args.substract+"AllInterpolatedx"+str(nx)+"z"+str(nz)) and not args.reset: # If the interpolation has already been done and written
        if args.verbose:
            print("Interpolated substracted field file has been found. Loading...")
        intEnergyiSubstract = np.load(args.substract+"AllInterpolatedx"+str(nx)+"z"+str(nz))
        if args.verbose:
            print("Done")
    else:
        # Interpolate using delaunay triangularization:
        if args.verbose:
            print("Interpolation of substracted file...")
        intEnergyiSubstract = mlab.griddata(xSubstract,zSubstract,intEnergySubstract,xi,zi,interp="linear")
        if args.verbose:
            print("Done")
        if args.writeInterpolatedField:
            if args.verbose:
                print("Writing the interpolated substracted field to file..."+args.substract+"AllInterpolatedx"+str(nx)+"z"+str(nz))
            intEnergyiSubstract.dump(args.substract+"AllInterpolatedx"+str(nx)+"z"+str(nz))
            if args.verbose:
                print("Done")

    intEnergyi = abs(intEnergyi - intEnergyiSubstract) # Substract field with given other field

##########################

intEnergyi = np.log10(intEnergyi)

# Normalize
#if "max" in args.name_of_files or "integrated" in args.name_of_files:
#    if not args.substract: # TODO maybe not needed but could create problem
#        if args.verbose:
#            print("Normalizing...")
#        valueAtRef = interpolateValue(intEnergyi,xil,zil,float(args.ref[0]),float(args.ref[1]))
#        if args.verbose:
#            print("Value at reference point (",args.ref[0],",",args.ref[1],") is ",valueAtRef)
#        intEnergyi = intEnergyi - valueAtRef
#        valueAtRef = interpolateValue(intEnergyi,xil,zil,float(args.ref[0]),float(args.ref[1]))
#        if args.verbose:
#            print("Value at reference point (",args.ref[0],",",args.ref[1],") is ",valueAtRef)
#        if args.verbose:
#            print("Done")

#Plot:
if not args.noplot:
    if args.verbose:
        print("Plots...")
    plt.figure(1,figsize=(15,6))
    plt.pcolormesh(xi,zi,intEnergyi,shading='gouraud',cmap=cmap)
    plt.colorbar()
    plt.axis([xmin, xmax, zmin, zmax])
    if args.clim:
        plt.clim(climMin,climMax)
    plt.title(args.title)
    font = {'family' : 'serif','size':fontsize}
    plt.rc('font', **font)
    plt.xlabel("Range (m)",fontsize=fontsize+3)
    plt.ylabel("Depth (m)",fontsize=fontsize+3)
    #plt.rc('text', usetex=True)
    if args.verbose:
        print("Done")

if args.profiles:
    # Plot energy profiles:
    if args.verbose:
        print("Print profiles...")

    cmap2 = plt.get_cmap('prism')
    if not args.noplot and nzProfiles > 1:
        plt.figure(2)
    zVector=np.linspace(zminProfiles-zmaxProfiles,zmaxProfiles,nzProfiles) # z coordinates of horizontal profiles
    #zVector=np.arange(zmin/(nzProfiles + 1),zmin,zmin/(nzProfiles + 1)) # z coordinates of horizontal profiles
    #print("zVector:",zVector,"zmin:",zmin,"nzProfiles:",nzProfiles)
    colors = [cmap2(i) for i in np.linspace(0, 1, len(zVector))] # Color vector
    xvect=np.linspace(xmin,xmax,num) # x vector

    for i,zz in enumerate(zVector): # loop on the depths, plot all horizontal profiles in a figure (figure 2)
        x0,z0=xmin,zz
        x1,z1=xmax,zz

        idxX0,idxZ0 = find_index(x0,z0,xil,zil) # indices of the closest point in the 2D grid
        idxX1,idxZ1 = find_index(x1,z1,xil,zil) # indices of the closest point in the 2D grid
        if not args.noplot and nzProfiles > 1:
            plt.figure(1)
            plt.hold(True)
        xLine, zLine = np.linspace(idxX0,idxX1, num), np.linspace(idxZ0, idxZ1, num) # vector containing the indices
        if args.verbose:
            print("Profile 1 to be saved: (x0,z0) = (",x0,",",z0,")   (x1,z1) = (",x1,",",z1,")")
        #print("xmin:",xmin,"xmax:",xmax,"zz:",zz)
        zi = intEnergyi[zLine.astype(np.int),xLine.astype(np.int)] # If you have got an error here try to choose a lower z1 or a bigger z0! 1

        # Extract the values along the line, using cubic interpolation
        if type(intEnergyi) == np.ndarray:
            zi = scipy.ndimage.map_coordinates(intEnergyi, np.vstack((zLine,xLine)),order=1)
        else:
            zi = scipy.ndimage.map_coordinates(intEnergyi.filled(), np.vstack((zLine,xLine)),order=1)
        if not args.noplot and nzProfiles > 1:
            plt.plot([x0, x1], [z0, z1], 'o-',color=colors[i])
            plt.figure(2)
            plt.plot(xvect,zi,color=colors[i])
            plt.xlabel("Range (m)",fontsize=fontsize+3)
            if not args.nolog:
                plt.ylabel("Log of integrated energy",fontsize=fontsize+3)
            else:
                plt.ylabel("Integrated energy",fontsize=fontsize+3)
            plt.title(args.title)
    if not args.noplot and nzProfiles > 1:
        plt.xlim([xmin,xmax])

    if not args.noplot and nxProfiles > 1:
        plt.figure(3)

    xVector=np.arange(xmax/(nxProfiles + 1),xmax,xmax/(nxProfiles + 1))
    colors = [cmap2(i) for i in np.linspace(0, 1, len(xVector))]
    z0=zminProfiles
    z1=zmaxProfiles # Be careful! This point can't be too close to zmax!
    zvect=np.linspace(z0,z1,num)
    depthIntegratedEnergy=np.zeros(len(xVector))
    #depthIntegratedEnergy2=np.zeros(len(xVector))

    for i,xx in enumerate(xVector): # Loop on the ranges, plot all vertical profiles in a figure.
        x0=xx
        x1=xx

        idxX0,idxZ0 = find_index(x0,z0,xil,zil) # indices of the closest point in the 2D grid
        idxX1,idxZ1 = find_index(x1,z1,xil,zil) # indices of the closest point in the 2D grid
        if not args.noplot and nxProfiles > 1:
            plt.figure(1)
            plt.hold(True)
        if args.verbose:
            print("Profile 2 to be saved: (x0,z0) = (",x0,",",z0,")   (x1,z1) = (",x1,",",z1,")")
        xLine, zLine = np.linspace(idxX0,idxX1, num), np.linspace(idxZ0, idxZ1, num)
        #print("xx:",xmin,"xil:",xil,"zil:",zil)
        zi = intEnergyi[zLine.astype(np.int),xLine.astype(np.int)] # If you have got an error here try to choose a lower z1 or a bigger z0! 2


        # Extract the values along the line, using cubic interpolation
        if type(intEnergyi) == np.ndarray:
            zi = scipy.ndimage.interpolation.map_coordinates(intEnergyi, np.vstack((zLine,xLine)),order=1)
        else:
             zi = scipy.ndimage.interpolation.map_coordinates(intEnergyi.filled(), np.vstack((zLine,xLine)),order=1)
        #depthIntegratedEnergy[i]=zi2.sum()
        #depthIntegratedEnergy2[i]=zi.sum()

        if not args.nolog:
            depthIntegratedEnergy[i]=10*np.log10(np.power(10,zi/10.0).sum())
        else:
            depthIntegratedEnergy[i]=np.power(10,zi/10.0).sum()
        if not args.noplot and nxProfiles > 1:
            plt.plot([x0, x1], [z0, z1], 'o-',color=colors[i])
            plt.figure(3)
            plt.plot(zi,zvect,color=colors[i]) # Without filtering
            plt.xlabel("Depth (m)",fontsize=fontsize+3)
            if not args.nolog:
                plt.ylabel("Log of integrated energy",fontsize=fontsize+3)
            else:
                plt.ylabel("Integrated energy",fontsize=fontsize+3)
            #plt.plot(zi2,zvect,color=colors[i])
            plt.ylim([z0,z1])
            plt.title(args.title)
    if not args.noplot and nxProfiles > 1:
        plt.figure(4)
        plt.plot(xVector,depthIntegratedEnergy,'o-')
        #plt.plot(xVector,depthIntegratedEnergy2,'o-')
        plt.xlabel("Range (m)",fontsize=fontsize+3)
        if not args.nolog:
            plt.ylabel("Log 10 of total energy in water",fontsize=fontsize+3)
        else:
            plt.ylabel("Total energy in water",fontsize=fontsize+3)

        plt.title(args.title)
        if args.verbose:
            print("Done")

    if args.writeInterpolatedField:
        if args.verbose:
            print("Saving energy vs range...")
        np.savetxt(directory+args.name_of_files+"_energy_vs_range",np.dstack((xVector,depthIntegratedEnergy))[0])
        print("File ",directory+args.name_of_files+"_energy_vs_range has been written")

### SAVE ONE PROFILE ###

if args.saveOneProfile:
    if args.verbose:
        print("Saving one profile...")

    # properties of profile to be saved (if option --sp given):
    x0profile = float(args.saveOneProfile[0])
    x1profile = float(args.saveOneProfile[1])
    z0profile = float(args.saveOneProfile[2])
    z1profile = float(args.saveOneProfile[3])

    x0,z0 = x0profile,z0profile
    x1,z1 = x1profile,z1profile
    if z0 == z1:
        vect = np.linspace(x0,x1,num)
        xlabel = "Range (m)"
    elif x0 == x1:
        vect = np.linspace(z0,z1,num)
        xlabel = "Depth (m)"
    else:
        sys.exit("Tilted profiles are not handled for now!")
    if args.verbose:
        print("Profile to be saved: (x0,z0) = (",x0,",",z0,")   (x1,z1) = (",x1,",",z1,")")
    idxX0,idxZ0 = find_index(x0,z0,xil,zil)
    idxX1,idxZ1 = find_index(x1,z1,xil,zil)
    if not args.noplot:
        plt.figure(1)
        plt.hold(True)
    xLine, zLine = np.linspace(idxX0,idxX1, num), np.linspace(idxZ0, idxZ1, num)
    # Extract the values along the line, using cubic interpolation

    zi1 = intEnergyi[zLine.astype(np.int),xLine.astype(np.int)] # If you have got an error here try to choose a lower z1 or a bigger z0! 3

    if type(intEnergyi) == np.ndarray:
        #zi2 = scipy.ndimage.map_coordinates(np.transpose(intEnergyi), np.vstack((xLine,zLine)),order=1)
        sp = interpolate.RectBivariateSpline(zil,xil,intEnergyi, kx=3, ky=3, s=7)
    else:
        #zi2 = scipy.ndimage.map_coordinates(np.transpose(intEnergyi).filled(), np.vstack((xLine,zLine)),order=1)
        sp = interpolate.RectBivariateSpline(zil,xil,intEnergyi, kx=3, ky=3, s=7)
    if x0 == x1:
        zi = [float(sp([vect[i]],[x0])) for i in range(num)]
    if z0 == z1:
        zi = [float(sp([z0],[vect[i]])) for i in range(num)]
    #print(zi2,sp([140000.0],[-2000]),sp([140000.0],[-2500]))
    #depthIntegratedEnergy2[i]=zi.sum()
    #if not args.nolog:
    #    depthIntegratedEnergy[i]=np.power(10,zi/10.0).sum()
    #else:
    #    depthIntegratedEnergy[i]=zi.sum()
    if not args.noplot:
        plt.hold(True)
        plt.plot([x0, x1], [z0, z1], 'o-',color="black",linewidth=3)
        plt.figure(2)
        plt.plot(vect,zi,color="black",linewidth=3)
        plt.plot(vect,zi1,color="green",linewidth=3)
        #plt.plot(vect,zi2,color="red",linewidth=3)
        plt.xlabel(xlabel,fontsize=fontsize+3)
        if not args.nolog:
            plt.ylabel("Log of integrated energy",fontsize=fontsize+3)
        else:
            plt.ylabel("Integrated energy",fontsize=fontsize+3)
        if nzProfiles == 1:
           plt.xlim([x0,x1])
        plt.title(args.title)
    i = 0
    #while os.path.isfile(args.input_directory+args.name_of_files+"_profile_"+str(i)): # If the a profile of this name has already been written
    #    i = i+1
    np.savetxt(directory+args.name_of_files+"_profile_"+str(i),np.dstack((vect,zi))[0])
    print("File ",directory+args.name_of_files+"_profile_"+str(i)," has been written")
    if args.verbose:
        print("Done")

if not args.noplot:
    plt.show()

