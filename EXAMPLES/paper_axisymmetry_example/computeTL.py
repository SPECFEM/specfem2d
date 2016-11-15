#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Created on Wed Apr  1 19:06:20 2015
Python code to process data obtained with specfem2d
Processing of the data obtained with a slope bottom

@author: bottero
"""

from __future__ import (absolute_import, division, print_function)
import numpy as np
from numpy.fft import rfftfreq, rfft
import matplotlib.pyplot as plt
#from scipy.special import hankel1
import os.path
import glob
import argparse # To deal with arguments :
# https://docs.python.org/2/library/argparse.html

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """Filter the signal... Useless here..."""
    from math import factorial
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

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
                        raise ValueError('Incorrect value of nt read')
                if 'DT=' in line.replace(" ", ""):
                    self.dt=line.split(" = ")[1].split("#")[0].strip()
                    self.dt=self.dt.replace("d","e").replace("D","E") # Convert to real scientific fomat
                    if representsFloat(self.dt): # Verify that the string extracted is a int
                        self.dt = float(self.dt)
                    else:
                        raise ValueError('Incorrect value of deltat read')

def ReadSource(path):
    source=np.loadtxt(path)
    return source[:,0],source[:,1] #,source[:,2]

def ReadSourceCharacteristics(path):
    with open(path) as SOURCE_file:
        for line in SOURCE_file:
            if 'xs=' in line.replace(" ", ""):
                xs=line.split(" = ")[1].split("#")[0].strip()
                xs=xs.replace("d","e").replace("D","E") # Convert to real scientific fomat
                if representsFloat(xs): # Verify that the string extracted is a int
                    xs = float(xs)
                else:
                    raise ValueError('Incorrect value of xs read')
            if 'zs=' in line.replace(" ", ""):
                zs=line.split(" = ")[1].split("#")[0].strip()
                zs=zs.replace("d","e").replace("D","E") # Convert to real scientific fomat
                if representsFloat(zs): # Verify that the string extracted is a int
                    zs = float(zs)
                else:
                    raise ValueError('Incorrect value of zs read')
    return xs,zs

def ReadStations(path):
    stats=np.loadtxt(path,usecols=(2,3))
    return stats[:,0],stats[:,1]

def ReadSeismo(path):
    data=np.loadtxt(path)
    return data[:,0],data[:,1]

def zeroPad(t,ft,nZeros):
    """ From t and f(t) add points to t and zeros to f(t)
    dt must be constant : t[1]-t[0]=t[2]-t[1]=...=t[n]-t[n-1]
    """
    dt=t[1]-t[0]
    for i in np.arange(nZeros):
        t=np.append(t,t[-1]+dt)
        ft=np.append(ft,0.0)
    return t,ft

DEBUG_PLOT=False

####################### PARSE ARGUMENTS #######################
# Here we read the argument given and we check them

parser = argparse.ArgumentParser(
    description="This script creates calculate transmission losses as a function of range \
                 from specfem2d OUTPUT_FILES")
parser.add_argument("path",type=str,
                    help="path : path to the OUTPUT_FILE directory")
parser.add_argument('--NdownSampled', type=int, default = 1,
    help='To read one point over NdownSampled')
parser.add_argument('--pad', type=int, default = 1,
    help='How much we 0 pad (must be 1,2,4,8 ...) !! Costly !!')
parser.add_argument('--freq', type=float, default = 1.0,
    help='The frequency at which we want to compute transmission losses')
parser.add_argument('--outputPath',type=str, default = 'losses.txt',
    help='Path where we store the transmission Losses Files')
parser.add_argument('--noplot', action='store_true',
    help='do not display any curves')
parser.add_argument('--refAmpl', type=float, default = -1.0,
    help='Reference amplitude at the frequency given')
parser.add_argument('-v','--verbose', action='store_true',
    help='display more information')

args = parser.parse_args()

###################### OPTIONS #######################
NdownSampled=args.NdownSampled # Downsampling
pad=args.pad # How much we 0 pad (must be 1,2,4,8 ...) !! Costly !!
freq=args.freq # The frequency at which we want to compute transmission losses

##################### READ FILES #####################
# Directory where the Specfem2D files are
Dir=args.path
par_file=ParFile(Dir+'Par_file')
xs,zs = ReadSourceCharacteristics(Dir+'SOURCE') # Read SOURCE file
X, Z = ReadStations(Dir+'STATIONS')   # Read STATIONS file
X = (X - xs)/1000 # Replace X by the horizontal distance to the source and convert to kilometers
nStats=len(X)
Ndown=int(par_file.nt/NdownSampled) # We will read one point over NdownSampled
N2=2**(Ndown-1).bit_length()*pad # Smallest power of 2 greater than Ndown * pad

if args.refAmpl != -1:
    ref_amplitude = args.refAmpl
else: # If the reference amplitude is not given: calculate one from source.txt
    t_source, aS = ReadSource(Dir+'plot_source_time_function.txt') #source.txt') #'AA.S0001.PRE.semp')   # Read source.txt file
    if DEBUG_PLOT:
        plt.figure()
        plt.hold(True)
        plt.plot(t_source,aS,'+-b',label='Source brute')
    aS=savitzky_golay(aS, 51, 3) # Filtering of this signal: totally useless, I don't know why I keep it
    if DEBUG_PLOT:
        plt.plot(t_source,aS,'+-r',label='Source filtree')
    assert abs((t_source[1]-t_source[0])/par_file.dt-1) < 0.01 # Check that t_source[1]-t_source[0] ~ par_file.dt
    assert len(t_source) == par_file.nt # Check that len(t_source) ~ par_file.nt
    ############### Down sample the signal ###############
    t_source=t_source[::NdownSampled]
    aS=aS[::NdownSampled]
    if DEBUG_PLOT:
        plt.plot(t_source,aS,'+-g',label='Source filtree downsamplee')
    ########### Zero pad up to next power of 2 ###########
    t_source,aS=zeroPad(t_source,aS,N2-Ndown) # Zero-pad the source
    if DEBUG_PLOT:
        plt.plot(t_source,aS,'+-k',label='Source filtree downsamplee zeropadee')
        plt.legend()
    ############# Compute amplitude spectrum #############
    Sf = abs(rfft(aS,N2)/Ndown)#**2
    freq_source = rfftfreq(N2,d=par_file.dt*NdownSampled)
    if DEBUG_PLOT:
        plt.figure()
        plt.plot(freq_source,Sf,'o-')
        plt.title('Source spectrum')
    ##### Interpolate to calculate amplitude at freq #####
    idxInf=int(np.floor(freq/(freq_source[1]-freq_source[0])))
    idxSup=int(idxInf+1)
    ref_amplitude=((Sf[idxSup]-Sf[idxInf])/(freq_source[idxSup]-freq_source[idxInf])) \
            * freq + (Sf[idxSup]*freq_source[idxInf]-Sf[idxInf]*freq_source[idxSup])/(freq_source[idxInf]-freq_source[idxSup])
    if args.verbose:
        print("Calculated reference amplitude at",freq,'Hz:',ref_amplitude)
    if DEBUG_PLOT:
        plt.show()

#ref2=np.abs(1j/4.*hankel1(0,2*np.pi*freq/c))

############# Loop on files #############
amplitudes=np.zeros(nStats) # To store amplitudes
i = 0
nSeismos = len(sorted(glob.glob(Dir+'AA.S*')))
print("Computing "+str(nSeismos)+" FFTs ...")
for path_to_seismo in sorted(glob.glob(Dir+'AA.S*')): # All seismograms
    print("  "+path_to_seismo+" (total number: "+str(nSeismos)+")")
    ###################### READ FILE #####################
    t_seismo,pres_seismo=ReadSeismo(path_to_seismo)
    if DEBUG_PLOT:
        plt.figure()
        plt.hold(True)
        plt.plot(t_seismo,pres_seismo,'+-b',label='Seismo brut')
        plt.title(path_to_seismo)
    ############### Down sample the signal ###############
    t_seismo=t_seismo[::NdownSampled]
    pres_seismo=pres_seismo[::NdownSampled]
    #print("NdownSampled",NdownSampled," ",len(pres_seismo))
    #print("DT: ",t_seismo[1]-t_seismo[0])
    #print("Ndown: ",Ndown)
    #print("N2: ",N2)
    if DEBUG_PLOT:
        plt.plot(t_seismo,pres_seismo,'+-r',label='Seismo downsample')
    ########### Zero pad up to next power of 2 ###########
    t_seismo,pres_seismo=zeroPad(t_seismo,pres_seismo,N2-Ndown)
    if DEBUG_PLOT:
        plt.plot(t_seismo,pres_seismo,'+-g',label='Seismo downsample zeropadde')
    # (The number of points is the same than for the source)
    ############# Compute amplitude spectrum #############
    s_seismo = abs(rfft(pres_seismo,N2)/Ndown)**2
    freq_seismo = rfftfreq(N2,d=par_file.dt*NdownSampled) #t_source.shape[-1])
    #print("rfftfreq(",N2,",",par_file.dt*NdownSampled,") = ",freq_seismo[0],freq_seismo[1],freq_seismo[2]," ...")
    if DEBUG_PLOT:
        plt.figure()
        plt.plot(freq_seismo,s_seismo,'o-')
        plt.title(path_to_seismo+' spectrum')
    ###### Interpolate to calculate amplitude at freq #####
    idxInf=int(np.floor(freq/(freq_seismo[1]-freq_seismo[0])))
    idxSup=int(idxInf+1)
    amplitudes[i]=((s_seismo[idxSup]-s_seismo[idxInf])/(freq_seismo[idxSup]-freq_seismo[idxInf])) \
            * freq + (s_seismo[idxSup]*freq_seismo[idxInf]-s_seismo[idxInf]*freq_seismo[idxSup])/(freq_seismo[idxInf]-freq_seismo[idxSup])
    #print("df : ",freq_seismo[1]-freq_seismo[0])
    #print("IdxInf :",s_seismo[idxInf]," freq :",freq_seismo[idxInf])
    #print("IdxSup :",s_seismo[idxSup]," freq :",freq_seismo[idxSup])
    #print("amplitude :",amplitudes[i])
    if args.verbose:
        print("  Calculated amplitude at",freq,'Hz:',amplitudes[i])
    if DEBUG_PLOT:
        plt.show()
    i=i+1

#from scipy.io import loadmat
#fic = Dir+'/elast_s2b2sd590rd30f5_fem_dimitri_20150225.fig'
#d = loadmat(fic,squeeze_me=True, struct_as_record=False)
#print d
#
#ax1 = d['hgS_070000'].children
#if np.size(ax1) > 1:
#    ax1 = ax1[0]
#counter = 0
#for line in ax1.children:
#    if line.type == 'graph2d.lineseries':
#        x = line.properties.XData*1000
#        y = line.properties.YData
#    elif line.type == 'text':
#        if counter < 1:
#            plt.xlabel("%s" % line.properties.String,fontsize =16)
#            counter += 1
#        elif counter < 2:
#            plt.ylabel("%s" % line.properties.String,fontsize = 16)
#            counter += 1

if args.verbose:
    print("ref_amplitude:",ref_amplitude)
    print("amplitudes:",amplitudes)

loss=-10*np.log10(amplitudes/ref_amplitude)

if not args.noplot:
    plt.figure()
    plt.plot(X,loss,'b') #,marker,linestyle=linestyle,color = (r,g,b),markersize=marker_size)
    plt.gca().invert_yaxis()
    plt.show()
np.savetxt(args.outputPath,np.dstack((X,loss))[0])
np.savetxt(args.outputPath+"ampl",np.dstack((X,amplitudes))[0])
