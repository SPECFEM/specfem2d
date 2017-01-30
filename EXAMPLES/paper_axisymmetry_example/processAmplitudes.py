#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 15:05:59 2015

Small script to process transmission losses data

@author: bottero
"""

from __future__ import (absolute_import, division, print_function)
import argparse # To deal with arguments :
# https://docs.python.org/2/library/argparse.html
import numpy as np              # NumPy (multidimensional arrays, linear algebra, ...)
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    description='Small script to process transmission losses data')
parser.add_argument('files', nargs='+', type=argparse.FileType('r'),
    help='files to be plotted')
parser.add_argument('--ref',  type=float, default=1.0,
    help='reference amplitude')
parser.add_argument('--source_vel',  type=float, default=1500.0,
    help='source medium velocity')
parser.add_argument('--noplot', action='store_true',
    help='do not display any curves')
parser.add_argument('-v','--verbose', action='store_true',
    help='show more details')

args = parser.parse_args()

for seismo in args.files:      # Loop on the files given
    print(seismo.name)
    data = np.loadtxt(seismo)  # Load the seismogram
    X = data[:, 0]             # First column is range
    amplitudes = data[:, 1]    # Second column is amplitude
    TL=-10*np.log10(amplitudes/args.ref) # (2.0*args.source_vel**2)**2
    if not args.noplot:
        plt.figure()
        plt.plot(X,TL,'b') #,marker,linestyle=linestyle,color = (r,g,b),markersize=marker_size)
        plt.gca().invert_yaxis()
        plt.show()
    np.savetxt(seismo.name+"Losses",np.dstack((X,TL))[0])
