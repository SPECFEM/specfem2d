# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:39:16 2012

@author: Paul Cristini, CNRS, LMA Marseille, France.

This script creates the sound profile displayed in figure_1_of_sound_profile_generated_by_the_Python_script_of_Paul.png
"""
from __future__ import print_function
import numpy as np
import pylab as pl


zc = 30.   # position of
Cdeb = 1490.  # minimum sound speed
Cinc = 40.
cc = 5.
ymax = 100.
aa = 5.   # Add a linear variation with depth

#------------------------------------------------------
y = np.linspace(0., ymax,100)
zz = cc * (np.abs(y) - zc) / zc
# Analytical sound speed variation law
vp = Cdeb + Cinc/2. * (1-np.tanh(zz)) + aa/ymax*y

# save profile for future use
y1, vp1, vp2 = y[:-1], vp[1:], vp[:-1]
print(np.shape(y1), np.shape(vp1), np.shape(vp2))
np.savetxt('SW.out',np.transpose((y1,vp1,vp2)), fmt='%4.2f')

pl.figure(figsize=(10,14))
pl.plot(vp,y,lw=1.5)
ax=pl.gca()
ax.set_ylim(ax.get_ylim()[::-1])
mm = np.array(ax.get_xlim())+[-20.,20]
ax.set_xlim(mm)
pl.ylabel('Depth (m)',fontsize=24)
pl.xlabel('Sound speed (m/s)',fontsize=24)
fontsize=20
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
ax.set_xticks(ax.get_xticks()[::2])
pl.grid()
pl.show()
