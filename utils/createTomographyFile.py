#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Example of external velocity model
Build 1D velocity model and write it in ASCII tomo file that can be read by specfem 2d.
This is the velocity profile set:

                 vmin
               o----------------> velocity (m/s)
               |   |   /
               |   |  /
               |   | /  gradient : -d1
          ZMIN |___|/
               |    \
               |      -
               |        \  gradient : d2
               |          -
               |            \
     Depth (m) v

How to use external velocity model in specfem (03/20/2017) :

To set an arbitrary velocity model you have to use the option TOMOGRAPHY_FILE that I have implemented and that you have to set like that in the Par_file when you define the velocity model:

2 -1 9999 9999 9999 9999 9999 9999 9999 0 0 0 0 0 0

# external tomography file
TOMOGRAPHY_FILE                 = ./profile.xyz

The program understand that the velocity model number 2 has to be read in the file profile.xyz as a regular grid it will then deal with the interpolation. You can set what you want in place of the 9999, it does not matter.
The file profile.xyz has to be written under the following format:

  ! The xyz file TOMOGRAPHY_FILE that describe the tomography should be located in the TOMOGRAPHY_PATH
  ! directory, set in the Par_file. The format of the file, as read from define_external_model_from_xyz_file.f90 looks like :
  !
  ! ORIGIN_X ORIGIN_Z END_X END_Z
  ! SPACING_X SPACING_Z
  ! NX NZ
  ! VP_MIN VP_MAX VS_MIN VS_MAX RHO_MIN RHO_MAX
  ! x(1) z(1) vp vs rho
  ! x(2) z(1) vp vs rho
  ! ...
  ! x(NX) z(1) vp vs rho
  ! x(1) z(2) vp vs rho
  ! x(2) z(2) vp vs rho
  ! ...
  ! x(NX) z(2) vp vs rho
  ! x(1) z(3) vp vs rho
  ! ...
  ! ...
  ! x(NX) z(NZ) vp vs rho
  !
  ! Where :
  ! _x and z must be increasing
  ! _ORIGIN_X, END_X are, respectively, the coordinates of the initial and final tomographic
  !  grid points along the x direction (in meters)
  ! _ORIGIN_Z, END_Z are, respectively, the coordinates of the initial and final tomographic
  !  grid points along the z direction (in meters)
  ! _SPACING_X, SPACING_Z are the spacing between the tomographic grid points along the x
  !  and z directions, respectively (in meters)
  ! _NX, NZ are the number of grid points along the spatial directions x and z,
  !  respectively; NX is given by [(END_X - ORIGIN_X)/SPACING_X]+1; NZ is the same as NX, but
  !  for z direction.
  ! _VP_MIN, VP_MAX, VS_MIN, VS_MAX, RHO_MIN, RHO_MAX are the minimum and maximum values of
  !  the wave speed vp and vs (in m.s-1) and of the density rho (in kg.m-3); these values
  !  could be the actual limits of the tomographic parameters in the grid or the minimum
  !  and maximum values to which we force the cut of velocity and density in the model.
  ! _After these first four lines, in the file file_name the tomographic grid points are
  !  listed with the corresponding values of vp, vs and rho, scanning the grid along the x
  !  coordinate (from ORIGIN_X to END_X with step of SPACING_X) for each given z (from ORIGIN_Z
  !  to END_Z, with step of SPACING_Z).

This is a working script writing such file for a 1D profile (velocity depending just on depth). Read carefully all the comments but in particular the lines containing "!!WARNING!!" that I have added for
people who need range dependent models.

run that script ./createTomographyFile.py and look at the file profile.xyz created.

In the case of range dependent models vp, vs and rho will be 2D arrays, you have to modify that.

For the moment just one model can be read from an external file but it is not very difficult to implement that. Likewise viscoelastic parameters can not vary with position for now.

@author: alexis bottero (alexis dot bottero at gmail dot com)
"""

import numpy as np
import matplotlib.pyplot as plt

orig_x = 0     # Min x coordinate (take a small margin)
orig_z = -5000 # Min z coordinate
end_x  = 5000  # Max x coordinate (take a small margin)
end_z  = 0     # Max z coordinate
d1   = 0.5     # m.s-1.m-1 velocity gradient
d2   = 0.2     # m.s-1.m-1
vmin = 1500    # m.s-1

# Regular sampling (the characteristic of GLL points found in between these points will be interpolated)
nx = 20  # Number of sampling points in x direction
nz = 200 # Number of sampling points in z direction

DENSITY = 1000  # Density of water in kg.m-3. We suppose it homogeneous
ZMIN    = -1500 # Depth of minimum velocity in meters

##################### Set your velocity model below #####################

vs  = np.zeros(nz) # In that case a fluid will be described : vs = 0 m/s everywhere
rho = np.zeros(nz) + DENSITY # Homogeneous density
vp  = np.zeros(nz) # We will fill that array. !!WARNING!! if you want 2D variations you have to replace that by a 2D array.

# Fill vp array : !!WARNING!! if you want 2D variations you have to add a loop over x. vp[ix,iz] = ...
for i,zi in enumerate(z):
    if zi > ZMIN:
        vp[i]=d1*(zi-ZMIN)+vmin
    else:
        vp[i]=-d2*(zi-ZMIN)+vmin

############## Writing tomography file under specfem format : ##############
############## This format is described in specfem user manual #############

x = np.linspace(orig_x,end_x,nx)
z = np.linspace(orig_z,end_z,nz)

# Plot vp as a function of depth:
plt.plot(vp,z,'o-')

spacing_x = (end_x-orig_x)/(nx-1)
spacing_z = (end_z-orig_z)/(nz-1)

# Open a file in write mode
fo = open("profile.xyz", "w+") # Name of the file it has to be set in the Par_file
print "Name of the file: ", fo.name
line1 = str(orig_x)+" "+str(orig_z)+" "+str(end_x)+" "+str(end_z)+"\n"
line2 = str(spacing_x)+" "+str(spacing_z)+"\n"
line3 = str(nx)+" "+str(nz)+"\n"
# line4 : vpMin vpMax vsMin vsMax rhoMin rhoMax
line4 = str(min(vp))+" "+str(max(vp))+" "+str(min(vs))+" "+str(max(vs))+" "+str(min(rho))+" "+str(max(rho))+"\n"
# Write a line at the end of the file.
fo.write(line1)
fo.write(line2)
fo.write(line3)
fo.write(line4)
for iz in np.arange(nz):
    for ix in np.arange(nx):
        #!!WARNING!! if you want 2D variations you have to add a dependence of vp on x (vp[ix,iz]).
        lineToWrite= str(x[ix])+" "+str(z[iz])+" "+str(vp[iz])+" "+str(vs[iz])+" "+str(rho[iz])+"\n" # x z vp vs rho
        fo.write(lineToWrite)

plt.show()
# Close opened file
fo.close()

