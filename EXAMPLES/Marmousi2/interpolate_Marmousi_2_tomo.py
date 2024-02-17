#!/usr/bin/env python
#
# this script downloads Marmousi2 velocity model and interpolates them onto a sparser, regular grid;
# the regular grid can then be read in by SPECFEM2D when the DATA/Par_file uses MODEL = tomo
#
# Marmousi2 reference:
#   Martin, G. S., Wiley, R., & Marfurt, K. J. (2006). Marmousi2: An elastic upgrade for Marmousi. The leading edge, 25(2), 156-166.
#
# Marmousi2 download info:
#   https://wiki.seg.org/wiki/Open_data#AGL_Elastic_Marmousi
# with a download link:
#   https://s3.amazonaws.com/open.source.geoscience/open_data/elastic-marmousi/elastic-marmousi-model.tar.gz
# the package provides gridded Marmousi2 model files:
#   - MODEL_DENSITY_1.25m.segy.tar.gz
#   - MODEL_P-WAVE_VELOCITY_1.25m.segy.tar.gz
#   - MODEL_S-WAVE_VELOCITY_1.25m.segy.tar.gz
# these model files can also be downloaded directly from the Allied Geophysical Laboratories (AGL) homepage server,
# using these links:
#   - http://www.agl.uh.edu/downloads/density_marmousi-ii.segy.gz
#   - http://www.agl.uh.edu/downloads/vp_marmousi-ii.segy.gz
#   - http://www.agl.uh.edu/downloads/vs_marmousi-ii.segy.gz
# this script uses the latter to get the model files and stores them in a folder model_raw/.
#
###########################################################################

import sys
import os
import datetime
import array

try:
    import obspy
except:
    print("Importing obspy failed, please install obspy: pip install -U obspy")
    sys.exit(1)

import matplotlib.pyplot as plt
import numpy as np

from scipy.ndimage import gaussian_filter
from scipy.interpolate import RectBivariateSpline   # replaces: interp2d

###########################################################################
## USER Parameters

# Marmousi2 model
# Marmousi2 original domain size
MODEL_MARMOUSI2_LENGTH_X = 17000.0  # horizontal in m, original size: 17 km
MODEL_MARMOUSI2_LENGTH_Z =  3500.0  # depth in m, original size: 3.5 km

## modifications
# uses only solid model part
# (stripping off the ocean water layer from Marmousi2)
use_solid_domain_only = False

# instead of water layer, convert it to a solid
use_solid_domain_for_water_layer = False

# uses Gaussian smoothing according to given minimum frequency
use_Gaussian_smoothing = False
# maximum frequency for Gaussian smoothing
Gaussian_max_freq = 0.0

# resampling grid dimensions
# (original Marmousi2 files grid has dimensions 17km x 3.5km and a grid shape [13601, 2801], leading to a 1.25m spacing)
#
# target number of sampling points:
#   nx x nz == 3401 x 701   - to get a  5m spacing
#   nx x nz == 1701 x 351   - to get a 10m spacing
#   nx x nz ==  851 x 176   - to get a 20m spacing
#   nx x nz ==  341 x  71   - to get a 50m spacing
#
target_grid_nx = 1701   # horizontal
target_grid_nz =  351   # depth

# minimum number of grid points per wavelength
MINIMUM_PTS_PER_WAVELENGTH = 5    # for spectral-elements

# tomography model format
use_tomography_file_ascii_format = True

# plotting shows figures
show_plot = False

###########################################################################


def determine_index_fluid_solid(vs_array):
    """
    determines the last index of the fluid domain layer and the first index of the solid domain
    """
    # gets last index of the fluid domain (where vs==0)
    # this assumes that all the top layers are fluid (contiguous range of 0 to 360 for Marmousi2)
    index_fluid = -1
    for k in range(vs_array.shape[1]):
        # check if this horizontal layer has solid vs velocities
        if np.any(vs_array[:,k] > 0.0) == False: index_fluid = k

    # sets the first index of the solid domain, assuming that it is the next following index to the fluid one
    if index_fluid >= 0:
        index_solid = index_fluid + 1
    else:
        index_solid = 0

    return index_fluid,index_solid

#
#----------------------------------------------------------------------------
#

def load_marmousi2():
    """
    downloads Marmousi2 files if needed and loads in the SEGY model as numpy arrays
    """
    print("loading Marmousi2:")
    print("")

    dir = './model_raw/'
    os.system("mkdir -p "+dir)

    if not os.path.exists(os.path.join(dir, "vp_marmousi-ii.segy.gz")):
        print("  downloading marmousi-ii files...")
        os.system("wget {} -P {}".format("http://www.agl.uh.edu/downloads/vp_marmousi-ii.segy.gz", dir))
        os.system("wget {} -P {}".format("http://www.agl.uh.edu/downloads/vs_marmousi-ii.segy.gz", dir))
        os.system("wget {} -P {}".format("http://www.agl.uh.edu/downloads/density_marmousi-ii.segy.gz", dir))

    # read in data using obspy
    model_vs = obspy.read(os.path.join(dir, "vs_marmousi-ii.segy.gz"), format='segy')
    data = [tr.data for tr in model_vs ]  # segy file contains streams for each line of data
    vs = np.array(data)

    model_vp = obspy.read(os.path.join(dir, "vp_marmousi-ii.segy.gz"), format='segy')
    data = [tr.data for tr in model_vp ]  # segy file contains streams for each line of data
    vp = np.array(data)

    model_rho = obspy.read(os.path.join(dir, "density_marmousi-ii.segy.gz"), format='segy')
    data = [tr.data for tr in model_rho ]  # segy file contains streams for each line of data
    rho = np.array(data)

    # original model range
    x_range = [0.0, MODEL_MARMOUSI2_LENGTH_X]   # x in m
    z_range = [0.0, MODEL_MARMOUSI2_LENGTH_Z]   # depth in m

    print("  model: marmousi-ii")
    print("    data arrays shape = ",vs.shape)

    # convert units
    vs *= 1000.0       # from km/s to m/s
    vp *= 1000.0
    rho *= 1000.0      # from g/cm^3 to kg/m^3 

    # determine grid spacing
    nx, nz = vp.shape

    #debug
    #print("  array shape: ",vp.shape)

    x = np.linspace(x_range[0], x_range[1], nx)
    z = np.linspace(z_range[0], z_range[1], nz)

    dx = x[1] - x[0]
    dz = z[1] - z[0]

    model = {'vp': vp, 'vs': vs, 'rho': rho,
              'x': x, 'z': z, 'dx': dx, 'dz': dz}

    # gets minimum velocity for elastic domain (vs != 0)
    v_min = vs[np.where(vs > 0.0)].min()   # Vs
    if v_min == 0.0: v_min = vp.min()
    # maximum grid spacing
    d_max = max(dx,dz)
    # maximum frequency resolved (using a minimum number of points per wavelength)
    f_max = v_min / (MINIMUM_PTS_PER_WAVELENGTH * d_max)

    # stats
    print("    vp  min/max = {:8.2f} / {:8.2f} (m/s)".format(vp.min(),vp.max()))
    print("    vs  min/max = {:8.2f} / {:8.2f} (m/s)".format(vs.min(),vs.max()))
    print("    rho min/max = {:8.2f} / {:8.2f} (kg/m^3)".format(rho.min(),rho.max()))
    print("")
    print("  grid: x     range min/max = {:.2f} / {:.2f} (m)".format(x.min(),x.max()))
    print("        depth range min/max = {:.2f} / {:.2f} (m)".format(z.min(),z.max()))
    print("        nx / nz             = {} / {}".format(nx,nz))
    print("        dx / dz             = {:.2f} / {:.2f} (m)".format(dx,dz))
    print("")
    print("        resolved maximum frequency = ",f_max)
    print("")

    # plot model
    plt.figure(figsize=(10,6),dpi=150)
    plt.title("Marmousi2 - Vs original")
    plt.pcolormesh(model['x'][::5], model['z'][::5], model['vs'][::5,::5].T, shading='auto', cmap="RdYlBu_r")
    #plt.axis("scaled")
    plt.gca().invert_yaxis()
    plt.colorbar(shrink=0.3)

    # saves as JPEG file
    filename = "tmp_marmousi2_original" + ".jpg"
    plt.savefig(filename)
    print("  plotted as ",filename)
    print("")
    if show_plot:
        plt.show()

    return model

#
#----------------------------------------------------------------------------
#

def modify_original_model(model):
    """
    modifies original Marmousi2 model by taking off or replacing water layer
    """
    print("modifying original Marmousi2 model:")
    print("")

    # original model
    vs = model['vs']
    vp = model['vp']
    rho = model['rho']

    print("  original model:")
    print("    data arrays shape = ",vs.shape)

    # Find the z-index above which all vs values are non-zero
    # (the original Marmousi2 should have the top 361 layers all fluid)
    index_fluid,index_solid = determine_index_fluid_solid(vs)

    if index_fluid >= 0:
        print("    fluid layers      = {:4d} to {:4d}".format(0,index_fluid))
    print("    solid layers      = {:4d} to {:4d}".format(index_solid,vs.shape[1]))
    print("")

    # stripping off the ocean layer
    if use_solid_domain_only:
        print("  using only solid domain part of Marmousi2")

        # fluid ocean layer
        # Marmousi2 has a water layer of 450 m depth
        # (with a grid sampling of 1.25 m, this leads to a fluid layer with index between 0 to 360)
        if index_fluid >= 0:
            fluid_layer_depth = index_fluid * 1.25 # original grid sampling uses 1.25 m
            # takes off the fluid layer
            vs = vs[:,index_solid::]
            vp = vp[:,index_solid::]
            rho = rho[:,index_solid::]
        else:
            fluid_layer_depth = 0.0

        print("  fluid layer depth = {} (m)".format(fluid_layer_depth))
        print("")
        print("  taking off fluid layer: new model depth = ",MODEL_MARMOUSI2_LENGTH_Z - fluid_layer_depth)
        print("")

        # solid model range
        x_range = [0.0, MODEL_MARMOUSI2_LENGTH_X]  # x in m
        z_range = [0.0, MODEL_MARMOUSI2_LENGTH_Z - fluid_layer_depth]   # depth in m

    elif use_solid_domain_for_water_layer:
        print("  using solid domain for water layer in Marmousi2")

        # fluid ocean layer
        # Marmousi2 has a water layer of 450 m depth
        # (with a grid sampling of 1.25 m, this leads to a fluid layer with index between 0 to 360)
        if index_fluid >= 0:
            fluid_layer_depth = index_fluid * 1.25 # original grid sampling uses 1.25 m
            print("  fluid layer depth = {} (m)".format(fluid_layer_depth))
            print("")
            
            # replaces fluid layer with extended solid domain velocities
            # takes top layer from solid domain
            vs_top_solid = vs[:,index_solid]
            vp_top_solid = vp[:,index_solid]
            rho_top_solid = rho[:,index_solid]

            print("  solid top layer velocities: vs  min/max = {:8.2f} / {:8.2f}".format(vs_top_solid.min(),vs_top_solid.max()))
            print("                              vp  min/max = {:8.2f} / {:8.2f}".format(vp_top_solid.min(),vp_top_solid.max()))
            print("                              rho min/max = {:8.2f} / {:8.2f}".format(rho_top_solid.min(),rho_top_solid.max()))
            print("")

            # replaces fluid layer with top velocities
            print("  replacing fluid with solid top velocities...")
            print("")
            for i in range(index_fluid+1):
                vs[:,i] = vs_top_solid
                vp[:,i] = vp_top_solid
                rho[:,i] = rho_top_solid
        else:
            fluid_layer_depth = 0.0

        # stays at full model range
        x_range = [0.0, MODEL_MARMOUSI2_LENGTH_X]   # x in m
        z_range = [0.0, MODEL_MARMOUSI2_LENGTH_Z]   # depth in m


    else:
        print("  using original Marmousi2")

        # full model range
        x_range = [0.0, MODEL_MARMOUSI2_LENGTH_X]   # x in m
        z_range = [0.0, MODEL_MARMOUSI2_LENGTH_Z]   # depth in m

    # determine new grid spacing
    print("  new array shape: ",vp.shape)
    print("")
    nx, nz = vp.shape

    x = np.linspace(x_range[0], x_range[1], nx)
    z = np.linspace(z_range[0], z_range[1], nz)

    dx = x[1] - x[0]
    dz = z[1] - z[0]

    # replace with new model
    model = {'vp': vp, 'vs': vs, 'rho': rho,
              'x': x, 'z': z, 'dx': dx, 'dz': dz}

    # gets minimum velocity for elastic domain (vs != 0)
    v_min = vs[np.where(vs > 0.0)].min()   # Vs
    if v_min == 0.0: v_min = vp.min()
    # maximum grid spacing
    d_max = max(dx,dz)
    # maximum frequency resolved (using a minimum number of points per wavelength)
    f_max = v_min / (MINIMUM_PTS_PER_WAVELENGTH * d_max)

    # stats
    print("  modified model:")
    print("    vp  min/max = {:8.2f} / {:8.2f} (m/s)".format(vp.min(),vp.max()))
    print("    vs  min/max = {:8.2f} / {:8.2f} (m/s)".format(vs.min(),vs.max()))
    print("    rho min/max = {:8.2f} / {:8.2f} (kg/m^3)".format(rho.min(),rho.max()))
    print("")
    print("    grid: x     range min/max = {:.2f} / {:.2f} (m)".format(x.min(),x.max()))
    print("          depth range min/max = {:.2f} / {:.2f} (m)".format(z.min(),z.max()))
    print("          nx / nz             = {} / {}".format(nx,nz))
    print("          dx / dz             = {:.2f} / {:.2f} (m)".format(dx,dz))
    print("")
    print("          resolved maximum frequency = ",f_max)
    print("")

    # plot model
    plt.figure(figsize=(10,6),dpi=150)
    plt.title("Marmousi2 - Vs modified")
    plt.pcolormesh(model['x'][::5], model['z'][::5], model['vs'][::5,::5].T, shading='auto', cmap="RdYlBu_r")
    #plt.axis("scaled")
    plt.gca().invert_yaxis()
    plt.colorbar(shrink=0.3)

    # saves as JPEG file
    filename = "tmp_marmousi2_modified" + ".jpg"
    plt.savefig(filename)
    print("  plotted as ",filename)
    print("")
    if show_plot:
        plt.show()

    return model

#
#----------------------------------------------------------------------------
#


def smooth_model(max_freq,model):
    """
    uses a Gaussian to smooth the input model, with a Gaussian width depending on the maximum frequency
    """
    global Gaussian_max_freq

    print("smoothing model:")
    print("  target maximum frequency         : ",max_freq)
    print("")
    Gaussian_max_freq = max_freq

    dim_x = model['x'].max() - model['x'].min()
    dim_z = model['z'].max() - model['z'].min()

    print("  model dimensions: x / z = {} / {}".format(dim_x,dim_z))
    print("")

    # note: Marmousi2 has an ocean layer on top, using a Gaussian to smooth the model (or resampling) will smooth out
    #       the fluid/solid interface as well. we thus separate first the solid model and only smooth in this solid domain.
    #       (and similar to resampling the model)
    #
    # Find the z-index above which all vs values are non-zero
    # (the original Marmousi2 should have the top 361 layers all fluid)
    # Find the z-index above which all vs values are non-zero
    # (the original Marmousi2 should have the top 361 layers all fluid)
    index_fluid,index_solid = determine_index_fluid_solid(model['vs'])

    # w/ fluid domain
    if index_fluid >= 0:
        vs_fluid = model['vs'][:,0:index_solid]
        vp_fluid = model['vp'][:,0:index_solid]
        rho_fluid = model['rho'][:,0:index_solid]

    # solid domain only
    vs_solid = model['vs'][:,index_solid::]
    vp_solid = model['vp'][:,index_solid::]
    rho_solid = model['rho'][:,index_solid::]

    if index_fluid >= 0:
        print("  original fluid domain shape   = ",vs_fluid.shape)
    print("  original solid domain shape   = ",vs_solid.shape)
    print("")

    # gets original maximum velocity from solid domain
    vp_min = vp_solid.min()
    vp_max = vp_solid.max()

    vs_min = vs_solid.min()
    vs_max = vs_solid.max()

    print("  model vp solid min/max = {:.2f} / {:.2f}".format(vp_min,vp_max))
    print("        vs solid min/max = {:.2f} / {:.2f}".format(vs_min,vs_max))
    print("")

    # debug: show solid model
    if 1 == 0:
        plt.figure(figsize=(10,6),dpi=150)
        plt.title("Marmousi2 - Vs solid")
        plt.pcolormesh(model['x'][::5], model['z'][index_solid::5], vs_solid[::5,::5].T, shading='auto', cmap="RdYlBu_r")
        #plt.axis("scaled")
        plt.gca().invert_yaxis()
        plt.colorbar(shrink=0.3)
        plt.show()

    # gaussian filter to smooth model
    #
    # we'll give some frequency informations like minimum wavelengths and sampling interval,
    # then compute the sampling ratios to determine an appropriate Gaussian lengths in X/Z-direction
    # gets minimum/ maximum wavelength
    lambda_max = 1.0/Gaussian_max_freq * vp_max
    lambda_min = 1.0/Gaussian_max_freq * vs_min

    print("  Gaussian smoothing:")
    print("    target frequency: ",Gaussian_max_freq)
    print("    wavelengths     : min / max = {:.2f} / {:.2f} (m)".format(lambda_min,lambda_max))
    print("")

    # sampling increment for minimum wavelength
    sampling_incr = lambda_min / MINIMUM_PTS_PER_WAVELENGTH

    ratio_x = sampling_incr / model['dx']
    ratio_z = sampling_incr / model['dz']

    print("    minimum points per wavelength: ",MINIMUM_PTS_PER_WAVELENGTH)
    print("    sampling length              : increment = {} (m)".format(sampling_incr))
    print("    sampling ratios              : x   = ",ratio_x)
    print("                                   z   = ",ratio_z)
    print("")

    if sampling_incr < dim_x/target_grid_nx:
        print("Error: sampling in x-direction needs a higher target grid shape : please set target_grid_nx > ",int(dim_x/sampling_incr))
        print("")
        sys.exit(1)

    if sampling_incr < dim_z/target_grid_nz:
        print("Error: sampling in z-direction needs a higher target grid shape : please set target_grid_nz > ",int(dim_z/sampling_incr))
        print("")
        sys.exit(1)

    # note: we use a pre-defined grid sampling by target_grid_nx/target_grid_nz;
    #       instead, we could also use a grid sampling defined by the minimum sampling length required for the target frequency.
    #       however, this could lead to a very sparse grid for low-frequency models and thus exhibit gridding artifacts.
    #
    # according to sampling increment
    #x_new = np.arange(model['x'].min(), model['x'].max(), sampling_incr)
    #z_new = np.arange(model['z'].min(), model['z'].max(), sampling_incr)


    #debug
    #print("debug: solid index   = ",index_solid)

    # determines Gaussian width
    # scalelength: approximately S ~ sigma * sqrt(8.0) for a Gaussian smoothing
    # sets width such that scalelength becomes 1/2 of minimum (or dominant) wavelength
    # standard deviation for Gaussian kernel
    sigma_x = 0.5 * lambda_min / np.sqrt(8.0)
    sigma_z = 0.5 * lambda_min / np.sqrt(8.0)

    print("    using Gaussian widths equivalent to one half of the minimum wavelength")
    print("    scalelengths: x / z = {:.2f} / {:.2f} (m)".format(sigma_x * np.sqrt(8.0),sigma_z * np.sqrt(8.0)))
    print("    sigma       : x / z = {:.2f} / {:.2f}".format(sigma_x,sigma_z))
    print("")

    # Gaussian smoothing
    # only smooths the solid domain part
    vs_smooth  = gaussian_filter(vs_solid, [sigma_x, sigma_z], mode='reflect')
    vp_smooth  = gaussian_filter(vp_solid, [sigma_x, sigma_z], mode='reflect')
    rho_smooth = gaussian_filter(rho_solid, [sigma_x, sigma_z], mode='reflect')

    # re-combines fluid/solid domains
    vp = model['vp']
    vp[:,index_solid::] = vp_smooth
    vs = model['vs']
    vs[:,index_solid::] = vs_smooth
    rho = model['rho']
    rho[:,index_solid::] = rho_smooth

    print("  Gaussian filtered model: vp shape   = ",vp.shape)
    print("                              min/max = {} / {}".format(vp.min(),vp.max()))
    print("")

    # replace with smooth velocities
    model['vp'] = vp
    model['vs'] = vs
    model['rho'] = rho

    # plot model
    plt.figure(figsize=(10,6),dpi=150)
    plt.title("Marmousi2 - Vs smoothed")
    plt.pcolormesh(model['x'][::5], model['z'][::5], model['vs'][::5,::5].T, shading='auto', cmap="RdYlBu_r")
    #plt.axis("scaled")
    plt.gca().invert_yaxis()
    plt.colorbar(shrink=0.3)

    # saves as JPEG file
    filename = "tmp_marmousi2_smoothed" + ".jpg"
    plt.savefig(filename)
    print("  plotted as ",filename)
    print("")
    if show_plot:
        plt.show()

    return model

#
#----------------------------------------------------------------------------
#

def resample_model(model):
    """
    resamples the input model onto a coarser grid as specified by the user parameters
    """
    print("resampling model:")
    print("  target grid dimensions   = {} / {}".format(target_grid_nx,target_grid_nz))
    print("")

    # gets maximum velocity
    vp_max = model['vp'].max()
    vp_min = model['vp'].min()

    vs_max = model['vs'].max()
    vs_min = model['vs'][np.where(model['vs'] > 0.0)].min()

    print("  model vp min/max = {:.2f} / {:.2f}".format(vp_min,vp_max))
    print("        vs min/max = {:.2f} / {:.2f}".format(vs_min,vs_max))
    print("")

    # target gridding
    x_new = np.linspace(model['x'].min(), model['x'].max(), target_grid_nx)
    z_new = np.linspace(model['z'].min(), model['z'].max(), target_grid_nz)

    print("  new grid length x  = ",x_new.size)
    print("           length z  = ",z_new.size)
    print("")
    print("  input model            : vp shape   = ",model['vp'].shape)
    print("                              min/max = {} / {}".format(model['vp'].min(),model['vp'].max()))
    print("")

    # note: Marmousi2 has an ocean layer on top, using a Gaussian to smooth the model (or resampling) will smooth out
    #       the fluid/solid interface as well. we thus separate first the solid model and only smooth in this solid domain.
    #       (and similar to resampling the model)
    #
    # Find the z-index above which all vs values are non-zero
    # (the original Marmousi2 should have the top 361 layers all fluid)
    index_fluid,index_solid = determine_index_fluid_solid(model['vs'])

    # grid interpolation
    print("  resampling:")
    if index_fluid >= 0:
        # separates fluid and solid domain parts
        # fluid domain
        vs_fluid = model['vs'][:,0:index_solid]
        vp_fluid = model['vp'][:,0:index_solid]
        rho_fluid = model['rho'][:,0:index_solid]

        # solid domain
        vs_solid = model['vs'][:,index_solid::]
        vp_solid = model['vp'][:,index_solid::]
        rho_solid = model['rho'][:,index_solid::]

        print("    input fluid domain shape   = ",vs_fluid.shape)
        print("    input solid domain shape   = ",vs_solid.shape)
        print("")

        # debug: show solid model
        if 1 == 0:
            plt.figure(figsize=(10,6),dpi=150)
            plt.title("Marmousi2 - Vs solid input")
            plt.pcolormesh(model['x'][::5], model['z'][index_solid::5], vs_solid[::5,::5].T, shading='auto', cmap="RdYlBu_r")
            #plt.axis("scaled")
            plt.gca().invert_yaxis()
            plt.colorbar(shrink=0.3)
            plt.show()

        # separates resampling for fluid and elastic parts
        nx_org,nz_org = model['vs'].shape
        index_fluid_resampled = int( index_fluid * (target_grid_nz - 1.0) / (nz_org - 1.0) )
        index_solid_resampled = index_fluid_resampled + 1

        print("    resampled fluid index = ",index_fluid_resampled)
        print("    resampled solid index = ",index_solid_resampled)
        print("")

        if index_fluid_resampled <= 0:
            print("Error: resampled fluid/solid domain parts have invalid fluid index ",index_fluid_resampled)
            print("       index for resampled fluid domain must be > 0, please modify size of target_grid_nz == ",target_grid_nz)
            sys.exit(1)

        print("    fluid size x / z   = {} / {}".format(model['x'][:].size,model['z'][0:index_solid].size))
        print("    fluid domain shape = ",vp_fluid.shape)
        print("")
        print("    solid size x / z   = {} / {}".format(model['x'][:].size,model['z'][index_solid::].size))
        print("    solid domain shape = ",vp_solid.shape)
        print("")

        # fluid domain part
        # deprecated: vp = interp2d(model['z'], model['x'], vp, kind='cubic')(z_new, x_new)
        vp_fluid_new  = RectBivariateSpline(model['z'][0:index_solid], model['x'][:], vp_fluid.T)(z_new[0:index_solid_resampled], x_new[:])
        vs_fluid_new  = RectBivariateSpline(model['z'][0:index_solid], model['x'][:], vs_fluid.T)(z_new[0:index_solid_resampled], x_new[:])
        rho_fluid_new = RectBivariateSpline(model['z'][0:index_solid], model['x'][:], rho_fluid.T)(z_new[0:index_solid_resampled], x_new[:])

        # transpose result
        vp_fluid_new = vp_fluid_new.T
        vs_fluid_new = vs_fluid_new.T
        rho_fluid_new = rho_fluid_new.T

        # solid domain part
        vp_solid_new  = RectBivariateSpline(model['z'][index_solid::], model['x'][:], vp_solid.T)(z_new[index_solid_resampled::], x_new[:])
        vs_solid_new  = RectBivariateSpline(model['z'][index_solid::], model['x'][:], vs_solid.T)(z_new[index_solid_resampled::], x_new[:])
        rho_solid_new = RectBivariateSpline(model['z'][index_solid::], model['x'][:], rho_solid.T)(z_new[index_solid_resampled::], x_new[:])

        # transpose result
        vp_solid_new = vp_solid_new.T
        vs_solid_new = vs_solid_new.T
        rho_solid_new = rho_solid_new.T

        # combines model domains
        vp = np.concatenate((vp_fluid_new, vp_solid_new), axis=1)
        vs = np.concatenate((vs_fluid_new, vs_solid_new), axis=1)
        rho = np.concatenate((rho_fluid_new, rho_solid_new), axis=1)

    else:
        # solid model
        # deprecated: vp = interp2d(model['z'], model['x'], vp, kind='cubic')(z_new, x_new)
        vp  = RectBivariateSpline(model['z'], model['x'], model['vp'].T)(z_new, x_new)
        vs  = RectBivariateSpline(model['z'], model['x'], model['vs'].T)(z_new, x_new)
        rho = RectBivariateSpline(model['z'], model['x'], model['rho'].T)(z_new, x_new)

        # transpose result
        vp = vp.T
        vs = vs.T
        rho = rho.T

    print("  interpolated model     : vp shape   = ",vp.shape)
    print("                              min/max = {} / {}".format(vp.min(),vp.max()))
    print("")

    # Find the z-index above which all vs values are non-zero
    index_fluid,index_solid = determine_index_fluid_solid(vs)

    # w/ fluid domain
    if index_fluid >= 0:
        vs_fluid = vs[:,0:index_solid]
        vp_fluid = vp[:,0:index_solid]
        rho_fluid = rho[:,0:index_solid]

        # solid domain only
        vs_solid = vs[:,index_solid::]
        vp_solid = vp[:,index_solid::]
        rho_solid = rho[:,index_solid::]

        print("  resampled fluid domain shape  = ",vs_fluid.shape)
        print("  resampled solid domain shape  = ",vs_solid.shape)
        print("")

    # return updated model
    dx = x_new[1] - x_new[0]
    dz = z_new[1] - z_new[0]
    nx, nz = vp.shape

    new_model = {'vp': vp, 'vs': vs, 'rho': rho,
                 'x': x_new, 'z': z_new, 'dx': dx, 'dz': dz}

    # gets minimum velocity for elastic domain (vs != 0)
    v_min = vs[np.where(vs > 0.0)].min()   # Vs
    if v_min == 0.0: v_min = vp.min()
    # maximum grid spacing
    d_max = max(dx,dz)
    # maximum frequency resolved (using a minimum number of points per wavelength)
    f_max = v_min / (MINIMUM_PTS_PER_WAVELENGTH * d_max)

    print("  resampled grid: x     range min/max = {:.2f} / {:.2f} (m)".format(x_new.min(),x_new.max()))
    print("                  depth range min/max = {:.2f} / {:.2f} (m)".format(z_new.min(),z_new.max()))
    print("                  nx / nz  = {} / {}".format(nx,nz))
    print("                  dx / dz  = {} / {} (m)".format(dx,dz))
    print("")
    print("                  resolved maximum frequency = ",f_max)
    print("")

    # plot model
    plt.figure(figsize=(10,6),dpi=150)
    plt.title("Marmousi2 - Vs resampled")
    plt.pcolormesh(new_model['x'], new_model['z'], new_model['vs'].T, cmap="RdYlBu_r")
    #plt.axis("scaled")
    plt.gca().invert_yaxis()
    plt.colorbar(shrink=0.3)

    # saves as JPEG file
    filename = "tmp_marmousi2_resampled" + ".jpg"
    plt.savefig(filename)
    print("  plotted as ",filename)
    print("")
    if show_plot:
        plt.show()

    return new_model

#
#----------------------------------------------------------------------------
#

def create_tomography_file(model):
    """
    creates an ascii tomography file in SPECFEM2D format
    """
    global Gaussian_max_freq
    print("creating tomography file:")

    # initializes header variables
    header_origin_x = sys.float_info.max
    header_origin_z = sys.float_info.max

    header_end_x = -sys.float_info.max
    header_end_z = -sys.float_info.max

    header_vp_min = sys.float_info.max
    header_vp_max = -sys.float_info.max

    header_vs_min = sys.float_info.max
    header_vs_max = -sys.float_info.max

    header_rho_min = sys.float_info.max
    header_rho_max = -sys.float_info.max

    # model dimensions
    nx, nz = model['vp'].shape

    # header infos
    header_dx = model['dx']        # increment in x-direction
    header_dz = model['dz']        # increment in x-direction

    header_nx = nx   # x-direction
    header_nz = nz   # z-direction

    print("  dimension:")
    print("    nx / nz = {} / {}".format(header_nx,header_nz))
    print("    dx / dz = {} / {} (m)".format(header_dx,header_dz))

    # note: the current model arrays have x/z coordinates with z being depth.
    #       for SPECFEM2D, we need instead Z-coordinates where z-coordinates would point positive upwards.
    #       thus, we will reverse the looping direction and loop to go bottom-up, instead of top-down.
    ztop_model = model['z'].max()

    # collect model data
    if use_tomography_file_ascii_format:
        output_data = list()
    else:
        output_data_array = np.empty([5,nx*nz],dtype='f')
        iline = 0

    # loops over z
    for j in range(nz):
        # loops over x (must have inner loop over x)
        for i in range(nx):
            # reverse indexing direction to bottom-up
            k = nz - 1 - j

            # position
            x_val = model['x'][i]
            depth_val = model['z'][k]

            # will need to flip z-values to have z-direction positive up
            # from:     0.0 == top-water   , 3500.0 == bottom-layer of model
            # to  :  3500.0 == top-layer   ,    0.0 == bottom-layer of model
            z_val = ztop_model - depth_val

            # model parameters
            vp_val = model['vp'][i,k]
            vs_val = model['vs'][i,k]
            rho_val = model['rho'][i,k]

            # data line format:
            #x #z #vp #vs #density (see in src/specfem2d/define_external_model_from_tomo_file.f90)
            if use_tomography_file_ascii_format:
                # ascii format output
                output_data.append("{} {} {} {} {}\n".format(x_val,z_val,vp_val,vs_val,rho_val))
            else:
                # binary format
                output_data_array[:,iline] = [x_val,z_val,vp_val,vs_val,rho_val]
                iline += 1

            # header stats
            header_origin_x = min(header_origin_x,x_val)
            header_origin_z = min(header_origin_z,z_val)

            header_end_x = max(header_end_x,x_val)
            header_end_z = max(header_end_z,z_val)

            header_vp_min = min(header_vp_min,vp_val)
            header_vp_max = max(header_vp_max,vp_val)

            header_vs_min = min(header_vs_min,vs_val)
            header_vs_max = max(header_vs_max,vs_val)

            header_rho_min = min(header_rho_min,rho_val)
            header_rho_max = max(header_rho_max,rho_val)

    print("    x origin / end = {} / {} (m)".format(header_origin_x,header_end_x))
    print("    z origin / end = {} / {} (m)".format(header_origin_z,header_end_z))
    print("")

    print("  tomographic model statistics:")
    print("    vp  min/max : {:8.2f} / {:8.2f} (m/s)".format(header_vp_min,header_vp_max))
    print("    vs  min/max : {:8.2f} / {:8.2f} (m/s)".format(header_vs_min,header_vs_max))
    print("    rho min/max : {:8.2f} / {:8.2f} (kg/m3)".format(header_rho_min,header_rho_max))
    print("")

    if header_vp_min <= 0.0:
        print("WARNING: Vp has invalid entries with a minimum of zero!")
        print("         The provided output model is likely invalid.")
        print("         Please check with your inputs...")
        print("")

    if header_vs_min <= 0.0:
        print("WARNING: Vs has entries with a minimum of zero!")
        print("         The provided output model is likely invalid.")
        print("         Please check with your inputs...")
        print("")

    if header_rho_min <= 0.0:
        print("WARNING: Density has entries with a minimum of zero!")
        print("         The provided output model is likely invalid.")
        print("         Please check with your inputs...")
        print("")

    # data header
    data_header = list()
    data_header.append("# tomography model - converted using script interpolate_Marmousi_2_tomo.py\n")
    data_header.append("#\n")

    # providence
    data_header.append("# providence\n")
    data_header.append("# created             : {}\n".format(str(datetime.datetime.now())))
    data_header.append("# command             : {}\n".format(" ".join(sys.argv)))
    data_header.append("#\n")

    # tomographic model format
    data_header.append("# model format\n")
    data_header.append("# model type          : Marmousi2\n")
    if use_Gaussian_smoothing:
        data_header.append("#                       Gaussian smoothed for maximum frequency: {} (Hz)\n".format(Gaussian_max_freq))
    if use_solid_domain_only:
        data_header.append("#                       using only solid domain part of Marmousi2\n")
    if use_solid_domain_for_water_layer:
        data_header.append("#                       using solid domain for water layer in Marmousi2\n")
    data_header.append("# coordinate format   : x / z  # z-direction (positive up)\n")
    data_header.append("#\n")

    # tomography model header infos
    #origin_x #origin_y #origin_z #end_x #end_y #end_z          - start-end dimensions
    data_header.append("#origin_x #origin_z #end_x #end_z\n")
    data_header.append("{} {} {} {}\n".format(header_origin_x,header_origin_z,header_end_x,header_end_z))
    #dx #dy #dz                                                 - increments
    data_header.append("#dx #dz\n")
    data_header.append("{} {}\n".format(header_dx,header_dz))
    #nx #ny #nz                                                 - number of models entries
    data_header.append("#nx #nz\n")
    data_header.append("{} {}\n".format(header_nx,header_nz))
    #vp_min #vp_max #vs_min #vs_max #density_min #density_max   - min/max stats
    data_header.append("#vp_min #vp_max #vs_min #vs_max #density_min #density_max\n")
    data_header.append("{} {} {} {} {} {}\n".format(header_vp_min,header_vp_max,header_vs_min,header_vs_max,header_rho_min,header_rho_max))
    data_header.append("#\n")

    # data record info
    data_header.append("# data records - format:\n")
    # format: x #z #vp #vs #density
    data_header.append("#x #z #vp #vs #density\n")

    ## SPECFEM2D tomographic model format
    if use_tomography_file_ascii_format:
        ## file header
        print("  using ascii file format")
        print("")

        ## writes output file
        os.system("mkdir -p DATA")
        filename = "./DATA/tomography_model_marmousi2.xyz"
        with open(filename, "w") as fp:
            fp.write(''.join(data_header))
            fp.write(''.join(output_data))
            fp.close()

    else:
        # binary format
        print("  using binary file format")
        print("")
        
        # collects header info
        output_header_data = np.array([ header_origin_x, header_origin_z, header_end_x, header_end_z, \
                                        header_dx, header_dz, \
                                        header_nx, header_nz, \
                                        header_vp_min, header_vp_max, header_vs_min, header_vs_max, header_rho_min, header_rho_max ])

        ## writes output file (in binary format)
        os.system("mkdir -p DATA")
        filename = "./DATA/tomography_model_marmousi2.xyz.bin"
        with open(filename, "wb") as fp:
            # header
            print("    header:")
            write_binary_file_custom_real_array(fp,output_header_data)
            # model data
            print("    model data:")
            write_binary_file_custom_real_array(fp,output_data_array)
            fp.close()

        # use data header as additional metadata infos
        data_header.append("#\n")
        data_header.append("# data file: {}\n".format(filename))
        data_header.append("#\n")

        filename_info = "./DATA/tomography_model_marmousi2.xyz.bin.info"
        with open(filename_info, "w") as fp:
            fp.write(''.join(data_header))
            fp.close()

    print("tomography model written to: ",filename)
    print("")


#
#----------------------------------------------------------------------------
#


def write_binary_file_custom_real_array(file,data):
    """
    writes data array to binary file
    """
    # defines float 'f' or double precision 'd' for binary values
    custom_type = 'f'

    # gets array length in bytes
    # marker
    binlength = array.array('i')
    num_points = data.size
    if custom_type == 'f':
        # float (4 bytes) for single precision
        binlength.fromlist([num_points * 4])
    else:
        # double precision
        binlength.fromlist([num_points * 8])

    # user output
    print("    array length = ",binlength," Bytes")
    print("    number of points in array = ",num_points)
    print("    memory required: in (kB) {:.4f} / in (MB): {:.4f}".format(binlength[0] / 1024.0, binlength[0] / 1024.0 / 1024.0))
    print("")

    # writes array data
    binvalues = array.array(custom_type)

    data = np.reshape(data, (num_points), order='F') # fortran-like index ordering
    #print("debug: data ",data.tolist())

    binvalues.fromlist(data.tolist()) #.tolist())

    # fortran binary file: file starts with array-size again
    binlength.tofile(file)
    # data
    binvalues.tofile(file)
    # fortran binary file: file ends with array-size again
    binlength.tofile(file)

#
#----------------------------------------------------------------------------
#

def interpolate_marmousi(max_freq):
    """
    downloads Marmousi2 model and interpolates velocities onto a regular grid,
    with a grid spacing determined by the user parameters. Gaussian smoothing considers the target maximum frequency.
    """
    print("*****************************")
    print("Marmousi2 model interpolation")
    print("*****************************")

    # load marmousi2 model
    model = load_marmousi2()

    # modify original model
    if use_solid_domain_only or use_solid_domain_for_water_layer:
        model = modify_original_model(model)

    # smooth model
    if use_Gaussian_smoothing:
        model = smooth_model(max_freq,model)

    # resample to target grid
    new_model = resample_model(model)

    # creates tomography file(s)
    create_tomography_file(new_model)

#
#----------------------------------------------------------------------------
#

def usage():
    print("usage:")
    print("    ./interpolate_Marmousi_2_tomo.py [--smooth=maximum_frequency] [--no-water or --replace-water] [--binary]")
    print("  with")
    print("    --smooth=maximum_frequency  - (optional) maximum frequency (in Hz) to determine Gaussian smoothing width")
    print("    --no-water                  - (optional) removes top water layer from original model (shinks model size)")
    print("    --replace-water             - (optional) replaces top water layer from original model with solid velocities (keeps original model size)")
    print("    --binary                    - (optional) use a binary format for the tomography model output file")
    sys.exit(1)
#
#----------------------------------------------------------------------------
#

if __name__ == '__main__':
    max_freq = 0.0

    # reads arguments
    if len(sys.argv) <= 0: usage()
    i = 0
    for arg in sys.argv:
        i += 1
        #print("argument "+str(i)+": " + arg)
        # get arguments
        if "--help" in arg:
            usage()
        elif "--smooth=" in arg:
            str_array = arg.split('=')[1]
            use_Gaussian_smoothing = True
            max_freq = float(arg.split('=')[1])
        elif "--no-water" in arg:
            use_solid_domain_only = True
        elif "--replace-water" in arg:
            use_solid_domain_for_water_layer = True
        elif "--binary" in arg:
            use_tomography_file_ascii_format = False
        elif i >= 2:
            print("argument not recognized: ",arg)
            sys.exit(1)

    # main routine
    interpolate_marmousi(max_freq)


