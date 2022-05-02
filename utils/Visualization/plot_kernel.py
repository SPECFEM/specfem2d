#!/usr/bin/env python
# SPECFEM2D
#
# script plots kernels and saves image as PNG-file
#
# usage:
#   ./plot_kernel.py file [1 == show figure / 0 == just plot file]
#
# for example:
#   ./plot_kernel.py OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat 1
#
#
# requirements:
#   - numpy
#   - matplotlib
#
from __future__ import print_function
import os.path
import sys
import numpy as np

try:
    import matplotlib.pyplot as plt
except:
    print("Error importing pyplot from matplotlib, please install matplotlib package first...")
    sys.tracebacklimit=0
    raise Exception("Importing matplotlib failed")

##################################################

## set default colormap
#colormap = 'jet'
colormap = 'RdBu'

##################################################


def grid(x, y, z, resX=100, resY=100):
    """
    Converts 3 column data to matplotlib grid
    """
    #griddata deprecated since matplotlib 3.1: from matplotlib.mlab import griddata
    from scipy.interpolate import griddata

    xi = np.linspace(min(x), max(x), resX)
    yi = np.linspace(min(y), max(y), resY)

    # mlab version
    #Z = griddata(x, y, z, xi, yi, interp='linear')
    # scipy version
    Z = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

    X, Y = np.meshgrid(xi, yi)
    return X, Y, Z

def plot_kernels(filename,show=False,max_val=None):
    """
    plots ASCII kernel file
    """
    global colormap

    print("plotting kernel file: ",filename)
    print("")

    data = np.loadtxt(filename)

    # checks data
    if data.ndim != 2:
        print("Error: wrong data dimension for kernel file",data.ndim)
        sys.tracebacklimit=0
        raise Exception("Invalid data dimension")

    # number of columns
    num_col = len(data[1,:])
    # number of kernels
    num_kernels = num_col - 2
    # number of points
    num_p = len(data[:,1])

    # checks array
    if num_kernels != 2 and num_kernels != 3:
        print("data shape  : ",data.shape)
        print("data lengths: ",len(data[:,1]),len(data[1,:]))
        print("Error: wrong data format for kernel file",data.shape)
        sys.tracebacklimit=0
        raise Exception("Invalid data format")

    # splits up data
    x = data[:,0]
    y = data[:,1]

    print("dimensions:")
    print("  x-range min/max = %f / %f" % (x.min(), x.max()))
    print("  y-range min/max = %f / %f" % (y.min(), y.max()))
    print("")

    z1 = np.zeros(num_p,dtype='float')
    z2 = np.zeros(num_p,dtype='float')
    z3 = np.zeros(num_p,dtype='float')

    if num_kernels == 2:
        z1 = data[:,2] # e.g. rho
        z2 = data[:,3] # e.g. bulk_c
    else:
        z1 = data[:,2] # e.g. rho
        z2 = data[:,3] # e.g. alpha
        z3 = data[:,4] # e.g. beta

    # names like
    #   proc000000_rhop_c_kernel.dat
    # or
    #   proc000000_rhop_alpha_beta_kernel.dat
    # or
    #   rhop_alpha_beta_kernel.dat
    name = os.path.basename(file)

    name_kernels = str.split(name,"_")
    num_k = len(name_kernels)
    istart = 0
    if "proc" in name_kernels[0]:
        istart = 1
        num_k = num_k - 1
    if "kernel.dat" in name_kernels[len(name_kernels)-1]:
        num_k = num_k - 1

    # check
    if num_k != num_kernels:
        print("warning: name has ",num_k,",but file columns provide ",num_kernels,"kernels")

    if num_k == 2:
        # rhop_c
        kernel1 = 'K_' + name_kernels[istart] # rhop
        kernel2 = 'K_' + name_kernels[istart+1] # c
        kernel3 = 'K_zero'
    elif num_k == 3:
        # rhop_alpha_beta
        kernel1 = 'K_' + name_kernels[istart] # rhop
        kernel2 = 'K_' + name_kernels[istart+1] # alpha
        kernel3 = 'K_' + name_kernels[istart+2] # beta
    else:
        kernel1 = 'K_1'
        kernel2 = 'K_2'
        kernel3 = 'K_3'

    print("statistics:")
    print("  %12s : min/max = %e / %e" % (kernel1,z1.min(),z1.max()))
    print("  %12s : min/max = %e / %e" % (kernel2,z2.min(),z2.max()))
    if num_kernels == 3:
        print("  %12s : min/max = %e / %e" % (kernel3,z3.min(),z3.max()))
    print("")

    total_max = abs(np.concatenate((z1,z2,z3))).max()
    print("  data max = ",total_max)
    print("")

    # plotting limits
    if not max_val is None:
        # user defined
        total_max = max_val
    else:
        # limit size to data range
        if total_max < 1.e-5:
            total_max = 1.0 * 10**(int(np.log10(total_max))-1)  # example: 2.73e-11 limits to 1.e-11
        else:
            total_max = 1.e-5
    print("plot: color scale max = ",total_max)
    print("")

    # determine orientation
    plot_horizontal = False

    x_range =  x.max() - x.min()
    y_range =  y.max() - y.min()
    if abs(y_range) > abs(x_range): plot_horizontal = True

    # setup figure
    if num_kernels == 2:
        # with 2 subplots
        if plot_horizontal:
            fig, axes = plt.subplots(nrows=1, ncols=2)
        else:
            fig, axes = plt.subplots(nrows=2, ncols=1)
    else:
        # with 3 subplots
        if plot_horizontal:
            fig, axes = plt.subplots(nrows=1, ncols=3)
        else:
            fig, axes = plt.subplots(nrows=3, ncols=1)

    for i,ax in enumerate(axes.flat,start=1):
        # top
        if i == 1:
            X, Y, Z = grid(x,y,z1)
            ax.set_title("Kernels")
            ax.set_ylabel(kernel1)
        elif i == 2:
            X, Y, Z = grid(x,y,z2)
            ax.set_ylabel(kernel2)
        elif i == 3:
            X, Y, Z = grid(x,y,z3)
            ax.set_ylabel(kernel3)

        im = ax.imshow(Z,vmax=total_max, vmin=-total_max,
                       extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap=colormap)

    # moves plots together
    if plot_horizontal:
        fig.subplots_adjust(wspace=0)
        plt.tight_layout()
    else:
        fig.subplots_adjust(hspace=0)

    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    # colorbar
    if plot_horizontal:
        fig.colorbar(im, ax=axes.ravel().tolist(),orientation='horizontal')
    else:
        fig.colorbar(im, ax=axes.ravel().tolist())

    # show the figure
    if show:
        plt.show()

    # saves kernel figure as file
    dir = os.path.dirname(file)
    name_without_ending = str.split(name,".")[0]
    outfile = dir + "/" + name_without_ending + ".png"
    fig.savefig(outfile, format="png")

    print("*****")
    print("plotted file: ",outfile)
    print("*****")
    print("")


def usage():
    print("usage: ./plot_kernel.py file [1 == show figure / 0 == just plot file] [max_val]")
    print("   where")
    print("       file    - ASCII kernel file, e.g. OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat")
    print("       1 or 0  - (optional) to show figure use 1, or 0 to just plot file (default to just plot file)")
    print("       max_val - (optional) to use max_val for kernel plot limit")

if __name__ == '__main__':
    # initializes
    show = False
    max_val = None

    # gets arguments
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
    else:
        file = sys.argv[1]

    if len(sys.argv) > 2:
        if sys.argv[2] == "show" or sys.argv[2] == "1": show = True

    if len(sys.argv) == 4:
        max_val = float(sys.argv[3])

    plot_kernels(file,show,max_val)

