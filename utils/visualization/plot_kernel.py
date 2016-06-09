#!/usr/bin/env python

import sys
import numpy as np

try:
    import matplotlib.pyplot as plt
except:
    print("Error importing pyplot from matplotlib, please install matplotlib package first...")
    sys.tracebacklimit=0
    raise Exception("Importing matplotlib failed")

def grid(x, y, z, resX=100, resY=100):
    """
    Converts 3 column data to matplotlib grid
    """
    from matplotlib.mlab import griddata
    #from scipy.interpolate import griddata

    xi = np.linspace(min(x), max(x), resX)
    yi = np.linspace(min(y), max(y), resY)

    # mlab version
    Z = griddata(x, y, z, xi, yi, interp='linear')
    # scipy version
    #Z = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

    X, Y = np.meshgrid(xi, yi)
    return X, Y, Z

def plot_kernels(filename):
    """
    plots ASCII kernel file
    """
    print "plotting kernel file: ",filename

    data = np.loadtxt(filename)

    # checks data
    if data.ndim != 2:
        print "Error: wrong data dimension for kernel file",data.ndim
        sys.tracebacklimit=0
        raise Exception("Invalid data dimension")

    # checks array
    if len(data[1,:]) != 5:
        print "data shape  : ",data.shape
        print "data lengths: ",len(data[:,1]),len(data[1,:])
        print "Error: wrong data format for kernel file",data.shape
        sys.tracebacklimit=0
        raise Exception("Invalid data format")

    # splits up data
    x = data[:,0]
    y = data[:,1]

    print "x-range min/max = %f / %f" % (x.min(), x.max())
    print "y-range min/max = %f / %f" % (y.min(), y.max())

    z1 = data[:,2] # e.g. rho
    z2 = data[:,3] # e.g. alpha
    z3 = data[:,4] # e.g. beta

    print "data 1: min/max = %e / %e" % (z1.min(),z1.max())
    print "data 2: min/max = %e / %e" % (z2.min(),z2.max())
    print "data 3: min/max = %e / %e" % (z3.min(),z3.max())

    total_max = abs(np.concatenate((z1,z2,z3))).max()
    print "data max = ",total_max

    total_max = 1.e-8

    # setup figure (with 3 subplots)
    fig, axes = plt.subplots(nrows=3, ncols=1)

    for i,ax in enumerate(axes.flat,start=1):
        # top
        if i == 1:
            X, Y, Z = grid(x,y,z1)
            ax.set_title("Kernels")
            ax.set_ylabel('K_rho')
        elif i == 2:
            X, Y, Z = grid(x,y,z2)
            ax.set_ylabel('K_alpha')
        elif i == 3:
            X, Y, Z = grid(x,y,z3)
            ax.set_ylabel('K_beta')

        im = ax.imshow(Z,vmax=total_max, vmin=-total_max,
                       extent=[x.min(), x.max(), y.min(), y.max()],cmap='jet') # cmap='RdBu'

    # moves plots together
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    # colorbar
    fig.colorbar(im, ax=axes.ravel().tolist())
    #fig.colorbar(im, ax=axes.ravel().tolist(),orientation='horizontal')

    # show the images
    plt.show()


def usage():
    print "usage: ./plot_kernel.py file"
    print "   where"
    print "       file - ASCII kernel file, e.g. OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat"

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
    else:
        file = sys.argv[1]

    plot_kernels(file)

