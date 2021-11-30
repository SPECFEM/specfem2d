#!/usr/bin/env python
#
# plots a model defined as GLL or binary file
# e.g. DATA/proc000000_vp.bin
#
# requires x,z and model file
# e.g. DATA/proc000000_x.bin
#      DATA/proc000000_z.bin
#      DATA/proc000000_vp.bin
#
#########################################
from __future__ import print_function

import sys
import array

import numpy as np
import matplotlib.pyplot as plt

##########################################################

## default setup
IN_DATA_FILES = './DATA/'

# defines float 'f' or double precision 'd' for binary values
custom_type = 'f'

# verbosity
verbose = False

##########################################################


def read_marker(file,verbose=False):
    """
    reads array marker from fortran binary file
    (each time a fortran write(*)-routine call is executed, a marker at beginning and end is placed)
    """
    binlength = array.array('i')
    binlength.fromfile(file,1)

    if verbose:
        print("marker length = ",binlength[0])

    return binlength[0]

#
#------------------------------------------------------------------------------------------
#

def read_binary_file_custom_real_array(file,verbose=False):
    """
    reads data array from file
    """
    global custom_type

    # gets array length in bytes
    binlength = read_marker(file)

    if custom_type == 'f':
        # float (4 bytes) for single precision
        num_points = int(binlength / 4)
    else:
        # double precision
        num_points = int(binlength / 8)

    # user output
    if verbose:
        print("  array length = ",binlength," Bytes")
        print("  number of points in array = ",num_points)
        print("")

    # reads in array data
    binvalues = array.array(custom_type)
    binvalues.fromfile(file,num_points)

    # fortran binary file: file ends with array-size again
    # checks markers
    binlength_end = read_marker(file)
    if binlength_end != binlength:
        print("Error array markers in fortran binary file:",binlength,binlength_end)
        print("start array length = ",binlength," Bytes")
        print("final array length = ",binlength_end," Bytes")
        print("array lengths should match, please check your file")
        raise Exception('array markers invalid')

    data = list()
    data.append(binvalues)

    # returns data from list-output (first item)
    # converts list to numpy array
    data_array = np.array(data[0],dtype=custom_type)

    return data_array

#
#------------------------------------------------------------------------------------------
#

def read_binary_SEM_file(filename,verbose=False):
    # reads proc000000_x.bin,.. arrays
    global custom_type

    if verbose: print("binary file: ",filename)

    with open(filename,'rb') as f:
        # fortran file format: binary or gll
        # for example by:
        #   write(172) x_save
        data = read_binary_file_custom_real_array(f,verbose)

    return data

#
#------------------------------------------------------------------------------------------
#

def plot_model_image(xstore,zstore,array,name,plot_kernel=False,verbose=False,show=False):
    """
    plots model or kernel image
    """
    # 2d array
    # define grid.
    xi = np.linspace(xstore.min(), xstore.max(), 150)
    yi = np.linspace(zstore.min(), zstore.max(), 150)
    # grid data
    x = xstore.reshape((xstore.size))
    y = zstore.reshape((zstore.size))
    z = array.reshape((array.size))
    # Linearly interpolate the data (x, y) on a grid defined by (xi, yi).
    if 1 == 1:
        import matplotlib.tri as tri
        triang = tri.Triangulation(x, y)
        interpolator = tri.LinearTriInterpolator(triang, z)
        Xi, Yi = np.meshgrid(xi, yi)
        data = interpolator(Xi, Yi)
    else:
        # scipy
        from scipy.interpolate import griddata
        data = griddata((x,y), z, (xi[None, :], yi[:, None]), method='linear')

    #non-interpolated
    #dimx = int(xstore.max())
    #dimz = int(zstore.max())
    #data = np.zeros((dimx+1,dimz+1))
    #for ispec  in range(0, nspec):
    #    for j in range(0, NGLLZ):
    #        for i in range(0, NGLLX):
    #            # GLL point location (given in m: dimension 2640 m x 2640 x x 1440 m)
    #            iglob = ibool[i,j,ispec]
    #            x = xstore[i,j,ispec]
    #            z = zstore[i,j,ispec]
    #            ix = int(x)
    #            iz = int(z)
    #            data[ix,iz] = vpstore[i,j,ispec]

    extent = [x.min(),x.max(),y.min(),y.max()]

    # limit size
    total_max = abs(data).max()
    if verbose: print("plot: data max = ",total_max)

    plt.clf()
    plt.title(name)

    if plot_kernel:
        if total_max < 1.e-8:
            if total_max != 0.0:
                total_max = 1.0 * 10**(int(np.log10(total_max))-1)  # example: 2.73e-11 limits to 1.e-11
            else:
                total_max = 0.0
        if verbose: print("plot: color scale max = ",total_max)

        plt.imshow(data, extent=extent, origin='lower', interpolation='none', cmap="seismic",
                   vmax=total_max, vmin=-total_max)  # colormaps: "RdBu_r", "twilight"
    else:
        plt.imshow(data, extent=extent, origin='lower', interpolation='none', cmap="seismic")  # colormaps: "RdBu_r", "twilight"

    plt.colorbar()  # draw colorbar

    # save as JPEG file
    filename = name + ".jpg"
    plt.savefig(filename)
    print("  plotted as ",filename)

    # show plot
    if show: plt.show()

#
#------------------------------------------------------------------------------------------
#

def plot_model_from_gll_or_binary_file(filename,show=False,verbose=False):
    """
    reads binary/gll file together with positions x/z and plots model as pdf-image
    """
    global custom_type

    # note: simple binary files like proc000000_vp.bin, .. only store a single real array.
    #       reading this kind of file with numpy's fromfile() works just fine.
    #       for a more custom like reading of binary fortran files, we can also use routine read_binary_SEM_file()
    #
    # by default, we use the custom read routine as numpy's fromfile might misinterprete some initial marker values
    # coming from fortran-binary file formats
    USE_CUSTOM_READ = True

    # user output
    print("input file: ",filename)
    print("")

    # gets basename
    # e.g. proc000000_vp.bin -> proc000000_vp
    names = filename.split(".")
    # checks ending
    if names[-1] != "bin":
        print("must have a binary file **.bin as input")
        sys.exit(1)

    basename = ''.join(names[:len(names)-1])

    # gets processors name from input filename
    # e.g. proc000000_vp.bin -> proc000000_
    names = filename.split("_")
    if "proc" in names[0]:
        prname = names[0] + "_"
    else:
        prname = ''.join(names[:len(names)-1]) + "_"

    # reads in model file
    # model rho,vp,vs
    # for example: DATA/proc000000_vp.bin
    print("reading in file: ",filename)

    # simple
    if USE_CUSTOM_READ:
        model_array = read_binary_SEM_file(filename,verbose)
    else:
        model_array = np.fromfile(filename,dtype=custom_type)

    print("")
    print("model: ",basename)
    print("  values min/max = ",model_array.min(),"/",model_array.max())
    print("")

    # gets coordinates
    # for example: DATA/proc000000_x.bin
    #              DATA/proc000000_z.bin
    print("reading in position files...")

    # global point arrays
    #  x_save(NGLLX,NGLLZ,nspec)
    filename = prname + "x.bin"
    print(filename)
    if USE_CUSTOM_READ:
        xstore = read_binary_SEM_file(filename,verbose)
    else:
        xstore = np.fromfile(filename,dtype=custom_type)

    filename = prname + "z.bin"
    print(filename)
    if USE_CUSTOM_READ:
        zstore = read_binary_SEM_file(filename,verbose)
    else:
        zstore = np.fromfile(filename,dtype='float32')

    print("")
    print("mesh dimensions: x min/max = ",xstore.min(),"/",xstore.max())
    print("                 z min/max = ",zstore.min(),"/",zstore.max())
    print("")

    # plot image
    if "kernel" in basename:
        # kernel plot
        print("plotting kernel...")
        plot_model_image(xstore,zstore,model_array,basename,plot_kernel=True,show=show)
    else:
        # model plot
        print("plotting model...")
        plot_model_image(xstore,zstore,model_array,basename,show=show)

    print("")
    print("done")


#
#------------------------------------------------------------------------------------------
#

def usage():
    print("Usage: ./plot_model_from_gll_or_binary_file.py filename [show]")
    print("  with")
    print("    filename - file path, e.g. DATA/proc000000_vp.bin")
    print("               (requires to have DATA/proc000000_x.bin,DATA/proc000000_z.bin available as well)")
    print("    show     - (optional) show matplot image, otherwise will just be safed as .pdf file")
    sys.exit(1)

#
#------------------------------------------------------------------------------------------
#

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) < 2:
        usage()

    ## input parameters
    filename = sys.argv[1]

    if len(sys.argv) >= 3:
        show = True
    else:
        show = False

    plot_model_from_gll_or_binary_file(filename,show=show,verbose=verbose)
