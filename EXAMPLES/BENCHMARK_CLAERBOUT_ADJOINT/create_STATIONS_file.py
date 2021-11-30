#!/usr/bin/env python
#
# creates adjoint sources f^adj = (syn - data)
#
#########################################
from __future__ import print_function

import sys
import os.path
import math

#########################################

# model dimensions from Par_file
dim_x = 1000.0
dim_z = 1000.0

# start position
x0 = 700.0
z0 = 500.0

#########################################


def create_STATIONS_file(nlinesx,nlinesz):
    """
    creates a STATIONS files in DATA/
    """
    global dim_x,dim_z,x0,z0

    # user output
    print("creating STATIONS file:")
    print("  nlinesx = ",nlinesx)
    print("  nlinesz = ",nlinesz)
    print("")

    # checks if DATA/ is there
    if not os.path.isdir("DATA/"):
        print("DATA/ directory does not exist, please check...")
        sys.exit(1)

    # opens file
    filename = "DATA/STATIONS"
    try:
        f = open(filename,'w')
    except:
        print("Error opening file ",filename)
        sys.tracebacklimit=0
        raise Exception("file does not open: %s" % filename)

    # position increments
    if nlinesx > 1:
        dx = dim_x / (nlinesx-1)
    else:
        dx = 0.0

    if nlinesz > 1:
        dz = dim_z / (nlinesz-1)
    else:
        dz = 0.0

    # string formatting
    idxlength = int(math.ceil(math.log10(nlinesx*nlinesz)))

    # creates nlinesx * nlinez receivers
    irec = 0
    for iz in range(nlinesz):
        for ix in range(nlinesx):
            # receiver counter
            irec += 1

            # position
            x = x0 + ix * dx
            z = z0 + iz * dz

            # name
            sta = "S{0:0{1}d}".format(irec,idxlength)
            #sta = "S{:06d}".format(irec)

            # network
            net = "AA"

            # STATIONS
            # format:
            #   #name #network #x #z #elevation #burial_depth
            f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sta,net,x,z,0.0,0.0))

    f.close()

    # user output
    print("  total number of stations: ",nlinesx * nlinesz)
    print("")
    print("  stations written to: ",filename)
    print("")

#
#------------------------------------------------------------------------------------------
#

def usage():
    #print("Usage: ./adj_seismogram.py NSTEP DT NPROC type[acoustic/elastic]")
    print("Usage: ./create_STATIONS_file.py nlinesx nlinesz")
    print("")
    sys.exit(1)


#
#------------------------------------------------------------------------------------------
#

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) < 3:
        usage()

    nlinesx = int(sys.argv[1])
    nlinesz = int(sys.argv[2])

    create_STATIONS_file(nlinesx,nlinesz)
