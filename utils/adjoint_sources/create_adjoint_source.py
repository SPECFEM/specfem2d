#!/usr/bin/env python
# SPECFEM2D
#
# creates adjoint source for an amplitude kernel
# according to Tromp et al. (2005), eq. (75)
#
# requirements:
#   - numpy
#   - scipy
#   - matplotlib
#
from __future__ import absolute_import, division, print_function

import os
import sys

############################################################
## user parameters

# select time-domain taper:
#   cosine == 1 - default
#   boxcar == 2
#   Welch  == 3
TAPER_TYPE = 1

############################################################

# initializes global variables
type = 0              # adjoint source type
file_syn = ""         # synthetics file
file_dat = ""         # data file
show_figures = False  # show matplotlib figures


try:
    import numpy as np
except:
    print("Error importing numpy, please install numpy package first...")
    sys.tracebacklimit=0
    raise Exception("Importing numpy failed")

try:
    import matplotlib.pyplot as plt
except:
    print("Error importing pyplot from matplotlib, please install matplotlib package first...")
    sys.tracebacklimit=0
    raise Exception("Importing matplotlib failed")

try:
    import scipy
except:
    print("Error importing scipy, please install scipy package first...")
    sys.tracebacklimit=0
    raise Exception("Importing scipy failed")


def print_trace_stats(syn,dat,dt=None):
    """
    print stats about the traces: integrated waveform difference and power

    synthetics and data are 2-d arrays, with trace[:,0] time and trace[:,1] data values
    """
    length = len(syn[:,1])
    if length <= 1:
        print("Error integral measures: traces too short %i" % length)
        sys.exit(1)

    # time step
    if dt == None:
        dt = syn[1,0] - syn[0,0]

    # integrals
    m_syn = 0.5 * np.sum( syn[:,1]**2 ) * dt

    if not dat is None:
      m_dat = 0.5 * np.sum( dat[:,1]**2 ) * dt
      m_diff = 0.5 * np.sum( (dat[:,1]-syn[:,1])**2 ) * dt

    print("stats:")
    print("  trace length: %i" % len(syn[:,1]))
    print("  time step dt: %f" % dt)
    print("  syn min/max   = {} / {}".format(syn[:,1].min(),syn[:,1].max()))
    print("  syn integral  = {}".format(m_syn))

    if not dat is None:
        print("  dat min/max   = {} / {}".format(dat[:,1].min(),dat[:,1].max()))
        print("  dat integral = {}".format(m_dat))
        print("  difference integrated = {}".format(m_diff))
    print("")


def get_time_window_taper(trace,t_start,t_end):
    """
    defines the time window W(t) between [t_start,t_end]

    trace is a 2-d array, with trace[:,0] time and trace[:,1] data values
    """
    global TAPER_TYPE,show_figures

    # checks start/end
    if t_start > t_end:
        tmp = t_start
        t_start = t_end
        t_end = tmp

    # time window length
    length = len(trace[:,0])
    for i in range(length):
        if trace[i,0] < t_start: i_start = i
        if trace[i,0] <= t_end: i_end = i

    length_window = i_end - i_start + 1

    print("Taper:")
    print("  time start/end = %f / %f " % (t_start,t_end))
    print("  trace length   = %i" % length)
    print("  window length  = %i" % length_window)
    print("")

    TW = np.zeros(length)

    # creates time window
    if TAPER_TYPE == 1:
        # cosine taper
        # for time-domain cosine taper: 1 - [cos(t)]^(ipwr)
        power_taper = 10

        for i in range(i_start,i_end+1):
            TW[i] = 1.0 - np.cos(np.pi*(i-i_start)/length_window)**power_taper

    elif TAPER_TYPE == 2:
        # boxcar taper
        for i in range(i_start,i_end+1):
            TW[i] = 1.0

    elif TAPER_TYPE == 3:
        # Welch taper
        factor = (2.0/length_window)**2

        for i in range(i_start,i_end+1):
            TW[i] = 1.0 - factor*((i-i_start) - length_window/2.0)**2

    else:
        print("Error invalid TAPER_TYPE: %i - must be 1,2 or 3" % TAPER_TYPE)
        sys.exit(1)

    if show_figures:
        #plt.figure()
        plt.subplots(nrows=2, ncols=1)
        #plt.hold(True)
        plt.subplot(2,1,1)
        plt.title('Trace')
        plt.plot(trace[:,0],trace[:,1],'-g',label='Trace')

        plt.subplot(2,1,2)
        plt.title('Taper')
        plt.plot(trace[:,0],TW,'-b',label='Window taper')
        plt.legend()

        plt.tight_layout()
        plt.show()

    return TW


def create_adjoint_source(type,file_syn,file_dat):
    """
    takes two traces (data and synthetics) and computes the adjoint source for an attenuation kernel
    """
    global show_figures

    # initializes
    syn = None
    dat = None

    # user output
    print("adjoint source:")
    if type == 1:
        print("  type: ",type," - for cross-correlation traveltime kernel")
    elif type == 2:
        print("  type: ",type," - for amplitude (attenuation Q) kernel")
    else:
        print("Invalid adjoint source type, must be 1 or 2")
        sys.exit(1)

    print("  synthetics file : %s\n" % file_syn)
    if type == 2: print("  data file       : %s" % file_dat)

    # makes sure files are available
    if not os.path.isfile(file_syn):
        print("  file " + file_syn + " not found")
        sys.exit(1)
    if type == 2:
        if not os.path.isfile(file_dat):
            print("  file " + file_dat + " not found")
            sys.exit(1)

    # loads traces
    # file format:
    # #time #value
    # ..    ..
    #
    # synthetics file
    syn = np.loadtxt(file_syn)

    # data file
    if type == 2:
        dat = np.loadtxt(file_ref)

    # time
    syn_t = syn[:, 0]
    t0 = syn_t[0] # start
    t1 = syn_t[len(syn_t)-1] # end
    dt = syn_t[1] - syn_t[0]

    print("synthetics file:")
    print("  trace length = %i" % len(syn_t))
    print("  time step : %f" % dt)
    print("")
    print("  start time: %f" % t0)
    print("    end time: %f" % t1)
    print("")

    # more stats
    print_trace_stats(syn,dat,dt)

    # taper for windowing
    t_start = t0 + 0.1 * len(syn_t) * dt
    t_end   = t1 - 0.1 * len(syn_t) * dt

    taper = get_time_window_taper(syn,t_start,t_end)

    # tappered synthetic trace (data values)
    trace = taper[:] * syn[:,1]

    # adjoint source
    if type == 1:
        # traveltime adjoint source
        adj = np.zeros(len(trace))

        # gets time derivative to obtain velocity
        for i in range(1,len(trace)-1):
            # central finite-differences
            val = (trace[i+1] - trace[i-1]) / (2.0 * dt)
            # adding negative sign
            adj[i] = - val

        # 2nd-order scheme
        #for i in range(1,len(trace)-1):
        #    # central finite difference (2nd-order scheme)
        #    val = (trace[i+1] - 2.0 * trace[i] + trace[i-1]) / (dt**2)
        #    # adding negative sign
        #    adj[i] = - val

        # 4th-order scheme
        #for i in range(2,len(adj_trace)-2):
        #    # central finite difference (4th-order scheme)
        #    val = ( - 1.0/12.0 * trace[i+2] + 4.0/3.0 * trace[i+1] - 5.0/2.0 * trace[i]
        #            + 4.0/3.0 * trace[i-1] - 1.0/12.0 * trace[i-2] ) / (dt**2)
        #    # adding negative sign
        #    adj[i] = - val

    elif type == 2:
        # amplitude adjoint source
        # kernels following Tromp et al. (2005) eq.67
        adj = np.zeros(len(trace))

        # waveform
        adj = trace[:]

    else:
        print("Invalid adjoint source type, must be 1 or 2")
        sys.exit(1)

    # normalization factor
    print("normalization:")
    norm = dt * np.sum( adj[:] * adj[:] )
    print("  norm = ", norm)
    print("")

    # normalizes
    if np.abs(norm) > 0.0:
        adj /= norm
    else:
        # zero trace
        adj[:] = 0.0

    # plotting
    if show_figures:
        plt.title('Adjoint source')
        plt.plot(syn_t,adj,'-g',label='adjoint source')
        plt.plot(syn_t,trace,'-b',label='tappered trace')
        plt.legend()
        plt.show()

    # work directory
    path = os.getcwd()

    # output in SEM/ directory
    if not os.path.exists("SEM"): os.mkdir("SEM")
    os.chdir("SEM")

    # filename
    fname = os.path.basename(file_syn)
    names = fname.split(".")

    # remove file ending .semd and add .adj instead
    name = '.'.join(names[0:-1])
    filename = name + ".adj"

    # file output
    f = open(filename, "w")

    # file header
    #f.write("# adjoint source: %s\n" % (filename))

    # format: #time #adjoint_source_value
    for ii in range(len(adj)):
        time = syn_t[ii]
        val = adj[ii]
        f.write("{}\t{}\n".format(time,val))

    # numpy save
    #xy = np.empty(shape=(0,2))
    #for ii in range(len(adj)):
    #    time = syn_t[ii]
    #    xy = np.append(xy,[[time,adj[ii]]],axis=0)
    ## saves as ascii
    #np.savetxt(f, xy, fmt="%f\t%f")

    # closes file
    f.close()

    print("written to: ", "SEM" + "/" + filename)
    print("")

    # changes back to work directory
    os.chdir(path)

def usage():
    print("usage: ./create_adjoint_source.py type[1==Cross-Correlation/2==Amplitude] synthetics [data] [show]")
    print("   where")
    print("       type       - 1 == cross-correlation adjoint source / 2 == amplitude adjoint source")
    print("       synthetics - synthetics file, e.g. OUTPUT_FILES/AA.S0001.BXX.semd")
    print("       data       - (optional) reference data trace, e.g. REF_DATA/AA.S0001.BXX.semd, needed for amplitude adjoint sources")


if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)
    else:
        type = int(sys.argv[1])
        file_syn = sys.argv[2]

    # optional for amplitude adjoint sources
    if len(sys.argv) == 4:
        if sys.argv[3] == "show":
            show_figures = True
            file_dat = None
        else:
            file_dat = sys.argv[3]

    # optional for showing figures
    if len(sys.argv) == 5:
        show_figures = (sys.argv[4] == "show")

    # creates adjoint source files
    create_adjoint_source(type,file_syn,file_dat)


