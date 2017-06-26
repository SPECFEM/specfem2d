#!/usr/bin/env python
#
# plot the cross-correlation and L2-norm between reference and output seismograms
#
import sys
import glob
import os
import numpy as np

# tolerance values
TOL_CORR = 0.8
TOL_ERR = 0.01

def plot_correlations(syn_file,ref_file):
    """
    plots correlation and L2-norm values between reference and output file
    """
    print('comparing kernels')
    print('  reference file: %s' % ref_file)
    print('  file          : %s\n' % syn_file)

    # makes sure files are both available
    if not os.path.isfile(ref_file):
        print "  file " + ref_file + " not found"
        sys.exit(1)
    if not os.path.isfile(syn_file):
        print "  file " + syn_file + " not found"
        sys.exit(1)

    corr_min = 1.0
    err_max = 0.0
    shift_max = 0.0

    # gets x-coordinates
    syn_x = np.loadtxt(syn_file)[:, 0]
    dx = syn_x[1] - syn_x[0]
    print "  dx size = ",dx
    print ""

    # gets y-coordinates
    #syn_y = np.loadtxt(syn_file)[:, 1]
    #dy = syn_y[1] - syn_y[0]
    #print "  dy size = ",dy

    # outputs table header
    print("|%-30s| %13s| %13s|" % ('kernel name', 'corr', 'err'))
    print("|------------------------------------------------------------|")

    # counter
    n = 0
    for i in range(1,4):
        # build reference and synthetics file names
        # specfem file: proc******_rhop_alpha_beta_kernel.dat
        fname = os.path.basename(syn_file)
        names = str.split(fname,"_")

        # trace
        kernel = names[i]

        # numpy: reads in file data
        # data format: e.g. #x-coord #y-coord #rhop #alpha #beta
        iker = i + 1
        ref = np.loadtxt(ref_file)[:, iker]
        syn = np.loadtxt(syn_file)[:, iker]

        # length warning
        if len(ref) != len(syn):
          print("Mismatch of file length in both files syn/ref = %d / %d" %(len(syn),len(ref)))
          sys.exit(1)

        # dx step size in reference file
        ref_x = np.loadtxt(ref_file)[:, 0]
        dx_ref = ref_x[1] - ref_x[0]
        # mismatch warning
        if abs(dx - dx_ref)/dx > 1.e-5:
          print("Mismatch of dx size in both files syn/ref = %e / %e" %(dx,dx_ref))
          sys.exit(1)

        # least square test
        norm = np.linalg.norm
        sqrt = np.sqrt

        # normalized by power in reference solution
        fac_norm = norm(ref)
        # or normalized by power in (ref*syn)
        #fac_norm = sqrt(norm(ref)*norm(syn))
        if fac_norm > 0.0:
            err = norm(ref-syn)/fac_norm
        else:
            err = norm(ref-syn)
        #print('norm syn = %e norm ref = %e' % (norm(syn),fac_norm))

        # correlation test
        # total length
        if fac_norm > 0.0:
            corr_mat = np.corrcoef(ref, syn)
        else:
            corr_mat = np.cov(ref-syn)
        corr = np.min(corr_mat)

        # statistics
        corr_min = min(corr, corr_min)
        err_max = max(err, err_max)

        # info string
        info = ""
        if corr < TOL_CORR:   info += "  poor correlation"
        if err > TOL_ERR:     info += "      poor match"

        # print results to screen
        print("|%-30s| %13.5f| %13.5le| %s" % (kernel, corr, err, info))

        # counter
        n += 1

    # print min(coor) max(err)
    print("|------------------------------------------------------------|")
    print("|%30s| %13.5f| %13.5le|" % ('min/max', corr_min, err_max))

    # output summary
    print("\nsummary:")
    print("%d kernels compared" % n)
    print("correlations: values 1.0 perfect, < %.1f poor correlation" % TOL_CORR)
    if corr_min < TOL_CORR:
        print("              poor correlation seismograms found")
    else:
        print("              no poor correlations found")
    print ""

    print("L2-error    : values 0.0 perfect, > %.2f poor match" % TOL_ERR)
    if err_max > TOL_ERR:
        print("              poor matching seismograms found")
    else:
        print("              no poor matches found")
    print ""

def usage():
    print "usage: ./compare_kernel_correlations.py kernel-file1 kernel-file2"
    print "  with"
    print "     kernel-file1 - ASCII kernel file,"
    print "                    OUTPUT_FILES/proc000000_rhop_alpha_beta_kernel.dat"
    print "     kernel-file1 - ASCII kernel file for reference"

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) != 3:
        usage()
        sys.exit(1)
    else:
        out_kernel = sys.argv[1]
        ref_kernel = sys.argv[2]

    plot_correlations(out_kernel,ref_kernel)

