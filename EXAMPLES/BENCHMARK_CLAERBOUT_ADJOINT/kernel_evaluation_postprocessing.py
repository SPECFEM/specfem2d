#!/usr/bin/env python
#
# perturbs model with a Gaussian perturbation
#
#########################################
from __future__ import print_function

import sys
import numpy as np
import array

from helper_functions import helper

##########################################################

## default setup
IN_DATA_FILES = './DATA/'
MODEL_FILES   = './MODELS/'
KERNEL_FILES  = './KERNELS/'

# misfits
use_data_sources    = True   # using data/synthetic traces for misfit calculation
use_adjoint_sources = False  # using adjoint traces for misfit calculation

# plotting kernel/model figures
do_plot_figures = False

##########################################################


def postprocessing(NSTEP,DT,NPROC,SIM_TYPE):
    """
    computes Claerbout's test < F* F dm, dm> =?= < F dm, F dm >
    """
    global do_plot_figures

    # user output
    print("postprocessing:")
    print("  NSTEP = ",NSTEP)
    print("  DT    = ",DT)
    print("  NPROC = ",NPROC)
    print("  type  = ",SIM_TYPE,"(1 == acoustic / 2 == elastic P-SV / 3 == elastic SH / 4 == coupled acoustic-elastic)")
    print("")

    # selects components
    if SIM_TYPE == 1:
        # acoustic
        comp_total = [0] # single trace pressure/potential
    elif SIM_TYPE == 2:
        # elastic
        # P-SV
        comp_total = [1,3]  # Ux/Uz
    elif SIM_TYPE == 3:
        # elastic
        # SH
        comp_total = [2]  # Uy
    elif SIM_TYPE == 4:
        # coupled acoustic-elastic (P-SV like with receiver(s) in elastic domain)
        comp_total = [1,3]
    else:
        print("Invalid SIM_TYPE {} for selecting component".format(SIM_TYPE))
        sys.exit(1)

    # helper functions
    hf = helper()

    MODELSPACE = 0.0
    DATASPACE  = 0.0
    DATASPACE_NEW = 0.0

    itrace_total = 0

    for myrank in range(0, NPROC):
        print("rank: ",myrank)

        # processors name
        prname = "proc{:06d}_".format(myrank)

        # coordinates
        print("")
        print("  reading model coordinates...")
        # global point arrays
        #  x_save(NGLLX,NGLLZ,nspec)
        filename = IN_DATA_FILES + prname + "x.bin"
        xstore = hf.read_binary_SEM_file(filename)
        #xstore = np.reshape(xstore, (NGLLX, NGLLZ, nspec), order='F') # fortran-like index ordering

        filename = IN_DATA_FILES + prname + "z.bin"
        zstore = hf.read_binary_SEM_file(filename)

        print("")
        print("  model dimensions: x min/max = ",xstore.min(),xstore.max())
        print("                    z min/max = ",zstore.min(),zstore.max())
        print("")

        # kernels
        print("")
        print("  reading kernels...")

        # K_rho_prime
        if SIM_TYPE == 1 or SIM_TYPE == 4:
            filename = KERNEL_FILES + prname + "rhop_acoustic_kernel.bin"
        else:
            filename = KERNEL_FILES + prname + "rhop_kernel.bin"
        krhop = hf.read_binary_SEM_file(filename)

        # coupled acoustic-elastic
        if SIM_TYPE == 4:
            # gets additional elastic kernel
            filename = KERNEL_FILES + prname + "rhop_kernel.bin"
            krhop_el = hf.read_binary_SEM_file(filename)
            # combines elastic & acoustic kernels for model update
            krhop += krhop_el

        # plotting
        if do_plot_figures:
            hf.plot_model_image(xstore,zstore,krhop,filename,plot_kernel=True)

        # K_alpha
        if SIM_TYPE == 1 or SIM_TYPE == 4:
            filename = KERNEL_FILES + prname + "c_acoustic_kernel.bin"
        else:
            filename = KERNEL_FILES + prname + "alpha_kernel.bin"
        kalpha = hf.read_binary_SEM_file(filename)

        # coupled acoustic-elastic
        if SIM_TYPE == 4:
            # gets additional elastic kernel
            filename = KERNEL_FILES + prname + "alpha_kernel.bin"
            kalpha_el = hf.read_binary_SEM_file(filename)
            # combines elastic & acoustic kernels for model update
            kalpha += kalpha_el

        # plotting
        if do_plot_figures:
            hf.plot_model_image(xstore,zstore,kalpha,filename,plot_kernel=True)

        # K_beta
        if SIM_TYPE == 1:
            # acoustic, zero shear for acoustic simulation
            kbeta = np.zeros((kalpha.size))
        else:
            filename = KERNEL_FILES + prname + "beta_kernel.bin"
            kbeta = hf.read_binary_SEM_file(filename)

        # plotting
        if do_plot_figures:
            hf.plot_model_image(xstore,zstore,kbeta,filename,plot_kernel=True)

        print("")
        print("  kernel rhop : min/max = ",krhop.min(),krhop.max())
        print("  kernel alpha: min/max = ",kalpha.min(),kalpha.max())
        if SIM_TYPE !=  1: print("  kernel beta : min/max = ",kbeta.min(),kbeta.max())
        print("")

        # integration weights
        print("  reading integration weights...")
        filename = KERNEL_FILES + prname + "weights_kernel.bin"
        weights = hf.read_binary_SEM_file(filename)

        # model
        print("")
        print("reading model...")
        print("")
        # target
        # rho
        filename = MODEL_FILES + 'target_model/' + prname + 'rho.bin'
        rho = hf.read_binary_SEM_file(filename)
        # vp
        filename = MODEL_FILES + 'target_model/' + prname + 'vp.bin'
        vp = hf.read_binary_SEM_file(filename)
        # vs
        if (SIM_TYPE == 1) :
            # acoustic, zero shear velocity
            vs = np.zeros((vp.size))
        else:
            filename = MODEL_FILES + 'target_model/' + prname + 'vs.bin'
            vs = hf.read_binary_SEM_file(filename)
        print("")
        print("  target model : rho min/max = ",rho.min(),rho.max())
        print("                 vp  min/max = ",vp.min(),vp.max())
        if SIM_TYPE != 1: print("                 vs  min/max = ",vs.min(),vs.max())
        print("")

        # initial
        # rho
        filename = MODEL_FILES + 'initial_model/' + prname + 'rho.bin'
        rho0 = hf.read_binary_SEM_file(filename)
        # vp
        filename = MODEL_FILES + 'initial_model/' + prname + 'vp.bin'
        vp0 = hf.read_binary_SEM_file(filename)
        # vs
        if (SIM_TYPE == 1) :
            # acoustic, zero shear velocity
            vs0 = np.zeros((vp.size))
        else:
            filename = MODEL_FILES + 'initial_model/' + prname + 'vs.bin'
            vs0 = hf.read_binary_SEM_file(filename)
        print("")
        print("  initial model: rho min/max = ",rho0.min(),rho0.max())
        print("                 vp  min/max = ",vp0.min(),vp0.max())
        if SIM_TYPE != 1: print("                 vs  min/max = ",vs0.min(),vs0.max())
        print("")

        # relative perturbations
        dlnrho = (rho - rho0)/rho0
        dlnvp = (vp - vp0)/vp0
        if SIM_TYPE == 1:
            # acoustic, has zero shear perturbation
            dlnvs = np.zeros(vs.size)
        elif SIM_TYPE == 4:
            # coupled acoustic-elastic vs models have zero values in acoustic domain
            dlnvs = np.zeros(vs.size)
            # takes only non-zero part (in elastic domain)
            idx = np.nonzero(vs0)
            dlnvs[idx] = (vs[idx] - vs0[idx])/vs0[idx]
        else:
            # relative shear perturbations
            dlnvs = (vs - vs0)/vs0

        print("")
        print("  relative model perturbation (rho - rho0)/rho0 : min/max = ",dlnrho.min(),dlnrho.max())
        print("  relative model perturbation (vp - vp0)/vp0    : min/max = ",dlnvp.min(),dlnvp.max())
        if SIM_TYPE != 1:
          print("  relative model perturbation (vs - vs0)/vs0    : min/max = ",dlnvs.min(),dlnvs.max())
        print("")

        ### calculate inner product in model space --- < F* F dm, dm>
        print("")
        print("  calculating inner product in model space: < F* F dm, dm>")
        print("")

        # F dm = S(m+dm) - S(m)
        # note: original statement
        #       we backpropogate syn-dat (see adj_seismogram.f90) => we have to add a minus sign in front of kernels
        if SIM_TYPE == 1:
            # acoustic
            MS_rhol = np.sum(weights * (-krhop) * dlnrho)
            MS_vpl  = np.sum(weights * (-kalpha) * dlnvp)
            MODELSPACE = MODELSPACE + MS_rhol + MS_vpl
            print("  model space contributions: Mrho = {:e} Mvp = {:e}".format(MS_rhol,MS_vpl))
        else:
            # elastic
            MS_rhol = np.sum(weights * (-krhop) * dlnrho)
            MS_vpl  = np.sum(weights * (-kalpha) * dlnvp)
            MS_vsl  = np.sum(weights * (-kbeta) * dlnvs)
            MODELSPACE = MODELSPACE + MS_rhol + MS_vpl + MS_vsl
            print("  model space contributions: Mrho = {:e} Mvp = {:e} Mvs = {:e}".format(MS_rhol,MS_vpl,MS_vsl))

        print("  total model space = {:e}".format(MODELSPACE))
        print("")

        ### calculate inner product in data space --- < F dm, F dm>
        print("")
        print("  calculating inner product in data space: < F dm, F dm>")
        print("")

        # data
        if use_data_sources:
            print("  data sources: ")
            for icomp  in comp_total:
                name = "file_single"
                if icomp == 0:
                    # acoustic
                    filename = "Up_" + name + "_p.su"     # pressure
                elif icomp == 1:
                    # elastic x-component
                    filename = "Ux_" + name + "_d.su"    # displacement
                elif icomp == 2:
                    # y-component (SH-case)
                    filename = "Uy_" + name + "_d.su"    # displacement
                elif icomp == 3:
                    # z-component (P-SV case)
                    filename = "Uz_" + name + "_d.su"    # displacement
                else:
                  print("Invalid component ",icomp)
                  sys.exit(1)

                # data d
                filename_dat = "OUTPUT_FILES.dat.forward/" + filename
                print("  ",filename_dat)
                print("")
                # reads in all traces
                dat = hf.read_SU_file(filename_dat)

                # synthetics s
                filename_syn = "OUTPUT_FILES.syn.forward/" + filename
                print("  ",filename_syn)
                print("")
                # reads in all traces
                syn = hf.read_SU_file(filename_syn)

                # updated synthetics s
                filename_syn_new = "OUTPUT_FILES.syn.updated/" + filename
                print("  ",filename_syn_new)
                print("")
                # reads in all traces
                syn_new = hf.read_SU_file(filename_syn_new)

                num_receivers = len(dat)
                num_samples = len(dat[0])

                print("")
                print("  num receivers = ",num_receivers)
                print("  samples       = ",num_samples)
                print("")

                # misfit (s - d)
                diff = syn - dat

                DS_l = 0.0
                for irec in range(num_receivers):
                    # single receiver trace
                    trace = diff[irec]

                    # inner product
                    DS_l = DS_l + np.sum(trace * trace) * DT

                    # total counter
                    itrace_total += 1

                DATASPACE = DATASPACE + DS_l
                print("  data space contribution = {:e}".format(DS_l))
                print("         total data space = {:e}".format(DATASPACE))
                print("")

                # updates misfit (s_new - d)
                diff_new = syn_new - dat
                DS_l = 0.0
                for irec in range(num_receivers):
                    # single receiver trace
                    trace = diff_new[irec]

                    # inner product
                    DS_l = DS_l + np.sum(trace * trace) * DT

                DATASPACE_NEW = DATASPACE_NEW + DS_l
                print("  data space contribution (updated) = {:e}".format(DS_l))
                print("         total data space (updated) = {:e}".format(DATASPACE_NEW))
                print("")


        # adjoint sources
        if use_adjoint_sources:
            print("  adjoint sources: ")
            for icomp  in comp_total:
                name = "file_single"
                if icomp == 0:
                    # acoustic
                    filename = "Up_" + name + ".su.adj"  # adjoint trace assumes Up_file_single.su name for acoustic simulations (pressure/potential type)
                elif icomp == 1:
                    # elastic
                    filename = "Ux_" + name + ".su.adj" # adjoint traces assumes Ux_file_single.su.adj name
                elif icomp == 2:
                    # y-component (SH-case)
                    filename = "Uy_" + name + ".su.adj" # adjoint traces assumes Uy_file_single.su.adj
                elif icomp == 3:
                    # z-component (P-SV case)
                    filename = "Ux_" + name + ".su.adj" # adjoint traces assumes Uz_file_single.su.adj
                else:
                  print("Invalid component ",icomp)
                  sys.exit(1)

                filename_adj = "SEM/" + filename
                print("  ",filename_adj)
                print("")

                # reads in all adjoint traces
                adj = hf.read_SU_file(filename_adj)

                num_receivers = len(adj)
                num_samples = len(adj[0])

                print("")
                print("  num receivers = ",num_receivers)
                print("  samples       = ",num_samples)
                print("")

                DS_l = 0.0
                for irec in range(num_receivers):
                    # single receiver trace
                    adj_trace = adj[irec]

                    # inner product
                    DS_l = DS_l + np.sum(adj_trace * adj_trace) * DT

                    # total counter
                    itrace_total += 1

                DATASPACE = DATASPACE + DS_l
                print("  data space contribution = ",DS_l)
                print("         total data space = ",DATASPACE)
                print("")

    print("")
    print("total number of traces considered: ",itrace_total)
    print("")
    print("DATASPACE      = {:e}".format(DATASPACE))
    print("DATASPACE_NEW  = {:e}".format(DATASPACE_NEW))
    print("MODELSPACE     = {:e}".format(MODELSPACE))
    print("")

    # data misfit decrease
    # relative error
    if DATASPACE != 0.0:
        error = (DATASPACE - DATASPACE_NEW) / DATASPACE
    else:
        error = DATASPACE_NEW
    print("relative data misfit (dataspace - dataspace_new)/dataspace = ",error)

    # checks decrease
    if DATASPACE_NEW > DATASPACE:
        print("data misfit increases for update model, please check kernels...")
        # exit with an error
        sys.exit(1)

    # at least a 20 % decrease should be achievable
    if abs(error) > 0.2:
        print("data misfit decrease looks ok")
    else:
        print("relative data misfit decrease should be larger, please check kernels...")
        # exit with an error
        sys.exit(1)

    # model vs. data space
    # relative error
    if DATASPACE != 0.0:
        error = (DATASPACE - MODELSPACE) / DATASPACE
    else:
        error = abs(MODELSPACE)
    print("")
    print("relative error       (dataspace - modelspace)/dataspace    = ",error)
    print("")
    print("total relative error should be small enough, < 0.15 should be OK")
    if abs(error) < 0.15:
        print("looks ok")
        print("")
    else:
        print("error seems too large, please check...")
        # exit with an error
        sys.exit(1)



#
#------------------------------------------------------------------------------------------
#

def usage():
    print("Usage: ./kernel_evaluation_postprocessing.py NSTEP DT NPROC SIM_TYPE(1 == acoustic / 2 == elastic P-SV / 3 == elastic SH / 4 == coupled acoustic-elastic)")
    print("")
    sys.exit(1)

#
#------------------------------------------------------------------------------------------
#

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) < 5:
        usage()

    ## input parameters
    NSTEP = int(sys.argv[1])
    DT = float(sys.argv[2])
    NPROC = int(sys.argv[3])
    SIM_TYPE = int(sys.argv[4])

    postprocessing(NSTEP,DT,NPROC,SIM_TYPE)
