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
KERNEL_FILES  = './KERNELS/'

# plotting model figures
do_plot_figures = True

##########################################################

#
#------------------------------------------------------------------------------------------
#

def model_update(NPROC,SIM_TYPE,percent):
    """
    update model using kernels
    """
    global do_plot_figures

    # helper functions
    hf = helper()

    # gradient step length factor
    step_factor = percent

    # loop over processes
    for myrank in range(0, NPROC):
        print("")
        print("reading in model files for rank {}...".format(myrank))
        print("")

        # processors name
        prname = IN_DATA_FILES + "proc{:06d}_".format(myrank)

        # note: we could deal with adding the Gaussian perturbation without having nspec and ibool arrays,
        #       and use just 1D-arrays read in from data files (since all positional xstore/zstore and
        #       model vpstore/.. arrays have the same dimensions (NGLLX,NGLLZ,nspec) -> 1D (NGLLX*NGLLZ*nspec) arrays)
        #       thus, we could just compute position and add Gaussian perturbations for each 1D-array entries.
        #
        #       nevertheless, having nspec & ibool might become useful for future handlings, so we'll use it
        #       here as an example use case.

        # mesh index
        #filename = prname + "NSPEC_ibool.bin"
        #nspec,ibool = hf.read_binary_NSPEC_ibool_file(filename=filename)
        #ibool = np.reshape(ibool, (NGLLX, NGLLZ, nspec), order='F') # fortran-like index ordering

        # coordinates
        # global point arrays
        #  x_save(NGLLX,NGLLZ,nspec)
        filename = prname + "x.bin"
        xstore = hf.read_binary_SEM_file(filename)

        filename = prname + "z.bin"
        zstore = hf.read_binary_SEM_file(filename)

        # user output
        print("model input directory: ",IN_DATA_FILES)
        print("")
        print("mesh:")
        print("  model dimension: x min/max = ",xstore.min(),xstore.max())
        print("                   z min/max = ",zstore.min(),zstore.max())
        print("")

        # model files
        print("initial model:")
        print("  directory: ",IN_DATA_FILES)

        # model rho,vp,vs
        filename = prname + "rho.bin"
        rhostore = hf.read_binary_SEM_file(filename)
        filename = prname + "vp.bin"
        vpstore = hf.read_binary_SEM_file(filename)
        filename = prname + "vs.bin"
        vsstore = hf.read_binary_SEM_file(filename)

        print("")
        print("  rho min/max = ",rhostore.min(),rhostore.max())
        print("  vp  min/max = ",vpstore.min(),vpstore.max())
        print("  vs  min/max = ",vsstore.min(),vsstore.max())

        # plot images
        if do_plot_figures:
            hf.plot_model_image(xstore,zstore,rhostore,prname + "rho",verbose=True)
            hf.plot_model_image(xstore,zstore,vpstore,prname + "vp",verbose=True)
            hf.plot_model_image(xstore,zstore,vsstore,prname + "vs",verbose=True)

        # kernels
        print("")
        print("kernels:")
        print("  input directory: ",KERNEL_FILES)

        prname = KERNEL_FILES + "proc{:06d}_".format(myrank)

        # K_rhop
        if SIM_TYPE == 1 or SIM_TYPE == 4:
            # acoustic domain kernel
            filename = prname + "rhop_acoustic_kernel.bin"
        else:
            # elastic domain kernel
            filename = prname + "rhop_kernel.bin"
        krhop = hf.read_binary_SEM_file(filename)

        # coupled acoustic-elastic
        if SIM_TYPE == 4:
            # gets additional elastic kernel
            filename = prname + "rhop_kernel.bin"
            krhop_el = hf.read_binary_SEM_file(filename)
            # combines elastic & acoustic kernels for model update
            krhop += krhop_el

        # plotting
        if do_plot_figures:
            hf.plot_model_image(xstore,zstore,krhop,filename,plot_kernel=True,verbose=True)


        # K_alpha
        if SIM_TYPE == 1 or SIM_TYPE == 4:
            # acoustic domain kernel
            filename = prname + "c_acoustic_kernel.bin"
        else:
            # elastic domain kernel
            filename = prname + "alpha_kernel.bin"
        kalpha = hf.read_binary_SEM_file(filename)

        # coupled acoustic-elastic
        if SIM_TYPE == 4:
            # gets additional elastic kernel
            filename = prname + "alpha_kernel.bin"
            kalpha_el = hf.read_binary_SEM_file(filename)
            # combines elastic & acoustic kernels for model update
            kalpha += kalpha_el

        # plotting
        if do_plot_figures:
            hf.plot_model_image(xstore,zstore,kalpha,filename,plot_kernel=True,verbose=True)

        # K_beta
        if SIM_TYPE == 1:
            # acoustic domain, zero shear for acoustic simulation
            kbeta = np.zeros((kalpha.size))
            print("  acoustic simulation: setting zero shear kernel")
        else:
            # elastic domain kernel
            filename = prname + "beta_kernel.bin"
            kbeta = hf.read_binary_SEM_file(filename)
            # plotting
            if do_plot_figures:
                hf.plot_model_image(xstore,zstore,kbeta,filename,plot_kernel=True,verbose=True)

        print("")
        print("  kernel rhop : min/max = ",krhop.min(),krhop.max())
        print("  kernel alpha: min/max = ",kalpha.min(),kalpha.max())
        print("  kernel beta : min/max = ",kbeta.min(),kbeta.max())
        print("")
        print("  norm of rhop kernel  = ",np.sum(krhop * krhop))
        print("  norm of alpha kernel = ",np.sum(kalpha * kalpha))
        print("  norm of beta kernel  = ",np.sum(kbeta * kbeta))

        # model update
        print("")
        print("model update:")

        # note: model update is using a gradient in the negative direction to decrease misfits
        #       this is a steepest descent update, using a step size of getting a maximum 1 percent update
        krhop_max = np.abs(krhop).max()
        kalpha_max = np.abs(kalpha).max()
        kbeta_max = np.abs(kbeta).max()

        max_val = max(krhop_max,kalpha_max,kbeta_max)
        if max_val != 0.0:
            step_size = step_factor / max_val
        else:
            print("Error: kernels are all zero, please check...")
            sys.exit(1)

        print("  kernel maximum value            = ",max_val)
        print("  maximum gradient step dln(m)    = ",step_factor)
        print("")
        print("  resulting step size             = ",step_size)
        print("")

        # gradients in negative direction
        drho   = - step_size * krhop
        dalpha = - step_size * kalpha
        dbeta  = - step_size * kbeta

        # update
        # kernels are for relative perturbations dln(m)
        #
        # dln(m) = G -> m_new = m * exp( dln(m) ) = m * exp( G )
        rhostore = rhostore * np.exp(drho)
        vpstore  = vpstore * np.exp(dalpha)
        vsstore  = vsstore * np.exp(dbeta)
        # or
        # to first order: dln(m) = dm/m = G -> dm = G * m
        #                                      m_new = m + dm = m + G * m = m (1 + G)
        #
        #rhostore = rhostore * (1.0 + drho)
        #vpstore  = vpstore * (1.0 + dalpha)
        #vsstore  = vsstore * (1.0 + dbeta)

        print("updated model:")
        print("  rho min/max = ",rhostore.min(),rhostore.max())
        print("  vp  min/max = ",vpstore.min(),vpstore.max())
        print("  vs  min/max = ",vsstore.min(),vsstore.max())
        print("")

        # stores perturbed model files
        print("storing new files...")
        prname = IN_DATA_FILES + "proc{:06d}_".format(myrank)
        # rho
        filename = prname + "rho_new.bin"
        hf.write_binary_file_custom_real_array(filename,rhostore)
        # vp
        filename = prname + "vp_new.bin"
        hf.write_binary_file_custom_real_array(filename,vpstore)
        # vs
        filename = prname + "vs_new.bin"
        hf.write_binary_file_custom_real_array(filename,vsstore)

        print("  updated model files written to:")
        print("  ",prname + "rho_new.bin")
        print("  ",prname + "vp_new.bin")
        print("  ",prname + "vs_new.bin")
        print("")

        # plot images
        if do_plot_figures:
            hf.plot_model_image(xstore,zstore,rhostore,prname + "rho_new",verbose=True)
            hf.plot_model_image(xstore,zstore,vpstore,prname + "vp_new",verbose=True)
            hf.plot_model_image(xstore,zstore,vsstore,prname + "vs_new",verbose=True)


#
#------------------------------------------------------------------------------------------
#

def usage():
    print("Usage: ./model_update.py NPROC SIM_TYPE percent")
    print("")
    sys.exit(1)

#
#------------------------------------------------------------------------------------------
#

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) < 4:
        usage()

    ## input parameters
    NPROC = int(sys.argv[1])
    SIM_TYPE = int(sys.argv[2])
    percent = float(sys.argv[3])

    print("")
    print("model update:")
    print("  NPROC        = ",NPROC)
    print("  SIM_TYPE     = ",SIM_TYPE)
    print("  percent      = ",percent)
    print("")

    # checks
    if percent <= 0.0:
        print("Invalid percent {} for step_factor, must be strictly positive".format(percent))
        sys.exit(1)

    model_update(NPROC,SIM_TYPE,percent)

    print("")
    print("all done")
    print("")
