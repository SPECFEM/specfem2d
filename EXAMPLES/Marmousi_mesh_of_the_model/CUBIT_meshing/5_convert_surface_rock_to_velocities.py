#!/usr/bin/env python
#
# reads in surface rock definitions from Capdeville's original file marmousi2.mat,
# using format:
#          #type #rock #k #v0 #vs #rho #unused #unused
#  with
#          type  - S or P for surface or point definition
#          rock  - rock type
#          k     - compaction gradient
#          v0    - compressional velocity, initial velocity (vp)
#          vs    - shear velocity, if negative then derived from table 2.3
#          rho   - density, if negative then derived from table 2.3
#
# rock types can be 1 = sand, 2 = shale, 3 = limestone (not used in table 2.8), and 4 = none (vp, vs, rho given directly)
#
# -----------------
# Marmousi 2 model:
# -----------------
#
# see paper by G. Martin et al., "Marmousi2: An elastic upgrade for Marmousi", The leading edge, 2006.
#
# note:
# - uses compaction gradient k to calculate velocity:
#     v = v_0_new + k * z,
#   with z being depth
#
# - initial velocity v_0_new is the "corrected" value from the original Marmousi model:
#      v_0_new = v_0_old - (k * 468)
#   with v_0_old being the original Marmousi vp value
#
# Table 2.3 - rock types and conversion formula
#  water    : vp = 1500         vs =    0                              rho = 1.01
#  sand     : vp from marmousi  vs = 0.804 * vp - 856                  rho = 0.2736 * vp**(0.261)
#  shale    : vp from marmousi  vs = 0.770 * vp - 867                  rho = 0.2806 * vp**(0.265)
#  limestone: vp from marmousi  vs = 1.017 * vp - 0.055* vp**2 - 1030  rho = 0.3170 * vp**(0.225)
#  salt     : vp = 4500         vs = 2600                              rho = 2.14
#
# Table 2.8 - Marmousi2 horizons, layers, elastic properties
#   Marmousi 2 specifies 192 layers, between horizons (top/bottom),
#   and corresponding rock type, fluid, compaction gradient, initial velocity v_0, vs, rho
#
# The cubit mesh from Capdeville/Hom Nath contains 435 surfaces.
# For each of these surfaces, this script will assign material properties taken from Capdeville's file marmousi2.mat.
#
from __future__ import print_function
import sys
import os

###############################################################################
# USER PARAMETER

# material_file
data_file = "../original_mesh_from_Yann_Capdeville_2009_using_Gmsh/marmousi2.mat"

## modifications
# note: the Marmousi2 model shown in Capdeville et al. (2010) replaces the original water layer with top solid velocities.
#       we can choose here to either use the original water layer or replace it with Capdeville's elastic part.
#
# instead of water layer, convert it to a solid
use_solid_domain_for_water_layer = False


###############################################################################

def convert_to_velocity():
    global data_file

    print("convert Marmousi2 surface rock definitions to velocities")

    # opens file
    print("data file: ",data_file,"\n")
    if not os.path.isfile(data_file):
        print("File not found, please check file: ",data_file)
        sys.exit(1)

    try:
        f_in = open(data_file,'r')
    except:
        print("error opening file ",data_file)
        sys.tracebacklimit=0
        raise Exception('file does not open: %s' % data_file)

    # get header
    txt = f_in.readline()
    print("header            : ",txt.strip())

    # get number of surfaces
    txt = f_in.readline()
    ns = int(txt)
    print("number of surfaces: ",ns)
    print("")

    # check
    if ns != 435:
        print("Error wrong number of surfaces, must be 435! please check...")
        sys.exit(1)

    # writes out to nummaterial_velocity_file_marmousi2
    filename = "./MESH/nummaterial_velocity_file_marmousi2"

    print("creating output file: ",filename)
    print("")

    os.system('mkdir -p MESH')
    try:
        f_out = open(filename,'w')
    except:
        print("error opening file ",filename)
        sys.tracebacklimit=0
        raise Exception('file does not open: %s' % filename)

    iline = 0
    for i in range(435):
        iline += 1

        # get data line
        txt = f_in.readline()
        dat = txt.split() # splits on whitespace

        # checks number of items
        if len(dat) != 8:
            print("Error line format: ",txt)
            sys.exit(1)

        # get items
        id = dat[0]
        rock_type = int(dat[1])    # 1==sand / 2==shale / 3==limestone / 4==user_defined
        k = float(dat[2])          # compaction gradient
        v0 = float(dat[3])         # P-velocity
        v1 = float(dat[4])         # Vs      for user_defined rock
        v2 = float(dat[5])         # Density for user_defined
        #dummy = dat[6]
        #dummy = dat[7]

        # identifier
        if id != "S":
            print("Error wrong surface type: ",txt)
            sys.exit(1)

        # determines elastic properties
        if rock_type == 1:
            # sand
            # (see Table 2.3 from Martin MSc thesis, page 9)
            vp = v0
            vs = 0.804 * vp - 856.0
            rho = 0.2736 * vp**(0.261); rho = 1000.0 * rho
            # check
            if vs <= 0.0:
                # note: v0 is corrected by v0 = v0_old - k * 468
                #       this can lead to negative vs if depth is not considered (v = v0 + k * z)
                #
                # since we have no depth information of the surface, we then take the old value for vp, vs, etc.
                # that is: v0 = v0_old - k * 468 -> v0_old = v0 + k * 468
                vp = v0 + k * 468.0
                vs = 0.804 * vp - 856.0
                rho = 0.2736 * vp**(0.261); rho = 1000.0 * rho

        elif rock_type == 2:
            # shale
            # (see Table 2.3 from Martin MSc thesis, page 9)
            vp = v0
            vs = 0.770 * vp - 867.0
            rho = 0.2806 * vp**(0.265); rho = 1000.0 * rho
            # check
            if vs < 0.0:
                # note: v0 is corrected by v0 = v0_old - k * 468
                #       this can lead to negative vs if depth is not considered (v = v0 + k * z)
                #
                # since we have no depth information of the surface, we then take the old value for vp, vs, etc.
                # that is: v0 = v0_old - k * 468 -> v0_old = v0 + k * 468
                vp = v0 + k * 468.0
                vs = 0.770 * vp - 867.0
                rho = 0.2806 * vp**(0.265); rho = 1000.0 * rho

        elif rock_type == 3:
            # limestone
            print("Error limestone not defined in Table 2.8, please check data line:",iline," -- ",txt)
            sys.exit(1)

        elif rock_type == 4:
            # as defined by entry
            vp = v0
            vs = v1
            rho = v2

        else:
            print("Error rock type undefined: ",txt)
            sys.exit(1)

        # check
        if vs < 0.0:
            print("Error invalid shear velocity ",txt," --- vp/vs/rho = ",vp,vs,rho)
            sys.exit(1)

        # note: the original Marmousi2 model gets modified on top by
        #
        #       Capdeville et al. (2010),
        #       2-D non-periodic homogenization to upscale elastic media for P–SV waves,
        #       GJI, 182, 903-922.
        #
        #       "For the original Marmousi and Marmousi2 models, the top layer is a water layer corresponding to the ocean.
        #       We replace this layer by an elastic layer with the same P-wave velocity but a non-zero S-wave velocity.
        #       The reason for this modification is to avoid the occurrence of a solid–fluid interface and the associated boundary layer
        #       from the point of view of homogenization which we shall present later. (..) "


        # surface 240 is water layer
        if iline == 240:
            if use_solid_domain_for_water_layer:
                # keep using the elastic velocities read in from Capdeville's file
                print("  replacing water layer with solid velocities (from Capdeville et al. 2010)")
                print("")
                pass
            else:
                # sets layer to water again (see table 2.3 from Martin MSc thesis, page 9)
                print("  using water layer (from original Marmousi2 model)")
                print("")
                vs = 0.0
                vp = 1500.0
                k = 0.0
                rho = 1.01; rho = 1000.0 * rho     # converts g/cm^3 to kg/m^3

        # determine acoustic==1 or elastic==2 domain
        domain_id = 1 if vs == 0.0 else 2

        # outputs properties
        #
        # format of nummaterial_velocity_file must be:
        #     #(1)material_domain_id #(2)material_id  #(3)rho  #(4)vp   #(5)vs   #(6)Q_kappa   #(7)Q_mu  #(8)anisotropy_flag
        # with material_domain_id ==1 for acoustic or ==2 for elastic materials
        #
        # available material types
        #   acoustic:    1 mat_id rho Vp 0  QKappa Qmu 0
        #   elastic:     2 mat_id rho Vp Vs QKappa Qmu 0
        #
        if iline == 1:
            # header
            f_out.write("# nummaterial_velocity_file for Marmousi2 - created by script convert_surface_rock_to_velocities.py\n")
            f_out.write("#\n")
            f_out.write("# velocity and density models:\n")
            f_out.write("#   nbmodels  = {}\n".format(ns))
            if use_solid_domain_for_water_layer:
                f_out.write("#   replacing water layer with solid velocities (material_id==240)\n".format(ns))
            else:
                f_out.write("#   using water layer (material_id==240)\n".format(ns))
            f_out.write("#\n")
            f_out.write("# format:\n")
            f_out.write("# (instead of anisotropy flag, added compaction gradient k)\n")
            f_out.write("#domain_id #material_id #rho #Vp #Vs #QKappa #Qmu #k(instead of ani==0)\n")
            f_out.write("#\n")
        # entry
        f_out.write("%d %3d %10.4f %10.4f %10.4f 9999.0 9999.0 %f\n" % (domain_id, iline, rho, vp, vs, k))

    f_in.close()
    f_out.close()

    print("written to: ",filename)
    print("")
    print("now, consider replacing the default `MESH/nummaterial_velocity_file` with this marmousi2 file,")
    print("or in the `Par_file` change the filename entry `nummaterial_velocity_file = ..` accordingly." )
    print("")
    print("all done")
    print("")

#
#----------------------------------------------------------------------------
#

def usage():
    print("usage:")
    print("    ./convert_surface_rock_to_velocities.py [--with-water / --without-water]")
    print("  with")
    print("    --with-water             - uses top water layer from original model (default)")
    print("    --without-water          - replaces top water layer from original model with solid velocities (as in Capdeville et al. 2010)")
    sys.exit(1)
#
#----------------------------------------------------------------------------
#

if __name__ == '__main__':

    # reads arguments
    if len(sys.argv) <= 0: usage()
    i = 0
    for arg in sys.argv:
        i += 1
        #print("argument "+str(i)+": " + arg)
        # get arguments
        if "--help" in arg:
            usage()
        elif "--with-water" in arg:
            use_solid_domain_for_water_layer = False
        elif "--without-water" in arg:
            use_solid_domain_for_water_layer = True
        elif i >= 2:
            print("argument not recognized: ",arg)
            sys.exit(1)

    # main routine
    convert_to_velocity()


