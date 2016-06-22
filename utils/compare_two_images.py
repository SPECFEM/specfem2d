#!/usr/bin/env python
#
# compares two images

# taken from http://www.pyimagesearch.com/2014/09/15/python-compare-two-images in April 2015

# computes the Mean Squared Error (MSE) and the Structural Similarity Index
# (SSIM) between two images (having the exact same size)

# Slightly modified by Dimitri Komatitsch and Elliott Sales de Andrade, April 2015

#########################################################
# requires the "python-skimage" package to be installed
#########################################################

# In 3D comparing kernels i.e. volumetric stuff would be a bit difficult,
# however in 2D we can see them as images if we draw them and then we can use
# Python routines that compute a similarity index between two images.
#
# I took one from
# http://www.pyimagesearch.com/2014/09/15/python-compare-two-images/ and
# slightly modified it. It works well. You then just need to check if the
# Structural Similarity Index (SSIM) is greater than 0.99.
#
# For the 2D code we can apply this to EXAMPLES/LuoYang_fluid_solid_kernel for
# instance in the nightly builds.
#
# For the 3D code for now we could just draw pictures along a couple of cut
# planes and use that script as well. Better than nothing.

# import the necessary packages
from __future__ import (absolute_import, division, print_function)

import sys

try:
    import numpy as np
except:
    print("Error importing numpy, please install numpy package first...")
    sys.tracebacklimit=0
    raise Exception("Importing numpy failed")

try:
    from skimage.color import rgb2grey
except:
    print("Error importing rgb2grey from skimage.color, please install skimage package first...")
    sys.tracebacklimit=0
    raise Exception("Importing skimage.color failed")

try:
    from skimage.data import imread
except:
    print("Error importing imread from skimage.data, please install skimage package first...")
    sys.tracebacklimit=0
    raise Exception("Importing skimage.data failed")

try:
    import skimage
    version = skimage.__version__
    print("skimage version: ",version)
    v = version.split('.')
    # determines skimage version function
    new_version = 0
    if len(v) > 2:
        print("skimage structural similarity support: %i.%i.x" % (int(v[0]),int(v[1])))
        if int(v[1]) > 11:
            new_version = 1
        else:
            new_version = 0
    # imports structural similarity function
    if new_version == 1:
        # deprecated: from skimage.measure import structural_similarity as ssim
        print("importing new function compare_ssim as ssim")
        from skimage.measure import compare_ssim as ssim
    else:
        # older version <= 0.11.x
        print("importing old function structural_similarity as ssim")
        from skimage.measure import structural_similarity as ssim
except:
    print("Error importing structural similarity from skimage.measure, please install skimage package first...")
    sys.tracebacklimit=0
    raise Exception("Importing skimage.measure failed")

#####################################################################
# USER PARAMETERS

# tolerance values
TOL_SIM = 0.99

# verbosity
VERBOSE = False

#####################################################################

def mse(imageA, imageB):
    # the 'Mean Squared Error' between the two images is the
    # sum of the squared difference between the two images;
    # NOTE: the two images must have the same dimension
    err = np.sum((imageA.astype("float") - imageB.astype("float")) ** 2)
    err /= float(imageA.shape[0] * imageA.shape[1])

    # return the MSE, the lower the error, the more "similar"
    # the two images are
    return err


def compare_images(imageA, imageB, title, show_plot=True):
    """
    computes the mean squared error and structural similarity
    """

    # index values for mean squared error
    if VERBOSE: print("comparing mean squared error...")
    m = mse(imageA, imageB)

    # convert the images to grayscale
    if VERBOSE: print("converting to greyscale...")
    imageA_grey = rgb2grey(imageA)
    imageB_grey = rgb2grey(imageB)

    # uses image copies to avoid runtime warning for ssim computation
    img1_grey = np.copy(imageA_grey)
    img2_grey = np.copy(imageB_grey)

    # index values for structural similarity
    if VERBOSE: print("comparing structural similarity...")
    s = ssim(img1_grey, img2_grey)

    if show_plot:
        if VERBOSE: print("plotting images...")
        try:
            import matplotlib.pyplot as plt
        except:
            print("Error importing pyplot from matplotlib, please install matplotlib package first...")
            sys.tracebacklimit=0
            raise Exception("Importing matplotlib failed")

        # setup the figure
        fig, ax = plt.subplots(2, 2)
        fig.suptitle("%s\nMSE: %.5f, SSIM: %.5f" % (title, m, s))

        ax[0][0].text(-10, -10, 'MSE: %.5f' %(m))

        # show first image
        ax[0][0].imshow(imageA, cmap=plt.cm.gray)
        ax[0][0].axis("off")

        # show the second image
        ax[0][1].imshow(imageB, cmap=plt.cm.gray)
        ax[0][1].axis("off")

        ax[1][0].text(-10, -10, 'SSIM: %.5f' %(s))

        # show first grey image
        ax[1][0].imshow(img1_grey, cmap=plt.cm.gray)
        ax[1][0].axis("off")

        # show the second grey image
        ax[1][1].imshow(img2_grey, cmap=plt.cm.gray)
        ax[1][1].axis("off")

        # show the images
        plt.show()

    return m, s


def plot_image_comparison(image1,image2,show):
    """
    plots comparison between two images
    """
    print("comparing images:")
    print("  image 1 = %s" % image1)
    print("  image 2 = %s" % image2)
    print("")
    if show:
        print("  with image plots")
    else:
        print("  without image plots")
    print("")

    # load the images
    if VERBOSE: print("loading images...")
    imageA = imread(image1)
    imageB = imread(image2)

    # compare the images
    m, s = compare_images(imageA, imageB, "Image 1 vs. Image 2", show_plot=show)

    # user output
    print("")
    print("mean squared error   : values 0.0 perfect match")
    print("structural similarity: values 1.0 perfect, < %.2f poor similarity" % TOL_SIM)
    print("")
    print("mean squared error    = %f" % m)
    print("structural similarity = %f" % s)
    print("")
    print("result:")
    if s < TOL_SIM:
        # Failure
        print("  poor image similarity found")
        print("")
        sys.exit(1)
    else:
        # Success
        print("  good image similarity found")
        print("")
        sys.exit(0)


def usage():
    print("usage: ./compare_two_images.py image1 image2 (show)")
    print("  with")
    print("     image1,image2   - path to images (jpg,png) for comparison")
    print("     (optional) show - set to 1 to show image plots, otherwise only outputs comparison values")

if __name__ == '__main__':
    # initialize
    image1 = ''
    image2 = ''
    show = False

    # gets arguments
    if len(sys.argv) < 2 or len(sys.argv) > 4:
        usage()
        sys.exit(1)
    else:
        image1 = sys.argv[1]
        image2 = sys.argv[2]
        if len(sys.argv) == 4:
            if int(sys.argv[3]) == 1: show = True

    # compares images
    plot_image_comparison(image1,image2,show)


