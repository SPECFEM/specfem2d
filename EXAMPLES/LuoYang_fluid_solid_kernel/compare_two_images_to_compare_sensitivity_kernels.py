#! /usr/bin/python

from __future__ import (absolute_import, division, print_function)

# taken from http://www.pyimagesearch.com/2014/09/15/python-compare-two-images/
# in April 2015

# computes the Mean Squared Error (MSE) and the Structural Similarity Index
# (SSIM) between two images (having the exact same size)

# Slightly modified by Dimitri Komatitsch, April 2015

############################################################################
# requires the "python-skimage" and "python-opencv" packages to be installed
############################################################################

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
from skimage.measure import structural_similarity as ssim
import matplotlib.pyplot as plt
import numpy as np
import cv2


def mse(imageA, imageB):
    # the 'Mean Squared Error' between the two images is the
    # sum of the squared difference between the two images;
    # NOTE: the two images must have the same dimension
    err = np.sum((imageA.astype("float") - imageB.astype("float")) ** 2)
    err /= float(imageA.shape[0] * imageA.shape[1])

    # return the MSE, the lower the error, the more "similar"
    # the two images are
    return err


def compare_images(imageA, imageB, title):
    # compute the mean squared error and structural similarity
    # index for the images
    m = mse(imageA, imageB)
    s = ssim(imageA, imageB)

    # setup the figure
    fig = plt.figure(title)
    plt.suptitle("MSE: %.5f, SSIM: %.5f" % (m, s))

    # show first image
    ax = fig.add_subplot(1, 2, 1)
    plt.imshow(imageA, cmap=plt.cm.gray)
    plt.axis("off")

    # show the second image
    ax = fig.add_subplot(1, 2, 2)
    plt.imshow(imageB, cmap=plt.cm.gray)
    plt.axis("off")

    # show the images
    plt.show()


# load the images
original_reference = cv2.imread(
    "rho_and_kappa_kernels_acoustic_and_elastic_only_no_caption_reference.png")
image_from_new_calculation = cv2.imread(
    "rho_and_kappa_kernels_acoustic_and_elastic_only_no_caption.png")

# convert the images to grayscale
original_reference = cv2.cvtColor(original_reference, cv2.COLOR_BGR2GRAY)
image_from_new_calculation = cv2.cvtColor(image_from_new_calculation,
                                          cv2.COLOR_BGR2GRAY)

# initialize the figure
# fig = plt.figure("Images")
images = [("Original reference", original_reference),
          ("Image from new calculation", image_from_new_calculation)]

# loop over the images
# for (i, (name, image)) in enumerate(images):
#     # show the image
#     ax = fig.add_subplot(1, 3, i + 1)
#     ax.set_title(name)
#     plt.imshow(image, cmap = plt.cm.gray)
#     plt.axis("off")

# show the figure
# plt.show()

# compare the images
compare_images(original_reference, image_from_new_calculation,
               "Original reference vs. image from new calculation")
