#!/usr/bin/env python
#
# compares kernel image against reference image
# (needs script utils/compare_two_images.py)

from __future__ import (absolute_import, division, print_function)

import sys

try:
    from compare_two_images import plot_image_comparison
except:
    print("Error importing python file compare_two_images.py (see in utils/ directory), please make sure it is available/linked in this working directory...")
    sys.tracebacklimit=0
    raise Exception("Importing compare_two_images.py failed")


# two images (produced by: gnuplot plot_kernel.gnu)
image_original_reference = "REF_KERNEL/rho_and_kappa_kernels_acoustic_and_elastic_only_no_caption.png"
image_from_new_calculation = "OUTPUT_FILES/rho_and_kappa_kernels_acoustic_and_elastic_only_no_caption.png"

# show a plot
show_plot = "-p" in sys.argv

plot_image_comparison(image_original_reference,image_from_new_calculation,show_plot)
