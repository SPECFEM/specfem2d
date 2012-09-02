#!/bin/csh

montage -gravity Center -page A4 -geometry +0+0 -tile 2x2 -strip -depth 8 image*.jpg final_image.jpg

