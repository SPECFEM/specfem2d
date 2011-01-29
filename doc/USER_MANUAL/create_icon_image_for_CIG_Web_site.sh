#!/bin/sh

# Dimitri Komatitsch, November 2010: use "convert" from the Imagemagick package

convert -quality 100 -resize 38% -colorspace RGB figures/specfem_2d-cover.pdf jpg:specfem_2d-cover_small.jpg

