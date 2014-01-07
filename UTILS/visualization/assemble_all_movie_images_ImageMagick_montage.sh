#!/bin/sh

rm movie_image*.jpg

# montage -gravity Center -page A4 -border 2 -geometry +0+0 -tile 4x4 -resize 33% -strip -depth 8 image*.jpg movie_image.jpg
montage -border 2 -geometry +0+0 -resize 33% -strip -depth 8 image*.jpg movie_image.jpg

