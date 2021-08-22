#!/bin/bash

file=$1

if [ "$file" == "" ]; then echo "usage: ./plot_kernel.sh kernel-file [e.g. OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat]"; exit 1; fi

if [ ! -f $file ]; then echo "file not found: $file"; exit 1; fi

# kernel visualization
echo "plotting kernel: $file"

# determines image name
# example proc000000_rhop_alpha_beta_kernel.dat  -> rhop
filename=`basename $file`
kernelrho=`echo "$filename" | cut -d '_' -f 2`
# sets image name
if [ "$kernelrho" == "rhop" ]; then
  imagename=image_rhop_alpha_beta_kernels.png
else
  imagename=image_rho_kappa_mu_kernels.png
fi

# creates temporary plot script
sed "s:OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat:$file:g" plot_kernel.gnu > tmp.gnu
sed -i "s:image_rho_kappa_mu_kernels.png:$imagename:g" tmp.gnu

gnuplot tmp.gnu

echo
echo "see image: $imagename"
echo



