#!/bin/bash

file=$1

if [ "$file" == "" ]; then echo "usage: ./plot_kernel.sh kernel-file [e.g. OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat]"; exit 1; fi

if [ ! -f $file ]; then echo "file not found: $file"; exit 1; fi

# kernel visualization
echo "plotting kernel: $file"

sed "s:OUTPUT_FILES/proc000000_rho_kappa_mu_kernel.dat:$file:g" plot_kernel.gnu > tmp.gnu

gnuplot tmp.gnu

echo
echo "see image: image_rho_kappa_mu_kernels.png"
echo



