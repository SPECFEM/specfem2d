#!/bin/bash

codefile=../../src/specfem2D/noise_tomography.f90

cd masks
select maskfile in `find -type f`
do
case $maskfile in
    *) echo $maskfile; break;;
esac
done
cd ..


line1=`grep -n '^[ \t]*subroutine create_mask_noise' $codefile | cut -d':' -f1`
line2=`grep -n '^[ \t]*end subroutine create_mask_noise' $codefile | cut -d':' -f1`

if [ `echo $line1 | wc -l` -ne 1 ]; then echo "Error reading noise_tomography.f90"; exit 1; fi
if [ `echo $line2 | wc -l` -ne 1 ]; then echo "Error reading noise_tomography.f90"; exit 2; fi
if [ ! -n "$line1" ]; then echo "Error reading noise_tomography.f90"; exit 3; fi
if [ ! -n "$line2" ]; then echo "Error reading noise_tomography.f90"; exit 4; fi

awk -vline1=$line1 '{if (NR < line1) print $0}' < $codefile > temp1
awk -vline2=$line2 '{if (NR > line2) print $0}' < $codefile > temp2

cat temp1 masks/$maskfile temp2 > $codefile
rm temp?
