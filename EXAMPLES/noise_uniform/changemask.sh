#!/bin/bash

mask=$1
if [ "$mask" == "" ]; then echo "./changemask [uniform/nonuniform]"; exit 1; fi

echo
echo "mask: $mask"
echo
if [ "$mask" == "uniform" ]; then
echo "setting uniform mask ..."
echo "0" > NOISE_TOMOGRAPHY/use_non_uniform_noise_distribution

elif [ "$mask" == "nonuniform" ]; then
echo "setting nonuniform mask ..."
echo "1" > NOISE_TOMOGRAPHY/use_non_uniform_noise_distribution

else
echo "mask: $mask not recognized, please use uniform or nonuniform"
exit 1
fi

echo
echo "done, see file: NOISE_TOMOGRAPHY/use_non_uniform_noise_distribution"
echo
