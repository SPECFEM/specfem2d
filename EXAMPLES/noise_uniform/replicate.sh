#!/bin/bash

if [ $# -ne 2 ]; then echo "Usage: ./replicate.sh  par_file  dir"; exit 1; fi

if [ ! -e $1 ];  then echo "File not found: $1"; exit 1; fi
if [ ! -d $2 ];  then echo "Directory not found: $2"; exit 1; fi

FILE=$1
DIR=$2

FILE_BAK=$DIR/replicate.bak
cp $FILE $FILE_BAK
cd $DIR

Line1=`grep ^SIMULATION_TYPE $FILE_BAK`
Line2=`grep ^NOISE_TOMOGRAPHY $FILE_BAK`
Line3=`grep ^SAVE_FORWARD $FILE_BAK`

if [ `echo $Line1 | wc -l` -ne 1 ]; then echo "Error reading SIMULATION_TYPE."; exit 1; fi
if [ `echo $Line2 | wc -l` -ne 1 ]; then echo "Error reading NOISE_TOMOGRAPHY."; exit 1; fi
if [ `echo $Line3 | wc -l` -ne 1 ]; then echo "Error reading SAVE_FORWARD."; exit 1; fi


#PAR_FILE_NOISE_1
cp $FILE_BAK Par_file_noise_1
sed -i -e "s:${Line1}:SIMULATION_TYPE                 = 1:"  Par_file_noise_1
sed -i -e "s:${Line2}:NOISE_TOMOGRAPHY                = 1:"  Par_file_noise_1
sed -i -e "s:${Line3}:SAVE_FORWARD                    = .false.:"  Par_file_noise_1

#PAR_FILE_NOISE_2
cp $FILE_BAK Par_file_noise_2
sed -i -e "s:${Line1}:SIMULATION_TYPE                 = 1:"  Par_file_noise_2
sed -i -e "s:${Line2}:NOISE_TOMOGRAPHY                = 2:"  Par_file_noise_2
sed -i -e "s:${Line3}:SAVE_FORWARD                    = .true.:"  Par_file_noise_2

#PAR_FILE_NOISE_3
cp $FILE_BAK Par_file_noise_3
sed -i -e "s:${Line1}:SIMULATION_TYPE                 = 2:"  Par_file_noise_3
sed -i -e "s:${Line2}:NOISE_TOMOGRAPHY                = 3:"  Par_file_noise_3
sed -i -e "s:${Line3}:SAVE_FORWARD                    = .false.:"  Par_file_noise_3

