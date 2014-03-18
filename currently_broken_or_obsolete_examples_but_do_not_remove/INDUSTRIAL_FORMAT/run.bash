#!/bin/bash -eu

DIR_SEM="../../"
FILE="$DIR_SEM/DATA/Par_file"
sed -e "s#^nt.*#nt    =  10000 #g"         < $FILE > temp; mv temp $FILE
sed -e "s#^assign_external_model.*#assign_external_model    =  .true. #g"         < $FILE > temp; mv temp $FILE
sed -e "s#^READ_EXTERNAL_SEP_FILE.*#READ_EXTERNAL_SEP_FILE   =  .true. #g"        < $FILE > temp; mv temp $FILE

cd $DIR_SEM
./xmeshfem2D
./xspecfem2D
