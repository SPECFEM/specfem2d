#!/bin/sh

USAGE="USAGE:  submit  par_file_directory"
if [ $# -eq 0 ]; then PAR_DIR=uniform; fi
if [ $# -eq 1 ]; then PAR_DIR=$1; fi
if [ $# -ne 0 ] && [ $# -ne 1 ]; then echo "$USAGE"; exit 1; fi


PAR_DIR_FULL=$PWD/$PAR_DIR
RUN_DIR=../..
cd $RUN_DIR


# prepare directories
rm -rf   SEM NOISE_TOMOGRAPHY OUTPUT_FILES OUTPUT_ALL
mkdir -p SEM NOISE_TOMOGRAPHY OUTPUT_FILES OUTPUT_ALL


# prepare files
cp $PAR_DIR_FULL/SOURCE_noise           DATA/SOURCE
cp $PAR_DIR_FULL/STATIONS_target_noise  DATA/STATIONS_target
#cp $PAR_DIR_FULL/S_squared              NOISE_TOMOGRAPHY
echo 1 >                                NOISE_TOMOGRAPHY/irec_master


#simulation 1
cp $PAR_DIR_FULL/Par_file_noise_1  DATA/Par_file
make; bin/xmeshfem2D; bin/xspecfem2D
mkdir OUTPUT_ALL/step_1
mv OUTPUT_FILES/image*  OUTPUT_ALL/step_1
mv OUTPUT_FILES/*.semd  OUTPUT_ALL/step_1
mv DATA/Par_file  OUTPUT_ALL/step_1


#simulation 2
cp $PAR_DIR_FULL/Par_file_noise_2  DATA/Par_file
bin/xmeshfem2D; bin/xspecfem2D
mkdir OUTPUT_ALL/step_2

ADJ_CODE=$PAR_DIR_FULL/adj_cc.f90
gfortran $ADJ_CODE -o adj_cc
cp OUTPUT_FILES/*.semd  SEM
./adj_cc SEM/S0003.AA.BXY.semd

cd SEM
rename '.semd' '' *.adj
awk '{printf(" %20.10f %20.10f\n",$1,0.)}' < S0003.AA.BXY.adj > ZEROS
cp ZEROS S0001.AA.BXX.adj; cp ZEROS S0001.AA.BXY.adj; cp ZEROS S0001.AA.BXZ.adj
cp ZEROS S0002.AA.BXX.adj; cp ZEROS S0002.AA.BXY.adj; cp ZEROS S0002.AA.BXZ.adj
cp ZEROS S0003.AA.BXX.adj; cp ZEROS S0003.AA.BXZ.adj

cd ..

mv OUTPUT_FILES/image*  OUTPUT_ALL/step_2
mv OUTPUT_FILES/*.semd  OUTPUT_ALL/step_2
mv DATA/Par_file  OUTPUT_ALL/step_2


#simulation 3
cp $PAR_DIR_FULL/Par_file_noise_3  DATA/Par_file
bin/xmeshfem2D; bin/xspecfem2D
mkdir OUTPUT_ALL/step_3
mv OUTPUT_FILES/image*  OUTPUT_ALL/step_3
mv OUTPUT_FILES/*.semd  OUTPUT_ALL/step_3
mv OUTPUT_FILES/proc*   OUTPUT_ALL/step_3
mv DATA/Par_file  OUTPUT_ALL/step_3

