#!/bin/sh

RUN_DIR=$PWD

# prepare directories
rm -rf   DATA SEM NOISE_TOMOGRAPHY OUTPUT_FILES OUTPUT_ALL
mkdir -p DATA SEM NOISE_TOMOGRAPHY OUTPUT_FILES OUTPUT_ALL


# prepare files
cp SOURCE_noise           DATA/SOURCE
cp STATIONS_target_noise  DATA/STATIONS_target
cp S_squared              NOISE_TOMOGRAPHY
cp uniform.dat            DATA/
echo 1 >                  NOISE_TOMOGRAPHY/irec_master
cd ../..
make
cd $RUN_DIR
ln -s ../../bin/xmeshfem2D .
ln -s ../../bin/xspecfem2D .

#simulation 1
cp Par_file_noise_1  DATA/Par_file
./xmeshfem2D; ./xspecfem2D
mkdir OUTPUT_ALL/step_1
mv OUTPUT_FILES/image*  OUTPUT_ALL/step_1
mv OUTPUT_FILES/*.semd  OUTPUT_ALL/step_1
mv DATA/Par_file  OUTPUT_ALL/step_1


#simulation 2
cp Par_file_noise_2  DATA/Par_file
./xmeshfem2D; ./xspecfem2D
mkdir OUTPUT_ALL/step_2

ADJ_CODE=adj_cc.f90
gfortran $ADJ_CODE -o adj_run
cp OUTPUT_FILES/*.semd  SEM
./adj_run SEM/S0003.AA.BXY.semd

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
cp Par_file_noise_3  DATA/Par_file
./xmeshfem2D; ./xspecfem2D
mkdir OUTPUT_ALL/step_3
mv OUTPUT_FILES/image*  OUTPUT_ALL/step_3
mv SEM/*Y.adj           OUTPUT_ALL/step_3
mv OUTPUT_FILES/proc*   OUTPUT_ALL/step_3
mv DATA/Par_file        OUTPUT_ALL/step_3

