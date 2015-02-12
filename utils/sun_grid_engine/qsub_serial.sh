#!/bin/sh

# shell to be used
#$ -S /bin/sh

# name of the job
#$ -N SPECFEM2D

# this is to make the current directory the working directory
#$ -cwd


# SGE script to use on the UPPA cluster, usually called by submit_sun_grid_engine_UPPA.sh
# It dispatches the Database to the corresponding node in its local storage devices (assuming /scratch/$user already exists).
# The host of the scp from iplmas014 should be changed according to the user's lab (can be ipigps014 for MIGP).
# Many unnecessary copies should be disposed of as soon as FS sync issues have been addressed.




CURRENT_DIR=$PWD
echo $CURRENT_DIR

mkdir $JOB_NAME$JOB_ID
mkdir $JOB_NAME$JOB_ID/OUTPUT_FILES
mkdir $JOB_NAME$JOB_ID/DATA

cd $JOB_NAME$JOB_ID

scp iplmas014.univ-pau.fr:$CURRENT_DIR/OUTPUT_FILES/Database* ./OUTPUT_FILES/
scp iplmas014.univ-pau.fr:$CURRENT_DIR/DATA/STATIONS ./DATA/
scp iplmas014.univ-pau.fr:$CURRENT_DIR/xspecfem2D ./

cp -r ../$JOB_NAME$JOB_ID /scratch/$USER/

cd $CURRENT_DIR

rm -r $JOB_NAME$JOB_ID/

rm ./OUTPUT_FILES/Database*

cd /scratch/$USER/$JOB_NAME$JOB_ID

./xspecfem2D

cp /scratch/$USER/$JOB_NAME$JOB_ID/OUTPUT_FILES/* $CURRENT_DIR/OUTPUT_FILES/

rm -r /scratch/$USER/$JOB_NAME$JOB_ID
