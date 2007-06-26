#!/bin/sh

# shell to be used
#$ -S /bin/sh

# name of the job
#$ -N SPECFEM2D

# this is to make the current directory the working directory
#$ -cwd


# SGE script to use on the UPPA cluster, usually called by submit_sun_grid_engine_UPPA.sh
# It dispatches the Databases to the corresponding nodes in their local storage devices (assuming /scratch/$user already exists).
# The current limit on the number of processes is 99 for this script.
# The host of the scp from iplmas014 should be changed according to the user's lab (can be ipigps014 for MIGP).
# Many unnecessary copies should be disposed of as soon as FS sync issues have been addressed.
# The call to clean_scratch_UPPA.sh should be replaced when an alternate way to do this becomes available, like "rsh" 
# or "shmux" for instance (according to OpenMPI manual it can be used, but beware of other MPI libraries, and it is in fact ugly).





PeHostfile2MachineFile()
{
    cat $1
    
    j=0
    
    host0=`head -n 1 $1 | cut -f1 -d" "|cut -f1 -d"."`

    cat $1 | while read line; do
      # echo $line
      host=`echo $line|cut -f1 -d" "|cut -f1 -d"."`
      nslots=`echo $line|cut -f2 -d" "`
      i=1
      
       scp -r ../$3$4 $host.univ-pau.fr:/scratch/$2

      while [ $i -le $nslots ]; do
         echo $host
         i=`expr $i + 1`

	 if [ $j -lt 10 ] 
	     then
	     scp ./OUTPUT_FILES/Database0000$j $host.univ-pau.fr:/scratch/$2/$3$4/OUTPUT_FILES/
	 else
	     if [ $j -lt 100 ]
		 then
		 scp ./OUTPUT_FILES/Database000$j $host.univ-pau.fr:/scratch/$2$3$4/OUTPUT_FILES/
	     fi  
	 fi

	 j=`expr $j + 1`	 
	 

      done
    done
}


CURRENT_DIR=$PWD
echo $CURRENT_DIR

mkdir $JOB_NAME$JOB_ID
mkdir $JOB_NAME$JOB_ID/OUTPUT_FILES
mkdir $JOB_NAME$JOB_ID/DATA

cd $JOB_NAME$JOB_ID

scp iplmas014.univ-pau.fr:$CURRENT_DIR/OUTPUT_FILES/Database* ./OUTPUT_FILES/
scp iplmas014.univ-pau.fr:$CURRENT_DIR/DATA/STATIONS ./DATA/
scp iplmas014.univ-pau.fr:$CURRENT_DIR/xspecfem2D ./
scp iplmas014.univ-pau.fr:$CURRENT_DIR/clean_scratch.sh ./

PeHostfile2MachineFile $PE_HOSTFILE $USER $JOB_NAME $JOB_ID

cd $CURRENT_DIR

rm -r $JOB_NAME$JOB_ID/

rm ./OUTPUT_FILES/Database*

cd /scratch/$USER/$JOB_NAME$JOB_ID


LD_LIBRARY_PATH=/opt/openmpi-1.2.1/gfortran64/lib  /opt/openmpi-1.2.1/gfortran64/bin/mpirun -np 8 ./xspecfem2D

scp $host0.univ-pau.fr:/scratch/$USER/$JOB_NAME$JOB_ID/OUTPUT_FILES/* $CURRENT_DIR/OUTPUT_FILES/

LD_LIBRARY_PATH=/opt/openmpi-1.2.1/gfortran64/lib  /opt/openmpi-1.2.1/gfortran64/bin/mpirun -np 8 ./clean_scratch_UPPA.sh $USER $JOB_NAME $JOB_ID

