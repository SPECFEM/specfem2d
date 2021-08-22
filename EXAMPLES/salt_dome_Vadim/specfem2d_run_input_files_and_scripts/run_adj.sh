#!/bin/bash
#OAR -l nodes=4,walltime=20:00:00
#OAR -n test_48_cores

module load intel/13.0
module load openmpi/intel/1.6.3
module list

ulimit -s unlimited
echo ${OAR_NODEFILE}

SEM=/home/vmonteil/progs/specfem2d/bin/

MPIRUN=mpirun
OPTION_MPI=" -np 48  -machinefile "${OAR_NODEFILE}" -bysocket -bind-to-core"

declare -i i n

i=1
n=19

while [ ${i} -le ${n} ]; do

mkdir OUTPUT_FILES
mkdir OUTPUT_FILES${i}
mkdir SEM
#mv OUTPUT_FILES${i} OUTPUT_FILES # forward already computed

cp DATA/Par_file_fwd DATA/Par_file
cp ../2D_salt_simu/DATA/SOURCE${i} DATA/SOURCE

# direct run with save snapshot
$MPIRUN $OPTION_MPI  $SEM/xmeshfem2D > output_mesher1_${i}.txt
$MPIRUN $OPTION_MPI $SEM/xspecfem2D > output_solver1_${i}.txt


cp DATA/Par_file_adj DATA/Par_file

./compute_adjoint_source.x  ../2D_salt_simu/OUTPUT_FILES${i}/  ./OUTPUT_FILES/ ./SEM/

# adjoint run and retored direct run from snapshots
$MPIRUN $OPTION_MPI $SEM/xmeshfem2D > output_mesher2_${i}.txt
$MPIRUN $OPTION_MPI $SEM/xspecfem2D > output_solver2_${i}.txt

mv OUTPUT_FILES/proc*ker*  OUTPUT_FILES${i}/.
mv SEM SEM${i}

i="$i+1"

done;
