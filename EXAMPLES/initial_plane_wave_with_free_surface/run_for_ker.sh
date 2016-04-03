#!/bin/bash
#
# script runs mesher and solver (in serial)
# using this example setup
#

if [ $# != 2 ]; then
  echo 'process.sh forward kernel'; exit
fi
forward=$1
kernel=$2
echo ' '
date
if [ $forward -eq 1 ] ; then
currentdir=`pwd`
mkdir -p OUTPUT_FILES
mkdir -p DATA
rm -rf OUTPUT_FILES/* OUTPUT_FILES_FORWARD OUTPUT_FILES_KERNEL

cd ../../
make > tmp.log
cd $currentdir
rm -f xmeshfem2D xspecfem2D
ln -s ../../bin/xmeshfem2D
ln -s ../../bin/xspecfem2D

# sets up local DATA/ directory
cd DATA/
rm -f Par_file SOURCE
ln -s ../Par_file_Slave_for Par_file
ln -s ../SOURCE_Slave SOURCE
cd ../

echo 'Forward simulation ...'
echo "  running mesher..."
./xmeshfem2D > OUTPUT_FILES/output_mesher.txt
echo "  running solver..."
./xspecfem2D > OUTPUT_FILES/output_solver.txt
cp DATA/SOURCE_xz.dat DATA/STATIONS_* DATA/Par_file DATA/SOURCE OUTPUT_FILES/
grep Receiver OUTPUT_FILES/output_mesher.txt | awk '{print $4,$5}' > OUTPUT_FILES/receiver.dat
mv OUTPUT_FILES OUTPUT_FILES_FORWARD
echo "see results in directory: OUTPUT_FILES_FORWARD/"
fi


if [ $kernel -eq 1 ]; then
# run select_adj in matlab
matlab.sh "select_adj"

echo 'Kernel simulation ...'
mkdir -p OUTPUT_FILES/
cd SEM; rm -f *; cd ..
cd OUTPUT_FILES_FORWARD; cp lastframe* absorb* ../OUTPUT_FILES; cp *.adj ../SEM; cd ..
cd SEM; rename 's/semd.adj/adj/' *.adj; cd ..
cd DATA; rm -f Par_file; ln -s ../Par_file_Slave_kernel Par_file; cd ..
echo "  running mesher..."
./xmeshfem2D > OUTPUT_FILES/output_mesher.txt
echo "  running solver..."
./xspecfem2D > OUTPUT_FILES/output_solver.txt
cp DATA/Par_file OUTPUT_FILES/
mv OUTPUT_FILES OUTPUT_FILES_KERNEL
echo "see results in directory: OUTPUT_FILES_KERNEL/"

matlab.sh "plot_kernel"
fi


echo "done"
date
