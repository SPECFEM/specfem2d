#!/bin/bash
###################################################

# example directory name
NAME="simple_topography_and_also_a_simple_fluid_layer"

# relative location of SPECFEM2D EXAMPLES/ directory for tests directories (e.g. SPECFEM2D/tests/compilations)
EXAMPLES="../../EXAMPLES/"

###################################################

# bash function for checking seismogram output with reference solutions
my_test(){
  echo "testing seismograms:";
  ../../utils/compare_seismogram_correlations.py OUTPUT_FILES/ REF_SEIS/;
  ../../utils/compare_seismogram_correlations.py OUTPUT_FILES/ REF_SEIS/ | grep min/max | cut -d \| -f 3 | awk '{print "correlation:",$1; if ($1 < 0.9 ){print $1,"failed"; exit 1;}else{ print $1,"good"; exit 0;}}';
}

testdir=`pwd`

# title
echo >> $testdir/results.log
echo "$NAME in: $testdir" >> $testdir/results.log
echo >> $testdir/results.log

# checks if compilation done
if [ ! -e ./bin/xspecfem2D ]; then
  echo "no binaries found, compilation must be done before, please check..." >> $testdir/results.log
  exit 1
fi

# checks if example directory exists
if [ ! -e $EXAMPLES/$NAME/run_this_example.sh ]; then
  echo "run-script in directory $EXAMPLES/$NAME not found, please check..." >> $testdir/results.log
  exit 1
fi

#cleanup output
rm -rf ./OUTPUT_FILES ./DATA ./run_this_example.sh ./REF_SEIS
mkdir -p OUTPUT_FILES

# setup
cp -rp $EXAMPLES/$NAME/DATA .
cp -p $EXAMPLES/$NAME/run_this_example.sh .
ln -s $EXAMPLES/$NAME/REF_SEIS

sed -i "s:^NPROC .*:NPROC    = 4:" DATA/Par_file

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "setup failed, please check..." >> $testdir/results.log
  exit 1
fi

./run_this_example.sh >> $testdir/results.log 2>&1

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "run failed, please check..." >> $testdir/results.log
  exit 1
fi

# test seismograms
my_test >> $testdir/results.log

# checks exit code
if [[ $? -ne 0 ]]; then
  echo >> $testdir/results.log
  echo "seismograms failed, please check..." >> $testdir/results.log
  exit 1
fi

# cleanup
rm -rf ./OUTPUT_FILES ./DATA ./run_this_example.sh ./REF_SEIS

echo "successful compilation" >> $testdir/results.log

