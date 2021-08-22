#!/bin/bash -eu
#
#reproduce results from Tromp et al. 2010
#
# Note that one noise sensitivity kernel contains two contributions,
# both the 1st and the 2nd contributions may be obtained through THREE steps.
# Each contribution requires a distinct 'main' receiver, as shown in Tromp et al., 2010, GJI,
# each step requires slightly different Par_file, as documented in the Manual
#
# Please see the paper & Manual for more details on noise simulations
date
echo "running directory: `pwd`"
echo

## negative branch
date
echo "*************************************"
echo "1. contribution..."
echo "*************************************"
echo

./run_kernel_branch.sh 3 1 0
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# store contribution
rm -rf contribution_1_BA_NEG
mv -v OUTPUT_ALL contribution_1_BA_NEG
echo

## positive branch
date
echo "*************************************"
echo "2. contribution..."
echo "*************************************"
echo

./run_kernel_branch.sh 1 3 1
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

rm -rf contribution_2_AB_POS
mv -v OUTPUT_ALL contribution_2_AB_POS
echo

# put all in a new output directory:
rm -rf OUTPUT_FILES
mkdir -p OUTPUT_FILES
mv -v contribution_1_BA_NEG OUTPUT_FILES/
mv -v contribution_2_AB_POS OUTPUT_FILES/

echo
echo `date`
echo "see all output in directory: OUTPUT_FILES/"
echo "done"
echo


