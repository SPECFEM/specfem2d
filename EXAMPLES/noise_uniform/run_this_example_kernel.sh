#!/bin/bash -eu

#reproduce results from Tromp et al. 2010

# negative branch
./run_kernel_branch.sh 3 1 0

rm -rf OUTPUT_Contribution_BA_NEG
mv OUTPUT_ALL OUTPUT_Contribution_BA_NEG

# positive branch
./run_kernel_branch.sh 1 3 1

rm -rf OUTPUT_Contribution_AB_POS
mv OUTPUT_ALL OUTPUT_Contribution_AB_POS

