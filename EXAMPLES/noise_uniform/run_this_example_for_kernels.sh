#!/bin/bash

#reproduce results from Tromp et al. 2010

# negative branch
./run_kernel_branch.sh 3 1 0
mv OUTPUT_ALL OUTPUT_CBA_NEG

# positive branch
./run_kernel_branch.sh 1 3 1
mv OUTPUT_ALL OUTPUT_CAB_POS

