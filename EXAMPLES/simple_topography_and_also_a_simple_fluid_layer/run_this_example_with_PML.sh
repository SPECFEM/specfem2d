#!/bin/bash

# setup
sed -i "s:^NPROC .*:NPROC   = 4:" DATA/Par_file_with_PML
cp -v DATA/Par_file_with_PML DATA/Par_file

./run_this_example.sh


