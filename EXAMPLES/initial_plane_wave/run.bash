#!/bin/bash

for dir in + -; do
 for ang in 23 20 18; do
   for nrec in 10 20 30; do
     set_source.bash $dir $ang $nrec
     run_for_ker.sh 1 1
     mkdir -p ang-$ang; mkdir -p ang-$ang/nrec-$nrec-dir-$dir
     mv OUTPUT_FILES_FORWARD OUTPUT_FILES_KERNEL ang-$ang/nrec-$nrec-dir-$dir
   done
  done
done
