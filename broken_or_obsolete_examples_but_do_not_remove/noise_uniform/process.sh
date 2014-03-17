#!/bin/bash

#reproduce results from Tromp et al. 2010

./use_negative_branch 3 1
mv OUTPUT_ALL CBA_NEG

./use_positive_branch 1 3
mv OUTPUT_ALL CAB_POS

