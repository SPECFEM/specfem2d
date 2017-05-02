#!/bin/bash

# extract left edge information (flag value == 4) from the absorbing surface file

grep ' 4$' modelY1_absorbing_surface_file | gawk '{print $3  "\n" $4}' | sort -n | uniq > left_edge

