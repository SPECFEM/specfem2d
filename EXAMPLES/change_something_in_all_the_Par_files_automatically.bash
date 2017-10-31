#!/bin/bash

#
# If you want to change the name of a parameter or anything else in all the Par_files in an automatic way
# you can also use a "find" command combined with a "sed" command, below is an example replacing
# ATTENUATION_VISCOELASTIC_SOLID by ATTENUATION_VISCOELASTIC in all Par_files
#
# However this will *not* make the change in ../DATA/Par_file, which you will need to do separately.
#

#  [ ]*  below is very useful, it matches any number of white spaces, thus makes the search insensitive to where the = sign is located in the initial line

find . -type f -iname "Par_file*" -exec sed -i 's/UNDO_ATTENUATION[ ]*=/UNDO_ATTENUATION_AND_OR_PML     =/g' {} \;
#find . -type f -iname "Par_file*" -exec sed -i 's/ATTENUATION_FLUID/ATTENUATION_VISCOACOUSTIC/g' {} \;

