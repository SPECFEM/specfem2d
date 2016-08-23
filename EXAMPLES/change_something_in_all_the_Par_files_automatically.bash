#!/bin/bash

#
# If you want to change the name of a parameter or anything else in all the Par_files in an automatic way
# you can also use a "find" command combined with a "sed" command, below is an example.
#
# However this will *not* make the change in ../DATA/Par_file, which you will need to do separately.
#

find . -type f -name "*Par_file*" -exec sed -i 's/ATTENUATION_VISCOELASTIC_SOLID/ATTENUATION_VISCOELASTIC/g' {} \;

