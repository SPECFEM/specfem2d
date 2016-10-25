#!/bin/bash
#
# runs all three noise tomography simulations
#

echo "Cleaning up."

rm -rf SEM OUTPUT_FILES OUTPUT_ALL SNAPSHOTS &> /dev/null
rm xmeshfem2D xspecfem2D &> /dev/null
rm adj_run &> /dev/null
rm log.* &> /dev/null


