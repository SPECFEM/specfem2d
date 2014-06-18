#!/bin/bash

./make_specific_mesher_for_2D_Earth.csh

mkdir -p OUTPUT_FILES
../../bin/xmeshfem2D
../../bin/xspecfem2D

