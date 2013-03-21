#!/bin/bash

nohup mpirun -np 20 -machinefile /home/cluster_maintenance/mymachines --byslot ./bin/xspecfem2D > output_night.txt 2> output_night.err &
