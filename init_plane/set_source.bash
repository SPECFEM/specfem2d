#!/bin/bash
if [ $# != 3 ]; then
  echo 'set_source +/- angle nrec'; exit
fi
sign=$1
angle=$2
nrec=$3

sed -i "s/^angleforce .*$/angleforce                      = $sign$angle/"   SOURCE_Slave
if [ $sign = "+" ]; then
    sed -i "s/^xs .*$/xs                              = -150000/" SOURCE_Slave
    sed -i "s/^xdeb .*$/xdeb                            = 100000./" Par_file_Slave_for Par_file_Slave_kernel
    sed -i "s/^xfin .*$/xfin                            = 500000./" Par_file_Slave_for Par_file_Slave_kernel
    sed -i "s/^jpeg_header.*$/jpeg_header='left-kernel-';/" plot_kernel.m
elif [ $sign = "-" ]; then
   sed -i "s/^xs .*$/xs                              = 750000/" SOURCE_Slave
   sed -i "s/^xdeb .*$/xdeb                            = 500000./" Par_file_Slave_for Par_file_Slave_kernel
   sed -i "s/^xfin .*$/xfin                            = 100000./" Par_file_Slave_for Par_file_Slave_kernel
   sed -i "s/^jpeg_header.*$/jpeg_header='right-kernel-';/" plot_kernel.m
fi

sed -i "s/^nrec .*$/nrec                            = $nrec/" Par_file_Slave_for
sed -i "s/^nrec .*$/nrec                            = $nrec/" Par_file_Slave_kernel

# set select_adj.m and plot_kernel.m
# X0 is hard-wired into selected_adj.m
sed -i "s/^angle_inc=.*$/angle_inc=$angle;/"   select_adj.m
sed -i "s/^nsta=.*$/nsta=$nrec;/" select_adj.m
sed -i "s/^angle_inc=.*$/angle_inc=$angle;/"   plot_kernel.m
