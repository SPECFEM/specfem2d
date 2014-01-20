#!/bin/sh

for file in imag*.gif; do
  jfile=`echo $file | sed 's/gif/jpg/'`
  convert -quality 100 $file $jfile
  echo convert -quality 100 $file $jfile
done
mencoder "mf://imag*.jpg" -mf fps=1 -o output.avi -ovc lavc -lavcopts vcodec=mpeg4
