#!/usr/bin/perl

#
#  Script to clean a given file
#
#  Author : Dimitri Komatitsch, EPS - Harvard University, May 1998
#

      $name  = $ARGV[0];
      $nametoprint  = $ARGV[1];

# change tabs to white spaces
            system("expand -6 < $name > _____testzzzXXXyyy_____");
            $f90name = $name;
            print STDOUT "Cleaning file $nametoprint ...\n";

            open(FILEF77,"<_____testzzzXXXyyy_____");
            open(FILEF90,">$f90name");

# read the input file
      while($line = <FILEF77>) {

# suppress trailing white spaces and carriage return
      $line =~ s/\s*$//;

# write the output line
      print FILEF90 "$line\n";

  }

            close(FILEF77);
            close(FILEF90);

            system("rm -f _____testzzzXXXyyy_____");

