#!/usr/bin/perl

#
# read and clean all Fortran files in the current directory and subdirectories
#

## DK DK April 2018: convert all "stop" statements in the code to calls to a subroutine called "stop_the_code", which itself calls exit_MPI()
## DK DK April 2018: this avoids wasting CPU hours on clusters when a thread exits with a stop statement, which may not kill the whole run properly, and the run may keep hanging.

## DK DK April 2018: two known minor issues for this script:
## - it replaces the word " stop " even if it appears in a sentence, for instance in a print statement; to avoid problems with that, in the SPECFEM codes we use "stopping" instead of "stop" in print statements
## - it has a minor issue if the stop is in a statement that extends over more than a single line, i.e. a line split using the & continuation symbol of Fortran;
##   if so, the script puts the closing parenthesis at the end of the first line instead of at the end of the group of lines (because the script analyzes the code
##   line by line for simplicity, instead of globally). Fixing this would be a lot of work, and this does not happen very often, thus instead I use this script,
##   then compile the Fortran codes and manually fix the very few lines for which the parenthesis is misplaced (the compiler will stop there, thus easy to locate them).
## Because of these two known issues, I do not put this script in the automatic cleaning script of Buildbot (specfem3d/utils/clean_listings_specfem.pl)

# DK DK only do this in the "src" directory, otherwise independent programs in other directories such as "utils" will not have access to the "stop_the_code()" subroutine
# DK DK also prune the "src/auxiliaries" directory, which contains independent programs rather than subroutines of the main code.

# when using this "find" command from Perl we need to use \\ instead of \ below otherwise Perl tries to interpret it
      @objects = `find 'src' -name '.git' -prune -o -name 'm4' -prune -o -path './src/auxiliaries' -prune -o -path './auxiliaries' -prune -o -path 'src/auxiliaries' -prune -o -path 'auxiliaries' -prune -o -type f -regextype posix-extended -regex '.*\\.(fh|f90|F90|h\\.in|fh\\.in)' -print`;

      foreach $name (@objects) {
            chop $name;
# change tabs to white spaces
            system("expand -2 < $name > _____temp08_____");
            $f90name = $name;
            print STDOUT "Cleaning $f90name ...\n";

            open(FILE_INPUT,"<_____temp08_____");
            open(FILEF90,">$f90name");

# open the input f90 file
      while($line = <FILE_INPUT>) {

# suppress trailing white spaces and carriage return
      $line =~ s/\s*$//;

      $linewithnospaceatall = $line;
      $linewithnospaceatall =~ s# ##ogi;
      $first_letter = substr(($linewithnospaceatall),0,1);
# do not make replacements in comments
      if($first_letter ne '!') {

# add a space if "stop" starts at the first column, to avoid having to test for that particular case in what follows
        $line =~ s#^stop # stop #ogi;

# add a dummy error message if "stop" ends at the last column, to avoid having a stop statement with an empty message
        $line =~ s# stop$# stop 'error: stopping the code'#ogi;

# test if the line contains a "stop" statement
        if (index($line, " stop ") != -1) {

# get the first part of the line, before the "stop" statement
        $first_part_of_line = substr($line, 0, index($line, ' stop '));

# get the last part of the line, after the "stop" statement
# the + 6 corresponds to the length of " stop " with two spaces
        $last_part_of_line = substr($line, index($line, ' stop ') + 6, length($line));

        $line = $first_part_of_line . " call stop_the_code(" . $last_part_of_line . ")";

# if the new line created is longer than the standard
        if (length($line) >= 132) {
          # split the line using a & line continuation symbol
            $line = $first_part_of_line . " call stop_the_code( &\n" . $last_part_of_line . ")";
          }

# suppress trailing white spaces, just in case we have added any in the above processing
          $line =~ s/\s*$//;

        }
      }

      print FILEF90 "$line\n";

      }

      close(FILE_INPUT);
      close(FILEF90);

      }

      system("rm -f _____temp08_____");

