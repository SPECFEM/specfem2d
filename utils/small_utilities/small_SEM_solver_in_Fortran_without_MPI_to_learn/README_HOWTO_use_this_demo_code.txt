
Demo code written by Dimitri Komatitsch, CNRS, Marseille, France, September 2015
--------------------------------------------------------------------------------

Useful to learn how the spectral-element method works, and how to write or modify a code to implement it.
Also useful to test new ideas by modifying this simple code to run some tests.

Here is how to compile and use it:

1/ edit make_Fortran_2D_plane_strain.csh if needed to change the compiler or compiler options

2/ ./make_Fortran_2D_plane_strain.csh

3/ ./xspecfem2D

and that's it. Then type:

  gs vect*.ps

to see the results displayed in PostScript format, and type:

  gnuplot plotall.gnu

to see the seismogram.

