
How to use SPECFEM2D:
---------------------

See file "todo_list_please_dont_remove.txt" for a list of known bugs, problems, or missing options.

To use the code:

- edit the Makefile. There are several options available : -DUSE_MPI compiles with use of an MPI library. -DUSE_METIS enables use of graph partitioner METIS, the same goes for -DUSE_SCOTCH for SCOTCH.

- type "make all"

- edit the input file "DATA/Par_file" which describes the simulation. It contains comments and should be almost self-explanatory, if you need more details we do not have a manual for the 2D version but you can find useful information in the manuals of the 3D versions, since many parameters and the general philosophy is similar. They are available at http://geodynamics.org/wsvn/cig/seismo/3D in subdirectories USER_MANUAL. To create acoustic (fluid) regions, just set the S wave speed to zero and the code will see that these elements are fluid and switch to the right equations there automatically, and automatically match them with the solid regions

- if you are using an external mesher (like GID or CUBIT), you should set "read_external_mesh" to true.
     "mesh_file" is the file describing the mesh : first line is the number of elements, then a list of 4 nodes (quadrilaterals only) forming each elements on each line.
     "nodes_coords_file" is the file containing the coordinates (x and z) of each nodes : number of nodes on the first line, then coordinates x and z on each line.
     "materials_file" is the number of the material for every elements : an integer ranging from 1 to nbmodels on each line.
     "free_surface_file" is the file describing the edges forming the acoustic free surface : number of edges on the first line, then on each line number of the element, number of nodes forming the free surface (1 for a point, 2 for an edge), the nodes forming the free surface for this element. If you do not want free surface, jusr put 0 on the first line.
     "absorbing_surface_file" is the file describing the edges forming the absorbing boundaries : the format is the same as the "free_surface_file".
     "tangential_detection_curve_file" contains points describing the envelope, used for source_normal_to_surface and rec_normal_to_surface. Should be fine grained, and ordained clockwise. Number of points on the first line, then (x,z) coordinates on each line.

- if you have compiled with MPI, you can specify the number of processes, and the partitioning method used to dispatch the elements on the processes. See the manual of METIS and SCOTCH for more informations on the partitioning strategies.

- then type xmeshfem2D to create the mesh (which will be stored in directory OUTPUT_FILES). xmeshfem2D is serial; it will output several files called Databasexxxxx, one for each process.

- then type xspecfem2D to run the main solver (use mpirun or equivalent if you compiled ). This will output the seismograms and snapshots of the wave fronts at different time steps in directory OUTPUT_FILES. To visualize them, type "gs OUTPUT_FILES/vect*.ps" to see the Postscript files (in which the wave field is represented with small arrows, fluid/solid matching interfaces with a thick pink line, and absorbing edges with a thick green line) and "gimp OUTPUT_FILES/image*.gif" to see the color snapshot showing a pixelized image of one of the two components of the wave field (or pressure, depending on what you have selected for the output in DATA/Par_file).

- the DATA/Par_file given with the code works fine, you can use it without any modification to test the code

- the seismograms OUTPUT_FILES/*.sem* are simple ASCII files with two columns: time in the first colum and amplitude in the second, therefore they can be visualized with any tool you like, for instance "gnuplot"

- if you set flag "assign_external_model" to .true. in DATA/Par_file, the velocity and density model that is given at the end of DATA/Par_file is then ignored and overwritten by the external velocity and density model that you define yourself in define_external_model.f90

- when compiling with Intel ifort, use " -assume byterecl " option to create binary PNM images displaying the wave field

- you can convolve them with any source time function in postprocessing later using "convolve_source_timefunction.csh" and "convolve_source_timefunction.f90", see the manual of the 3D code for details on how to do this

- we do not have PML absorbing conditions implemented in the fluid/solid code yet. We use (older and less efficient) paraxial Clayton-Engquist or Sommerfeld equations instead. This is only by lack of time, we have a developer who is currently implementing PML but the code is not fully ready. For now, since the paraxial conditions are less efficient, please use a larger model

- there are a few useful scripts and Fortran routines in directory UTILS

- if you find bugs (or if you have comments or suggestions) please send an email to cig-seismo AT geodynamics.org and the developers will try to fix them and send you an updated version

--------------------------

Regarding the structure of some of the database files:

Question: Can anyone tell me what the columns of the SPECFEM2D-5.2.2 boundary
condition files in SPECFEM2D-5.2.2/DATA/Mesh_canyon are?

SPECFEM2D-5.2.2/DATA/Mesh_canyoncanyon_absorbing_surface_file
SPECFEM2D-5.2.2/DATA/Mesh_canyoncanyon_free_surface_file

Answer: "canyon_absorbing_surface_file" refers to parameters related to the
absorbing conditions:
The first number (180) is the number of absorbing elements (nelemabs in the
code).
Then the columns are:
column 1 = the element number
column 2 = the number of nodes of this element that form the absorbing surface
column 3 =  the first node
column 4 = the second node

"canyon_free_surface_file" refers to the elements of the free surface
(relevant for enforcing free surface condition for acoustic media):
The first number (160) is the number of  elements of the free surface.
Then the columns are (similar to the absorbing case):
column 1 = the element number
column 2 = the number of nodes of this element that form the absorbing surface
column 3 =  the first node
column 4 = the second node

Concerning the free surface description file, nodes/edges pertaining to
elastic elements are discarded when the file is read (if for whatever
reason it was simpler to include all the nodes/edges on one side of a
studied area and that there are among them some elements that are
elastic elements, only the nodes/edges of acoustic elements are kept).

These files are opened and read in meshfem2D.F90 using subroutines
read_abs_surface and read_acoustic_surface, which are in part_unstruct.F90

