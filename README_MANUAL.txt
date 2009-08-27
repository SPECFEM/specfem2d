
SPECFEM2D 6.0: SPECFEM2D facilitates 2D simulations of
        acoustic, (an)elastic, and poroelastic seismic wave propagation.
        With version 6.0, the 2D spectral-element solver accommodates
        regular and unstructured meshes, generated for example by Cubit
        (cubit.sandia.gov). The solver has adjoint capabilities and can
        calculate finite-frequency sensitivity kernels for acoustic,
        (an)elastic, and poroelastic media. Finally, the solver can run
        both in serial and in parallel. See SPECFEM2D
        <http://www.geodynamics.org/cig/software/packages/seismo/specfem2d>
        for the source code.

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

- you can find a Fortran code to compute the analytical solution for simple media that we use as a reference in benchmarks in many of our articles at http://www.spice-rtn.org/library/software/EX2DDIR . That code is described in:

@INCOLLECTION{BeIfNiSk94,
  author = {P. Berg and F. If and P. Nielsen and O. Skovegaard},
  title = {Analytic reference solutions},
  booktitle = {Modeling the Earth for oil exploration, Final report of the CEC's GEOSCIENCE~I Program 1990-1993},
  publisher = {Pergamon Press, Oxford, United Kingdom},
  year = {1994},
  editor = {K. Helbig},
  pages = {421-427}}

How to use poroelasticity:
--------------------------

---------------------------------
      NEW INPUTS IN Par_file
---------------------------------
In section "# geometry of model and mesh description":
TURN_VISCATTENUATION_ON, Q0, and FREQ0 deal with viscous damping in a poroelastic medium.
Q0 is the quality factor set at the central frequency FREQ0. For more details
see Morency & Tromp, GJI 2008.

In section "# time step parameters":
ISOLVER defines the type of simulations
(1) forward simulation
(2) adjoint method and kernels calculation

In section "# source parameters":
The code now support multi sources.
NSOURCE is the number of source.
Parameters of the sources are displayed in the file CMTSOLUTION, which must be
in the directory DATA.

In section "# receiver line parameters for seismograms":
SAVE_FORWARD determines if the last frame of a forward simulation is saved (.true.) or not (.false)

In section "# define models....":
Three types of models:
I: (model_number,1,rho,Vp,Vs,0,0,Qp,Qs,0,0,0,0,0,0), for isotropic elastic/acoustic
material
or II: (model_number,2,rho,c11,c13,c33,c44,Qp,Qs,0,0,0,0,0,0), for anisotropic material
or III: (model_number,3,rhos,rhof,phi,c,kxx,kxz,kzz,Ks,Kf,Kfr,etaf,mufr,Qs),
for isotropic poroelastic material

rho_s = solid density
rho_f = fluid density
phi = porosity
tort = tortuosity
permxx = xx component of permeability tensor
permxz = xz,zx components of permeability tensor
permzz = zz component of permeability tensor
kappa_s = solid bulk modulus
kappa_f= fluid bulk modulus
kappa_fr= frame bulk modulus
eta_f = fluid viscosity
mu_fr = frame shear modulus
Qs = shear quality factor

Note: for the poroelastic case, mu_s is irrelevant.
For details on the poroelastic theory see Morency and Tromp, GJI 2008.

--------------------------------------------------
     HOW TO OBTAIN FINITE SENSITIVITY KERNELS
--------------------------------------------------

First: run a forward simulation
=> isolver = 1
=> save_forward = .true.
=> seismotype = 1 (we need to save the displacement fields to later on derive the
adjoint source. Note: if the user forgets it, the program corrects it when reading the proper
isolver/save_forward combination and a warning message appears in the ouput
file)

Important output files (for example, for the elastic case)
absorb_elastic_bottom*****.bin
absorb_elastic_left*****.bin
absorb_elastic_right*****.bin
absorb_elastic_top*****.bin
lastframe_elastic*****.bin
S****.AA.BHX.semd
S****.AA.BHZ.semd

Second: define the adjoint source
Use adj_seismogram.f90 which is in UTILS/adjoint.
Edit to update NSTEP, nrec, t0, deltat, and the position of the cut to pic
any given phase if needed (tstart,tend), add the right number of stations, and
put one component of the source to zero if needed.
The ouput files of adj_seismogram.f90 are S****.AA.BHX.adj and S****.AA.BHZ.adj. They need to be
kept in the OUTPUT_FILES directory together with the absorb_elastic_****.bin
and lastframe_elastic.bin files to be read when running the adjoint simulation.

Third: run the adjoint simulation
Make sure that the adjoint source files absorbing boundaries and last frame files are
in the OUTPUT_FILES directory.
=> isolver = 2
=> save_forward = .false.

Output_files (for example for the elastic case)
snapshot_rho_kappa_mu*****
snapshot_rhop_alpha_beta*****
which are the primary moduli kernels and the phase velocities kernels respectively.

Note: At the moment, adjoint simulations do not support anisotropy, attenuation, and viscous damping.

--------------------------------------------------
               COUPLED SIMULATIONS
--------------------------------------------------

The code support acoustic/elastic, acoustic/poroelastic, elastic/poroelastic,
and acoustic,elastic/poroelastic simulations.

Elastic/poroelastic coupling supports anisotropy, but not attenuation for the
elastic material.

--------------------------
--------------------------
--------------------------

Regarding the structure of some of the database files:

Question: Can anyone tell me what the columns of the SPECFEM2D boundary
condition files in SPECFEM2D/DATA/Mesh_canyon are?

SPECFEM2D/DATA/Mesh_canyoncanyon_absorbing_surface_file
SPECFEM2D/DATA/Mesh_canyoncanyon_free_surface_file

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

