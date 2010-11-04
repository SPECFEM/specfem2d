
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

!====================================================================================
!
!   An explicit 2D parallel MPI spectral element solver
!   for the anelastic anisotropic or poroelastic wave equation.
!
!====================================================================================

! If you use this code for your own research, please cite at least one article
! written by the developers of the package, for instance:
!
! @ARTICLE{TrKoLi08,
! author = {Jeroen Tromp and Dimitri Komatitsch and Qinya Liu},
! title = {Spectral-Element and Adjoint Methods in Seismology},
! journal = {Communications in Computational Physics},
! year = {2008},
! volume = {3},
! pages = {1-32},
! number = {1}}
!
! or
!
! @ARTICLE{VaCaSaKoVi99,
! author = {R. Vai and J. M. Castillo-Covarrubias and F. J. S\'anchez-Sesma and
! D. Komatitsch and J. P. Vilotte},
! title = {Elastic wave propagation in an irregularly layered medium},
! journal = {Soil Dynamics and Earthquake Engineering},
! year = {1999},
! volume = {18},
! pages = {11-18},
! number = {1},
! doi = {10.1016/S0267-7261(98)00027-X}}
!
! @ARTICLE{LeChKoHuTr09,
! author = {Shiann Jong Lee and Yu Chang Chan and Dimitri Komatitsch and Bor
! Shouh Huang and Jeroen Tromp},
! title = {Effects of realistic surface topography on seismic ground motion
! in the {Y}angminshan region of {T}aiwan based upon the spectral-element
! method and {LiDAR DTM}},
! journal = {Bull. Seismol. Soc. Am.},
! year = {2009},
! volume = {99},
! pages = {681-693},
! number = {2A},
! doi = {10.1785/0120080264}}
!
! @ARTICLE{LeChLiKoHuTr08,
! author = {Shiann Jong Lee and How Wei Chen and Qinya Liu and Dimitri Komatitsch
! and Bor Shouh Huang and Jeroen Tromp},
! title = {Three-Dimensional Simulations of Seismic Wave Propagation in the
! {T}aipei Basin with Realistic Topography Based upon the Spectral-Element Method},
! journal = {Bull. Seismol. Soc. Am.},
! year = {2008},
! volume = {98},
! pages = {253-264},
! number = {1},
! doi = {10.1785/0120070033}}
!
! @ARTICLE{LeKoHuTr09,
! author = {S. J. Lee and Dimitri Komatitsch and B. S. Huang and J. Tromp},
! title = {Effects of topography on seismic wave propagation: An example from
! northern {T}aiwan},
! journal = {Bull. Seismol. Soc. Am.},
! year = {2009},
! volume = {99},
! pages = {314-325},
! number = {1},
! doi = {10.1785/0120080020}}
!
! @ARTICLE{KoErGoMi10,
! author = {Dimitri Komatitsch and Gordon Erlebacher and Dominik G\"oddeke and
! David Mich\'ea},
! title = {High-order finite-element seismic wave propagation modeling with
! {MPI} on a large {GPU} cluster},
! journal = {J. Comput. Phys.},
! year = {2010},
! volume = {229},
! pages = {7692-7714},
! number = {20},
! doi = {10.1016/j.jcp.2010.06.024}}
!
! @ARTICLE{KoGoErMi10,
! author = {Dimitri Komatitsch and Dominik G\"oddeke and Gordon Erlebacher and
! David Mich\'ea},
! title = {Modeling the propagation of elastic waves using spectral elements
! on a cluster of 192 {GPU}s},
! journal = {Computer Science Research and Development},
! year = {2010},
! volume = {25},
! pages = {75-82},
! number = {1-2},
! doi = {10.1007/s00450-010-0109-1}}
!
! @ARTICLE{KoMiEr09,
! author = {Dimitri Komatitsch and David Mich\'ea and Gordon Erlebacher},
! title = {Porting a high-order finite-element earthquake modeling application
! to {NVIDIA} graphics cards using {CUDA}},
! journal = {Journal of Parallel and Distributed Computing},
! year = {2009},
! volume = {69},
! pages = {451-460},
! number = {5},
! doi = {10.1016/j.jpdc.2009.01.006}}
!
! @ARTICLE{LiPoKoTr04,
! author = {Qinya Liu and Jascha Polet and Dimitri Komatitsch and Jeroen Tromp},
! title = {Spectral-element moment tensor inversions for earthquakes in {S}outhern {C}alifornia},
! journal={Bull. Seismol. Soc. Am.},
! year = {2004},
! volume = {94},
! pages = {1748-1761},
! number = {5},
! doi = {10.1785/012004038}}
!
! @INCOLLECTION{ChKoViCaVaFe07,
! author = {Emmanuel Chaljub and Dimitri Komatitsch and Jean-Pierre Vilotte and
! Yann Capdeville and Bernard Valette and Gaetano Festa},
! title = {Spectral Element Analysis in Seismology},
! booktitle = {Advances in Wave Propagation in Heterogeneous Media},
! publisher = {Elsevier - Academic Press},
! year = {2007},
! editor = {Ru-Shan Wu and Val\'erie Maupin},
! volume = {48},
! series = {Advances in Geophysics},
! pages = {365-419}}
!
! @ARTICLE{KoVi98,
! author={D. Komatitsch and J. P. Vilotte},
! title={The spectral-element method: an efficient tool to simulate the seismic response of 2{D} and 3{D} geological structures},
! journal={Bull. Seismol. Soc. Am.},
! year=1998,
! volume=88,
! number=2,
! pages={368-392}}
!
! @ARTICLE{KoTr99,
! author={D. Komatitsch and J. Tromp},
! year=1999,
! title={Introduction to the spectral-element method for 3-{D} seismic wave propagation},
! journal={Geophys. J. Int.},
! volume=139,
! number=3,
! pages={806-822},
! doi={10.1046/j.1365-246x.1999.00967.x}}
!
! @ARTICLE{KoLiTrSuStSh04,
! author={Dimitri Komatitsch and Qinya Liu and Jeroen Tromp and Peter S\"{u}ss
!   and Christiane Stidham and John H. Shaw},
! year=2004,
! title={Simulations of Ground Motion in the {L}os {A}ngeles {B}asin
!   based upon the Spectral-Element Method},
! journal={Bull. Seism. Soc. Am.},
! volume=94,
! number=1,
! pages={187-206}}
!
! @ARTICLE{MoTr08,
! author={C. Morency and J. Tromp},
! title={Spectral-element simulations of wave propagation in poroelastic media},
! journal={Geophys. J. Int.},
! year=2008,
! volume=175,
! pages={301-345}}
!
! and/or other articles from http://web.univ-pau.fr/~dkomati1/publications.html
!
! If you use the kernel capabilities of the code, please cite at least one article
! written by the developers of the package, for instance:
!
! @ARTICLE{TrKoLi08,
! author = {Jeroen Tromp and Dimitri Komatitsch and Qinya Liu},
! title = {Spectral-Element and Adjoint Methods in Seismology},
! journal = {Communications in Computational Physics},
! year = {2008},
! volume = {3},
! pages = {1-32},
! number = {1}}
!
! or
!
! @ARTICLE{LiTr06,
! author={Qinya Liu and Jeroen Tromp},
! title={Finite-frequency kernels based on adjoint methods},
! journal={Bull. Seismol. Soc. Am.},
! year=2006,
! volume=96,
! number=6,
! pages={2383-2397},
! doi={10.1785/0120060041}}
!
! @ARTICLE{MoLuTr09,
! author={C. Morency and Y. Luo and J. Tromp},
! title={Finite-frequency kernels for wave propagation in porous media based upon adjoint methods},
! year=2009,
! journal={Geophys. J. Int.},
! doi={10.1111/j.1365-246X.2009.04332}}
!
! If you use the SCOTCH / CUBIT non-structured capabilities, please also cite:
!
! @ARTICLE{MaKoBlLe08,
! author = {R. Martin and D. Komatitsch and C. Blitz and N. {Le Goff}},
! title = {Simulation of seismic wave propagation in an asteroid based upon
! an unstructured {MPI} spectral-element method: blocking and non-blocking
! communication strategies},
! journal = {Lecture Notes in Computer Science},
! year = {2008},
! volume = {5336},
! pages = {350-363}}
!


SPECFEM2D 6.1: SPECFEM2D facilitates 2D simulations of
        acoustic, (an)elastic, and poroelastic seismic wave propagation.
        With version 6.1, the 2D spectral-element solver accommodates
        regular and unstructured meshes, generated for example by Cubit
        (http://cubit.sandia.gov), Gmsh (http://geuz.org/gmsh)
        or GiD (http://www.gid.cimne.upc.es).
        The solver has adjoint capabilities and can
        calculate finite-frequency sensitivity kernels for acoustic,
        (an)elastic, and poroelastic media. The package also considers 2D SH
        and P-SV wave propagation. Finally, the solver can run
        both in serial and in parallel. See SPECFEM2D
        <http://www.geodynamics.org/cig/software/packages/seismo/specfem2d>
        for the source code.

The SEM is a continuous Galerkin technique, which can easily be made discontinous;
it is then a particular case of the discontinuous Galerkin technique,
with optimized efficiency because of its tensorized basis functions.
Note that in many (most?) geological models in the context of seismic wave propagation studies
a discontinous mesh is not needed because material property contrasts are not drastic and thus a continuous formulation is sufficient.

Caution:
--------

- the units for the components of a moment tensor source are different
in SPECFEM2D and in SPECFEM3D:
 - in SPECFEM3D the moment tensor components are in dyne*cm
 - in SPECFEM2D the moment tensor components are in N*m

How to use SPECFEM2D:
---------------------

See file "todo_list_please_dont_remove.txt" for a list of known bugs, problems, or missing options.

To use the code:

- edit the Makefile. There are several options available : -DUSE_MPI compiles with use of an MPI library. -DUSE_SCOTCH enables use of graph partitioner SCOTCH.

- if you want to run in parallel, i.e., using more than one processor core, then you need to compile the SCOTCH graph partitioner. Go to subdirectory scotch_5.1.10b and read INSTALL.txt. You may want to download more recent versions of SCOTCH in the future from http://www.labri.fr/perso/pelegrin/scotch/scotch_en.html . Support for the METIS graph partitioner has been discontinued because SCOTCH is more recent and performs better.

- then go back to the main directory of SPECFEM2D and type "make all"

- edit the input file "DATA/Par_file" which describes the simulation. It contains comments and should be almost self-explanatory, if you need more details we do not have a manual for the 2D version but you can find useful information in the manuals of the 3D versions, since many parameters and the general philosophy is similar. They are available at http://geodynamics.org/wsvn/cig/seismo/3D in subdirectories USER_MANUAL. To create acoustic (fluid) regions, just set the S wave speed to zero and the code will see that these elements are fluid and switch to the right equations there automatically, and automatically match them with the solid regions

- if you are using an external mesher (like GID or CUBIT), you should set "read_external_mesh" to true.
     "mesh_file" is the file describing the mesh : first line is the number of elements, then a list of 4 nodes (quadrilaterals only) forming each elements on each line.
     "nodes_coords_file" is the file containing the coordinates (x and z) of each nodes : number of nodes on the first line, then coordinates x and z on each line.
     "materials_file" is the number of the material for every elements : an integer ranging from 1 to nbmodels on each line.
     "free_surface_file" is the file describing the edges forming the acoustic free surface : number of edges on the first line, then on each line number of the element, number of nodes forming the free surface (1 for a point, 2 for an edge), the nodes forming the free surface for this element. If you do not want free surface, jusr put 0 on the first line.
     "absorbing_surface_file" is the file describing the edges forming the absorbing boundaries : the format is the same as the "free_surface_file".
     "tangential_detection_curve_file" contains points describing the envelope, used for source_normal_to_surface and rec_normal_to_surface. Should be fine grained, and ordained clockwise. Number of points on the first line, then (x,z) coordinates on each line.

- if you have compiled with MPI, you must specify the number of processes.

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
SIMULATION_TYPE defines the type of simulations
(1) forward simulation
(2) adjoint method and kernels calculation

In section "# source parameters":
The code now support multi sources.
NSOURCE is the number of source.
Parameters of the sources are displayed in the file SOURCE, which must be
in the directory DATA. The components of a moment tensor source must be given in N.m,
not in dyne.cm as in the DATA/CMTSOLUTION source file of the 3D version of the code.

In section "# receiver line parameters for seismograms":
SAVE_FORWARD determines if the last frame of a forward simulation is saved (.true.) or not (.false)

In section "# define models....":
Three types of models:
I: (model_number,1,rho,Vp,Vs,0,0,Qp,Qs,0,0,0,0,0,0), for isotropic elastic/acoustic
material
or II: (model_number,2,rho,Vp,Vs,c11,c13,c15,c33,c35,c55,Qp,Qs,0,0), for anisotropic material
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

get_poroelastic_velocities.f90 allows to compute cpI, cpII, and cs function of
the source dominant frequency. Notice that for this calculation we use permxx
and the dominant frequency of the first source , f0(1). Caution if you use
several sources with different frequencies and if you consider anistropic
permeability.

--------------------------------------------------
     HOW TO OBTAIN FINITE SENSITIVITY KERNELS
--------------------------------------------------

First: run a forward simulation
=> SIMULATION_TYPE = 1
=> SAVE_FORWARD = .true.
=> seismotype = 1 (we need to save the displacement fields to later on derive the
adjoint source. Note: if the user forgets it, the program corrects it when reading the proper
SIMULATION_TYPE/SAVE_FORWARD combination and a warning message appears in the ouput
file)

Important output files (for example, for the elastic case, P-SV waves)
absorb_elastic_bottom*****.bin
absorb_elastic_left*****.bin
absorb_elastic_right*****.bin
absorb_elastic_top*****.bin
lastframe_elastic*****.bin
S****.AA.BHX.semd
S****.AA.BHZ.semd

Second: define the adjoint source
Use adj_seismogram.f90
Edit to update NSTEP, nrec, t0, deltat, and the position of the cut to pic
any given phase if needed (tstart,tend), add the right number of stations, and
put one component of the source to zero if needed.
The ouput files of adj_seismogram.f90 are S****.AA.BHX.adj and S****.AA.BHZ.adj, for P-SV waves (and
S****.AA.BHY.adj, for SH (membrane) waves). Note that you will need these three
files (S****.AA.BHX.adj, S****.AA.BHY.adj and S****.AA.BHZ.adj) to be present in the OUTPUT_FILES directory
together with the absorb_elastic_****.bin and lastframe_elastic.bin files to be read
when running the adjoint simulation.

Third: run the adjoint simulation
Make sure that the adjoint source files absorbing boundaries and last frame files are
in the OUTPUT_FILES directory.
=> SIMULATION_TYPE = 2
=> SAVE_FORWARD = .false.

Output_files (for example for the elastic case)
snapshot_rho_kappa_mu*****
snapshot_rhop_alpha_beta*****
which are the primary moduli kernels and the phase velocities kernels respectively, in ascii format
and at the local level, that is as "kernels(i,j,ispec)".

Note1: At the moment, adjoint simulations do not support anisotropy, attenuation, and viscous damping.
Note2: You will need S****.AA.BHX.adj, S****.AA.BHY.adj and S****.AA.BHZ.adj
to be present in OUTPUT_FILES even if you are just running an acoustic or
poroelastic adjoint simulation.
S****.AA.BHX.adj is the only relevant component for an acoustic case.
S****.AA.BHX.adj and S****.AA.BHZ.adj are the only relevant components for a
poroelastic case.

--------------------------------------------------
               COUPLED SIMULATIONS
--------------------------------------------------

The code supports acoustic/elastic, acoustic/poroelastic, elastic/poroelastic,
and acoustic,elastic/poroelastic simulations.

Elastic/poroelastic coupling supports anisotropy, but not attenuation for the
elastic material.


How to run P-SV or SH (membrane) wave simulations :
---------------------------------------------------

To run a P-SV waves calculation propagating in the x-z plane,
set p_sv = .true. in the Par_file.
To run a SH (membrane) waves calculation traveling in the x-z plane with a
y-component of motion, set p_sv = .false.

This feature is only implemented for elastic materials and sensitivity kernels
can be calculated (see Tape, Liu & Tromp, GJI 2006 for details on membrane
surface waves).


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

