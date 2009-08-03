Extra README: 
addresses the modifications to the code to run adjoint and poroelastic simulations. [09/10/08]
updated [07/30/09]

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

elastic/poroelastic coupling support anisotropy, but not attenuation for the
elastic material.



