Extra README: 
addresses the modifications to the code to run adjoint and poroelastic simulations. [09/10/08]

---------------------------------
      NEW INPUTS IN Par_file
---------------------------------
In section "# time step parameters":
ISOLVER defines the type of simulations 
(1) forward simulation
(2) adjoint method and kernels calculation

In section "# receiver line parameters for seismograms":
SAVE_FORWARD determines if the last frame of a forward simulation is saved (.true.) or not (.false)

In section "# define models....":
Contrary to the previous version, we don't define density and velocity, but:
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
mu_s = solid shear modulus
eta_f = fluid viscosity
mu_fr = frame shear modulus

Set the porosity phi to 1 to make a given model acoustic [then edit rho_f,
kappa_f], and to 0 to make it elastic [then edit rho_s, kappa_s, mu_s]
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
absorb_elastic_bottom.bin
absorb_elastic_left.bin
absorb_elastic_right.bin
absorb_elastic_top.bin
lastframe_elastic.bin
S****.AA.BHX.semd
S****.AA.BHZ.semd

Second: define the adjoint source
Use adj_seismogram.f90 which is in UTILS/adjoint.
Edit to update NSTEP, nrec, t0, deltat, and the position of the cut to pic
any given phase if needed (tstart,tend), add the right number of stations, and
put one component of the source to zero if needed.

The ouput files are S****.AA.BHX.adj and S****.AA.BHZ.adj. They need to be
kept in the OUTPUT_FILES directory together with the absorb_elastic_****.bin
files to be read when running the "adjoint" simulation.

Third: run the "adjoint" simulation
Make sure that the adjoint source files and the absorbing boundaries files are
in the OUTPUT_FILES directory.
=> isolver = 2
=> save_forward = .false.

Output_files (for example for the elastic case)
snapshot_rho_kappa_mu*****
snapshot_rhop_alpha_beta*****
which are the moduli kernels and the phase velocities kernels respectively.
Edit and use plot_snapshot.csh located in UTILS/adjoint to generate kernels
plot.







