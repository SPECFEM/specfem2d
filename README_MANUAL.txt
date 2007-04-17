
To use the code:

- edit the Makefile

- type "make all"

- edit the input file "DATA/Par_file" which describes the simulation. It contains comments and should be almost self-explanatory, if you need more details we do not have a manual for the 2D version but you can find useful information in the manuals of the 3D versions, since many parameters and the general philosophy is similar. They are available at http://www.univ-pau.fr/~dkomati1/published_papers/manual_SPECFEM3D_GLOBE.pdf and http://www.univ-pau.fr/~dkomati1/published_papers/manual_SPECFEM3D_BASIN.pdf. To create acoustic (fluid) regions, just set the S wave speed to zero and the code will see that these elements are fluid and switch to the right equations there automatically, and automatically match them with the solid regions

- then type xmeshfem2D to create the mesh (which will be stored in directory OUTPUT_FILES)

- then type xspecfem2D to run the main solver. This will output the seismograms and snapshots of the wave fronts at different time steps in directory OUTPUT_FILES. To visualize them, type "gs OUTPUT_FILES/vect*.ps" to see the Postscript files (in which the wave field is represented with small arrows, fluid/solid matching interfaces with a thick pink line, and absorbing edges with a thick green line) and "gimp OUTPUT_FILES/image*.gif" to see the color snapshot showing a pixelized image of one of the two components of the wave field (or pressure, depending on what you have selected for the output in DATA/Par_file)

- the DATA/Par_file given with the code works fine, you can use it without any modification to test the code

- the seismograms OUTPUT_FILES/*.sem* are simple ASCII files with two columns: time in the first colum and amplitude in the second, therefore they can be visualized with any tool you like, for instance "gnuplot"

- you can convolve them with any source time function in postprocessing later using "convolve_source_timefunction.csh" and "convolve_source_timefunction.f90", see the manual of the 3D code for details on how to do this

- we do not have PML absorbing conditions implemented in the fluid/solid code yet. We use (older and less efficient) paraxial Clayton-Engquist or Sommerfeld equations instead. This is only by lack of time, I have a student who is currently implementing PML but the code is not fully ready. I will send it to you when it is. (We already have PML in the purely elastic code, see http://www.univ-pau.fr/~dkomati1/published_papers/pml_2nd_order_GJI_typos_fixed.pdf for details, therefore it is only a matter of cutting/pasting the routines). For now, since the paraxial conditions are less efficient, please use a larger model until we send you the code with PML

- I have made many modifications in the source code recently, I have carefully tested them but if you find bugs please send me an email and I will fix them and send you the patch

