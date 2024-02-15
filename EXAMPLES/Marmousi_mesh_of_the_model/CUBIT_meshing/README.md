README
======

**Scripts to mesh Marmousi model:**

The scripts here can be used with CUBIT/Trelis to create a mesh for the Marmousi model.
Just run in consecutive order, or open in CUBIT/Trelis -> run script:

  ```
  ./1_create_surfaces.py
  ```
  this will create a CUBIT  model of the Marmousi model using all vertices from the horizons to create mesh surfaces.
  ```
  ./2_mesh_marmousi.py
  ```
  this will mesh the surfaces. due to some small edges for a few surfaces, this initial mesh needs to be tweaked a bit in a subsequent step.
  ```
  ./3_smooth_mesh.py
  ```
  this step will try to improve the mesh by running a smoothing operation on the mesh surface - takes a while to finish...
  ```
  ./4_export_mesh.py
  ```
  as it says, it will export the CUBIT mesh into SPECFEM readible files stored in folder MESH/


  For assigning velocities to the surfaces, run script:
  ```
  ./5_convert_surface_rock_to_velocities.py
  ```
  to create a file nummaterial_velocity_file_marmousi2. This will assign the original water layer to the mesh. 
  To replace the water layer with solid velocities as in Capdeville et al. 2010, use
  ```
  ./5_convert_surface_rock_to_velocities.py --without-water
  ```

**Installation:**

To run these python scripts in the command line, instead of opening the Cubit application and use the "run script"-Button,
you will need to make sure that within the python environment you can load the "cubit" module.

The most recent Coreform-Cubit version is 2023.11. It uses internally a python3.10 version. 
The details below describe how such a installation setup could look like for MacOS.
 
* for MacOS:
   it requires to run this python script with the shipped Cubit version python
   (due to different architecture compilations).

   Cubit Coreform version 2023.11:

   - the cubit python version is here
     `/Applications/Coreform-Cubit-2023.11.app/Contents/lib/python3/Python.framework/Versions/3.10/python3.10`


   - either create an alias which can be set in `~/.alias` or `~/.bashrc`:
     ```
     export CUBITDIR=/Applications/Coreform-Cubit-2023.11.app
     function python-cubit() { ${CUBITDIR}/Contents/lib/python3/Python.framework/Versions/3.10/python3.10 "$@" ; }
     export -f python-cubit
     ```

     or create a symbolic link to the python version in `/usr/local/bin`:
     ``` 
     sudo ln -s /Applications/Coreform-Cubit-2023.11.app/Contents/lib/python3/Python.framework/Versions/3.10/python3.10 python-cubit
     ```

     I ended up using this latter approach, as the export option didn't work properly.

   - to import the cubit module, this also requires to add the Cubit library path to the `PYTHONPATH` environment variable, in `~/.bashrc` use:
     ```
     export PYTHONPATH=$PYTHONPATH:${CUBITDIR}/Contents/lib/
     ``` 


