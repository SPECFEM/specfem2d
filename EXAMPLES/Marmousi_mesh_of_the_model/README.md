README
======

**Marmousi2 with CUBIT/Trelis mesh of the model**

This examples holds several subfolders with the idea of providing a CUBIT/Trelis mesh example of Marmousi2.
A default mesh to run the example is provided in folder `MESH-default/` that was created using the scripts in the subfolder `CUBIT_meshing/`.

To run the example, type:
```
./run_this_example.sh
```

* mesh improvements:

  For more recent meshes with manual tweaking to improve the element Jacobians, for example in subfolder
  `improved_mesh_from_Hom_Nath_Gharti_2023_using_CUBIT_higher_resolution/`, one can create a symbolic link to the provided improved mesh
  as `MESH-improved/` by:
  ```
  > cd improved_mesh_from_Hom_Nath_Gharti_2023_using_CUBIT_higher_resolution/
  > tar -xvjf Marmousi_mesh_HomNath_Gharti_higher_resolution_nov2023.tar.bz2
  > cd ../
  > ln -s improved_mesh_from_Hom_Nath_Gharti_2023_using_CUBIT_higher_resolution/Marmousi_mesh_HomNath_Gharti_higher_resolution_nov2023/
  ```
  and replace the default `Par_file` & `SOURCE` files in `DATA/` with the corresponding `Par_file.mesh-improved` & `SOURCE.mesh-improved` files.

* Marmousi2 water-layer modification:

  One might want to apply a main modification to the original Marmousi2 model, consisting of replacing the water layer with a solid layer. To do so, one can create a new `nummaterial_velocity_file_marmousi2` using the script provided in folder `CUBIT_meshing/`:
  ```
  > cd CUBIT_meshing/
  > ./5_convert_surface_rock_to_velocities.py --without-water
  ```
  The newly created `nummaterial_velocity_file_marmousi2` file replaces the waterlayer with the solid as done in Capdeville et al. (2010), *2-D non-periodic homogenization to upscale elastic media for P–SV waves*, GJI, 182, p. 903-922.

  For convenience, this file is provided in the default `MESH-default/` folder as `nummaterial_velocity_file_marmousi2_nowater`. You can adapt the `Par_file` accordingly to switch to this one.
