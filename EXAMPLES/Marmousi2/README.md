Readme
======

This example is demonstrating how to use Marmousi2 based on a tomography file.
The meshing is done by the internal mesher, using a regular mesh with an element size of ~50m.
You can change this mesh by modifying the corresponding section in the `DATA/Par_file`.


## Tomography file

  We provide a Marmousi2 tomography file in `DATA/tomography_model_marmousi2.xyz.bin` that was created
  using the script `interpolate_Marmousi_2_tomo.py`:
  ```
  ./interpolate_Marmousi_2_tomo.py --smooth=5.0 --replace-water --binary
  ```

  This script downloads the original, gridded Marmousi2 SEGY-files and together with the above options, it smooths the original model
  with a Gaussian kernel corresponding to a 5 Hz signal wavelength. Furthermore, it replaces the top water layer with the solid velocities
  of the upper-most elastic layer. This avoids having problems with reading in the tomography file in SPECFEM2D as well as smearing out
  the fluid-solid interface by the interpolation of the tomography model velocities onto the mesh.


## Mesh setup

  The default mesh setup in the `DATA/Par_file` creates two layers, the top water layer and the solid elastic layer below.
  For the water layer, the velocities are set to the original water velocities from Marmousi2.
  The solid, elastic layer below defines a tomographic material. This will interpolate the tomography file velocities from Marmousi2
  onto all the GLL points in this solid layer.

  In case you want the water layer removed, you can choose:
  - define the whole mesh as a single solid layer with a tomographic material. This will interpolate the provided tomography file, which has
    the water replaced by the top-most elastic velocities, onto the whole region.
    Thus, you would only need to modify in `DATA/Par_file` the `nbmodels` section.

  - adapt the mesh position such that it falls within the solid domain. Marmousi2 has a top water layer with a depth of 450 m.
    Thus, you can specify the mesh to lie within a vertical range of 0 to 3050 m to be in the solid part only.
    For this, you would modify in `DATA/Par_file` the `nbmodels` and `nbregions` sections, as well as changing
    the `DATA/interfaces.dat` file accordingly.

There are many ways how to run external models, and in particular how to interpolate external model velocities onto a SPECFEM2D grid.
This example demonstrates the use of a tomography file to achieve this, but might not be the best option for your own use cases.

Please consider contributing more examples to this package!


## References

Martin, G. S., 2004,
*The Marmousi2 model, elastic synthetic data, and an analysis of imaging and AVO in a structurally complex environment*:
Master’s thesis. University of Houston. Retrieved from www.agl.uh.edu/pdf/theses/martin.pdf

Martin, G. S., Marfurt, K. J., and Larsen, S., 2005,
*Marmousi‐2: An updated model for the investigation of AVO in structurally complex areas*:
SEG Technical Program Expanded Abstracts 2002, 1979–1982. doi:10.1190/1.1817083

Martin, G. S., Wiley, R., and Marfurt, K. J., 2006,
*Marmousi2: An elastic upgrade for Marmousi*:
The Leading Edge, 25, 156–166. doi:10.1190/1.2172306
