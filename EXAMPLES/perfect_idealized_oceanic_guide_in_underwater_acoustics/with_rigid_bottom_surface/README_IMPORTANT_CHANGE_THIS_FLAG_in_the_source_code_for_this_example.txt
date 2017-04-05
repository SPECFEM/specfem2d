
IMPORTANT: this example requires changing the following flag:

 logical, parameter :: ENFORCE_RIGID_SURFACE_BOTTOM = .false.

from .false. to .true. in file src/specfem2D/enforce_acoustic_free_surface.f90
and then recompiling the code ("make clean ; make all") for this example to work,
otherwise it will run but will give incorrect results.

If you have any questions feel free to contact the author of this example: Laurent dot Guillon at ecole-navale dot fr

