
rm -f xspecfem2D* *.o *__genmod.*

# Intel ifort compiler
# do NOT suppress -ftz, which is important for performance
#ifort -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -vec-report0 -std03 -implicitnone -gen-interfaces -warn all -O3 -check nobounds -o xspecfem2D specfem2D_plane_strain.f90 plot_post.f90 createnum_slow.f90 define_derivation_matrices.f90 recompute_jacobian.f90 define_shape_functions.f90 gll_library.f90 lagrange_poly.f90

# with full debugging turned on
#ifort -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -vec-report0 -std03 -implicitnone -gen-interfaces -warn all -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv -o xspecfem2D specfem2D_plane_strain.f90 plot_post.f90 createnum_slow.f90 define_derivation_matrices.f90 recompute_jacobian.f90 define_shape_functions.f90 gll_library.f90 lagrange_poly.f90

# other compilers

gfortran -O3 -std=f2003 -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused -Werror -o xspecfem2D specfem2D_plane_strain.f90 plot_post.f90 createnum_slow.f90 define_derivation_matrices.f90 recompute_jacobian.f90 define_shape_functions.f90 gll_library.f90 lagrange_poly.f90

# with full debugging turned on
#gfortran -std=f2003 -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused -Werror -g -O0 -ggdb -fbacktrace -fbounds-check -o xspecfem2D specfem2D_plane_strain.f90 plot_post.f90 createnum_slow.f90 define_derivation_matrices.f90 recompute_jacobian.f90 define_shape_functions.f90 gll_library.f90 lagrange_poly.f90

#pgf90 -fast -Mnobounds -Mneginfo -Mdclchk -Knoieee -Minform=warn -Mstandard -fastsse -tp amd64e -o xspecfem2D specfem2D_plane_strain.f90 plot_post.f90 createnum_slow.f90 define_derivation_matrices.f90 recompute_jacobian.f90 define_shape_functions.f90 gll_library.f90 lagrange_poly.f90

#xlf_r -O3 -qsave -qstrict -qtune=ppc970 -qarch=ppc64v -qcache=auto -qfree=f90 -Q -qflttrap=en:ov:zero:inv -o xspecfem2D specfem2D_plane_strain.f90 plot_post.f90 createnum_slow.f90 define_derivation_matrices.f90 recompute_jacobian.f90 define_shape_functions.f90 gll_library.f90 lagrange_poly.f90

rm -f *.o *__genmod.*

