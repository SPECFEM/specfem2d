
rm -f xspecfem2D* *.o libjpeg/*.o *__genmod.*

# Intel ifort compiler

#cd libjpeg
#icc -c -O3 *.c
#cd ..
#icc -c -O3 -o write_jpeg_image.o write_jpeg_image.c

# with full optimization on
# do NOT suppress -ftz, which is important for performance
#ifort -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -vec-report0 -std03 -implicitnone -gen-interfaces -warn all -O3 -check nobounds -o xspecfem2D specfem2D_plane_strain.f90 plot_post.f90 createnum_slow.f90 define_derivation_matrices.f90 recompute_jacobian.f90 define_shape_functions.f90 gll_library.f90 lagrange_poly.f90 create_color_image.f90 write_color_image_snaphot.f90 prepare_color_image.f90 is_in_convex_quadrilateral.f90 write_jpeg_image.o libjpeg/*.o

# with full debugging turned on
#ifort -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence -vec-report0 -std03 -implicitnone -gen-interfaces -warn all -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv -o xspecfem2D specfem2D_plane_strain.f90 plot_post.f90 createnum_slow.f90 define_derivation_matrices.f90 recompute_jacobian.f90 define_shape_functions.f90 gll_library.f90 lagrange_poly.f90 create_color_image.f90 write_color_image_snaphot.f90 prepare_color_image.f90 is_in_convex_quadrilateral.f90 write_jpeg_image.o libjpeg/*.o

# other compilers

cd libjpeg
gcc -c -O3 *.c
cd ..
gcc -c -O3 -o write_jpeg_image.o write_jpeg_image.c

# with full optimization on
gfortran -O3 -std=f2003 -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused -Werror -o xspecfem2D specfem2D_plane_strain.f90 plot_post.f90 createnum_slow.f90 define_derivation_matrices.f90 recompute_jacobian.f90 define_shape_functions.f90 gll_library.f90 lagrange_poly.f90 create_color_image.f90 write_color_image_snaphot.f90 prepare_color_image.f90 is_in_convex_quadrilateral.f90 write_jpeg_image.o libjpeg/*.o

# with full debugging turned on
#gfortran -std=f2003 -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused -Werror -g -O0 -ggdb -fbacktrace -fbounds-check -o xspecfem2D specfem2D_plane_strain.f90 plot_post.f90 createnum_slow.f90 define_derivation_matrices.f90 recompute_jacobian.f90 define_shape_functions.f90 gll_library.f90 lagrange_poly.f90 create_color_image.f90 write_color_image_snaphot.f90 prepare_color_image.f90 is_in_convex_quadrilateral.f90 write_jpeg_image.o libjpeg/*.o

rm -f *.o libjpeg/*.o *__genmod.*

