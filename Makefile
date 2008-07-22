
#========================================================================
#
#                   S P E C F E M 2 D  Version 5.2
#                   ------------------------------
#
# Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
# Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
#               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
#               Roland Martin, roland DOT martin aT univ-pau DOT fr
#
# This software is a computer program whose purpose is to solve
# the two-dimensional viscoelastic anisotropic wave equation
# using a spectral-element method (SEM).
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The full text of the license is available in file "LICENSE".
#
#========================================================================

SHELL=/bin/sh

# uncomment this to generate ParaVer traces on MareNostrum in Barcelona
#MPITRACE_HOME = /gpfs/apps/CEPBATOOLS/mpitrace-devel/64
#PAPI_HOME = /gpfs/apps/PAPI/3.2.1-970mp/64
#PERFCTR_HOME  = /gpfs/apps/PAPI/papi-3.2.1-970mp/64

O = obj

# Portland
#F90 = pgf90
#F90 = /opt/openmpi-1.2.2/pgi64/bin/mpif90 -DUSE_MPI -DUSE_METIS -DUSE_SCOTCH
#CC = pgcc
#FLAGS_NOCHECK=-fast -Mnobounds -Minline -Mneginfo -Mdclchk -Knoieee -Minform=warn -fastsse -tp amd64e -Msmart
#FLAGS_CHECK=-fast -Mbounds -Mneginfo -Mdclchk -Minform=warn

# Intel
# NOTE FOR USERS OF IFORT 10.0 AND ABOVE :
# Use of option -heap-arrays <size> can be required, depending on the size of the simulation. 
# Another workaround can be to increase your stack size (ulimit -s).
#F90 = ifort
#CC = gcc
#FLAGS_NOCHECK=-O0 -implicitnone -warn stderrors -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -assume byterecl -check nobounds
#FLAGS_CHECK = $(FLAGS_NOCHECK)

# GNU gfortran
F90 = gfortran
#F90 = mpif90 -DUSE_MPI -DUSE_METIS -DUSE_SCOTCH
#F90 = /opt/openmpi-1.2.1/gfortran64/bin/mpif90 -DUSE_MPI -DUSE_METIS -DUSE_SCOTCH
CC = gcc
#FLAGS_NOCHECK = -O3 -march=opteron -m64 -mfpmath=sse,387
FLAGS_NOCHECK = -std=gnu -fimplicit-none -frange-check -O3 -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -fno-trapping-math # -mcmodel=medium
FLAGS_CHECK = $(FLAGS_NOCHECK) -fbounds-check

# IBM
#####F90 = xlf_r
#F90 = mpif90 -WF,-DUSE_MPI,-DUSE_METIS
#CC = xlc -g -q64
# uncomment this to generate ParaVer traces on MareNostrum in Barcelona
#FLAGS_NOCHECK_ADD = -L$(MPITRACE_HOME)/lib -lmpitracef -lxml2 -L${PAPI_HOME}/lib -lpapi -lperfctr
#FLAGS_NOCHECK = $(FLAGS_NOCHECK_ADD) -qextname=attenuation_compute_param -O3 -qstrict -q64 -qtune=ppc970 -qarch=ppc970 -qcache=auto -qfree=f90 -Q -qsuffix=f=f90 -qhalt=w -qflttrap=overflow:zerodivide:invalid:enable -qfullpath  
#####FLAGS_NOCHECK = $(FLAGS_NOCHECK_ADD) -qextname=attenuation_compute_param -O3 -qstrict -q64 -qtune=ppc970 -qarch=ppc970 -qcache=auto -qfree=f90 -qsuffix=f=f90 -qhalt=w -qflttrap=overflow:zerodivide:invalid:enable -qinitauto=7FBFFFFF -C # -qlanglvl=2003pure
#####FLAGS_NOCHECK = $(FLAGS_NOCHECK_ADD) -qextname=attenuation_compute_param -O0 -q64 -qtune=ppc970 -qarch=ppc970 -qcache=auto -qfree=f90 -qsuffix=f=f90 -qhalt=w -qflttrap=overflow:zerodivide:invalid:enable -qinitauto=7FBFFFFF -C -g -qfullpath -qlinedebug
#FLAGS_CHECK = $(FLAGS_NOCHECK)

#LIB = /opt/metis-4.0/gcc64/lib/libmetis.a /opt/scotch-4.0/gcc64/lib/libscotch.a  /opt/scotch-4.0/gcc64/lib/libscotcherr.a
# uncomment this to use Metis on MareNostrum in Barcelona
#LIB = /home/hpce08/hpce08548/utils/metis-4.0/libmetis.a

LINK = $(F90)

OBJS_MESHFEM2D = $O/part_unstruct.o $O/meshfem2D.o $O/read_value_parameters.o $O/spline_routines.o

OBJS_SPECFEM2D = $O/checkgrid.o $O/datim.o $O/enforce_acoustic_free_surface.o\
        $O/compute_forces_acoustic.o $O/compute_forces_elastic.o\
        $O/lagrange_poly.o $O/gmat01.o $O/gll_library.o $O/plotgll.o $O/define_derivation_matrices.o\
        $O/plotpost.o $O/locate_receivers.o $O/locate_source_force.o $O/compute_gradient_attenuation.o\
        $O/specfem2D.o $O/write_seismograms.o $O/define_external_model.o $O/createnum_fast.o $O/createnum_slow.o\
        $O/define_shape_functions.o $O/attenuation_model.o $O/create_color_image.o $O/compute_vector_field.o $O/compute_pressure.o\
        $O/recompute_jacobian.o $O/compute_arrays_source.o $O/locate_source_moment_tensor.o $O/netlib_specfun_erf.o\
        $O/construct_acoustic_surface.o $O/assemble_MPI.o $O/compute_energy.o $O/compute_curl_one_element.o\
        $O/attenuation_compute_param.o $O/compute_Bielak_conditions.o $O/paco_beyond_critical.o\
        $O/paco_convolve_fft.o $O/is_in_convex_quadrilateral.o $O/get_perm_cuthill_mckee.o

default: clean meshfem2D specfem2D convolve_source_timefunction

all: default

clean:
	/bin/rm -r -f xmeshfem2D xmeshfem2D.trace xspecfem2D xspecfem2D.trace $O/*.o *.o $O/*.il *.mod core xconvolve_source_timefunction *.oo *.ipo

meshfem2D: $(OBJS_MESHFEM2D)
	$(LINK) $(FLAGS_CHECK) -o xmeshfem2D $(OBJS_MESHFEM2D) $(LIB)

### use optimized compilation option for solver only
specfem2D: $(OBJS_SPECFEM2D)
	$(LINK) $(FLAGS_NOCHECK) -o xspecfem2D $(OBJS_SPECFEM2D)

convolve_source_timefunction: $O/convolve_source_timefunction.o
	${F90} $(FLAGS_CHECK) -o xconvolve_source_timefunction $O/convolve_source_timefunction.o

$O/checkgrid.o: checkgrid.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/checkgrid.o checkgrid.F90
    
$O/meshfem2D.o: meshfem2D.F90
	${F90} $(FLAGS_CHECK) -c -o $O/meshfem2D.o meshfem2D.F90

$O/createnum_fast.o: createnum_fast.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/createnum_fast.o createnum_fast.f90
    
$O/createnum_slow.o: createnum_slow.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/createnum_slow.o createnum_slow.f90
    
$O/convolve_source_timefunction.o: convolve_source_timefunction.f90
	${F90} $(FLAGS_CHECK) -c -o $O/convolve_source_timefunction.o convolve_source_timefunction.f90

$O/read_value_parameters.o: read_value_parameters.f90
	${F90} $(FLAGS_CHECK) -c -o $O/read_value_parameters.o read_value_parameters.f90

$O/datim.o: datim.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/datim.o datim.f90
    
$O/lagrange_poly.o: lagrange_poly.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/lagrange_poly.o lagrange_poly.f90
    
$O/gmat01.o: gmat01.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/gmat01.o gmat01.f90
    
$O/gll_library.o: gll_library.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/gll_library.o gll_library.f90
    
$O/define_derivation_matrices.o: define_derivation_matrices.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/define_derivation_matrices.o define_derivation_matrices.f90
    
$O/plotgll.o: plotgll.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/plotgll.o plotgll.f90
    
$O/plotpost.o: plotpost.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/plotpost.o plotpost.F90
    
$O/locate_receivers.o: locate_receivers.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/locate_receivers.o locate_receivers.F90
    
$O/recompute_jacobian.o: recompute_jacobian.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/recompute_jacobian.o recompute_jacobian.f90
    
$O/locate_source_force.o: locate_source_force.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/locate_source_force.o locate_source_force.F90
    
$O/locate_source_moment_tensor.o: locate_source_moment_tensor.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/locate_source_moment_tensor.o locate_source_moment_tensor.F90
    
$O/define_shape_functions.o: define_shape_functions.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/define_shape_functions.o define_shape_functions.f90
    
$O/attenuation_model.o: attenuation_model.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/attenuation_model.o attenuation_model.f90
    
### use optimized compilation option for solver only
$O/specfem2D.o: specfem2D.F90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/specfem2D.o specfem2D.F90
    
### use optimized compilation option for solver only
$O/enforce_acoustic_free_surface.o: enforce_acoustic_free_surface.f90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/enforce_acoustic_free_surface.o enforce_acoustic_free_surface.f90
    
### use optimized compilation option for solver only
$O/compute_forces_acoustic.o: compute_forces_acoustic.f90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/compute_forces_acoustic.o compute_forces_acoustic.f90
    
### use optimized compilation option for solver only
$O/compute_forces_elastic.o: compute_forces_elastic.f90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/compute_forces_elastic.o compute_forces_elastic.f90
    
### use optimized compilation option for solver only
$O/compute_gradient_attenuation.o: compute_gradient_attenuation.f90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/compute_gradient_attenuation.o compute_gradient_attenuation.f90
    
$O/compute_energy.o: compute_energy.f90 constants.h
	${F90} $(FLAGS_NOCHECK) -c -o $O/compute_energy.o compute_energy.f90
    
$O/compute_vector_field.o: compute_vector_field.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/compute_vector_field.o compute_vector_field.f90
    
$O/compute_pressure.o: compute_pressure.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/compute_pressure.o compute_pressure.f90
    
$O/compute_curl_one_element.o: compute_curl_one_element.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/compute_curl_one_element.o compute_curl_one_element.f90
    
$O/compute_Bielak_conditions.o: compute_Bielak_conditions.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/compute_Bielak_conditions.o compute_Bielak_conditions.f90
    
$O/compute_arrays_source.o: compute_arrays_source.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/compute_arrays_source.o compute_arrays_source.f90
    
$O/create_color_image.o: create_color_image.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/create_color_image.o create_color_image.f90
    
$O/spline_routines.o: spline_routines.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/spline_routines.o spline_routines.f90
    
$O/netlib_specfun_erf.o: netlib_specfun_erf.f90
	${F90} $(FLAGS_CHECK) -c -o $O/netlib_specfun_erf.o netlib_specfun_erf.f90
    
$O/define_external_model.o: define_external_model.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/define_external_model.o define_external_model.f90
    
$O/write_seismograms.o: write_seismograms.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/write_seismograms.o write_seismograms.F90
    
$O/part_unstruct.o: part_unstruct.F90 constants.h 
	${F90} $(FLAGS_CHECK) -c -o $O/part_unstruct.o part_unstruct.F90

$O/construct_acoustic_surface.o: construct_acoustic_surface.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/construct_acoustic_surface.o construct_acoustic_surface.f90

$O/assemble_MPI.o: assemble_MPI.F90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/assemble_MPI.o assemble_MPI.F90

$O/attenuation_compute_param.o: attenuation_compute_param.c
	${CC} -c -o $O/attenuation_compute_param.o attenuation_compute_param.c

$O/paco_beyond_critical.o: paco_beyond_critical.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/paco_beyond_critical.o paco_beyond_critical.f90

$O/paco_convolve_fft.o: paco_convolve_fft.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/paco_convolve_fft.o paco_convolve_fft.f90

$O/is_in_convex_quadrilateral.o: is_in_convex_quadrilateral.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/is_in_convex_quadrilateral.o is_in_convex_quadrilateral.f90

$O/get_perm_cuthill_mckee.o: get_perm_cuthill_mckee.f90 constants.h
	${F90} $(FLAGS_CHECK) -c -o $O/get_perm_cuthill_mckee.o get_perm_cuthill_mckee.f90

