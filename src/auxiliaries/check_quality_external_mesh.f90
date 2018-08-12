!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

! read an external 2D mesh file and display statistics about mesh quality;
! and create an OpenDX file showing a given range of elements or a single element

! Dimitri Komatitsch, University of Toulouse, France, January 2011.
! (adapted from the version that is available in our 3D code, SPECFEM3D)

!! DK DK
!! DK DK this routine could be improved by computing the mean in addition to min and max of ratios
!! DK DK

  program check_quality_external_mesh

  implicit none

  include "constants.h"

  integer, parameter :: NGNOD = 4     ! quadrangles

  integer :: NPOIN                    ! number of nodes
  integer :: NSPEC                    ! number of elements

  double precision, dimension(:), allocatable :: x,y,z

  integer, dimension(:,:), allocatable :: ibool

  integer :: i,ispec,iformat,ispec_min_edge_length,ispec_max_edge_length, &
             ispec_begin,ispec_end,ispec_to_output,ispec_equiangle_skewness_max

! for quality of mesh
  double precision :: equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio
  double precision :: equiangle_skewness_min,edge_aspect_ratio_min,diagonal_aspect_ratio_min
  double precision :: equiangle_skewness_max,edge_aspect_ratio_max,diagonal_aspect_ratio_max
  double precision :: skewness_AVS_DX_min,skewness_AVS_DX_max,distance_min,distance_max
  double precision :: distmin,distmax

! for histogram
  integer, parameter :: NCLASSES = 20
  integer classes_skewness(0:NCLASSES-1)
  integer :: iclass
  double precision :: current_percent,total_percent

! to export elements that have a certain skewness range to OpenDX
  integer :: ntotspecAVS_DX
  logical :: USE_OPENDX

! parameters for external mesh
  logical  :: read_external_mesh
  character(len=256)  :: mesh_file, nodes_coords_file

  integer :: NPOIN_unique_needed
  integer, dimension(:), allocatable :: ibool_reduced
  logical, dimension(:), allocatable :: mask_ibool
  integer,external :: err_occurred

! to check if any element with a negative Jacobian is found

! 2D shape functions and their derivatives at receiver
  double precision shape2D(ngnod)
  double precision dershape2D(NDIM,ngnod)

  double precision xxi,zxi,xgamma,zgamma,xelm,zelm
  double precision xi,gamma,jacobian

  integer ia

  if (NGNOD /= 4) stop 'NGNOD must be 4'

  ! ***
  ! *** read the parameter file
  ! ***

  print *,'Reading the parameter file ... '
  print *

  call open_parameter_file()

  ! read info about external mesh
  call read_value_logical_p(read_external_mesh, 'mesher.read_external_mesh')
  if (err_occurred() /= 0) stop 'error reading parameter read_external_mesh in Par_file'

  if (.not. read_external_mesh) stop 'this program is designed for read_external_mesh = .true.'

  call read_value_string_p(mesh_file, 'mesher.mesh_file')
  if (err_occurred() /= 0) stop 'error reading parameter mesh_file in Par_file'

  call read_value_string_p(nodes_coords_file, 'mesher.nodes_coords_file')
  if (err_occurred() /= 0) stop 'error reading parameter nodes_coords_file in Par_file'

  call close_parameter_file()


  print *
  print *,'1 = output elements above a certain skewness threshold in OpenDX format'
  print *,'2 = output a given element in OpenDX format'
  print *,'3 = do not output any OpenDX file'
  print *
  print *,'enter value:'
  read(5,*) iformat

  if (iformat < 1 .or. iformat > 3) stop 'exiting...'

  if (iformat == 1 .or. iformat == 2) then
    USE_OPENDX = .true.
  else
    USE_OPENDX = .false.
  endif

! read the nodes
  print *
  print *,'start reading the external node file: ',nodes_coords_file(1:len_trim(nodes_coords_file))
  open(unit=10,file=nodes_coords_file,status='unknown',action='read')

! read the header
  read(10,*) NPOIN

! read the mesh
  print *,'start reading the external mesh file: ',mesh_file(1:len_trim(mesh_file))
  open(unit=11,file=mesh_file,status='unknown',action='read')

! read the header
  read(11,*) NSPEC

  allocate(x(NPOIN))
  allocate(y(NPOIN))
  allocate(z(NPOIN))

  allocate(ibool(NGNOD,NSPEC))

  if (USE_OPENDX) then

  if (iformat == 1) then

! read range of skewness used for elements
  print *
  print *,'enter minimum skewness for OpenDX (between 0. and 0.99):'
  read(5,*) skewness_AVS_DX_min
  if (skewness_AVS_DX_min < 0.d0) skewness_AVS_DX_min = 0.d0
  if (skewness_AVS_DX_min > 0.99999d0) skewness_AVS_DX_min = 0.99999d0

!!!!!!!!  print *,'enter maximum skewness for OpenDX (between 0. and 1.):'
!!!!!!!!!!!!!  read(5,*) skewness_AVS_DX_max
  skewness_AVS_DX_max = 0.99999d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (skewness_AVS_DX_max < 0.d0) skewness_AVS_DX_max = 0.d0
  if (skewness_AVS_DX_max > 0.99999d0) skewness_AVS_DX_max = 0.99999d0

  if (skewness_AVS_DX_min > skewness_AVS_DX_max) stop 'incorrect skewness range'

  else
    print *,'enter the element number to output in OpenDX format between 1 and ',NSPEC
    read(5,*) ispec_to_output
    if (ispec_to_output < 1 .or. ispec_to_output > NSPEC) stop 'incorrect element number to output'
  endif

  endif

! read the points
  print *,'NPOIN = ',NPOIN
  do i = 1,NPOIN
    read(10,*) x(i),y(i)
! the 2D mesh is flat, therefore the third coordinate is zero
    z(i) = 0
  enddo
  close(10)

! read the elements
  print *,'NSPEC = ',NSPEC
  do i = 1,NSPEC
    read(11,*) ibool(1,i),ibool(2,i),ibool(3,i),ibool(4,i)
  enddo
  close(11)

  print *,'done reading the external files'
  print *

! ************* compute min and max of skewness and ratios ******************

  print *,'start checking if any element with a negative Jacobian is found'

  do i = 1,NSPEC

! compute jacobian matrix
  xxi = ZERO
  zxi = ZERO
  xgamma = ZERO
  zgamma = ZERO

  do ia = 1,ngnod
    xelm = x(ibool(ia,i))
    zelm = y(ibool(ia,i))

! create the 2D shape functions
    if (ia == 1) then
      xi    = -1.d0
      gamma = -1.d0
    else if (ia == 2) then
      xi    = +1.d0
      gamma = -1.d0
    else if (ia == 3) then
      xi    = +1.d0
      gamma = +1.d0
    else if (ia == 4) then
      xi    = -1.d0
      gamma = +1.d0
    else
      stop 'ia must be between 1 and NGNOD = 4'
    endif

    call define_shape_functions(shape2D,dershape2D,xi,gamma,ngnod)

    xxi = xxi + dershape2D(1,ia)*xelm
    zxi = zxi + dershape2D(1,ia)*zelm
    xgamma = xgamma + dershape2D(2,ia)*xelm
    zgamma = zgamma + dershape2D(2,ia)*zelm
  enddo

  jacobian = xxi*zgamma - xgamma*zxi

! the Jacobian is negative, this means that there is an error in the mesh
  if (jacobian <= ZERO) then
    print *,'element ',i,' has a negative Jacobian'
    stop 'negative Jacobian found'
  endif

  enddo

  print *,'OK, no element with a negative Jacobian found in the mesh'
  print *

! ************* compute min and max of skewness and ratios ******************

  print *,'start computing the minimum and maximum edge size'

! erase minimum and maximum of quality numbers
  equiangle_skewness_min = + HUGEVAL
  edge_aspect_ratio_min = + HUGEVAL
  diagonal_aspect_ratio_min = + HUGEVAL
  distance_min = + HUGEVAL

  equiangle_skewness_max = - HUGEVAL
  edge_aspect_ratio_max = - HUGEVAL
  diagonal_aspect_ratio_max = - HUGEVAL
  distance_max = - HUGEVAL

  ispec_min_edge_length = -1
  ispec_max_edge_length = -1

! loop on all the elements
  do ispec = 1,NSPEC

    if (mod(ispec,100000) == 0) print *,'processed ',ispec,' elements out of ',NSPEC

      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax)

! store element number in which the edge of minimum or maximum length is located
    if (distmin < distance_min) ispec_min_edge_length = ispec
    if (distmax > distance_max) ispec_max_edge_length = ispec

! compute minimum and maximum of quality numbers
    equiangle_skewness_min = min(equiangle_skewness_min,equiangle_skewness)
    edge_aspect_ratio_min = min(edge_aspect_ratio_min,edge_aspect_ratio)
    diagonal_aspect_ratio_min = min(diagonal_aspect_ratio_min,diagonal_aspect_ratio)
    distance_min = min(distance_min,distmin)

    if (equiangle_skewness > equiangle_skewness_max) ispec_equiangle_skewness_max = ispec
    equiangle_skewness_max = max(equiangle_skewness_max,equiangle_skewness)
    edge_aspect_ratio_max = max(edge_aspect_ratio_max,edge_aspect_ratio)
    diagonal_aspect_ratio_max = max(diagonal_aspect_ratio_max,diagonal_aspect_ratio)
    distance_max = max(distance_max,distmax)

  enddo
  print *,'done processing ',NSPEC,' elements out of ',NSPEC

  print *
  print *,'------------'
  print *,'mesh quality parameter definitions:'
  print *
  print *,'equiangle skewness: 0. perfect,  1. bad'
  print *,'skewness max deviation angle: 0. perfect,  90. bad'
  print *,'edge aspect ratio: 1. perfect,  above 1. gives stretching factor'
  print *,'diagonal aspect ratio: 1. perfect,  above 1. gives stretching factor'
  print *,'------------'

  print *
  print *,'minimum length of an edge in the whole mesh (m) = ',distance_min,' in element ',ispec_min_edge_length
  print *
  print *,'maximum length of an edge in the whole mesh (m) = ',distance_max,' in element ',ispec_max_edge_length
  print *
  print *,'max equiangle skewness = ',equiangle_skewness_max
  print *,'in element ',ispec_equiangle_skewness_max
! print *,'min equiangle skewness = ',equiangle_skewness_min
  print *
  print *,'max deviation angle from a right angle (90 degrees) is therefore = ',90.*equiangle_skewness_max
  print *
  print *,'worst angle in the mesh is therefore either ',90.*(1. - equiangle_skewness_max)
  print *,'or ',180. - 90.*(1. - equiangle_skewness_max),' degrees'
  print *
  print *,'max edge aspect ratio = ',edge_aspect_ratio_max
! print *,'min edge aspect ratio = ',edge_aspect_ratio_min
  print *
  print *,'max diagonal aspect ratio = ',diagonal_aspect_ratio_max
! print *,'min diagonal aspect ratio = ',diagonal_aspect_ratio_min
  print *

! create statistics about mesh quality
  print *,'creating histogram and statistics of mesh quality'

! erase histogram of skewness
  classes_skewness(:) = 0

! loop on all the elements
  do ispec = 1,NSPEC

      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax)

! store skewness in histogram
    iclass = int(equiangle_skewness * dble(NCLASSES))
    if (iclass < 0) iclass = 0
    if (iclass > NCLASSES-1) iclass = NCLASSES-1
    classes_skewness(iclass) = classes_skewness(iclass) + 1

  enddo

! create histogram of skewness and save in Gnuplot file
  print *
  print *,'histogram of skewness (0. good - 1. bad):'
  print *
  total_percent = 0.
  open(unit=14,file='OUTPUT_FILES/mesh_quality_histogram.txt',status='unknown')
  do iclass = 0,NCLASSES-1
    current_percent = 100.*dble(classes_skewness(iclass))/dble(NSPEC)
    total_percent = total_percent + current_percent
    print *,real(iclass/dble(NCLASSES)),' - ',real((iclass+1)/dble(NCLASSES)),classes_skewness(iclass),' ', &
       sngl(current_percent),' %'
    write(14,*) 0.5*(real(iclass/dble(NCLASSES)) + real((iclass+1)/dble(NCLASSES))),' ',sngl(current_percent)
  enddo
  close(14)

! create script for Gnuplot histogram file
  open(unit=14,file='OUTPUT_FILES/plot_mesh_quality_histogram.gnu',status='unknown')
  write(14,*) 'set term wxt'
  write(14,*) '#set term gif'
  write(14,*) '#set output "mesh_quality_histogram.gif"'
  write(14,*)
  write(14,*) 'set xrange [0:1]'
  write(14,*) 'set xtics 0,0.1,1'
  write(14,*) 'set boxwidth ',1./real(NCLASSES)
  write(14,*) 'set xlabel "Skewness range"'
  write(14,*) 'set ylabel "Percentage of elements (%)"'
  write(14,*) 'set loadpath "./OUTPUT_FILES"'
  write(14,*) 'plot "mesh_quality_histogram.txt" with boxes'
  write(14,*) 'pause -1 "hit any key..."'
  close(14)

  print *
  print *,'total number of elements = ',NSPEC
  print *

! display warning if maximum skewness is too high
  if (equiangle_skewness_max >= 0.75d0) then
    print *
    print *,'*********************************************'
    print *,'*********************************************'
    print *,' WARNING, mesh is bad (max skewness >= 0.75)'
    print *,'*********************************************'
    print *,'*********************************************'
    print *
  endif

  if (total_percent < 99.9d0 .or. total_percent > 100.1d0) then
    print *,'total percentage = ',total_percent,' %'
    stop 'total percentage should be 100%'
  endif

! ************* create OpenDX file with elements in a certain range of skewness

  if (USE_OPENDX) then

  print *
  if (iformat == 1) then
    print *,'creating OpenDX file with subset of elements in skewness range'
    print *,'between ',skewness_AVS_DX_min,' and ',skewness_AVS_DX_max
  else
    print *,'creating OpenDX file with element #',ispec_to_output
  endif
  print *

! ************* count number of elements in skewness range *************

! erase number of elements belonging to skewness range for AVS_DX
  ntotspecAVS_DX = 0

! loop on all the elements
  if (iformat == 1) then

  do ispec = 1,NSPEC

      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax)

! check if element belongs to requested skewness range
    if (equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max) &
        ntotspecAVS_DX = ntotspecAVS_DX + 1

  enddo

  else
! outputing a single element
    ntotspecAVS_DX = 1
  endif

  if (ntotspecAVS_DX == 0) then
    stop 'no elements in skewness range, no file created'
  else if (iformat == 1) then
    print *
    print *,'there are ',ntotspecAVS_DX,' elements in AVS or DX skewness range ',skewness_AVS_DX_min,skewness_AVS_DX_max
    print *
  endif

  open(unit=11,file='DX_mesh_quality.dx',status='unknown')

! generate the subset of points that are needed

! count the number of unique points
  NPOIN_unique_needed = 0
  allocate(mask_ibool(NPOIN))
  mask_ibool(:) = .false.

! loop on all the elements
  if (iformat == 1) then
    ispec_begin = 1
    ispec_end = NSPEC
  else
    ispec_begin = ispec_to_output
    ispec_end = ispec_to_output
  endif

  do ispec = ispec_begin,ispec_end

      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax)

! check if element needs to be output
    if (iformat == 2 .or. (iformat == 1 .and. &
       equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max)) then
! create point for first corner of the element
       if (.not. mask_ibool(ibool(1,ispec))) then
         mask_ibool(ibool(1,ispec)) = .true.
         NPOIN_unique_needed = NPOIN_unique_needed + 1
       endif

! create point for second corner of the element
       if (.not. mask_ibool(ibool(2,ispec))) then
         mask_ibool(ibool(2,ispec)) = .true.
         NPOIN_unique_needed = NPOIN_unique_needed + 1
       endif

! create point for third corner of the element
       if (.not. mask_ibool(ibool(3,ispec))) then
         mask_ibool(ibool(3,ispec)) = .true.
         NPOIN_unique_needed = NPOIN_unique_needed + 1
       endif

! create point for fourth corner of the element
       if (.not. mask_ibool(ibool(4,ispec))) then
         mask_ibool(ibool(4,ispec)) = .true.
         NPOIN_unique_needed = NPOIN_unique_needed + 1
       endif

    endif

  enddo


! ************* generate points ******************

! write OpenDX header
  write(11,*) 'object 1 class array type float rank 1 shape 3 items ',NPOIN_unique_needed,' data follows'

  allocate(ibool_reduced(NPOIN))

! count the number of unique points
  NPOIN_unique_needed = 0
  mask_ibool(:) = .false.

! loop on all the elements
  if (iformat == 1) then
    ispec_begin = 1
    ispec_end = NSPEC
  else
    ispec_begin = ispec_to_output
    ispec_end = ispec_to_output
  endif

  do ispec = ispec_begin,ispec_end

      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax)

! check if element needs to be output
    if (iformat == 2 .or. (iformat == 1 .and. &
       equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max)) then
! create point for first corner of the element
       if (.not. mask_ibool(ibool(1,ispec))) then
         mask_ibool(ibool(1,ispec)) = .true.
         ibool_reduced(ibool(1,ispec)) = NPOIN_unique_needed
         write(11,*) sngl(x(ibool(1,ispec))),sngl(y(ibool(1,ispec))),sngl(z(ibool(1,ispec)))
         NPOIN_unique_needed = NPOIN_unique_needed + 1
       endif

! create point for second corner of the element
       if (.not. mask_ibool(ibool(2,ispec))) then
         mask_ibool(ibool(2,ispec)) = .true.
         ibool_reduced(ibool(2,ispec)) = NPOIN_unique_needed
         write(11,*) sngl(x(ibool(2,ispec))),sngl(y(ibool(2,ispec))),sngl(z(ibool(2,ispec)))
         NPOIN_unique_needed = NPOIN_unique_needed + 1
       endif

! create point for third corner of the element
       if (.not. mask_ibool(ibool(3,ispec))) then
         mask_ibool(ibool(3,ispec)) = .true.
         ibool_reduced(ibool(3,ispec)) = NPOIN_unique_needed
         write(11,*) sngl(x(ibool(3,ispec))),sngl(y(ibool(3,ispec))),sngl(z(ibool(3,ispec)))
         NPOIN_unique_needed = NPOIN_unique_needed + 1
       endif

! create point for fourth corner of the element
       if (.not. mask_ibool(ibool(4,ispec))) then
         mask_ibool(ibool(4,ispec)) = .true.
         ibool_reduced(ibool(4,ispec)) = NPOIN_unique_needed
         write(11,*) sngl(x(ibool(4,ispec))),sngl(y(ibool(4,ispec))),sngl(z(ibool(4,ispec)))
         NPOIN_unique_needed = NPOIN_unique_needed + 1
       endif

    endif

  enddo

  deallocate(mask_ibool)

! ************* generate elements ******************

  write(11,*) 'object 2 class array type int rank 1 shape ',NGNOD,' items ',ntotspecAVS_DX,' data follows'

! loop on all the elements
  if (iformat == 1) then
    ispec_begin = 1
    ispec_end = NSPEC
  else
    ispec_begin = ispec_to_output
    ispec_end = ispec_to_output
  endif

  do ispec = ispec_begin,ispec_end

      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax)

! check if element needs to be output
    if (iformat == 2 .or. (iformat == 1 .and. &
       equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max)) then
! point order in OpenDX in 2D is 1,4,2,3 *not* 1,2,3,4 as in AVS
! point order in OpenDX in 3D is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
! in the case of OpenDX, node numbers start at zero
      write(11,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") &
            ibool_reduced(ibool(1,ispec)), ibool_reduced(ibool(4,ispec)), &
            ibool_reduced(ibool(2,ispec)), ibool_reduced(ibool(3,ispec))
      if (iformat == 1) print *,'element ',ispec,' belongs to the range and has skewness = ',sngl(equiangle_skewness)
    endif

  enddo

! ************* generate element data values ******************

! output OpenDX header for data
  write(11,*) 'attribute "element type" string "quads"'
  write(11,*) 'attribute "ref" string "positions"'
  write(11,*) 'object 3 class array type float rank 0 items ',ntotspecAVS_DX,' data follows'

! loop on all the elements
  do ispec = ispec_begin,ispec_end

      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax)

! check if element needs to be output
    if (iformat == 2 .or. (iformat == 1 .and. &
       equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max)) &
    write(11,*) sngl(equiangle_skewness)

  enddo

! define OpenDX field
  write(11,*) 'attribute "dep" string "connections"'
  write(11,*) 'object "irregular positions irregular connections" class field'
  write(11,*) 'component "positions" value 1'
  write(11,*) 'component "connections" value 2'
  write(11,*) 'component "data" value 3'
  write(11,*) 'end'

! close OpenDX file
  close(11)

  endif

  end program check_quality_external_mesh

!
!=====================================================================
!

! create mesh quality data for a given 2D spectral element

  subroutine create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax)

  implicit none

  include "constants.h"

  integer :: icorner,ispec,NSPEC,NPOIN,NGNOD,i

  double precision, dimension(NPOIN) :: x,y,z

  integer, dimension(NGNOD,NSPEC) :: ibool

  double precision, dimension(NGNOD) :: xelm,yelm,zelm

  double precision vectorA_x,vectorA_y,vectorA_z
  double precision vectorB_x,vectorB_y,vectorB_z
  double precision norm_A,norm_B,angle_vectors
  double precision distmin,distmax,dist,dist1,dist2
  double precision equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio

! topology of faces of cube for skewness
! only one face in 2D
  integer faces_topo(6)

! store the corners of this element for the skewness routine
  do i = 1,NGNOD
    xelm(i) = x(ibool(i,ispec))
    yelm(i) = y(ibool(i,ispec))
    zelm(i) = z(ibool(i,ispec))
  enddo

! define topology of faces of cube for skewness

! only one face in 2D
  faces_topo(1) = 1
  faces_topo(2) = 2
  faces_topo(3) = 3
  faces_topo(4) = 4

! define wraparound for angles for skewness calculation
  faces_topo(5) = faces_topo(1)
  faces_topo(6) = faces_topo(2)

! compute equiangle skewness (as defined in Fluent/Gambit manual)
! and compute edge aspect ratio using the corners of the element
     distmin = + HUGEVAL
     distmax = - HUGEVAL
     equiangle_skewness = - HUGEVAL

     do icorner = 1,4

! first vector of angle
       vectorA_x = xelm(faces_topo(icorner)) - xelm(faces_topo(icorner+1))
       vectorA_y = yelm(faces_topo(icorner)) - yelm(faces_topo(icorner+1))
       vectorA_z = zelm(faces_topo(icorner)) - zelm(faces_topo(icorner+1))

! second vector of angle
       vectorB_x = xelm(faces_topo(icorner+2)) - xelm(faces_topo(icorner+1))
       vectorB_y = yelm(faces_topo(icorner+2)) - yelm(faces_topo(icorner+1))
       vectorB_z = zelm(faces_topo(icorner+2)) - zelm(faces_topo(icorner+1))

! norm of vectors A and B
       norm_A = sqrt(vectorA_x**2 + vectorA_y**2 + vectorA_z**2)
       norm_B = sqrt(vectorB_x**2 + vectorB_y**2 + vectorB_z**2)

! angle formed by the two vectors
       angle_vectors = dacos((vectorA_x*vectorB_x + vectorA_y*vectorB_y + vectorA_z*vectorB_z) / (norm_A * norm_B))

! compute equiangle skewness
       equiangle_skewness = max(equiangle_skewness,dabs(2.d0 * angle_vectors - PI) / PI)

! compute min and max size of an edge
       dist = sqrt(vectorA_x**2 + vectorA_y**2 + vectorA_z**2)

       distmin = min(distmin,dist)
       distmax = max(distmax,dist)

     enddo

! compute edge aspect ratio
   edge_aspect_ratio = distmax / distmin

! compute diagonal aspect ratio
   dist1 = sqrt((xelm(1) - xelm(3))**2 + (yelm(1) - yelm(3))**2 + (zelm(1) - zelm(3))**2)
   dist2 = sqrt((xelm(2) - xelm(4))**2 + (yelm(2) - yelm(4))**2 + (zelm(2) - zelm(4))**2)
   diagonal_aspect_ratio = max(dist1,dist2) / min(dist1,dist2)

  end subroutine create_mesh_quality_data_2D

