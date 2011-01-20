
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

! read a 2D Gmsh mesh file and display statistics about mesh quality;
! and create an OpenDX file showing a given range of elements or a single element

! Dimitri Komatitsch, University of Toulouse, France, January 2011.
! (adapted from the version that is available in our 3D code, SPECFEM3D)

!! DK DK
!! DK DK this routine could be improved by computing the mean in addition to min and max of ratios
!! DK DK

  program check_mesh_quality_Gmsh

  implicit none

  include "constants.h"

  integer, parameter :: NGNOD = 4                       ! quadrangles

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
  integer, parameter :: NCLASS = 20
  integer classes_skewness(0:NCLASS-1)
  integer :: iclass
  double precision :: current_percent,total_percent

! to export elements that have a certain skewness range to OpenDX
  integer :: ntotspecAVS_DX
  logical :: USE_OPENDX

  character(len=50) interfacesfile,title

  ! flag to save the last frame for kernels calculation purpose and type of simulation
  logical :: SAVE_FORWARD
  integer :: SIMULATION_TYPE

  ! parameters for external mesh
  logical  :: read_external_mesh
  character(len=256)  :: mesh_file, nodes_coords_file

  ! ignore variable name field (junk) at the beginning of each input line
  logical, parameter :: IGNORE_JUNK = .true.

  if(NGNOD /= 4) stop 'NGNOD must be 4'

  ! ***
  ! *** read the parameter file
  ! ***

  print *,'Reading the parameter file ... '
  print *

  open(unit=IIN,file='DATA/Par_file',status='old')

  ! read and ignore file names and path for output
  call read_value_string(IIN,IGNORE_JUNK,title)
  call read_value_string(IIN,IGNORE_JUNK,interfacesfile)

  ! read and ignore type of simulation
  call read_value_integer(IIN,IGNORE_JUNK,SIMULATION_TYPE)
  call read_value_logical(IIN,IGNORE_JUNK,SAVE_FORWARD)

  ! read info about external mesh
  call read_value_logical(IIN,IGNORE_JUNK,read_external_mesh)
  call read_value_string(IIN,IGNORE_JUNK,mesh_file)
  call read_value_string(IIN,IGNORE_JUNK,nodes_coords_file)

  print *
  print *,'1 = output elements above a certain skewness threshold in OpenDX format'
  print *,'2 = output a given element in OpenDX format'
  print *,'3 = do not output any OpenDX file'
  print *
  print *,'enter value:'
  read(5,*) iformat

  if(iformat < 1 .or. iformat > 3) stop 'exiting...'

  if(iformat == 1 .or. iformat == 2) then
    USE_OPENDX = .true.
  else
    USE_OPENDX = .false.
  endif

! read the nodes
  print *
  print *,'start reading the Gmsh node file: ',nodes_coords_file(1:len_trim(nodes_coords_file))
  open(unit=10,file=nodes_coords_file,status='unknown',action='read')

! read the header
  read(10,*) NPOIN

! read the mesh
  print *,'start reading the Gmsh mesh file: ',mesh_file(1:len_trim(mesh_file))
  open(unit=11,file=mesh_file,status='unknown',action='read')

! read the header
  read(11,*) NSPEC

  allocate(x(NPOIN))
  allocate(y(NPOIN))
  allocate(z(NPOIN))

  allocate(ibool(NGNOD,NSPEC))

  if(USE_OPENDX) then

  if(iformat == 1) then

! read range of skewness used for elements
  print *,'enter minimum skewness for OpenDX (between 0. and 0.99):'
  read(5,*) skewness_AVS_DX_min
  if(skewness_AVS_DX_min < 0.d0) skewness_AVS_DX_min = 0.d0
  if(skewness_AVS_DX_min > 0.99999d0) skewness_AVS_DX_min = 0.99999d0

!!!!!!!!  print *,'enter maximum skewness for OpenDX (between 0. and 1.):'
!!!!!!!!!!!!!  read(5,*) skewness_AVS_DX_max
  skewness_AVS_DX_max = 0.99999d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(skewness_AVS_DX_max < 0.d0) skewness_AVS_DX_max = 0.d0
  if(skewness_AVS_DX_max > 0.99999d0) skewness_AVS_DX_max = 0.99999d0

  if(skewness_AVS_DX_min > skewness_AVS_DX_max) stop 'incorrect skewness range'

  else
    print *,'enter the element number to output in OpenDX format between 1 and ',NSPEC
    read(5,*) ispec_to_output
    if(ispec_to_output < 1 .or. ispec_to_output > NSPEC) stop 'incorrect element number to output'
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

  print *,'done reading the Gmsh files'
  print *

  print *,'start computing the minimum and maximum edge size'

! ************* compute min and max of skewness and ratios ******************

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

    if(mod(ispec,100000) == 0) print *,'processed ',ispec,' elements out of ',NSPEC

      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax)

! store element number in which the edge of minimum or maximum length is located
    if(distmin < distance_min) ispec_min_edge_length = ispec
    if(distmax > distance_max) ispec_max_edge_length = ispec

! compute minimum and maximum of quality numbers
    equiangle_skewness_min = min(equiangle_skewness_min,equiangle_skewness)
    edge_aspect_ratio_min = min(edge_aspect_ratio_min,edge_aspect_ratio)
    diagonal_aspect_ratio_min = min(diagonal_aspect_ratio_min,diagonal_aspect_ratio)
    distance_min = min(distance_min,distmin)

    if(equiangle_skewness > equiangle_skewness_max) ispec_equiangle_skewness_max = ispec
    equiangle_skewness_max = max(equiangle_skewness_max,equiangle_skewness)
    edge_aspect_ratio_max = max(edge_aspect_ratio_max,edge_aspect_ratio)
    diagonal_aspect_ratio_max = max(diagonal_aspect_ratio_max,diagonal_aspect_ratio)
    distance_max = max(distance_max,distmax)

  enddo
  print *,'done processing ',NSPEC,' elements out of ',NSPEC

  print *
  print *,'------------'
  print *,'mesh quality parameter definitions'
  print *
  print *,'equiangle skewness: 0. perfect  1. bad'
  print *,'skewness max deviation angle: 0. perfect  90. bad'
  print *,'edge aspect ratio: 1. perfect  above 1. gives stretching factor'
  print *,'diagonal aspect ratio: 1. perfect  above 1. gives stretching factor'
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
  print *,'worst angle in the mesh is therefore ',90.*(1. - equiangle_skewness_max)
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
    iclass = int(equiangle_skewness * dble(NCLASS))
    if(iclass < 0) iclass = 0
    if(iclass > NCLASS-1) iclass = NCLASS-1
    classes_skewness(iclass) = classes_skewness(iclass) + 1

  enddo

! create histogram of skewness and save in Gnuplot file
  print *
  print *,'histogram of skewness (0. good - 1. bad):'
  print *
  total_percent = 0.
  open(unit=14,file='mesh_quality_histogram.txt',status='unknown')
  do iclass = 0,NCLASS-1
    current_percent = 100.*dble(classes_skewness(iclass))/dble(NSPEC)
    total_percent = total_percent + current_percent
    print *,real(iclass/dble(NCLASS)),' - ',real((iclass+1)/dble(NCLASS)),classes_skewness(iclass),' ',sngl(current_percent),' %'
    write(14,*) 0.5*(real(iclass/dble(NCLASS)) + real((iclass+1)/dble(NCLASS))),' ',sngl(current_percent)
  enddo
  close(14)

! create script for Gnuplot histogram file
  open(unit=14,file='plot_mesh_quality_histogram.gnu',status='unknown')
  write(14,*) 'set term x11'
  write(14,*) '#set term gif'
  write(14,*) '#set output "mesh_quality_histogram.gif"'
  write(14,*)
  write(14,*) 'set xrange [0:1]'
  write(14,*) 'set xtics 0,0.1,1'
  write(14,*) 'set boxwidth ',1./real(NCLASS)
  write(14,*) 'set xlabel "Skewness range"'
  write(14,*) 'set ylabel "Percentage of elements (%)"'
  write(14,*) 'plot "mesh_quality_histogram.txt" with boxes'
  write(14,*) 'pause -1 "hit any key..."'
  close(14)

  print *
  print *,'total number of elements = ',NSPEC
  print *

! display warning if maximum skewness is too high
  if(equiangle_skewness_max >= 0.75d0) then
    print *
    print *,'*********************************************'
    print *,'*********************************************'
    print *,' WARNING, mesh is bad (max skewness >= 0.75)'
    print *,'*********************************************'
    print *,'*********************************************'
    print *
  endif

  if(total_percent < 99.9d0 .or. total_percent > 100.1d0) then
    print *,'total percentage = ',total_percent,' %'
    stop 'total percentage should be 100%'
  endif

! ************* create OpenDX file with elements in a certain range of skewness

  if(USE_OPENDX) then

  print *
  if(iformat == 1) then
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
  if(iformat == 1) then

  do ispec = 1,NSPEC

      call create_mesh_quality_data_2D(x,y,z,ibool,ispec,NSPEC,NPOIN,NGNOD, &
               equiangle_skewness,edge_aspect_ratio,diagonal_aspect_ratio,distmin,distmax)

! check if element belongs to requested skewness range
    if(equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max) &
        ntotspecAVS_DX = ntotspecAVS_DX + 1

  enddo

  else
! outputing a single element
    ntotspecAVS_DX = 1
  endif

  if(ntotspecAVS_DX == 0) then
    stop 'no elements in skewness range, no file created'
  else if(iformat == 1) then
    print *
    print *,'there are ',ntotspecAVS_DX,' elements in AVS or DX skewness range ',skewness_AVS_DX_min,skewness_AVS_DX_max
    print *
  endif

  open(unit=11,file='DX_mesh_quality.dx',status='unknown')

! ************* generate points ******************

! write OpenDX header
  write(11,*) 'object 1 class array type float rank 1 shape 3 items ',NPOIN,' data follows'

! write all the points
  do i = 1,NPOIN
    write(11,*) sngl(x(i)),sngl(y(i)),sngl(z(i))
  enddo

! ************* generate elements ******************

  write(11,*) 'object 2 class array type int rank 1 shape ',NGNOD,' items ',ntotspecAVS_DX,' data follows'

! loop on all the elements
  if(iformat == 1) then
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
    if(iformat == 2 .or. (iformat == 1 .and. &
       equiangle_skewness >= skewness_AVS_DX_min .and. equiangle_skewness <= skewness_AVS_DX_max)) then
! point order in OpenDX in 2D is 1,4,2,3 *not* 1,2,3,4 as in AVS
! point order in OpenDX in 3D is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
! in the case of OpenDX, node numbers start at zero
      write(11,"(i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9,1x,i9)") &
            ibool(1,ispec)-1, ibool(4,ispec)-1, ibool(2,ispec)-1, ibool(3,ispec)-1
      if(iformat == 1) print *,'element ',ispec,' belongs to the range and has skewness = ',sngl(equiangle_skewness)
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
    if(iformat == 2 .or. (iformat == 1 .and. &
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

  end program check_mesh_quality_Gmsh

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

