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

  program add_CPML_layers_to_a_given_mesh

! Dimitri Komatitsch, CNRS Marseille, France, March 2016 and January 2017

! add PML layers around an existing mesh (i.e. create new elements and new points)

  implicit none

! this code is for QUAD4 only for now
  integer, parameter :: NGNOD = 4

! we are in 2D
  integer, parameter :: NDIM = 2

! number of GLL points in each direction, to check for negative Jacobians
  integer, parameter :: NGLLX = 5,NGLLZ = NGLLX

! number of PML and non-PML layers to add on each side of the mesh
  integer :: NUMBER_OF_PML_LAYERS_TO_ADD,NUMBER_OF_TRANSITION_LAYERS_TO_ADD,TOTAL_NUMBER_OF_LAYERS_TO_ADD

! size of each PML element to add when they are added on the Xmin and Xmax faces, Zmin and/or Zmax faces
  double precision :: SIZE_OF_X_ELEMENT_TO_ADD,SIZE_OF_Z_ELEMENT_TO_ADD
  double precision :: SIZE_OF_XMIN_ELEMENT_TO_ADD,SIZE_OF_ZMIN_ELEMENT_TO_ADD
  double precision :: SIZE_OF_XMAX_ELEMENT_TO_ADD,SIZE_OF_ZMAX_ELEMENT_TO_ADD

  integer :: nspec,npoin,npoin_new_max,npoin_new_real,nspec_new,count_elem_faces_to_extend,iextend
  integer :: factor_x,factor_z
  integer :: i,k,ispec,ipoin,iloop_on_X_Z_faces,iloop_on_min_face_then_max_face
  integer :: i1,i2,i3,i4,elem_counter,ia,iflag,icompute_size
  integer :: p1,p2

  double precision :: x_value_to_create,z_value_to_create
  logical :: point_already_exists

  double precision, dimension(:), allocatable, target :: x,z
  double precision, dimension(:), allocatable :: x_new,z_new
  double precision, dimension(:), pointer :: coord_to_use1,coord_to_use2

  integer, dimension(:), allocatable :: imaterial,imaterial_new

  integer, dimension(:,:), allocatable :: ibool,ibool_new

! Gauss-Lobatto-Legendre points of integration, to check for negative Jacobians
  double precision xigll(NGLLX)
  double precision zigll(NGLLZ)

! 2D shape function derivatives, to check for negative Jacobians
  double precision shape2D(NGNOD,NGLLX,NGLLZ)
  double precision dershape2D(NDIM,NGNOD,NGLLX,NGLLZ)

  double precision, dimension(NGNOD) :: xelm,zelm

  double precision :: xmin,xmax,zmin,zmax,limit,xsize,zsize
  double precision :: value_min,value_max,value_size,sum_of_distances,mean_distance,distance_squared,very_small_distance_squared

  logical :: ADD_ON_THE_XMIN_SURFACE,ADD_ON_THE_XMAX_SURFACE
  logical :: ADD_ON_THE_ZMIN_SURFACE,ADD_ON_THE_ZMAX_SURFACE

  logical :: need_to_extend_this_element,found_a_negative_Jacobian1,found_a_negative_Jacobian2

  double precision, parameter :: SMALL_RELATIVE_VALUE = 1.d-5

  print *
  print *,'IMPORTANT: it is your responsibility to make sure that the input mesh'
  print *,'that this code will read in SPECFEM2D format from files "nodes_coords_file" and "mesh_file"'
  print *,'has flat outer edges aligned with the coordinate grid axes (X and/or Z),'
  print *,'so that this code can external CPML elements to it.'
  print *,'This code does NOT check that (because it cannot, in any easy way).'
  print *,'The mesh does not need to be structured nor regular, any non-structured'
  print *,'mesh is fine as long as it has flat outer faces, parallel to the axes.'
  print *

  print *,'enter the number of PML layers to add on each side of the mesh (usually 3, can also be 4):'
  read(*,*) NUMBER_OF_PML_LAYERS_TO_ADD
  if (NUMBER_OF_PML_LAYERS_TO_ADD < 1) stop 'NUMBER_OF_PML_LAYERS_TO_ADD must be >= 1'
  print *

  print *,'enter the number of non-PML transition layers to add between each side of the mesh and the new PMLs &
                &that will be created (usually 0, 1 or 2) (if you have no idea, just enter 0):'
  read(*,*) NUMBER_OF_TRANSITION_LAYERS_TO_ADD
  if (NUMBER_OF_TRANSITION_LAYERS_TO_ADD < 0) stop 'NUMBER_OF_TRANSITION_LAYERS_TO_ADD must be >= 0'
  print *

  TOTAL_NUMBER_OF_LAYERS_TO_ADD = NUMBER_OF_TRANSITION_LAYERS_TO_ADD + NUMBER_OF_PML_LAYERS_TO_ADD

  ADD_ON_THE_XMIN_SURFACE = .true.
  ADD_ON_THE_XMAX_SURFACE = .true.
  ADD_ON_THE_ZMIN_SURFACE = .true.
  ADD_ON_THE_ZMAX_SURFACE = .true.

  print *,'1 = use a free surface at the top of the mesh i.e. do not add a CPML layer at the top (most classical option)'
  print *,'2 = add a CPML absorbing layer at the top of the mesh (less classical option)'
  print *,'3 = exit'
  read(*,*) iflag
  if (iflag /= 1 .and. iflag /= 2) stop 'exiting...'
  if (iflag == 1) then
    ADD_ON_THE_ZMAX_SURFACE = .false.
  else
    ADD_ON_THE_ZMAX_SURFACE = .true.
  endif
  print *

  print *,'1 = compute the size of the PML elements to add automatically using the average size of edge elements'
  print *,'2 = enter the size of the PML elements to add manually'
  print *,'3 = exit'
  read(*,*) icompute_size
  if (icompute_size /= 1 .and. icompute_size /= 2) stop 'exiting...'

  if (icompute_size == 2) then

  print *,'enter the X size (in meters) of each CPML element to add on the Xmin face (enter -1 to turn that PML layer off):'
  read(*,*) SIZE_OF_XMIN_ELEMENT_TO_ADD
  if (SIZE_OF_XMIN_ELEMENT_TO_ADD <= 0.d0) then
    SIZE_OF_XMIN_ELEMENT_TO_ADD = 0.d0
    ADD_ON_THE_XMIN_SURFACE = .false.
  endif
  print *

  print *,'enter the X size (in meters) of each CPML element to add on the Xmax face (enter -1 to turn that PML layer off):'
  read(*,*) SIZE_OF_XMAX_ELEMENT_TO_ADD
  if (SIZE_OF_XMAX_ELEMENT_TO_ADD <= 0.d0) then
    SIZE_OF_XMAX_ELEMENT_TO_ADD = 0.d0
    ADD_ON_THE_XMAX_SURFACE = .false.
  endif
  print *

  print *,'enter the Z size (in meters) of each CPML element to add on the Zmin faces (enter -1 to turn that PML layer off):'
  read(*,*) SIZE_OF_ZMIN_ELEMENT_TO_ADD
  if (SIZE_OF_ZMIN_ELEMENT_TO_ADD <= 0.d0) then
    SIZE_OF_ZMIN_ELEMENT_TO_ADD = 0.d0
    ADD_ON_THE_ZMIN_SURFACE = .false.
  endif
  print *

  if (ADD_ON_THE_ZMAX_SURFACE) then
    print *,'enter the Z size (in meters) of each CPML element to add on the Zmax faces (enter -1 to turn that PML layer off):'
    read(*,*) SIZE_OF_ZMAX_ELEMENT_TO_ADD
    if (SIZE_OF_ZMAX_ELEMENT_TO_ADD <= 0.d0) then
      SIZE_OF_ZMAX_ELEMENT_TO_ADD = 0.d0
      ADD_ON_THE_ZMAX_SURFACE = .false.
    endif
    print *
  endif

  endif

! check that we need to add at least one PML, otherwise this code is useless
  if (.not. ADD_ON_THE_XMIN_SURFACE .and. .not. ADD_ON_THE_XMAX_SURFACE &
      .and. .not. ADD_ON_THE_ZMIN_SURFACE .and. .not. ADD_ON_THE_ZMAX_SURFACE) &
    stop 'Error: the purpose of this code is to add at least one PML, but you have added none'

! hardwire GLL point location values to avoid having to link with a long library to compute them
  xigll(:) = (/ -1.d0 , -0.654653670707977d0 , 0.d0 , 0.654653670707977d0 , 1.d0 /)
  zigll(:) = xigll(:)

! shape arrays
  do k = 1,NGLLZ
    do i = 1,NGLLX
      call define_shape_functions(shape2D(:,i,k),dershape2D(:,:,i,k),xigll(i),zigll(k),NGNOD)
    enddo
  enddo

! open SPECFEM2D mesh file to read the points
  open(unit=23,file='nodes_coords_file',status='old',action='read')
  read(23,*) npoin
  allocate(x(npoin))
  allocate(z(npoin))
  do ipoin = 1,npoin
    read(23,*) x(ipoin),z(ipoin)
  enddo
  close(23)

! ************* read mesh elements *************

! open SPECFEM2D topology file to read the mesh elements
  open(unit=23,file='mesh_file',status='old',action='read')
  read(23,*) nspec

  allocate(ibool(NGNOD,nspec))

! loop on the whole mesh
  do ispec = 1,nspec
    read(23,*) ibool(1,ispec),ibool(2,ispec),ibool(3,ispec),ibool(4,ispec)
  enddo

  close(23)

  print *,'Total number of elements in the mesh read = ',nspec

! read the materials file
  allocate(imaterial(nspec))
  open(unit=23,file='materials_file',status='old',action='read')
! loop on the whole mesh
  do ispec = 1,nspec
    read(23,*) imaterial(ispec)
  enddo
  close(23)

! we need to extend/extrude the existing mesh by adding CPML elements
! along the X faces, then along the Z faces.

! loop on the three sets of faces to first add CPML elements along X, then along Z
  do iloop_on_X_Z_faces = 1,NDIM

    do iloop_on_min_face_then_max_face = 1,2  ! 1 is min face and 2 is max face (Xmin then Xmax, or Zmin then Zmax)

! do not add a CPML layer on a given surface if the user asked not to
      if (iloop_on_X_Z_faces == 1 .and. iloop_on_min_face_then_max_face == 1 .and. .not. ADD_ON_THE_XMIN_SURFACE) cycle
      if (iloop_on_X_Z_faces == 1 .and. iloop_on_min_face_then_max_face == 2 .and. .not. ADD_ON_THE_XMAX_SURFACE) cycle
      if (iloop_on_X_Z_faces == 2 .and. iloop_on_min_face_then_max_face == 1 .and. .not. ADD_ON_THE_ZMIN_SURFACE) cycle
      if (iloop_on_X_Z_faces == 2 .and. iloop_on_min_face_then_max_face == 2 .and. .not. ADD_ON_THE_ZMAX_SURFACE) cycle

    print *
    print *,'********************************************************************'
    if (iloop_on_X_Z_faces == 1) then
      print *,'adding CPML elements along one of the two X faces of the existing mesh'
    else if (iloop_on_X_Z_faces == 2) then
      print *,'adding CPML elements along one of the two Z faces of the existing mesh'
    else
      stop 'wrong index in loop on faces'
    endif
    print *,'********************************************************************'
    print *

! compute the min and max values of each coordinate
  xmin = minval(x)
  xmax = maxval(x)

  zmin = minval(z)
  zmax = maxval(z)

  xsize = xmax - xmin
  zsize = zmax - zmin

  print *,'Xmin and Xmax of the mesh read = ',xmin,xmax
  print *,'Zmin and Zmax of the mesh read = ',zmin,zmax
  print *

  print *,'Size of the mesh read along X = ',xsize
  print *,'Size of the mesh read along Z = ',zsize
  print *

  count_elem_faces_to_extend = 0

    if (iloop_on_X_Z_faces == 1) then ! Xmin or Xmax face
      value_min = xmin
      value_max = xmax
      value_size = xsize
      coord_to_use1 => x  ! make coordinate array to use point to array x()
      coord_to_use2 => z
    else if (iloop_on_X_Z_faces == 2) then ! Zmin or Zmax face
      value_min = zmin
      value_max = zmax
      value_size = zsize
      coord_to_use1 => z  ! make coordinate array to use point to array z()
      coord_to_use2 => x
    else
      stop 'wrong index in loop on faces'
    endif

  if (minval(ibool) /= 1) stop 'error in minval(ibool)'

  sum_of_distances = 0.d0

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)

      if (iloop_on_min_face_then_max_face == 1) then ! min face

! detect elements belonging to the min face
    limit = value_min + value_size * SMALL_RELATIVE_VALUE

! test face 1 (bottom)
    if (coord_to_use1(i1) < limit .and. coord_to_use1(i2) < limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + abs(coord_to_use2(i1) - coord_to_use2(i2))
    endif

! test face 2 (top)
    if (coord_to_use1(i3) < limit .and. coord_to_use1(i4) < limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + abs(coord_to_use2(i3) - coord_to_use2(i4))
    endif

! test face 3 (left)
    if (coord_to_use1(i1) < limit .and. coord_to_use1(i4) < limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + abs(coord_to_use2(i1) - coord_to_use2(i4))
    endif

! test face 4 (right)
    if (coord_to_use1(i2) < limit .and. coord_to_use1(i3) < limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + abs(coord_to_use2(i2) - coord_to_use2(i3))
    endif

      else ! max face

! detect elements belonging to the max face
    limit = value_max - value_size * SMALL_RELATIVE_VALUE

! test face 1 (bottom)
    if (coord_to_use1(i1) > limit .and. coord_to_use1(i2) > limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + abs(coord_to_use2(i1) - coord_to_use2(i2))
    endif

! test face 2 (top)
    if (coord_to_use1(i3) > limit .and. coord_to_use1(i4) > limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + abs(coord_to_use2(i3) - coord_to_use2(i4))
    endif

! test face 3 (left)
    if (coord_to_use1(i1) > limit .and. coord_to_use1(i4) > limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + abs(coord_to_use2(i1) - coord_to_use2(i4))
    endif

! test face 4 (right)
    if (coord_to_use1(i2) > limit .and. coord_to_use1(i3) > limit) then
      count_elem_faces_to_extend = count_elem_faces_to_extend + 1
      sum_of_distances = sum_of_distances + abs(coord_to_use2(i2) - coord_to_use2(i3))
    endif

      endif

  enddo

  print *,'Total number of elements in the mesh before extension = ',nspec
  print *,'Number of element faces to extend  = ',count_elem_faces_to_extend
  if (count_elem_faces_to_extend == 0) stop 'error: number of element faces to extend detected is zero!'

! we will add TOTAL_NUMBER_OF_LAYERS_TO_ADD to each of the element faces detected that need to be extended
  nspec_new = nspec + count_elem_faces_to_extend * TOTAL_NUMBER_OF_LAYERS_TO_ADD
  print *,'Total number of elements in the mesh after extension = ',nspec_new

! and each of these elements will have NGNOD points
! (some of them shared with other elements, but we will remove the multiples below, thus here it is a maximum
  npoin_new_max = npoin + count_elem_faces_to_extend * TOTAL_NUMBER_OF_LAYERS_TO_ADD * NGNOD

  mean_distance = sum_of_distances / dble(count_elem_faces_to_extend)
  very_small_distance_squared = (mean_distance / 10000.d0)**2
  if (icompute_size == 1) print *,'Computed mean size of the elements to extend = ',mean_distance
  print *

! allocate a new set of elements, i.e. a new ibool()
! allocate arrays with the right size of the future extended mesh
  allocate(imaterial_new(nspec_new))
  allocate(ibool_new(NGNOD,nspec_new))

! clear the arrays
  imaterial_new(:) = 0
  ibool_new(:,:) = 0

! copy the original arrays into the first part i.e. the non-extended part of the new ones
  ibool_new(:,1:nspec) = ibool(:,1:nspec)
  imaterial_new(1:nspec) = imaterial(1:nspec)

  if (minval(ibool) /= 1) stop 'error in minval(ibool)'

! allocate a new set of points, with multiples
  allocate(x_new(npoin_new_max))
  allocate(z_new(npoin_new_max))

! copy the original points into the new set
  x_new(1:npoin) = x(1:npoin)
  z_new(1:npoin) = z(1:npoin)

! position after which to start to create the new elements
  elem_counter = nspec
  npoin_new_real = npoin

! loop on the whole original mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)

! reset flag
    need_to_extend_this_element = .false.

      if (iloop_on_min_face_then_max_face == 1) then ! min face

! detect elements belonging to the min face
    limit = value_min + value_size * SMALL_RELATIVE_VALUE

! test face 1 (bottom)
    if (coord_to_use1(i1) < limit .and. coord_to_use1(i2) < limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i2
    endif

! test face 2 (top)
    if (coord_to_use1(i3) < limit .and. coord_to_use1(i4) < limit) then
      need_to_extend_this_element = .true.
      p1 = i3
      p2 = i4
    endif

! test face 3 (left)
    if (coord_to_use1(i1) < limit .and. coord_to_use1(i4) < limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i4
    endif

! test face 4 (right)
    if (coord_to_use1(i2) < limit .and. coord_to_use1(i3) < limit) then
      need_to_extend_this_element = .true.
      p1 = i2
      p2 = i3
    endif

      else ! max face

! detect elements belonging to the max face
    limit = value_max - value_size * SMALL_RELATIVE_VALUE

! test face 1 (bottom)
    if (coord_to_use1(i1) > limit .and. coord_to_use1(i2) > limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i2
    endif

! test face 2 (top)
    if (coord_to_use1(i3) > limit .and. coord_to_use1(i4) > limit) then
      need_to_extend_this_element = .true.
      p1 = i3
      p2 = i4
    endif

! test face 3 (left)
    if (coord_to_use1(i1) > limit .and. coord_to_use1(i4) > limit) then
      need_to_extend_this_element = .true.
      p1 = i1
      p2 = i4
    endif

! test face 4 (right)
    if (coord_to_use1(i2) > limit .and. coord_to_use1(i3) > limit) then
      need_to_extend_this_element = .true.
      p1 = i2
      p2 = i3
    endif

      endif

    if (need_to_extend_this_element) then

! create the TOTAL_NUMBER_OF_LAYERS_TO_ADD new elements

    factor_x = 0
    factor_z = 0

    SIZE_OF_X_ELEMENT_TO_ADD = 0.d0
    SIZE_OF_Z_ELEMENT_TO_ADD = 0.d0

    if (iloop_on_X_Z_faces == 1) then  ! Xmin or Xmax
      if (iloop_on_min_face_then_max_face == 1) then ! min face
        factor_x = -1
        if (icompute_size == 1) SIZE_OF_XMIN_ELEMENT_TO_ADD = mean_distance
        SIZE_OF_X_ELEMENT_TO_ADD = SIZE_OF_XMIN_ELEMENT_TO_ADD
      else ! max face
        factor_x = +1
        if (icompute_size == 1) SIZE_OF_XMAX_ELEMENT_TO_ADD = mean_distance
        SIZE_OF_X_ELEMENT_TO_ADD = SIZE_OF_XMAX_ELEMENT_TO_ADD
      endif
    else if (iloop_on_X_Z_faces == 2) then
      if (iloop_on_min_face_then_max_face == 1) then ! min face
        factor_z = -1
        if (icompute_size == 1) SIZE_OF_ZMIN_ELEMENT_TO_ADD = mean_distance
        SIZE_OF_Z_ELEMENT_TO_ADD = SIZE_OF_ZMIN_ELEMENT_TO_ADD
      else ! max face
        factor_z = +1
        if (icompute_size == 1) SIZE_OF_ZMAX_ELEMENT_TO_ADD = mean_distance
        SIZE_OF_Z_ELEMENT_TO_ADD = SIZE_OF_ZMAX_ELEMENT_TO_ADD
      endif
    else
      stop 'wrong index in loop on faces'
    endif

      do iextend = 1,TOTAL_NUMBER_OF_LAYERS_TO_ADD

        ! create a new element
        elem_counter = elem_counter + 1

        ! use the same material property for the extended elements as for the element being extended
        imaterial_new(elem_counter) = imaterial(ispec)

        ! create a new point if it does not exist yet, otherwise use the existing one to avoid creating multiples
        x_value_to_create = x(p1) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
        z_value_to_create = z(p1) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)
        point_already_exists = .false.
        ! this loop to avoid creating multiples is slow because it is a large loop for each point
        ! inside an existing loop on all the new points to create; it is usually OK for 2D meshes though.
        ! for a much faster version based on a quick sorting routine, see utils/CPML/add_CPML_layers_to_an_existing_mesh.f90
        ! in the 3D package (SPECFEM3D).
        do ipoin = 1,npoin_new_real
          distance_squared = (x_new(ipoin) - x_value_to_create)**2 + (z_new(ipoin) - z_value_to_create)**2
          if (distance_squared < very_small_distance_squared) then
            point_already_exists = .true.
            exit
          endif
        enddo
        if (point_already_exists) then
          ibool_new(1,elem_counter) = ipoin
        else
          npoin_new_real = npoin_new_real + 1
          ibool_new(1,elem_counter) = npoin_new_real
          x_new(npoin_new_real) = x_value_to_create
          z_new(npoin_new_real) = z_value_to_create
        endif

        ! create a new point if it does not exist yet, otherwise use the existing one to avoid creating multiples
        x_value_to_create = x(p1) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
        z_value_to_create = z(p1) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend
         point_already_exists = .false.
        ! this loop to avoid creating multiples is slow because it is a large loop for each point
        ! inside an existing loop on all the new points to create; it is usually OK for 2D meshes though.
        ! for a much faster version based on a quick sorting routine, see utils/CPML/add_CPML_layers_to_an_existing_mesh.f90
        ! in the 3D package (SPECFEM3D).
        do ipoin = 1,npoin_new_real
          distance_squared = (x_new(ipoin) - x_value_to_create)**2 + (z_new(ipoin) - z_value_to_create)**2
          if (distance_squared < very_small_distance_squared) then
            point_already_exists = .true.
            exit
          endif
        enddo
        if (point_already_exists) then
          ibool_new(2,elem_counter) = ipoin
        else
          npoin_new_real = npoin_new_real + 1
          ibool_new(2,elem_counter) = npoin_new_real
          x_new(npoin_new_real) = x_value_to_create
          z_new(npoin_new_real) = z_value_to_create
        endif

        ! create a new point if it does not exist yet, otherwise use the existing one to avoid creating multiples
        x_value_to_create = x(p2) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*iextend
        z_value_to_create = z(p2) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*iextend
        point_already_exists = .false.
        ! this loop to avoid creating multiples is slow because it is a large loop for each point
        ! inside an existing loop on all the new points to create; it is usually OK for 2D meshes though.
        ! for a much faster version based on a quick sorting routine, see utils/CPML/add_CPML_layers_to_an_existing_mesh.f90
        ! in the 3D package (SPECFEM3D).
        do ipoin = 1,npoin_new_real
          distance_squared = (x_new(ipoin) - x_value_to_create)**2 + (z_new(ipoin) - z_value_to_create)**2
          if (distance_squared < very_small_distance_squared) then
            point_already_exists = .true.
            exit
          endif
        enddo
        if (point_already_exists) then
          ibool_new(3,elem_counter) = ipoin
        else
          npoin_new_real = npoin_new_real + 1
          ibool_new(3,elem_counter) = npoin_new_real
          x_new(npoin_new_real) = x_value_to_create
          z_new(npoin_new_real) = z_value_to_create
        endif

        ! create a new point if it does not exist yet, otherwise use the existing one to avoid creating multiples
        x_value_to_create = x(p2) + factor_x*SIZE_OF_X_ELEMENT_TO_ADD*(iextend-1)
        z_value_to_create = z(p2) + factor_z*SIZE_OF_Z_ELEMENT_TO_ADD*(iextend-1)
        point_already_exists = .false.
        ! this loop to avoid creating multiples is slow because it is a large loop for each point
        ! inside an existing loop on all the new points to create; it is usually OK for 2D meshes though.
        ! for a much faster version based on a quick sorting routine, see utils/CPML/add_CPML_layers_to_an_existing_mesh.f90
        ! in the 3D package (SPECFEM3D).
        do ipoin = 1,npoin_new_real
          distance_squared = (x_new(ipoin) - x_value_to_create)**2 + (z_new(ipoin) - z_value_to_create)**2
          if (distance_squared < very_small_distance_squared) then
            point_already_exists = .true.
            exit
          endif
        enddo
        if (point_already_exists) then
          ibool_new(4,elem_counter) = ipoin
        else
          npoin_new_real = npoin_new_real + 1
          ibool_new(4,elem_counter) = npoin_new_real
          x_new(npoin_new_real) = x_value_to_create
          z_new(npoin_new_real) = z_value_to_create
        endif

! now we need to test if the element created is flipped i.e. it has a negative Jacobian,
! and if so we will use the mirrored version of that element, which will then have a positive Jacobian

! check the element for a negative Jacobian
      xelm(1) = x_new(ibool_new(1,elem_counter))
      xelm(2) = x_new(ibool_new(2,elem_counter))
      xelm(3) = x_new(ibool_new(3,elem_counter))
      xelm(4) = x_new(ibool_new(4,elem_counter))

      zelm(1) = z_new(ibool_new(1,elem_counter))
      zelm(2) = z_new(ibool_new(2,elem_counter))
      zelm(3) = z_new(ibool_new(3,elem_counter))
      zelm(4) = z_new(ibool_new(4,elem_counter))

      call calc_jacobian(xelm,zelm,dershape2D,found_a_negative_jacobian1,NDIM,NGNOD,NGLLX,NGLLZ)

! check the mirrored (i.e. flipped/swapped) element for a negative Jacobian
! either this one or the non-mirrored one above should be OK, and thus we will select it
      xelm(1) = x_new(ibool_new(2,elem_counter))
      xelm(2) = x_new(ibool_new(1,elem_counter))
      xelm(3) = x_new(ibool_new(4,elem_counter))
      xelm(4) = x_new(ibool_new(3,elem_counter))

      zelm(1) = z_new(ibool_new(2,elem_counter))
      zelm(2) = z_new(ibool_new(1,elem_counter))
      zelm(3) = z_new(ibool_new(4,elem_counter))
      zelm(4) = z_new(ibool_new(3,elem_counter))

      call calc_jacobian(xelm,zelm,dershape2D,found_a_negative_jacobian2,NDIM,NGNOD,NGLLX,NGLLZ)

! this should never happen, it is just a safety test
      if (found_a_negative_jacobian1 .and. found_a_negative_jacobian2) &
        stop 'error: found a negative Jacobian that could not be fixed'

! this should also never happen, it is just a second safety test
      if (.not. found_a_negative_jacobian1 .and. .not. found_a_negative_jacobian2) &
        stop 'strange error: both the element created and its mirrored version have a positive Jacobian!'

! if we have found that the original element has a negative Jacobian and its mirrored element is fine,
! swap the points so that we use that mirrored element in the final mesh saved to disk instead of the original one
      if (found_a_negative_jacobian1) then
        i1 = ibool_new(2,elem_counter)
        i2 = ibool_new(1,elem_counter)
        i3 = ibool_new(4,elem_counter)
        i4 = ibool_new(3,elem_counter)

        ibool_new(1,elem_counter) = i1
        ibool_new(2,elem_counter) = i2
        ibool_new(3,elem_counter) = i3
        ibool_new(4,elem_counter) = i4
      endif

      enddo
    endif

  enddo

  if (minval(ibool_new) /= 1) stop 'error in minval(ibool_new)'
  if (maxval(ibool_new) > npoin_new_max) stop 'error in maxval(ibool_new)'

! deallocate the original arrays
  deallocate(x,z)
  deallocate(ibool)
  deallocate(imaterial)

! reallocate them with the new size
  allocate(x(npoin_new_real))
  allocate(z(npoin_new_real))
  allocate(imaterial(nspec_new))
  allocate(ibool(NGNOD,nspec_new))

! make the new ones become the old ones, to prepare for the next iteration of the two nested loops we are in,
! i.e. to make sure the next loop will extend the mesh from the new arrays rather than from the old ones
  x(:) = x_new(1:npoin_new_real)
  z(:) = z_new(1:npoin_new_real)
  imaterial(:) = imaterial_new(:)
  ibool(:,:) = ibool_new(:,:)

! the new number of elements and points becomes the old one, for the same reason
  nspec = nspec_new
  npoin = npoin_new_real

! deallocate the new ones, to make sure they can be allocated again in the next iteration of the nested loops we are in
  deallocate(x_new,z_new)
  deallocate(ibool_new)
  deallocate(imaterial_new)

  if (minval(ibool) /= 1) stop 'error in minval(ibool)'
  if (maxval(ibool) > npoin) stop 'error in maxval(ibool)'

    enddo ! of iloop_on_min_face_then_max_face loop on Xmin then Xmax, or Zmin then Zmax

! end of loop on the three sets of faces to first add CPML elements along X, then along Z
  enddo

! write the new points (overwrite the old file)
  open(unit=23,file='nodes_coords_file',status='old',action='write')
  write(23,*) npoin
  do ipoin = 1,npoin
    write(23,*) sngl(x(ipoin)),sngl(z(ipoin))
  enddo
  close(23)

! write the new mesh elements (overwrite the old file)
  open(unit=23,file='mesh_file',status='old',action='write')
  write(23,*) nspec
! loop on the whole mesh
  do ispec = 1,nspec
    write(23,"(i9,1x,i9,1x,i9,1x,i9)") (ibool(ia,ispec), ia = 1,NGNOD)
  enddo
  close(23)

! write the new material properties (overwrite the old file)
  open(unit=23,file='materials_file',status='old',action='write')
! loop on the whole mesh
  do ispec = 1,nspec
    write(23,*) imaterial(ispec)
  enddo
  close(23)

! output information for the next code (xconvert_external_layers_of_a_given_mesh_to_CPML_layers)
  print *
  print *,'Here are the values to use as input in the next code, xconvert_external_layers_of_a_given_mesh_to_CPML_layers:'
  print *,'  (also saved in file values_to_use_for_convert_external_layers_of_a_given_mesh_to_CPML_layers.txt)'
  print *
  if (ADD_ON_THE_XMIN_SURFACE) print *,'THICKNESS_OF_XMIN_PML = ',sngl(SIZE_OF_XMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  if (ADD_ON_THE_XMAX_SURFACE) print *,'THICKNESS_OF_XMAX_PML = ',sngl(SIZE_OF_XMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  if (ADD_ON_THE_ZMIN_SURFACE) print *,'THICKNESS_OF_ZMIN_PML = ',sngl(SIZE_OF_ZMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  if (ADD_ON_THE_ZMAX_SURFACE) print *,'THICKNESS_OF_ZMAX_PML = ',sngl(SIZE_OF_ZMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD)
  print *

! save the thickness values to a text file
  open(unit=23,file='values_to_use_for_convert_external_layers_of_a_given_mesh_to_CPML_layers.txt',status='unknown',action='write')

  if (ADD_ON_THE_XMIN_SURFACE) then
    write(23,*) SIZE_OF_XMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  else
! convention (negative value) to say that this absorbing edge is turned off
    write(23,*) -1
  endif

  if (ADD_ON_THE_XMAX_SURFACE) then
    write(23,*) SIZE_OF_XMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  else
! convention (negative value) to say that this absorbing edge is turned off
    write(23,*) -1
  endif

  if (ADD_ON_THE_ZMIN_SURFACE) then
    write(23,*) SIZE_OF_ZMIN_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  else
! convention (negative value) to say that this absorbing edge is turned off
    write(23,*) -1
  endif

  if (ADD_ON_THE_ZMAX_SURFACE) then
    write(23,*) SIZE_OF_ZMAX_ELEMENT_TO_ADD * NUMBER_OF_PML_LAYERS_TO_ADD
  else
! convention (negative value) to say that this absorbing edge is turned off
    write(23,*) -1
  endif

  close(23)

  end program add_CPML_layers_to_a_given_mesh

!
!=====================================================================
!

  subroutine calc_jacobian(xelm,zelm,dershape2D,found_a_negative_jacobian,NDIM,NGNOD,NGLLX,NGLLZ)

  implicit none

  integer :: NDIM,NGNOD,NGLLX,NGLLZ

  logical :: found_a_negative_jacobian

  double precision, dimension(NGNOD) :: xelm,zelm
  double precision dershape2D(NDIM,NGNOD,NGLLX,NGLLZ)

  integer i,k,ia
  double precision xxi,xgamma,zxi,zgamma
  double precision jacobian

  double precision, parameter :: ZERO = 0.d0

  found_a_negative_jacobian = .false.

! do k = 1,NGLLZ
!   do i = 1,NGLLX
! for this CPML mesh extrusion routine it is sufficient to test the 4 corners of each element to reduce the cost
! because we just want to detect if the element is flipped or not, and if so flip it back
! do k = 1,NGLLZ,NGLLZ-1
!   do i = 1,NGLLX,NGLLX-1
! it is even sufficient to test a single corner
  do k = 1,1
    do i = 1,1

      xxi = ZERO
      zxi = ZERO
      xgamma = ZERO
      zgamma = ZERO

      do ia=1,NGNOD
        xxi = xxi + dershape2D(1,ia,i,k)*xelm(ia)
        zxi = zxi + dershape2D(1,ia,i,k)*zelm(ia)
        xgamma = xgamma + dershape2D(2,ia,i,k)*xelm(ia)
        zgamma = zgamma + dershape2D(2,ia,i,k)*zelm(ia)
      enddo

      jacobian = xxi*zgamma - xgamma*zxi

! check that the Jacobian transform is invertible, i.e. that the Jacobian never becomes negative or null
      if (jacobian <= ZERO) found_a_negative_jacobian = .true.

    enddo
  enddo

  end subroutine calc_jacobian

!
!=====================================================================
!

  subroutine define_shape_functions(shape2D,dershape2D,xi,gamma,NGNOD)

!=======================================================================
!
!  Set up the shape functions for the superparametric transformation
! (i.e. the geometry is defined with lower-order functions than the unknowns
!  of the problem; see for instance Chapter 16 of the finite-element course of
!  Carlos Felippa in Colorado for a discussion about this).
!  The routine can handle 4 or 9 control nodes defined as follows:
!
!                               4 . . . . 7 . . . . 3
!                               .                   .
!                               .         gamma     .
!                               .                   .
!                               8         9  xi     6
!                               .                   .
!                               .                   .
!                               .                   .
!                               1 . . . . 5 . . . . 2
!
!                           Local coordinate system : s,t
!
!=======================================================================

  implicit none

! we are in 2D
  integer, parameter :: NDIM = 2

! a few useful constants
  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0
  double precision, parameter :: HALF = 0.5d0,TWO = 2.d0,QUARTER = 0.25d0

! very large and very small values
  double precision, parameter :: TINYVAL = 1.d-9

  integer,intent(in) :: NGNOD

  double precision,intent(out) :: shape2D(NGNOD)
  double precision,intent(out) :: dershape2D(NDIM,NGNOD)
  double precision,intent(in) :: xi,gamma

  ! local parameters
  double precision :: s,t,sp,sm,tp,tm,s2,t2,ss,tt,st

!
!---- set up the shape functions and their local derivatives
!
  s  = xi
  t  = gamma

!----    4-node element
  if (NGNOD == 4) then
       sp = s + ONE
       sm = s - ONE
       tp = t + ONE
       tm = t - ONE

!----  corner nodes
       shape2D(1) = QUARTER * sm * tm
       shape2D(2) = - QUARTER * sp * tm
       shape2D(3) = QUARTER * sp * tp
       shape2D(4) = - QUARTER * sm * tp

       dershape2D(1,1) = QUARTER * tm
       dershape2D(1,2) = - QUARTER * tm
       dershape2D(1,3) =  QUARTER * tp
       dershape2D(1,4) = - QUARTER * tp

       dershape2D(2,1) = QUARTER * sm
       dershape2D(2,2) = - QUARTER * sp
       dershape2D(2,3) =  QUARTER * sp
       dershape2D(2,4) = - QUARTER * sm

!----    9-node element
  else if (NGNOD == 9) then

       sp = s + ONE
       sm = s - ONE
       tp = t + ONE
       tm = t - ONE
       s2 = s * TWO
       t2 = t * TWO
       ss = s * s
       tt = t * t
       st = s * t

!----  corner nodes
       shape2D(1) = QUARTER * sm * st * tm
       shape2D(2) = QUARTER * sp * st * tm
       shape2D(3) = QUARTER * sp * st * tp
       shape2D(4) = QUARTER * sm * st * tp

       dershape2D(1,1) = QUARTER * tm * t * (s2 - ONE)
       dershape2D(1,2) = QUARTER * tm * t * (s2 + ONE)
       dershape2D(1,3) = QUARTER * tp * t * (s2 + ONE)
       dershape2D(1,4) = QUARTER * tp * t * (s2 - ONE)

       dershape2D(2,1) = QUARTER * sm * s * (t2 - ONE)
       dershape2D(2,2) = QUARTER * sp * s * (t2 - ONE)
       dershape2D(2,3) = QUARTER * sp * s * (t2 + ONE)
       dershape2D(2,4) = QUARTER * sm * s * (t2 + ONE)

!----  midside nodes
       shape2D(5) = HALF * tm * t * (ONE - ss)
       shape2D(6) = HALF * sp * s * (ONE - tt)
       shape2D(7) = HALF * tp * t * (ONE - ss)
       shape2D(8) = HALF * sm * s * (ONE - tt)

       dershape2D(1,5) = -ONE  * st * tm
       dershape2D(1,6) =  HALF * (ONE - tt) * (s2 + ONE)
       dershape2D(1,7) = -ONE  * st * tp
       dershape2D(1,8) =  HALF * (ONE - tt) * (s2 - ONE)

       dershape2D(2,5) =  HALF * (ONE - ss) * (t2 - ONE)
       dershape2D(2,6) = -ONE  * st * sp
       dershape2D(2,7) =  HALF * (ONE - ss) * (t2 + ONE)
       dershape2D(2,8) = -ONE  * st * sm

!----  center node
       shape2D(9) = (ONE - ss) * (ONE - tt)

       dershape2D(1,9) = -ONE * s2 * (ONE - tt)
       dershape2D(2,9) = -ONE * t2 * (ONE - ss)

  else
     stop 'Error: wrong number of control nodes'
  endif

!--- check the shape functions and their derivatives
  ! sum of shape functions should be one
  if (abs(sum(shape2D)-ONE) > TINYVAL) stop 'Error shape functions'

  ! sum of derivatives of shape functions should be zero
  if (abs(sum(dershape2D(1,:))) > TINYVAL) stop 'Error deriv xi shape functions'
  if (abs(sum(dershape2D(2,:))) > TINYVAL) stop 'Error deriv gamma shape functions'

  end subroutine define_shape_functions

