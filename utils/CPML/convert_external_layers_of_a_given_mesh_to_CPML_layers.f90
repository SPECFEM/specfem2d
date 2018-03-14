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

  program convert_mesh_to_CPML

! Dimitri Komatitsch, CNRS Marseille, France, January 2017

! convert outer layers of an existing CUBIT mesh file stored in SPECFEM2D format to CPML layers,
! i.e., create a 'absorbing_cpml_file' file for that existing mesh

  implicit none

! this is for QUAD4; the code below also works for QUAD9, it then just uses the first 4 points of an element
! to determine if it belongs to a CPML layer
  integer, parameter :: NGNOD = 4

  integer :: nspec,npoin
  integer :: ispec,ipoin
  integer :: iflag,iread
  integer :: i1,i2,i3,i4,number_of_CPML_elements,count_faces_found

  double precision, dimension(:), allocatable :: x,z

  logical, dimension(:), allocatable :: is_X_CPML,is_Z_CPML

  double precision :: xmin,xmax,zmin,zmax,limit,size_of_model

  logical :: ADD_ON_THE_XMIN_SURFACE,ADD_ON_THE_XMAX_SURFACE
  logical :: ADD_ON_THE_ZMIN_SURFACE,ADD_ON_THE_ZMAX_SURFACE

  logical :: already_found_a_face

  double precision :: THICKNESS_OF_XMIN_PML,THICKNESS_OF_ZMIN_PML
  double precision :: THICKNESS_OF_XMAX_PML,THICKNESS_OF_ZMAX_PML

  integer, dimension(:,:), allocatable :: ibool

! to make sure coordinate roundoff problems do not occur, use a tolerance of 0.5%
  double precision, parameter :: SMALL_PERCENTAGE_TOLERANCE = 1.005d0

  double precision, parameter :: SMALL_RELATIVE_VALUE = 0.5d-3

! flags for CPML absorbing boundaries
  integer, parameter :: CPML_X_ONLY = 1
  integer, parameter :: CPML_Z_ONLY = 2
  integer, parameter :: CPML_XZ = 3

  print *
  print *,'IMPORTANT: it is your responsibility to make sure that in the input CUBIT (or similar) mesh'
  print *,'that this code will read in SPECFEM3D format from files "nodes_coords_file" and "mesh_file"'
  print *,'you have created layers of elements that constitute a layer of constant thickness aligned'
  print *,'with the coordinate grid axes (X and/or Z), so that this code can assign CPML flags to them.'
  print *,'This code does NOT check that (because it cannot, in any easy way).'
  print *,'The mesh inside these CPML layers does not need to be structured nor regular, any non-structured'
  print *,'mesh is fine as long as it has flat PML inner and outer faces, parallel to the axes, and thus'
  print *,'of a constant thickness.'
  print *,'The thickness can be different for the X and Z sides. But for X it must not vary, and for Z it must not vary.'
  print *,'If you do not know the exact thickness, you can use a slightly LARGER value'
  print *,'in this code (say 2% to 5% more) and this code will fix that and will adjust it;'
  print *,'never use a SMALLER value otherwise this code will miss some CPML elements.'
  print *
  print *,'Note: in a future release we will remove the constraint of having CPML layers aligned with the'
  print *,'coordinate axes; we will allow for meshes that are titled by any constant angle in the horizontal plane.'
  print *,'However this is not implemented yet.'
  print *

  print *,'1 = enter the CPML thickness values to create manually'
  print *,'2 = read them from a file created by the previous code, xadd_CPML_layers_to_an_existing_mesh'
  print *,'3 = exit'
  read(*,*) iread
  if (iread /= 1 .and. iread /= 2) stop 'exiting...'

  if (iread /= 2) then
    ADD_ON_THE_XMIN_SURFACE = .true.
    ADD_ON_THE_XMAX_SURFACE = .true.
    ADD_ON_THE_ZMIN_SURFACE = .true.
    ADD_ON_THE_ZMAX_SURFACE = .true.

    print *
    print *,'1 = use a free surface at the top of the mesh (most classical option)'
    print *,'2 = use a CPML absorbing layer at the top of the mesh (less classical option)'
    print *,'3 = exit'
    read(*,*) iflag
    if (iflag /= 1 .and. iflag /= 2) stop 'exiting...'
    if (iflag == 1) then
      ADD_ON_THE_ZMAX_SURFACE = .false.
    else
      ADD_ON_THE_ZMAX_SURFACE = .true.
    endif
  endif

  print *

! open SPECFEM2D mesh file to read the points
  open(unit=23,file='nodes_coords_file',status='old',action='read')
  read(23,*) npoin
  allocate(x(npoin))
  allocate(z(npoin))
  do ipoin = 1,npoin
    read(23,*) x(ipoin),z(ipoin)
  enddo
  close(23)

! compute the min and max values of each coordinate
  xmin = minval(x)
  xmax = maxval(x)

  zmin = minval(z)
  zmax = maxval(z)

  print *,'Xmin and Xmax of the mesh read = ',xmin,xmax
  print *,'Zmin and Zmax of the mesh read = ',zmin,zmax
  print *

  if (iread == 2) then

! read the thickness values from an existing text file
  open(unit=23,file='values_to_use_for_convert_external_layers_of_a_given_mesh_to_CPML_layers.txt',status='old',action='read')
  read(23,*) THICKNESS_OF_XMIN_PML
  read(23,*) THICKNESS_OF_XMAX_PML
  read(23,*) THICKNESS_OF_ZMIN_PML
  read(23,*) THICKNESS_OF_ZMAX_PML
  close(23)

! use the convention that a negative value means that that PML is turned off
  ADD_ON_THE_XMIN_SURFACE = .true.
  ADD_ON_THE_XMAX_SURFACE = .true.
  ADD_ON_THE_ZMIN_SURFACE = .true.
  ADD_ON_THE_ZMAX_SURFACE = .true.
  if (THICKNESS_OF_XMIN_PML <= 0) ADD_ON_THE_XMIN_SURFACE = .false.
  if (THICKNESS_OF_XMAX_PML <= 0) ADD_ON_THE_XMAX_SURFACE = .false.
  if (THICKNESS_OF_ZMIN_PML <= 0) ADD_ON_THE_ZMIN_SURFACE = .false.
  if (THICKNESS_OF_ZMAX_PML <= 0) ADD_ON_THE_ZMAX_SURFACE = .false.

! check convention (negative value) that says that this absorbing edge is turned off
  if (ADD_ON_THE_XMIN_SURFACE .and. THICKNESS_OF_XMIN_PML <= 0) &
    stop 'negative thickness is not allowed; ADD_ON_THE_XMIN_SURFACE is maybe inconsistent with the previous code; exiting...'
  if (.not. ADD_ON_THE_XMIN_SURFACE .and. THICKNESS_OF_XMIN_PML > 0) &
    stop 'ADD_ON_THE_XMIN_SURFACE seems inconsistent with the previous code; exiting...'

! check convention (negative value) that says that this absorbing edge is turned off
  if (ADD_ON_THE_XMAX_SURFACE .and. THICKNESS_OF_XMAX_PML <= 0) &
    stop 'negative thickness is not allowed; ADD_ON_THE_XMAX_SURFACE is maybe inconsistent with the previous code; exiting...'
  if (.not. ADD_ON_THE_XMAX_SURFACE .and. THICKNESS_OF_XMAX_PML > 0) &
    stop 'ADD_ON_THE_XMAX_SURFACE seems inconsistent with the previous code; exiting...'

! check convention (negative value) that says that this absorbing edge is turned off
  if (ADD_ON_THE_ZMIN_SURFACE .and. THICKNESS_OF_ZMIN_PML <= 0) &
    stop 'negative thickness is not allowed; ADD_ON_THE_ZMIN_SURFACE is maybe inconsistent with the previous code; exiting...'
  if (.not. ADD_ON_THE_ZMIN_SURFACE .and. THICKNESS_OF_ZMIN_PML > 0) &
    stop 'ADD_ON_THE_ZMIN_SURFACE seems inconsistent with the previous code; exiting...'

! check convention (negative value) that says that this absorbing edge is turned off
  if (ADD_ON_THE_ZMAX_SURFACE .and. THICKNESS_OF_ZMAX_PML <= 0) &
    stop 'negative thickness is not allowed; ADD_ON_THE_ZMAX_SURFACE is maybe inconsistent with the previous code; exiting...'
  if (.not. ADD_ON_THE_ZMAX_SURFACE .and. THICKNESS_OF_ZMAX_PML > 0) &
    stop 'ADD_ON_THE_ZMAX_SURFACE seems inconsistent with the previous code; exiting...'

  else

  print *,'What is the exact thickness of the PML layer that you want (enter -1 to turn that PML layer off)'
  print *,'on the Xmin face of your mesh? (it needs to correspond exactly'
  print *,'to the flat layers you created in your input CUBIT mesh, as mentioned in'
  print *,'the comment printed above; if you think you have roundoff issues or very'
  print *,'slightly varying thickness, give 2% or 5% more here, but never less'
  read(*,*) THICKNESS_OF_XMIN_PML
  if (THICKNESS_OF_XMIN_PML <= 0) then
    THICKNESS_OF_XMIN_PML = 0
    ADD_ON_THE_XMIN_SURFACE = .false.
  else if (THICKNESS_OF_XMIN_PML > 0.30*(xmax - xmin)) then
    stop 'thickness of each CPML layer greater than 30% of the size of the mesh is not a good idea; exiting...'
  endif
  print *

  print *,'What is the exact thickness of the PML layer that you want (enter -1 to turn that PML layer off)'
  print *,'on the Xmax face of your mesh?'
  read(*,*) THICKNESS_OF_XMAX_PML
  if (THICKNESS_OF_XMAX_PML <= 0) then
    THICKNESS_OF_XMAX_PML = 0
    ADD_ON_THE_XMAX_SURFACE = .false.
  else if (THICKNESS_OF_XMAX_PML > 0.30*(xmax - xmin)) then
    stop 'thickness of each CPML layer greater than 30% of the size of the mesh is not a good idea; exiting...'
  endif
  print *

  print *,'What is the exact thickness of the PML layer that you want (enter -1 to turn that PML layer off)'
  print *,'on the Zmin face of your mesh?'
  read(*,*) THICKNESS_OF_ZMIN_PML
  if (THICKNESS_OF_ZMIN_PML <= 0) then
    THICKNESS_OF_ZMIN_PML = 0
    ADD_ON_THE_ZMIN_SURFACE = .false.
  else if (THICKNESS_OF_ZMIN_PML > 0.30*(zmax - zmin)) then
    stop 'thickness of each CPML layer greater than 30% of the size of the mesh is not a good idea; exiting...'
  endif
  print *

  if (ADD_ON_THE_ZMAX_SURFACE) then
    print *,'What is the exact thickness of the PML layer that you want (enter -1 to turn that PML layer off)'
    print *,'on the Zmax face of your mesh?'
    read(*,*) THICKNESS_OF_ZMAX_PML
    if (THICKNESS_OF_ZMAX_PML <= 0) then
      THICKNESS_OF_ZMAX_PML = 0
      ADD_ON_THE_ZMAX_SURFACE = .false.
    else if (THICKNESS_OF_ZMAX_PML > 0.30*(zmax - zmin)) then
      stop 'thickness of each CPML layer greater than 30% of the size of the mesh is not a good idea; exiting...'
    endif
    print *
  endif

  endif

! check that we need to create at least one PML, otherwise this code is useless
  if (.not. ADD_ON_THE_XMIN_SURFACE .and. .not. ADD_ON_THE_XMAX_SURFACE &
      .and. .not. ADD_ON_THE_ZMIN_SURFACE .and. .not. ADD_ON_THE_ZMAX_SURFACE) &
    stop 'Error: the purpose of this code is to create at least one PML, but you are creating none'

! ************* read mesh elements and generate CPML flags *************

! open SPECFEM2D topology file to read the mesh elements
  open(unit=23,file='mesh_file',status='old',action='read')
  read(23,*) nspec

  allocate(is_X_CPML(nspec))
  allocate(is_Z_CPML(nspec))

  is_X_CPML(:) = .false.
  is_Z_CPML(:) = .false.

  allocate(ibool(NGNOD,nspec))

! loop on the whole mesh
  do ispec = 1,nspec
    read(23,*) ibool(1,ispec),ibool(2,ispec),ibool(3,ispec),ibool(4,ispec)
  enddo

  close(23)

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)

! Xmin CPML
  if (ADD_ON_THE_XMIN_SURFACE) then
    limit = xmin + THICKNESS_OF_XMIN_PML * SMALL_PERCENTAGE_TOLERANCE
    if (x(i1) < limit .and. x(i2) < limit .and. x(i3) < limit .and. x(i4) < limit) is_X_CPML(ispec) = .true.
  endif

! Xmax CPML
  if (ADD_ON_THE_XMAX_SURFACE) then
    limit = xmax - THICKNESS_OF_XMAX_PML * SMALL_PERCENTAGE_TOLERANCE
    if (x(i1) > limit .and. x(i2) > limit .and. x(i3) > limit .and. x(i4) > limit) is_X_CPML(ispec) = .true.
  endif

! Zmin CPML
  if (ADD_ON_THE_ZMIN_SURFACE) then
    limit = zmin + THICKNESS_OF_ZMIN_PML * SMALL_PERCENTAGE_TOLERANCE
    if (z(i1) < limit .and. z(i2) < limit .and. z(i3) < limit .and. z(i4) < limit) is_Z_CPML(ispec) = .true.
  endif

! Zmax CPML
  if (ADD_ON_THE_ZMAX_SURFACE) then
    limit = zmax - THICKNESS_OF_ZMAX_PML * SMALL_PERCENTAGE_TOLERANCE
    if (z(i1) > limit .and. z(i2) > limit .and. z(i3) > limit .and. z(i4) > limit) is_Z_CPML(ispec) = .true.
  endif

  enddo

  print *,'Total number of elements in the mesh read = ',nspec
  print *
  print *,'Found ',count(is_X_CPML),' X_CPML elements'
  print *,'Found ',count(is_Z_CPML),' Z_CPML elements'
  number_of_CPML_elements = 0
  do ispec = 1,nspec
    if (is_X_CPML(ispec) .and. is_Z_CPML(ispec)) number_of_CPML_elements = number_of_CPML_elements + 1
  enddo
  print *,'of which ',number_of_CPML_elements,' are both X_CPML and Z_CPML elements'
  print *
  if (ADD_ON_THE_XMIN_SURFACE) print *,'    (converted the Xmin surface from free surface to CPML)'
  if (ADD_ON_THE_XMAX_SURFACE) print *,'    (converted the Xmax surface from free surface to CPML)'
  if (ADD_ON_THE_ZMIN_SURFACE) print *,'    (converted the Zmin surface from free surface to CPML)'
  if (ADD_ON_THE_ZMAX_SURFACE) print *,'    (converted the Zmax surface from free surface to CPML)'
  print *

  if (count(is_X_CPML) == 0 .and. count(is_Z_CPML) == 0) stop 'error: no CPML elements detected on any of the sides!'

  number_of_CPML_elements = 0
  do ispec = 1,nspec
    if (is_X_CPML(ispec) .or. is_Z_CPML(ispec)) number_of_CPML_elements = number_of_CPML_elements + 1
  enddo
  print *,'Created a total of ',number_of_CPML_elements,' unique CPML elements'
  print *,'   (i.e., ',100.*number_of_CPML_elements/real(nspec),'% of the mesh)'

! ************* generate the CPML database file *************

  open(unit=24,file='absorbing_cpml_file',status='unknown',action='write')

! write the total number of unique CPML elements
  write(24,*) number_of_CPML_elements

! write the CPML flag for each CPML element
  do ispec=1,nspec
    if (is_X_CPML(ispec) .and. is_Z_CPML(ispec)) then
      write(24,*) ispec,CPML_XZ

    else if (is_X_CPML(ispec)) then
      write(24,*) ispec,CPML_X_ONLY

    else if (is_Z_CPML(ispec)) then
      write(24,*) ispec,CPML_Z_ONLY
    endif

  enddo
  close(24)

  print *
  print *,'CPML absorbing layer file "absorbing_cpml_file" has been successfully created'
  print *

  !-----------------------------------------------
  ! Generate outer absorbing surface.
  !-----------------------------------------------

! first count the number of faces that are along that edge

  count_faces_found = 0

! ************* generate "absorbing_surface_file" for Xmin *************

! Xmin CPML
  size_of_model = xmax - xmin
  limit = xmin + SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)

    if (is_X_CPML(ispec)) then

      already_found_a_face = .false.

! test face 1 (bottom)
      if (x(i1) < limit .and. x(i2) < limit) then
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 2 (top)
      if (x(i3) < limit .and. x(i4) < limit) then
        if (already_found_a_face) then
          print *,'problem detected in the mesh near coordinates x,z = ',x(i3),z(i3)
          stop 'error: element with two faces on the same PML edge found!'
        endif
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 3 (left)
      if (x(i1) < limit .and. x(i4) < limit) then
        if (already_found_a_face) then
          print *,'problem detected in the mesh near coordinates x,z = ',x(i1),z(i1)
          stop 'error: element with two faces on the same PML edge found!'
        endif
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 4 (right)
      if (x(i2) < limit .and. x(i3) < limit) then
        if (already_found_a_face) then
          print *,'problem detected in the mesh near coordinates x,z = ',x(i2),z(i2)
          stop 'error: element with two faces on the same PML edge found!'
        endif
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

    endif

  enddo

! ************* generate "absorbing_surface_file" for Xmax *************

! Xmax CPML
  size_of_model = xmax - xmin
  limit = xmax - SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)

    if (is_X_CPML(ispec)) then

      already_found_a_face = .false.

! test face 1 (bottom)
      if (x(i1) > limit .and. x(i2) > limit) then
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 2 (top)
      if (x(i3) > limit .and. x(i4) > limit) then
        if (already_found_a_face) then
          print *,'problem detected in the mesh near coordinates x,z = ',x(i3),z(i3)
          stop 'error: element with two faces on the same PML edge found!'
        endif
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 3 (left)
      if (x(i1) > limit .and. x(i4) > limit) then
        if (already_found_a_face) then
          print *,'problem detected in the mesh near coordinates x,z = ',x(i1),z(i1)
          stop 'error: element with two faces on the same PML edge found!'
        endif
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 4 (right)
      if (x(i2) > limit .and. x(i3) > limit) then
        if (already_found_a_face) then
          print *,'problem detected in the mesh near coordinates x,z = ',x(i2),z(i2)
          stop 'error: element with two faces on the same PML edge found!'
        endif
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

    endif

  enddo

! ************* generate "absorbing_surface_file" for bottom *************

! Zmin CPML
  size_of_model = zmax - zmin
  limit = zmin + SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)

    if (is_Z_CPML(ispec)) then

      already_found_a_face = .false.

! test face 1 (bottom)
      if (z(i1) < limit .and. z(i2) < limit) then
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 2 (top)
      if (z(i3) < limit .and. z(i4) < limit) then
        if (already_found_a_face) then
          print *,'problem detected in the mesh near coordinates x,z = ',x(i3),z(i3)
          stop 'error: element with two faces on the same PML edge found!'
        endif
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 3 (left)
      if (z(i1) < limit .and. z(i4) < limit) then
        if (already_found_a_face) then
          print *,'problem detected in the mesh near coordinates x,z = ',x(i1),z(i1)
          stop 'error: element with two faces on the same PML edge found!'
        endif
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 4 (right)
      if (z(i2) < limit .and. z(i3) < limit) then
        if (already_found_a_face) then
          print *,'problem detected in the mesh near coordinates x,z = ',x(i2),z(i2)
          stop 'error: element with two faces on the same PML edge found!'
        endif
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

    endif

  enddo

! ************* generate "absorbing_surface_file" for top *************

! Zmax CPML
  size_of_model = zmax - zmin
  limit = zmax - SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)

    if (is_Z_CPML(ispec)) then

      already_found_a_face = .false.

! test face 1 (bottom)
      if (z(i1) > limit .and. z(i2) > limit) then
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 2 (top)
      if (z(i3) > limit .and. z(i4) > limit) then
        if (already_found_a_face) then
          print *,'problem detected in the mesh near coordinates x,z = ',x(i3),z(i3)
          stop 'error: element with two faces on the same PML edge found!'
        endif
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 3 (left)
      if (z(i1) > limit .and. z(i4) > limit) then
        if (already_found_a_face) then
          print *,'problem detected in the mesh near coordinates x,z = ',x(i1),z(i1)
          stop 'error: element with two faces on the same PML edge found!'
        endif
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

! test face 4 (right)
      if (z(i2) > limit .and. z(i3) > limit) then
        if (already_found_a_face) then
          print *,'problem detected in the mesh near coordinates x,z = ',x(i2),z(i2)
          stop 'error: element with two faces on the same PML edge found!'
        endif
        count_faces_found = count_faces_found + 1
        already_found_a_face = .true.
      endif

    endif

  enddo

  print *,'found ',count_faces_found,' full faces on all PML faces'

  !-----------------------------------------------
  ! Write the absorbing surface created.
  ! 'abs_surface' contains 1/ element number, 2/ number of nodes that form the absorbing edge
  ! (which currently must always be equal to 2),
  ! 3/ first node on the abs surface, 4/ second node on the abs surface
  ! 5/ 1=IBOTTOM, 2=IRIGHT, 3=ITOP, 4=ILEFT
  !-----------------------------------------------

! open SPECFEM2D topology file to read the mesh elements
  open(unit=24,file='absorbing_surface_file',status='unknown',action='write')

! write the total number of face elements
  write(24,*) count_faces_found

! ************* generate "absorbing_surface_file" for Xmin *************

! Xmin CPML
  size_of_model = xmax - xmin
  limit = xmin + SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)

    if (is_X_CPML(ispec)) then

! for the four faces below it is important to make sure we write the four points
! in an order for which the normal to the face points outwards

! test face 1 (bottom)
      if (x(i1) < limit .and. x(i2) < limit) write(24,*) ispec,' 2 ',i1,i2,' 4'

! test face 2 (top)
      if (x(i3) < limit .and. x(i4) < limit) write(24,*) ispec,' 2 ',i3,i4,' 4'

! test face 3 (left)
      if (x(i1) < limit .and. x(i4) < limit) write(24,*) ispec,' 2 ',i4,i1,' 4'

! test face 4 (right)
      if (x(i2) < limit .and. x(i3) < limit) write(24,*) ispec,' 2 ',i2,i3,' 4'

    endif

  enddo

! ************* generate "absorbing_surface_file" for Xmax *************

! Xmax CPML
  size_of_model = xmax - xmin
  limit = xmax - SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)

    if (is_X_CPML(ispec)) then

! for the four faces below it is important to make sure we write the four points
! in an order for which the normal to the face points outwards

! test face 1 (bottom)
      if (x(i1) > limit .and. x(i2) > limit) write(24,*) ispec,' 2 ',i1,i2,' 2'

! test face 2 (top)
      if (x(i3) > limit .and. x(i4) > limit) write(24,*) ispec,' 2 ',i3,i4,' 2'

! test face 3 (left)
      if (x(i1) > limit .and. x(i4) > limit) write(24,*) ispec,' 2 ',i4,i1,' 2'

! test face 4 (right)
      if (x(i2) > limit .and. x(i3) > limit) write(24,*) ispec,' 2 ',i2,i3,' 2'

    endif

  enddo

! ************* generate "absorbing_surface_file" for bottom *************

! Zmin CPML
  size_of_model = zmax - zmin
  limit = zmin + SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)

    if (is_Z_CPML(ispec)) then

! for the four faces below it is important to make sure we write the four points
! in an order for which the normal to the face points outwards

! test face 1 (bottom)
      if (z(i1) < limit .and. z(i2) < limit) write(24,*) ispec,' 2 ',i1,i2,' 1'

! test face 2 (top)
      if (z(i3) < limit .and. z(i4) < limit) write(24,*) ispec,' 2 ',i3,i4,' 1'

! test face 3 (left)
      if (z(i1) < limit .and. z(i4) < limit) write(24,*) ispec,' 2 ',i4,i1,' 1'

! test face 4 (right)
      if (z(i2) < limit .and. z(i3) < limit) write(24,*) ispec,' 2 ',i2,i3,' 1'

    endif

  enddo

! ************* generate "absorbing_surface_file" for top *************

! Zmax CPML
  size_of_model = zmax - zmin
  limit = zmax - SMALL_RELATIVE_VALUE*size_of_model

! loop on the whole mesh
  do ispec = 1,nspec

    i1 = ibool(1,ispec)
    i2 = ibool(2,ispec)
    i3 = ibool(3,ispec)
    i4 = ibool(4,ispec)

    if (is_Z_CPML(ispec)) then

! for the four faces below it is important to make sure we write the four points
! in an order for which the normal to the face points outwards

! test face 1 (bottom)
      if (z(i1) > limit .and. z(i2) > limit) write(24,*) ispec,' 2 ',i1,i2,' 3'

! test face 2 (top)
      if (z(i3) > limit .and. z(i4) > limit) write(24,*) ispec,' 2 ',i3,i4,' 3'

! test face 3 (left)
      if (z(i1) > limit .and. z(i4) > limit) write(24,*) ispec,' 2 ',i4,i1,' 3'

! test face 4 (right)
      if (z(i2) > limit .and. z(i3) > limit) write(24,*) ispec,' 2 ',i2,i3,' 3'

    endif

  enddo

  close(24)

  print *
  print *,'CPML file "absorbing_surface_file" has been successfully created'
  print *

! ************* create an empty "free_surface_file" file if all edges are PML
  if (ADD_ON_THE_XMIN_SURFACE .and. ADD_ON_THE_XMAX_SURFACE .and. ADD_ON_THE_ZMIN_SURFACE .and. ADD_ON_THE_ZMAX_SURFACE) then
    open(unit=24,file='free_surface_file',status='unknown',action='write')
    write(24,*) '0'
    close(24)
  endif

  end program convert_mesh_to_CPML

