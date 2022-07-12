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

! in the case of very large meshes, this option can be useful to switch from ASCII to binary for the mesh files
! #define USE_BINARY_FOR_EXTERNAL_MESH_DATABASE

  !-----------------------------------------------
  ! Read the mesh and storing it in array 'elmnts' (which is allocated here).
  ! 'remove_min_to_start_at_zero' is used to have the numbering of the nodes starting at '0'.
  ! 'nelmnts' is the number of elements, 'nnodes' is the number of nodes in the mesh.
  !-----------------------------------------------

  subroutine read_external_mesh_file(filename, remove_min_to_start_at_zero, NGNOD)

  use constants, only: MAX_STRING_LEN,IMAIN,myrank
  use part_unstruct_par, only: elmnts,nelmnts,nnodes

  implicit none

  character(len=MAX_STRING_LEN), intent(in)  :: filename
  integer, intent(out)  :: remove_min_to_start_at_zero
  integer, intent(in)  :: NGNOD

  integer  :: i,ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Reading data from external mesh file: ',trim(filename)
    write(IMAIN,*) '    NGNOD = ',NGNOD
    call flush_IMAIN()
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  ! note: adding access='stream' would further decrease file size
  open(unit=990, file=trim(filename), form='unformatted' , status='old', action='read',iostat=ier)
#else
  open(unit=990, file=trim(filename), form='formatted' , status='old', action='read',iostat=ier)
#endif
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    call stop_the_code('Error read external mesh file')
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  read(990) nelmnts
#else
  read(990,*) nelmnts
#endif

  allocate(elmnts(0:NGNOD*nelmnts-1),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating elmnts array')
  elmnts(:) = -1

  do i = 0, nelmnts-1
    if (NGNOD == 4) then
      ! linear elements
      ! 4 corner nodal points
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
      read(990) elmnts(i*NGNOD), elmnts(i*NGNOD+1), elmnts(i*NGNOD+2), elmnts(i*NGNOD+3)
#else
      read(990,*) elmnts(i*NGNOD), elmnts(i*NGNOD+1), elmnts(i*NGNOD+2), elmnts(i*NGNOD+3)
#endif
    else if (NGNOD == 9) then
      ! quadratic elements
      ! 4 corners + 4 edge mid-points + 1 center nodal point
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
      read(990) elmnts(i*NGNOD), elmnts(i*NGNOD+1), elmnts(i*NGNOD+2), elmnts(i*NGNOD+3), &
                  elmnts(i*NGNOD+4), elmnts(i*NGNOD+5), elmnts(i*NGNOD+6), elmnts(i*NGNOD+7), elmnts(i*NGNOD+8)
#else
      read(990,*) elmnts(i*NGNOD), elmnts(i*NGNOD+1), elmnts(i*NGNOD+2), elmnts(i*NGNOD+3), &
                  elmnts(i*NGNOD+4), elmnts(i*NGNOD+5), elmnts(i*NGNOD+6), elmnts(i*NGNOD+7), elmnts(i*NGNOD+8)
#endif
    else
      call stop_the_code('error, NGNOD should be either 4 or 9 for external meshes')
    endif
  enddo

  close(990)

  remove_min_to_start_at_zero = minval(elmnts(:))

  ! checks if we missed some elements
  if (remove_min_to_start_at_zero == -1) then
    print *,'Error: reading elements ',nelmnts,': mesh has missing nodal points. please check NGNOD setting...'
    call stop_the_code('Invalid mesh with missing nodal points')
  endif

  elmnts(:) = elmnts(:) - remove_min_to_start_at_zero

  nnodes = maxval(elmnts) + 1

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    Total number of spectral elements   :',nelmnts
    write(IMAIN,*) '    Total number of nodal points        :',nnodes
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_external_mesh_file

!
!---------------------------------------------------------------------------------------
!

  !-----------------------------------------------
  ! Read the material for each element and storing it in array 'num_materials'
  !-----------------------------------------------

  subroutine read_external_materials_file(filename)

  use constants, only: MAX_STRING_LEN,IMAIN,myrank

  use part_unstruct_par, only: nelmnts
  use shared_parameters, only: num_material

  implicit none

  character(len=MAX_STRING_LEN), intent(in)  :: filename

  ! local parameters
  integer  :: i,ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Reading materials from external mesh file: ',trim(filename)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! assigns materials to mesh elements
  allocate(num_material(nelmnts),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating num_material array')
  num_material(:) = 0

  ! file input
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  ! note: adding access='stream' would further decrease file size
  open(unit=992, file=trim(filename), form='unformatted' , status='old', action='read',iostat=ier)
#else
  open(unit=992, file=trim(filename), form='formatted' , status='old', action='read',iostat=ier)
#endif
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    call stop_the_code('Error read external mat file')
  endif

  do i = 1, nelmnts
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
     read(992) num_material(i)
#else
     read(992,*) num_material(i)
#endif
  enddo
  close(992)

  ! quick check
  if (any(num_material(:) == 0)) call stop_the_code('Error reading material array, some elements have zero material index')

  end subroutine read_external_materials_file

!
!---------------------------------------------------------------------------------------
!

  !-----------------------------------------------
  ! Read the PML elements, storing them in array 'region_pml_external_mesh'
  !-----------------------------------------------

  subroutine read_external_pml_element(filename, region_pml_external_mesh, nspec_cpml)

  use constants, only: MAX_STRING_LEN,CPML_X_ONLY,CPML_Z_ONLY,CPML_XZ,IMAIN,myrank
  use part_unstruct_par, only: nelmnts
  use compute_elements_load_par, only: is_pml

  implicit none

  character(len=MAX_STRING_LEN), intent(in)  :: filename
  integer, dimension(1:nelmnts), intent(out)  :: region_pml_external_mesh
  integer, intent(out)  :: nspec_cpml

  integer  :: i,ier,ispec,pml_flag

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Reading PML elements from external mesh file: ',trim(filename)
    call flush_IMAIN()
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  ! note: adding access='stream' would further decrease file size
  open(unit=992, file=trim(filename), form='unformatted' , status='old', action='read',iostat=ier)
#else
  open(unit=992, file=trim(filename), form='formatted' , status='old', action='read',iostat=ier)
#endif
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    call stop_the_code('Error read external absorbing_cpml_file')
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  read(992) nspec_cpml
#else
  read(992,*) nspec_cpml
#endif

  do i = 1, nspec_cpml
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
    read(992) ispec, pml_flag
#else
    read(992,*) ispec, pml_flag
#endif
    if (pml_flag /= CPML_X_ONLY .and. pml_flag /= CPML_Z_ONLY .and. pml_flag /= CPML_XZ) &
      call stop_the_code('error: incorrect CPML element flag found, should be CPML_X_ONLY or CPML_Z_ONLY or CPML_XZ only')

    ! stores element
    region_pml_external_mesh(ispec) = pml_flag
    is_pml(ispec-1) = .true.
  enddo

  close(992)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    Total number of PML elements: ',nspec_cpml
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_external_pml_element

!
!---------------------------------------------------------------------------------------
!

  !-----------------------------------------------
  ! read the node coordinates and store them in array 'nodes_coords'
  !-----------------------------------------------

  subroutine read_external_mesh_nodes_coords(filename)

  use constants, only: MAX_STRING_LEN,IMAIN,myrank
  use part_unstruct_par, only: nodes_coords,nnodes

  implicit none

  character(len=MAX_STRING_LEN), intent(in)  :: filename

  integer  :: i,ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Reading node coords from external mesh file: ',trim(filename)
    call flush_IMAIN()
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  ! note: adding access='stream' would further decrease file size
  open(unit=991, file=trim(filename), form='unformatted' , status='old', action='read', iostat=ier)
#else
  open(unit=991, file=trim(filename), form='formatted' , status='old', action='read', iostat=ier)
#endif
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    call stop_the_code('Error read external nodes coords file')
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  read(991) nnodes
#else
  read(991,*) nnodes
#endif

  allocate(nodes_coords(2,nnodes))
  nodes_coords(:,:) = 0.d0

  do i = 1, nnodes
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
     read(991) nodes_coords(1,i), nodes_coords(2,i)
#else
     read(991,*) nodes_coords(1,i), nodes_coords(2,i)
#endif
  enddo
  close(991)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    node coords: X min/max = ',minval(nodes_coords(1,:)),"/",maxval(nodes_coords(1,:))
    write(IMAIN,*) '                 Z min/max = ',minval(nodes_coords(2,:)),"/",maxval(nodes_coords(2,:))
    write(IMAIN,*)
    call flush_IMAIN()
  endif


  end subroutine read_external_mesh_nodes_coords

!
!---------------------------------------------------------------------------------------
!

  !-----------------------------------------------
  ! Read free surface.
  ! Edges from elastic elements are discarded.
  ! 'acoustic_surface' contains 1/ element number, 2/ number of nodes that form the free surface,
  ! 3/ first node on the free surface, 4/ second node on the free surface, if relevant (if 2/ is equal to 2)
  !-----------------------------------------------

  subroutine read_external_acoustic_surface(filename, num_material, &
                                            nbmodels, icodemat, phi_material, remove_min_to_start_at_zero)

  use constants, only: MAX_STRING_LEN,ANISOTROPIC_MATERIAL,IMAIN,myrank
  use part_unstruct_par, only: nelmnts,nelem_acoustic_surface,acoustic_surface

  implicit none

  character(len=MAX_STRING_LEN), intent(in)  :: filename
  integer, dimension(0:nelmnts-1)  :: num_material
  integer, intent(in)  :: nbmodels
  integer, dimension(1:nbmodels), intent(in)  :: icodemat
  double precision, dimension(1:nbmodels), intent(in)  :: phi_material
  integer, intent(in)  :: remove_min_to_start_at_zero

  ! local parameters
  integer, dimension(:,:), allocatable  :: acoustic_surface_tmp
  integer  :: nelmnts_surface
  integer  :: i,ier
  integer  :: imaterial_number

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Reading acoustic surface from external mesh file: ',trim(filename)
    call flush_IMAIN()
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  ! note: adding access='stream' would further decrease file size
  open(unit=993, file=trim(filename), form='unformatted' , status='old', action='read', iostat=ier)
#else
  open(unit=993, file=trim(filename), form='formatted' , status='old', action='read', iostat=ier)
#endif
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    call stop_the_code('Error read acoustic surface file')
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  read(993) nelmnts_surface
#else
  read(993,*) nelmnts_surface
#endif

  allocate(acoustic_surface_tmp(4,nelmnts_surface))

  do i = 1, nelmnts_surface
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
     read(993) acoustic_surface_tmp(1,i), acoustic_surface_tmp(2,i), acoustic_surface_tmp(3,i), acoustic_surface_tmp(4,i)
#else
     read(993,*) acoustic_surface_tmp(1,i), acoustic_surface_tmp(2,i), acoustic_surface_tmp(3,i), acoustic_surface_tmp(4,i)
#endif
  enddo

  close(993)

  acoustic_surface_tmp(1,:) = acoustic_surface_tmp(1,:) - remove_min_to_start_at_zero
  acoustic_surface_tmp(3,:) = acoustic_surface_tmp(3,:) - remove_min_to_start_at_zero
  acoustic_surface_tmp(4,:) = acoustic_surface_tmp(4,:) - remove_min_to_start_at_zero

  nelem_acoustic_surface = 0
  do i = 1, nelmnts_surface
     imaterial_number = num_material(acoustic_surface_tmp(1,i))
     if (icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. phi_material(imaterial_number) >= 1.d0) then
        nelem_acoustic_surface = nelem_acoustic_surface + 1
     endif
  enddo

  allocate(acoustic_surface(4,nelem_acoustic_surface),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating acoustic_surface array')

  nelem_acoustic_surface = 0
  do i = 1, nelmnts_surface
     imaterial_number = num_material(acoustic_surface_tmp(1,i))
     if (icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. phi_material(imaterial_number) >= 1.d0) then
        nelem_acoustic_surface = nelem_acoustic_surface + 1
        acoustic_surface(:,nelem_acoustic_surface) = acoustic_surface_tmp(:,i)
     endif
  enddo

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Total number of surface elements         : ',nelmnts_surface
    write(IMAIN,*) '  Total number of acoustic surface elements: ',nelem_acoustic_surface
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_external_acoustic_surface

!
!---------------------------------------------------------------------------------------
!

  !-----------------------------------------------
  ! Read absorbing surface.
  ! 'abs_surface' contains 1/ element number, 2/ number of nodes that form the absorbing edge
  ! (which currently must always be equal to 2),
  ! 3/ first node on the abs surface, 4/ second node on the abs surface
  ! 5/ 1=IBOTTOM, 2=IRIGHT, 3=ITOP, 4=ILEFT
  !-----------------------------------------------
  subroutine read_external_abs_surface(filename, remove_min_to_start_at_zero)

  use constants, only: MAX_STRING_LEN,IMAIN,myrank
  use part_unstruct_par, only: abs_surface,nelemabs

  implicit none

  character(len=MAX_STRING_LEN), intent(in)  :: filename
  integer, intent(in)  :: remove_min_to_start_at_zero

  integer  :: i,ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Reading absorbing surface from external mesh file: ',trim(filename)
    call flush_IMAIN()
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  ! note: adding access='stream' would further decrease file size
  open(unit=994, file=trim(filename), form='unformatted' , status='old', action='read', iostat=ier)
#else
  open(unit=994, file=trim(filename), form='formatted' , status='old', action='read', iostat=ier)
#endif
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    call stop_the_code('Error read absorbing surface file')
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  read(994) nelemabs
#else
  read(994,*) nelemabs
#endif

  allocate(abs_surface(5,nelemabs))

  do i = 1, nelemabs
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
    read(994) abs_surface(1,i), abs_surface(2,i), abs_surface(3,i), abs_surface(4,i), abs_surface(5,i)
#else
    read(994,*) abs_surface(1,i), abs_surface(2,i), abs_surface(3,i), abs_surface(4,i), abs_surface(5,i)
#endif

    if (abs_surface(2,i) /= 2) then
      print *,'Only two nodes per absorbing element can be listed.'
      print *,'If one of your elements has more than one edge along a given absorbing contour'
      print *,'(e.g., if that contour has a corner) then list it twice,'
      print *,'putting the first edge on the first line and the second edge on the second line.'
      print *,'If one of your elements has a single point along the absording contour rather than a full edge, do NOT list it'
      print *,'(it would have no weight in the contour integral anyway because it would consist of a single point).'
      print *,'If you use 9-node elements, list only the first and last points of the edge and not the intermediate point'
      print *,'located around the middle of the edge; the right 9-node curvature will be restored automatically by the code.'

      call stop_the_code('only two nodes per element should be listed for absorbing edges')
    endif

    if (abs_surface(5,i) < 1 .or. abs_surface(5,i) > 4) then
      call stop_the_code('absorbing element type must be between 1 (IBOTTOM) and 4 (ILEFT)')
    endif

  enddo

  close(994)

  abs_surface(1,:) = abs_surface(1,:) - remove_min_to_start_at_zero
  abs_surface(3,:) = abs_surface(3,:) - remove_min_to_start_at_zero
  abs_surface(4,:) = abs_surface(4,:) - remove_min_to_start_at_zero

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    Total number of absorbing surface elements: ',nelemabs
    write(IMAIN,*)
    call flush_IMAIN()
  endif


  end subroutine read_external_abs_surface

!
!---------------------------------------------------------------------------------------
!

  !-----------------------------------------------
  ! Read acoustic forcing surface.
  ! 'acforcing_surface' contains 1/ element number, 2/ number of nodes that form the acoustic forcing edge
  ! (which currently must always be equal to 2),
  ! 3/ first node on the acforcing surface, 4/ second node on the acforcing surface
  ! 5/ 1=IBOTTOME, 2=IRIGHT, 3=ITOP, 4=ILEFT
  !-----------------------------------------------

  subroutine read_external_acoustic_forcing_surface(filename, remove_min_to_start_at_zero)

  use constants, only: MAX_STRING_LEN,IMAIN,myrank
  use part_unstruct_par, only: acforcing_surface,nelemacforcing

  implicit none

  character(len=MAX_STRING_LEN), intent(in)  :: filename
  integer, intent(in)  :: remove_min_to_start_at_zero

  integer  :: i,ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Reading acoustic forcing surface from external mesh file: ',trim(filename)
    call flush_IMAIN()
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  ! note: adding access='stream' would further decrease file size
  open(unit=995, file=trim(filename), form='unformatted' , status='old', action='read', iostat=ier)
#else
  open(unit=995, file=trim(filename), form='formatted' , status='old', action='read', iostat=ier)
#endif
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    call stop_the_code('Error read acoustic forcing surface file')
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  read(995) nelemacforcing
#else
  read(995,*) nelemacforcing
#endif

  allocate(acforcing_surface(5,nelemacforcing))

  do i = 1, nelemacforcing
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
    read(995) acforcing_surface(1,i), acforcing_surface(2,i), acforcing_surface(3,i), acforcing_surface(4,i), &
                acforcing_surface(5,i)
#else
    read(995,*) acforcing_surface(1,i), acforcing_surface(2,i), acforcing_surface(3,i), acforcing_surface(4,i), &
                acforcing_surface(5,i)
#endif

    if (acforcing_surface(2,i) /= 2) then
      print *,'Only two nodes per acoustic forcing element can be listed.'
      print *,'If one of your elements has more than one edge along a given acoustic forcing contour'
      print *,'(e.g., if that contour has a corner) then list it twice,'
      print *,'putting the first edge on the first line and the second edge on the second line.'
      print *,'If one of your elements has a single point along the acoustic forcing contour rather than a full edge, do NOT'
      print *,'list it (it would have no weight in the contour integral anyway because it would consist of a single point).'
      print *,'If you use 9-node elements, list only the first and last points of the edge and not the intermediate point'
      print *,'located around the middle of the edge; the right 9-node curvature will be restored automatically by the code.'

      call stop_the_code('only two nodes per element should be listed for absorbing edges')
    endif

    if (acforcing_surface(5,i) < 1 .or. acforcing_surface(5,i) > 4) then
      call stop_the_code('absorbing element type must be between 1 (IBOTTOM) and 4 (ILEFT)')
    endif

  enddo

  close(995)

  acforcing_surface(1,:) = acforcing_surface(1,:) - remove_min_to_start_at_zero
  acforcing_surface(3,:) = acforcing_surface(3,:) - remove_min_to_start_at_zero
  acforcing_surface(4,:) = acforcing_surface(4,:) - remove_min_to_start_at_zero

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    Total number of acoustic forcing surface elements: ',nelemacforcing
    write(IMAIN,*)
    call flush_IMAIN()
  endif


  end subroutine read_external_acoustic_forcing_surface

!
!---------------------------------------------------------------------------------------
!

  !-----------------------------------------------
  ! Read axial elements file.
  ! The first line of this file must be the total number of axial elements then it contains
  ! the ispec of each of these elements.
  ! 'axial_elements' contains the list of the ispec corresponding to axial elements
  !-----------------------------------------------

  subroutine read_external_axial_elements_file(filename,remove_min_to_start_at_zero)

  use constants, only: MAX_STRING_LEN,IMAIN,myrank
  use part_unstruct_par, only: ispec_of_axial_elements,nelem_on_the_axis,inode1_axial_elements,inode2_axial_elements

  implicit none

  character(len=MAX_STRING_LEN), intent(in)  :: filename
  integer, intent(in)  :: remove_min_to_start_at_zero

  integer :: i,j,ier
  integer :: dump

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Reading axial elements from external mesh file: ',trim(filename)
    call flush_IMAIN()
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  ! note: adding access='stream' would further decrease file size
  open(unit=994, file=trim(filename), form='unformatted' , status='old', action='read', iostat=ier)
#else
  open(unit=994, file=trim(filename), form='formatted' , status='old', action='read', iostat=ier)
#endif
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    call stop_the_code('Error read axial elements file')
  endif

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  read(994) nelem_on_the_axis
#else
  read(994,*) nelem_on_the_axis
#endif

  ! user output
  write(IMAIN,*) "Number of elements on the axis: ", nelem_on_the_axis

  allocate(ispec_of_axial_elements(nelem_on_the_axis),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating array ispec_of_axial_elements')

  ! needed for rotation
  allocate(inode1_axial_elements(nelem_on_the_axis), &
           inode2_axial_elements(nelem_on_the_axis),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating array inode**_axial_elements')

  do i = 1, nelem_on_the_axis ! Dump is always 2 (old convention from absorbing surfaces)
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
    read(994) ispec_of_axial_elements(i), dump, inode1_axial_elements(i), inode2_axial_elements(i)
#else
    read(994,*) ispec_of_axial_elements(i), dump, inode1_axial_elements(i), inode2_axial_elements(i)
#endif

  enddo

  close(994)

  ! this has a cost of O(nelem_on_the_axis^2), could be reduced by using a quicksort algorithm
  ! and checking the sorted list
  write(IMAIN,*) 'testing for duplicates in the axial element input file...'
  do i = 1,nelem_on_the_axis
    do j = i+1,nelem_on_the_axis
      if (ispec_of_axial_elements(i) == ispec_of_axial_elements(j)) then
        call stop_the_code('At least one element appears twice in the axial element file')
      endif
    enddo
  enddo
  write(IMAIN,*) 'done testing for duplicates'
  write(IMAIN,*)

  ! axisym TODO : Test if the informations supplied are compatible with axisym

  ispec_of_axial_elements(:) = ispec_of_axial_elements(:) - remove_min_to_start_at_zero
  inode1_axial_elements(:) = inode1_axial_elements(:) - remove_min_to_start_at_zero
  inode2_axial_elements(:) = inode2_axial_elements(:) - remove_min_to_start_at_zero

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    Total number of axial elements: ',nelem_on_the_axis
    write(IMAIN,*)
    call flush_IMAIN()
  endif


  end subroutine read_external_axial_elements_file

!
!---------------------------------------------------------------------------------------
!

  subroutine read_external_tangential_curve_file(filename)

! reads in tangential detection curve file

  use constants, only: MAX_STRING_LEN,IIN,IMAIN,myrank

  use part_unstruct_par, only: nnodes_tangential_curve,nodes_tangential_curve

  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: filename

  ! local parameters
  integer :: i,ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Reading tangential curve from external mesh file: ',trim(filename)
    call flush_IMAIN()
  endif

  ! reads in specified external file
  open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    call stop_the_code('Error read tangential curve file')
  endif

  read(IIN,*) nnodes_tangential_curve

  allocate(nodes_tangential_curve(2,nnodes_tangential_curve),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating tangential array')

  do i = 1, nnodes_tangential_curve
    read(IIN,*) nodes_tangential_curve(1,i), nodes_tangential_curve(2,i)
  enddo
  close(IIN)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    Total number of tangential curve nodes: ',nnodes_tangential_curve
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_external_tangential_curve_file

