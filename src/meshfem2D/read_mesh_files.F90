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

  subroutine read_mesh_tangential_curve_file()

! reads in tangential detection

  use part_unstruct_par, only: nnodes_tangential_curve

  use shared_parameters, only: force_normal_to_surface,rec_normal_to_surface,read_external_mesh

  implicit none

  ! initializes
  nnodes_tangential_curve = 0

  if (force_normal_to_surface .or. rec_normal_to_surface) then
    ! reads in file
    if (read_external_mesh) then
      ! reads in specified external file
      call read_external_tangential_curve_file()
    else
      ! safety stop
      call stop_the_code('Error read_external_mesh must be set to .true. to use external tangential_dectection_curve_file')
    endif
  endif

  end subroutine read_mesh_tangential_curve_file

!
!---------------------------------------------------------------------------------------
!

  subroutine read_mesh_nodes_coords_from_interfaces()

  use constants, only: IMAIN,IIN_INTERFACES,DONT_IGNORE_JUNK,MAX_STRING_LEN,mygroup,IN_DATA_FILES, &
                       ADD_RANDOM_PERTURBATION_TO_THE_MESH,RANDOM_AMPLITUDE_ALONG_X,RANDOM_AMPLITUDE_ALONG_Z, &
                       ADD_PERTURBATION_AROUND_SOURCE_ONLY,RADIUS_TO_USE_AROUND_THE_SOURCE

  use part_unstruct_par, only: nodes_coords,max_npoints_interface,number_of_interfaces, &
    npoints_interface_top,xinterface_top,zinterface_top,coefs_interface_top, &
    nx,nz,nxread,nzread,nz_layer,number_of_layers,nnodes,grid_point_x,grid_point_z

  use source_file_par, only: source_surf,xs,zs

  use shared_parameters, only: ngnod,interfacesfile,NSOURCES,xmin_param,xmax_param,NUMBER_OF_SIMULTANEOUS_RUNS, &
                                     PML_BOUNDARY_CONDITIONS,NELEM_PML_THICKNESS

  implicit none

  ! local parameters
  integer :: ilayer,ipoint_current
  integer  :: num_node

  integer :: npoints_interface_bottom
  double precision, dimension(:), allocatable :: xinterface_bottom,zinterface_bottom,coefs_interface_bottom

  ! to compute the coordinate transformation
  integer :: ioffset
  double precision :: gamma,absx,a00,a01,bot0,top0,random_value,radius_squared
  double precision :: tang1,tangN

  integer :: i_source
  integer :: ix,iz
  integer :: i,j

  logical :: need_to_add_the_perturbation

  ! external functions
  integer, external :: num_4, num_9
  double precision, external :: value_spline
  character(len=MAX_STRING_LEN) :: interfaces_filename,path_to_add

  interfaces_filename = trim(IN_DATA_FILES)//trim(interfacesfile)

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    interfaces_filename = path_to_add(1:len_trim(path_to_add))//interfaces_filename(1:len_trim(interfaces_filename))
  endif

  ! get interface data from external file
  write(IMAIN,*) 'Reading interface data from file: '//trim(IN_DATA_FILES)//interfacesfile(1:len_trim(interfacesfile))
  open(unit=IIN_INTERFACES,file=trim(interfaces_filename),status='old')

  ! allocate arrays for the grid
  allocate(grid_point_x(0:nx,0:nz))
  allocate(grid_point_z(0:nx,0:nz))

  grid_point_x(:,:) = 0.d0
  grid_point_z(:,:) = 0.d0

  ! top interface arrays (required also for station locations)
  allocate(xinterface_top(max_npoints_interface))
  allocate(zinterface_top(max_npoints_interface))
  allocate(coefs_interface_top(max_npoints_interface))

  ! temporary arrays
  allocate(xinterface_bottom(max_npoints_interface))
  allocate(zinterface_bottom(max_npoints_interface))
  allocate(coefs_interface_bottom(max_npoints_interface))

  ! read number of interfaces
  call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,number_of_interfaces)

  ! read bottom interface
  call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface_bottom)

  ! loop on all the points describing this interface
  do ipoint_current = 1,npoints_interface_bottom
    call read_two_interface_points(IIN_INTERFACES,DONT_IGNORE_JUNK, &
                                   xinterface_bottom(ipoint_current),zinterface_bottom(ipoint_current))
  enddo

  ! loop on all the layers
  do ilayer = 1,number_of_layers

    ! read top interface
    call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface_top)

    ! loop on all the points describing this interface
    do ipoint_current = 1,npoints_interface_top
       call read_two_interface_points(IIN_INTERFACES,DONT_IGNORE_JUNK, &
                                      xinterface_top(ipoint_current),zinterface_top(ipoint_current))
    enddo

    ! compute the spline for the bottom interface, impose the tangent on both edges
    tang1 = (zinterface_bottom(2)-zinterface_bottom(1)) / (xinterface_bottom(2)-xinterface_bottom(1))
    tangN = (zinterface_bottom(npoints_interface_bottom)-zinterface_bottom(npoints_interface_bottom-1)) / &
         (xinterface_bottom(npoints_interface_bottom)-xinterface_bottom(npoints_interface_bottom-1))

    call spline_construction(xinterface_bottom,zinterface_bottom,npoints_interface_bottom, &
                             tang1,tangN,coefs_interface_bottom)

    ! compute the spline for the top interface, impose the tangent on both edges
    tang1 = (zinterface_top(2)-zinterface_top(1)) / (xinterface_top(2)-xinterface_top(1))
    tangN = (zinterface_top(npoints_interface_top)-zinterface_top(npoints_interface_top-1)) / &
         (xinterface_top(npoints_interface_top)-xinterface_top(npoints_interface_top-1))

    call spline_construction(xinterface_top,zinterface_top,npoints_interface_top, &
                             tang1,tangN,coefs_interface_top)

    ! check if we are in the last layer, which contains topography,
    ! and modify the position of the source accordingly if it is located exactly at the surface
    do i_source= 1,NSOURCES
       if (source_surf(i_source) .and. ilayer == number_of_layers) then
          ! user output
          write(IMAIN,*) 'source ', i_source
          write(IMAIN,*) '  target (input) z: ', zs(i_source)

          zs(i_source) = value_spline(xs(i_source),xinterface_top,zinterface_top, &
                                      coefs_interface_top,npoints_interface_top)

          write(IMAIN,*) '  surface (actual) z: ', zs(i_source)
       endif
    enddo

    ! compute the offset of this layer in terms of number of spectral elements below along Z
    if (ilayer > 1) then
       ioffset = sum(nz_layer(1:ilayer-1))
    else
       ioffset = 0
    endif

    !--- definition of the mesh

    do ix = 0,nx

       ! evenly spaced points along X
       absx = xmin_param + (xmax_param - xmin_param) * dble(ix) / dble(nx)

       ! value of the bottom and top splines
       bot0 = value_spline(absx,xinterface_bottom,zinterface_bottom,coefs_interface_bottom,npoints_interface_bottom)
       top0 = value_spline(absx,xinterface_top,zinterface_top,coefs_interface_top,npoints_interface_top)

       do iz = 0,nz_layer(ilayer)

          ! linear interpolation between bottom and top
          gamma = dble(iz) / dble(nz_layer(ilayer))
          a00 = 1.d0 - gamma
          a01 = gamma

          ! coordinates of the grid points
          grid_point_x(ix,iz + ioffset) = absx
          grid_point_z(ix,iz + ioffset) = a00*bot0 + a01*top0

          ! add a random perturbation to the mesh if needed
          if (ADD_RANDOM_PERTURBATION_TO_THE_MESH) then

            need_to_add_the_perturbation = .true.

            ! do not make any modification on the outer edges of the mesh (fictitious outer edges, and free surface)
            if ((ilayer == 1 .and. iz == 0) .or. (ilayer == number_of_layers .and. iz == nz_layer(ilayer)) &
                           .or. ix == 0 .or. ix == nx) need_to_add_the_perturbation = .false.

            ! do not make any modification inside the PML layers
            if (PML_BOUNDARY_CONDITIONS) then
              ! this works only if the whole bottom PML is comprised inside the bottom layer
              if (ilayer == 1 .and. iz <= NELEM_PML_THICKNESS) need_to_add_the_perturbation = .false.
              if (ix <= NELEM_PML_THICKNESS) need_to_add_the_perturbation = .false.
              if (ix > nx - NELEM_PML_THICKNESS - 1) need_to_add_the_perturbation = .false.
            endif

            ! apply the perturbation in a disk around the source only
            ! for simplicity here we assume that there is a single source
            if (ADD_PERTURBATION_AROUND_SOURCE_ONLY) then
              radius_squared = (grid_point_x(ix,iz + ioffset) - xs(1))**2 + (grid_point_z(ix,iz + ioffset) - zs(1))**2
              if (radius_squared > RADIUS_TO_USE_AROUND_THE_SOURCE**2) need_to_add_the_perturbation = .false.
            endif

            if (need_to_add_the_perturbation) then
              ! get a random number between 0. and 1.
              call random_number(random_value)
              ! map this random number to between -1. and +1., so that the point can be moved to the left or to the right
              grid_point_x(ix,iz + ioffset) = grid_point_x(ix,iz + ioffset) + &
                                                  RANDOM_AMPLITUDE_ALONG_X * 2.d0 * (random_value - 0.5d0)
              call random_number(random_value)
              grid_point_z(ix,iz + ioffset) = grid_point_z(ix,iz + ioffset) + &
                                                  RANDOM_AMPLITUDE_ALONG_Z * 2.d0 * (random_value - 0.5d0)
            endif

          endif

       enddo

    enddo

    ! the top interface becomes the bottom interface before switching to the next layer
    npoints_interface_bottom = npoints_interface_top
    xinterface_bottom(:) = xinterface_top(:)
    zinterface_bottom(:) = zinterface_top(:)

  enddo

  close(IIN_INTERFACES)

  deallocate(xinterface_bottom,zinterface_bottom,coefs_interface_bottom)

  ! sets node coordinates
  nnodes = (nz+1)*(nx+1)

  allocate(nodes_coords(2,nnodes))
  nodes_coords(:,:) = 0.d0

  if (ngnod == 4) then
    do j = 0, nz
      do i = 0, nx
        num_node = num_4(i,j,nxread)
        nodes_coords(1, num_node) = grid_point_x(i,j)
        nodes_coords(2, num_node) = grid_point_z(i,j)
      enddo
    enddo
  else if (ngnod == 9) then
    do j = 0, nz
      do i = 0, nx
        num_node = num_9(i,j,nxread,nzread)
        nodes_coords(1, num_node) = grid_point_x(i,j)
        nodes_coords(2, num_node) = grid_point_z(i,j)
      enddo
    enddo
  else
    call stop_the_code('ngnod should be either 4 or 9')
  endif

  end subroutine read_mesh_nodes_coords_from_interfaces





