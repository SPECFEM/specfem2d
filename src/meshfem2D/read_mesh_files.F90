!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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
      stop 'Error read_external_mesh must be set to .true. to use external tangential_dectection_curve_file'
    endif
  endif

  end subroutine read_mesh_tangential_curve_file

!
!---------------------------------------------------------------------------------------
!

  subroutine read_mesh_nodes_coords_from_interfaces()

  use constants, only: IMAIN,IIN_INTERFACES,DONT_IGNORE_JUNK

  use part_unstruct_par, only: nodes_coords,max_npoints_interface,number_of_interfaces, &
    npoints_interface_top,xinterface_top,zinterface_top,coefs_interface_top, &
    nx,nz,nxread,nzread,nz_layer,number_of_layers,nnodes,grid_point_x,grid_point_z

  use source_file_par, only: source_surf,xs,zs

  use shared_parameters, only: ngnod,interfacesfile,NSOURCES,xmin_param,xmax_param

  implicit none

  ! local parameters
  integer :: ilayer,ipoint_current
  integer  :: num_node

  integer :: npoints_interface_bottom
  double precision, dimension(:), allocatable :: xinterface_bottom,zinterface_bottom,coefs_interface_bottom

  ! to compute the coordinate transformation
  integer :: ioffset
  double precision :: gamma,absx,a00,a01,bot0,top0
  double precision :: tang1,tangN

  integer :: i_source
  integer :: ix,iz
  integer :: i,j

  ! external functions
  integer, external :: num_4, num_9
  double precision, external :: value_spline

  ! get interface data from external file
  write(IMAIN,*) 'Reading interface data from file: ','DATA/' // interfacesfile(1:len_trim(interfacesfile))
  open(unit=IIN_INTERFACES,file='DATA/'//interfacesfile,status='old')

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
    stop 'ngnod should be either 4 or 9'
  endif

  end subroutine read_mesh_nodes_coords_from_interfaces





