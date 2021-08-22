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

  subroutine read_interfaces_file()

  use constants, only: IMAIN,IIN_INTERFACES,DONT_IGNORE_JUNK,HUGEVAL,MAX_STRING_LEN,IN_DATA_FILES

  use shared_parameters, only: interfacesfile,nx_param, &
    nz_layer,number_of_layers, &
    max_npoints_interface,number_of_interfaces,npoints_of_interfaces, &
    xinterface_coords,zinterface_coords, &
    nxread,nzread

  ! see comment below for uncommenting this feature
  !use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS
  !use constants, only: mygroup
  implicit none

  ! local parameters
  integer :: ier,interface_current,ipoint_current,ilayer,idummy
  integer :: npoints_interface
  double precision :: xval,zval,xval_previous
  character(len=MAX_STRING_LEN) :: interfaces_filename
  !character(len=MAX_STRING_LEN) :: path_to_add

  interfaces_filename = trim(IN_DATA_FILES)//trim(interfacesfile) ! by default: DATA/..

! in case we want to allow for different mesh interfaces between different run setups
!
! note: this is a somewhat dangerous feature and differs from 3D version behavior.
!       simultaneous runs had the original idea of using the same mesh, but allow for different source/station setups.
!       thus, saving time in having to run the mesher multiple times for multiple runs.
!
!       allowing to read in different interfaces here defeats this original idea for simultaneous runs
!       we'll disable it for now and require the mesh setup (interface,tomography, external mesh files)
!       be read from a main reference folder DATA/, not run****/DATA/..
!
! please uncomment in case needed and make sure you know how to run your simultaneous runs properly :)

  !if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
  !  write(path_to_add,"('run',i4.4,'/')") mygroup + 1
  !  interfaces_filename = path_to_add(1:len_trim(path_to_add))//interfaces_filename(1:len_trim(interfaces_filename))
  !endif

  ! user output
  write(IMAIN,*) 'Reading interface data from file: '//trim(interfaces_filename)
  call flush_IMAIN()

  ! get interface data from external file to count the spectral elements along Z
  open(unit=IIN_INTERFACES,file=trim(interfaces_filename),status='old',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: '//trim(interfaces_filename)
    call exit_MPI(0,'Error read interface file in meshfem2D')
  endif

  max_npoints_interface = -1

  ! read number of interfaces
  call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,number_of_interfaces)

  ! safety check
  if (number_of_interfaces < 2) call stop_the_code('not enough interfaces (minimum is 2)')

  ! loop on all the interfaces
  do interface_current = 1,number_of_interfaces

    call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface)

    ! check
    if (npoints_interface < 2) call stop_the_code('not enough interface points (minimum is 2)')

    max_npoints_interface = max(npoints_interface,max_npoints_interface)

    ! user output
    write(IMAIN,*) 'Reading ',npoints_interface,' points for interface ',interface_current
    call flush_IMAIN()

    ! loop on all the points describing this interface
    xval_previous = -HUGEVAL

    do ipoint_current = 1,npoints_interface

      call read_two_interface_points(IIN_INTERFACES,DONT_IGNORE_JUNK,xval,zval)

      ! safety check
      if (ipoint_current > 1 .and. xval <= xval_previous) &
        call stop_the_code('interface points must be sorted in increasing X')

      xval_previous = xval
    enddo
  enddo

  ! allocates interface points for mesher
  allocate(xinterface_coords(max_npoints_interface,number_of_interfaces), &
           zinterface_coords(max_npoints_interface,number_of_interfaces),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating interface arrays')

  xinterface_coords(:,:) = 0.0
  zinterface_coords(:,:) = 0.0

  allocate(npoints_of_interfaces(number_of_interfaces),stat=ier)
  if (ier /= 0) call exit_MPI(0,'Error allocating npoints_interfaces array')
  npoints_of_interfaces(:) = 0

  ! start reading from beginning to store all interface points
  !
  ! note: this will avoid reading this same file again later in the mesher.
  !       reading all interface information is done here once and for all.
  rewind(IIN_INTERFACES)

  ! skip
  call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,idummy)  ! number_of_interfaces

  ! loop on all the interfaces
  do interface_current = 1,number_of_interfaces
    ! reads number of interface points
    call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface)

    ! check
    if (npoints_interface < 2) call stop_the_code('not enough interface points (minimum is 2) to store values')

    ! stores number of points
    npoints_of_interfaces(interface_current) = npoints_interface

    ! loop on all the points describing this interface
    xval_previous = -HUGEVAL
    do ipoint_current = 1,npoints_interface
      call read_two_interface_points(IIN_INTERFACES,DONT_IGNORE_JUNK,xval,zval)

      ! safety check
      if (ipoint_current > 1 .and. xval <= xval_previous) &
        call stop_the_code('interface points must be sorted in increasing X to store values')

      ! stores interface point
      xinterface_coords(ipoint_current,interface_current) = xval
      zinterface_coords(ipoint_current,interface_current) = zval
    enddo
  enddo

  ! define number of layers
  number_of_layers = number_of_interfaces - 1

  allocate(nz_layer(number_of_layers),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating array nz_layer')
  nz_layer(:) = 0

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'Total number of layers in z direction = ', number_of_layers
  call flush_IMAIN()

  ! loop on all the layers
  do ilayer = 1,number_of_layers

    ! read number of spectral elements in vertical direction in this layer
    call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,nz_layer(ilayer))

    if (nz_layer(ilayer) < 1) call stop_the_code('not enough spectral elements along Z in layer (minimum is 1)')

    ! user output
    write(IMAIN,*) 'There are ',nz_layer(ilayer),' spectral elements along Z in layer ',ilayer
    call flush_IMAIN()
  enddo

  ! all done reading
  close(IIN_INTERFACES)

  ! compute total number of spectral elements in vertical direction
  nzread = sum(nz_layer)

  ! sets number of elements along X for internal mesher
  nxread = nx_param

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'Total number of spectral elements along X = ',nxread
  write(IMAIN,*) 'Total number of spectral elements along Z = ',nzread
  write(IMAIN,*)
  call flush_IMAIN()

  end subroutine read_interfaces_file
