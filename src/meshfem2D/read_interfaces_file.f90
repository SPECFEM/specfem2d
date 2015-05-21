
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
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
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

module interfaces_file

  ! note: we use this module definition only to be able to allocate
  !          arrays for receiverlines and materials in this subroutine rather than in the main
  !          routine in meshfem2D.F90

  ! note 2: the filename ending is .F90 to have pre-compilation with pragmas
  !            (like #ifndef USE_MPI) working properly

  implicit none

contains

  subroutine read_interfaces_file(interfacesfile,max_npoints_interface, &
                                number_of_interfaces,npoints_interface_bottom, &
                                number_of_layers,nz_layer,nx,nz,nxread,nzread,ngnod, &
                                nelmnts,elmnts)
  implicit none
  include "constants.h"

  character(len=100) :: interfacesfile

  integer :: max_npoints_interface,number_of_interfaces,npoints_interface_bottom
  integer :: number_of_layers,nx,nz,nxread,nzread,ngnod
  integer :: nelmnts
  integer, dimension(:), pointer :: nz_layer
  integer, dimension(:), pointer  :: elmnts

  ! local parameters
  integer :: ios,interface_current,ipoint_current,ilayer,i,j,num_elmnt
  double precision :: xinterface_dummy,zinterface_dummy,xinterface_dummy_previous

  ! get interface data from external file to count the spectral elements along Z
  print *,'Reading interface data from file DATA/',interfacesfile(1:len_trim(interfacesfile)),' to count the spectral elements'
  open(unit=IIN_INTERFACES,file='DATA/'//interfacesfile,status='old',iostat=ios)
  if( ios /= 0 ) then
    print *,'error opening file: ',trim('DATA/'//interfacesfile)
    stop 'error read interface file in meshfem2D'
  endif

  max_npoints_interface = -1

  ! read number of interfaces
  call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,number_of_interfaces)
  if(number_of_interfaces < 2) stop 'not enough interfaces (minimum is 2)'

  ! loop on all the interfaces
  do interface_current = 1,number_of_interfaces

    call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface_bottom)
    if(npoints_interface_bottom < 2) stop 'not enough interface points (minimum is 2)'
    max_npoints_interface = max(npoints_interface_bottom,max_npoints_interface)
    print *,'Reading ',npoints_interface_bottom,' points for interface ',interface_current

    ! loop on all the points describing this interface
    xinterface_dummy_previous = -HUGEVAL
    do ipoint_current = 1,npoints_interface_bottom
       call read_two_interface_points(IIN_INTERFACES,DONT_IGNORE_JUNK,xinterface_dummy,zinterface_dummy)
       if(ipoint_current > 1 .and. xinterface_dummy <= xinterface_dummy_previous) &
            stop 'interface points must be sorted in increasing X'
       xinterface_dummy_previous = xinterface_dummy
    enddo
  enddo

  ! define number of layers
  number_of_layers = number_of_interfaces - 1

  allocate(nz_layer(number_of_layers))

  print *, 'Total number of layers in z direction = ', number_of_layers

  ! loop on all the layers
  do ilayer = 1,number_of_layers

    ! read number of spectral elements in vertical direction in this layer
    call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,nz_layer(ilayer))
    if(nz_layer(ilayer) < 1) stop 'not enough spectral elements along Z in layer (minimum is 1)'
    print *,'There are ',nz_layer(ilayer),' spectral elements along Z in layer ',ilayer

  enddo

  close(IIN_INTERFACES)

  ! compute total number of spectral elements in vertical direction
  nz = sum(nz_layer)

  print *
  print *,'Total number of spectral elements along Z = ',nz
  print *

  nxread = nx
  nzread = nz

  ! multiply by 2 if elements have 9 nodes
  if(ngnod == 9) then
    nx = nx * 2
    nz = nz * 2
    nz_layer(:) = nz_layer(:) * 2
  endif

  nelmnts = nxread * nzread
  allocate(elmnts(0:ngnod*nelmnts-1))

  if ( ngnod == 4 ) then
    num_elmnt = 0
    do j = 1, nzread
       do i = 1, nxread
          elmnts(num_elmnt*ngnod)   = (j-1)*(nxread+1) + (i-1)
          elmnts(num_elmnt*ngnod+1) = (j-1)*(nxread+1) + (i-1) + 1
          elmnts(num_elmnt*ngnod+2) = j*(nxread+1) + (i-1) + 1
          elmnts(num_elmnt*ngnod+3) = j*(nxread+1) + (i-1)
          num_elmnt = num_elmnt + 1
       enddo
    enddo
  else
    num_elmnt = 0
    do j = 1, nzread
       do i = 1, nxread
          elmnts(num_elmnt*ngnod)   = (j-1)*(nxread+1) + (i-1)
          elmnts(num_elmnt*ngnod+1) = (j-1)*(nxread+1) + (i-1) + 1
          elmnts(num_elmnt*ngnod+2) = j*(nxread+1) + (i-1) + 1
          elmnts(num_elmnt*ngnod+3) = j*(nxread+1) + (i-1)
          elmnts(num_elmnt*ngnod+4) = (nxread+1)*(nzread+1) + (j-1)*nxread + (i-1)
          elmnts(num_elmnt*ngnod+5) = (nxread+1)*(nzread+1) + nxread*(nzread+1) + (j-1)*(nxread*2+1) + (i-1)*2 + 2
          elmnts(num_elmnt*ngnod+6) = (nxread+1)*(nzread+1) + j*nxread + (i-1)
          elmnts(num_elmnt*ngnod+7) = (nxread+1)*(nzread+1) + nxread*(nzread+1) + (j-1)*(nxread*2+1) + (i-1)*2
          elmnts(num_elmnt*ngnod+8) = (nxread+1)*(nzread+1) + nxread*(nzread+1) + (j-1)*(nxread*2+1) + (i-1)*2 + 1
          num_elmnt = num_elmnt + 1
       enddo
    enddo

  endif

  end subroutine read_interfaces_file

end module interfaces_file
