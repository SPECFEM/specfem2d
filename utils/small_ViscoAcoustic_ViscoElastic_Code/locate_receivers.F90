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

!----
!---- locate_receivers finds the correct position of the receivers
!----

  subroutine locate_receivers(ibool,coord,nspec,nglob,xigll,zigll, &
                              nrec,st_xval,st_zval,ispec_selected_rec, &
                              xi_receiver,gamma_receiver, &
                              x_source,z_source, &
                              coorg,knods,ngnod,npgeo, &
                              x_final_receiver, z_final_receiver,NDIM,NGLLX,NGLLZ)

  implicit none

  integer :: nrec,nspec,nglob,ngnod,npgeo
  integer :: NDIM,NGLLX,NGLLZ
  integer :: knods(ngnod,nspec)
  double precision :: coorg(NDIM,npgeo)

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

! array containing coordinates of the points
  double precision :: coord(NDIM,nglob)

  integer :: irec,i,j,ispec,iglob,iter_loop

  double precision :: x_source,z_source,dist_squared
  double precision :: xi,gamma,dx,dz,dxi,dgamma

! Gauss-Lobatto-Legendre points of integration
  double precision :: xigll(NGLLX)
  double precision :: zigll(NGLLZ)

  double precision :: x,z,xix,xiz,gammax,gammaz,jacobian

! use dynamic allocation
  double precision :: distmin_squared
  double precision, dimension(:), allocatable :: final_distance

! receiver information
  integer, dimension(nrec) :: ispec_selected_rec
  double precision, dimension(nrec) :: xi_receiver,gamma_receiver,distance_receiver

  double precision, dimension(nrec) :: st_xval,st_zval

! tangential detection
  double precision, dimension(nrec)  :: x_final_receiver, z_final_receiver
  integer, dimension(nrec) :: ix_initial_guess,iz_initial_guess

  double precision :: HUGEVAL = 1.0d30
  integer :: NUM_ITER = 10

  ! user output
    write(*,*)
    write(*,*) '********************'
    write(*,*) ' locating receivers'
    write(*,*) '********************'
    write(*,*)
    write(*,*) 'reading receiver information from the DATA/STATIONS file'
    write(*,*)

  open(unit= 1,file='DATA/STATIONS',status='old',action='read')

! allocate memory for arrays using number of stations
  allocate(final_distance(nrec))

! loop on all the stations
  do irec= 1,nrec

    ! set distance to huge initial value
    distmin_squared = HUGEVAL

    read(1,*) st_xval(irec),st_zval(irec)
    ! compute distance between source and receiver
    distance_receiver(irec) = sqrt((st_zval(irec)-z_source)**2 + (st_xval(irec)-x_source)**2)

    do ispec= 1,nspec

      ! loop only on points inside the element
      ! exclude edges to ensure this point is not shared with other elements
      do j = 2,NGLLZ-1
        do i = 2,NGLLX-1

          iglob = ibool(i,j,ispec)

          !  we compare squared distances instead of distances themselves to significantly speed up calculations
          dist_squared = (st_xval(irec)-dble(coord(1,iglob)))**2 + (st_zval(irec)-dble(coord(2,iglob)))**2

          ! keep this point if it is closer to the receiver
          if (dist_squared < distmin_squared) then
            distmin_squared = dist_squared
            ispec_selected_rec(irec) = ispec
            ix_initial_guess(irec) = i
            iz_initial_guess(irec) = j
          endif

        enddo
      enddo

    ! end of loop on all the spectral elements
    enddo


    ! ****************************************
    ! find the best (xi,gamma) for each receiver
    ! ****************************************
    ! use initial guess in xi and gamma

    xi = xigll(ix_initial_guess(irec))
    gamma = zigll(iz_initial_guess(irec))

    ! iterate to solve the non linear system
    do iter_loop = 1,NUM_ITER
      ! compute coordinates of the new point and derivatives dxi/dx, dxi/dz

      call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                  coorg,knods,ispec_selected_rec(irec),ngnod,nspec,npgeo, &
                  NDIM)

      ! compute distance to target location
      dx = - (x - st_xval(irec))
      dz = - (z - st_zval(irec))

      ! compute increments
      dxi  = xix*dx + xiz*dz
      dgamma = gammax*dx + gammaz*dz

      ! update values
      xi = xi + dxi
      gamma = gamma + dgamma

      ! impose that we stay in that element
      ! (useful if user gives a receiver outside the mesh for instance)
      ! we can go slightly outside the [1,1] segment since with finite elements
      ! the polynomial solution is defined everywhere
      ! this can be useful for convergence of itertive scheme with distorted elements
      if (xi > 1.10d0) xi = 1.10d0
      if (xi < -1.10d0) xi = -1.10d0
      if (gamma > 1.10d0) gamma = 1.10d0
      if (gamma < -1.10d0) gamma = -1.10d0

    ! end of non linear iterations
    enddo

    ! compute final coordinates of point found
    call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                coorg,knods,ispec_selected_rec(irec),ngnod,nspec,npgeo, &
                NDIM)

    ! store xi,gamma of point found
    xi_receiver(irec) = xi
    gamma_receiver(irec) = gamma

    ! compute final distance between asked and found
    final_distance(irec) = sqrt((st_xval(irec)-x)**2 + (st_zval(irec)-z)**2)

    x_final_receiver(irec) = x
    z_final_receiver(irec) = z

  enddo

  ! close receiver file
  close(1)

    do irec = 1, nrec
      write(*,*)
      write(*,*) 'Station # ',irec,'    '

      if (final_distance(irec) == HUGEVAL) &
        stop('Error locating receiver')

      write(*,*) '            original x: ',sngl(st_xval(irec))
      write(*,*) '            original z: ',sngl(st_zval(irec))
      write(*,*) '  distance from source: ',sngl(distance_receiver(irec))
      write(*,*) 'closest estimate found: ',sngl(final_distance(irec)), &
                    ' m away'
      write(*,*) ' in element ',ispec_selected_rec(irec)
      write(*,*) ' closest point :',ix_initial_guess(irec),iz_initial_guess(irec)
      write(*,*) ' at xi,gamma coordinates = ',xi_receiver(irec), &
                                  gamma_receiver(irec)
      write(*,*)
    enddo

    write(*,*)
    write(*,*) 'end of receiver detection'
    write(*,*)

  ! deallocate arrays
  deallocate(final_distance)

  end subroutine locate_receivers

