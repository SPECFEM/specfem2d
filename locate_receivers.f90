
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!                         Dimitri Komatitsch
!                     University of Pau, France
!
!                          (c) April 2007
!
!========================================================================

!----
!---- locate_receivers finds the correct position of the receivers
!----

  subroutine locate_receivers(ibool,coord,nspec,npoin,xigll,zigll,nrec,st_xval,st_zval,ispec_selected_rec, &
                 xi_receiver,gamma_receiver,station_name,network_name,x_source,z_source,coorg,knods,ngnod,npgeo)

  implicit none

  include "constants.h"

  integer nrec,nspec,npoin,ngnod,npgeo

  integer knods(ngnod,nspec)
  double precision coorg(NDIM,npgeo)

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

! array containing coordinates of the points
  double precision coord(NDIM,npoin)

  integer nrec_dummy,irec,i,j,ispec,iglob,iter_loop,ix_initial_guess,iz_initial_guess

  double precision x_source,z_source,dist,stele,stbur,distance_receiver
  double precision xi,gamma,dx,dz,dxi,dgamma

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision zigll(NGLLZ)

  double precision x,z,xix,xiz,gammax,gammaz,jacobian

! use dynamic allocation
  double precision distmin
  double precision, dimension(:), allocatable :: final_distance

! receiver information
  integer, dimension(nrec) :: ispec_selected_rec
  double precision, dimension(nrec) :: xi_receiver,gamma_receiver

! station information for writing the seismograms
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  double precision, dimension(nrec) :: st_xval,st_zval

! **************

  write(IOUT,*)
  write(IOUT,*) '********************'
  write(IOUT,*) ' locating receivers'
  write(IOUT,*) '********************'
  write(IOUT,*)
  write(IOUT,*) 'reading receiver information from the DATA/STATIONS file'
  write(IOUT,*)

! get number of stations from receiver file
  open(unit=1,file='DATA/STATIONS',status='old')
  read(1,*) nrec_dummy

  if(nrec_dummy /= nrec) stop 'problem with number of receivers'

! allocate memory for arrays using number of stations
  allocate(final_distance(nrec))

! loop on all the stations
  do irec=1,nrec

! set distance to huge initial value
  distmin=HUGEVAL

    read(1,*) station_name(irec),network_name(irec),st_xval(irec),st_zval(irec),stele,stbur

! check that station is not buried, burial is not implemented in current code
    if(abs(stbur) > TINYVAL) stop 'stations with non-zero burial not implemented yet'

! compute distance between source and receiver
      distance_receiver = sqrt((st_zval(irec)-z_source)**2 + (st_xval(irec)-x_source)**2)

      do ispec=1,nspec

! loop only on points inside the element
! exclude edges to ensure this point is not shared with other elements
        do j=2,NGLLZ-1
          do i=2,NGLLX-1

            iglob = ibool(i,j,ispec)
            dist = sqrt((st_xval(irec)-dble(coord(1,iglob)))**2 + (st_zval(irec)-dble(coord(2,iglob)))**2)

!           keep this point if it is closer to the receiver
            if(dist < distmin) then
              distmin = dist
              ispec_selected_rec(irec) = ispec
              ix_initial_guess = i
              iz_initial_guess = j
            endif

          enddo
        enddo

! end of loop on all the spectral elements
      enddo

! ****************************************
! find the best (xi,gamma) for each receiver
! ****************************************

! use initial guess in xi and gamma
        xi = xigll(ix_initial_guess)
        gamma = zigll(iz_initial_guess)

! iterate to solve the non linear system
  do iter_loop = 1,NUM_ITER

! recompute jacobian for the new point
    call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian,coorg,knods,ispec_selected_rec(irec),ngnod,nspec,npgeo)

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
    call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian,coorg,knods,ispec_selected_rec(irec),ngnod,nspec,npgeo)

! store xi,gamma of point found
  xi_receiver(irec) = xi
  gamma_receiver(irec) = gamma

! compute final distance between asked and found
  final_distance(irec) = sqrt((st_xval(irec)-x)**2 + (st_zval(irec)-z)**2)

    write(IOUT,*)
    write(IOUT,*) 'Station # ',irec,'    ',station_name(irec),network_name(irec)

    if(final_distance(irec) == HUGEVAL) stop 'error locating receiver'

    write(IOUT,*) '            original x: ',sngl(st_xval(irec))
    write(IOUT,*) '            original z: ',sngl(st_zval(irec))
    write(IOUT,*) '  distance from source: ',sngl(distance_receiver)
    write(IOUT,*) 'closest estimate found: ',sngl(final_distance(irec)),' m away'
    write(IOUT,*) ' in element ',ispec_selected_rec(irec)
    write(IOUT,*) ' at xi,gamma coordinates = ',xi_receiver(irec),gamma_receiver(irec)
    write(IOUT,*)

  enddo

! close receiver file
  close(1)

! display maximum error for all the receivers
  write(IOUT,*) 'maximum error in location of all the receivers: ',sngl(maxval(final_distance(:))),' m'

  write(IOUT,*)
  write(IOUT,*) 'end of receiver detection'
  write(IOUT,*)

! deallocate arrays
  deallocate(final_distance)

  end subroutine locate_receivers

