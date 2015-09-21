
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

!----
!---- locate_receivers finds the correct position of the receivers
!----

  subroutine locate_receivers(ibool,coord,nspec,nglob,xigll,zigll, &
                          nrec,nrecloc,recloc,which_proc_receiver,nproc,myrank, &
                          st_xval,st_zval,ispec_selected_rec, &
                          xi_receiver,gamma_receiver,station_name,network_name, &
                          x_source,z_source, &
                          coorg,knods,ngnod,npgeo, &
                          x_final_receiver, z_final_receiver)

use specfem_par, only : AXISYM,is_on_the_axis,xiglj,gather_ispec_selected_rec,acoustic,USE_TRICK_FOR_BETTER_PRESSURE

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  include "constants.h"

  integer nrec,nspec,nglob,ngnod,npgeo
  integer, intent(in)  :: nproc, myrank

  integer knods(ngnod,nspec)
  double precision coorg(NDIM,npgeo)

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

! array containing coordinates of the points
  double precision coord(NDIM,nglob)

  integer irec,i,j,ispec,iglob,iter_loop,ix_initial_guess,iz_initial_guess

  double precision x_source,z_source,dist_squared,stele,stbur
  double precision, dimension(nrec)  :: distance_receiver
  double precision xi,gamma,dx,dz,dxi,dgamma

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision zigll(NGLLZ)

  double precision x,z,xix,xiz,gammax,gammaz,jacobian

! use dynamic allocation
  double precision distmin_squared
  double precision, dimension(:), allocatable :: final_distance

! receiver information
  integer  :: nrecloc
  integer, dimension(nrec) :: ispec_selected_rec, recloc
  double precision, dimension(nrec) :: xi_receiver,gamma_receiver

! station information for writing the seismograms
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  double precision, dimension(nrec) :: st_xval,st_zval

! tangential detection
  double precision, dimension(nrec)  :: x_final_receiver, z_final_receiver

  double precision, dimension(nrec,nproc)  :: gather_final_distance
  double precision, dimension(nrec,nproc)  :: gather_xi_receiver, gather_gamma_receiver

  integer, dimension(nrec), intent(inout)  :: which_proc_receiver
  integer  :: ierror

  allocate(gather_ispec_selected_rec(nrec,nproc))
  ierror = 0

! **************

  if (myrank == 0) then
    write(IOUT,*)
    write(IOUT,*) '********************'
    write(IOUT,*) ' locating receivers'
    write(IOUT,*) '********************'
    write(IOUT,*)
    write(IOUT,*) 'reading receiver information from the DATA/STATIONS file'
    write(IOUT,*)
  endif

  open(unit=1,file='DATA/STATIONS',status='old',action='read')

! allocate memory for arrays using number of stations
  allocate(final_distance(nrec))

! loop on all the stations
  do irec=1,nrec

    ! set distance to huge initial value
    distmin_squared = HUGEVAL

    read(1,*) station_name(irec),network_name(irec),st_xval(irec),st_zval(irec),stele,stbur

    ! check that station is not buried, burial is not implemented in current code
    if(abs(stbur) > TINYVAL) call exit_MPI('stations with non-zero burial not implemented yet')

    ! compute distance between source and receiver
    distance_receiver(irec) = sqrt((st_zval(irec)-z_source)**2 + (st_xval(irec)-x_source)**2)

    do ispec=1,nspec

      ! loop only on points inside the element
      ! exclude edges to ensure this point is not shared with other elements
      do j=2,NGLLZ-1
        do i=2,NGLLX-1

          iglob = ibool(i,j,ispec)

          !  we compare squared distances instead of distances themselves to significantly speed up calculations
          dist_squared = (st_xval(irec)-dble(coord(1,iglob)))**2 + (st_zval(irec)-dble(coord(2,iglob)))**2

          ! keep this point if it is closer to the receiver
          if(dist_squared < distmin_squared) then
            distmin_squared = dist_squared
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

    if (AXISYM) then ! TODO
      if (is_on_the_axis(ispec_selected_rec(irec))) then
        xi = xiglj(ix_initial_guess)
      else
        xi = xigll(ix_initial_guess)
      endif
    else
      xi = xigll(ix_initial_guess)
    endif
    gamma = zigll(iz_initial_guess)

    ! iterate to solve the non linear system
    do iter_loop = 1,NUM_ITER
      ! compute coordinates of the new point and derivatives dxi/dx, dxi/dz
      call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                  coorg,knods,ispec_selected_rec(irec),ngnod,nspec,npgeo, &
                  .true.)

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
                .true.)

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

! select one mesh slice for each receiver
#ifdef USE_MPI
  call MPI_GATHER(final_distance(1),nrec,MPI_DOUBLE_PRECISION,&
        gather_final_distance(1,1),nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
  call MPI_GATHER(xi_receiver(1),nrec,MPI_DOUBLE_PRECISION,&
        gather_xi_receiver(1,1),nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
  call MPI_GATHER(gamma_receiver(1),nrec,MPI_DOUBLE_PRECISION,&
        gather_gamma_receiver(1,1),nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
  call MPI_GATHER(ispec_selected_rec(1),nrec,MPI_INTEGER,&
        gather_ispec_selected_rec(1,1),nrec,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

  if ( myrank == 0 ) then
    do irec = 1, nrec
      which_proc_receiver(irec:irec) = minloc(gather_final_distance(irec,:)) - 1
    enddo
  endif

  call MPI_BCAST(which_proc_receiver(1),nrec,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

#else

  gather_final_distance(:,1) = final_distance(:)

  gather_xi_receiver(:,1) = xi_receiver(:)
  gather_gamma_receiver(:,1) = gamma_receiver(:)
  gather_ispec_selected_rec(:,1) = ispec_selected_rec(:)
  which_proc_receiver(:) = 0

#endif

  if (USE_TRICK_FOR_BETTER_PRESSURE) then
    do irec=1,nrec
      if (which_proc_receiver(irec) == myrank) then
        if (.not. acoustic(ispec_selected_rec(irec))) then
          call exit_MPI('USE_TRICK_FOR_BETTER_PRESSURE : receivers must be in acoustic elements')
        endif
      endif
    enddo
  endif

  nrecloc = 0
  do irec = 1, nrec
    if ( which_proc_receiver(irec) == myrank ) then
      nrecloc = nrecloc + 1
      recloc(nrecloc) = irec
    endif
  enddo

  if (myrank == 0) then

    do irec = 1, nrec
      write(IOUT,*)
      write(IOUT,*) 'Station # ',irec,'    ',network_name(irec),station_name(irec)

      if(gather_final_distance(irec,which_proc_receiver(irec)+1) == HUGEVAL) &
        call exit_MPI('error locating receiver')

      write(IOUT,*) '            original x: ',sngl(st_xval(irec))
      write(IOUT,*) '            original z: ',sngl(st_zval(irec))
      write(IOUT,*) '  distance from source: ',sngl(distance_receiver(irec))
      write(IOUT,*) 'closest estimate found: ',sngl(gather_final_distance(irec,which_proc_receiver(irec)+1)), &
                    ' m away'
      write(IOUT,*) ' in element ',gather_ispec_selected_rec(irec,which_proc_receiver(irec)+1)
      write(IOUT,*) ' at process ', which_proc_receiver(irec)
      write(IOUT,*) ' at xi,gamma coordinates = ',gather_xi_receiver(irec,which_proc_receiver(irec)+1),&
                                  gather_gamma_receiver(irec,which_proc_receiver(irec)+1)
      write(IOUT,*)
    enddo

    write(IOUT,*)
    write(IOUT,*) 'end of receiver detection'
    write(IOUT,*)

  endif

  ! deallocate arrays
  deallocate(final_distance)

  end subroutine locate_receivers

