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

!----
!---- locate_receivers finds the correct position of the receivers
!----

  subroutine locate_receivers(ibool,coord,nspec,nglob,xigll,zigll, &
                              nrec,nrecloc,recloc,islice_selected_rec,NPROC,myrank, &
                              st_xval,st_zval,ispec_selected_rec, &
                              xi_receiver,gamma_receiver,station_name,network_name, &
                              x_source,z_source, &
                              coorg,knods,ngnod,npgeo, &
                              x_final_receiver, z_final_receiver)

  use constants, only: NDIM,NGLLX,NGLLZ,MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME, &
    IMAIN,HUGEVAL,TINYVAL,NUM_ITER,mygroup,MAX_STRING_LEN,IN_DATA_FILES

  use specfem_par, only: AXISYM,is_on_the_axis,xiglj,ispec_is_acoustic,USE_TRICK_FOR_BETTER_PRESSURE,NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

  integer :: nrec,nspec,nglob,ngnod,npgeo
  integer, intent(in)  :: NPROC, myrank

  integer :: knods(ngnod,nspec)
  double precision :: coorg(NDIM,npgeo)

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

! array containing coordinates of the points
  double precision :: coord(NDIM,nglob)

  integer :: irec,i,j,ispec,iglob,iter_loop,ix_initial_guess,iz_initial_guess

  double precision :: x_source,z_source,dist_squared,stele,stbur
  double precision, dimension(nrec)  :: distance_receiver
  double precision :: xi,gamma,dx,dz,dxi,dgamma

! Gauss-Lobatto-Legendre points of integration
  double precision :: xigll(NGLLX)
  double precision :: zigll(NGLLZ)

  double precision :: x,z,xix,xiz,gammax,gammaz,jacobian

! use dynamic allocation
  double precision :: distmin_squared
  double precision, dimension(:), allocatable :: final_distance
  double precision :: final_distance_this_element

! receiver information
  integer  :: nrecloc
  integer, dimension(nrec) :: ispec_selected_rec, recloc
  double precision, dimension(nrec) :: xi_receiver,gamma_receiver

!! DK DK dec 2017: also loop on all the elements in contact with the initial guess element to improve accuracy of estimate
  logical, dimension(nglob) :: flag_topological
  integer :: number_of_mesh_elements_for_the_initial_guess
  integer, dimension(:), allocatable :: array_of_all_elements_of_ispec_selected_rec

! station information for writing the seismograms
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  double precision, dimension(nrec) :: st_xval,st_zval

! tangential detection
  double precision, dimension(nrec)  :: x_final_receiver, z_final_receiver

  double precision, dimension(nrec,NPROC)  :: gather_final_distance
  double precision, dimension(nrec,NPROC)  :: gather_xi_receiver, gather_gamma_receiver

  integer, dimension(nrec), intent(inout)  :: islice_selected_rec
  integer, dimension(:,:), allocatable  :: gather_ispec_selected_rec
  integer  :: ier
  character(len=MAX_STRING_LEN) :: stations_filename,path_to_add

  stations_filename = trim(IN_DATA_FILES)//'STATIONS'

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    stations_filename = path_to_add(1:len_trim(path_to_add))//stations_filename(1:len_trim(stations_filename))
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating receivers'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    write(IMAIN,*) 'reading receiver information from the '//trim(stations_filename)//' file'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  open(unit= 1,file=trim(stations_filename),status='old',action='read')

! allocate memory for arrays using number of stations
  allocate(final_distance(nrec))

! loop on all the stations
  do irec = 1,nrec

    ! set distance to huge initial value
    distmin_squared = HUGEVAL

    read(1,*) station_name(irec),network_name(irec),st_xval(irec),st_zval(irec),stele,stbur

    ! check that station is not buried, burial is not implemented in current code
    if (abs(stbur) > TINYVAL) call exit_MPI(myrank,'stations with non-zero burial not implemented yet')

    ! compute distance between source and receiver
    distance_receiver(irec) = sqrt((st_zval(irec)-z_source)**2 + (st_xval(irec)-x_source)**2)

    do ispec = 1,nspec

      ! loop only on points inside the element
      ! exclude edges to ensure this point is not shared with other elements
      do j = 1,NGLLZ
        do i = 1,NGLLX

          iglob = ibool(i,j,ispec)

          !  we compare squared distances instead of distances themselves to significantly speed up calculations
          dist_squared = (st_xval(irec)-coord(1,iglob))**2 + (st_zval(irec)-coord(2,iglob))**2

          ! keep this point if it is closer to the receiver
          if (dist_squared < distmin_squared) then

            ! this statement is useless, it is there because of a bug in some releases of the Intel ifort compiler
            ! in which at optimization level -O3 the "if" statement above undergoes heavy optimization and
            ! something goes wrong in the compiler, resulting in a comparison in the "if" statement that can
            ! be performed wrong; adding this dummy call in the "if" means that the compiler cannot optimize
            ! as aggressively, and the compiler bug is not triggered;
            ! thus, please do *NOT* remove this statement (it took us a while to discover this nasty compiler problem).
            ! GNU gfortran and all other compilers we have tested are fine and do not have any problem even if this
            ! statement is removed. Some releases of Intel ifort are also OK.
            call dummy_routine()

            distmin_squared = dist_squared
            ispec_selected_rec(irec) = ispec
            ix_initial_guess = i
            iz_initial_guess = j
          endif

        enddo
      enddo

    ! end of loop on all the spectral elements
    enddo

!! DK DK dec 2017: also loop on all the elements in contact with the initial guess element to improve accuracy of estimate
  flag_topological(:) = .false.

! mark the four corners of the initial guess element
  flag_topological(ibool(1,1,ispec_selected_rec(irec))) = .true.
  flag_topological(ibool(NGLLX,1,ispec_selected_rec(irec))) = .true.
  flag_topological(ibool(NGLLX,NGLLZ,ispec_selected_rec(irec))) = .true.
  flag_topological(ibool(1,NGLLZ,ispec_selected_rec(irec))) = .true.

! loop on all the elements to count how many are shared with the initial guess
  number_of_mesh_elements_for_the_initial_guess = 1
  do ispec = 1,nspec
    if (ispec == ispec_selected_rec(irec)) cycle
    ! loop on the four corners only, no need to loop on the rest since we just want to detect adjacency
    do j = 1,NGLLZ,NGLLZ-1
      do i = 1,NGLLX,NGLLX-1
        if (flag_topological(ibool(i,j,ispec))) then
          ! this element is in contact with the initial guess
          number_of_mesh_elements_for_the_initial_guess = number_of_mesh_elements_for_the_initial_guess + 1
          ! let us not count it more than once, it may have a full edge in contact with it and would then be counted twice
          goto 700
        endif
      enddo
    enddo
    700 continue
  enddo

! now that we know the number of elements, we can allocate the list of elements and create it
  allocate(array_of_all_elements_of_ispec_selected_rec(number_of_mesh_elements_for_the_initial_guess))

! first store the initial guess itself
  number_of_mesh_elements_for_the_initial_guess = 1
  array_of_all_elements_of_ispec_selected_rec(number_of_mesh_elements_for_the_initial_guess) = ispec_selected_rec(irec)

! then store all the others
  do ispec = 1,nspec
    if (ispec == ispec_selected_rec(irec)) cycle
    ! loop on the four corners only, no need to loop on the rest since we just want to detect adjacency
    do j = 1,NGLLZ,NGLLZ-1
      do i = 1,NGLLX,NGLLX-1
        if (flag_topological(ibool(i,j,ispec))) then
          ! this element is in contact with the initial guess
          number_of_mesh_elements_for_the_initial_guess = number_of_mesh_elements_for_the_initial_guess + 1
          array_of_all_elements_of_ispec_selected_rec(number_of_mesh_elements_for_the_initial_guess) = ispec
          ! let us not count it more than once, it may have a full edge in contact with it and would then be counted twice
          goto 800
        endif
      enddo
    enddo
    800 continue
  enddo

!! DK DK dec 2017
  final_distance(irec) = HUGEVAL

  do i = 1,number_of_mesh_elements_for_the_initial_guess

!! DK DK dec 2017 set initial guess in the middle of the element, since we computed the true one only for the true initial guess
!! DK DK dec 2017 the nonlinear process below will converge anyway
    if (i > 1) then
      ix_initial_guess = NGLLX / 2
      iz_initial_guess = NGLLZ / 2
    endif

    ispec = array_of_all_elements_of_ispec_selected_rec(i)

    ! ****************************************
    ! find the best (xi,gamma) for each receiver
    ! ****************************************

    ! use initial guess in xi and gamma
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        xi = xiglj(ix_initial_guess)
      else
        xi = xigll(ix_initial_guess)
      endif
    else
      xi = xigll(ix_initial_guess)
    endif
    gamma = zigll(iz_initial_guess)

    ! iterate to solve the nonlinear system
    do iter_loop = 1,NUM_ITER

      ! compute coordinates of the new point and derivatives dxi/dx, dxi/dz
      call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                  coorg,knods,ispec,ngnod,nspec,npgeo,.true.)

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
      if (xi > 1.01d0) xi = 1.01d0
      if (xi < -1.01d0) xi = -1.01d0
      if (gamma > 1.01d0) gamma = 1.01d0
      if (gamma < -1.01d0) gamma = -1.01d0

    ! end of nonlinear iterations
    enddo

    ! compute final coordinates of point found
    call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                coorg,knods,ispec,ngnod,nspec,npgeo,.true.)

! compute final distance between asked and found
    final_distance_this_element = sqrt((st_xval(irec)-x)**2 + (st_zval(irec)-z)**2)

! if we have found an element that gives a shorter distance
  if (final_distance_this_element < final_distance(irec)) then
!   store element number found
    ispec_selected_rec(irec) = ispec

    ! store xi,gamma of point found
    xi_receiver(irec) = xi
    gamma_receiver(irec) = gamma

    x_final_receiver(irec) = x
    z_final_receiver(irec) = z

!   store final distance between asked and found
    final_distance(irec) = final_distance_this_element
  endif

!! DK DK dec 2017
  enddo

!! DK DK dec 2017: also loop on all the elements in contact with the initial guess element to improve accuracy of estimate
    deallocate(array_of_all_elements_of_ispec_selected_rec)

  enddo ! of loop on all the receivers

  ! close receiver file
  close(1)

  ! select one mesh slice for each receiver
  allocate(gather_ispec_selected_rec(nrec,NPROC),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating gather array')

  ! gathers infos onto master process
  call gather_all_dp(final_distance(1),nrec,gather_final_distance(1,1),nrec,NPROC)
  call gather_all_dp(xi_receiver(1),nrec,gather_xi_receiver(1,1),nrec,NPROC)
  call gather_all_dp(gamma_receiver(1),nrec,gather_gamma_receiver(1,1),nrec,NPROC)

  call gather_all_i(ispec_selected_rec(1),nrec,gather_ispec_selected_rec(1,1),nrec, NPROC)

  if (myrank == 0) then
    ! selects best slice which minimum distance to receiver location
    do irec = 1, nrec
      islice_selected_rec(irec:irec) = minloc(gather_final_distance(irec,:)) - 1
    enddo
  endif
  call bcast_all_i(islice_selected_rec(1),nrec)

  if (USE_TRICK_FOR_BETTER_PRESSURE) then
    do irec= 1,nrec
      if (myrank == islice_selected_rec(irec)) then
        if (.not. ispec_is_acoustic(ispec_selected_rec(irec))) then
          call exit_MPI(myrank,'USE_TRICK_FOR_BETTER_PRESSURE : receivers must be in acoustic elements')
        endif
      endif
    enddo
  endif

  ! counts local receivers in this slice
  nrecloc = 0
  do irec = 1, nrec
    if (myrank == islice_selected_rec(irec)) then
      nrecloc = nrecloc + 1
      recloc(nrecloc) = irec
    endif
  enddo

  ! user output
  if (myrank == 0) then
    do irec = 1, nrec
      write(IMAIN,*)
      write(IMAIN,*) 'Station # ',irec,'    ',network_name(irec),station_name(irec)

      if (gather_final_distance(irec,islice_selected_rec(irec)+1) == HUGEVAL) &
        call exit_MPI(myrank,'Error locating receiver')

      write(IMAIN,*) '            original x: ',sngl(st_xval(irec))
      write(IMAIN,*) '            original z: ',sngl(st_zval(irec))
      write(IMAIN,*) 'Closest estimate found: ',sngl(gather_final_distance(irec,islice_selected_rec(irec)+1)), &
                    ' m away'
      write(IMAIN,*) ' in element ',gather_ispec_selected_rec(irec,islice_selected_rec(irec)+1)
      write(IMAIN,*) ' in rank ', islice_selected_rec(irec)
      write(IMAIN,*) ' at xi,gamma coordinates = ',gather_xi_receiver(irec,islice_selected_rec(irec)+1), &
                                  gather_gamma_receiver(irec,islice_selected_rec(irec)+1)
      write(IMAIN,*) 'Distance from source: ',sngl(distance_receiver(irec)),' m'
      write(IMAIN,*)
    enddo

    write(IMAIN,*)
    write(IMAIN,*) 'end of receiver detection'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! deallocate arrays
  deallocate(final_distance)
  deallocate(gather_ispec_selected_rec)

  end subroutine locate_receivers

