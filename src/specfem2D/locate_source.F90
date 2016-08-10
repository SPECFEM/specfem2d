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
!---- locate_source finds the correct position of the (point force/momen-tensor) source
!----

  subroutine locate_source(ibool,coord,nspec,nglob,xigll,zigll, &
                           x_source,z_source, &
                           ispec_selected_source,islice_selected_source, &
                           NPROC,myrank, &
                           xi_source,gamma_source,coorg,knods,ngnod,npgeo,iglob_source,is_force_source)

  use constants, only: NDIM,NGLLX,NGLLZ,IMAIN,HUGEVAL,TINYVAL,NUM_ITER,USE_BEST_LOCATION_FOR_SOURCE,SOURCE_IS_MOVING, &
    IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC,IDOMAIN_GRAVITOACOUSTIC

  use specfem_par, only: AXISYM,is_on_the_axis,xiglj, &
    ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic,ispec_is_gravitoacoustic

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer,intent(in) :: nspec,nglob

  integer, dimension(NGLLX,NGLLZ,nspec),intent(in) :: ibool

  ! array containing coordinates of the points
  double precision,intent(in) :: coord(NDIM,nglob)

  ! Gauss-Lobatto-Legendre points of integration
  double precision,intent(in) :: xigll(NGLLX)
  double precision,intent(in) :: zigll(NGLLZ)

  double precision,intent(in) :: x_source,z_source

  integer,intent(in) :: ngnod,npgeo

  integer,intent(in) :: knods(ngnod,nspec)
  double precision,intent(in) :: coorg(NDIM,npgeo)

  ! source information
  integer,intent(out) :: ispec_selected_source,islice_selected_source,iglob_source

  integer,intent(in)  :: NPROC, myrank
  double precision,intent(out) :: xi_source,gamma_source

  logical,intent(in) :: is_force_source

  ! local parameters
  integer :: i,j,ispec,iglob,iter_loop,ix_initial_guess,iz_initial_guess,number_of_iterations
  integer :: imin,imax,jmin,jmax
  integer :: idomain
  double precision :: dist_squared
  double precision :: xi,gamma,dx,dz,dxi,dgamma
  double precision :: x,z,xix,xiz,gammax,gammaz,jacobian
  double precision :: distmin_squared,dist_glob_squared
  double precision :: final_distance

  integer :: is_proc_source
  integer, dimension(1:NPROC)  :: allgather_is_proc_source
  integer, dimension(1)  :: locate_is_proc_source

  ! user output
  if (myrank == 0) then
    if (is_force_source .and. (.not. SOURCE_IS_MOVING)) then ! TODO
      write(IMAIN,*)
      write(IMAIN,*) '*******************************'
      write(IMAIN,*) ' locating force source'
      write(IMAIN,*) '*******************************'
      write(IMAIN,*)
      call flush_IMAIN()
    else
      write(IMAIN,*)
      write(IMAIN,*) '*******************************'
      write(IMAIN,*) ' locating moment-tensor source'
      write(IMAIN,*) '*******************************'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! initializes slice and element which hold source
  is_proc_source = 0
  ispec_selected_source = 0
  islice_selected_source = 0
  iglob_source = 0 ! closest global point

  ! set distance to huge initial value
  distmin_squared = HUGEVAL

  ! determines search range
  if (is_force_source) then
    !  if (USE_BEST_LOCATION_FOR_SOURCE) then
    if (.false.) then
      ! only search best element with inner points (exact location will be interpolated afterwards)
      imin = 2
      imax = NGLLX - 1
      jmin = 2
      jmax = NGLLZ - 1
    else
      ! search best GLL location
      imin = 1
      imax = NGLLX
      jmin = 1
      jmax = NGLLZ
    endif
  else
    ! moment-tensor
    imin = 2
    imax = NGLLX - 1
    jmin = 2
    jmax = NGLLZ - 1
  endif

  ! loops over all elements
  do ispec = 1,nspec
    ! loop only on points inside the element
    ! exclude edges to ensure this point is not shared with other elements
    do j = jmin,jmax
      do i = imin,imax
        iglob = ibool(i,j,ispec)

        !  we compare squared distances instead of distances themselves to significantly speed up calculations
        dist_squared = (x_source-dble(coord(1,iglob)))**2 + (z_source-dble(coord(2,iglob)))**2

        ! keep this point if it is closer to the source
        if (dist_squared < distmin_squared) then
          iglob_source = iglob
          distmin_squared = dist_squared
          ispec_selected_source = ispec
          ix_initial_guess = i
          iz_initial_guess = j
          ! determines domain for outputting element type
          if (ispec_is_acoustic(ispec)) then
            idomain = IDOMAIN_ACOUSTIC
          else if (ispec_is_elastic(ispec)) then
            idomain = IDOMAIN_ELASTIC
          else if (ispec_is_poroelastic(ispec)) then
            idomain = IDOMAIN_POROELASTIC
          else if (ispec_is_gravitoacoustic(ispec)) then
            idomain = IDOMAIN_GRAVITOACOUSTIC
          else
            stop 'Invalid element type in locating source found!'
          endif
        endif
      enddo
    enddo
  enddo

  ! global minimum distance computed over all processes
  call min_all_all_dp(distmin_squared, dist_glob_squared)

  ! check if this process contains the source
  if (abs(sqrt(dist_glob_squared) - sqrt(distmin_squared)) < TINYVAL ) is_proc_source = 1

  ! master collects info
  call gather_all_singlei(is_proc_source,allgather_is_proc_source,NPROC)
  if (myrank == 0) then
    ! select slice with maximum rank which contains source
    locate_is_proc_source = maxloc(allgather_is_proc_source) - 1
    islice_selected_source = locate_is_proc_source(1)
  endif
  ! selects slice which holds source
  call bcast_all_singlei(islice_selected_source)

! ****************************************
! find the best (xi,gamma) for each source
! ****************************************

! use initial guess in xi and gamma
  if (AXISYM) then
    if (is_on_the_axis(ispec_selected_source)) then
      xi = xiglj(ix_initial_guess)
    else
      xi = xigll(ix_initial_guess)
    endif
  else
    xi = xigll(ix_initial_guess)
  endif
  gamma = zigll(iz_initial_guess)

! iterate to solve the non linear system
  if (USE_BEST_LOCATION_FOR_SOURCE) then
    number_of_iterations = NUM_ITER
  else
    number_of_iterations = 0 ! this means that the loop below will not be executed, i.e. we will not iterate
  endif

  do iter_loop = 1,number_of_iterations

    ! recompute jacobian for the new point
    call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                  coorg,knods,ispec_selected_source,ngnod,nspec,npgeo, &
                  .true.)

    ! compute distance to target location
    dx = - (x - x_source)
    dz = - (z - z_source)

    ! compute increments
    dxi  = xix*dx + xiz*dz
    dgamma = gammax*dx + gammaz*dz

    ! update values
    xi = xi + dxi
    gamma = gamma + dgamma

    ! impose that we stay in that element
    ! (useful if user gives a source outside the mesh for instance)
    ! we can go slightly outside the [1,1] segment since with finite elements
    ! the polynomial solution is defined everywhere
    ! this can be useful for convergence of itertive scheme with distorted elements
    if (is_force_source) then
      if (xi > 1.010d0) xi = 1.010d0
      if (xi < -1.010d0) xi = -1.010d0
      if (gamma > 1.010d0) gamma = 1.010d0
      if (gamma < -1.010d0) gamma = -1.010d0
    else
      if (xi > 1.10d0) xi = 1.10d0
      if (xi < -1.10d0) xi = -1.10d0
      if (gamma > 1.10d0) gamma = 1.10d0
      if (gamma < -1.10d0) gamma = -1.10d0
    endif
! end of non linear iterations
  enddo

! compute final coordinates of point found
  call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian, &
                          coorg,knods,ispec_selected_source,ngnod,nspec,npgeo, &
                          .true.)

! store xi,gamma of point found
  xi_source = xi
  gamma_source = gamma

! compute final distance between asked and found
  final_distance = sqrt((x_source-x)**2 + (z_source-z)**2)

#ifdef USE_MPI
  ! for MPI version, gather information from all the nodes
  if (islice_selected_source /= 0) then
    ! source is in another slice than the master process
    if (myrank == islice_selected_source) then
      ! send information from process holding source
      call send_singlei(ispec_selected_source,0,0)
      call send_singlei(idomain,0,1)
      call send_singledp(xi_source,0,2)
      call send_singledp(gamma_source,0,3)
      call send_singledp(final_distance,0,4)
    else if (myrank == 0) then
      ! master collects
      call recv_singlei(ispec_selected_source,islice_selected_source,0)
      call recv_singlei(idomain,islice_selected_source,1)
      call recv_singledp(xi_source,islice_selected_source,2)
      call recv_singledp(gamma_source,islice_selected_source,3)
      call recv_singledp(final_distance,islice_selected_source,4)
    endif
  endif
#endif

  ! user output
  if (myrank == 0) then
    if ((.not. is_force_source) .or. (is_force_source .and. (.not. SOURCE_IS_MOVING))) then
      write(IMAIN,*)
      if (is_force_source) then
        write(IMAIN,*) 'Force source:'
      else
        write(IMAIN,*) 'Moment-tensor source:'
      endif
      write(IMAIN,*) '            original x: ',sngl(x_source)
      write(IMAIN,*) '            original z: ',sngl(z_source)
      write(IMAIN,*) 'closest estimate found: ',sngl(final_distance),' m away'
      write(IMAIN,*) ' in rank ',islice_selected_source
      write(IMAIN,*) ' in element ',ispec_selected_source
      if (idomain == IDOMAIN_ACOUSTIC) then
        write(IMAIN,*) ' in acoustic domain'
      else if (idomain == IDOMAIN_ELASTIC) then
        write(IMAIN,*) ' in elastic domain'
      else if (idomain == IDOMAIN_POROELASTIC) then
        write(IMAIN,*) ' in poroelastic domain'
      else if (idomain == IDOMAIN_GRAVITOACOUSTIC) then
        write(IMAIN,*) ' in gravitoacoustic domain'
      else
        write(IMAIN,*) ' in unknown domain'
      endif
      write(IMAIN,*) ' at xi,gamma coordinates = ',xi_source,gamma_source
      write(IMAIN,*)

      write(IMAIN,*)
      if (is_force_source) then
        write(IMAIN,*) 'end of force source detection'
      else
        write(IMAIN,*) 'end of moment-tensor source detection'
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! check position
    if (final_distance == HUGEVAL) call exit_MPI(myrank,'Error locating source')
  endif

  end subroutine locate_source

