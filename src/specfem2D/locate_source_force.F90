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
!---- locate_source_force finds the correct position of the point force source
!----

  subroutine locate_source_force(ibool,coord,nspec,nglob,xigll,zigll, &
                                 x_source,z_source, &
                                 ispec_selected_source,is_proc_source,nb_proc_source,NPROC,myrank, &
                                 xi_source,gamma_source,coorg,knods,ngnod,npgeo,iglob_source)

  use constants,only: NDIM,NGLLX,NGLLZ,IMAIN,HUGEVAL,TINYVAL,NUM_ITER,USE_BEST_LOCATION_FOR_SOURCE

  use specfem_par, only : AXISYM,is_on_the_axis,xiglj

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  integer :: nspec,nglob,ngnod,npgeo

  integer :: knods(ngnod,nspec)
  double precision :: coorg(NDIM,npgeo)

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

! array containing coordinates of the points
  double precision :: coord(NDIM,nglob)

  integer :: i,j,ispec,iglob,iter_loop,ix_initial_guess,iz_initial_guess,number_of_iterations

  double precision :: x_source,z_source,dist_squared
  double precision :: xi,gamma,dx,dz,dxi,dgamma

! Gauss-Lobatto-Legendre points of integration
  double precision :: xigll(NGLLX)
  double precision :: zigll(NGLLZ)

  double precision :: x,z,xix,xiz,gammax,gammaz,jacobian
  double precision :: distmin_squared,final_distance,dist_glob_squared

! source information
  integer :: ispec_selected_source,is_proc_source,nb_proc_source,iglob_source
  integer, intent(in)  :: NPROC, myrank
  double precision :: xi_source,gamma_source

#ifdef USE_MPI
  integer, dimension(1:NPROC)  :: allgather_is_proc_source
  integer, dimension(1)  :: locate_is_proc_source
  integer  :: ierror
#endif



! **************
  if (myrank == 0 .or. NPROC == 1) then
    write(IMAIN,*)
    write(IMAIN,*) '*******************************'
    write(IMAIN,*) ' locating force source'
    write(IMAIN,*) '*******************************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

! set distance to huge initial value
  distmin_squared = HUGEVAL

  is_proc_source = 0

  do ispec = 1,nspec

! loop only on points inside the element
! exclude edges to ensure this point is not shared with other elements
     do j = 2,NGLLZ-1
        do i = 2,NGLLX-1

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
           endif

        enddo
     enddo

! end of loop on all the spectral elements
  enddo

  ! global minimum distance computed over all processes
  call min_all_all_dp(distmin_squared, dist_glob_squared)

  ! check if this process contains the source
  if (abs(sqrt(dist_glob_squared) - sqrt(distmin_squared)) < TINYVAL ) is_proc_source = 1

  ! determining the number of processes that contain the source
  ! (useful when the source is located on an interface)
  call sum_all_all_i(is_proc_source, nb_proc_source)


#ifdef USE_MPI
  ! when several processes contain the source, we elect one of them (minimum rank).
  if (nb_proc_source > 1) then

     call MPI_ALLGATHER(is_proc_source, 1, MPI_INTEGER, allgather_is_proc_source(1), &
                        1, MPI_INTEGER, MPI_COMM_WORLD, ierror)
     locate_is_proc_source = maxloc(allgather_is_proc_source) - 1

     if (myrank /= locate_is_proc_source(1)) then
        is_proc_source = 0
     endif
     nb_proc_source = 1

  endif

#endif

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

  do iter_loop = 1,NUM_ITER

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
    if (xi > 1.10d0) xi = 1.10d0
    if (xi < -1.10d0) xi = -1.10d0
    if (gamma > 1.10d0) gamma = 1.10d0
    if (gamma < -1.10d0) gamma = -1.10d0

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

  if (is_proc_source == 1) then
     write(IMAIN,*)
     write(IMAIN,*) 'Force source:'

     if (final_distance == HUGEVAL) call exit_MPI(myrank,'Error locating force source')

     write(IMAIN,*) '            original x: ',sngl(x_source)
     write(IMAIN,*) '            original z: ',sngl(z_source)
     write(IMAIN,*) 'closest estimate found: ',sngl(final_distance),' m away'
#ifdef USE_MPI
     write(IMAIN,*) ' in rank ',myrank
#endif
     write(IMAIN,*) ' in element ',ispec_selected_source
     write(IMAIN,*) ' at xi,gamma coordinates = ',xi_source,gamma_source
     write(IMAIN,*)

     write(IMAIN,*)
     write(IMAIN,*) 'end of force source detection'
     write(IMAIN,*)
     call flush_IMAIN()
  endif

  end subroutine locate_source_force

