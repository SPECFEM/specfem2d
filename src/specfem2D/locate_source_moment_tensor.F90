
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
!---- locate_source_moment_tensor finds the correct position of the moment-tensor source
!----

  subroutine locate_source_moment_tensor(ibool,coord,nspec,nglob, &
               xigll,zigll,x_source,z_source, &
               ispec_selected_source,is_proc_source,nb_proc_source,nproc,myrank, &
               xi_source,gamma_source,coorg,knods,ngnod,npgeo)

  use specfem_par, only : AXISYM,is_on_the_axis,xiglj

#ifdef USE_MPI
  use mpi
#endif

  implicit none

  include "constants.h"

  integer nspec,nglob,ngnod,npgeo

  integer knods(ngnod,nspec)
  double precision coorg(NDIM,npgeo)

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

! array containing coordinates of the points
  double precision coord(NDIM,nglob)

  integer i,j,ispec,iglob,iter_loop,ix_initial_guess,iz_initial_guess

  double precision x_source,z_source,dist_squared
  double precision xi,gamma,dx,dz,dxi,dgamma

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision zigll(NGLLZ)

  double precision x,z,xix,xiz,gammax,gammaz,jacobian
  double precision distmin_squared,final_distance,dist_glob_squared

! source information
  integer ispec_selected_source,is_proc_source,nb_proc_source
  integer, intent(in)  :: nproc, myrank
  double precision xi_source,gamma_source

#ifdef USE_MPI
  integer, dimension(1:nproc)  :: allgather_is_proc_source
  integer, dimension(1)  :: locate_is_proc_source
  integer  :: ierror
#endif

! **************

  if (myrank == 0 .or. nproc == 1) then
    write(IOUT,*)
    write(IOUT,*) '*******************************'
    write(IOUT,*) ' locating moment-tensor source'
    write(IOUT,*) '*******************************'
    write(IOUT,*)
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
           if(dist_squared < distmin_squared) then
              distmin_squared = dist_squared
              ispec_selected_source = ispec
              ix_initial_guess = i
              iz_initial_guess = j
           endif

        enddo
     enddo

! end of loop on all the spectral elements
  enddo

#ifdef USE_MPI
  ! global minimum distance computed over all processes
  call MPI_ALLREDUCE (distmin_squared, dist_glob_squared, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
#else
  dist_glob_squared = distmin_squared
#endif

! check if this process contains the source
  if ( abs(sqrt(dist_glob_squared) - sqrt(distmin_squared)) < TINYVAL ) is_proc_source = 1

#ifdef USE_MPI
  ! determining the number of processes that contain the source
  ! (useful when the source is located on an interface)
  call MPI_ALLREDUCE (is_proc_source, nb_proc_source, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
#else
  nb_proc_source = is_proc_source
#endif


#ifdef USE_MPI
  ! when several processes contain the source, we elect one of them (minimum rank).
  if ( nb_proc_source > 1 ) then

     call MPI_ALLGATHER(is_proc_source, 1, MPI_INTEGER, allgather_is_proc_source(1), &
                        1, MPI_INTEGER, MPI_COMM_WORLD, ierror)
     locate_is_proc_source = maxloc(allgather_is_proc_source) - 1

     if ( myrank /= locate_is_proc_source(1) ) then
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
     write(IOUT,*)
     write(IOUT,*) 'Moment-tensor source:'

     if(final_distance == HUGEVAL) call exit_MPI('error locating moment-tensor source')

     write(IOUT,*) '            original x: ',sngl(x_source)
     write(IOUT,*) '            original z: ',sngl(z_source)
     write(IOUT,*) 'closest estimate found: ',sngl(final_distance),' m away'
#ifdef USE_MPI
     write(IOUT,*) ' in rank ',myrank
#endif
     write(IOUT,*) ' in element ',ispec_selected_source
     write(IOUT,*) ' at xi,gamma coordinates = ',xi_source,gamma_source
     write(IOUT,*)

     write(IOUT,*)
     write(IOUT,*) 'end of moment-tensor source detection'
     write(IOUT,*)
  endif

  end subroutine locate_source_moment_tensor

