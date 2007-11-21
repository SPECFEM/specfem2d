
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!  Main authors: Dimitri Komatitsch, Nicolas Le Goff and Roland Martin
!                 University of Pau and CNRS, France
!
!        (c) University of Pau and CNRS, France, November 2007
!
!========================================================================

!----
!---- locate_source_moment_tensor finds the correct position of the moment-tensor source
!----

  subroutine locate_source_moment_tensor(ibool,coord,nspec,npoin,xigll,zigll,x_source,z_source, &
               ispec_selected_source,is_proc_source,nb_proc_source,nproc,myrank,xi_source,gamma_source,coorg,knods,ngnod,npgeo)

  implicit none

  include "constants.h"
#ifdef USE_MPI
  include "mpif.h"
#endif

  integer nspec,npoin,ngnod,npgeo

  integer knods(ngnod,nspec)
  double precision coorg(NDIM,npgeo)

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

! array containing coordinates of the points
  double precision coord(NDIM,npoin)

  integer i,j,ispec,iglob,iter_loop,ix_initial_guess,iz_initial_guess

  double precision x_source,z_source,dist
  double precision xi,gamma,dx,dz,dxi,dgamma

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX)
  double precision zigll(NGLLZ)

  double precision x,z,xix,xiz,gammax,gammaz,jacobian
  double precision distmin,final_distance,dist_glob

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
  if ( myrank == 0 .or. nproc == 1 ) then
  write(IOUT,*)
  write(IOUT,*) '*******************************'
  write(IOUT,*) ' locating moment-tensor source'
  write(IOUT,*) '*******************************'
  write(IOUT,*)
  end if

! set distance to huge initial value
  distmin=HUGEVAL

  is_proc_source = 0

  do ispec=1,nspec

! loop only on points inside the element
! exclude edges to ensure this point is not shared with other elements
     do j=2,NGLLZ-1
        do i=2,NGLLX-1

           iglob = ibool(i,j,ispec)
           dist = sqrt((x_source-dble(coord(1,iglob)))**2 + (z_source-dble(coord(2,iglob)))**2)

!          keep this point if it is closer to the source
           if(dist < distmin) then
              distmin = dist
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
  call MPI_ALLREDUCE (distmin, dist_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)

#else
  dist_glob = distmin

#endif


  ! check if this process contains the source
  if ( dist_glob == distmin ) then
     is_proc_source = 1
  end if


#ifdef USE_MPI
  ! determining the number of processes that contain the source (useful when the source is located on an interface)
  call MPI_ALLREDUCE (is_proc_source, nb_proc_source, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)

#else
  nb_proc_source = is_proc_source

#endif


#ifdef USE_MPI
  ! when several processes contain the source, we elect one of them (minimum rank).
  if ( nb_proc_source > 1 ) then

     call MPI_ALLGATHER(is_proc_source, 1, MPI_INTEGER, allgather_is_proc_source(1), 1, MPI_INTEGER, MPI_COMM_WORLD, ierror)
     locate_is_proc_source = maxloc(allgather_is_proc_source) - 1

     if ( myrank /= locate_is_proc_source(1) ) then
        is_proc_source = 0
     end if
     nb_proc_source = 1

  end if

#endif

! ****************************************
! find the best (xi,gamma) for each source
! ****************************************

! use initial guess in xi and gamma
        xi = xigll(ix_initial_guess)
        gamma = zigll(iz_initial_guess)

! iterate to solve the non linear system
  do iter_loop = 1,NUM_ITER

! recompute jacobian for the new point
    call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian,coorg,knods,ispec_selected_source,ngnod,nspec,npgeo)

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
    call recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian,coorg,knods,ispec_selected_source,ngnod,nspec,npgeo)

! store xi,gamma of point found
  xi_source = xi
  gamma_source = gamma

! compute final distance between asked and found
  final_distance = sqrt((x_source-x)**2 + (z_source-z)**2)

  if ( is_proc_source == 1 ) then
     write(IOUT,*)
     write(IOUT,*) 'Moment-tensor source:'

     if(final_distance == HUGEVAL) call exit_MPI('error locating moment-tensor source')

     write(IOUT,*) '            original x: ',sngl(x_source)
     write(IOUT,*) '            original z: ',sngl(z_source)
     write(IOUT,*) 'closest estimate found: ',sngl(final_distance),' m away'
     write(IOUT,*) ' in element ',ispec_selected_source
     write(IOUT,*) ' at xi,gamma coordinates = ',xi_source,gamma_source
     write(IOUT,*)

     write(IOUT,*)
     write(IOUT,*) 'end of moment-tensor source detection'
     write(IOUT,*)
  end if

#ifdef USE_MPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
#endif

  end subroutine locate_source_moment_tensor

