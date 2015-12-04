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


  subroutine setup_mesh()

! creates mesh related properties, local to global mesh numbering and node locations

  implicit none

  ! generate the global numbering
  call setup_mesh_numbering()

  ! sets point coordinates
  call setup_mesh_coordinates()

  ! sets material properties on node points
  call setup_mesh_properties()

  end subroutine setup_mesh

!
!-----------------------------------------------------------------------------------
!

  subroutine setup_mesh_numbering()

  use specfem_par

  implicit none

  ! local parameters
  integer :: ier
  ! to count the number of degrees of freedom
  integer :: count_nspec_acoustic_total,nspec_total,nglob_total
  integer :: nb_acoustic_DOFs,nb_elastic_DOFs
  double precision :: ratio_1DOF,ratio_2DOFs

  ! "slow and clean" or "quick and dirty" version
  if (FAST_NUMBERING) then
    call createnum_fast()
  else
    call createnum_slow()
  endif

  ! gets total numbers for all slices
  call sum_all_i(count_nspec_acoustic,count_nspec_acoustic_total)
  call sum_all_i(nspec,nspec_total)
  call sum_all_i(nglob,nglob_total)

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of elements: ',nspec_total
    write(IMAIN,*) 'decomposed as follows:'
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of elastic/visco/poro elements: ',nspec_total - count_nspec_acoustic_total
    write(IMAIN,*) 'Total number of acoustic elements: ',count_nspec_acoustic_total
    write(IMAIN,*)
#ifdef USE_MPI
    write(IMAIN,*) 'Approximate total number of grid points in the mesh'
    write(IMAIN,*) '(with a few duplicates coming from MPI buffers): ',nglob_total
#else
    write(IMAIN,*) 'Exact total number of grid points in the mesh: ',nglob_total
#endif

    ! percentage of elements with 2 degrees of freedom per point
    ratio_2DOFs = (nspec_total - count_nspec_acoustic_total) / dble(nspec_total)
    ratio_1DOF  = count_nspec_acoustic_total / dble(nspec_total)

    nb_acoustic_DOFs = nint(nglob_total*ratio_1DOF)

    ! elastic elements have two degrees of freedom per point
    nb_elastic_DOFs  = nint(nglob_total*ratio_2DOFs*2)

    if (p_sv) then
      write(IMAIN,*)
      write(IMAIN,*) 'Approximate number of acoustic degrees of freedom in the mesh: ',nb_acoustic_DOFs
      write(IMAIN,*) 'Approximate number of elastic degrees of freedom in the mesh: ',nb_elastic_DOFs
      write(IMAIN,*) '  (there are 2 degrees of freedom per point for elastic elements)'
      write(IMAIN,*)
      write(IMAIN,*) 'Approximate total number of degrees of freedom in the mesh'
      write(IMAIN,*) '(sum of the two values above): ',nb_acoustic_DOFs + nb_elastic_DOFs
      write(IMAIN,*)
      write(IMAIN,*) ' (for simplicity viscoelastic or poroelastic elements, if any,'
      write(IMAIN,*) '  are counted as elastic in the above three estimates;'
      write(IMAIN,*) '  in reality they have more degrees of freedom)'
      write(IMAIN,*)
    endif
    call flush_IMAIN()
  endif

  ! allocate temporary arrays
  allocate(integer_mask_ibool(nglob),stat=ier)
  if (ier /= 0 ) stop 'error allocating integer_mask_ibool'
  allocate(copy_ibool_ori(NGLLX,NGLLZ,nspec),stat=ier)
  if (ier /= 0 ) stop 'error allocating copy_ibool_ori'

  ! reduce cache misses by sorting the global numbering in the order in which it is accessed in the time loop.
  ! this speeds up the calculations significantly on modern processors
  call get_global()

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_mesh_numbering

!
!-----------------------------------------------------------------------------------
!


  subroutine setup_mesh_coordinates()

  use specfem_par

  implicit none

  ! local parameters
  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl
  double precision :: xi,gamma,x,z

  integer :: i,j,ispec,iglob,ier

  ! to help locate elements with a negative Jacobian using OpenDX
  logical :: found_a_negative_jacobian

  ! allocate other global arrays
  allocate(coord(NDIM,nglob),stat=ier)
  if (ier /= 0) stop 'Error allocating coord array'

  ! sets the coordinates of the points of the global grid
  found_a_negative_jacobian = .false.
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            xi = xiglj(i)
          else
            xi = xigll(i)
          endif
        else
          xi = xigll(i)
        endif
        gamma = zigll(j)

        call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl, &
                                jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo, &
                                .false.)

        if (jacobianl <= ZERO) found_a_negative_jacobian = .true.

        ! coordinates of global nodes
        iglob = ibool(i,j,ispec)
        coord(1,iglob) = x
        coord(2,iglob) = z

        xix(i,j,ispec) = xixl
        xiz(i,j,ispec) = xizl
        gammax(i,j,ispec) = gammaxl
        gammaz(i,j,ispec) = gammazl
        jacobian(i,j,ispec) = jacobianl

      enddo
    enddo
  enddo

! create an OpenDX file containing all the negative elements displayed in red, if any
! this allows users to locate problems in a mesh based on the OpenDX file created at the second iteration
! do not create OpenDX files if no negative Jacobian has been found, or if we are running in parallel
! (because writing OpenDX routines is much easier in serial)
  if (found_a_negative_jacobian .and. nproc == 1) then
    call save_openDX_jacobian(nspec,npgeo,ngnod,knods,coorg,xigll,zigll,AXISYM,is_on_the_axis,xiglj)
  endif

  ! stop the code at the first negative element found, because such a mesh cannot be computed
  if (found_a_negative_jacobian) then
    do ispec = 1,nspec
      do j = 1,NGLLZ
        do i = 1,NGLLX
          if (AXISYM) then
            if (is_on_the_axis(ispec)) then
              xi = xiglj(i)
            else
              xi = xigll(i)
            endif
          else
            xi = xigll(i)
          endif
          gamma = zigll(j)

          call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl, &
                                  jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo, &
                                  .true.)
        enddo
      enddo
    enddo
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_mesh_coordinates

!
!-----------------------------------------------------------------------------------
!

  subroutine setup_mesh_properties()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  double precision :: xmin,xmax,zmin,zmax
  double precision :: xmin_local,xmax_local,zmin_local,zmax_local
  integer :: i,n,ier

  ! determines mesh dimensions
  xmin_local = minval(coord(1,:))
  xmax_local = maxval(coord(1,:))
  zmin_local = minval(coord(2,:))
  zmax_local = maxval(coord(2,:))

  ! collect min/max
#ifdef USE_MPI
  call MPI_ALLREDUCE(xmin_local, xmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE(xmax_local, xmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE(zmin_local, zmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE(zmax_local, zmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
#else
  xmin = xmin_local
  xmax = xmax_local
  zmin = zmin_local
  zmax = zmax_local
#endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Xmin,Xmax of the whole mesh = ',xmin,xmax
    write(IMAIN,*) 'Zmin,Zmax of the whole mesh = ',zmin,zmax
    write(IMAIN,*)
  endif

  ! checks that no source is located outside the mesh
  if (myrank == 0) then
    do i = 1,NSOURCES
      if (x_source(i) < xmin) stop 'error: at least one source has x < xmin of the mesh'
      if (x_source(i) > xmax) stop 'error: at least one source has x > xmax of the mesh'

      if (z_source(i) < zmin) stop 'error: at least one source has z < zmin of the mesh'
      if (z_source(i) > zmax) stop 'error: at least one source has z > zmax of the mesh'
    enddo
  endif

  ! use a spring to improve the stability of the Stacey condition
  if(STACEY_BOUNDARY_CONDITIONS .and. ADD_SPRING_TO_STACEY) then
    x_center_spring = (xmax + xmin)/2.d0
    z_center_spring = (zmax + zmin)/2.d0
  endif

  ! saves the grid of points in a file
  if (output_grid_ASCII .and. myrank == 0) then
     write(IMAIN,*)
     write(IMAIN,*) 'Saving the grid in an ASCII text file...'
     write(IMAIN,*)
     open(unit=55,file='OUTPUT_FILES/ASCII_dump_of_grid_points.txt',status='unknown')
     write(55,*) nglob
     do n = 1,nglob
        write(55,*) (coord(i,n), i = 1,NDIM)
     enddo
     close(55)
  endif

  ! plots the GLL mesh in a Gnuplot file
  if (output_grid_Gnuplot .and. myrank == 0) then
    call plot_gll()
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_mesh_properties

