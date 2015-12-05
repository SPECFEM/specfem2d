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

  use specfem_par

  implicit none

  ! generate the global numbering
  call setup_mesh_numbering()

  ! sets point coordinates
  call setup_mesh_coordinates()

  ! sets material properties on node points
  call setup_mesh_properties()

  ! for periodic edges
  call setup_mesh_periodic_edges()

  ! for acoustic forcing
  call setup_mesh_acoustic_forcing_edges()

  ! reads in external models and re-assigns material properties
  call setup_mesh_external_models()

  ! synchronizes all processes
  call synchronize_all()

  ! performs basic checks on parameters read
  all_anisotropic = .false.
  if (count(ispec_is_anisotropic(:) .eqv. .true.) == nspec) all_anisotropic = .true.

  if (all_anisotropic .and. anyabs) &
    call exit_MPI(myrank,'Cannot put absorbing boundaries if anisotropic materials along edges')

  if (ATTENUATION_VISCOELASTIC_SOLID .and. all_anisotropic) then
    call exit_MPI(myrank,'Cannot turn attenuation on in anisotropic materials')
  endif

  ! synchronizes all processes
  call synchronize_all()

  ! global domain flags
  ! (sets global flag for all slices)
  call any_all_l(any_elastic, ELASTIC_SIMULATION)
  call any_all_l(any_poroelastic, POROELASTIC_SIMULATION)
  call any_all_l(any_acoustic, ACOUSTIC_SIMULATION)
  call any_all_l(any_gravitoacoustic, GRAVITOACOUSTIC_SIMULATION)

  ! check for acoustic
  if (ATTENUATION_VISCOELASTIC_SOLID .and. .not. ELASTIC_SIMULATION) &
    call exit_MPI(myrank,'currently cannot have attenuation if acoustic/poroelastic simulation only')

  ! sets up domain coupling, i.e. edge detection for domain coupling
  call get_coupling_edges()

  ! synchronizes all processes
  call synchronize_all()

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
  integer :: i,n

  ! determines mesh dimensions
  xmin_local = minval(coord(1,:))
  xmax_local = maxval(coord(1,:))
  zmin_local = minval(coord(2,:))
  zmax_local = maxval(coord(2,:))

  ! collect min/max
  call min_all_all_dp(xmin_local, xmin)
  call max_all_all_dp(xmax_local, xmax)
  call min_all_all_dp(zmin_local, zmin)
  call max_all_all_dp(zmax_local, zmax)

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

!
!-----------------------------------------------------------------------------------
!

  subroutine setup_mesh_periodic_edges()

  use constants,only: IMAIN,NGLLX,NGLLZ,HUGEVAL
  use specfem_par

  implicit none

  ! local parameters
  integer :: ispec,i,j,iglob,iglob2,ier
  double precision :: xmaxval,xminval,ymaxval,yminval,xtol,xtypdist
  integer :: counter

! allocate an array to make sure that an acoustic free surface is not enforced on periodic edges
  allocate(this_ibool_is_a_periodic_edge(NGLOB),stat=ier)
  if (ier /= 0) stop 'Error allocating periodic edge array'

  this_ibool_is_a_periodic_edge(:) = .false.

! periodic conditions: detect common points between left and right edges and replace one of them with the other
  if (ADD_PERIODIC_CONDITIONS) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'implementing periodic boundary conditions'
      write(IMAIN,*) 'in the horizontal direction with a periodicity distance of ',PERIODIC_HORIZ_DIST,' m'
      if (PERIODIC_HORIZ_DIST <= 0.d0) stop 'PERIODIC_HORIZ_DIST should be greater than zero when using ADD_PERIODIC_CONDITIONS'
      write(IMAIN,*)
      write(IMAIN,*) '*****************************************************************'
      write(IMAIN,*) '*****************************************************************'
      write(IMAIN,*) '**** BEWARE: because of periodic conditions, values computed ****'
      write(IMAIN,*) '****         by check_grid() below will not be reliable       ****'
      write(IMAIN,*) '*****************************************************************'
      write(IMAIN,*) '*****************************************************************'
      write(IMAIN,*)
    endif

    ! set up a local geometric tolerance
    xtypdist = +HUGEVAL

    do ispec = 1,nspec

      xminval = +HUGEVAL
      yminval = +HUGEVAL
      xmaxval = -HUGEVAL
      ymaxval = -HUGEVAL

      ! only loop on the four corners of each element to get a typical size
      do j = 1,NGLLZ,NGLLZ-1
        do i = 1,NGLLX,NGLLX-1
          iglob = ibool(i,j,ispec)
          xmaxval = max(coord(1,iglob),xmaxval)
          xminval = min(coord(1,iglob),xminval)
          ymaxval = max(coord(2,iglob),ymaxval)
          yminval = min(coord(2,iglob),yminval)
        enddo
      enddo

      ! compute the minimum typical "size" of an element in the mesh
      xtypdist = min(xtypdist,xmaxval-xminval)
      xtypdist = min(xtypdist,ymaxval-yminval)

    enddo

    ! define a tolerance, small with respect to the minimum size
    xtol = 1.d-4 * xtypdist

! detect the points that are on the same horizontal line (i.e. at the same height Z)
! and that have a value of the horizontal coordinate X that differs by exactly the periodicity length;
! if so, make them all have the same global number, which will then implement periodic boundary conditions automatically.
! We select the smallest value of iglob and assign it to all the points that are the same due to periodicity,
! this way the maximum value of the ibool() array will remain as small as possible.
!
! *** IMPORTANT: this simple algorithm will be slow for large meshes because it has a cost of NGLOB^2 / 2
! (where NGLOB is the number of points per MPI slice, not of the whole mesh though). This could be
! reduced to O(NGLOB log(NGLOB)) by using a quicksort algorithm on the coordinates of the points to detect the multiples
! (as implemented in routine createnum_fast() elsewhere in the code). This could be done one day if needed instead
! of the very simple double loop below.
    if (myrank == 0) then
      write(IMAIN,*) 'start detecting points for periodic boundary conditions '// &
                     '(the current algorithm can be slow and could be improved)...'
    endif

    counter = 0
    do iglob = 1,NGLOB-1
      do iglob2 = iglob + 1,NGLOB
        ! check if the two points have the exact same Z coordinate
        if (abs(coord(2,iglob2) - coord(2,iglob)) < xtol) then
          ! if so, check if their X coordinate differs by exactly the periodicity distance
          if (abs(abs(coord(1,iglob2) - coord(1,iglob)) - PERIODIC_HORIZ_DIST) < xtol) then
            ! if so, they are the same point, thus replace the highest value of ibool with the lowest
            ! to make them the same global point and thus implement periodicity automatically
            counter = counter + 1
            this_ibool_is_a_periodic_edge(iglob) = .true.
            this_ibool_is_a_periodic_edge(iglob2) = .true.
            do ispec = 1,nspec
              do j = 1,NGLLZ
                do i = 1,NGLLX
                  if (ibool(i,j,ispec) == iglob2) ibool(i,j,ispec) = iglob
                enddo
              enddo
            enddo
          endif
        endif
      enddo
    enddo

    if (myrank == 0) write(IMAIN,*) 'done detecting points for periodic boundary conditions.'

    if (counter > 0) write(IMAIN,*) 'implemented periodic conditions on ',counter,' grid points on proc ',myrank

  endif ! of if (ADD_PERIODIC_CONDITIONS)

  end subroutine setup_mesh_periodic_edges

!
!-----------------------------------------------------------------------------------
!

  subroutine setup_mesh_acoustic_forcing_edges()

! acoustic forcing edge detection

  use specfem_par

  implicit none

  ! local parameters
  integer :: ipoin1D

  ! acoustic forcing edge detection
  ! the elements forming an edge are already known (computed in meshfem2D),
  ! the common nodes forming the edge are computed here
  if (ACOUSTIC_FORCING) then

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Acoustic forcing simulation'
      write(IMAIN,*)
      write(IMAIN,*) 'Beginning of acoustic forcing edge detection'
      call flush_IMAIN()
    endif

    ! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,IBOTTOM) = ipoin1D
      jvalue(ipoin1D,IBOTTOM) = NGLLZ
      jvalue_inverse(ipoin1D,IBOTTOM) = NGLLZ

      ivalue(ipoin1D,IRIGHT) = 1
      ivalue_inverse(ipoin1D,IRIGHT) = 1
      jvalue(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,IRIGHT) = ipoin1D

      ivalue(ipoin1D,ITOP) = ipoin1D
      ivalue_inverse(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,ITOP) = 1
      jvalue_inverse(ipoin1D,ITOP) = 1

      ivalue(ipoin1D,ILEFT) = NGLLX
      ivalue_inverse(ipoin1D,ILEFT) = NGLLX
      jvalue(ipoin1D,ILEFT) = ipoin1D
      jvalue_inverse(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1

    enddo

  endif ! if (ACOUSTIC_FORCING)

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_mesh_acoustic_forcing_edges


!
!-----------------------------------------------------------------------------------
!

  subroutine setup_mesh_external_models()

! external models

  use specfem_par

  implicit none

  ! local parameters
  integer :: nspec_ext,ier

  ! allocates material arrays
  if (assign_external_model) then
    nspec_ext = nspec
  else
    ! dummy allocations
    nspec_ext = 1
  endif

  allocate(vpext(NGLLX,NGLLZ,nspec_ext), &
           vsext(NGLLX,NGLLZ,nspec_ext), &
           rhoext(NGLLX,NGLLZ,nspec_ext), &
           gravityext(NGLLX,NGLLZ,nspec_ext), &
           Nsqext(NGLLX,NGLLZ,nspec_ext), &
           QKappa_attenuationext(NGLLX,NGLLZ,nspec_ext), &
           Qmu_attenuationext(NGLLX,NGLLZ,nspec_ext), &
           c11ext(NGLLX,NGLLZ,nspec_ext), &
           c13ext(NGLLX,NGLLZ,nspec_ext), &
           c15ext(NGLLX,NGLLZ,nspec_ext), &
           c33ext(NGLLX,NGLLZ,nspec_ext), &
           c35ext(NGLLX,NGLLZ,nspec_ext), &
           c55ext(NGLLX,NGLLZ,nspec_ext), &
           c12ext(NGLLX,NGLLZ,nspec_ext), &
           c23ext(NGLLX,NGLLZ,nspec_ext), &
           c25ext(NGLLX,NGLLZ,nspec_ext),stat=ier)
  if (ier /= 0) stop 'Error allocating external model arrays'

  ! reads in external models
  if (assign_external_model) then
    if (myrank == 0) then
      write(IMAIN,*) 'Assigning an external velocity and density model'
      call flush_IMAIN()
    endif
    call read_external_model()
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine setup_mesh_external_models

