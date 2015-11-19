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


subroutine prepare_timerun()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

  integer :: i,j,ispec,k,iglob,irec,i_source,ispecabs,irecloc,ier

#ifdef USE_MPI
  include "precision.h"
#endif

  ! reads sources, stations, and mesh from database
  call prepare_timerun_read()


!
!---- compute shape functions and their derivatives for SEM grid
!

! set up Gauss-Lobatto-Legendre points, weights and also derivation matrices
  call define_derivation_matrices()

  if (AXISYM) then
    ! set up Gauss-Lobatto-Jacobi points, weights and also derivation matrices
    call define_GLJ_derivation_matrix()
  endif

  do j = 1,NGLLZ
    do i = 1,NGLLX
      call define_shape_functions(shape2D(:,i,j),dershape2D(:,:,i,j),xigll(i),zigll(j),ngnod)
    enddo
  enddo

!
!---- generate the global numbering
!

! "slow and clean" or "quick and dirty" version
  if (FAST_NUMBERING) then
    call createnum_fast()
  else
    call createnum_slow()
  endif

#ifdef USE_MPI
  call MPI_REDUCE(count_nspec_acoustic, count_nspec_acoustic_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
  call MPI_REDUCE(nspec, nspec_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
  call MPI_REDUCE(nglob, nglob_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
#else
  count_nspec_acoustic_total = count_nspec_acoustic
  nspec_total = nspec
  nglob_total = nglob
#endif
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
    if (ier /= 0 ) stop 'error allocating mask_ibool'
    allocate(copy_ibool_ori(NGLLX,NGLLZ,nspec),stat=ier)
    if (ier /= 0 ) stop 'error allocating copy_ibool_ori'

    ! reduce cache misses by sorting the global numbering in the order in which it is accessed in the time loop.
    ! this speeds up the calculations significantly on modern processors
    call get_global()

!---- compute shape functions and their derivatives for regular interpolated display grid
  do j = 1,pointsdisp
    do i = 1,pointsdisp
      xirec  = 2.d0*dble(i-1)/dble(pointsdisp-1) - 1.d0
      gammarec  = 2.d0*dble(j-1)/dble(pointsdisp-1) - 1.d0
      call define_shape_functions(shape2D_display(:,i,j),dershape2D_display(:,:,i,j),xirec,gammarec,ngnod)
    enddo
  enddo

!---- compute Lagrange interpolants on a regular interpolated grid in (xi,gamma)
!---- for display (assumes NGLLX = NGLLZ)
  do j = 1,NGLLX
    do i = 1,pointsdisp
      xirec  = 2.d0*dble(i-1)/dble(pointsdisp-1) - 1.d0
      flagrange(j,i) = hgll(j-1,xirec,xigll,NGLLX)
      if (AXISYM) flagrange_GLJ(j,i) = hglj(j-1,xirec,xiglj,NGLJ)
    enddo
  enddo

! get number of stations from receiver file
  open(unit=IIN,file='DATA/STATIONS',iostat=ios,status='old',action='read')
  nrec = 0
  do while(ios == 0)
    read(IIN,"(a)",iostat=ios) dummystring
    if (ios == 0) nrec = nrec + 1
  enddo
  close(IIN)

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of receivers = ',nrec
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (nrec < 1) call exit_MPI('need at least one receiver')

! receiver information
    allocate(ispec_selected_rec(nrec))
    allocate(st_xval(nrec))
    allocate(st_zval(nrec))
    allocate(xi_receiver(nrec))
    allocate(gamma_receiver(nrec))
    allocate(station_name(nrec))
    allocate(network_name(nrec))
    allocate(recloc(nrec))
    allocate(which_proc_receiver(nrec))
    allocate(x_final_receiver(nrec))
    allocate(z_final_receiver(nrec))
    allocate(ix_image_color_receiver(nrec))
    allocate(iy_image_color_receiver(nrec))

! allocate 1-D Lagrange interpolators and derivatives
    allocate(hxir(NGLLX))
    allocate(hxis(NGLLX))
    allocate(hpxir(NGLLX))
    allocate(hpxis(NGLLX))
    allocate(hgammar(NGLLZ))
    allocate(hgammas(NGLLZ))
    allocate(hpgammar(NGLLZ))
    allocate(hpgammas(NGLLZ))

! allocate Lagrange interpolators for receivers
    allocate(hxir_store(nrec,NGLLX))
    allocate(hgammar_store(nrec,NGLLZ))

! allocate Lagrange interpolators for sources
    allocate(hxis_store(NSOURCES,NGLLX))
    allocate(hgammas_store(NSOURCES,NGLLZ))

! allocate other global arrays
    allocate(coord(NDIM,nglob))

! to display the whole vector field (it needs to be computed from the potential in acoustic elements,
! thus it does not exist as a whole it case of simulations that contain some acoustic elements
! and it thus needs to be computed specifically for display purposes)
    allocate(vector_field_display(3,nglob))

! when periodic boundary conditions are on, some global degrees of freedom are going to be removed,
! thus we need to set this array to zero otherwise some of its locations may contain random values
! if the memory is not cleaned
    vector_field_display(:,:) = 0.d0

    if (assign_external_model) then
      allocate(vpext(NGLLX,NGLLZ,nspec))
      allocate(vsext(NGLLX,NGLLZ,nspec))
      allocate(rhoext(NGLLX,NGLLZ,nspec))
      allocate(gravityext(NGLLX,NGLLZ,nspec))
      allocate(Nsqext(NGLLX,NGLLZ,nspec))
      allocate(QKappa_attenuationext(NGLLX,NGLLZ,nspec))
      allocate(Qmu_attenuationext(NGLLX,NGLLZ,nspec))
      allocate(c11ext(NGLLX,NGLLZ,nspec))
      allocate(c13ext(NGLLX,NGLLZ,nspec))
      allocate(c15ext(NGLLX,NGLLZ,nspec))
      allocate(c33ext(NGLLX,NGLLZ,nspec))
      allocate(c35ext(NGLLX,NGLLZ,nspec))
      allocate(c55ext(NGLLX,NGLLZ,nspec))
      allocate(c12ext(NGLLX,NGLLZ,nspec))
      allocate(c23ext(NGLLX,NGLLZ,nspec))
      allocate(c25ext(NGLLX,NGLLZ,nspec))
    else
      allocate(vpext(1,1,1))
      allocate(vsext(1,1,1))
      allocate(rhoext(1,1,1))
      allocate(gravityext(1,1,1))
      allocate(Nsqext(1,1,1))
      allocate(QKappa_attenuationext(1,1,1))
      allocate(Qmu_attenuationext(1,1,1))
      allocate(c11ext(1,1,1))
      allocate(c13ext(1,1,1))
      allocate(c15ext(1,1,1))
      allocate(c33ext(1,1,1))
      allocate(c35ext(1,1,1))
      allocate(c55ext(1,1,1))
      allocate(c12ext(1,1,1))
      allocate(c23ext(1,1,1))
      allocate(c25ext(1,1,1))
    endif

!
!----  set the coordinates of the points of the global grid
!
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

          coord(1,ibool(i,j,ispec)) = x
          coord(2,ibool(i,j,ispec)) = z

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

  xmin_local = minval(coord(1,:))
  xmax_local = maxval(coord(1,:))
  zmin_local = minval(coord(2,:))
  zmax_local = maxval(coord(2,:))

#ifdef USE_MPI
  call MPI_REDUCE(xmin_local, xmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ier)
  call MPI_REDUCE(xmax_local, xmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ier)
  call MPI_REDUCE(zmin_local, zmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ier)
  call MPI_REDUCE(zmax_local, zmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ier)
#else
  xmin = xmin_local
  xmax = xmax_local
  zmin = zmin_local
  zmax = zmax_local
#endif

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Xmin,Xmax of the whole mesh = ',xmin,xmax
    write(IMAIN,*) 'Zmin,Zmax of the whole mesh = ',zmin,zmax
    write(IMAIN,*)

! check that no source is located outside the mesh
    do i = 1,NSOURCES
      if (x_source(i) < xmin) stop 'error: at least one source has x < xmin of the mesh'
      if (x_source(i) > xmax) stop 'error: at least one source has x > xmax of the mesh'

      if (z_source(i) < zmin) stop 'error: at least one source has z < zmin of the mesh'
      if (z_source(i) > zmax) stop 'error: at least one source has z > zmax of the mesh'
    enddo

  endif

! use a spring to improve the stability of the Stacey condition
  x_center_spring = (xmax + xmin)/2.d0
  z_center_spring = (zmax + zmin)/2.d0

! allocate an array to make sure that an acoustic free surface is not enforced on periodic edges
  allocate(this_ibool_is_a_periodic_edge(NGLOB))
  this_ibool_is_a_periodic_edge(:) = .false.

! periodic conditions: detect common points between left and right edges and replace one of them with the other
    if (ADD_PERIODIC_CONDITIONS) then

      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'implementing periodic boundary conditions'
        write(IMAIN,*) 'in the horizontal direction with a periodicity distance of ',PERIODIC_HORIZ_DIST,' m'
        if (PERIODIC_HORIZ_DIST <= 0.d0) stop 'PERIODIC_HORIZ_DIST should be greater than zero when using ADD_PERIODIC_CONDITIONS'
        write(IMAIN,*)
        write(IMAIN,*) '*****************************************************************'
        write(IMAIN,*) '*****************************************************************'
        write(IMAIN,*) '**** BEWARE: because of periodic conditions, values computed ****'
        write(IMAIN,*) '****         by checkgrid() below will not be reliable       ****'
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
      if (myrank == 0) write(IMAIN,*) &
        'start detecting points for periodic boundary conditions (the current algorithm can be slow and could be improved)...'
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

!
!--- save the grid of points in a file
!
  if (output_grid_ASCII .and. myrank == 0) then
     write(IMAIN,*)
     write(IMAIN,*) 'Saving the grid in an ASCII text file...'
     write(IMAIN,*)
     open(unit=55,file='OUTPUT_FILES/ASCII_dump_of_grid_points.txt',status='unknown')
     write(55,*) nglob
     do n = 1,nglob
        write(55,*) (coord(i,n), i= 1,NDIM)
     enddo
     close(55)
  endif

!
!-----   plot the GLL mesh in a Gnuplot file
!
  if (output_grid_Gnuplot .and. myrank == 0)  &
    call plotgll()

  if (assign_external_model) then
    if (myrank == 0) write(IMAIN,*) 'Assigning an external velocity and density model...'
    call read_external_model()
  endif

!
!----  perform basic checks on parameters read
!
  all_anisotropic = .false.
  if (count(anisotropic(:) .eqv. .true.) == nspec) all_anisotropic = .true.

  if (all_anisotropic .and. anyabs) &
    call exit_MPI('Cannot put absorbing boundaries if anisotropic materials along edges')

  if (ATTENUATION_VISCOELASTIC_SOLID .and. all_anisotropic) then
    call exit_MPI('Cannot turn attenuation on in anisotropic materials')
  endif

  ! global domain flags
  any_elastic_glob = any_elastic
#ifdef USE_MPI
  call MPI_ALLREDUCE(any_elastic, any_elastic_glob, 1, MPI_LOGICAL, &
                    MPI_LOR, MPI_COMM_WORLD, ier)
#endif

  any_poroelastic_glob = any_poroelastic
#ifdef USE_MPI
  call MPI_ALLREDUCE(any_poroelastic, any_poroelastic_glob, 1, MPI_LOGICAL, &
                    MPI_LOR, MPI_COMM_WORLD, ier)
#endif

  any_acoustic_glob = any_acoustic
#ifdef USE_MPI
  call MPI_ALLREDUCE(any_acoustic, any_acoustic_glob, 1, MPI_LOGICAL, &
                    MPI_LOR, MPI_COMM_WORLD, ier)
#endif

   any_gravitoacoustic_glob = any_gravitoacoustic
#ifdef USE_MPI
  call MPI_ALLREDUCE(any_gravitoacoustic, any_gravitoacoustic_glob, 1, MPI_LOGICAL, &
                    MPI_LOR, MPI_COMM_WORLD, ier)
#endif

  ! for acoustic
  if (ATTENUATION_VISCOELASTIC_SOLID .and. .not. any_elastic_glob) &
    call exit_MPI('currently cannot have attenuation if acoustic/poroelastic simulation only')

!
!----   define coefficients of the Newmark time scheme
!
  deltatover2 = HALF*deltat
  deltatsquareover2 = HALF*deltat*deltat

  if (SIMULATION_TYPE == 3) then
!  define coefficients of the Newmark time scheme for the backward wavefield
    b_deltat = - deltat
    b_deltatover2 = HALF*b_deltat
    b_deltatsquareover2 = HALF*b_deltat*b_deltat
  endif

!---- define actual location of source and receivers

  call setup_sources_receivers()

! compute source array for adjoint source
  nadj_rec_local = 0
  if (SIMULATION_TYPE == 3) then  ! adjoint calculation

  allocate(source_adjointe(nrecloc,NSTEP,2))

    do irec = 1,nrec
      if (myrank == which_proc_receiver(irec)) then
        ! check that the source proc number is okay
        if (which_proc_receiver(irec) < 0 .or. which_proc_receiver(irec) > NPROC-1) &
              call exit_MPI('something is wrong with the source proc number in adjoint simulation')
        nadj_rec_local = nadj_rec_local + 1
      endif
    enddo
    allocate(adj_sourcearray(NSTEP,3,NGLLX,NGLLZ))
    if (nadj_rec_local > 0)  then
      allocate(adj_sourcearrays(nadj_rec_local,NSTEP,3,NGLLX,NGLLZ))
    else
      allocate(adj_sourcearrays(1,1,1,1,1))
    endif

  if (seismotype == 1 .or. seismotype == 2 .or. seismotype == 3) then

    if (.not. SU_FORMAT) then
       irec_local = 0
       allocate(adj_src_s(NSTEP,3))
       do irec = 1, nrec
         ! compute only adjoint source arrays in the local proc
         if (myrank == which_proc_receiver(irec)) then
           irec_local = irec_local + 1
           adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
           call compute_arrays_adj_source(xi_receiver(irec), gamma_receiver(irec))
           adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
         endif
       enddo
    else ! (SU_FORMAT)
        call add_adjoint_sources_SU()
    endif

  else if (seismotype == 4 .or. seismotype == 6) then

    if (.not. SU_FORMAT) then
       irec_local = 0
       allocate(adj_src_s(NSTEP,3))
       do irec = 1, nrec
         ! compute only adjoint source arrays in the local proc
         if (myrank == which_proc_receiver(irec)) then
           irec_local = irec_local + 1
           adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
           call compute_arrays_adj_source(xi_receiver(irec), gamma_receiver(irec))
           adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
         endif
       enddo
    else
       irec_local = 0
       write(filename, "('./SEM/Up_file_single.su.adj')")
       open(111,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ios)
               if (ios /= 0) call exit_MPI(' file '//trim(filename)//'does not exist')
       allocate(adj_src_s(NSTEP,3))

       do irec = 1, nrec
         if (myrank == which_proc_receiver(irec)) then
          irec_local = irec_local + 1
          adj_sourcearray(:,:,:,:) = 0.0
          read(111,rec=irec,iostat=ios) r4head, adj_src_s(:,1)
               if (ios /= 0) call exit_MPI(' file '//trim(filename)//' read error')
          if (irec==1) print *, r4head(1),r4head(19),r4head(20),r4head(21),r4head(22),header2(2)

          if (AXISYM) then
            if (is_on_the_axis(ispec_selected_rec(irec))) then
              call lagrange_any(xi_receiver(irec),NGLJ,xiglj,hxir,hpxir)
              !do j = 1,NGLJ ! Same result with that loop
              !  hxir(j) = hglj(j-1,xi_receiver(irec),xiglj,NGLJ)
              !enddo
            else
              call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
            endif
          else
            call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
          endif

          call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
          source_adjointe(irec_local,:,1) = adj_src_s(:,1)

          if (.not. GPU_MODE) then
            do k = 1, NGLLZ
              do i = 1, NGLLX
                adj_sourcearray(:,:,i,k) = hxir(i) * hgammar(k) * adj_src_s(:,:)
              enddo
            enddo
            adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
          endif

         endif !  if (myrank == which_proc_receiver(irec))
       enddo ! irec
       close(111)
       deallocate(adj_src_s)

    endif ! SU_format

  endif ! seismotype

  else
     allocate(adj_sourcearrays(1,1,1,1,1))
  endif ! SIMULATION_TYPE

    if (nrecloc > 0) then
      allocate(anglerec_irec(nrecloc))
      allocate(cosrot_irec(nrecloc))
      allocate(sinrot_irec(nrecloc))
      allocate(rec_tangential_detection_curve(nrecloc))
    else
      allocate(anglerec_irec(1))
      allocate(cosrot_irec(1))
      allocate(sinrot_irec(1))
      allocate(rec_tangential_detection_curve(1))
    endif

    if (rec_normal_to_surface .and. abs(anglerec) > 1.d-6) &
      stop 'anglerec should be zero when receivers are normal to the topography'

    anglerec_irec(:) = anglerec * pi / 180.d0
    cosrot_irec(:) = cos(anglerec_irec(:))
    sinrot_irec(:) = sin(anglerec_irec(:))

!
!--- tangential computation
!

! for receivers
    if (rec_normal_to_surface) then
      irecloc = 0
      do irec = 1, nrec
        if (which_proc_receiver(irec) == myrank) then
          irecloc = irecloc + 1
          distmin = HUGEVAL
          do i = 1, nnodes_tangential_curve
            dist_current = sqrt((x_final_receiver(irec)-nodes_tangential_curve(1,i))**2 + &
               (z_final_receiver(irec)-nodes_tangential_curve(2,i))**2)
            if (dist_current < distmin) then
              n1_tangential_detection_curve = i
              distmin = dist_current
            endif
         enddo

         rec_tangential_detection_curve(irecloc) = n1_tangential_detection_curve
         call tri_quad(n_tangential_detection_curve, n1_tangential_detection_curve, &
                      nnodes_tangential_curve)

         call compute_normal_vector( anglerec_irec(irecloc), &
           nodes_tangential_curve(1,n_tangential_detection_curve(1)), &
           nodes_tangential_curve(1,n_tangential_detection_curve(2)), &
           nodes_tangential_curve(1,n_tangential_detection_curve(3)), &
           nodes_tangential_curve(1,n_tangential_detection_curve(4)), &
           nodes_tangential_curve(2,n_tangential_detection_curve(1)), &
           nodes_tangential_curve(2,n_tangential_detection_curve(2)), &
           nodes_tangential_curve(2,n_tangential_detection_curve(3)), &
           nodes_tangential_curve(2,n_tangential_detection_curve(4)) )
       endif

      enddo
      cosrot_irec(:) = cos(anglerec_irec(:))
      sinrot_irec(:) = sin(anglerec_irec(:))
    endif

! for the source
    if (force_normal_to_surface) then

      do i_source= 1,NSOURCES
        if (is_proc_source(i_source) == 1) then
          distmin = HUGEVAL
          do i = 1, nnodes_tangential_curve
            dist_current = sqrt((coord(1,iglob_source(i_source))-nodes_tangential_curve(1,i))**2 + &
                                (coord(2,iglob_source(i_source))-nodes_tangential_curve(2,i))**2)
            if (dist_current < distmin) then
              n1_tangential_detection_curve = i
              distmin = dist_current

            endif
          enddo

          call tri_quad(n_tangential_detection_curve, n1_tangential_detection_curve, &
                       nnodes_tangential_curve)

          ! in the case of a source force vector
          ! users can give an angle with respect to the normal to the topography surface,
          ! in which case we must compute the normal to the topography
          ! and add it the existing rotation angle
          call compute_normal_vector( anglesource(i_source), &
                            nodes_tangential_curve(1,n_tangential_detection_curve(1)), &
                            nodes_tangential_curve(1,n_tangential_detection_curve(2)), &
                            nodes_tangential_curve(1,n_tangential_detection_curve(3)), &
                            nodes_tangential_curve(1,n_tangential_detection_curve(4)), &
                            nodes_tangential_curve(2,n_tangential_detection_curve(1)), &
                            nodes_tangential_curve(2,n_tangential_detection_curve(2)), &
                            nodes_tangential_curve(2,n_tangential_detection_curve(3)), &
                            nodes_tangential_curve(2,n_tangential_detection_curve(4)) )

          source_courbe_eros(i_source) = n1_tangential_detection_curve
          if (myrank == 0 .and. is_proc_source(i_source) == 1 .and. nb_proc_source(i_source) == 1) then
            source_courbe_eros(i_source) = n1_tangential_detection_curve
            anglesource_recv = anglesource(i_source)
#ifdef USE_MPI
          else if (myrank == 0) then
            do i = 1, nb_proc_source(i_source) - is_proc_source(i_source)
              call MPI_recv(source_courbe_eros(i_source),1,MPI_INTEGER, &
                          MPI_ANY_SOURCE,42,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
              call MPI_recv(anglesource_recv,1,MPI_DOUBLE_PRECISION, &
                          MPI_ANY_SOURCE,43,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
            enddo
          else if (is_proc_source(i_source) == 1) then
            call MPI_send(n1_tangential_detection_curve,1,MPI_INTEGER,0,42,MPI_COMM_WORLD,ier)
            call MPI_send(anglesource(i_source),1,MPI_DOUBLE_PRECISION,0,43,MPI_COMM_WORLD,ier)
#endif
          endif

#ifdef USE_MPI
          call MPI_bcast(anglesource_recv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
          anglesource(i_source) = anglesource_recv
#endif
        endif !  if (is_proc_source(i_source) == 1)
      enddo ! do i_source= 1,NSOURCES
    endif !  if (force_normal_to_surface)

! CHRIS --- how to deal with multiple source. Use first source now. ---
! compute distance from source to receivers following the curve
    if (force_normal_to_surface .and. rec_normal_to_surface) then
      dist_tangential_detection_curve(source_courbe_eros(1)) = 0
      do i = source_courbe_eros(1)+1, nnodes_tangential_curve
        dist_tangential_detection_curve(i) = dist_tangential_detection_curve(i-1) + &
            sqrt((nodes_tangential_curve(1,i)-nodes_tangential_curve(1,i-1))**2 + &
            (nodes_tangential_curve(2,i)-nodes_tangential_curve(2,i-1))**2)
      enddo
      dist_tangential_detection_curve(1) = dist_tangential_detection_curve(nnodes_tangential_curve) + &
           sqrt((nodes_tangential_curve(1,1)-nodes_tangential_curve(1,nnodes_tangential_curve))**2 + &
           (nodes_tangential_curve(2,1)-nodes_tangential_curve(2,nnodes_tangential_curve))**2)
      do i = 2, source_courbe_eros(1)-1
        dist_tangential_detection_curve(i) = dist_tangential_detection_curve(i-1) + &
            sqrt((nodes_tangential_curve(1,i)-nodes_tangential_curve(1,i-1))**2 + &
            (nodes_tangential_curve(2,i)-nodes_tangential_curve(2,i-1))**2)
      enddo
      do i = source_courbe_eros(1)-1, 1, -1
        dist_current = dist_tangential_detection_curve(i+1) + &
            sqrt((nodes_tangential_curve(1,i)-nodes_tangential_curve(1,i+1))**2 + &
            (nodes_tangential_curve(2,i)-nodes_tangential_curve(2,i+1))**2)
        if (dist_current < dist_tangential_detection_curve(i)) then
          dist_tangential_detection_curve(i) = dist_current
        endif
      enddo
      dist_current = dist_tangential_detection_curve(1) + &
         sqrt((nodes_tangential_curve(1,1)-nodes_tangential_curve(1,nnodes_tangential_curve))**2 + &
         (nodes_tangential_curve(2,1)-nodes_tangential_curve(2,nnodes_tangential_curve))**2)
      if (dist_current < dist_tangential_detection_curve(nnodes_tangential_curve)) then
        dist_tangential_detection_curve(nnodes_tangential_curve) = dist_current
      endif
      do i = nnodes_tangential_curve-1, source_courbe_eros(1)+1, -1
        dist_current = dist_tangential_detection_curve(i+1) + &
            sqrt((nodes_tangential_curve(1,i)-nodes_tangential_curve(1,i+1))**2 + &
            (nodes_tangential_curve(2,i)-nodes_tangential_curve(2,i+1))**2)
        if (dist_current < dist_tangential_detection_curve(i)) then
          dist_tangential_detection_curve(i) = dist_current
        endif
      enddo

      if (myrank == 0) then
        open(unit=11,file='OUTPUT_FILES/dist_rec_tangential_detection_curve', &
              form='formatted', status='unknown')
      endif
      irecloc = 0
      do irec = 1,nrec

        if (myrank == 0) then
          if (which_proc_receiver(irec) == myrank) then
            irecloc = irecloc + 1
            n1_tangential_detection_curve = rec_tangential_detection_curve(irecloc)
            x_final_receiver_dummy = x_final_receiver(irec)
            z_final_receiver_dummy = z_final_receiver(irec)
#ifdef USE_MPI
          else

            call MPI_RECV(n1_tangential_detection_curve,1,MPI_INTEGER,&
               which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
            call MPI_RECV(x_final_receiver_dummy,1,MPI_DOUBLE_PRECISION,&
               which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)
            call MPI_RECV(z_final_receiver_dummy,1,MPI_DOUBLE_PRECISION,&
               which_proc_receiver(irec),irec,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

#endif
          endif

#ifdef USE_MPI
        else
          if (which_proc_receiver(irec) == myrank) then
            irecloc = irecloc + 1
            call MPI_SEND(rec_tangential_detection_curve(irecloc),1,MPI_INTEGER,0,irec,MPI_COMM_WORLD,ier)
            call MPI_SEND(x_final_receiver(irec),1,MPI_DOUBLE_PRECISION,0,irec,MPI_COMM_WORLD,ier)
            call MPI_SEND(z_final_receiver(irec),1,MPI_DOUBLE_PRECISION,0,irec,MPI_COMM_WORLD,ier)
          endif
#endif

        endif
        if (myrank == 0) then
          write(11,*) dist_tangential_detection_curve(n1_tangential_detection_curve)
          write(12,*) x_final_receiver_dummy
          write(13,*) z_final_receiver_dummy
        endif
      enddo

      if (myrank == 0) then
        close(11)
        close(12)
        close(13)
      endif

    endif ! force_normal_to_surface

!
!---
!

! allocate seismogram arrays
  allocate(sisux(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrecloc))
  allocate(sisuz(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrecloc))
  allocate(siscurl(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrecloc))

! check if acoustic receiver is exactly on the free surface because pressure is zero there
  do ispec_acoustic_surface = 1,nelem_acoustic_surface
    ispec = acoustic_surface(1,ispec_acoustic_surface)
    ixmin = acoustic_surface(2,ispec_acoustic_surface)
    ixmax = acoustic_surface(3,ispec_acoustic_surface)
    izmin = acoustic_surface(4,ispec_acoustic_surface)
    izmax = acoustic_surface(5,ispec_acoustic_surface)
    do irecloc = 1,nrecloc
      irec = recloc(irecloc)
      if (acoustic(ispec) .and. ispec == ispec_selected_rec(irec)) then
        if ((izmin==1 .and. izmax==1 .and. ixmin==1 .and. ixmax==NGLLX .and. &
        gamma_receiver(irec) < -0.99d0) .or.&
        (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==NGLLX .and. &
        gamma_receiver(irec) > 0.99d0) .or.&
        (izmin==1 .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==1 .and. &
        xi_receiver(irec) < -0.99d0) .or.&
        (izmin==1 .and. izmax==NGLLZ .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
        xi_receiver(irec) > 0.99d0) .or.&
        (izmin==1 .and. izmax==1 .and. ixmin==1 .and. ixmax==1 .and. &
        gamma_receiver(irec) < -0.99d0 .and. xi_receiver(irec) < -0.99d0) .or.&
        (izmin==1 .and. izmax==1 .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
        gamma_receiver(irec) < -0.99d0 .and. xi_receiver(irec) > 0.99d0) .or.&
        (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==1 .and. &
        gamma_receiver(irec) > 0.99d0 .and. xi_receiver(irec) < -0.99d0) .or.&
        (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
        gamma_receiver(irec) > 0.99d0 .and. xi_receiver(irec) > 0.99d0)) then
          if (seismotype == 4) then
            call exit_MPI('an acoustic pressure receiver cannot be located exactly '// &
                            'on the free surface because pressure is zero there')
          else
            print *, '**********************************************************************'
            print *, '*** Warning: acoustic receiver located exactly on the free surface ***'
            print *, '*** Warning: tangential component will be zero there               ***'
            print *, '**********************************************************************'
            print *
          endif
        endif
      endif
    enddo
  enddo



  allocate(xir_store_loc(nrecloc,NGLLX))
  allocate(gammar_store_loc(nrecloc,NGLLX))

! define and store Lagrange interpolators at all the receivers
  irec_local=0
  do irec = 1,nrec

    if (AXISYM) then
      if (is_on_the_axis(ispec_selected_rec(irec)) .and. myrank == which_proc_receiver(irec)) then
        call lagrange_any(xi_receiver(irec),NGLJ,xiglj,hxir,hpxir)
        !do j = 1,NGLJ ! AB AB Same result with that loop
        !  hxir(j) = hglj(j-1,xi_receiver(irec),xiglj,NGLJ)
        !enddo
      else
        call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
      endif
    else
      call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
    endif

    call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
    hxir_store(irec,:) = hxir(:)
    hgammar_store(irec,:) = hgammar(:)

    if (myrank == which_proc_receiver(irec)) then
     irec_local = irec_local + 1
      do i = 1, NGLLX
        xir_store_loc(irec_local,i)    = sngl(hxir(i))
        gammar_store_loc(irec_local,i) = sngl(hgammar(i))
      enddo
    endif
  enddo

  ! define and store Lagrange interpolators at all the sources
  do i = 1,NSOURCES
    if (AXISYM) then
      if (is_on_the_axis(ispec_selected_source(i)) .and. is_proc_source(i) == 1) then
        call lagrange_any(xi_source(i),NGLJ,xiglj,hxis,hpxis)
        !do j = 1,NGLJ ! AB AB same result with that loop
        !  hxis(j) = hglj(j-1,xi_source(i),xiglj,NGLJ)
        !enddo
      else
        call lagrange_any(xi_source(i),NGLLX,xigll,hxis,hpxis)
      endif
    else
      call lagrange_any(xi_source(i),NGLLX,xigll,hxis,hpxis)
    endif

    call lagrange_any(gamma_source(i),NGLLZ,zigll,hgammas,hpgammas)
    hxis_store(i,:) = hxis(:)
    hgammas_store(i,:) = hgammas(:)
  enddo

! displacement, velocity, acceleration and inverse of the mass matrix for elastic elements
  if (myrank == 0) then
    write(IMAIN,*) "preparing array allocations..."
    call flush_IMAIN()
  endif

    if (any_elastic) then
      nglob_elastic = nglob
    else
      ! allocate unused arrays with fictitious size
      nglob_elastic = 1
    endif
    allocate(displ_elastic(3,nglob_elastic))
    allocate(displ_elastic_old(3,nglob_elastic))
    allocate(veloc_elastic(3,nglob_elastic))
    allocate(accel_elastic(3,nglob_elastic))
    allocate(accel_elastic_adj_coupling(3,nglob_elastic))
    allocate(accel_elastic_adj_coupling2(3,nglob_elastic))

    allocate(rmass_inverse_elastic_one(nglob_elastic))
    allocate(rmass_inverse_elastic_three(nglob_elastic))

    if (time_stepping_scheme==2) then
      allocate(displ_elastic_LDDRK(3,nglob_elastic))
      allocate(veloc_elastic_LDDRK(3,nglob_elastic))
      allocate(veloc_elastic_LDDRK_temp(3,nglob_elastic))
    endif

    if (time_stepping_scheme == 3) then
      allocate(accel_elastic_rk(3,nglob_elastic,stage_time_scheme))
      allocate(veloc_elastic_rk(3,nglob_elastic,stage_time_scheme))
      allocate(veloc_elastic_initial_rk(3,nglob_elastic))
      allocate(displ_elastic_initial_rk(3,nglob_elastic))
    endif

    ! extra array if adjoint and kernels calculation
    if (SIMULATION_TYPE == 3 .and. any_elastic) then
      allocate(b_displ_elastic(3,nglob))
      allocate(b_displ_elastic_old(3,nglob))
      allocate(b_veloc_elastic(3,nglob))
      allocate(b_accel_elastic(3,nglob))
      allocate(rho_kl(NGLLX,NGLLZ,nspec))
      allocate(rho_k(nglob))
      allocate(rhol_global(nglob))
      allocate(mu_kl(NGLLX,NGLLZ,nspec))
      allocate(mu_k(nglob))
      allocate(mul_global(nglob))
      allocate(kappa_kl(NGLLX,NGLLZ,nspec))
      allocate(kappa_k(nglob))
      allocate(kappal_global(nglob))
      allocate(rhop_kl(NGLLX,NGLLZ,nspec))
      allocate(alpha_kl(NGLLX,NGLLZ,nspec))
      allocate(beta_kl(NGLLX,NGLLZ,nspec))
      allocate(bulk_c_kl(NGLLX,NGLLZ,nspec))
      allocate(bulk_beta_kl(NGLLX,NGLLZ,nspec))
      if (APPROXIMATE_HESS_KL) then
        allocate(rhorho_el_hessian_final2(NGLLX,NGLLZ,nspec))
        allocate(rhorho_el_hessian_temp2(nglob))
        allocate(rhorho_el_hessian_final1(NGLLX,NGLLZ,nspec))
        allocate(rhorho_el_hessian_temp1(nglob))
      endif
    else
      allocate(b_displ_elastic(1,1))
      allocate(b_displ_elastic_old(1,1))
      allocate(b_veloc_elastic(1,1))
      allocate(b_accel_elastic(1,1))
      allocate(rho_kl(1,1,1))
      allocate(rho_k(1))
      allocate(rhol_global(1))
      allocate(mu_kl(1,1,1))
      allocate(mu_k(1))
      allocate(mul_global(1))
      allocate(kappa_kl(1,1,1))
      allocate(kappa_k(1))
      allocate(kappal_global(1))
      allocate(rhop_kl(1,1,1))
      allocate(alpha_kl(1,1,1))
      allocate(beta_kl(1,1,1))
      if (APPROXIMATE_HESS_KL) then
        allocate(rhorho_el_hessian_final2(1,1,1))
        allocate(rhorho_el_hessian_temp2(1))
        allocate(rhorho_el_hessian_final1(1,1,1))
        allocate(rhorho_el_hessian_temp1(1))
      endif
    endif

    if (any_poroelastic) then
      nglob_poroelastic = nglob
    else
      ! allocate unused arrays with fictitious size
      nglob_poroelastic = 1
    endif
    allocate(displs_poroelastic(NDIM,nglob_poroelastic))
    allocate(displs_poroelastic_old(NDIM,nglob_poroelastic))
    allocate(velocs_poroelastic(NDIM,nglob_poroelastic))
    allocate(accels_poroelastic(NDIM,nglob_poroelastic))
    allocate(accels_poroelastic_adj_coupling(NDIM,nglob_poroelastic))
    allocate(rmass_s_inverse_poroelastic(nglob_poroelastic))
    allocate(displw_poroelastic(NDIM,nglob_poroelastic))
    allocate(velocw_poroelastic(NDIM,nglob_poroelastic))
    allocate(accelw_poroelastic(NDIM,nglob_poroelastic))
    allocate(accelw_poroelastic_adj_coupling(NDIM,nglob_poroelastic))
    allocate(rmass_w_inverse_poroelastic(nglob_poroelastic))

    if (time_stepping_scheme == 2) then
      allocate(displs_poroelastic_LDDRK(NDIM,nglob_poroelastic))
      allocate(velocs_poroelastic_LDDRK(NDIM,nglob_poroelastic))
      allocate(displw_poroelastic_LDDRK(NDIM,nglob_poroelastic))
      allocate(velocw_poroelastic_LDDRK(NDIM,nglob_poroelastic))
    endif

    if (time_stepping_scheme == 3) then
      allocate(accels_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
      allocate(velocs_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
      allocate(accelw_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
      allocate(velocw_poroelastic_rk(NDIM,nglob_poroelastic,stage_time_scheme))
      allocate(displs_poroelastic_initial_rk(NDIM,nglob_poroelastic))
      allocate(velocs_poroelastic_initial_rk(NDIM,nglob_poroelastic))
      allocate(displw_poroelastic_initial_rk(NDIM,nglob_poroelastic))
      allocate(velocw_poroelastic_initial_rk(NDIM,nglob_poroelastic))
    endif

    ! extra array if adjoint and kernels calculation
    if (SIMULATION_TYPE == 3 .and. any_poroelastic) then
      allocate(b_displs_poroelastic(NDIM,nglob))
      allocate(b_velocs_poroelastic(NDIM,nglob))
      allocate(b_accels_poroelastic(NDIM,nglob))
      allocate(b_displw_poroelastic(NDIM,nglob))
      allocate(b_velocw_poroelastic(NDIM,nglob))
      allocate(b_accelw_poroelastic(NDIM,nglob))
      allocate(rhot_kl(NGLLX,NGLLZ,nspec))
      allocate(rhot_k(nglob))
      allocate(rhof_kl(NGLLX,NGLLZ,nspec))
      allocate(rhof_k(nglob))
      allocate(sm_kl(NGLLX,NGLLZ,nspec))
      allocate(sm_k(nglob))
      allocate(eta_kl(NGLLX,NGLLZ,nspec))
      allocate(eta_k(nglob))
      allocate(mufr_kl(NGLLX,NGLLZ,nspec))
      allocate(mufr_k(nglob))
      allocate(B_kl(NGLLX,NGLLZ,nspec))
      allocate(B_k(nglob))
      allocate(C_kl(NGLLX,NGLLZ,nspec))
      allocate(C_k(nglob))
      allocate(M_kl(NGLLX,NGLLZ,nspec))
      allocate(M_k(nglob))
      allocate(rhob_kl(NGLLX,NGLLZ,nspec))
      allocate(rhofb_kl(NGLLX,NGLLZ,nspec))
      allocate(phi_kl(NGLLX,NGLLZ,nspec))
      allocate(Bb_kl(NGLLX,NGLLZ,nspec))
      allocate(Cb_kl(NGLLX,NGLLZ,nspec))
      allocate(Mb_kl(NGLLX,NGLLZ,nspec))
      allocate(mufrb_kl(NGLLX,NGLLZ,nspec))
      allocate(rhobb_kl(NGLLX,NGLLZ,nspec))
      allocate(rhofbb_kl(NGLLX,NGLLZ,nspec))
      allocate(phib_kl(NGLLX,NGLLZ,nspec))
      allocate(cpI_kl(NGLLX,NGLLZ,nspec))
      allocate(cpII_kl(NGLLX,NGLLZ,nspec))
      allocate(cs_kl(NGLLX,NGLLZ,nspec))
      allocate(ratio_kl(NGLLX,NGLLZ,nspec))
      allocate(phil_global(nglob))
      allocate(mulfr_global(nglob))
      allocate(etal_f_global(nglob))
      allocate(rhol_s_global(nglob))
      allocate(rhol_f_global(nglob))
      allocate(rhol_bar_global(nglob))
      allocate(tortl_global(nglob))
      allocate(permlxx_global(nglob))
      allocate(permlxz_global(nglob))
      allocate(permlzz_global(nglob))
    else
      allocate(b_displs_poroelastic(1,1))
      allocate(b_velocs_poroelastic(1,1))
      allocate(b_accels_poroelastic(1,1))
      allocate(b_displw_poroelastic(1,1))
      allocate(b_velocw_poroelastic(1,1))
      allocate(b_accelw_poroelastic(1,1))
      allocate(rhot_kl(1,1,1))
      allocate(rhot_k(1))
      allocate(rhof_kl(1,1,1))
      allocate(rhof_k(1))
      allocate(sm_kl(1,1,1))
      allocate(sm_k(1))
      allocate(eta_kl(1,1,1))
      allocate(eta_k(1))
      allocate(mufr_kl(1,1,1))
      allocate(mufr_k(1))
      allocate(B_kl(1,1,1))
      allocate(B_k(1))
      allocate(C_kl(1,1,1))
      allocate(C_k(1))
      allocate(M_kl(1,1,1))
      allocate(M_k(1))
      allocate(rhob_kl(1,1,1))
      allocate(rhofb_kl(1,1,1))
      allocate(phi_kl(1,1,1))
      allocate(Bb_kl(1,1,1))
      allocate(Cb_kl(1,1,1))
      allocate(Mb_kl(1,1,1))
      allocate(mufrb_kl(1,1,1))
      allocate(rhobb_kl(1,1,1))
      allocate(rhofbb_kl(1,1,1))
      allocate(phib_kl(1,1,1))
      allocate(cpI_kl(1,1,1))
      allocate(cpII_kl(1,1,1))
      allocate(cs_kl(1,1,1))
      allocate(ratio_kl(1,1,1))
      allocate(phil_global(1))
      allocate(mulfr_global(1))
      allocate(etal_f_global(1))
      allocate(rhol_s_global(1))
      allocate(rhol_f_global(1))
      allocate(rhol_bar_global(1))
      allocate(tortl_global(1))
      allocate(permlxx_global(1))
      allocate(permlxz_global(1))
      allocate(permlzz_global(1))
    endif

    if (any_poroelastic .and. any_elastic) then
      allocate(icount(nglob))
    else
      allocate(icount(1))
    endif

    ! potential, its first and second derivative, and inverse of the mass matrix for acoustic elements
    if (any_acoustic) then
      nglob_acoustic = nglob
    else
      ! allocate unused arrays with fictitious size
      nglob_acoustic = 1
    endif
    allocate(potential_acoustic(nglob_acoustic))
    allocate(potential_acoustic_old(nglob_acoustic))
    allocate(potential_acoustic_adj_coupling(nglob_acoustic))
    allocate(potential_dot_acoustic(nglob_acoustic))
    allocate(potential_dot_dot_acoustic(nglob_acoustic))
    allocate(rmass_inverse_acoustic(nglob_acoustic))
    if (time_stepping_scheme == 2) then
      allocate(potential_acoustic_LDDRK(nglob_acoustic))
      allocate(potential_dot_acoustic_LDDRK(nglob_acoustic))
      allocate(potential_dot_acoustic_temp(nglob_acoustic))
    endif

    if (time_stepping_scheme == 3) then
      allocate(potential_acoustic_init_rk(nglob_acoustic))
      allocate(potential_dot_acoustic_init_rk(nglob_acoustic))
      allocate(potential_dot_dot_acoustic_rk(nglob_acoustic,stage_time_scheme))
      allocate(potential_dot_acoustic_rk(nglob_acoustic,stage_time_scheme))
    endif

    if (SIMULATION_TYPE == 3 .and. any_acoustic) then
      allocate(b_potential_acoustic(nglob))
      allocate(b_potential_acoustic_old(nglob))
      allocate(b_potential_dot_acoustic(nglob))
      allocate(b_potential_dot_dot_acoustic(nglob))
      allocate(b_displ_ac(2,nglob))
      allocate(b_accel_ac(2,nglob))
      allocate(accel_ac(2,nglob))
      allocate(rho_ac_kl(NGLLX,NGLLZ,nspec))
      allocate(rhol_ac_global(nglob))
      allocate(kappa_ac_kl(NGLLX,NGLLZ,nspec))
      allocate(kappal_ac_global(nglob))
      allocate(rhop_ac_kl(NGLLX,NGLLZ,nspec))
      allocate(alpha_ac_kl(NGLLX,NGLLZ,nspec))
      if (APPROXIMATE_HESS_KL) then
        allocate(rhorho_ac_hessian_final2(NGLLX,NGLLZ,nspec))
        allocate(rhorho_ac_hessian_final1(NGLLX,NGLLZ,nspec))
      endif
    else
    ! allocate unused arrays with fictitious size
      allocate(b_potential_acoustic(1))
      allocate(b_potential_acoustic_old(1))
      allocate(b_potential_dot_acoustic(1))
      allocate(b_potential_dot_dot_acoustic(1))
      allocate(b_displ_ac(1,1))
      allocate(b_accel_ac(1,1))
      allocate(accel_ac(1,1))
      allocate(rho_ac_kl(1,1,1))
      allocate(rhol_ac_global(1))
      allocate(kappa_ac_kl(1,1,1))
      allocate(kappal_ac_global(1))
      allocate(rhop_ac_kl(1,1,1))
      allocate(alpha_ac_kl(1,1,1))
      if (APPROXIMATE_HESS_KL) then
        allocate(rhorho_ac_hessian_final2(1,1,1))
        allocate(rhorho_ac_hessian_final1(1,1,1))
      endif
    endif

    ! potential, its first and second derivative, and inverse of the mass matrix for gravitoacoustic elements
    if (any_gravitoacoustic) then
      nglob_gravitoacoustic = nglob
    else
      ! allocate unused arrays with fictitious size
      nglob_gravitoacoustic = 1
    endif
    allocate(potential_gravitoacoustic(nglob_gravitoacoustic))
    allocate(potential_dot_gravitoacoustic(nglob_gravitoacoustic))
    allocate(potential_dot_dot_gravitoacoustic(nglob_gravitoacoustic))
    allocate(rmass_inverse_gravitoacoustic(nglob_gravitoacoustic))
    allocate(potential_gravito(nglob_gravitoacoustic))
    allocate(potential_dot_gravito(nglob_gravitoacoustic))
    allocate(potential_dot_dot_gravito(nglob_gravitoacoustic))
    allocate(rmass_inverse_gravito(nglob_gravitoacoustic))



  if (myrank == 0) then
    write(IMAIN,*) "preparing PML..."
    call flush_IMAIN()
  endif
  call prepare_timerun_pml()


! Test compatibility with axisymmetric formulation
  if (AXISYM) call check_compatibility_axisym()

  if (myrank == 0) then
    write(IMAIN,*) "preparing mass matrices..."
    call flush_IMAIN()
  endif
  call prepare_timerun_mass_matrix()

  if (myrank == 0) then
    write(IMAIN,*) "preparing image coloring..."
    call flush_IMAIN()
  endif
  call prepare_timerun_image_coloring()
!
!---- initialize seismograms
!
  if (myrank == 0) then
    write(IMAIN,*) "preparing array initializations..."
    call flush_IMAIN()
  endif

  sisux = ZERO ! double precision zero
  sisuz = ZERO

! initialize arrays to zero
  displ_elastic = 0._CUSTOM_REAL
  displ_elastic_old = 0._CUSTOM_REAL
  veloc_elastic = 0._CUSTOM_REAL
  accel_elastic = 0._CUSTOM_REAL

    if (SIMULATION_TYPE == 3 .and. any_elastic) then
      b_displ_elastic_old = 0._CUSTOM_REAL
      b_displ_elastic = 0._CUSTOM_REAL
      b_veloc_elastic = 0._CUSTOM_REAL
      b_accel_elastic = 0._CUSTOM_REAL
    endif

  if (time_stepping_scheme == 2) then
  displ_elastic_LDDRK = 0._CUSTOM_REAL
  veloc_elastic_LDDRK = 0._CUSTOM_REAL
  veloc_elastic_LDDRK_temp = 0._CUSTOM_REAL
  endif

  if (time_stepping_scheme == 3) then
   accel_elastic_rk = 0._CUSTOM_REAL
   veloc_elastic_rk = 0._CUSTOM_REAL
   veloc_elastic_initial_rk = 0._CUSTOM_REAL
   displ_elastic_initial_rk = 0._CUSTOM_REAL
  endif

  displs_poroelastic = 0._CUSTOM_REAL
  displs_poroelastic_old = 0._CUSTOM_REAL
  velocs_poroelastic = 0._CUSTOM_REAL
  accels_poroelastic = 0._CUSTOM_REAL
  displw_poroelastic = 0._CUSTOM_REAL
  velocw_poroelastic = 0._CUSTOM_REAL
  accelw_poroelastic = 0._CUSTOM_REAL

  if (time_stepping_scheme == 2) then
    displs_poroelastic_LDDRK = 0._CUSTOM_REAL
    velocs_poroelastic_LDDRK = 0._CUSTOM_REAL
    displw_poroelastic_LDDRK = 0._CUSTOM_REAL
    velocw_poroelastic_LDDRK = 0._CUSTOM_REAL
  endif

  if (time_stepping_scheme == 3) then
    accels_poroelastic_rk = 0._CUSTOM_REAL
    velocs_poroelastic_rk = 0._CUSTOM_REAL

    accelw_poroelastic_rk = 0._CUSTOM_REAL
    velocw_poroelastic_rk = 0._CUSTOM_REAL

    velocs_poroelastic_initial_rk = 0._CUSTOM_REAL
    displs_poroelastic_initial_rk = 0._CUSTOM_REAL

    velocw_poroelastic_initial_rk = 0._CUSTOM_REAL
    displw_poroelastic_initial_rk = 0._CUSTOM_REAL

  endif

  potential_acoustic = 0._CUSTOM_REAL
  potential_acoustic_old = 0._CUSTOM_REAL
  potential_dot_acoustic = 0._CUSTOM_REAL
  potential_dot_dot_acoustic = 0._CUSTOM_REAL

  if (time_stepping_scheme == 2) then
    potential_acoustic_LDDRK = 0._CUSTOM_REAL
    potential_dot_acoustic_LDDRK = 0._CUSTOM_REAL
    potential_dot_acoustic_temp = 0._CUSTOM_REAL
  endif

  if (time_stepping_scheme == 3) then
    potential_acoustic_init_rk = 0._CUSTOM_REAL
    potential_dot_acoustic_init_rk = 0._CUSTOM_REAL
    potential_dot_dot_acoustic_rk = 0._CUSTOM_REAL
    potential_dot_acoustic_rk = 0._CUSTOM_REAL
  endif

  potential_gravitoacoustic = 0._CUSTOM_REAL
  potential_dot_gravitoacoustic = 0._CUSTOM_REAL
  potential_dot_dot_gravitoacoustic = 0._CUSTOM_REAL
  potential_gravito = 0._CUSTOM_REAL
  potential_dot_gravito = 0._CUSTOM_REAL
  potential_dot_dot_gravito = 0._CUSTOM_REAL

!
!----- Files where viscous damping are saved during forward wavefield calculation
!
  if (any_poroelastic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3)) then
    allocate(b_viscodampx(nglob))
    allocate(b_viscodampz(nglob))
    write(outputname,'(a,i6.6,a)') 'viscodampingx',myrank,'.bin'
    write(outputname2,'(a,i6.6,a)') 'viscodampingz',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      reclen = CUSTOM_REAL * nglob
      open(unit=23,file='OUTPUT_FILES/'//outputname,status='old',&
            action='read',form='unformatted',access='direct',&
            recl=reclen)
      open(unit=24,file='OUTPUT_FILES/'//outputname2,status='old',&
            action='read',form='unformatted',access='direct',&
            recl=reclen)
    else
      reclen = CUSTOM_REAL * nglob
      open(unit=23,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',access='direct',&
            recl=reclen)
      open(unit=24,file='OUTPUT_FILES/'//outputname2,status='unknown',&
            form='unformatted',access='direct',&
            recl=reclen)
    endif
  else
    allocate(b_viscodampx(1))
    allocate(b_viscodampz(1))
  endif

!
!----- Read last frame for forward wavefield reconstruction
!

  if (((SAVE_FORWARD .and. SIMULATION_TYPE ==1) .or. SIMULATION_TYPE == 3) .and. anyabs &
      .and. (.not. PML_BOUNDARY_CONDITIONS)) then
    ! opens files for absorbing boundary data
    call prepare_absorb_files()
  endif

  if (anyabs .and. SIMULATION_TYPE == 3 .and. (.not. PML_BOUNDARY_CONDITIONS)) then

    ! reads in absorbing boundary data

    if (any_elastic) call prepare_absorb_elastic()

    if (any_poroelastic) call prepare_absorb_poroelastic()

    if (any_acoustic) call prepare_absorb_acoustic()

  endif ! if (anyabs .and. SIMULATION_TYPE == 3)

  if (myrank == 0) then
    write(IMAIN,*) "preparing kernels..."
    call flush_IMAIN()
  endif
  call prepare_timerun_kernel()

!
!----  read initial fields from external file if needed
!

! if we are looking a plane wave beyond critical angle we use other method
  over_critical_angle = .false.

  if (initialfield) then

    ! Calculation of the initial field for a plane wave
    if (any_elastic) then
      call prepare_initialfield()
    endif

    if (over_critical_angle) then

      allocate(left_bound(nelemabs*NGLLX))
      allocate(right_bound(nelemabs*NGLLX))
      allocate(bot_bound(nelemabs*NGLLZ))

      call prepare_initialfield_paco()

      allocate(v0x_left(count_left,NSTEP))
      allocate(v0z_left(count_left,NSTEP))
      allocate(t0x_left(count_left,NSTEP))
      allocate(t0z_left(count_left,NSTEP))

      allocate(v0x_right(count_right,NSTEP))
      allocate(v0z_right(count_right,NSTEP))
      allocate(t0x_right(count_right,NSTEP))
      allocate(t0z_right(count_right,NSTEP))

      allocate(v0x_bot(count_bottom,NSTEP))
      allocate(v0z_bot(count_bottom,NSTEP))
      allocate(t0x_bot(count_bottom,NSTEP))
      allocate(t0z_bot(count_bottom,NSTEP))

      ! call Paco's routine to compute in frequency and convert to time by Fourier transform
      call paco_beyond_critical(anglesource(1),&
                                f0(1),QKappa_attenuation(1),source_type(1),left_bound(1:count_left),&
                                right_bound(1:count_right),bot_bound(1:count_bottom), &
                                count_left,count_right,count_bottom,x_source(1))

      deallocate(left_bound)
      deallocate(right_bound)
      deallocate(bot_bound)

      if (myrank == 0) then
        write(IMAIN,*)  '***********'
        write(IMAIN,*)  'done calculating the initial wave field'
        write(IMAIN,*)  '***********'
      endif

    endif ! beyond critical angle

    write(IMAIN,*) 'Max norm of initial elastic displacement = ', &
      maxval(sqrt(displ_elastic(1,:)**2 + displ_elastic(3,:)**2))

  endif ! initialfield

! compute the source time function and store it in a text file
  if (.not. initialfield) then

    allocate(source_time_function(NSOURCES,NSTEP,stage_time_scheme))
    source_time_function(:,:,:) = 0._CUSTOM_REAL

    ! computes source time function array
    call prepare_source_time_function()
  else
    ! uses an initialfield
    ! dummy allocation
    allocate(source_time_function(1,1,1))
  endif

! acoustic forcing edge detection
! the elements forming an edge are already known (computed in meshfem2D),
! the common nodes forming the edge are computed here
  if (ACOUSTIC_FORCING) then

    if (myrank == 0) then
      print *
      print *,'Acoustic forcing simulation'
      print *
      print *,'Beginning of acoustic forcing edge detection'
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


! determine if coupled fluid-solid simulation
  coupled_acoustic_elastic = any_acoustic .and. any_elastic
  coupled_acoustic_poro = any_acoustic .and. any_poroelastic

! fluid/solid (elastic) edge detection
! the two elements (fluid and solid) forming an edge are already known (computed in meshfem2D),
! the common nodes forming the edge are computed here
  if (coupled_acoustic_elastic) then

    if (myrank == 0) then
      print *
      print *,'Mixed acoustic/elastic simulation'
      print *
      print *,'Beginning of fluid/solid edge detection'
    endif

! define the edges of a given element
    i_begin(IBOTTOM) = 1
    j_begin(IBOTTOM) = 1
    i_end(IBOTTOM) = NGLLX
    j_end(IBOTTOM) = 1

    i_begin(IRIGHT) = NGLLX
    j_begin(IRIGHT) = 1
    i_end(IRIGHT) = NGLLX
    j_end(IRIGHT) = NGLLZ

    i_begin(ITOP) = NGLLX
    j_begin(ITOP) = NGLLZ
    i_end(ITOP) = 1
    j_end(ITOP) = NGLLZ

    i_begin(ILEFT) = 1
    j_begin(ILEFT) = NGLLZ
    i_end(ILEFT) = 1
    j_end(ILEFT) = 1

! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = ipoin1D
      ivalue_inverse(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,IBOTTOM) = 1
      jvalue_inverse(ipoin1D,IBOTTOM) = 1

      ivalue(ipoin1D,IRIGHT) = NGLLX
      ivalue_inverse(ipoin1D,IRIGHT) = NGLLX
      jvalue(ipoin1D,IRIGHT) = ipoin1D
      jvalue_inverse(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1

      ivalue(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,ITOP) = ipoin1D
      jvalue(ipoin1D,ITOP) = NGLLZ
      jvalue_inverse(ipoin1D,ITOP) = NGLLZ

      ivalue(ipoin1D,ILEFT) = 1
      ivalue_inverse(ipoin1D,ILEFT) = 1
      jvalue(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,ILEFT) = ipoin1D

    enddo

    do inum = 1, num_fluid_solid_edges
       ispec_acoustic =  fluid_solid_acoustic_ispec(inum)
       ispec_elastic =  fluid_solid_elastic_ispec(inum)

! one element must be acoustic and the other must be elastic
        if (ispec_acoustic /= ispec_elastic .and. .not. elastic(ispec_acoustic) .and. &
             .not. poroelastic(ispec_acoustic) .and. elastic(ispec_elastic)) then

! loop on the four edges of the two elements
          do iedge_acoustic = 1,NEDGES
            do iedge_elastic = 1,NEDGES
! store the matching topology if the two edges match in inverse order
              if (ibool(i_begin(iedge_acoustic),j_begin(iedge_acoustic),ispec_acoustic) == &
                  ibool(i_end(iedge_elastic),j_end(iedge_elastic),ispec_elastic) .and. &
                  ibool(i_end(iedge_acoustic),j_end(iedge_acoustic),ispec_acoustic) == &
                  ibool(i_begin(iedge_elastic),j_begin(iedge_elastic),ispec_elastic)) then
                 fluid_solid_acoustic_iedge(inum) = iedge_acoustic
                 fluid_solid_elastic_iedge(inum) = iedge_elastic
!                  print *,'edge ',iedge_acoustic,' of acoustic element ',ispec_acoustic, &
!                          ' is in contact with edge ',iedge_elastic,' of elastic element ',ispec_elastic
              endif

            enddo
          enddo

       endif

    enddo

! make sure fluid/solid (elastic) matching has been perfectly detected: check that the grid points
! have the same physical coordinates
! loop on all the coupling edges

    if (myrank == 0) print *,'Checking fluid/solid edge topology...'

    do inum = 1,num_fluid_solid_edges

! get the edge of the acoustic element
      ispec_acoustic = fluid_solid_acoustic_ispec(inum)
      iedge_acoustic = fluid_solid_acoustic_iedge(inum)

! get the corresponding edge of the elastic element
      ispec_elastic = fluid_solid_elastic_ispec(inum)
      iedge_elastic = fluid_solid_elastic_iedge(inum)

! implement 1D coupling along the edge
      do ipoin1D = 1,NGLLX

! get point values for the elastic side, which matches our side in the inverse direction
        i = ivalue_inverse(ipoin1D,iedge_elastic)
        j = jvalue_inverse(ipoin1D,iedge_elastic)
        iglob = ibool(i,j,ispec_elastic)

! get point values for the acoustic side
        i = ivalue(ipoin1D,iedge_acoustic)
        j = jvalue(ipoin1D,iedge_acoustic)
        iglob2 = ibool(i,j,ispec_acoustic)

! if distance between the two points is not negligible, there is an error, since it should be zero
        if (sqrt((coord(1,iglob) - coord(1,iglob2))**2 + (coord(2,iglob) - coord(2,iglob2))**2) > TINYVAL) &
            call exit_MPI( 'error in fluid/solid coupling buffer')

      enddo

    enddo

    if (myrank == 0) then
      print *,'End of fluid/solid edge detection'
      print *
    endif

  endif

! fluid/solid (poroelastic) edge detection
! the two elements (fluid and solid) forming an edge are already known (computed in meshfem2D),
! the common nodes forming the edge are computed here
  if (coupled_acoustic_poro) then
    if (myrank == 0) then
    print *
    print *,'Mixed acoustic/poroelastic simulation'
    print *
    print *,'Beginning of fluid/solid (poroelastic) edge detection'
    endif

! define the edges of a given element
    i_begin(IBOTTOM) = 1
    j_begin(IBOTTOM) = 1
    i_end(IBOTTOM) = NGLLX
    j_end(IBOTTOM) = 1

    i_begin(IRIGHT) = NGLLX
    j_begin(IRIGHT) = 1
    i_end(IRIGHT) = NGLLX
    j_end(IRIGHT) = NGLLZ

    i_begin(ITOP) = NGLLX
    j_begin(ITOP) = NGLLZ
    i_end(ITOP) = 1
    j_end(ITOP) = NGLLZ

    i_begin(ILEFT) = 1
    j_begin(ILEFT) = NGLLZ
    i_end(ILEFT) = 1
    j_end(ILEFT) = 1

! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = ipoin1D
      ivalue_inverse(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,IBOTTOM) = 1
      jvalue_inverse(ipoin1D,IBOTTOM) = 1

      ivalue(ipoin1D,IRIGHT) = NGLLX
      ivalue_inverse(ipoin1D,IRIGHT) = NGLLX
      jvalue(ipoin1D,IRIGHT) = ipoin1D
      jvalue_inverse(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1

      ivalue(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,ITOP) = ipoin1D
      jvalue(ipoin1D,ITOP) = NGLLZ
      jvalue_inverse(ipoin1D,ITOP) = NGLLZ

      ivalue(ipoin1D,ILEFT) = 1
      ivalue_inverse(ipoin1D,ILEFT) = 1
      jvalue(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,ILEFT) = ipoin1D

    enddo

    do inum = 1, num_fluid_poro_edges
       ispec_acoustic =  fluid_poro_acoustic_ispec(inum)
       ispec_poroelastic =  fluid_poro_poroelastic_ispec(inum)

! one element must be acoustic and the other must be poroelastic
        if (ispec_acoustic /= ispec_poroelastic .and. .not. poroelastic(ispec_acoustic) .and. &
                 .not. elastic(ispec_acoustic) .and. poroelastic(ispec_poroelastic)) then

! loop on the four edges of the two elements
          do iedge_acoustic = 1,NEDGES
            do iedge_poroelastic = 1,NEDGES

! store the matching topology if the two edges match in inverse order
              if (ibool(i_begin(iedge_acoustic),j_begin(iedge_acoustic),ispec_acoustic) == &
                   ibool(i_end(iedge_poroelastic),j_end(iedge_poroelastic),ispec_poroelastic) .and. &
                   ibool(i_end(iedge_acoustic),j_end(iedge_acoustic),ispec_acoustic) == &
                   ibool(i_begin(iedge_poroelastic),j_begin(iedge_poroelastic),ispec_poroelastic)) then
                 fluid_poro_acoustic_iedge(inum) = iedge_acoustic
                 fluid_poro_poroelastic_iedge(inum) = iedge_poroelastic
                endif

             enddo
          enddo

       endif

    enddo


! make sure fluid/solid (poroelastic) matching has been perfectly detected: check that the grid points
! have the same physical coordinates
! loop on all the coupling edges

    if (myrank == 0) then
    print *,'Checking fluid/solid (poroelastic) edge topology...'
    endif

    do inum = 1,num_fluid_poro_edges

! get the edge of the acoustic element
      ispec_acoustic = fluid_poro_acoustic_ispec(inum)
      iedge_acoustic = fluid_poro_acoustic_iedge(inum)

! get the corresponding edge of the poroelastic element
      ispec_poroelastic = fluid_poro_poroelastic_ispec(inum)
      iedge_poroelastic = fluid_poro_poroelastic_iedge(inum)

! implement 1D coupling along the edge
      do ipoin1D = 1,NGLLX

! get point values for the poroelastic side, which matches our side in the inverse direction
        i = ivalue_inverse(ipoin1D,iedge_poroelastic)
        j = jvalue_inverse(ipoin1D,iedge_poroelastic)
        iglob = ibool(i,j,ispec_poroelastic)

! get point values for the acoustic side
        i = ivalue(ipoin1D,iedge_acoustic)
        j = jvalue(ipoin1D,iedge_acoustic)
        iglob2 = ibool(i,j,ispec_acoustic)

! if distance between the two points is not negligible, there is an error, since it should be zero
        if (sqrt((coord(1,iglob) - coord(1,iglob2))**2 + (coord(2,iglob) - coord(2,iglob2))**2) > TINYVAL) &
            call exit_MPI( 'error in fluid/solid (poroelastic) coupling buffer')

      enddo

    enddo

    if (myrank == 0) then
    print *,'End of fluid/solid (poroelastic) edge detection'
    print *
    endif

  endif

! exclude common points between acoustic absorbing edges and acoustic/elastic matching interfaces
  if (coupled_acoustic_elastic .and. anyabs) then

    if (myrank == 0) &
      print *,'excluding common points between acoustic absorbing edges and acoustic/elastic matching interfaces, if any'

! loop on all the absorbing elements
    do ispecabs = 1,nelemabs

      ispec = numabs(ispecabs)

! loop on all the coupling edges
      do inum = 1,num_fluid_solid_edges

! get the edge of the acoustic element
        ispec_acoustic = fluid_solid_acoustic_ispec(inum)
        iedge_acoustic = fluid_solid_acoustic_iedge(inum)

! if acoustic absorbing element and acoustic/elastic coupled element is the same
        if (ispec_acoustic == ispec) then

          if (iedge_acoustic == IBOTTOM) then
            ibegin_edge4(ispecabs) = 2
            ibegin_edge2(ispecabs) = 2
          endif

          if (iedge_acoustic == ITOP) then
            iend_edge4(ispecabs) = NGLLZ - 1
            iend_edge2(ispecabs) = NGLLZ - 1
          endif

          if (iedge_acoustic == ILEFT) then
            ibegin_edge1(ispecabs) = 2
            ibegin_edge3(ispecabs) = 2
          endif

          if (iedge_acoustic == IRIGHT) then
            iend_edge1(ispecabs) = NGLLX - 1
            iend_edge3(ispecabs) = NGLLX - 1
          endif

        endif

      enddo

    enddo

  endif

! exclude common points between acoustic absorbing edges and acoustic/poroelastic matching interfaces
  if (coupled_acoustic_poro .and. anyabs) then

    if (myrank == 0) &
      print *,'excluding common points between acoustic absorbing edges and acoustic/poroelastic matching interfaces, if any'

! loop on all the absorbing elements
    do ispecabs = 1,nelemabs

      ispec = numabs(ispecabs)

! loop on all the coupling edges
      do inum = 1,num_fluid_poro_edges

! get the edge of the acoustic element
        ispec_acoustic = fluid_poro_acoustic_ispec(inum)
        iedge_acoustic = fluid_poro_acoustic_iedge(inum)

! if acoustic absorbing element and acoustic/poroelastic coupled element is the same
        if (ispec_acoustic == ispec) then

          if (iedge_acoustic == IBOTTOM) then
            ibegin_edge4(ispecabs) = 2
            ibegin_edge2(ispecabs) = 2
          endif

          if (iedge_acoustic == ITOP) then
            iend_edge4(ispecabs) = NGLLZ - 1
            iend_edge2(ispecabs) = NGLLZ - 1
          endif

          if (iedge_acoustic == ILEFT) then
            ibegin_edge1(ispecabs) = 2
            ibegin_edge3(ispecabs) = 2
          endif

          if (iedge_acoustic == IRIGHT) then
            iend_edge1(ispecabs) = NGLLX - 1
            iend_edge3(ispecabs) = NGLLX - 1
          endif

        endif

      enddo

    enddo

  endif


! determine if coupled elastic-poroelastic simulation
  coupled_elastic_poro = any_elastic .and. any_poroelastic

! solid/porous edge detection
! the two elements forming an edge are already known (computed in meshfem2D),
! the common nodes forming the edge are computed here

if (ATTENUATION_PORO_FLUID_PART .and. any_poroelastic .and. &
(time_stepping_scheme == 2.or. time_stepping_scheme == 3)) &
    stop 'RK and LDDRK time scheme not supported poroelastic simulations with attenuation'

if (coupled_elastic_poro) then

    if (ATTENUATION_VISCOELASTIC_SOLID .or. ATTENUATION_PORO_FLUID_PART) &
                   stop 'Attenuation not supported for mixed elastic/poroelastic simulations'

    if (time_stepping_scheme == 2.or. time_stepping_scheme == 3) &
                   stop 'RK and LDDRK time scheme not supported for mixed elastic/poroelastic simulations'

    if (myrank == 0) then
    print *
    print *,'Mixed elastic/poroelastic simulation'
    print *
    print *,'Beginning of solid/porous edge detection'
    endif

! define the edges of a given element
    i_begin(IBOTTOM) = 1
    j_begin(IBOTTOM) = 1
    i_end(IBOTTOM) = NGLLX
    j_end(IBOTTOM) = 1

    i_begin(IRIGHT) = NGLLX
    j_begin(IRIGHT) = 1
    i_end(IRIGHT) = NGLLX
    j_end(IRIGHT) = NGLLZ

    i_begin(ITOP) = NGLLX
    j_begin(ITOP) = NGLLZ
    i_end(ITOP) = 1
    j_end(ITOP) = NGLLZ

    i_begin(ILEFT) = 1
    j_begin(ILEFT) = NGLLZ
    i_end(ILEFT) = 1
    j_end(ILEFT) = 1

! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = ipoin1D
      ivalue_inverse(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,IBOTTOM) = 1
      jvalue_inverse(ipoin1D,IBOTTOM) = 1

      ivalue(ipoin1D,IRIGHT) = NGLLX
      ivalue_inverse(ipoin1D,IRIGHT) = NGLLX
      jvalue(ipoin1D,IRIGHT) = ipoin1D
      jvalue_inverse(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1

      ivalue(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,ITOP) = ipoin1D
      jvalue(ipoin1D,ITOP) = NGLLZ
      jvalue_inverse(ipoin1D,ITOP) = NGLLZ

      ivalue(ipoin1D,ILEFT) = 1
      ivalue_inverse(ipoin1D,ILEFT) = 1
      jvalue(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,ILEFT) = ipoin1D

    enddo


    do inum = 1, num_solid_poro_edges
       ispec_elastic =  solid_poro_elastic_ispec(inum)
       ispec_poroelastic =  solid_poro_poroelastic_ispec(inum)

! one element must be elastic and the other must be poroelastic
        if (ispec_elastic /= ispec_poroelastic .and. elastic(ispec_elastic) .and. &
                 poroelastic(ispec_poroelastic)) then

! loop on the four edges of the two elements
          do iedge_poroelastic = 1,NEDGES
            do iedge_elastic = 1,NEDGES

! store the matching topology if the two edges match in inverse order
              if (ibool(i_begin(iedge_poroelastic),j_begin(iedge_poroelastic),ispec_poroelastic) == &
                   ibool(i_end(iedge_elastic),j_end(iedge_elastic),ispec_elastic) .and. &
                   ibool(i_end(iedge_poroelastic),j_end(iedge_poroelastic),ispec_poroelastic) == &
                   ibool(i_begin(iedge_elastic),j_begin(iedge_elastic),ispec_elastic)) then
                 solid_poro_elastic_iedge(inum) = iedge_elastic
                 solid_poro_poroelastic_iedge(inum) = iedge_poroelastic
                endif

             enddo
          enddo

       endif

    enddo

! make sure solid/porous matching has been perfectly detected: check that the grid points
! have the same physical coordinates
! loop on all the coupling edges

    if (myrank == 0) then
    print *,'Checking solid/porous edge topology...'
    endif

    do inum = 1,num_solid_poro_edges

! get the edge of the elastic element
      ispec_elastic = solid_poro_elastic_ispec(inum)
      iedge_elastic = solid_poro_elastic_iedge(inum)

! get the corresponding edge of the poroelastic element
      ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
      iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

! implement 1D coupling along the edge
      do ipoin1D = 1,NGLLX

! get point values for the poroelastic side, which matches our side in the inverse direction
        i = ivalue_inverse(ipoin1D,iedge_elastic)
        j = jvalue_inverse(ipoin1D,iedge_elastic)
        iglob = ibool(i,j,ispec_elastic)

! get point values for the elastic side
        i = ivalue(ipoin1D,iedge_poroelastic)
        j = jvalue(ipoin1D,iedge_poroelastic)
        iglob2 = ibool(i,j,ispec_poroelastic)

! if distance between the two points is not negligible, there is an error, since it should be zero
        if (sqrt((coord(1,iglob) - coord(1,iglob2))**2 + (coord(2,iglob) - coord(2,iglob2))**2) > TINYVAL) &
            call exit_MPI( 'error in solid/porous coupling buffer')

      enddo

    enddo

    if (myrank == 0) then
    print *,'End of solid/porous edge detection'
    print *
    endif

  endif

! initiation
 if (any_poroelastic .and. anyabs) then
! loop on all the absorbing elements
    do ispecabs = 1,nelemabs
            ibegin_edge4_poro(ispecabs) = 1
            ibegin_edge2_poro(ispecabs) = 1

            iend_edge4_poro(ispecabs) = NGLLZ
            iend_edge2_poro(ispecabs) = NGLLZ

            ibegin_edge1_poro(ispecabs) = 1
            ibegin_edge3_poro(ispecabs) = 1

            iend_edge1_poro(ispecabs) = NGLLX
            iend_edge3_poro(ispecabs) = NGLLX
    enddo
 endif

! exclude common points between poroelastic absorbing edges and elastic/poroelastic matching interfaces
  if (coupled_elastic_poro .and. anyabs) then

    if (myrank == 0) &
      print *,'excluding common points between poroelastic absorbing edges and elastic/poroelastic matching interfaces, if any'

! loop on all the absorbing elements
    do ispecabs = 1,nelemabs

      ispec = numabs(ispecabs)

! loop on all the coupling edges
      do inum = 1,num_solid_poro_edges

! get the edge of the acoustic element
        ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
        iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

! if poroelastic absorbing element and elastic/poroelastic coupled element is the same
        if (ispec_poroelastic == ispec) then

          if (iedge_poroelastic == IBOTTOM) then
            ibegin_edge4_poro(ispecabs) = 2
            ibegin_edge2_poro(ispecabs) = 2
          endif

          if (iedge_poroelastic == ITOP) then
            iend_edge4_poro(ispecabs) = NGLLZ - 1
            iend_edge2_poro(ispecabs) = NGLLZ - 1
          endif

          if (iedge_poroelastic == ILEFT) then
            ibegin_edge1_poro(ispecabs) = 2
            ibegin_edge3_poro(ispecabs) = 2
          endif

          if (iedge_poroelastic == IRIGHT) then
            iend_edge1_poro(ispecabs) = NGLLX - 1
            iend_edge3_poro(ispecabs) = NGLLX - 1
          endif

        endif

      enddo

    enddo

  endif

!----  create a Gnuplot script to display the energy curve in log scale
  if (output_energy .and. myrank == 0) then
    close(IOUT_ENERGY)
    open(unit=IOUT_ENERGY,file='plot_energy.gnu',status='unknown',action='write')
    write(IOUT_ENERGY,*) 'set term wxt'
    write(IOUT_ENERGY,*) '#set term postscript landscape color solid "Helvetica" 22'
    write(IOUT_ENERGY,*) '#set output "energy.ps"'
    write(IOUT_ENERGY,*) '# set xrange [0:60]'
    write(IOUT_ENERGY,*) 'set logscale y'
    write(IOUT_ENERGY,*) 'set xlabel "Time (s)"'
    write(IOUT_ENERGY,*) 'set ylabel "Energy (J)"'
    write(IOUT_ENERGY,'(A)') &
      'plot "energy.dat" us 1:4 t ''Total Energy'' w l lc 1, "energy.dat" us 1:3 t ''Potential Energy'' w l lc 2'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    close(IOUT_ENERGY)
  endif

! open the file in which we will store the energy curve
  if (output_energy .and. myrank == 0) open(unit=IOUT_ENERGY,file='energy.dat',status='unknown',action='write')

  if (myrank == 0) then
    write(IMAIN,*) "preparing noise..."
    call flush_IMAIN()
  endif
  call prepare_timerun_noise()


  ! prepares image background
  if (output_color_image) then
    if (myrank == 0) then
      write(IMAIN,*) "preparing color image vp..."
      call flush_IMAIN()
    endif
    call prepare_color_image_vp()
  endif

! dummy allocation of plane wave arrays if they are unused (but still need to exist because
! they are used as arguments to subroutines)
  if (.not. over_critical_angle) then
    allocate(v0x_left(1,NSTEP))
    allocate(v0z_left(1,NSTEP))
    allocate(t0x_left(1,NSTEP))
    allocate(t0z_left(1,NSTEP))

    allocate(v0x_right(1,NSTEP))
    allocate(v0z_right(1,NSTEP))
    allocate(t0x_right(1,NSTEP))
    allocate(t0z_right(1,NSTEP))

    allocate(v0x_bot(1,NSTEP))
    allocate(v0z_bot(1,NSTEP))
    allocate(t0x_bot(1,NSTEP))
    allocate(t0z_bot(1,NSTEP))
  endif

! initialize variables for writing seismograms
  seismo_offset = 0
  seismo_current = 0

  if (myrank == 0) then
    write(IMAIN,*) "preparing attenuation..."
    call flush_IMAIN()
  endif
  call prepare_timerun_attenuation()



  allocate(kappastore(NGLLX,NGLLZ,nspec))
  allocate(mustore(NGLLX,NGLLZ,nspec))
  allocate(rhostore(NGLLX,NGLLZ,nspec))
  allocate(rho_vp(NGLLX,NGLLZ,nspec))
  allocate(rho_vs(NGLLX,NGLLZ,nspec))


if (assign_external_model) then

do ispec= 1,nspec
          do j = 1,NGLLZ
              do i = 1,NGLLX

              rhostore(i,j,ispec)    = rhoext(i,j,ispec)
              rho_vp(i,j,ispec)      = rhostore(i,j,ispec) * vpext(i,j,ispec)
              rho_vs(i,j,ispec)      = rhostore(i,j,ispec) * vsext(i,j,ispec)
              mustore(i,j,ispec)     = rho_vs(i,j,ispec) * vsext(i,j,ispec)
              kappastore(i,j,ispec)  = rho_vp(i,j,ispec) * vpext(i,j,ispec)-TWO*TWO*mustore(i,j,ispec)/3._CUSTOM_REAL

              enddo
          enddo
enddo

else ! Internal rho vp vs model

do ispec= 1,nspec
          do j = 1,NGLLZ
              do i = 1,NGLLX

              rhostore(i,j,ispec)       = density(1,kmato(ispec))
              lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
              mul_unrelaxed_elastic     = poroelastcoef(2,1,kmato(ispec))
              mustore(i,j,ispec)        = mul_unrelaxed_elastic
              kappastore(i,j,ispec)     = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic/3._CUSTOM_REAL
              rho_vp(i,j,ispec)         = density(1,kmato(ispec)) * sqrt((kappastore(i,j,ispec) + &
                                          4._CUSTOM_REAL*mul_unrelaxed_elastic/ &
                                          3._CUSTOM_REAL)/density(1,kmato(ispec)))
              rho_vs(i,j,ispec)         = density(1,kmato(ispec)) * sqrt(mul_unrelaxed_elastic/density(1,kmato(ispec)))
              enddo
        enddo
enddo

endif ! Internal/External model


  if (GPU_MODE) then

    call init_host_to_dev_variable()

    if (myrank == 0) then
      write(IMAIN,*) "preparing GPU..."
      call flush_IMAIN()
    endif
    call prepare_timerun_GPU()

  endif

  ! synchronizes all processes
  call synchronize_all()

  if (myrank == 0) then
    write(IMAIN,*) ""
    write(IMAIN,*) "done, preparation successful"
    write(IMAIN,*) ""
    call flush_IMAIN()
  endif

end subroutine prepare_timerun

