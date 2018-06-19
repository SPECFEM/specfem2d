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

  subroutine get_MPI()

! sets up MPI arrays

  use constants, only: IMAIN,USE_A_STRONG_FORMULATION_FOR_E1
  use shared_parameters, only: NPROC
  use specfem_par

  implicit none

  ! local parameters
  integer :: i,iinterface,imax_all,ier,n_sls_loc

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Setting up MPI communication arrays'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! initializes
  ninterface_acoustic = 0
  ninterface_elastic = 0
  ninterface_poroelastic = 0

  max_nibool_interfaces_ext_mesh = 0
  max_ibool_interfaces_size_ac = 0
  max_ibool_interfaces_size_el = 0
  max_ibool_interfaces_size_po = 0

  ! user output
  call max_all_i(ninterface,imax_all)
  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  maximum number of MPI interfaces (for a single slice) = ',imax_all
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! preparing for MPI communications
  call get_MPI_interfaces()

  if (NPROC > 1) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  master process:'
      write(IMAIN,*) '  number of MPI interfaces in acoustic domain    = ',ninterface_acoustic
      write(IMAIN,*) '  number of MPI interfaces in elastic domain     = ',ninterface_elastic
      write(IMAIN,*) '  number of MPI interfaces in poroelastic domain = ',ninterface_poroelastic
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! total MPI interfaces
    max_nibool_interfaces_ext_mesh = maxval(nibool_interfaces_ext_mesh(:))

    ! domain buffers for MPI communications
    max_ibool_interfaces_size_ac = maxval(nibool_interfaces_acoustic(:))
    max_ibool_interfaces_size_el = NDIM*maxval(nibool_interfaces_elastic(:))
    max_ibool_interfaces_size_po = NDIM*maxval(nibool_interfaces_poroelastic(:))

    if (ACOUSTIC_SIMULATION) then
      n_sls_loc = 0
      if (ATTENUATION_VISCOACOUSTIC .and. .not. USE_A_STRONG_FORMULATION_FOR_E1) n_sls_loc = N_SLS
      allocate(request_send_recv_acoustic(ninterface_acoustic*2),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array request_send_recv_acoustic')
      allocate(buffer_send_faces_vector_ac(max_ibool_interfaces_size_ac*(n_sls_loc+1),ninterface_acoustic),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array buffer_send_faces_vector_ac')
      allocate(buffer_recv_faces_vector_ac(max_ibool_interfaces_size_ac*(n_sls_loc+1),ninterface_acoustic),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array buffer_recv_faces_vector_ac')
    endif

    if (ELASTIC_SIMULATION) then
      allocate(request_send_recv_elastic(ninterface_elastic*2),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array request_send_recv_elastic')
      allocate(buffer_send_faces_vector_el(max_ibool_interfaces_size_el,ninterface_elastic),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array buffer_send_faces_vector_el')
      allocate(buffer_recv_faces_vector_el(max_ibool_interfaces_size_el,ninterface_elastic),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array buffer_recv_faces_vector_el')
    endif

    if (POROELASTIC_SIMULATION) then
      allocate(request_send_recv_poro(ninterface_poroelastic*4),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array request_send_recv_poro')
      allocate(buffer_send_faces_vector_pos(max_ibool_interfaces_size_po,ninterface_poroelastic),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array buffer_send_faces_vector_pos')
      allocate(buffer_recv_faces_vector_pos(max_ibool_interfaces_size_po,ninterface_poroelastic),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array buffer_recv_faces_vector_pos')
      allocate(buffer_send_faces_vector_pow(max_ibool_interfaces_size_po,ninterface_poroelastic),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array buffer_send_faces_vector_pow')
      allocate(buffer_recv_faces_vector_pow(max_ibool_interfaces_size_po,ninterface_poroelastic),stat=ier)
      if (ier /= 0) call stop_the_code('error in allocation of array buffer_recv_faces_vector_pow')
    endif

  else
    ! safety check
    if (myrank /= 0) call stop_the_code('Invalid myrank for serial simulation')

    ! user output
    write(IMAIN,*) '  This is a single process simulation, no need for MPI communication'
    write(IMAIN,*)
    call flush_IMAIN()
  endif ! end of test on whether there is more than one process (NPROC > 1)

  ! MPI interfaces arrays
  if (ninterface > 0) then

    allocate(ibool_interfaces_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface),stat=ier)
    if (ier /= 0) call stop_the_code('error in allocation of array ibool_interfaces_ext_mesh')
    ibool_interfaces_ext_mesh(:,:) = 0

    do iinterface = 1,ninterface
      do i = 1,nibool_interfaces_ext_mesh(iinterface)
        ibool_interfaces_ext_mesh(i,iinterface) = ibool_interfaces_ext_mesh_init(i,iinterface)
      enddo
    enddo
  else
    ! dummy array
    allocate(ibool_interfaces_ext_mesh(1,1))
  endif

  ! sets up inner and outer elements
  call get_MPI_setup_inner_outer_elements()

  ! setups up phase_inner arrays for looping over inner/outer elements in different domains (acoustic/elastic/..)
  ! note: must have flags ispec_is_elastic() etc. and flags ispec_is_inner() already set before calling this routine
  call get_MPI_phase_domains()

  end subroutine get_MPI

!
!-------------------------------------------------------------------------------------
!

  subroutine get_MPI_interfaces()

! sets up the MPI interface for communication between partitions

  use constants, only: CUSTOM_REAL,IMAIN,NGLLX,NGLLZ

  use shared_parameters, only: NPROC

  use specfem_par, only: nspec,ibool,nglob,ninterface,myrank,coord,ACOUSTIC_SIMULATION

  use specfem_par, only: ibool_interfaces_acoustic,ibool_interfaces_elastic, &
    ibool_interfaces_poroelastic, &
    nibool_interfaces_acoustic,nibool_interfaces_elastic, &
    nibool_interfaces_poroelastic, &
    ninterface_acoustic,ninterface_elastic,ninterface_poroelastic

  use specfem_par, only: buffer_send_faces_vector_ac,buffer_recv_faces_vector_ac,request_send_recv_acoustic

  implicit none

  !local parameters
  double precision, dimension(:), allocatable :: xp,zp
  double precision, dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: locval
  integer, dimension(:), allocatable :: nibool_interfaces_true

  ! for MPI buffers
  integer, dimension(:), allocatable :: reorder_interface,ind,ninseg,iwork
  integer, dimension(:), allocatable :: ibool_dummy
  logical, dimension(:), allocatable :: ifseg
  integer :: iinterface,ilocnum
  integer :: num_points1, num_points2

  ! assembly test
  integer :: i,j,ispec,iglob,inum,idomain
  integer :: max_nibool_interfaces,num_nibool,num_interface
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: test_flag_cr
  logical,dimension(:),allocatable :: mask_ibool

  ! checks if anything to do
  if (NPROC <= 1) return

  ! gets global indices for points on MPI interfaces
  ! (defined by my_interfaces) between different partitions
  ! and stores them in ibool_interfaces*** & nibool_interfaces*** (number of total points)
  call prepare_assemble_MPI()

  ! sorts ibool comm buffers lexicographically for all MPI interfaces
  num_points1 = 0
  num_points2 = 0
  allocate(nibool_interfaces_true(ninterface))

  ! acoustic/elastic/poroelastic domains
  do idomain = 1,3

    ! checks number of interface in this domain
    num_interface = 0
    if (idomain == 1) then
      num_interface = ninterface_acoustic
    else if (idomain == 2) then
      num_interface = ninterface_elastic
    else if (idomain == 3) then
      num_interface = ninterface_poroelastic
    endif
    if (num_interface == 0) cycle

    ! loops over interfaces
    do iinterface = 1, ninterface

      ! number of global points in this interface
      num_nibool = 0
      if (idomain == 1) then
        num_nibool = nibool_interfaces_acoustic(iinterface)
      else if (idomain == 2) then
        num_nibool = nibool_interfaces_elastic(iinterface)
      else if (idomain == 3) then
        num_nibool = nibool_interfaces_poroelastic(iinterface)
      endif
      ! checks if anything to sort
      if (num_nibool == 0) cycle

      allocate(xp(num_nibool))
      allocate(zp(num_nibool))
      allocate(locval(num_nibool))
      allocate(ifseg(num_nibool))
      allocate(reorder_interface(num_nibool))
      allocate(ibool_dummy(num_nibool))
      allocate(ind(num_nibool))
      allocate(ninseg(num_nibool))
      allocate(iwork(num_nibool))
      allocate(work(num_nibool))

      ! works with a copy of ibool array
      if (idomain == 1) then
        ibool_dummy(:) = ibool_interfaces_acoustic(1:num_nibool,iinterface)
      else if (idomain == 2) then
        ibool_dummy(:) = ibool_interfaces_elastic(1:num_nibool,iinterface)
      else if (idomain == 3) then
        ibool_dummy(:) = ibool_interfaces_poroelastic(1:num_nibool,iinterface)
      endif

      ! gets x,y,z coordinates of global points on MPI interface
      do ilocnum = 1, num_nibool
        iglob = ibool_dummy(ilocnum)
        xp(ilocnum) = coord(1,iglob)
        zp(ilocnum) = coord(2,iglob)
      enddo

      ! sorts (lexicographically?) ibool_interfaces and updates value
      ! of total number of points nibool_interfaces_true(iinterface)
      call sort_array_coordinates(num_nibool,xp,zp, &
                                  ibool_dummy, &
                                  reorder_interface,locval,ifseg, &
                                  nibool_interfaces_true(iinterface), &
                                  ind,ninseg,iwork,work)

      ! checks that number of MPI points are still the same
      num_points1 = num_points1 + num_nibool
      num_points2 = num_points2 + nibool_interfaces_true(iinterface)
      if (num_points1 /= num_points2) then
        write(IMAIN,*) 'Error sorting MPI interface points:',myrank
        write(IMAIN,*) '   domain:',idomain
        write(IMAIN,*) '   interface:',iinterface,num_points1,num_points2
        call exit_MPI(myrank,'error sorting MPI interface')
      endif

      ! stores new order of ibool array
      if (idomain == 1) then
        ibool_interfaces_acoustic(1:num_nibool,iinterface) = ibool_dummy(:)
      else if (idomain == 2) then
        ibool_interfaces_elastic(1:num_nibool,iinterface) = ibool_dummy(:)
      else if (idomain == 3) then
        ibool_interfaces_poroelastic(1:num_nibool,iinterface) = ibool_dummy(:)
      endif

      ! cleanup temporary arrays
      deallocate(xp)
      deallocate(zp)
      deallocate(locval)
      deallocate(ifseg)
      deallocate(reorder_interface)
      deallocate(ibool_dummy)
      deallocate(ind)
      deallocate(ninseg)
      deallocate(iwork)
      deallocate(work)
    enddo
  enddo

  ! cleanup
  deallocate(nibool_interfaces_true)

  ! outputs total number of MPI interface points
  call sum_all_i(num_points2,num_points1)
  if (myrank == 0) then
    write(IMAIN,*) '  total MPI interface points: ',num_points1
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! checks interfaces in acoustic domains
  if (ACOUSTIC_SIMULATION) then
    if (myrank == 0) then
      write(IMAIN,*) '  checking acoustic interfaces:'
      call flush_IMAIN()
    endif

    inum = 0
    if (ninterface_acoustic > 0) then
      ! checks with assembly of test fields
      allocate(test_flag_cr(nglob))
      test_flag_cr(:) = 0._CUSTOM_REAL

      ! sets identifier on all global nodes in this slice
      do ispec = 1, nspec
        ! sets flags on global points
        do j = 1, NGLLZ
          do i = 1, NGLLX
            ! global index
            iglob = ibool(i,j,ispec)
            ! sets identifier
            test_flag_cr(iglob) = myrank + 1.0
          enddo
        enddo
      enddo

      ! count global points on MPI interfaces to neighbors
      allocate(mask_ibool(nglob))
      mask_ibool(:) = .false.
      do iinterface = 1, ninterface
        do inum = 1,nibool_interfaces_acoustic(iinterface)
          iglob = ibool_interfaces_acoustic(inum,iinterface)
          if (.not. mask_ibool(iglob)) mask_ibool(iglob) = .true.
        enddo
      enddo
      inum = count(mask_ibool(:) .eqv. .true.)
      deallocate(mask_ibool)
    endif

    ! note: this MPI reduction awaits information from all processes.
    !          thus, avoid an MPI deadlock in case some of the paritions have no acoustic interface
    call sum_all_i(inum,num_points1)
    if (myrank == 0) then
      write(IMAIN,*) '  total number of global acoustic interface points: ',num_points1
      call flush_IMAIN()
    endif

    ! checks if assembly works
    inum = 0
    if (ninterface_acoustic > 0) then
      ! adds contributions from different partitions to flag arrays
      max_nibool_interfaces = maxval(nibool_interfaces_acoustic(:))

      allocate(request_send_recv_acoustic(ninterface_acoustic*2))
      allocate(buffer_send_faces_vector_ac(max_nibool_interfaces,ninterface_acoustic))
      allocate(buffer_recv_faces_vector_ac(max_nibool_interfaces,ninterface_acoustic))

      ! CUSTOM_REAL arrays
      call assemble_MPI_scalar_ac_blocking(test_flag_cr)

      ! checks number of interface points
      inum = 0
      do iglob = 1,nglob
        ! only counts flags with MPI contributions
        if (nint(test_flag_cr(iglob)) > myrank + 1 ) inum = inum + 1
      enddo

      deallocate(request_send_recv_acoustic)
      deallocate(buffer_send_faces_vector_ac)
      deallocate(buffer_recv_faces_vector_ac)
      deallocate(test_flag_cr)
    endif

    ! note: this MPI reduction awaits information from all processes.
    call sum_all_i(inum,num_points2)
    if (myrank == 0) then
      write(IMAIN,*) '  total number of global points assembled by acoustic MPI interfaces:',num_points2
      call flush_IMAIN()
    endif

    ! checks
    if (num_points1 /= num_points2) then
      call stop_the_code('Error acoustic MPI interfaces has invalid assembly')
    else
      if (myrank == 0) then
        write(IMAIN,*) '  interfaces okay'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
    endif

  endif

  end subroutine get_MPI_interfaces

!
!-------------------------------------------------------------------------------------
!

  subroutine get_MPI_setup_inner_outer_elements()

! sets up inner and outer elements for overlapping communication

  use constants, only: IMAIN,NGLLX,NGLLZ

  use specfem_par, only: nspec_inner,nspec_outer,ispec_is_inner, &
    ninterface,nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
    ibool,copy_ibool_ori,integer_mask_ibool, &
    myrank,nspec,nglob,NPROC

  implicit none

  ! local parameters
  integer :: ispec,i,j,iglob,iinterface
  ! inner/outer elements in the case of an MPI simulation
  integer :: ispec_inner,ispec_outer
  integer :: num_ispec_outer, num_ispec_inner

  logical,dimension(:),allocatable :: iglob_is_inner
  integer, dimension(:,:,:), allocatable :: ibool_outer,ibool_inner
  integer, dimension(:), allocatable  :: ispec_outer_to_glob, ispec_inner_to_glob

  real :: percentage_edge
  integer :: nglob_outer,nglob_inner

  ! initializes flags on elements (default to all inner elements)
  allocate(ispec_is_inner(nspec))
  ispec_is_inner(:) = .true.

  ! initializes flags on global nodes
  allocate(iglob_is_inner(nglob))
  iglob_is_inner(:) = .true.

  if (NPROC > 1) then
    ! sets flag on shared global nodes on MPI interfaces
    do iinterface = 1, ninterface
      do i = 1, nibool_interfaces_ext_mesh(iinterface)
        iglob = ibool_interfaces_ext_mesh(i,iinterface)
        iglob_is_inner(iglob) = .false.
      enddo
    enddo
  endif

  ! determines flags for inner elements (purely inside the partition) which contains an "outer" global node
  do ispec = 1, nspec
    do j = 1, NGLLZ
      do i = 1, NGLLX
        iglob = ibool(i,j,ispec)
        ispec_is_inner(ispec) = ( iglob_is_inner(iglob) .and. ispec_is_inner(ispec) )
      enddo
    enddo
  enddo

  ! counts number of outer elements (in this slice)
  nspec_inner = count(ispec_is_inner(:) .eqv. .true.)
  nspec_outer = nspec - nspec_inner

  ! building of corresponding arrays between inner/outer elements and their global number
  allocate(ispec_outer_to_glob(nspec_outer))
  allocate(ispec_inner_to_glob(nspec_inner))

  num_ispec_outer = 0
  num_ispec_inner = 0
  do ispec = 1, nspec
    if (ispec_is_inner(ispec)) then
      num_ispec_inner = num_ispec_inner + 1
      ispec_inner_to_glob(num_ispec_inner) = ispec
    else
      num_ispec_outer = num_ispec_outer + 1
      ispec_outer_to_glob(num_ispec_outer) = ispec
    endif
  enddo

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  number of outer elements  = ',nspec_outer
    write(IMAIN,*) '  number of inner elements  = ',nspec_inner
    write(IMAIN,*)

    percentage_edge = 100.*nspec_inner/real(nspec)
    write(IMAIN,*) '  percentage of outer elements ',100. - percentage_edge,'%'
    write(IMAIN,*) '  percentage of inner elements ',percentage_edge,'%'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! copies ibool into separate arrays for inner/outer elements
  allocate(ibool_outer(NGLLX,NGLLZ,nspec_outer))
  allocate(ibool_inner(NGLLX,NGLLZ,nspec_inner))

  ! loop over spectral elements
  do ispec_outer = 1,nspec_outer
    ! get global numbering for inner or outer elements
    ispec = ispec_outer_to_glob(ispec_outer)
    ibool_outer(:,:,ispec_outer) = ibool(:,:,ispec)
  enddo

  ! loop over spectral elements
  do ispec_inner = 1,nspec_inner
    ! get global numbering for inner or outer elements
    ispec = ispec_inner_to_glob(ispec_inner)
    ibool_inner(:,:,ispec_inner) = ibool(:,:,ispec)
  enddo

  ! reduces cache misses for outer elements
  if (nspec_outer > 0) then
    call get_global_indirect_addressing(nspec_outer,nglob,ibool_outer,copy_ibool_ori,integer_mask_ibool)

    ! the total number of points without multiples in this region is now known
    nglob_outer = maxval(ibool_outer)
  else
    nglob_outer = 0
  endif

  ! reduces cache misses for inner elements
  if (nspec_inner > 0) then
    call get_global_indirect_addressing(nspec_inner,nglob,ibool_inner,copy_ibool_ori,integer_mask_ibool)

    ! the total number of points without multiples in this region is now known
    nglob_inner = maxval(ibool_inner)
  else
    nglob_inner = 0
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  number of global nodes in outer elements  = ',nglob_outer
    write(IMAIN,*) '  number of global nodes in inner elements  = ',nglob_inner
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! frees temporary arrays
  deallocate(ibool_inner,ibool_outer)
  deallocate(ispec_inner_to_glob,ispec_outer_to_glob)
  deallocate(iglob_is_inner)

  end subroutine get_MPI_setup_inner_outer_elements


!
!-------------------------------------------------------------------------------------
!


  subroutine get_MPI_phase_domains()

  use constants, only: IMAIN
  use specfem_par

  implicit none

  ! local parameters
  integer :: ispec,ier
  integer :: ispec_inner,ispec_outer
  real :: percentage_edge

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  determining communication phases:'
    call flush_IMAIN()
  endif

  ! elastic domains
  nspec_inner_elastic = 0
  nspec_outer_elastic = 0

  ! only if this slice contains elastic elements
  if (any_elastic) then
    ! counts inner and outer elements
    do ispec = 1, nspec
      if (ispec_is_elastic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          nspec_inner_elastic = nspec_inner_elastic + 1
        else
          nspec_outer_elastic = nspec_outer_elastic + 1
        endif
      endif
    enddo

    ! stores indices of inner and outer elements
    num_phase_ispec_elastic = max(nspec_inner_elastic,nspec_outer_elastic)
    if (num_phase_ispec_elastic < 0 ) call stop_the_code('Error elastic simulation: num_phase_ispec_elastic is < zero')

    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if (ier /= 0 ) call stop_the_code('Error allocating array phase_ispec_inner_elastic')
    phase_ispec_inner_elastic(:,:) = 0

    ispec_inner = 0
    ispec_outer = 0
    do ispec = 1, nspec
      if (ispec_is_elastic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          ispec_inner = ispec_inner + 1
          phase_ispec_inner_elastic(ispec_inner,2) = ispec
        else
          ispec_outer = ispec_outer + 1
          phase_ispec_inner_elastic(ispec_outer,1) = ispec
        endif
      endif
    enddo
  else
    ! allocates dummy array
    num_phase_ispec_elastic = 0
    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if (ier /= 0 ) call stop_the_code('Error allocating dummy array phase_ispec_inner_elastic')
    phase_ispec_inner_elastic(:,:) = 0
  endif

  ! user output
  if (ELASTIC_SIMULATION) then
    call sum_all_i(nspec_outer_elastic,ispec_outer)
    call sum_all_i(nspec_inner_elastic,ispec_inner)
    if (myrank == 0) then
      ! check
      if (ispec_inner + ispec_outer == 0) call stop_the_code('Invalid total number of inner/outer elements for elastic simulation')
      ! ratio inner/outer
      percentage_edge = 100.0 * ispec_inner/real(ispec_inner + ispec_outer)
      ! output
      write(IMAIN,*) '  elastic domains:'
      write(IMAIN,*) '  total number of outer/inner elements = ',ispec_outer,ispec_inner
      write(IMAIN,*) '  total percentage of outer elements ',100. - percentage_edge,'%'
      write(IMAIN,*) '  total percentage of inner elements ',percentage_edge,'%'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! acoustic domains
  nspec_inner_acoustic = 0
  nspec_outer_acoustic = 0

  ! only if this slice contains acoustic elements
  if (any_acoustic) then

    ! counts inner and outer elements
    do ispec = 1, nspec
      if (ispec_is_acoustic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          nspec_inner_acoustic = nspec_inner_acoustic + 1
        else
          nspec_outer_acoustic = nspec_outer_acoustic + 1
        endif
      endif
    enddo

    ! stores indices of inner and outer elements for faster(?) computation
    num_phase_ispec_acoustic = max(nspec_inner_acoustic,nspec_outer_acoustic)
    if (num_phase_ispec_acoustic < 0 ) call stop_the_code('Error acoustic simulation: num_phase_ispec_acoustic is < zero')

    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if (ier /= 0 ) call stop_the_code('Error allocating array phase_ispec_inner_acoustic')
    phase_ispec_inner_acoustic(:,:) = 0

    ispec_inner = 0
    ispec_outer = 0
    do ispec = 1, nspec
      if (ispec_is_acoustic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          ispec_inner = ispec_inner + 1
          phase_ispec_inner_acoustic(ispec_inner,2) = ispec
        else
          ispec_outer = ispec_outer + 1
          phase_ispec_inner_acoustic(ispec_outer,1) = ispec
        endif
      endif
    enddo
  else
    ! allocates dummy array
    num_phase_ispec_acoustic = 0
    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if (ier /= 0 ) call stop_the_code('Error allocating dummy array phase_ispec_inner_acoustic')
    phase_ispec_inner_acoustic(:,:) = 0
  endif

  ! user output
  if (ACOUSTIC_SIMULATION) then
    ! master collects total
    call sum_all_i(nspec_outer_acoustic,ispec_outer)
    call sum_all_i(nspec_inner_acoustic,ispec_inner)
    if (myrank == 0) then
      ! check
      if (ispec_inner + ispec_outer == 0) call stop_the_code( &
'Invalid total number of inner/outer elements for acoustic simulation')
      ! ratio inner/outer
      percentage_edge = 100.0 * ispec_inner/real(ispec_inner + ispec_outer)
      ! output
      write(IMAIN,*) '  acoustic domains:'
      write(IMAIN,*) '  total number of outer/inner elements = ',ispec_outer,ispec_inner
      write(IMAIN,*) '  total percentage of outer elements ',100. - percentage_edge,'%'
      write(IMAIN,*) '  total percentage of inner elements ',percentage_edge,'%'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! poroelastic domains
  nspec_inner_poroelastic = 0
  nspec_outer_poroelastic = 0

  ! only if this slice contains elastic elements
  if (any_poroelastic) then
    ! counts inner and outer elements
    do ispec = 1, nspec
      if (ispec_is_poroelastic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          nspec_inner_poroelastic = nspec_inner_poroelastic + 1
        else
          nspec_outer_poroelastic = nspec_outer_poroelastic + 1
        endif
      endif
    enddo

    ! stores indices of inner and outer elements
    num_phase_ispec_poroelastic = max(nspec_inner_poroelastic,nspec_outer_poroelastic)
    if (num_phase_ispec_poroelastic < 0 ) call stop_the_code('Error poroelastic simulation: num_phase_ispec_poroelastic is < zero')

    allocate( phase_ispec_inner_poroelastic(num_phase_ispec_poroelastic,2),stat=ier)
    if (ier /= 0 ) call stop_the_code('Error allocating array phase_ispec_inner_poroelastic')
    phase_ispec_inner_poroelastic(:,:) = 0

    ispec_inner = 0
    ispec_outer = 0
    do ispec = 1, nspec
      if (ispec_is_poroelastic(ispec)) then
        if (ispec_is_inner(ispec) .eqv. .true.) then
          ispec_inner = ispec_inner + 1
          phase_ispec_inner_poroelastic(ispec_inner,2) = ispec
        else
          ispec_outer = ispec_outer + 1
          phase_ispec_inner_poroelastic(ispec_outer,1) = ispec
        endif
      endif
    enddo
  else
    ! allocates dummy array
    num_phase_ispec_poroelastic = 0
    allocate( phase_ispec_inner_poroelastic(num_phase_ispec_poroelastic,2),stat=ier)
    if (ier /= 0 ) call stop_the_code('Error allocating dummy array phase_ispec_inner_poroelastic')
    phase_ispec_inner_poroelastic(:,:) = 0
  endif

  ! user output
  if (POROELASTIC_SIMULATION) then
    call sum_all_i(nspec_outer_poroelastic,ispec_outer)
    call sum_all_i(nspec_inner_poroelastic,ispec_inner)
    if (myrank == 0) then
      ! check
      if (ispec_inner + ispec_outer == 0) call stop_the_code( &
'Invalid total number of inner/outer elements for poroelastic simulation')
      ! ratio inner/outer
      percentage_edge = 100.0 * ispec_inner/real(ispec_inner + ispec_outer)
      ! output
      write(IMAIN,*) '  poroelastic domains:'
      write(IMAIN,*) '  total number of outer/inner elements = ',ispec_outer,ispec_inner
      write(IMAIN,*) '  total percentage of outer elements ',100. - percentage_edge,'%'
      write(IMAIN,*) '  total percentage of inner elements ',percentage_edge,'%'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  end subroutine get_MPI_phase_domains
