
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

! XSMOOTH_SEM
!
! USAGE
!   mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR
!
!
! COMMAND LINE ARGUMENTS
!   SIGMA_H                - horizontal smoothing radius
!   SIGMA_V                - vertical smoothing radius
!   KERNEL_NAME            - kernel name, e.g. alpha_kernel
!   INPUT_DIR              - directory from which kernels are read
!   OUTPUT_DIR             - directory to which smoothed kernels are written
!
! DESCRIPTION
!   Smooths kernels by convolution with a Gaussian. Writes the resulting
!   smoothed kernels to OUTPUT_DIR.
!
!   Files written to OUTPUT_DIR have the suffix 'smooth' appended,
!   e.g. proc***alpha_kernel.bin becomes proc***alpha_kernel_smooth.bin
!
!   This program's primary use case is to smooth kernels. It can be used though on
!   any scalar field of dimension (NGLLX,NGLLZ,NSPEC).
!
!   This is a parallel program -- it must be invoked with mpirun or other
!   appropriate utility.  Operations are performed in embarassingly-parallel
!   fashion.


program smooth_sem

!! DK DK comments this out to avoid a compilation error  use mpi
  use postprocess_par

  implicit none
!! DK DK comments this out to avoid a compilation error  include "precision.h"
  integer, parameter :: NARGS = 5

  ! data must be of dimension: (NGLLX,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: dat
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dat_store,dat_smooth
!! DK DK comments this out to avoid a compilation error
!! DK DK  integer :: NGLOB_me, myrank,nproc, NGLOB_other, nspec_me, nspec_other
  integer :: myrank,nproc, NGLOB_other, nspec_me, nspec_other

!! DK DK comments this out to avoid a compilation error  integer :: i,j,iglob,ier,ispec2,ispec,inum,iker
  integer :: i,j,ier,ispec2,ispec,iker
!! DK DK comments this out to avoid a compilation error
!! DK DK  integer :: iproc, num_interfaces_ext_mesh, filesize, ninterface, max_interface_size
  integer :: iproc, ninterface, max_interface_size
  integer, dimension(:,:),allocatable :: ibool_interfaces_ext_mesh
  integer, dimension(:),allocatable :: nelmnts_neighbours, nibool_interfaces_ext_mesh

  character(len=MAX_STRING_LEN) :: arg(5)
  character(len=MAX_STRING_LEN) :: input_dir, output_dir
!! DK DK comments this out to avoid a compilation error  character(len=MAX_STRING_LEN) :: prname_lp, prname, local_path
  character(len=MAX_STRING_LEN) :: prname
!! DK DK comments this out to avoid a compilation error  character(len=MAX_STRING_LEN*2) :: local_data_file


  character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  integer :: nker

  ! smoothing parameters
  character(len=MAX_STRING_LEN*2) :: ks_file

  real(kind=CUSTOM_REAL) :: sigma_h, sigma_h2, sigma_h3, sigma_v, sigma_v2, sigma_v3
  real(kind=CUSTOM_REAL) :: norm_h, norm_v
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: norm, max_old, max_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: exp_val !,factor

  integer, dimension(:,:,:),allocatable :: ibool_other
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: tk
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: bk
!! DK DK comments this out to avoid a compilation error
!! DK DK  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: xl,  zl
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: xstore_me, zstore_me,xstore_other, zstore_other
  integer, dimension(:,:,:),allocatable :: imask

  real(kind=CUSTOM_REAL) :: dist_h,dist_v
  real(kind=CUSTOM_REAL) :: element_size
  REAL(kind=CUSTOM_REAL), PARAMETER :: PI = 3.1415927

!! DK DK comments this out to avoid a compilation error  call MPI_INIT(ier)
!! DK DK comments this out to avoid a compilation error  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ier)
!! DK DK comments this out to avoid a compilation error  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  if (myrank == 0) print *,"Running XSMOOTH_SEM"
!! DK DK comments this out to avoid a compilation error  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  stop 'there is a bug here thus DK DK adds this stop statement for now'

  ! parse command line arguments
  if (command_argument_count() /= NARGS) then
    if (myrank == 0) then
        print *, 'USAGE:  mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR'
      stop ' Please check command line arguments'
    endif
  endif
!! DK DK comments this out to avoid a compilation error  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
  enddo

  read(arg(1),*) sigma_h
  read(arg(2),*) sigma_v
  kernel_names_comma_delimited = arg(3)
  input_dir= arg(4)
  output_dir = arg(5)

  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  allocate(norm(nker),max_new(nker),max_old(nker))

!! DK DK comments this out to avoid a compilation error  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  ! check smoothing radii
  sigma_h2 = 2.0 * sigma_h ** 2  ! factor two for gaussian distribution with standard variance sigma
  sigma_v2 = 2.0 * sigma_v ** 2

  if (sigma_h2 < 1.e-18) stop 'Error sigma_h2 zero, must non-zero'
  if (sigma_v2 < 1.e-18) stop 'Error sigma_v2 zero, must non-zero'

  ! adds margin to search radius
  element_size = max(sigma_h,sigma_v) * 0.5

  ! search radius
  sigma_h3 = 3.0  * sigma_h + element_size
  sigma_v3 = 3.0  * sigma_v + element_size

  ! theoretic normal value
  ! (see integral over -inf to +inf of exp[- x*x/(2*sigma) ] = sigma * sqrt(2*pi) )
  ! note: smoothing is using a gaussian (ellipsoid for sigma_h /= sigma_v),
  norm_h = 2.0*PI*sigma_h**2
  norm_v = sqrt(2.0*PI) * sigma_v

  ! user output
  if (myrank == 0) then
    print *,"command line arguments:"
    print *,"  smoothing sigma_h , sigma_v                : ",sigma_h,sigma_v
    ! scalelength: approximately S ~ sigma * sqrt(8.0) for a gaussian smoothing
    print *,"  smoothing scalelengths horizontal, vertical: ",sigma_h*sqrt(8.0),sigma_v*sqrt(8.0)
    print *,"  input dir : ",trim(input_dir)
    print *,"  output dir: ",trim(output_dir)
    print *
  endif

  write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_NSPEC_ibool.bin'
  open(IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error opening smoothed kernel file'
  read(IIN) nspec_me
  close(IIN)
  allocate(xstore_me(NGLLX,NGLLZ,NSPEC_me),zstore_me(NGLLX,NGLLZ,NSPEC_me),stat=ier)


    write(prname, '(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_x.bin'
    ! gets the coordinate x of the points located in my slice
    open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      stop 'Error reading neighbors external mesh file'
    endif
    ! global point arrays
    read(IIN) xstore_me
    close(IIN)

    write(prname, '(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_z.bin'
    ! gets the coordinate z of the points located in my slice
    open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      stop 'Error reading neighbors external mesh file'
    endif
    ! global point arrays
    read(IIN) zstore_me
    close(IIN)

  ! synchronizes
!! DK DK comments this out to avoid a compilation error  call MPI_BARRIER(MPI_COMM_WORLD,ier)

! loops over slices
! each process reads in his own neighbor slices and gaussian filters the values
  allocate(tk(NGLLX,NGLLZ,NSPEC_ME,nker), bk(NGLLX,NGLLZ,NSPEC_ME),stat=ier)
  if (ier /= 0) stop 'Error allocating array tk and bk'

  tk = 0.0_CUSTOM_REAL
  bk = 0.0_CUSTOM_REAL




  do iproc = 0,NPROC-1
    ! slice database file
    write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',iproc,'_NSPEC_ibool.bin'
    open(IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) stop 'Error opening ibool file'
    read(IIN) nspec_other
    allocate(ibool_other(NGLLX,NGLLZ,nspec_other))
    read(IIN) ibool_other
    close(IIN)
    allocate(xstore_other(NGLLX,NGLLZ,NSPEC_other),zstore_other(NGLLX,NGLLZ,NSPEC_other),stat=ier)
    if (ier /= 0) stop 'Error allocating array xstore_other etc.'
    nglob_other = maxval(ibool_other(:,:,:))
    allocate(imask(NGLLX,NGLLZ,nglob_other))


     write(prname, '(a,i6.6,a)') trim(input_dir)//'/proc',iproc,'_x.bin'
    ! gets the coordinate x of the points located in the other slice
    open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      stop 'Error reading x coordinate'
    endif
    ! global point arrays
    read(IIN) xstore_other
    close(IIN)

    write(prname, '(a,i6.6,a)') trim(input_dir)//'/proc',iproc,'_z.bin'
    ! gets the coordinate z of the points located in the other slice
    open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      stop 'Error reading z coordinate'
    endif
    ! global point arrays
    read(IIN) zstore_other
    close(IIN)

    allocate(dat(NGLLX,NGLLZ,NSPEC_other),dat_store(NGLLX,NGLLZ,NSPEC_other,nker),stat=ier)
    if (ier /= 0) stop 'Error allocating dat array'


    write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',iproc,'_MPI_interfaces_info.bin'
    open(unit=172,file=prname,status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      stop 'Error reading MPI info'
    endif
    read(172) ninterface
    allocate(nelmnts_neighbours(ninterface),nibool_interfaces_ext_mesh(ninterface))
    read(172) nelmnts_neighbours
    read(172) nibool_interfaces_ext_mesh
    read(172) max_interface_size
    allocate(ibool_interfaces_ext_mesh(NGLLX*max_interface_size,ninterface))
    read(172) ibool_interfaces_ext_mesh
    close(172)



do iker= 1, nker
    ! data file
    write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',iproc,'_'//trim(kernel_names(iker))//'.bin'

    open(unit = IIN,file = trim(prname),status='old',action='read',form ='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening data file: ',trim(prname)
      stop 'Error opening data file'
    endif
    read(IIN) dat
    close(IIN)

dat_store(:,:,:,iker) = dat(:,:,:)

    if (iproc == myrank) max_old(iker) = maxval(abs(dat(:,:,:)))

enddo
    ! loop over elements to be smoothed in the current slice
    do ispec = 1, nspec_me
      ! --- only double loop over the elements in the search radius ---
      imask(:,:,:)= -1
      do ispec2 = 1, nspec_other

        ! calculates horizontal and vertical distance between two element centers
        call get_distance_square_vec(dist_h,dist_v,xstore_me(1,1,ispec),zstore_me(1,1,ispec),&
                          xstore_other(1,1,ispec2),zstore_other(1,1,ispec2))

        ! checks distance between centers of elements
        if (dist_h > sigma_h3**2 .or. dist_v > sigma_v3**2) cycle

    ! loop over GLL points of the elements in current slice (ispec)
        do j = 1, NGLLZ
            do i = 1, NGLLX

              ! calculate weights based on gaussian smoothing
              exp_val = 0.0_CUSTOM_REAL
              call smoothing_weights_vec(xstore_me(i,j,ispec),zstore_me(i,j,ispec),sigma_h2,sigma_v2,exp_val,&
                      xstore_other(:,:,ispec2),zstore_other(:,:,ispec2),nglob_other,imask,ispec2,nspec_other,ibool_other, i, j, &
                      ninterface,nelmnts_neighbours, nibool_interfaces_ext_mesh,max_interface_size,ibool_interfaces_ext_mesh, iproc)

              ! adds contribution of element ispec2 to smoothed kernel values
              do iker=1, nker
              tk(i,j,ispec,iker) = tk(i,j,ispec,iker) + sum(exp_val(:,:) * dat_store(:,:,ispec2,iker))
              enddo
              ! normalization, integrated values of gaussian smoothing function
              bk(i,j,ispec) = bk(i,j,ispec) + sum(exp_val(:,:))

          enddo
        enddo

      enddo ! iglob2
    enddo ! iglob

    ! frees arrays
    deallocate(dat,dat_store,xstore_other,zstore_other,imask,ibool_other,nelmnts_neighbours, &
               nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh)

  enddo ! iproc
  if (myrank == 0) print *

  ! normalizes/scaling factor
  if (myrank == 0) print *, 'Scaling values: min/max = ',minval(bk),maxval(bk)

  allocate(dat_smooth(NGLLX,NGLLZ,NSPEC_me,nker),stat=ier)
  if (ier /= 0) stop 'Error allocating array dat_smooth'

  dat_smooth(:,:,:,:) = 0.0_CUSTOM_REAL
  do ispec = 1, nspec_me
      do j = 1, NGLLZ
        do i = 1, NGLLX
    !      if (abs(bk(i,j,ispec)) < 1.e-18) then
    !        print *, 'Problem norm here --- ', ispec, i, j, bk(i,j,ispec), norm
     !     endif
          ! normalizes smoothed kernel values by integral value of gaussian weighting
          dat_smooth(i,j,ispec,:) = tk(i,j,ispec,:) / bk(i,j,ispec)
        enddo
      enddo
  enddo !  ispec
  deallocate(tk,bk)

do iker= 1, nker

  max_new(iker) = maxval(abs(dat_smooth(:,:,:,iker)))

  ! file output
  ! smoothed kernel file name


  write(ks_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//trim(kernel_names(iker))//'_smooth.bin'

  open(IOUT,file=trim(ks_file),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0) stop 'Error opening smoothed kernel file'
  write(IOUT) dat_smooth(:,:,:,iker)
  close(IOUT)
  if (myrank == 0) print *,'written: ',trim(ks_file)
enddo
  ! frees memory
  deallocate(dat_smooth)

  ! synchronizes
!! DK DK comments this out to avoid a compilation error  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  ! the maximum value for the smoothed kernel
  norm(:) = max_old(:)

!! DK DK comments this out to avoid a compilation error
!! DK DK call MPI_REDUCE(norm,max_old,nker,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  norm(:) = max_new(:)
!! DK DK comments this out to avoid a compilation error
!! DK DK call MPI_REDUCE(norm,max_new,nker,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
do iker= 1, nker
  if (myrank == 0) then
    print *
    print *,'  Maximum data value before smoothing = ', max_old(iker), 'for ', trim(kernel_names(iker))
    print *,'  Maximum data value after smoothing  = ', max_new(iker), 'for ', trim(kernel_names(iker))

  endif
enddo
   if (myrank == 0) close(IIN)

  ! stop all the processes and exit
!! DK DK comments this out to avoid a compilation error  call MPI_FINALIZE(ier)

end program smooth_sem

!
! -----------------------------------------------------------------------------
!
  subroutine smoothing_weights_vec(x0,z0,sigma_h2,sigma_v2,exp_val,&
                              xx_elem,zz_elem,nglob_other,imask,ispec,nspec,ibool,i,j,&
                      ninterface,nelmnts_neighbours, nibool_interfaces_ext_mesh,max_interface_size,ibool_interfaces_ext_mesh,iproc)

  use constants
  implicit none

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ),intent(out) :: exp_val
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ),intent(in) :: xx_elem, zz_elem
  real(kind=CUSTOM_REAL),intent(in) :: x0,z0,sigma_h2,sigma_v2
  integer :: nglob_other, ispec, nspec, i, j,max_interface_size, ninterface, iproc
  integer,dimension(NGLLX,NGLLZ,nglob_other) :: imask
  integer,dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer,dimension(ninterface) :: nelmnts_neighbours,  nibool_interfaces_ext_mesh
  integer,dimension(NGLLX*max_interface_size,ninterface) :: ibool_interfaces_ext_mesh

  ! local parameters
  integer :: ii,jj,iglob,k, iinterface, i_former_slice
  real(kind=CUSTOM_REAL) :: dist_h,dist_v
  real(kind=CUSTOM_REAL) :: sigma_h2_inv,sigma_v2_inv
  integer :: nb_former_neighbour_slice

  nb_former_neighbour_slice = 0

  sigma_h2_inv = 1.0_CUSTOM_REAL / sigma_h2
  sigma_v2_inv = 1.0_CUSTOM_REAL / sigma_v2

    do k=1, ninterface
    if (nelmnts_neighbours(k)<iproc)   nb_former_neighbour_slice = nb_former_neighbour_slice + 1
    enddo

    do jj = 1, NGLLZ
      do ii = 1, NGLLX

    iglob = ibool(ii,jj,ispec)

     ! Check if the point has not already be taken in account in an other element of the same slice
     if (imask(i,j,iglob)==1) cycle

    ! Check if the point has not already be taken in account in an other element of an other slice
    ! (located in an interface between MPI slices)
     if(nb_former_neighbour_slice>0) then
     iinterface=1
     do i_former_slice = 1,nb_former_neighbour_slice
       do while (nelmnts_neighbours(iinterface)<iproc)
          iinterface= iinterface+1
       enddo
       do k = 1, nibool_interfaces_ext_mesh(iinterface)
          if (iglob==ibool_interfaces_ext_mesh(k,i_former_slice)) imask(i,j,iglob)=1
       enddo
     enddo
       if(imask(i,j,iglob)==1) cycle
     endif
        ! point in second slice
        ! gets vertical and horizontal distance
        call get_distance_square_vec(dist_h,dist_v,x0,z0,xx_elem(ii,jj),zz_elem(ii,jj))
        ! gaussian function
        exp_val(ii,jj) = exp(- sigma_h2_inv*dist_h - sigma_v2_inv*dist_v)
        ! mark the current gll point as taken in account
        imask(i,j,iglob)=1
      enddo
    enddo

  end subroutine smoothing_weights_vec


!
! -----------------------------------------------------------------------------
!

  subroutine get_distance_square_vec(dist_h,dist_v,x0,z0,x1,z1)

! returns vector lengths as distances in radial and horizontal direction
! only for flat earth with z in vertical direction

  use constants
  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: dist_h,dist_v
  real(kind=CUSTOM_REAL),intent(in) :: x0,z0,x1,z1

  ! vertical distance
  dist_v =  (z0-z1)**2

  ! horizontal distance
  dist_h =  (x0-x1)**2

  end subroutine get_distance_square_vec
