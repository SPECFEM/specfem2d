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

! XSMOOTH_SEM
!
! USAGE
!   mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR GPU_MODE
!
!
! COMMAND LINE ARGUMENTS
!   SIGMA_H                - horizontal smoothing radius
!   SIGMA_V                - vertical smoothing radius
!   KERNEL_NAME            - kernel name, e.g. alpha_kernel
!   INPUT_DIR              - directory from which kernels are read
!   OUTPUT_DIR             - directory to which smoothed kernels are written
!   GPU_MODE               - use GPUs to process smoothing (logical T F)

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



#ifdef USE_MPI
  use mpi
#endif

  use postprocess_par

  implicit none

#ifdef USE_MPI
  include "precision.h"
#endif

  integer, parameter :: NARGS = 6

  ! data must be of dimension: (NGLLX,NGLLZ,NSPEC_AB)
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: dat
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: dat_store,dat_smooth
  integer :: NGLOB_me, nspec_me, nspec_other, ncuda_devices
  ! MPI
  integer :: myrank,NPROC
  integer(kind=8) :: Container

  integer :: i,j,ier,ispec2,ispec, iker, iglob
  integer :: iproc

  character(len=MAX_STRING_LEN) :: arg(6)
  character(len=MAX_STRING_LEN) :: input_dir, output_dir
  character(len=MAX_STRING_LEN) :: prname

  logical :: GPU_MODE

  character(len=MAX_STRING_LEN) :: kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  integer :: nker

  ! smoothing parameters
  character(len=MAX_STRING_LEN*2) :: ks_file

  real(kind=CUSTOM_REAL) :: sigma_h, sigma_h2_inv, sigma_h3_sq, sigma_v,sigma_v2_inv, sigma_v3_sq
  real(kind=CUSTOM_REAL) :: norm_h, norm_v
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: norm, max_old,max_new, min_old, min_new
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: exp_val, factor, wgll_sq
  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLX) :: xigll
  integer, dimension(:,:,:),allocatable :: ibool_me
  integer, dimension(:),allocatable :: imask
  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: tk
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: bk
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: xstore_me, zstore_me,xstore_other, zstore_other, jacobian


  real(kind=CUSTOM_REAL) :: dist_h,dist_v
  real(kind=CUSTOM_REAL) :: element_size
  real(kind=CUSTOM_REAL), PARAMETER :: PI = 3.1415927
  real t1,t2

  ! MPI initialization
  call init_mpi()
  call world_size(NPROC)
  call world_rank(myrank)

  if (myrank == 0) print *,"Running XSMOOTH_SEM on",NPROC,"processors"
  call cpu_time(t1)

  ! parse command line arguments
  if (command_argument_count() /= NARGS) then
    if (myrank == 0) then
        print *, 'USAGE:  mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR GPU_MODE'
      call stop_the_code(' Please check command line arguments')
    endif
  endif

  ! synchronizes all processes
  call synchronize_all()

  do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
  enddo

  read(arg(1),*) sigma_h
  read(arg(2),*) sigma_v
  kernel_names_comma_delimited = arg(3)
  input_dir= arg(4)
  output_dir = arg(5)
  read(arg(6),*) GPU_MODE

  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  allocate(norm(nker),max_new(nker),max_old(nker),min_new(nker),min_old(nker))

 if (GPU_MODE) call initialize_cuda_device(myrank,ncuda_devices)

  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  !We assume NGLLX=NGLLZ
  do j = 1,NGLLZ
    do i = 1,NGLLX
      wgll_sq(i,j) = real(wxgll(i)*wxgll(j),kind=CUSTOM_REAL)
    enddo
  enddo

  ! check smoothing radii
  sigma_h2_inv = ( 1.0 / (2.0 * (sigma_h ** 2)) ) ! factor two for Gaussian distribution with standard variance sigma
  sigma_v2_inv = ( 1.0 / (2.0 * (sigma_v ** 2)) )

  if ((1.0 / sigma_h2_inv) < 1.e-18) call stop_the_code('Error sigma_h2 zero, must non-zero')
  if ((1.0 / sigma_v2_inv) < 1.e-18) call stop_the_code('Error sigma_v2 zero, must non-zero')

  ! adds margin to search radius
  element_size = max(sigma_h,sigma_v) * 0.5

  ! search radius
  sigma_h3_sq = (3.0  * sigma_h + element_size)**2
  sigma_v3_sq = (3.0  * sigma_v + element_size)**2

  ! theoretic normal value
  ! (see integral over -inf to +inf of exp[- x*x/(2*sigma) ] = sigma * sqrt(2*pi) )
  ! note: smoothing is using a Gaussian (ellipsoid for sigma_h /= sigma_v),
  norm_h = 2.0*PI*sigma_h**2
  norm_v = sqrt(2.0*PI) * sigma_v

  ! user output
  if (myrank == 0) then
    print *,"command line arguments:"
    print *,"  smoothing sigma_h , sigma_v                : ",sigma_h,sigma_v
    ! scalelength: approximately S ~ sigma * sqrt(8.0) for a Gaussian smoothing
    print *,"  smoothing scalelengths horizontal, vertical: ",sigma_h*sqrt(8.0),sigma_v*sqrt(8.0)
    print *,"  input dir : ",trim(input_dir)
    print *,"  output dir: ",trim(output_dir)
    print *
  endif


  write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_NSPEC_ibool.bin'
  open(IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      call stop_the_code('Error opening _NSPEC_IBOOL file')
  endif
  read(IIN) nspec_me
  allocate(ibool_me(NGLLX,NGLLZ,nspec_me))
  read(IIN) ibool_me
  close(IIN)
  nglob_me = maxval(ibool_me(:,:,:))
  allocate(xstore_me(NGLLX,NGLLZ,NSPEC_me),zstore_me(NGLLX,NGLLZ,NSPEC_me),imask(nglob_me),stat=ier)

    write(prname, '(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_x.bin'
    ! gets the coordinate x of the points located in my slice
    open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      call stop_the_code('Error reading neighbors external mesh file')
    endif
    ! global point arrays
    read(IIN) xstore_me
    close(IIN)

    write(prname, '(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_z.bin'
    ! gets the coordinate z of the points located in my slice
    open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      call stop_the_code('Error reading neighbors external mesh file')
    endif
    ! global point arrays
    read(IIN) zstore_me
    close(IIN)


  ! synchronizes all processes
  call synchronize_all()

  ! GPU setup
  if (GPU_MODE) then
    call prepare_arrays_GPU(Container,xstore_me,zstore_me, &
                            sigma_h2_inv,sigma_v2_inv,sigma_h3_sq,sigma_v3_sq,nspec_me,nker,wgll_sq)

    ! synchronizes all processes
    call synchronize_all()
  endif


! loops over slices
! each process reads all the other slices and Gaussian filters the values
  allocate(tk(nglob_me,nker), bk(nglob_me),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating array tk and bk')

  tk = 0.0_CUSTOM_REAL
  bk = 0.0_CUSTOM_REAL

  do iproc = 0,NPROC-1
    ! slice database file
    write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',iproc,'_NSPEC_ibool.bin'
    open(IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening ibool file')
    read(IIN) nspec_other
    close(IIN)
    allocate(xstore_other(NGLLX,NGLLZ,NSPEC_other),zstore_other(NGLLX,NGLLZ,NSPEC_other), &
             jacobian(NGLLX,NGLLZ,NSPEC_other),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating array xstore_other etc.')

     write(prname, '(a,i6.6,a)') trim(input_dir)//'/proc',iproc,'_x.bin'
    ! gets the coordinate x of the points located in the other slice
    open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      call stop_the_code('Error reading x coordinate')
    endif
    ! global point arrays
    read(IIN) xstore_other
    close(IIN)

    write(prname, '(a,i6.6,a)') trim(input_dir)//'/proc',iproc,'_z.bin'
    ! gets the coordinate z of the points located in the other slice
    open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      call stop_the_code('Error reading z coordinate')
    endif
    ! global point arrays
    read(IIN) zstore_other
    close(IIN)

    write(prname, '(a,i6.6,a)') trim(input_dir)//'/proc',iproc,'_jacobian.bin'
    ! gets the jacobian of the points located in the other slice
    open(unit=IIN,file=trim(prname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: could not open database file: ',trim(prname)
      call stop_the_code('Error reading jacobian')
    endif
    ! global point arrays
    read(IIN) jacobian
    close(IIN)

    allocate(dat(NGLLX,NGLLZ,NSPEC_other),dat_store(NGLLX,NGLLZ,NSPEC_other,nker),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating dat array')

    do iker= 1, nker
      ! data file
      write(prname,'(a,i6.6,a)') trim(input_dir)//'/proc',iproc,'_'//trim(kernel_names(iker))//'.bin'

      open(unit = IIN,file = trim(prname),status='old',action='read',form ='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'Error opening data file: ',trim(prname)
        call stop_the_code('Error opening data file')
      endif
      read(IIN) dat
      close(IIN)

      dat_store(:,:,:,iker) = dat(:,:,:)

      if (iproc == myrank) then
        max_old(iker) = maxval(abs(dat(:,:,:)))
        min_old(iker) = minval(abs(dat(:,:,:)))
      endif
    enddo

    imask = -1

    if (GPU_MODE) then
      call compute_smooth(Container,jacobian,xstore_other,zstore_other,dat_store,nspec_other)
    else
      ! loop over elements to be smoothed in the current slice
      do ispec = 1, nspec_me
        ! --- only double loop over the elements in the search radius ---
        do ispec2 = 1, nspec_other

          ! calculates horizontal and vertical distance between two element centers
          call get_distance_square_vec(dist_h,dist_v,xstore_me(1,1,ispec),zstore_me(1,1,ispec), &
                            xstore_other(1,1,ispec2),zstore_other(1,1,ispec2))

          ! checks distance between centers of elements
     !     if (dist_h > 2*sigma_h3_sq .or. dist_v > 2*sigma_v3_sq) cycle



          factor(:,:) = jacobian(:,:,ispec2) * wgll_sq(:,:)

      ! loop over GLL points of the elements in current slice (ispec)
          do j = 1, NGLLZ
              do i = 1, NGLLX
                iglob = ibool_me(i,j,ispec)
                if (imask(iglob) == 1) cycle

                ! calculate weights based on Gaussian smoothing
                exp_val = 0.0_CUSTOM_REAL

                call smoothing_weights_vec(xstore_me(i,j,ispec),zstore_me(i,j,ispec),sigma_h2_inv,sigma_v2_inv,exp_val, &
                        xstore_other(:,:,ispec2),zstore_other(:,:,ispec2))

                exp_val(:,:) = exp_val(:,:) * factor(:,:)

                ! adds contribution of element ispec2 to smoothed kernel values
                do iker= 1, nker
                tk(iglob,iker) = tk(iglob,iker) + sum(exp_val(:,:) * dat_store(:,:,ispec2,iker))
                enddo
                ! normalization, integrated values of Gaussian smoothing function
                bk(iglob) = bk(iglob) + sum(exp_val(:,:))

              enddo
          enddo
        enddo ! ispec2

          do j = 1, NGLLZ
              do i = 1, NGLLX
                imask(ibool_me(i,j,ispec)) = 1
              enddo
          enddo

      enddo ! ispec

    endif !GPU_MODE

    ! frees arrays
    deallocate(dat,dat_store,xstore_other,zstore_other,jacobian)

  enddo ! iproc

  ! normalizes/scaling factor
  if (myrank == 0) print *
  if (myrank == 0) print *, 'Scaling values: min/max = ',minval(bk),maxval(bk)

  allocate(dat_smooth(NGLLX,NGLLZ,NSPEC_me,nker),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating array dat_smooth')

  dat_smooth(:,:,:,:) = 0.0_CUSTOM_REAL

  if (GPU_MODE) then
    call get_smooth(Container,dat_smooth)
  else
    do ispec = 1, nspec_me
        do j = 1, NGLLZ
          do i = 1, NGLLX
   !         if (abs(bk(i,j,ispec)) < 1.e-18) then
   !           print *, 'Problem norm here --- ', ispec, i, j, bk(i,j,ispec)
    !        endif
                iglob = ibool_me(i,j,ispec)
            ! normalizes smoothed kernel values by integral value of Gaussian weighting
            dat_smooth(i,j,ispec,:) = tk(iglob,:) / bk(iglob)
          enddo
        enddo
    enddo !  ispec
  endif

  deallocate(tk,bk)

  do iker= 1, nker
    max_new(iker) = maxval(abs(dat_smooth(:,:,:,iker)))
    min_new(iker) = minval(abs(dat_smooth(:,:,:,iker)))
    ! file output
    ! smoothed kernel file name
    write(ks_file,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//trim(kernel_names(iker))//'_smooth.bin'

    open(IOUT,file=trim(ks_file),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call stop_the_code('Error opening smoothed kernel file')
    write(IOUT) dat_smooth(:,:,:,iker)
    close(IOUT)
    if (myrank == 0) print *,'written: ',trim(ks_file)
  enddo

  ! frees memory
  deallocate(dat_smooth)

  ! synchronizes all processes
  call synchronize_all()

#ifdef USE_MPI
  if (NPROC > 1) then

    ! the maximum value for the smoothed kernel
    norm(:) = max_old(:)


    call MPI_REDUCE(norm,max_old,nker,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

    norm(:) = max_new(:)

    call MPI_REDUCE(norm,max_new,nker,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)
    norm(:) = min_old(:)


    call MPI_REDUCE(norm,min_old,nker,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

    norm(:) = min_new(:)

    call MPI_REDUCE(norm,min_new,nker,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  endif
#endif

  do iker= 1, nker
    if (myrank == 0) then
      print *
      print *,' Min / Max data value before smoothing = ', min_old(iker), max_old(iker), 'for ', trim(kernel_names(iker))
      print *,' Min / Max data value after smoothing  = ', min_new(iker), max_new(iker), 'for ', trim(kernel_names(iker))

    endif
  enddo

  call cpu_time(t2)

  if (GPU_Mode) then
    print *,'Computation time with GPU:',t2-t1
  else
    print *,'Computation time with CPU:',t2-t1
  endif

  if (myrank == 0) close(IIN)

  ! MPI finish
  call finalize_mpi()

end program smooth_sem

!
! -----------------------------------------------------------------------------
!
  subroutine smoothing_weights_vec(x0,z0,sigma_h2_inv,sigma_v2_inv,exp_val,xx_elem,zz_elem)

  use constants
  implicit none

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ),intent(out) :: exp_val
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ),intent(in) :: xx_elem, zz_elem
  real(kind=CUSTOM_REAL),intent(in) :: x0,z0,sigma_h2_inv,sigma_v2_inv


  ! local parameters
  integer :: ii,jj
  real(kind=CUSTOM_REAL) :: dist_h,dist_v


  do jj = 1, NGLLZ
    do ii = 1, NGLLX
      ! gets vertical and horizontal distance
      call get_distance_square_vec(dist_h,dist_v,x0,z0,xx_elem(ii,jj),zz_elem(ii,jj))
       ! Gaussian function
      exp_val(ii,jj) = exp(- sigma_h2_inv*dist_h - sigma_v2_inv*dist_v)
    enddo
  enddo

  end subroutine smoothing_weights_vec


!
! -----------------------------------------------------------------------------
!

  subroutine get_distance_square_vec(dist_h,dist_v,x0,z0,x1,z1)

! returns square lengths as distances in radial and horizontal direction
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
