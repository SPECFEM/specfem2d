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


! sum_kernels
!
! this program can be used for event kernel summation,
! where it sums up transverse isotropic kernel files:
!
!   - proc***_reg1_bulk_c_kernel.bin
!   - proc***_reg1_bulk_betav_kernel.bin
!   - proc***_reg1_bulk_betah_kernel.bin
!   - proc***_reg1_eta_kernel.bin
!
! input file: kernels_list.txt
!   lists all event kernel directories which should be summed together
!
! input directory:  INPUT_KERNELS/
!    contains links to all event kernel directories (listed in "kernels_list.txt")
!
! output directory: OUTPUT_SUM/
!    the resulting kernel files will be stored in this directory
!
!
! DEPRECATION WARNING: Eventually, all of the following routines, or at lesast
! some the subroutines, will be merged with src/tomography/xcombine_sem
!


program sum_kernels

  use tomography_par,only: MAX_STRING_LEN,MAX_KERNEL_PATHS,KERNEL_FILE_LIST,IIN, &
    myrank,sizeprocs, &
    USE_ALPHA_BETA_RHO,USE_ISO_KERNELS

  implicit none

  character(len=MAX_STRING_LEN) :: kernel_list(MAX_KERNEL_PATHS)
  character(len=MAX_STRING_LEN) :: sline, kernel_name
  integer :: nker
  integer :: ier
  integer :: NPROC

  ! ============ program starts here =====================

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)
  NPROC = sizeprocs

  if (myrank==0) then
    write(*,*) 'sum_kernels:'
    write(*,*)
    write(*,*) 'reading kernel list: '
  endif
  call synchronize_all()

  ! reads in event list
  nker=0
  open(unit = IIN, file = trim(KERNEL_FILE_LIST), status = 'old',iostat = ier)
  if (ier /= 0) then
     print *,'Error opening ',trim(KERNEL_FILE_LIST),myrank
     stop 1
  endif
  do while (1 == 1)
     read(IIN,'(a)',iostat=ier) sline
     if (ier /= 0) exit
     nker = nker+1
     if (nker > MAX_KERNEL_PATHS) stop 'Error number of kernels exceeds MAX_KERNEL_PATHS'
     kernel_list(nker) = sline
  enddo
  close(IIN)
  if (myrank == 0) then
    write(*,*) '  ',nker,' events'
    write(*,*)
  endif

  ! checks if number of MPI process as specified
  if (sizeprocs /= NPROC) then
    if (myrank == 0) then
      print *,''
      print *,'Error: run xsum_kernels with the same number of MPI processes '
      print *,'       as specified in Par_file by NPROC when slices were created'
      print *,''
      print *,'for example: mpirun -np ',NPROC,' ./xsum_kernels ...'
      print *,''
    endif
    call synchronize_all()
    stop 'Error total number of slices'
  endif
  call synchronize_all()

  ! reads mesh file
  !
  ! needs to get array dimensions (NSPEC)
  call read_model_nspec()

  ! user output
  if (myrank == 0) then
    print *,'summing kernels in INPUT_KERNELS/ directories:'
    print *,kernel_list(1:nker)
    print *
  endif

  ! synchronizes
  call synchronize_all()

  ! sums up kernels
  if (USE_ISO_KERNELS) then

    !  isotropic kernels
    if (myrank == 0) write(*,*) 'isotropic kernels: bulk_c, bulk_beta, rho'

    kernel_name = 'bulk_c_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'bulk_beta_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'rho_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

  else if (USE_ALPHA_BETA_RHO) then

    ! isotropic kernels
    if (myrank == 0) write(*,*) 'isotropic kernels: alpha, beta, rho'

    kernel_name = 'alpha_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'beta_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'rho_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

  else

    ! transverse isotropic kernels
    if (myrank == 0) write(*,*) 'transverse isotropic kernels: bulk_c, bulk_betav, bulk_betah,eta'

    kernel_name = 'bulk_c_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'bulk_betav_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'bulk_betah_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

    kernel_name = 'eta_kernel'
    call sum_kernel(kernel_name,kernel_list,nker)

  endif

  if (myrank==0) write(*,*) 'done writing all kernels, see directory OUTPUT_SUM/'

  ! stop all the processes, and exit
  call finalize_mpi()

end program sum_kernels

!
!-------------------------------------------------------------------------------------------------
!

subroutine sum_kernel(kernel_name,kernel_list,nker)

  use tomography_par

  implicit none

  character(len=MAX_STRING_LEN) :: kernel_name,kernel_list(MAX_KERNEL_PATHS)
  integer :: nker

  ! local parameters
  character(len=MAX_STRING_LEN*2) :: k_file
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: kernel,total_kernel
  double precision :: norm,norm_sum
  integer :: iker,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: mask_source

  ! initializes arrays
  allocate(kernel(NGLLX,NGLLZ,NSPEC), &
           total_kernel(NGLLX,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel arrays'

  if (USE_SOURCE_MASK) then
    allocate( mask_source(NGLLX,NGLLZ,NSPEC) )
    mask_source(:,:,:) = 1.0_CUSTOM_REAL
  endif

  ! loops over all event kernels
  total_kernel = 0._CUSTOM_REAL
  do iker = 1, nker
    ! user output
    if (myrank==0) then
      write(*,*) 'reading in event kernel for: ',trim(kernel_name)
      write(*,*) '    ',iker, ' out of ', nker
    endif

    ! sensitivity kernel / frechet derivative
    kernel = 0._CUSTOM_REAL
    write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                          //'/proc',myrank,trim(REG)//trim(kernel_name)//'.bin'

    open(IIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ier)
    if (ier /= 0) then
      write(*,*) '  kernel not found: ',trim(k_file)
      stop 'Error kernel file not found'
    endif
    read(IIN) kernel
    close(IIN)

    ! outputs norm of kernel
    norm = sum( kernel * kernel )
    call sum_all_dp(norm, norm_sum)
    if (myrank == 0) then
      print *,'  norm kernel: ',sqrt(norm_sum)
      print *
    endif

    ! source mask
    if (USE_SOURCE_MASK) then
      ! reads in mask
      write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                            //'/proc',myrank,trim(REG)//'mask_source.bin'
      open(IIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0) then
        write(*,*) '  file not found: ',trim(k_file)
        stop 'Error source mask file not found'
      endif
      read(IIN) mask_source
      close(IIN)

      ! masks source elements
      kernel = kernel * mask_source
    endif

    ! sums all kernels from each event
    total_kernel = total_kernel + kernel
  enddo

  ! stores summed kernels
  if (myrank==0) write(*,*) 'writing out summed kernel for: ',trim(kernel_name)

  write(k_file,'(a,i6.6,a)') 'OUTPUT_SUM/proc',myrank,trim(REG)//trim(kernel_name)//'.bin'

  open(IOUT,file=trim(k_file),form='unformatted',status='unknown',action='write',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'Error kernel not written: ',trim(k_file)
    write(*,*) 'Please check if directory OUTPUT_SUM/ exists...'
    stop 'Error kernel write'
  endif
  write(IOUT) total_kernel
  close(IOUT)

  if (myrank==0) write(*,*)

  ! frees memory
  deallocate(kernel,total_kernel)
  if (USE_SOURCE_MASK) deallocate(mask_source)

end subroutine sum_kernel

