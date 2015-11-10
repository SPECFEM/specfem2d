
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

! XCOMBINE_SEM
!
! USAGE
!   mpirun -n NPROC ./bin/xcombine_sem KERNEL_NAMES INPUT_FILE OUTPUT_DIR
!
!
! COMMAND LINE ARGUMENTS
!   KERNEL_NAMES           - kernel name, e.g. alpha_kernel
!   INPUT_FILE             - text file containing list of kernel directories
!   OUTPUT_PATH            - directory to which summed kernels are written
!
!
! DESCRIPTION
!   For each name in KERNEL_NAMES, sums kernels from directories specified in
!   INPUT_FILE. Writes the resulting sums to OUTPUT_DIR.
!
!   INPUT_FILE is a text file containing a list of absolute or relative paths to
!   kernel direcotires, one directoy per line.
!
!   This program's primary use case is to clip kernels. It can be used though on
!   any scalar field of dimension (NGLLX,NGLLZ,NSPEC).
!
!   This is a parallel program -- it must be invoked with mpirun or other
!   appropriate utility.  Operations are performed in embarassingly-parallel
!   fashion.


program combine_sem

#ifdef USE_MPI
  use mpi
#endif
  use postprocess_par, only: MAX_STRING_LEN, MAX_KERNEL_PATHS, MAX_KERNEL_NAMES, &
    CUSTOM_REAL, NGLLX, NGLLZ, IIN, IOUT

  implicit none


  character(len=MAX_STRING_LEN) :: kernel_paths(MAX_KERNEL_PATHS), kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: line,filename,output_dir,input_file
  character(len=MAX_STRING_LEN) :: arg(3)
  integer :: npath,nker,nspec
  integer :: i,ier,iker, myrank, nproc
  integer :: filesize

#ifdef USE_MPI
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)
#else
  nproc = 1
  myrank = 0
#endif
  if( ier /= 0 ) stop 'error MPI initialization'


  if (myrank==0) then
  print *, 'Running XCOMBINE_SEM'
  print *

  if (command_argument_count() /= 3) then
    print *, 'mpirun -n NPROC bin/xcombine_sem KERNEL_NAMES INPUT_FILE OUTPUT_DIR'
    print *, ''
    stop 'Please check command line arguments'
  endif
  endif
  do i = 1, 3
    call get_command_argument(i,arg(i), status=ier)
  enddo

  read(arg(1),'(a)') kernel_names_comma_delimited
  read(arg(2),'(a)') input_file
  read(arg(3),'(a)') output_dir

  ! parse names from KERNEL_NAMES
  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)

  ! parse paths from INPUT_FILE
  npath=0
  open(unit = IIN, file = trim(input_file), status = 'old',iostat = ier)
  if (ier /= 0) then
     print *,'Error opening ',trim(input_file)
     stop 'Please check command line argument: INPUT_FILE'
  endif
  do while (.true.)
     read(IIN,'(a)',iostat=ier) line
     if (ier /= 0) exit
     npath = npath+1
     if (npath > MAX_KERNEL_PATHS) stop 'Error number of paths exceeds MAX_KERNEL_PATHS'
     kernel_paths(npath) = line
  enddo
  close(IIN)

  ! Attempt to determine NSPEC directly from Fortran binary file.
  ! Advantage of this approach is that the utility doesn't have to be recompiled
  ! whenever mech changes, and avoids dealing with SPECFEM2D database system,
  ! which is a bit messy. Disadvantage of this approach is that it is a hack and
  ! possibly not portable.

  write(filename, '(a,i6.6,a)') '/proc',myrank,'_'//trim(kernel_names(1))//'.bin'
  open(IIN, file=trim(kernel_paths(1))//trim(filename))
  inquire(IIN, size=filesize)
  close(IIN)
  nspec=(filesize-8)/(CUSTOM_REAL*NGLLX*NGLLZ)

  do iker=1,nker
      call combine_sem_array(kernel_names(iker),kernel_paths,output_dir,npath,nspec,myrank)
  enddo

#ifdef USE_MPI
 call MPI_FINALIZE(ier)
#endif

end program combine_sem

!
!-------------------------------------------------------------------------------------------------
!

subroutine combine_sem_array(kernel_name,kernel_paths,output_dir,npath,nspec,myrank)

  use postprocess_par, only: MAX_STRING_LEN, MAX_KERNEL_PATHS, MAX_KERNEL_NAMES, &
    CUSTOM_REAL, NGLLX, NGLLZ, IIN, IOUT


  implicit none

  character(len=MAX_STRING_LEN) :: kernel_name,kernel_paths(MAX_KERNEL_PATHS)
  character(len=MAX_STRING_LEN) :: output_dir
  integer :: npath,nspec,myrank

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: array,sum_arrays
  integer :: iker,ier

  allocate(array(NGLLX,NGLLZ,NSPEC), &
           sum_arrays(NGLLX,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating array'

 ! loop over kernel paths
  sum_arrays = 0._CUSTOM_REAL
  do iker = 1, npath
  write(*,*) 'reading in array for: ',trim(kernel_name)
  write(*,*) '    ',iker, ' out of ', npath

    ! read array
    array = 0._CUSTOM_REAL
    write(filename,'(a,i6.6,a)') trim(kernel_paths(iker)) //'/proc',myrank,'_'//trim(kernel_name)//'.bin'
    open(IIN,file=trim(filename),status='old',form='unformatted',action='read',iostat=ier)
    if (ier /= 0) then
      write(*,*) '  array not found: ',trim(filename)
      stop 'Error array file not found'
    endif
    read(IIN) array
    close(IIN)

    ! keep track of sum
    sum_arrays = sum_arrays + array

  enddo

  ! write sum
  write(*,*) 'writing sum: ',trim(kernel_name)
  write(filename,'(a,i6.6,a)') trim(output_dir)//'/proc',myrank,'_'//trim(kernel_name)//'.bin'
  open(IOUT,file=trim(filename),form='unformatted',status='unknown',action='write',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'Error array not written:',trim(filename)
    stop 'Error array write'
  endif
  write(IOUT) sum_arrays
  close(IOUT)

  write(*,*)
  deallocate(array,sum_arrays)

end subroutine combine_sem_array


