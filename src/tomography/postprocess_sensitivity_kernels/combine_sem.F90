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

  use postprocess_par, only: MAX_STRING_LEN, MAX_KERNEL_PATHS, MAX_KERNEL_NAMES, &
    CUSTOM_REAL, NGLLX, NGLLZ, IIN, IOUT

  implicit none

  character(len=MAX_STRING_LEN) :: kernel_paths(MAX_KERNEL_PATHS), kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: line,filename,output_dir,input_file
  character(len=MAX_STRING_LEN) :: arg(3)
  integer :: npath,nker,nspec
  integer :: i,ier,iker
  integer :: filesize
  ! mpi
  integer :: myrank, NPROC

  ! MPI initialization
  call init_mpi()
  call world_size(NPROC)
  call world_rank(myrank)

  if (myrank == 0) then
    print *, 'Running XCOMBINE_SEM'
    print *

    if (command_argument_count() /= 3) then
      print *, 'mpirun -n NPROC bin/xcombine_sem KERNEL_NAMES INPUT_FILE OUTPUT_DIR'
      print *
      call stop_the_code('Please check command line arguments')
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
     call stop_the_code('Please check command line argument: INPUT_FILE')
  endif
  do while (.true.)
     read(IIN,'(a)',iostat=ier) line
     if (ier /= 0) exit
     npath = npath+1
     if (npath > MAX_KERNEL_PATHS) call stop_the_code('Error number of paths exceeds MAX_KERNEL_PATHS')
     kernel_paths(npath) = line
  enddo
  close(IIN)

  ! Attempt to determine NSPEC directly from Fortran binary file.
  ! Advantage of this approach is that the utility doesn't have to be recompiled
  ! whenever mesh changes, and avoids dealing with SPECFEM2D database system,
  ! which is a bit messy. Disadvantage of this approach is that it is a hack and
  ! possibly not portable.

  write(filename, '(a,i6.6,a)') '/proc',myrank,'_'//trim(kernel_names(1))//'.bin'
  open(IIN, file=trim(kernel_paths(1))//trim(filename))
  inquire(IIN, size=filesize)
  close(IIN)

  nspec = (filesize-8)/(CUSTOM_REAL*NGLLX*NGLLZ)

  do iker= 1,nker
      call combine_sem_array(kernel_names(iker),kernel_paths,output_dir,npath,nspec,myrank)
  enddo

  ! MPI finish
  call finalize_mpi()

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
  if (ier /= 0) call stop_the_code('Error allocating array')

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
      call stop_the_code('Error array file not found')
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
    call stop_the_code('Error array write')
  endif
  write(IOUT) sum_arrays
  close(IOUT)

  write(*,*)
  deallocate(array,sum_arrays)

  end subroutine combine_sem_array


