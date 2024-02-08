program filter_input_trace

! This program is used to filter signals with a Butterworth bandpass filter.
! It has been initially written by Vadim Monteiller and has been adapted by Alexis Bottero.
! To compile: ifort filter.f90 -o filter
!
! usage: filter [-h] [-v] [--norder FILTER_ORDER] [--output PATH_TO_FILTERED_SIGNAL]
!               PATH_TO_SIGNALS -f1 F1 -f2 F2
!
! Positional arguments :
!   PATH_TO_SIGNALS                  Paths to signal to filter (two columns file time,amplitude)
!   -f1 F1                           Inferior cutoff frequency
!   -f2 F2                           Superior cutoff frequency
!
! Optional arguments :
!   -h, --help                       Show this help message and exit
!   --norder FILTER_ORDER            Order of the filter (default 4)
!   --output PATH_TO_FILTERED_SIGNAL If one signal given: path where we store the filtered
!                                    signal (default PATH_TO_SIGNAL+Filt)
!   -v, --verbose                    Verbose mode

  implicit none

  !======= Variables for argument parsing ======!
  integer :: narg,cptArg ! number of arguments given by user & counter of arguments
  character(len=100) :: name,f1str,f2str,norderstr ! To store argument as character strings
  logical :: lookForF1 = .false.
  logical :: lookForF2 = .false.
  logical :: lookForNorder = .false.
  logical :: lookForOutput = .false.
  logical :: fileExist
  logical :: verbose = .false.

  !======= Local variables ======!
  character(len=100), dimension(:), allocatable :: filesToFilter,output_files
  real(kind=4), dimension(:), allocatable :: input_tr1,input_tr2, output_tr1,output_tr2

  !=============== Default parameters ===============!
  character(len=100) :: output_file_name = "default"
  real(kind=4) :: f1 = -9999
  real(kind=4) :: f2 = -9999
  real(kind=4) :: dt = -9999
  integer      :: norder = 4
  integer      :: nFiles = 0

  !======= Local variables ======!
  integer      :: nstep ! To store the number of lines in file
  integer      :: i = 1
  integer      :: ifile = 1
  integer      :: irek
  real(kind=4) :: x

  !=============== Help message ===============!
  200 format(/5x,'This program is used to filter signals with a Butterworth bandpass filter.',/5x, &
       'It has been initially written by Vadim Monteiller and has been adapted by Alexis Bottero.',/5x, &
       '',/5x, &
       'usage: filter [-h] [-v] [--norder FILTER_ORDER] [--output PATH_TO_FILTERED_SIGNAL]',/5x, &
       '              PATH_TO_SIGNALS -f1 F1 -f2 F2',/5x, &
       '',/5x, &
       'Positional arguments :',/5x, &
       '  PATH_TO_SIGNALS                  Paths to signal to filter (two columns file time,amplitude)',/5x, &
       '  -f1 F1                           Inferior cutoff frequency',/5x, &
       '  -f2 F2                           Superior cutoff frequency',/5x, &
       '',/5x, &
       'Optional arguments :',/5x, &
       '  -h, --help                       Show this help message and exit',/5x, &
       '  --norder FILTER_ORDER            Order of the filter (default 4)',/5x, &
       '  --output PATH_TO_FILTERED_SIGNAL If one signal given: path where we store the filtered signal',/5x, &
       '                                   (default PATH_TO_SIGNAL+Filt)',/5x, &
       '  -v, --verbose                    Verbose mode')

  !=============== Parse arguments ===============!
  ! Check if any arguments are found
  narg = command_argument_count()
  ! Loop over the arguments
  if (narg > 0) then
    ! Loop across options
    do cptArg = 1,narg
      call get_command_argument(cptArg,name)
      select case(adjustl(name))
        case("--help","-h")
          write(*,200)
          stop
        case("-f1")
          lookForF1 = .true.
        case("-f2")
          lookForF2 = .true.
        case("--norder")
          lookForNorder = .true.
        case("--output")
          lookForOutput = .true.
        case("--verbose","-v")
          verbose = .true.
        case default
          if (lookForF1) then
            f1str = adjustl(name)
            read(f1str,*) f1
            lookForF1 = .false.
          else if (lookForF2) then
            f2str = adjustl(name)
            read(f2str,*) f2
            lookForF2 = .false.
          else if (lookForNorder) then
            norderstr = adjustl(name)
            read(norderstr,*) norder
            lookForNorder = .false.
          else if (lookForOutput) then
            output_file_name = adjustl(name)
            lookForOutput = .false.
          else
            nFiles = nFiles + 1
          endif
      end select
    enddo

    if (nFiles > 0) then
      ! Allocate the array to store file names
      allocate(filesToFilter(nFiles),output_files(nFiles))
    else
      write(*,*) "At least one file must be given... Nothing has been done!"
      write(*,200)
      stop
    endif

    ! Second loop across options to store files
    do cptArg = 1,narg
      call get_command_argument(cptArg,name)
      select case(adjustl(name))
        case("-f1")
          lookForF1 = .true.
        case("-f2")
          lookForF2 = .true.
        case("--norder")
          lookForNorder = .true.
        case("--output")
          lookForOutput = .true.
        case("--verbose","-v")
          verbose = .true.
        case default
          if (lookForF1) then
            lookForF1 = .false.
          else if (lookForF2) then
            lookForF2 = .false.
          else if (lookForNorder) then
            lookForNorder = .false.
          else if (lookForOutput) then
            lookForOutput = .false.
          else
            filesToFilter(i) = adjustl(trim(name)) ! assign a value to filesToFilter(i)
            inquire(file=filesToFilter(i),exist=fileExist) ! check if the file exists
            if (.not. fileExist) then
              write(*,*)'File ',adjustl(trim(filesToFilter(i))),' not found.'
              write(*,*) "Nothing has been done!"
              stop
            endif
            i = i + 1
          endif
      end select
    enddo
  endif

  if (f1 == -9999 .or. f2 == -9999) then
    write(*,*) "Options -f1 and -f2 must be given... Nothing has been done!"
    write(*,200)
    stop
  endif

  if (output_file_name == "default") then
    do ifile = 1,nFiles
      output_files(ifile) = adjustl(trim(filesToFilter(ifile))) // "Filt"
    enddo
  else if (nFiles == 1) then
    output_files(1) = adjustl(trim(output_file_name))
  else
    write(*,*) "Option --output can only be used if just one file is suppose to be filtered"
    write(*,*) "Nothing has been done!"
    stop
  endif

  if (verbose) then
    if (nFiles == 1) then
      write(*,"(A,A)") "Name of the file containing the trace to filter: ",adjustl(trim(filesToFilter(1)))
      write(*,"(A,A)") "The filtered signal will be written in : ",adjustl(trim(output_files(1)))
    else
      write(*,"(A)") "Name of the files containing the trace to filter:"
      do ifile = 1,nFiles
        write(*,"(A,A)") "  ",adjustl(trim(filesToFilter(ifile)))
      enddo
      write(*,*)
      write(*,"(A)") "The filtered signals will be written in :"
      do ifile = 1,nFiles
        write(*,"(A,A)") "  ",adjustl(trim(output_files(ifile)))
      enddo
    endif
    write(*,*)
    write(*,"(A,F10.2,A)") "f1 : ",f1," Hz"
    write(*,"(A,F10.2,A)") "f2 : ",f2," Hz"
    write(*,"(A,I3)") "Filter order : ",norder
    write(*,*)
  endif

  do ifile = 1,nFiles
    ! Open file containing signal to filter
    open(10,file=filesToFilter(ifile))

    ! Count the number of values in file
    nstep = 0
    do
      read(10,*,end=99) x
      nstep = nstep + 1
    enddo
99  close(10)

    allocate(input_tr1(nstep),output_tr1(nstep),input_tr2(nstep),output_tr2(nstep))
    output_tr2(:) = 0.

    ! Read input trace
    open(10,file=filesToFilter(ifile))
    do i = 1,nstep
      read(10,*) input_tr1(i),input_tr2(i)
    enddo
    close(10)

    dt = input_tr1(2) - input_tr1(1)
    !if (verbose) then
    !  write(*,"(A,F5.2)") "dt :",dt
    !endif
    output_tr1(:) = input_tr1(:)

    irek = 1
    call bwfilt (input_tr2, output_tr2, dt, nstep, irek, norder, f1, f2)

    open(10,file=output_files(ifile))

    do i = 1,nstep
      write(10,*) output_tr1(i),output_tr2(i)
    enddo
    close(10)

    write(*,"(A,A,A)") "File ",adjustl(trim(output_files(ifile)))," has been written."
    deallocate(input_tr1,output_tr1,input_tr2,output_tr2)
  enddo

  deallocate(filesToFilter,output_files)

end program filter_input_trace


!#############################################################
subroutine bwfilt (x, y, dt, n, irek, norder, f1, f2)

  ! recursive filtering of data with butterworth filter
  ! x: input array
  ! y: output array
  ! dt: time increment
  ! n: number of data points

  ! irek=0: forward filtering only
  ! irek=1: forward and backward filtering

  ! norder: order of butterworth filter
  ! norder=0: only filtering, no determination of coefficients
  ! norder < 0: no starplots of transfer function and impulse response

  ! f1: low cutoff frequency (Hz)
  ! f1=0: low pass filter

  ! f2: high cutoff frequency (Hz)
  ! f2>0.5/dt: high pass filter

  implicit none

  real(kind=4), dimension(1)::x,y
  real(kind=4), dimension (10) ::  a, b1, b2
  real(kind=4) :: dt,f1,f2
  integer :: iunit, npoles,norder,irek,n,lx
  !real(kind(0d0)) :: x(n),y(n)

   iunit = 3

   if (norder /= 0) then
      npoles=iabs(norder)
      !determination of filter coefficients
      call bpcoeff(f1,f2,npoles, dt, a,b1, b2)
      if (norder >= 0) then
         !plot of transfer function and impuulse response
         lx = 100
         !filtering
      endif
   endif


   if (n /= 0) then
      call rekurs(x,y,n,a,b1,b2,npoles,irek)
   endif
   return
 end subroutine bwfilt

!---------------------------------------------------------------


subroutine rekurs(x,y,ndat,a,b1,b2,npoles,iflag)
  ! performs recursive filtering of data in array x of length ndat
  ! filtered output in y
  ! a, b1, b2 are the filtercoefficients previously determined in bwcoef
  ! npoles is the number of poles
  ! iflag=0: forward filtering only
  ! iflag /= 0: forward and backward filtering

  implicit none

  real(kind=4), dimension(10) :: z,z1,z2 ,a,b1,b2
  real(kind=4) ::  x1,x2
  integer :: ndat, npoles, iflag, n,i
  real(kind=4) :: x(ndat), y(ndat)

  !forward

  x1 = 0.d0
  x2 = 0.d0

  do i = 1, npoles
     z1(i) = 0.d0
     z2(i) = 0.d0
  enddo

  do n = 1, ndat
     z(1) = a(1)*(x(n)-x2) -b1(1)*z1(1) -b2(1)*z2(1)
     do i = 2, npoles
        z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
     enddo
     x2=x1
     x1=x(n)
     do i = 1, npoles
        z2(i) =z1(i)
        z1(i) =z(i)
     enddo
     y(n) = z(npoles)
  enddo

  if (iflag == 0) then
     return
  endif

  !backward

  x1 =0.d0
  x2 =0.d0

  do i = 1, npoles
     z1(i) = 0.d0
     z2(i) = 0.d0
  enddo

  do n = ndat, 1, -1
     z(1) = a(1)*(y(n)-x2)-b1(1)*z1(1)-b2(1)*z2(1)
     do i =2, npoles
        z(i) = a(i)*(z(i-1)-z2(i-1))-b1(i)*z1(i)-b2(i)*z2(i)
     enddo
     x2=x1
     x1=y(n)
     do i = 1,npoles
        z2(i)=z1(i)
        z1(i)=z(i)
     enddo
     y(n) = z(npoles)
  enddo
  return
end subroutine rekurs




!---------------------------------------------------------------


subroutine bpcoeff(f1,f2,npoles,dt,a,b1,b2)
  !determines filtercoefficients for recursive bandpassfilter

  real(kind=4),dimension(10) :: a,b1,b2
  complex(kind=4) :: s(20), t1,t2,p
  real(kind=4), parameter :: pi = 3.141592653589793d0
  real(kind=4) :: f1,f2,dt,d2,w0,w1,w2,ssum, sprod,fact1,fact2,fact3
  integer :: i,npol2,n,npoles


  if (npoles > 10) then
     stop ' npoles greater than 10: STOP '
  endif

  d2= 2.d0/dt
  w1=d2*tan(2.d0*pi*f1/d2)
  w2=d2*tan(2.d0*pi*f2/d2)
  w0=0.5*(w2-w1)

  i=1
  npol2=npoles/2+1
  do n =1,npoles
     p = cexp(cmplx(0.d0,dble(2*n-1+npoles)*pi/dble(2*npoles)))
     t1 = p*cmplx(w0,0.d0)
     t2 = sqrt(t1*t1-cmplx(w1*w2,0.d0))
     s(i)=t1+t2
     s(i+1)=t1-t2
     i=i+2
  enddo

  do n=1,npoles
     ssum=2*real(s(n))
     sprod=dble(s(n)*conjg(s(n)))
     fact1=d2*d2-d2*ssum+sprod
     fact2=2.d0*(sprod-d2*d2)
     fact3=d2*d2+d2*ssum+sprod
     a(n)=2.d0*d2*w0/fact1
     b1(n)=fact2/fact1
     b2(n)=fact3/fact1
  enddo
  return
end subroutine bpcoeff


!=====================================================================
