!
! simple program to create a time-windowed adjoint source trace
!
  program create_adjoint_source

  implicit none

  double precision :: time,displ_potential
  double precision :: displ_max

  integer :: nstep
  integer :: i,ier
  character(len=150) :: filename,filename_adj

  ! time window (starts at initial time)
  double precision,parameter :: TIME_WINDOW_END = 1.274

  ! default
  filename = 'OUTPUT_FILES/AA.S0001.POT.semx'
  filename_adj = 'SEM/AA.S0001.POT.adj'

  ! note: the forward simulation outputs the fluid potential values as trace (potential_acoustic),
  !       which will be used here directly as adjoint source
  !
  !       in general, for acoustic domains pressure would be defined as
  !         p = - \partial_t^2 chi
  !           = - potential_dot_dot_acoustic
  !
  !       and for adjoint simulations, the adjoint source in terms of pressure differences would be
  !         f^adj = - \partial_t^2 (p^syn - p^obs)
  !               = - \partial_t^2 (-potential_dot_dot_acoustic^syn + potential_dot_dot_acoustic^obs)
  !               = \partial_t^2 (potential_dot_dot_acoustic^syn - potential_dot_dot_acoustic^obs)
  !
  !       for simplicity, here a different adjoint source, simply defined by
  !         f^adj = chi
  !               = potential_acoustic
  !       will be used.

  ! user output
  print *,'creating adjoint source:'
  print *,'  input trace file           = ',trim(filename)
  print *,'  output adjoint source file = ',trim(filename_adj)

  ! open trace
  open(1,file=trim(filename),status="old",action="read",iostat=ier)
  if (ier /= 0) then
    print *,'Error could not open trace file for reading: ',trim(filename)
    print *,'Please check if file exists...'
    stop 'Error opening trace file'
  endif

  ! adjoint source file
  open(10,file=trim(filename_adj),status="unknown",action="write",iostat=ier)
  if (ier /= 0) then
    print *,'Error could not open adjoint source file for writing out: ',trim(filename_adj)
    print *,'Please check if directory SEM/ exists...'
    stop 'Error opening adjoint source file'
  endif

  ! gets number of time steps in trace file
  nstep = 0
  ier = 0
  do while (ier == 0)
    read(1,*,iostat=ier) time,displ_potential
    if (ier == 0) nstep = nstep + 1
  enddo
  rewind(1)

  print *
  print *,'number of time steps in trace file: ',nstep
  print *

  ! statistics
  displ_max = 0.d0

  do i = 1,nstep
    ! reads in trace values
    read(1,*) time,displ_potential

    ! zeroes out source after time window end
    if (time > TIME_WINDOW_END) displ_potential = 0.d0

    ! writes out as adjoint source
    write(10,*) time,displ_potential

    ! statistics
    if (abs(displ_potential) > displ_max) displ_max = abs(displ_potential)

  enddo

  ! closes files
  close(1)
  close(10)

  ! user output
  print *
  print *,'adjoint source: maximum amplitude = ',displ_max
  print *
  print *,'done, see adjoint source file: ',trim(filename_adj)


  end program
