!!!!! this subroutine is used for preparing adjoint sources in this example ONLY
!!!!! generally, you should make your own measuremnts and calculate corresponding adjoint sources
!!!!! FLEXWIN could be of help, or you can modify this code a little bit

program adj_traveltime

  implicit none

  !--------------------------------------------------------------
  ! USER PARAMETERS

  !choose whether to differentiate traces 0, 1, or 2 times
  integer, parameter :: differentiate = 1

  !choose whether to bandpass filter
  logical, parameter :: use_filtering = .true.

  ! cross-correlation branch for measurement (0 = negative / 1 = positive branch / 2 = custom window)
  integer,parameter :: branch_type = 1

  !choose whether to time reverse, carried out subsequent to all other processing
  logical, parameter :: time_reverse = .false.

  ! FILTERING PARAMETERS
  real, parameter :: freq_low  = 2.d-4
  real, parameter :: freq_high = 5.d-1

  ! WINDOW PARAMETERS
  real, parameter :: t_begin = 45.d0
  real, parameter :: t_end   = 65.d0

  ! Tukey window parameter
  real, parameter :: w_tukey = 0.4

  !--------------------------------------------------------------

  !see explanation below

  ! time variables
  integer :: it, nstep, nstep_half
  double precision :: dt

  ! data variables
  double precision, dimension(:), allocatable :: seismo_1, seismo_2, seismo_3, seismo_4, seismo_adj, t, w

  ! input/ output
  character(len=64) :: file_in
  integer :: ios

  ! miscellaneous
  double precision, parameter :: PI = 3.141592653589793d0
  double precision :: norm_adj

  integer :: it_off, it_wdt, it_begin, it_end, k
  integer :: ifreq, nfreq
  integer :: ISW
  real :: F1,F2,D(8),G,DELT
  real :: alpha, beta

  ! EXPLANATION OF WINDOW PARAMETERS
  !
  ! To select the desired branch of the cross-correlogram, we employ a Tukey window.
  ! A Tukey taper is just a variant of a cosine taper.
  !
  ! W_TUKEY is a number between 0 and 1, 0 being pure boxcar and 1 being pure cosine

  ! Get file info
  call get_command_argument(1,file_in)
  call getlen(file_in,nstep)
  call getinc(file_in,nstep,dt)
  nstep_half = (nstep+1)/2

  ! user output
  print *,'noise adjoint source'
  print *
  print *,'This routine works only for evenly sampled cross-correlograms.'
  print *
  if (time_reverse) print *,'using time reversed signal'
  print *,'using filter flag : ',use_filtering,'(false = no filter / true = bandpass)'
  print *,'using Tukey taper : ',w_tukey,'(between 0 = boxcar to 1 = cosine)'
  print *,'using cross-correlation branch : ',branch_type,'(0 = negative / 1 = positive branch / 2 = custom window)'
  print *

  print *,'Reading from file: '//trim(file_in)
  print *,'  nt: ', nstep
  print *,'  dt: ', dt

  ! Allocate, initialize
  allocate(t(nstep))
  allocate(w(nstep))
  allocate(seismo_1(nstep))
  allocate(seismo_2(nstep))
  allocate(seismo_3(nstep))
  allocate(seismo_4(nstep))
  allocate(seismo_adj(nstep))
  w(:)          = 0.0d0
  seismo_1(:)   = 0.0d0
  seismo_2(:)   = 0.0d0
  seismo_3(:)   = 0.0d0
  seismo_4(:)   = 0.0d0
  seismo_adj(:) = 0.0d0

  !!!!!!!!!! READ INPUT !!!!!!!!!!!!!!!!!!!!
  open(unit=1001,file=trim(file_in),status='old',action='read')
  do it = 1,nstep
    read(1001,*) t(it), seismo_1(it)
  enddo
  close(1001)

  !!!!!!!!!! DIFFERENTIATE !!!!!!!!!!!!!!!!!!!
  seismo_2(1) = 0.0
  seismo_2(nstep) = 0.0
  do it = 2,nstep-1
    ! displacement
    if (differentiate == 0) seismo_2(it) = seismo_1(it)
    ! displacement --> particle velocity
    if (differentiate == 1) seismo_2(it) = ( seismo_1(it+1) - seismo_1(it-1) ) / (2.d0*dt)
    ! displacement --> particle acceleration
    if (differentiate == 2) seismo_2(it) = ( seismo_1(it+1) - 2 * seismo_1(it) + seismo_1(it-1) ) / dt**2
  enddo

  !!!!!!!!!! FILTER !!!!!!!!!!!!!!!!!!!!
  seismo_3(:) = seismo_2(:)

  if (use_filtering) then
    ! THIS SECTION CALCULATES THE FILTER AND MUST BE CALLED BEFORE
    ! FILTER IS CALLED
    F1 = freq_low
    F2 = freq_high

    ! time step (in milliseconds)
    DELT = 1.0d3 * dt

    print *,'filtering:'
    print *,'  f_min = ',sngl(F1),'Hz , f_max = ',sngl(F2), 'Hz'
    print *,'  T_min = ',sngl(1.d0/F2),'(s) , T_max = ',sngl(1.d0/F1), '(s)'

    call BNDPAS(F1,F2,DELT,D,G,nstep,ISW)
    !    F1 = LOW FREQUENCY CUTOFF (6 DB DOWN)
    !    F2 = HIGH FREQUENCY CUTOFF (6 DB DOWN)
    !    DELT = SAMPLE INTERVAL IN MILLISECONDS
    !    D = WILL CONTAIN 8 Z DOMAIN COEFICIENTS OF RECURSIVE FILTER
    !    G = WILL CONTAIN THE GAIN OF THE FILTER,
    call FILTER(seismo_3,nstep,D,G,2,ISW)
    !    X = DATA VECTOR OF LENGTH N CONTAINING DATA TO BE FILTERED
    !    D = FILTER COEFFICIENTS CALCULATED BY BNDPAS
    !    G = FILTER GAIN
    !    IG = 1  one pass
    !    IG = 2  two passes
  endif

  !!!!!!!!!! WINDOW !!!!!!!!!!!!!!!!!!!!
  select case(branch_type)
  case (0)
    ! use negative branch
    print *,'Choosing negative branch'
    it_begin = 1
    it_end   = nstep_half
    if (it_begin < 1) it_begin = 1
    if (it_end > nstep_half) it_end = nstep_half

  case (1)
    ! use positive branch
    print *,'Choosing positive branch'
    it_begin = nstep_half+1
    it_end   = nstep
    if (it_begin < nstep_half) it_begin = nstep_half
    if (it_end > nstep) it_end = nstep

  case (2)
    ! custom window
    print *,'Choosing custom window: ',t_begin,'/',t_end,'(s)'
    it_begin = floor((t_begin - t(1))/dt)
    it_end   = ceiling((t_end - t(1))/dt)
    if (it_begin < 1) it_begin = 1
    if (it_end > nstep) it_end = nstep

  case default
    print *,'Must select one of the following: (0 = negative / 1 = positive branch / 2 = custom window)'
    stop 'Invalid branch type selection'

  end select

  write(*,'(a,2f10.3)') ' Time range: ', t(1), t(nstep)
  write(*,'(a,2f10.3)') ' Window:     ', t(it_begin), t(it_end)
  write(*,'(a,f10.3,f10.3)') ' Filtering:  ', 1./freq_high, 1./freq_low

  !! Tukey taper
  alpha = w_tukey
  k = 0
  do it = it_begin,it_end
    k = k+1
    beta = real(k-1)/(it_end-it_begin)

    if (beta < alpha/2.) then
      w(it) = 0.5*(1.+cos(2.*pi/alpha*(beta-alpha/2.)))

    else if (beta > alpha/2. .and. beta < 1.-alpha/2.) then
      w(it) = 1.0

    else
      w(it) = 0.5*(1.+cos(2*pi/w_tukey*(beta-1.+alpha/2.)))

    endif
  enddo
  seismo_4(:) = w(:) * seismo_3(:)

  !!!!!!!!!! NORMALIZE !!!!!!!!!!!!!!!!!!!!
  ! minus sign comes from integration by part
  norm_adj = - DOT_PRODUCT(seismo_4,seismo_4) * dt

  if (abs(norm_adj) > 0.d0) then
    seismo_adj(:) = seismo_4(:) / norm_adj
  else
    ! zero trace
    seismo_adj(:) = 0.d0
  endif

  print *
  print *,'adjoint source norm = ',sngl(abs(norm_adj))

  !!!!!!!!!! WRITE ADJOINT SOURCE !!!!!!!!!!!!!!!!!!!!
  open(unit=1002,file=trim(file_in)//'.adj',status='unknown',iostat=ios)
  if (ios /= 0) then
    print *,'Error opening output file ',trim(file_in)
    stop 'Error opening output file'
  endif
  print *
  print *,'Writing to file: '//trim(file_in)//'.adj'

  do it = 1,nstep
    if (.not. time_reverse) then
      ! non-reversed
      write(1002,*) t(it), seismo_adj(it)
    else
      ! time_reverse
      write(1002,*) t(it), seismo_adj(nstep-it+1)
    endif
  enddo
  close(1002)

  print *,'Finished writing to file.'
  print *

end program adj_traveltime


!=====================================================================

  subroutine getlen(filename,len)

  implicit none

  !input
  character(len=64) :: filename

  !output
  integer :: len

  !local
  integer :: i,ios
  real :: dummy1, dummy2

  open(unit=1001,file=trim(filename),status='old',action='read',iostat=ios)
  if (ios /= 0) then
    print *,'Error reading file ',trim(filename)
    stop 'Error reading input file'
  endif

  len = 0
  do while (ios == 0)
    read(1001,*,iostat=ios) dummy1, dummy2
    if (ios == 0) len = len+1
  enddo
  close(1001)

  end subroutine getlen


!=====================================================================

  subroutine getinc(filename,len,inc)

  implicit none

  !input
  character(len=64) :: filename
  integer :: len

  !output
  double precision :: inc

  !local
  integer :: it,ios
  double precision, dimension(len) :: t
  double precision :: sumdt
  real :: dummy

  open(unit=1001,file=trim(filename),status='old',action='read',iostat=ios)
  if (ios /= 0) then
    print *,'Error reading file ',trim(filename)
    stop 'Error reading input file'
  endif

  do it = 1,len
    read(1001,*) t(it), dummy
  enddo
  close(1001)

  sumdt = 0.0d0
  do it = 1,len-1
      sumdt = sumdt + t(it+1) - t(it)
  enddo
  inc = sumdt/(len-1)

  end subroutine getinc

!
!-----------------------------------------------------------------------------
!
! subroutine for bandpass filter
!
! slighlty modified version from reference in:
!
! E.R. Kanasewich, Time Sequence Analysis in Geophysics, 3rd edition, The University of Alberta Press, 1981
!
! see: (access July 2015)
! https://books.google.com.sa/books?id=k8SSLy-FYagC&pg=PA274&lpg=PA274&dq=bndpas.for&source=bl&ots=gsjFXkN1FZ&sig=W-qvA2kamMr5xIkEIzY_f2yciOI&hl=en&sa=X&redir_esc=y#v=onepage&q=bndpas.for&f=false
!
! or see: website (access ~2012)
! http://www-lgit.obs.ujf-grenoble.fr/users/jrevilla/seiscomp/patch/pack/plugins/seisan/LIB/bndpas.for


  subroutine BNDPAS(F1,F2,DELT,D,G,N,ISW)

! RECURSIVE BUTTERWORTH BAND PASS FILTER (KANASEWICH, TIME SERIES
! ANALYSIS IN GEOPHYSICS, UNIVERSITY OF ALBERTA PRESS, 1975; SHANKS,
! JOHN L, RECURSION FILTERS FOR DIGITAL PROCESSING, GEOPHYSICS, V32,
! FILTER.  THE FILTER WILL HAVE 8 POLES IN THE S PLANE AND IS
! APPLIED IN FORWARD AND REVERSE DIRECTIONS SO AS TO HAVE ZERO
! PHASE SHIFT.  THE GAIN AT THE TWO FREQUENCIES SPECIFIED AS
! CUTOFF FREQUENCIES WILL BE -6DB AND THE ROLLOFF WILL BE ABOUT
! THE FILTER TO PREVENT ALIASING PROBLEMS.

  implicit none

  real, intent(in) :: F1,F2,DELT
  integer, intent(in) :: N

  real, intent(out) :: D(8),G
  integer, intent(out) :: ISW

  ! local parameters
  COMPLEX :: P(4),S(8),Z1,Z2

  real :: DT,TDT,FDT,W1,W2,WW,HWID,A,B,C
  integer :: I

  double precision, parameter :: TWOPI = 2.d0 * 3.141592653589793d0

! THIS SECTION CALCULATES THE FILTER AND MUST BE CALLED BEFORE
! FILTER IS CALLED

!    F1 = LOW FREQUENCY CUTOFF (6 DB DOWN)
!    F2 = HIGH FREQUENCY CUTOFF (6 DB DOWN)
!    DELT = SAMPLE INTERVAL IN MILLISECONDS
!    D = WILL CONTAIN 8 Z DOMAIN COEFICIENTS OF RECURSIVE FILTER
!    G = WILL CONTAIN THE GAIN OF THE FILTER,

  DT = DELT/1000.0
  TDT = 2.0/DT
  FDT = 4.0/DT

  ISW = 1

  P(1) = CMPLX(-.3826834,.9238795)
  P(2) = CMPLX(-.3826834,-.9238795)
  P(3) = CMPLX(-.9238795,.3826834)
  P(4) = CMPLX(-.9238795,-.3826834)

  W1 = TWOPI*F1
  W2 = TWOPI*F2
  W1 = TDT*TAN(W1/TDT)
  W2 = TDT*TAN(W2/TDT)
  HWID = (W2-W1)/2.0
  WW = W1*W2

  do I = 1,4
    Z1 = P(I)*HWID
    Z2 = Z1*Z1-WW
    Z2 = CSQRT(Z2)
    S(I) = Z1+Z2
    S(I+4) = Z1-Z2
  enddo

  G = 0.5/HWID
  G = G*G
  G = G*G

  do I = 1,7,2
    B = -2.0*REAL(S(I))
    Z1 = S(I)*S(I+1)
    C = REAL(Z1)
    A = TDT+B+C/TDT
    G = G*A
    D(I) = (C*DT-FDT)/A
    D(I+1) = (A-2.0*B)/A
  enddo

  G = G*G

  !debug
  !print *,'debug: FILTER GAIN IS ',G

  end subroutine BNDPAS

  subroutine FILTER(X,N,D,G,IG,ISW)

!     X = DATA VECTOR OF LENGTH N CONTAINING DATA TO BE FILTERED
!     D = FILTER COEFFICIENTS CALCULATED BY BNDPAS
!     G = FILTER GAIN
!     IG = 1  one pass
!     IG = 2  two passes

  implicit none
  integer,intent(in) :: N
  double precision,intent(inout) :: X(N)

  real,intent(in) :: D(N),G
  integer,intent(in) :: IG

  integer, intent(in) :: ISW

  ! local parameters
  real :: XC(3),XD(3),XE(3)
  real :: XM,XM1,XM2,GG
  integer :: K,I,J,M,M1,M2

  ! check
  if (ISW /= 1) then
    print *,'1BNDPAS MUST BE CALLED BEFORE FILTER'
    stop 'Invalid ISW in FILTER() routine'
  endif

  ! APPLY FILTER IN FORWARD DIRECTION
  XM2 = X(1)
  XM1 = X(2)
  XM = X(3)

  XC(1) = XM2
  XC(2) = XM1-D(1)*XC(1)
  XC(3) = XM-XM2-D(1)*XC(2)-D(2)*XC(1)

  XD(1) = XC(1)
  XD(2) = XC(2)-D(3)*XD(1)
  XD(3) = XC(3)-XC(1)-D(3)*XD(2)-D(4)*XD(1)

  XE(1) = XD(1)
  XE(2) = XD(2)-D(5)*XE(1)
  XE(3) = XD(3)-XD(1)-D(5)*XE(2)-D(6)*XE(1)

  X(1) = XE(1)
  X(2) = XE(2)-D(7)*X(1)
  X(3) = XE(3)-XE(1)-D(7)*X(2)-D(8)*X(1)

  do I = 4,N
    XM2 = XM1
    XM1 = XM
    XM = X(I)

    K = I-((I-1)/3)*3

    select case(K)
    case (1)
      M = 1
      M1 = 3
      M2 = 2
    case (2)
      M = 2
      M1 = 1
      M2 = 3
    case (3)
      M = 3
      M1 = 2
      M2 = 1
    case default
      stop 'Invalid K value in FILTER'
    end select

    XC(M) = XM-XM2-D(1)*XC(M1)-D(2)*XC(M2)
    XD(M) = XC(M)-XC(M2)-D(3)*XD(M1)-D(4)*XD(M2)
    XE(M) = XD(M)-XD(M2)-D(5)*XE(M1)-D(6)*XE(M2)

    X(I) = XE(M)-XE(M2)-D(7)*X(I-1)-D(8)*X(I-2)
  enddo

  if (IG /= 1) then
    XM2 = X(N)
    XM1 = X(N-1)
    XM = X(N-2)

    XC(1) = XM2
    XC(2) = XM1-D(1)*XC(1)
    XC(3) = XM-XM2-D(1)*XC(2)-D(2)*XC(1)

    XD(1) = XC(1)
    XD(2) = XC(2)-D(3)*XD(1)
    XD(3) = XC(3)-XC(1)-D(3)*XD(2)-D(4)*XD(1)

    XE(1) = XD(1)
    XE(2) = XD(2)-D(5)*XE(1)
    XE(3) = XD(3)-XD(1)-D(5)*XE(2)-D(6)*XE(1)

    X(N) = XE(1)
    X(N-1) = XE(2)-D(7)*X(1)
    X(N-2) = XE(3)-XE(1)-D(7)*X(2)-D(8)*X(1)

    DO I = 4,N
      XM2 = XM1
      XM1 = XM

      J = N-I+1
      XM = X(J)

      K = I-((I-1)/3)*3

      select case (K)
      case (1)
        M = 1
        M1 = 3
        M2 = 2
      case (2)
        M = 2
        M1 = 1
        M2 = 3
      case (3)
        M = 3
        M1 = 2
        M2 = 1
      case default
        stop 'Invalid K value in FILTER'
      end select

      XC(M) = XM-XM2-D(1)*XC(M1)-D(2)*XC(M2)
      XD(M) = XC(M)-XC(M2)-D(3)*XD(M1)-D(4)*XD(M2)
      XE(M) = XD(M)-XD(M2)-D(5)*XE(M1)-D(6)*XE(M2)
      X(J) = XE(M)-XE(M2)-D(7)*X(J+1)-D(8)*X(J+2)
    enddo
  endif

  if (IG == 1) then
    GG = sqrt(G)   ! if only pass once, modify gain
  else
    GG = G
  endif

  do I = 1,N
    X(I) = X(I)/GG
  enddo

  end subroutine FILTER

