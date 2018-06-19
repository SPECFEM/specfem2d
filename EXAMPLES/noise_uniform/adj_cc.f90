program adj_cc

  implicit none

  !choose whether to differentiate traces 0, 1, or 2 times
  integer, parameter :: differentiate = 1

  !choose whether to bandpass filter
  logical, parameter :: use_filtering = .true.

  !choose exactly one of the following window options
  logical, parameter :: use_negative_branch = .true.
  logical, parameter :: use_positive_branch = .false.
  logical, parameter :: use_custom_window = .false.

  !choose whether to time reverse, carried out subsequent to all other processing
  logical, parameter :: time_reverse = .false.

  ! FILTERING PARAMETERS
  real  freq_low,freq_high
  data  freq_low  / 2.d-4 /
  data  freq_high / 5.d-1 /

  ! WINDOW PARAMETERS
  real :: t_begin, t_end, w_tukey
  data t_begin / 45.d0  /
  data t_end   / 65.d0 /
  data w_tukey / 0.4    /
  !see explanation below

  ! time variables
  integer :: it, nt, nthalf
  double precision :: dt

  ! data variables
  double precision, dimension(:), allocatable :: seismo_1, seismo_2, seismo_3, seismo_4, &
   seismo_adj, t, w

  ! input/ output
  character(len=64) :: file_in
  integer :: ios

  ! miscellaneous
  double precision, parameter :: PI = 3.141592653589793
  double precision :: norm_adj
  integer :: it_off, it_wdt, it_begin, it_end, k
  integer :: ifreq, nfreq
  real :: F1,F2,D(8),G,DELT
  real :: alpha, beta


  ! EXPLANATION OF WINDOW PARAMETERS

  !To select the desired branch of the cross-correlogram, we employ a Tukey window.  A Tukey taper is just a variant of a cosine taper.

  !W_TUKEY is a number between 0 and 1, 0 being pure boxcar and 1 being pure cosine

  ! Get file info
  call get_command_argument(1,file_in)
  call getlen(file_in,nt)
  call getinc(file_in,nt,dt)
  nthalf = (nt+1)/2

  print *
  write(*,*) 'This routine works only for evenly sampled cross-correlograms.'
  write(*,*) 'Reading from file: '//trim(file_in)
  write(*,'(a,i10)')   ' nt: ', nt
  write(*,'(a,f10.3)') ' dt: ', dt

  ! Allocate, initialize
  allocate(t(nt))
  allocate(w(nt))
  allocate(seismo_1(nt))
  allocate(seismo_2(nt))
  allocate(seismo_3(nt))
  allocate(seismo_4(nt))
  allocate(seismo_adj(nt))
  w(:)          = 0.0d0
  seismo_1(:)   = 0.0d0
  seismo_2(:)   = 0.0d0
  seismo_3(:)   = 0.0d0
  seismo_4(:)   = 0.0d0
  seismo_adj(:) = 0.0d0


  !!!!!!!!!! READ INPUT !!!!!!!!!!!!!!!!!!!!
  open(unit=1001,file=trim(file_in),status='old',action='read')
  do it = 1, nt
      read(1001,*) t(it), seismo_1(it)
  enddo
  close(1001)


  !!!!!!!!!! DIFFERENTIATE !!!!!!!!!!!!!!!!!!!
  seismo_2(1) = 0.0
  seismo_2(nt) = 0.0
  do it = 2, nt-1
      if (differentiate == 0) seismo_2(it) = seismo_1(it)
      if (differentiate == 1) seismo_2(it) = ( seismo_1(it+1) - seismo_1(it-1) ) / (2*dt)
      if (differentiate == 2) seismo_2(it) = ( seismo_1(it+1) - 2 * seismo_1(it) + seismo_1(it-1) ) / dt**2

  enddo


  !!!!!!!!!! FILTER !!!!!!!!!!!!!!!!!!!!
  seismo_3 = seismo_2
  if (use_filtering) then
  ! THIS SECTION CALCULATES THE FILTER AND MUST BE CALLED BEFORE
  ! FILTER IS CALLED
  DELT = 1.0d3 * dt
  F1=freq_low
  F2=freq_high
  call BNDPAS(F1,F2,DELT,D,G,nt)
  !    F1 = LOW FREQUENCY CUTOFF (6 DB DOWN)
  !    F2 = HIGH FREQUENCY CUTOFF (6 DB DOWN)
  !    DELT = SAMPLE INTERVAL IN MILLISECONDS
  !    D = WILL CONTAIN 8 Z DOMAIN COEFICIENTS OF RECURSIVE FILTER
  !    G = WILL CONTAIN THE GAIN OF THE FILTER,
  call FILTER(seismo_3,nt,D,G,2)
  !    X = DATA VECTOR OF LENGTH N CONTAINING DATA TO BE FILTERED
  !    D = FILTER COEFFICIENTS CALCULATED BY BNDPAS
  !    G = FILTER GAIN
  !    IG = 1  one pass
  !    IG = 2  two passes
  endif


  !!!!!!!!!! WINDOW !!!!!!!!!!!!!!!!!!!!
  if (use_custom_window) then
    it_begin = floor((t_begin - t(1))/dt)
    it_end   = ceiling((t_end - t(1))/dt)
    if (it_begin < 1) it_begin = 1
    if (it_end > nt) it_end = nt

  else if (use_positive_branch) then
    write(*,*) 'Choosing positive branch'
    it_begin = nthalf+1
    it_end   = nt
    if (it_begin < nthalf) it_begin = nthalf
    if (it_end > nt) it_end = nt

  else if (use_negative_branch) then
    write(*,*) 'Choosing negative branch'
    it_begin = 1
    it_end   = nthalf
    if (it_begin < 1) it_begin = 1
    if (it_end > nthalf) it_end = nthalf

  else
    write(*,*) 'Must select one of the following: positive_branch, &
                negative_branch, custom_window.'

  endif

  write(*,'(a,2f10.3)') ' Time range: ', t(1), t(nt)
  write(*,'(a,2f10.3)') ' Window:     ', t(it_begin), t(it_end)
  write(*,'(a,f10.3,f10.3)') ' Filtering:  ', 1./freq_high, 1./freq_low

  !! Tukey taper
  alpha = w_tukey
  k=0
  do it = it_begin,it_end
    k=k+1
    beta = real(k-1)/(it_end-it_begin)

    if (beta < alpha/2.) then
      w(it) = 0.5*(1.+cos(2.*pi/alpha*(beta-alpha/2.)))

    else if (beta > alpha/2. .and. beta < 1.-alpha/2.) then
      w(it) = 1.0

    else
      w(it) = 0.5*(1.+cos(2*pi/w_tukey*(beta-1.+alpha/2.)))

    endif
  enddo
  seismo_4 = w * seismo_3


  !!!!!!!!!! NORMALIZE !!!!!!!!!!!!!!!!!!!!
  norm_adj = DOT_PRODUCT(seismo_4,seismo_4) * dt

  if (norm_adj > 1.e-30) then
    seismo_adj(:) = - seismo_4(:) / norm_adj
  else
    ! zero trace
    seismo_adj(:) = 0.0
  endif

  print *
  write(*,*) 'adjoint source norm = ',sngl(norm_adj)

  !!!!!!!!!! WRITE ADJOINT SOURCE !!!!!!!!!!!!!!!!!!!!
  open(unit=1002,file=trim(file_in)//'.adj',status='unknown',iostat=ios)
  if (ios /= 0) write(*,*) 'Error opening output file.'

  print *
  write(*,*) 'Writing to file: '//trim(file_in)//'.adj'

  do it = 1,nt
      if (.not. time_reverse) write(1002,'(f16.12,1pe18.8)'), t(it), seismo_adj(it)
      if (time_reverse) write(1002,'(f16.12,1pe18.8)'), t(it), seismo_adj(nt-it+1)
  enddo
  close(1002)

  write(*,*) 'Finished writing to file.'
  print *


end program adj_cc


!=====================================================================

  subroutine getlen(filename,len)

  implicit none

  !input
  character(len=64) :: filename

  !output
  integer :: len

  !local
  integer, parameter :: IMAX = 1000000
  integer :: i,ios
  real :: dummy1, dummy2

  open(unit=1001,file=trim(filename),status='old',action='read')
  len=0
  do i=1,IMAX
      read(1001,*,iostat=ios) dummy1, dummy2
      if (ios == -1) exit
      len=len+1
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
  integer :: it
  double precision, dimension(len) :: t
  double precision :: sumdt
  real :: dummy

  open(unit=1001,file=trim(filename),status='old',action='read')
  do it=1,len
      read(1001,*) t(it), dummy
  enddo
  close(1001)

  sumdt = 0.0d0
  do it=1,len-1
      sumdt = sumdt + t(it+1) - t(it)
  enddo
  inc=sumdt/(len-1)

  end subroutine getinc


!=====================================================================

  subroutine BNDPAS(F1,F2,DELT,D,G,N)
! RECURSIVE BUTTERWORTH BAND PASS FILTER (KANASEWICH, TIME SERIES
! ANALYSIS IN GEOPHYSICS, UNIVERSITY OF ALBERTA PRESS, 1975; SHANKS,
! JOHN L, RECURSION FILTERS FOR DIGITAL PROCESSING, GEOPHYSICS, V32,
! FILTER.  THE FILTER WILL HAVE 8 POLES IN THE S PLANE AND IS
! APPLIED IN FORWARD AND REVERSE DIRECTIONS SO AS TO HAVE ZERO
! PHASE SHIFT.  THE GAIN AT THE TWO FREQUENCIES SPECIFIED AS
! CUTOFF FREQUENCIES WILL BE -6DB AND THE ROLLOFF WILL BE ABOUT
! THE FILTER TO PREVENT ALIASING PROBLEMS.
    COMPLEX P(4),S(8),Z1,Z2
    real D(8),XC(3),XD(3),XE(3)
    double precision :: X(N)
    DATA ISW/0/,TWOPI/6.2831853/
! THIS SECTION CALCULATES THE FILTER AND MUST BE CALLED BEFORE
! FILTER IS CALLED

!    F1 = LOW FREQUENCY CUTOFF (6 DB DOWN)
!    F2 = HIGH FREQUENCY CUTOFF (6 DB DOWN)
!    DELT = SAMPLE INTERVAL IN MILLISECONDS
!    D = WILL CONTAIN 8 Z DOMAIN COEFICIENTS OF RECURSIVE FILTER
!    G = WILL CONTAIN THE GAIN OF THE FILTER,

      DT=DELT/1000.0
      TDT=2.0/DT
      FDT=4.0/DT
      ISW=1
      P(1)=CMPLX(-.3826834,.9238795)
      P(2)=CMPLX(-.3826834,-.9238795)
      P(3)=CMPLX(-.9238795,.3826834)
      P(4)=CMPLX(-.9238795,-.3826834)
      W1=TWOPI*F1
      W2=TWOPI*F2
      W1=TDT*TAN(W1/TDT)
      W2=TDT*TAN(W2/TDT)
      HWID=(W2-W1)/2.0
      WW=W1*W2
      DO 19 I=1,4
      Z1=P(I)*HWID
      Z2=Z1*Z1-WW
      Z2=CSQRT(Z2)
      S(I)=Z1+Z2
   19 S(I+4)=Z1-Z2
      G=.5/HWID
      G=G*G
      G=G*G
      DO 29 I=1,7,2
      B=-2.0*REAL(S(I))
      Z1=S(I)*S(I+1)
      C=REAL(Z1)
      A=TDT+B+C/TDT
      G=G*A
      D(I)=(C*DT-FDT)/A
   29 D(I+1)=(A-2.0*B)/A
      G=G*G
    5 FORMAT ('-FILTER GAIN IS ', 9E12.6)
      return

      ENTRY FILTER(X,N,D,G,IG)

!     X = DATA VECTOR OF LENGTH N CONTAINING DATA TO BE FILTERED
!     D = FILTER COEFFICIENTS CALCULATED BY BNDPAS
!     G = FILTER GAIN
!     IG = 1  one pass
!     ig = 2  two passes

      if (ISW == 1) goto 31
      write(*,6)
    6 FORMAT ('1BNDPAS MUST BE CALLED BEFORE FILTER')
      return

!     APPLY FILTER IN FORWARD DIRECTION

   31 XM2=X(1)
      XM1=X(2)
      XM=X(3)
      XC(1)=XM2
      XC(2)=XM1-D(1)*XC(1)
      XC(3)=XM-XM2-D(1)*XC(2)-D(2)*XC(1)
      XD(1)=XC(1)
      XD(2)=XC(2)-D(3)*XD(1)
      XD(3)=XC(3)-XC(1)-D(3)*XD(2)-D(4)*XD(1)
      XE(1)=XD(1)
      XE(2)=XD(2)-D(5)*XE(1)
      XE(3)=XD(3)-XD(1)-D(5)*XE(2)-D(6)*XE(1)
      X(1)=XE(1)
      X(2)=XE(2)-D(7)*X(1)
      X(3)=XE(3)-XE(1)-D(7)*X(2)-D(8)*X(1)
      DO 39 I=4,N
      XM2=XM1
      XM1=XM
      XM=X(I)
      K=I-((I-1)/3)*3
      goto (34,35,36),K
   34 M=1
      M1=3
      M2=2
      goto 37
   35 M=2
      M1=1
      M2=3
      goto 37
   36 M=3
      M1=2
      M2=1
   37 XC(M)=XM-XM2-D(1)*XC(M1)-D(2)*XC(M2)
      XD(M)=XC(M)-XC(M2)-D(3)*XD(M1)-D(4)*XD(M2)
      XE(M)=XD(M)-XD(M2)-D(5)*XE(M1)-D(6)*XE(M2)
   39 X(I)=XE(M)-XE(M2)-D(7)*X(I-1)-D(8)*X(I-2)
!
!
      if (ig == 1) goto 3333
      XM2=X(N)
      XM1=X(N-1)
      XM=X(N-2)
      XC(1)=XM2
      XC(2)=XM1-D(1)*XC(1)
      XC(3)=XM-XM2-D(1)*XC(2)-D(2)*XC(1)
      XD(1)=XC(1)
      XD(2)=XC(2)-D(3)*XD(1)
      XD(3)=XC(3)-XC(1)-D(3)*XD(2)-D(4)*XD(1)
      XE(1)=XD(1)
      XE(2)=XD(2)-D(5)*XE(1)
      XE(3)=XD(3)-XD(1)-D(5)*XE(2)-D(6)*XE(1)
      X(N)=XE(1)
      X(N-1)=XE(2)-D(7)*X(1)
      X(N-2)=XE(3)-XE(1)-D(7)*X(2)-D(8)*X(1)
      DO 49 I=4,N
      XM2=XM1
      XM1=XM
      J=N-I+1
      XM=X(J)
      K=I-((I-1)/3)*3
      goto (44,45,46),K
   44 M=1
      M1=3
      M2=2
      goto 47
   45 M=2
      M1=1
      M2=3
      goto 47
   46 M=3
      M1=2
      M2=1
   47 XC(M)=XM-XM2-D(1)*XC(M1)-D(2)*XC(M2)
      XD(M)=XC(M)-XC(M2)-D(3)*XD(M1)-D(4)*XD(M2)
      XE(M)=XD(M)-XD(M2)-D(5)*XE(M1)-D(6)*XE(M2)
   49 X(J)=XE(M)-XE(M2)-D(7)*X(J+1)-D(8)*X(J+2)
 3333 continue
      if (ig == 1) then
        gg=sqrt(g)   ! if only pass once, modify gain
      else
        gg=g
      endif
      DO 59 I=1,N
   59 X(I)=X(I)/gg
      return
  END

