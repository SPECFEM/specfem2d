!
! This subroutine was written by Paco Sanchez-Sesma and his colleagues
! from the Autonomous University of Mexico (UNAM), Mexico City, Mexico
!
!     PROGRAMA PARA CALCULAR SISMOGRAMAS SINTETICOS DADA LA
!     FUNCION DE TRANSFERENCIA PARA COMPONENTES Ux, Uz, R2
!     Tx y Tz  SOLUCION DE CAMPO LIBRE   Caso P-SV, RAYLEIGH
!
! modified by Dimitri Komatitsch and Ronan Madec in March 2008
! in particular, converted to Fortran90 and to double precision

subroutine paco_convolve_fft(Field,label,NSTEP,dt,NFREC,output_field,tp,ts)

  implicit none

  integer :: NFREC,N,NSTEP

  complex(selected_real_kind(15,300)), dimension(NFREC+1) :: Field

  complex(selected_real_kind(15,300)) :: CR(2*NFREC)

  double precision, dimension(NSTEP) :: output_field

  integer :: J,label

  double precision :: AN,FUN,RAIZ,dt,tp,ts

  double precision, external :: RIC, deRIC, de2RIC

  N=2*NFREC

  AN  = N

!
! label=1 <=> champ U en entree =>convolution par un ricker pour U
! label=2 <=> champ U en entree =>convolution par la derivee de ricker pour V
! label=3 <=> champ U en entree =>convolution par la derivee seconde de ricker pour A
! label=4 <=> champ T en entree =>convolution par un ricker
!
! flag=0 on a besoin de U, V et A (pas T)
! flag/=0 on a besoin de T et V (pas U ni A)
!
! NSTEP==1 <=> FLAG==0 (flags: interior=0, left=1, right=2, bottom=3)
!

  do j=1,N
     if (label==1 .or. label==4) FUN=ric(j,tp,ts,dt)
     if (label==2) FUN=deric(j,tp,ts,dt)
     if (label==3) FUN=de2ric(j,tp,ts,dt)
     CR(j)=CMPLX(FUN,0.0d0)
  enddo

  CALL fourier_transform(N,CR,-1.0d0)

  RAIZ = SQRT(AN)

  CALL SINTER(Field,output_field,NSTEP,CR,RAIZ,NFREC,label,dt)

END subroutine paco_convolve_fft

SUBROUTINE SINTER(V,output_field,NSTEP,CR,RAIZ,NFREC,label,dt)

  implicit none

  integer NSTEP, j,jn,N,label,nfrec,mult,delay

  double precision :: RAIZ

  complex(selected_real_kind(15,300)) :: VC

  double precision VT(2*NFREC)

  double precision :: filt,dt

  double precision, dimension(NSTEP) :: output_field

  complex(selected_real_kind(15,300)), dimension(NFREC+1) :: V

  complex(selected_real_kind(15,300)) :: CY(2*NFREC),CR(2*NFREC)

  N=2*NFREC

  CY(1) = CR(1) * V(1) * RAIZ * dt

  DO J=2,N/2+1
     FILT = 1.0d0
     VC   = V(J)
     CY(J)= CR(J)*VC * RAIZ * dt/ FILT
     JN = N-J+2
     CY(JN)=CONJG(CY(J))
  enddo

  CALL fourier_transform(N,CY,1.0d0)

  if (label==1 .or. label==3 .or. (label==2 .and. NSTEP==1)) then
! coefficients to take time steps needed (t=0: first time step)
     mult=1
     delay=0
  else if(label==2 .and. NSTEP>1) then
! coefficients to take time steps needed (t=i*deltat+1/2: one step on two starting at 1/2)
     mult=2
     delay=0
  else if(label==4) then
! coefficients to take time steps needed (t=i*deltat+1: one step on two starting at 1)
     mult=2
     delay=1
  endif

  do J=1,NSTEP
     CY(mult*J+delay)=CY(mult*J+delay)/RAIZ/dt
     VT(mult*J+delay)=REAL(CY(mult*J+delay))
     output_field(J)=VT(mult*J+delay)
  enddo

END SUBROUTINE SINTER

!
! Ricker time function
!
FUNCTION RIC(J,tp,ts,dt)

  implicit none

  include "constants.h"

  double precision :: A,RIC,tp,ts,dt

  integer j

  A=PI*(dt*(J-1)-ts)/tp
  A=A*A
  RIC=0.0d0
  IF(A>30.0d0) RETURN
  RIC=(A-0.5)*EXP(-A)

END FUNCTION RIC

!
! first time derivative of Ricker time function
!
FUNCTION deRIC(J,tp,ts,dt)

  implicit none

  include "constants.h"

  double precision :: A,A_dot,deRIC,tp,ts,dt
  integer :: j

  A=PI*(dt*(J-1)-ts)/tp
  A=A*A
  A_dot=2*(PI/tp)**2*(dt*(J-1)-ts)
  deRIC=0.0d0
  IF(A>30.0d0) RETURN
  deRIC=A_dot*(1.5-A)*EXP(-A)

END FUNCTION deRIC

!
! second time derivative of Ricker time function
!
FUNCTION de2RIC(J,tp,ts,dt)

  implicit none

  include "constants.h"

  double precision :: A,A_dot,A_dot_dot,de2RIC,tp,ts,dt
  integer j

  A=PI*(dt*(J-1)-ts)/tp
  A=A*A
  A_dot=2*(PI/tp)**2*(dt*(J-1)-ts)
  A_dot_dot=2*(PI/tp)**2
  de2RIC=0.0d0
  IF(A>30.0d0) RETURN
  de2RIC=(A_dot_dot*(1.5-A)-A_dot*A_dot-A_dot*(1.5-A)*A_dot)*EXP(-A)

END FUNCTION de2RIC


! Fourier transform
SUBROUTINE fourier_transform(LX,CX,SIGNI)

  implicit none

  include "constants.h"

  integer LX,i,j,l,istep,m

  double precision SC

  complex(selected_real_kind(15,300)) :: CX(LX),CARG,CW,CTEMP

  double precision SIGNI

  J=1
  SC=SQRT(1.0d0/LX)
  DO I=1,LX
     IF (I<=J) then
        CTEMP=CX(J)*SC
        CX(J)=CX(I)*SC
        CX(I)=CTEMP
     endif
     M=LX/2
     do while (M>=1 .and. M<J)
        J=J-M
        M=M/2
     enddo
     J=J+M
  enddo
  L=1

  do while(L<LX)
     ISTEP=2*L
     DO  M=1,L
        CARG=(0.0d0,1.0d0)*(PI*SIGNI*(M-1))/L
        CW=EXP(CARG)
        DO  I=M,LX,ISTEP
           CTEMP=CW*CX(I+L)
           CX(I+L)=CX(I)-CTEMP
           CX(I)=CX(I)+CTEMP
        enddo
     enddo

     L=ISTEP
  enddo

END SUBROUTINE fourier_transform

