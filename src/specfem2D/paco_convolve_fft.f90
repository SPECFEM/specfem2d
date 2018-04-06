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

  integer,intent(in) :: NSTEP,NFREC,label

  complex(selected_real_kind(15,300)), dimension(NFREC+1),intent(in) :: Field

  double precision, dimension(NSTEP),intent(out) :: output_field
  double precision,intent(in) :: dt,tp,ts

  ! local parameters
  integer :: N,J
  double precision :: AN,FUN,RAIZ
  double precision, external :: RIC, deRIC, de2RIC
  complex(selected_real_kind(15,300)) :: CR(2*NFREC)

  N = 2 * NFREC
  AN = N

!
! label=1, displacement field U as input, convolution by a Ricker for displacement U
! label=2, displacement field U as input, convolution by the first derivative of a Ricker for velocity V
! label=3, displacement field U as input, convolution by the second derivative of a Ricker for acceleration A
! label=4, displacement field T as input, convolution by a Ricker
!
! flag = 0 on a besoin de U, V et A (pas T)
! flag /= 0 on a besoin de T et V (pas U ni A)
!
! NSTEP==1, FLAG==0 (flags: interior=0, left= 1, right= 2, bottom=3)
!

  do j = 1,N
    select case (label)
    case (1,4)
      FUN = ric(j,tp,ts,dt)
    case (2)
      FUN = deric(j,tp,ts,dt)
    case (3)
      FUN = de2ric(j,tp,ts,dt)
    case default
      call stop_the_code('Invalid label in paco_convolve_fft')
    end select

    CR(j) = CMPLX(FUN,0.0d0)
  enddo

  call fourier_transform(N,CR,-1.0d0)

  RAIZ = SQRT(AN)

  call SINTER(Field,output_field,NSTEP,CR,RAIZ,NFREC,label,dt)

  end subroutine paco_convolve_fft

!
!-------------------------------------------------------------------------------------
!

  subroutine SINTER(V,output_field,NSTEP,CR,RAIZ,NFREC,label,dt)

  implicit none

  integer,intent(in) :: NSTEP,NFREC
  integer,intent(in) :: label

  complex(selected_real_kind(15,300)), dimension(NFREC+1),intent(in) :: V
  double precision, dimension(NSTEP),intent(out) :: output_field

  complex(selected_real_kind(15,300)),intent(in) :: CR(2*NFREC)
  double precision,intent(in) :: RAIZ
  double precision,intent(in) :: dt

  ! local parameters
  integer :: j,jn,N,mult,delay
  double precision :: VT(2*NFREC)
  double precision :: filt
  complex(selected_real_kind(15,300)) :: VC
  complex(selected_real_kind(15,300)) :: CY(2*NFREC)

  N = 2 * NFREC

  CY(1) = CR(1) * V(1) * RAIZ * dt

  do J = 2,N/2+1
     FILT = 1.0d0
     VC = V(J)
     CY(J)= CR(J)*VC * RAIZ * dt/ FILT
     JN = N-J+2
     CY(JN) = CONJG(CY(J))
  enddo

  call fourier_transform(N,CY,1.0d0)

  if (label == 1 .or. label == 3 .or. (label == 2 .and. NSTEP == 1)) then
! coefficients to take time steps needed (t=0: first time step)
     mult = 1
     delay = 0
  else if (label == 2 .and. NSTEP > 1) then
! coefficients to take time steps needed (t=i*deltat+1/2: one step on two starting at 1/2)
     mult = 2
     delay = 0
  else if (label == 4) then
! coefficients to take time steps needed (t=i*deltat+1: one step on two starting at 1)
     mult = 2
     delay = 1
  endif

  do J = 1,NSTEP
     CY(mult*J+delay) = CY(mult*J+delay)/RAIZ/dt
     VT(mult*J+delay) = real(CY(mult*J+delay))
     output_field(J) = VT(mult*J+delay)
  enddo

  end subroutine SINTER

!
!-------------------------------------------------------------------------------------
!

!
! Ricker time function
!
  function RIC(J,tp,ts,dt)

  use constants, only: PI

  implicit none

  integer,intent(in) :: J
  double precision,intent(in) :: tp,ts,dt
  double precision :: RIC

  ! local parameters
  double precision :: A

  A = PI*(dt*(J-1)-ts)/tp
  A = A*A

  RIC = 0.0d0
  if (A > 30.0d0) return

  RIC = (A-0.5)*EXP(-A)

  end function RIC

!
!-------------------------------------------------------------------------------------
!

!
! first time derivative of Ricker time function
!
  function deRIC(J,tp,ts,dt)

  use constants, only: PI

  implicit none

  integer,intent(in) :: J
  double precision,intent(in) :: tp,ts,dt
  double precision :: deRIC

  ! local parameters
  double precision :: A,A_dot

  A = PI*(dt*(J-1)-ts)/tp
  A = A*A
  A_dot = 2*(PI/tp)**2*(dt*(J-1)-ts)

  deRIC = 0.0d0
  if (A > 30.0d0) return

  deRIC = A_dot*(1.5-A)*EXP(-A)

  end function deRIC

!
!-------------------------------------------------------------------------------------
!

!
! second time derivative of Ricker time function
!
  function de2RIC(J,tp,ts,dt)

  use constants, only: PI

  implicit none

  integer,intent(in) :: J
  double precision,intent(in) :: tp,ts,dt
  double precision :: de2RIC

  ! local parameters
  double precision :: A,A_dot,A_dot_dot

  A = PI*(dt*(J-1)-ts)/tp
  A = A*A
  A_dot = 2*(PI/tp)**2*(dt*(J-1)-ts)
  A_dot_dot = 2*(PI/tp)**2

  de2RIC = 0.0d0
  if (A > 30.0d0) return

  de2RIC = (A_dot_dot*(1.5-A)-A_dot*A_dot-A_dot*(1.5-A)*A_dot)*EXP(-A)

  end function de2RIC

!
!-------------------------------------------------------------------------------------
!

! Fourier transform
  subroutine fourier_transform(LX,CX,SIGNI)

! note: this routine is based on: http://pages.mtu.edu/~jdiehl/Potential_Fields/fork.f
!
! calculates the Fourier transform of a one-dimensional array.
! Algorithm from Claerbout, J.F.,
! "Fundamentals of Geophysical Data Processing with Applications to Petroleum Prospecting", McGraw-Hill, 1976.
!
! Input/output parameters:
!    Complex array cx of length lx is the input array.  Upon
!    return, cx contains the transformed array.  Length of
!    array must be a power of 2.
!    If signi = -1., then the forward calculation is performed;
!    if signi =  1., the inverse transform is performed.

  use constants, only: PI

  implicit none

  integer,intent(in) :: LX
  complex(selected_real_kind(15,300)),intent(inout) :: CX(LX)
  double precision,intent(in) :: SIGNI

  ! local parameters
  integer :: i,j,l,istep,m
  double precision :: SC
  complex(selected_real_kind(15,300)) :: CARG,CW,CTEMP

  J = 1
  SC = SQRT(1.0d0/LX)
  do I = 1,LX
    if (I <= J) then
      CTEMP = CX(J)*SC
      CX(J) = CX(I)*SC
      CX(I) = CTEMP
    endif
    M = LX/2
    do while (M >= 1 .and. M < J)
      J = J-M
      M = M/2
    enddo
    J = J+M
  enddo
  L = 1

  do while(L < LX)
    ISTEP = 2*L
    do  M = 1,L
      CARG = (0.0d0,1.0d0)*(PI*SIGNI*(M-1))/L
      CW = EXP(CARG)
      do  I = M,LX,ISTEP
        CTEMP = CW*CX(I+L)
        CX(I+L) = CX(I)-CTEMP
        CX(I) = CX(I)+CTEMP
      enddo
    enddo
    L = ISTEP
  enddo

  end subroutine fourier_transform

