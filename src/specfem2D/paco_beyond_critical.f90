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
! original name : DISTRAFF.f
!
!     CALCULO DE DESPLAZAMIENTOS (UX, UZ) y TRACCIONES (TX, TZ) DE CAMPO LIBRE
!     EN UN SEMIESPACIO ELASTICO Y EN LA VECINDAD DE LA SUPERFICIE
!
!     INCIDENCIA DE ONDAS P, SV Y DE RAYLEIGH
!
!     7 de febrero de 2007
!
! modified by Dimitri Komatitsch and Ronan Madec in March 2008
! in particular, converted to Fortran90 and to double precision

  subroutine paco_beyond_critical(anglesource,f0, &
                                  QD,source_type, &
                                  left_bound,right_bound,bot_bound, &
                                  nleft,nright,nbot, &
                                  x_source,cploc,csloc)

  use constants, only: PI,IMAIN

  use specfem_par, only: myrank, &
    coord,nglob,deltat,NSTEP,ATTENUATION_VISCOELASTIC, &
    displ_elastic,veloc_elastic,accel_elastic, &
    v0x_left,v0z_left,v0x_right,v0z_right,v0x_bot,v0z_bot, &
    t0x_left,t0z_left,t0x_right,t0z_right,t0x_bot,t0z_bot

  implicit none

  double precision,intent(in) :: anglesource,f0,QD
  integer,intent(in) :: source_type

  integer,intent(in) :: nleft,nright,nbot
  integer, dimension(nleft),intent(in) :: left_bound
  integer, dimension(nright),intent(in) :: right_bound
  integer, dimension(nbot),intent(in) :: bot_bound

  double precision,intent(in) :: x_source,cploc,csloc

  ! local parameters
  integer :: npt
  integer, dimension(:),allocatable :: local_pt

  double precision, dimension(:), allocatable :: temp_field
  double precision :: dt,TP,delta_in_period

  integer :: J, indice, NSTEP_local, FLAG, N, NFREC, NFREC1
  double precision :: ANU,BEALF,ALFBE,RLM,VNX,VNZ,A1,B1,TOTO,FJ,AKA,AQA,GAMR

  ! location of the point
  double precision :: X,Z,xmin,xmax,zmin,zmax,xmin_glob,xmax_glob,zmin_glob,zmax_glob
  integer :: inode

  double precision :: TS,offset
  complex(selected_real_kind(15,300)) :: CAKA,CAQA,UI,UR
  complex(selected_real_kind(15,300)) :: UX,UZ,SX,SZ,SXZ,A2,B2,AL,AK,AM
  complex(selected_real_kind(15,300)) :: TX,TZ
  complex(selected_real_kind(15,300)), dimension(:),allocatable :: Field_Ux,Field_Uz,Field_Tx,Field_Tz


  ! size of the model
  ! get minimum and maximum values of mesh coordinates
  xmin = minval(coord(1,:))
  xmax = maxval(coord(1,:))
  zmin = minval(coord(2,:))
  zmax = maxval(coord(2,:))
  ! collects min/max on all slices
  call min_all_all_dp(xmin, xmin_glob)
  call max_all_all_dp(xmax, xmax_glob)
  call min_all_all_dp(zmin, zmin_glob)
  call max_all_all_dp(zmax, zmax_glob)
  xmin = xmin_glob
  zmin = zmin_glob
  xmax = xmax_glob
  zmax = zmax_glob

  TS = 1.2d0/f0

! dominant period of the Ricker
  TP = 1.d0/f0

! offset to move the initial location of the source in the horizontal direction of the mesh
  offset = x_source

! find optimal period
! if period is too small, you should see several initial plane wave on your initial field
  delta_in_period = 2.d0
  do while(delta_in_period < 1.5*abs(xmax - xmin)/csloc)
     delta_in_period = 2.d0 * delta_in_period
  enddo

! test Deltat compatibility
  DT = 256.d0
  do while(DT > deltat)
     DT = DT/2.d0
  enddo
  if (abs(DT-deltat) > 1.0d-13) then
    print *, 'Error: Initial plane wave setting has invalid time step size deltat = ',deltat
    print *, 'You must take a deltat that is a power of two (power can be negative)'
    print *, 'For example, you can take ', DT
    call stop_the_code('Error in paco_beyond_critical routine cannot go further, restart with new deltat')
  endif

  DT = deltat/2.d0

  N = 2
  do while(N < 2*NSTEP+1)
     N = 2*N
  enddo

  do while(DT < (delta_in_period/N))
     N = 2*N
  enddo

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'N found to perform the frequency calculation =',N
    write(IMAIN,*) 'number of discrete frequencies               = ',N/2
    write(IMAIN,*) 'delta in period (seconds)                    = ',delta_in_period
    write(IMAIN,*) 'delta in frequency (Hz)                      = ',1.d0/delta_in_period
    write(IMAIN,*) 'dt (here we need deltat/2)                   = ', DT
  endif

  if (mod(N,2) /= 0) call stop_the_code('N must be a multiple of 2')

  NFREC = N/2
  NFREC1 = NFREC + 1


!
!     FDT:  FUNCION DE TRANSFERENCIA
!

! calculation of Poisson's ratio
  ANU = (cploc*cploc-2.d0*csloc*csloc)/(2.d0*(cploc*cploc-csloc*csloc))

  if (myrank == 0) then
    write(IMAIN,*) "Poisson's ratio                              = ",ANU
    call flush_IMAIN()
  endif

  UI = (0.0d0, 1.0d0) ! imaginary unit
  UR = (1.0d0, 0.0d0) ! real unit

! convert angle to radians
  GAMR = anglesource

  BEALF = SQRT((1.0d0-2.0d0*ANU)/(2.0d0*(1.0d0-ANU)))
  ALFBE = 1.0d0/BEALF
  RLM = ALFBE**2 - 2.0d0

! flags: interior=0, left= 1, right= 2, bottom=3
  do FLAG = 0,3

    ! sets points
    select case (FLAG)
    case (0)
      ! all points
      if (myrank == 0) write(IMAIN,*) "calculation of the initial field for every point of the mesh"
      npt = nglob
      allocate(local_pt(npt))
      do inode = 1,npt
         local_pt(inode) = inode
      enddo
      NSTEP_local = 1
    case (1)
      ! left boundary
      if (myrank == 0) write(IMAIN,*) "calculation of every time step on the left absorbing boundary"
      npt = nleft
      allocate(local_pt(npt))
      local_pt = left_bound
      NSTEP_local = NSTEP
    case (2)
      ! right boundary
      if (myrank == 0) write(IMAIN,*) "calculation of every time step on the right absorbing boundary"
      npt = nright
      allocate(local_pt(npt))
      local_pt = right_bound
      NSTEP_local = NSTEP
    case (3)
      ! bottom
      if (myrank == 0) write(IMAIN,*) "calculation of every time step on the bottom absorbing boundary"
      npt = nbot
      allocate(local_pt(npt))
      local_pt = bot_bound
      NSTEP_local = NSTEP
    case default
      call stop_the_code('Invalid flag')
    end select

    ! to distinguish all model case and boundary case
    allocate(temp_field(NSTEP_local))

    allocate(Field_Ux(NFREC1))
    allocate(Field_Uz(NFREC1))
    allocate(Field_Tx(NFREC1))
    allocate(Field_Tz(NFREC1))

! normal vector to the edge at this grid point
! therefore corners between two grid edges must be computed twice
! because the normal will change
    if (FLAG == 1) then
      VNZ = 0.d0
      VNX = 1.d0
    else if (FLAG == 2) then
      VNZ = 0.d0
      VNX = 1.d0
    else if (FLAG == 3) then
      VNZ = 1.d0
      VNX = 0.d0
    else
      VNZ = 0.d0
      VNX = 0.d0
    endif

    do indice = 1,npt

      if (FLAG == 0) then
        inode = indice
        X = coord(1,indice) - offset
! specfem coordinate axes are implemented from bottom to top whereas for this code
! we need from top to bottom
        Z = zmax - coord(2,indice)
      else
        inode = local_pt(indice)
        X = coord(1,inode) - offset
! specfem coordinate axes are implemented from bottom to top whereas for this code
! we need from top to bottom
        Z = zmax - coord(2,inode)
      endif

! user output
      if (myrank == 0 .and. mod(indice,500) == 0) then
        write(IMAIN,*) indice,"points have been computed out of ",npt
      endif

!
! first handle the particular case of zero frequency
!
      TOTO = 0.01d0
      if (source_type == 1) call ONDAS_P(GAMR,0.01d0*BEALF,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)
      if (source_type == 2) call ONDAS_S(GAMR,TOTO,0.01d0*BEALF,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)
      if (source_type == 3) call ONDAS_R(0.01d0*BEALF,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)


      TOTO = 0.0d0
      call DESFXY(TOTO,TOTO,source_type,UX,UZ,SX,SZ,SXZ,A1,B1,A2,B2,AL,AK,AM,RLM)

! write the frequency seismograms
      TX = SX *VNX + SXZ*VNZ
      TZ = SXZ*VNX + SZ *VNZ

      Field_Ux(1) = UX
      Field_Uz(1) = UZ
      if (FLAG /= 0) then
        Field_Tx(1) = TX
        Field_Tz(1) = TZ
      endif

!
! then loop on all the other discrete frequencies
!
      do J = 1,N/2
! compute the value of the frequency (= index * delta in frequency = index * 1/delta in period)
        FJ = dble(J) * 1.d0 / delta_in_period

! pulsation (= 2 * PI * frequency)
        AKA = 2.0d0*PI*FJ

        AQA = AKA*BEALF

! exclude attenuation completely if needed
        if (ATTENUATION_VISCOELASTIC) then
          CAKA = CMPLX(AKA,-AKA/(2.0d0*QD))
          CAQA = CMPLX(AQA,-AQA/(2.0d0*QD))
        else
          CAKA = CMPLX(AKA,0)
          CAQA = CMPLX(AQA,0)
        endif

        if (source_type == 1) then
          call ONDAS_P(GAMR,AQA,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)
        else if (source_type == 2) then
          call ONDAS_S(GAMR,AKA,AQA,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)
        else if (source_type == 3) then
          call ONDAS_R(AQA,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)
        endif

        call DESFXY(X,Z,source_type,UX,UZ,SX,SZ,SXZ,A1,B1,A2,B2,AL,AK,AM,RLM)

! write the frequency seismograms
        TX = SX *VNX + SXZ*VNZ
        TZ = SXZ*VNX + SZ *VNZ

        Field_Ux(J+1) = UX
        Field_Uz(J+1) = UZ
        if (FLAG /= 0) then
          Field_Tx(J+1) = TX
          Field_Tz(J+1) = TZ
        endif
      enddo

! to convert frequency field in time field
! (number at the end are unit numbers for writing in the good file,
! in the case of the traction we fill only one file per call)

! global model case for initial field
      select case (FLAG)
      case (0)
        call paco_convolve_fft(Field_Ux,1,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        displ_elastic(1,indice) = temp_field(1)

        call paco_convolve_fft(Field_Uz,1,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        displ_elastic(2,indice) = temp_field(1)

        call paco_convolve_fft(Field_Ux,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        veloc_elastic(1,indice) = temp_field(1)

        call paco_convolve_fft(Field_Uz,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        veloc_elastic(2,indice) = temp_field(1)

        call paco_convolve_fft(Field_Ux,3,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        accel_elastic(1,indice) = temp_field(1)

        call paco_convolve_fft(Field_Uz,3,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        accel_elastic(2,indice) = temp_field(1)

! absorbing boundaries

! left case
      case (1)
        call paco_convolve_fft(Field_Ux,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        v0x_left(indice,:) = temp_field(:)

        call paco_convolve_fft(Field_Uz,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        v0z_left(indice,:) = temp_field(:)

        call paco_convolve_fft(Field_Tx,4,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        t0x_left(indice,:) = temp_field(:)

        call paco_convolve_fft(Field_Tz,4,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        t0z_left(indice,:) = temp_field(:)

! right case
      case (2)
        call paco_convolve_fft(Field_Ux,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        v0x_right(indice,:) = temp_field(:)

        call paco_convolve_fft(Field_Uz,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        v0z_right(indice,:) = temp_field(:)

        call paco_convolve_fft(Field_Tx,4,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        t0x_right(indice,:) = temp_field(:)

        call paco_convolve_fft(Field_Tz,4,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        t0z_right(indice,:) = temp_field(:)

! bottom case
      case (3)
        call paco_convolve_fft(Field_Ux,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        v0x_bot(indice,:) = temp_field(:)

        call paco_convolve_fft(Field_Uz,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        v0z_bot(indice,:) = temp_field(:)

        call paco_convolve_fft(Field_Tx,4,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        t0x_bot(indice,:) = temp_field(:)

        call paco_convolve_fft(Field_Tz,4,NSTEP_local,dt,NFREC,temp_field,TP,TS)
        t0z_bot(indice,:) = temp_field(:)

      case default
        call stop_the_code('Invalid flag')
      end select
    enddo

    deallocate(temp_field)
    deallocate(local_pt)

    deallocate(Field_Ux)
    deallocate(Field_Uz)
    deallocate(Field_Tx)
    deallocate(Field_Tz)

  enddo

  end subroutine paco_beyond_critical

!
!-------------------------------------------------------------------------------------
!

  subroutine DESFXY(X,Z,ICAS,UX,UZ,SX,SZ,SXZ,A1,B1,A2,B2,AL,AK,AM,RLM)

  implicit none

  double precision,intent(in) :: A1,B1,RLM,X,Z
  integer,intent(in) :: ICAS
  complex(selected_real_kind(15,300)) :: UX,UZ,SX,SZ,SXZ,A2,B2,AL,AK,AM
  complex(selected_real_kind(15,300)) :: UI,FAC
  complex(selected_real_kind(15,300)) :: AUX1,AUX2,FI1,FI2,PS1,PS2

  UI = (0.0d0,1.0d0)
  if (A1 /= 0.0d0) then
     AUX1 = A1*EXP(UI*(AM*Z-AL*X))         ! campo P incidente
  else
     AUX1 = CMPLX(0.0d0)
  endif
  if (A2 /= 0.0d0) then
     AUX2 = A2*EXP(-UI*(AM*Z+AL*X)) *1.0d0      ! campo P reflejado
  else
     AUX2 = CMPLX(0.0d0)
  endif
  FI1 = AUX1+AUX2
  FI2 = AUX1-AUX2
  if (B1 /= 0.0d0) then
     AUX1 = B1*EXP(UI*(AK*Z-AL*X))            ! campo S incidente
  else
     AUX1 = CMPLX(0.0d0)
  endif
  if (B2 /= 0.0d0) then
     AUX2 = B2*EXP(-UI*(AK*Z+AL*X)) *1.0d0      ! campo S reflejado
  else
     AUX2 = CMPLX(0.0d0)
  endif
  PS1 = AUX1+AUX2
  PS2 = AUX1-AUX2

!
!     FAC ES PARA TENER CONSISTENCIA CON AKI & RICHARDS (1980)
!
  FAC = UI
  if (ICAS == 2) FAC = -UI

  UX = (-UI*AL*FI1+UI*AK*PS2)*FAC

  UZ = (UI*AM*FI2+UI*AL*PS1)*FAC
! Paco's convention for vertical coordinate axis is inverted
  UZ = - UZ

  AUX1 = AL*AL+AM*AM
  SX = (-RLM*AUX1*FI1-2.0d0*AL*(AL*FI1-AK*PS2))*FAC
  SZ = (-RLM*AUX1*FI1-2.0d0*(AM*AM*FI1+AK*AL*PS2))*FAC

  SXZ = (2.0d0*AM*AL*FI2+(AL*AL-AK*AK)*PS1)*FAC
! Paco's convention for vertical coordinate axis is inverted
  SXZ = - SXZ

  end subroutine DESFXY

!
!-------------------------------------------------------------------------------------
!

  subroutine FAFB(CA,CB,FA,FB)

  implicit none

  double precision,intent(in) :: CA,CB
  complex(selected_real_kind(15,300)),intent(out) :: FA,FB

  ! local parameters
  double precision :: A,B
  complex(selected_real_kind(15,300)),parameter :: ZER = (0.0d0,0.0d0), UI = (0.0d0,1.0d0)

  A = CA*CA - 1.0d0
  B = CB*CB - 1.0d0

  if (CA < 1.0d0) then
     FA = -UI * SQRT(-A)
  else
     FA = SQRT(A) + ZER
  endif

  if (CB < 1.0d0) then
     FB = -UI * SQRT(-B)
  else
     FB = CMPLX(SQRT(B),0.0d0)
  endif

  end subroutine FAFB

!
!-------------------------------------------------------------------------------------
!

  subroutine A2B2(FA,FB,A2,B2)

  implicit none

  complex(selected_real_kind(15,300)),intent(in) :: FA,FB
  complex(selected_real_kind(15,300)),intent(out) :: A2,B2

  ! local parameters
  complex(selected_real_kind(15,300)) :: DEN,AUX

  AUX = FB*FB - 1.0d0
  DEN = 4.0d0*FA*FB + AUX*AUX

  A2 = (4.0d0*FA*FB-AUX*AUX)/DEN
  B2 = 4.0d0*FA*AUX/DEN

  end subroutine A2B2

!
!-------------------------------------------------------------------------------------
!

! calculation of P waves
  subroutine ONDAS_P(GP,AQB,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)

  implicit none

  double precision,intent(in) :: GP,AQB,ANU
  double precision,intent(out) :: A1,B1,BEALF
  complex(selected_real_kind(15,300)),intent(out) :: A2,B2,AL,AK,AM

  ! local parameters
  double precision :: CA,CB
  complex(selected_real_kind(15,300)) :: FA,FB
  complex(selected_real_kind(15,300)),parameter :: ZER = (0.0d0,0.0d0)

  BEALF = SQRT((1.0d0-2.0d0*ANU)/2.0d0/(1.0d0-ANU))
  A1 = 1.0d0/AQB
  B1 = 0.0d0

  if (GP == 0.0d0) then
    AL = ZER
    AK = ZER
    AM = AQB+ZER
    A2 = (-1.0d0+ZER)/AQB
    B2 = ZER
    return
  endif

  CA = 1.0d0/SIN(GP)
  CB = CA/BEALF
  AL = AQB/CA+ZER
  CALL FAFB(CA,CB,FA,FB)
  AK = AL*FB
  AM = AL*FA
  CALL A2B2(FA,FB,A2,B2)
  A2 = A2/AQB
  B2 = B2/AQB

  end subroutine ONDAS_P

!
!-------------------------------------------------------------------------------------
!

! calculation of S waves
  subroutine ONDAS_S(GS,AKB,AQB,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)

  implicit none

  double precision,intent(in) :: GS,AKB,AQB,ANU
  double precision,intent(out) :: A1,B1,BEALF
  complex(selected_real_kind(15,300)),intent(out) :: A2,B2,AL,AK,AM

  ! local parameters
  double precision :: CA,CB
  complex(selected_real_kind(15,300)) :: FA,FB
  complex(selected_real_kind(15,300)),parameter :: ZER = (0.0d0,0.0d0)

  BEALF = SQRT((1.0d0-2.0d0*ANU)/2.0d0/(1.0d0-ANU))
  A1 = 0.0d0
  B1 = 1.0d0/AKB

  if (GS == 0.0d0) then
    AL = ZER
    AK = AKB+ZER
    AM = ZER
    A2 = ZER
    B2 = (-1.0d0+ZER)/AKB
    return
  endif

  CB = 1.0d0/SIN(GS)
  CA = CB*BEALF

!
! case of the critical angle
!
  if (CA == 1.d0) then
    AL = AQB+ZER
    AM = ZER
    call FAFB(CA,CB,FA,FB)
    AK = AL*FB
    B2 = -B1
    A2 = -4.0d0*COS(GS)*B1/(1./BEALF-2.*BEALF)

! case of an angle that is not critical
  else
    AL = AQB/CA+ZER
    call FAFB(CA,CB,FA,FB)
    AK = AL*FB
    AM = AL*FA
    call A2B2(FA,FB,B2,A2)
    A2 = -A2*FB/FA
    A2 = A2/AKB
    B2 = B2/AKB
  endif

  end subroutine ONDAS_S

!
!-------------------------------------------------------------------------------------
!

! calculation of Rayleigh waves
  subroutine ONDAS_R(AQB,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)

  implicit none

  double precision,intent(in) :: AQB,ANU
  double precision,intent(out) :: A1,B1,BEALF
  complex(selected_real_kind(15,300)),intent(out) :: A2,B2,AL,AK,AM

  ! local parameters
  double precision :: CA,CB,BA2
  complex(selected_real_kind(15,300)) :: FA,FB
  complex(selected_real_kind(15,300)),parameter :: ZER = (0.0d0,0.0d0)

  double precision, external :: CRB

  A1 = 0.0d0
  B1 = 0.0d0
  B2 = 1.0d0+ZER
  BEALF = SQRT((1.0d0-2.0d0*ANU)/2.0d0/(1.0d0-ANU))
  BA2 = BEALF*BEALF
  CB = CRB(BEALF)
  CA = CB*BEALF
  AL = AQB/CA+ZER

  CALL FAFB(CA,CB,FA,FB)

  AK = AL*FB
  AM = AL*FA
  A2 = 2.0d0*FB/(FB*FB-1.0d0)*B2
  B2 = B2/(AL*A2+AK)
  A2 = A2*B2

  end subroutine ONDAS_R

!
!-------------------------------------------------------------------------------------
!

  function CRB(BEALF)

  use constants, only: PI

  implicit none

  double precision,intent(in) :: BEALF
  double precision :: CRB

  ! local parameters
  double precision :: BA2,P,Q,FIND,F1,F2,F12,FACT

  double precision,parameter :: U3 = 1.0d0/3.0d0

  BA2 = BEALF*BEALF
  P = 8.0d0/3.0d0 - 16.0d0*BA2
  Q = 272.0d0/27.0d0 - 80.0d0/3.0d0*BA2

  FIND = Q*Q/4.0d0 + P*P*P/27.0d0

  if (FIND >= 0.0d0) then
    F1 = SQRT(FIND)-Q/2.0d0
    if (F1 > 0.0d0) then
      F1 = F1**U3
    else
      F1 = -(-F1)**U3
    endif
    F2 = -SQRT(FIND)-Q/2.0d0
    if (F2 > 0.0d0) then
      F2 = F2**U3
    else
      F2 = -(-F2)**U3
    endif
    FACT = F1+F2+8.0d0/3.0d0
    CRB = SQRT(FACT)
  else
    F1 = -27.0d0*Q*Q/(4.0d0*P*P*P)
    F1 = SQRT(F1)
    if (Q < 0.0d0) then
      F1 = COS((PI-ACOS(F1))/3.0d0)
    else
      F1 = COS(ACOS(F1)/3.0d0)
    endif
    F2 = -P/3.0d0
    F2 = SQRT(F2)
    F12 = -2.0d0*F1*F2+8.0d0/3.0d0
    CRB = SQRT(F12)
  endif

  end function CRB

