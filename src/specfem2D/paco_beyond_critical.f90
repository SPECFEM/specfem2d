!
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

subroutine paco_beyond_critical(anglesource,f0,QD,source_type,left_bound,right_bound,&
                                bot_bound,nleft,nright,nbot,x_source)

  use specfem_par,only : coord,nglob,deltat,NSTEP,cploc,csloc,ATTENUATION_VISCOELASTIC_SOLID,v0x_left,&
                         v0z_left,v0x_right,v0z_right,v0x_bot,v0z_bot,t0x_left,t0z_left,t0x_right,t0z_right,t0x_bot,t0z_bot,&
                         displ_elastic,veloc_elastic,accel_elastic


  implicit none

  include "constants.h"

  double precision :: f0,dt,TP,anglesource,QD,delta_in_period
  integer :: npt,source_type,nleft,nright,nbot

  integer, dimension(nleft) :: left_bound
  integer, dimension(nright) :: right_bound
  integer, dimension(nbot) :: bot_bound

  integer, dimension(:),allocatable :: local_pt

  double precision, dimension(:), allocatable :: temp_field

  integer :: J, indice, NSTEP_local, FLAG, N, NFREC, NFREC1

  double precision :: ANU,BEALF,ALFBE,RLM,VNX,VNZ,A1,B1,TOTO,FJ,AKA,AQA,GAMR

! location of the point
  double precision :: X, Z, xmin, xmax, zmin, zmax
  integer :: inode

  complex(selected_real_kind(15,300)) :: CAKA,CAQA,UI,UR
  complex(selected_real_kind(15,300)) :: UX,UZ,SX,SZ,SXZ,A2,B2,AL,AK,AM

  complex(selected_real_kind(15,300)) :: TX,TZ

  complex(selected_real_kind(15,300)), dimension(:),allocatable::Field_Ux,Field_Uz,Field_Tx,Field_Tz

  double precision :: TS,offset,x_source

! size of the model
  xmin=minval(coord(1,:))
  xmax=maxval(coord(1,:))
  zmin=minval(coord(2,:))
  zmax=maxval(coord(2,:))

  TS=1.2d0/f0

! dominant period of the Ricker
  TP=1.d0/f0

! offset to move the initial location of the source in the horizontal direction of the mesh
  offset = x_source

! find optimal period
! if period is too small, you should see several initial plane wave on your initial field
  delta_in_period=2.d0
  do while(delta_in_period<1.5*abs(xmax-xmin)/csloc)
     delta_in_period=2.d0*delta_in_period
  enddo

! test Deltat compatibility
  DT=256.d0
  do while(DT>deltat)
     DT=DT/2.d0
  enddo
  if (abs(DT-deltat)>1.0d-13) then
     print *, "you must take a deltat that is a power of two (power can be negative)"
     print *, "for example you can take", DT
     stop "cannot go further, restart with new deltat"
  endif

  DT=deltat/2.d0

  N=2
  do while(N<2*NSTEP+1)
     N=2*N
  enddo

  do while(DT<(delta_in_period/N))
     N=2*N
  enddo

  print *,'N found to perform the frequency calculation:',N
  print *,'number of discrete frequencies = ',N/2
  print *,'delta in period (seconds) = ',delta_in_period
  print *,'delta in frequency (Hz) = ',1.d0/delta_in_period
  print *,'dt (here we need deltat/2) = ', DT

  NFREC=N/2
  NFREC1=NFREC+1


!
!     FDT:  FUNCION DE TRANSFERENCIA
!

! calculation of Poisson's ratio
  ANU = (cploc*cploc-2.d0*csloc*csloc)/(2.d0*(cploc*cploc-csloc*csloc))
  print *,"Poisson's ratio = ",ANU

  UI=(0.0d0, 1.0d0)
  UR=(1.0d0, 0.0d0)

! convert angle to radians
  GAMR = anglesource

  BEALF=SQRT((1.0d0-2.0d0*ANU)/(2.0d0*(1.0d0-ANU)))
  ALFBE=1.0d0/BEALF
  RLM=ALFBE**2-2.0d0

! flags: interior=0, left=1, right=2, bottom=3
  do FLAG=0,3

     if (FLAG==0) then
        print *,"calculation of the initial field for every point of the mesh"
        npt=nglob
        allocate(local_pt(npt))
        do inode=1,npt
           local_pt(inode)=inode
        enddo
        NSTEP_local=1
     else if(FLAG==1) then
        print *,"calculation of every time step on the left absorbing boundary"
        npt=nleft
        allocate(local_pt(npt))
        local_pt=left_bound
        NSTEP_local=NSTEP
     else if(FLAG==2) then
        print *,"calculation of every time step on the right absorbing boundary"
        npt=nright
        allocate(local_pt(npt))
        local_pt=right_bound
        NSTEP_local=NSTEP
     else if(FLAG==3) then
        print *,"calculation of every time step on the bottom absorbing boundary"
        npt=nbot
        allocate(local_pt(npt))
        local_pt=bot_bound
        NSTEP_local=NSTEP
     endif

! to distinguish all model case and boundary case
     allocate(temp_field(NSTEP_local))

     allocate(Field_Ux(NFREC1))
     allocate(Field_Uz(NFREC1))
     allocate(Field_Tx(NFREC1))
     allocate(Field_Tz(NFREC1))


     if(mod(N,2) /= 0) stop 'N must be a multiple of 2'

! normal vector to the edge at this grid point
! therefore corners between two grid edges must be computed twice
! because the normal will change
     if (FLAG==1) then
        VNZ = 0.d0
        VNX = 1.d0
     else if (FLAG==2) then
        VNZ = 0.d0
        VNX = 1.d0
     else if (FLAG==3) then
        VNZ = 1.d0
        VNX = 0.d0
     else
        VNZ = 0.d0
        VNX = 0.d0
     endif


     do indice=1,npt

        if (FLAG==0) then
           inode=indice
           X=coord(1,indice) - offset
! specfem coordinate axes are implemented from bottom to top whereas for this code
! we need from top to bottom
           Z=zmax-coord(2,indice)
        else
           inode=local_pt(indice)
           X=coord(1,inode) - offset
! specfem coordinate axes are implemented from bottom to top whereas for this code
! we need from top to bottom
           Z=zmax-coord(2,inode)
        endif

        if (mod(indice,500)==0) then
           print *,indice,"points have been computed out of ",npt
        endif

!
! first handle the particular case of zero frequency
!
        TOTO=0.01d0
        IF (source_type==1) CALL ONDASP(GAMR,0.01d0*BEALF,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)
        IF (source_type==2) CALL ONDASS(GAMR,TOTO,0.01d0*BEALF,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)
        IF (source_type==3) CALL ONDASR(0.01d0*BEALF,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)


        TOTO=0.0d0
        CALL DESFXY(TOTO,TOTO,source_type,UX,UZ,SX,SZ,SXZ,A1,B1,A2,B2,AL,AK,AM,RLM)

! write the frequency seismograms
        TX = SX *VNX+SXZ*VNZ
        TZ = SXZ*VNX+SZ *VNZ

        Field_Ux(1)=UX
        Field_Uz(1)=UZ
        if (FLAG/=0) then
           Field_Tx(1)=TX
           Field_Tz(1)=TZ
        endif

!
! then loop on all the other discrete frequencies
!
        do J=1,N/2

! compute the value of the frequency (= index * delta in frequency = index * 1/delta in period)
           FJ = dble(J) * 1.d0 / delta_in_period

! pulsation (= 2 * PI * frequency)
           AKA=2.0d0*PI*FJ

           AQA=AKA*BEALF

! exclude attenuation completely if needed
           if(ATTENUATION_VISCOELASTIC_SOLID) then
              CAKA=CMPLX(AKA,-AKA/(2.0d0*QD))
              CAQA=CMPLX(AQA,-AQA/(2.0d0*QD))
           else
              CAKA=CMPLX(AKA,0)
              CAQA=CMPLX(AQA,0)
           endif

           IF (source_type==1) CALL ONDASP(GAMR,AQA,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)
           IF (source_type==2) CALL ONDASS(GAMR,AKA,AQA,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)
           IF (source_type==3) CALL ONDASR(AQA,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)

           CALL DESFXY(X,Z,source_type,UX,UZ,SX,SZ,SXZ,A1,B1,A2,B2,AL,AK,AM,RLM)

! write the frequency seismograms
           TX = SX *VNX+SXZ*VNZ
           TZ = SXZ*VNX+SZ *VNZ

           Field_Ux(J+1)=UX
           Field_Uz(J+1)=UZ
           if (FLAG/=0) then
              Field_Tx(J+1)=TX
              Field_Tz(J+1)=TZ
           endif

        enddo

! to convert frequency field in time field
! (number at the end are unit numbers for writing in the good file,
! in the case of the traction we fill only one file per call)

! global model case for initial field
        if (FLAG==0) then
           call paco_convolve_fft(Field_Ux,1,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           displ_elastic(1,indice)=temp_field(1)
           call paco_convolve_fft(Field_Uz,1,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           displ_elastic(3,indice)=temp_field(1)
           call paco_convolve_fft(Field_Ux,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           veloc_elastic(1,indice)=temp_field(1)
           call paco_convolve_fft(Field_Uz,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           veloc_elastic(3,indice)=temp_field(1)
           call paco_convolve_fft(Field_Ux,3,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           accel_elastic(1,indice)=temp_field(1)
           call paco_convolve_fft(Field_Uz,3,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           accel_elastic(3,indice)=temp_field(1)

! absorbing boundaries

! left case
        else if (FLAG==1) then
           call paco_convolve_fft(Field_Ux,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           v0x_left(indice,:)=temp_field(:)
           call paco_convolve_fft(Field_Uz,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           v0z_left(indice,:)=temp_field(:)
           call paco_convolve_fft(Field_Tx,4,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           t0x_left(indice,:)=temp_field(:)
           call paco_convolve_fft(Field_Tz,4,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           t0z_left(indice,:)=temp_field(:)

! right case
        else if (FLAG==2) then
           call paco_convolve_fft(Field_Ux,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           v0x_right(indice,:)=temp_field(:)
           call paco_convolve_fft(Field_Uz,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           v0z_right(indice,:)=temp_field(:)
           call paco_convolve_fft(Field_Tx,4,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           t0x_right(indice,:)=temp_field(:)
           call paco_convolve_fft(Field_Tz,4,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           t0z_right(indice,:)=temp_field(:)

! bottom case
        else if (FLAG==3) then
           call paco_convolve_fft(Field_Ux,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           v0x_bot(indice,:)=temp_field(:)
           call paco_convolve_fft(Field_Uz,2,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           v0z_bot(indice,:)=temp_field(:)
           call paco_convolve_fft(Field_Tx,4,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           t0x_bot(indice,:)=temp_field(:)
           call paco_convolve_fft(Field_Tz,4,NSTEP_local,dt,NFREC,temp_field,TP,TS)
           t0z_bot(indice,:)=temp_field(:)
        endif
     enddo

     deallocate(temp_field)
     deallocate(local_pt)

     deallocate(Field_Ux)
     deallocate(Field_Uz)
     deallocate(Field_Tx)
     deallocate(Field_Tz)

  enddo

end subroutine paco_beyond_critical

!---

SUBROUTINE DESFXY(X,Z,ICAS,UX,UZ,SX,SZ,SXZ,A1,B1,A2,B2,AL,AK,AM,RLM)

  implicit none

  double precision A1,B1,RLM,X,Z
  integer ICAS
  complex(selected_real_kind(15,300)) :: UX,UZ,SX,SZ,SXZ,A2,B2,AL,AK,AM
  complex(selected_real_kind(15,300)) :: UI,FAC
  complex(selected_real_kind(15,300)) :: AUX1,AUX2,FI1,FI2,PS1,PS2

  UI=(0.0d0,1.0d0)
  if (A1/=0.0d0) then
     AUX1=A1*EXP(UI*(AM*Z-AL*X))         ! campo P incidente
  else
     AUX1=CMPLX(0.0d0)
  endif
  if (A2/=0.0d0) then
     AUX2=A2*EXP(-UI*(AM*Z+AL*X)) *1.0d0      ! campo P reflejado
  else
     AUX2=CMPLX(0.0d0)
  endif
  FI1=AUX1+AUX2
  FI2=AUX1-AUX2
  if (B1/=0.0d0) then
     AUX1=B1*EXP(UI*(AK*Z-AL*X))            ! campo S incidente
  else
     AUX1=CMPLX(0.0d0)
  endif
  if (B2/=0.0d0) then
     AUX2=B2*EXP(-UI*(AK*Z+AL*X)) *1.0d0      ! campo S reflejado
  else
     AUX2=CMPLX(0.0d0)
  endif
  PS1=AUX1+AUX2
  PS2=AUX1-AUX2

!
!     FAC ES PARA TENER CONSISTENCIA CON AKI & RICHARDS (1980)
!
  FAC=UI
  IF (ICAS==2)FAC=-UI

  UX=(-UI*AL*FI1+UI*AK*PS2)*FAC

  UZ=(UI*AM*FI2+UI*AL*PS1)*FAC
! Paco's convention for vertical coordinate axis is inverted
  UZ = - UZ

  AUX1=AL*AL+AM*AM
  SX=(-RLM*AUX1*FI1-2.0d0*AL*(AL*FI1-AK*PS2))*FAC
  SZ=(-RLM*AUX1*FI1-2.0d0*(AM*AM*FI1+AK*AL*PS2))*FAC

  SXZ=(2.0d0*AM*AL*FI2+(AL*AL-AK*AK)*PS1)*FAC
! Paco's convention for vertical coordinate axis is inverted
  SXZ = - SXZ

END SUBROUTINE DESFXY

SUBROUTINE FAFB(CA,CB,FA,FB)

  implicit none

  double precision CA,CB,A,B
  complex(selected_real_kind(15,300)) :: FA,FB,ZER,UI

  ZER=(0.0d0,0.0d0)
  UI=(0.0d0,1.0d0)
  A=CA*CA-1.0d0
  B=CB*CB-1.0d0

  IF (CA<1.0d0) then
     FA=-UI*SQRT(-A)
  else
     FA=SQRT(A)+ZER
  endif

  IF (CB<1.0d0) then
     FB=-UI*SQRT(-B)
  else
     FB=CMPLX(SQRT(B),0.0d0)
  endif

END SUBROUTINE FAFB

SUBROUTINE A2B2(FA,FB,A2,B2)

  implicit none

  complex(selected_real_kind(15,300)) :: FA,FB,A2,B2,DEN,AUX

  AUX=FB*FB-1.0d0
  DEN=4.0d0*FA*FB+AUX*AUX
  A2=(4.0d0*FA*FB-AUX*AUX)/DEN
  B2=4.0d0*FA*AUX/DEN

END SUBROUTINE A2B2

! calculation of P waves
SUBROUTINE ONDASP(GP,AQB,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)

  implicit none

  double precision A1,B1,ANU,CA,CB,GP,AQB,BEALF
  complex(selected_real_kind(15,300)) :: A2,B2,FA,FB,ZER,AL,AK,AM

  ZER=(0.0d0,0.0d0)
  BEALF=SQRT((1.0d0-2.0d0*ANU)/2.0d0/(1.0d0-ANU))
  A1=1.0d0/AQB
  B1=0.0d0

  IF (GP==0.0d0) then
     AL=ZER
     AK=ZER
     AM=AQB+ZER
     A2=(-1.0d0+ZER)/AQB
     B2=ZER
     RETURN
  endif

  CA=1.0d0/SIN(GP)
  CB=CA/BEALF
  AL=AQB/CA+ZER
  CALL FAFB(CA,CB,FA,FB)
  AK=AL*FB
  AM=AL*FA
  CALL A2B2(FA,FB,A2,B2)
  A2=A2/AQB
  B2=B2/AQB

END SUBROUTINE ONDASP

! calculation of S waves
SUBROUTINE ONDASS(GS,AKB,AQB,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)

  implicit none

  double precision A1,B1,ANU,CA,CB,GS,AQB,BEALF,AKB
  complex(selected_real_kind(15,300)) :: A2,B2,FA,FB,ZER,AL,AK,AM

  ZER=(0.0d0,0.0d0)
  BEALF=SQRT((1.0d0-2.0d0*ANU)/2.0d0/(1.0d0-ANU))
  A1=0.0d0
  B1=1.0d0/AKB

  IF (GS==0.0d0) then
     AL=ZER
     AK=AKB+ZER
     AM=ZER
     A2=ZER
     B2=(-1.0d0+ZER)/AKB
     return
  endif

  CB=1.0d0/SIN(GS)
  CA=CB*BEALF

!
! case of the critical angle
!
  IF (CA==1.d0) then
    AL=AQB+ZER
    AM=ZER
    CALL FAFB(CA,CB,FA,FB)
    AK=AL*FB
    B2=-B1
    A2=-4.0d0*COS(GS)*B1/(1./BEALF-2.*BEALF)

! case of an angle that is not critical
  ELSE
    AL=AQB/CA+ZER
    CALL FAFB(CA,CB,FA,FB)
    AK=AL*FB
    AM=AL*FA
    CALL A2B2(FA,FB,B2,A2)
    A2=-A2*FB/FA
    A2=A2/AKB
    B2=B2/AKB
  endif

END SUBROUTINE ONDASS

! calculation of Rayleigh waves
SUBROUTINE ONDASR(AQB,A1,B1,A2,B2,AL,AK,AM,ANU,BEALF)

  implicit none

  double precision A1,B1,ANU,CA,CB,AQB,BEALF,ba2
  complex(selected_real_kind(15,300)) :: A2,B2,FA,FB,ZER,AL,AK,AM

  double precision, external :: crb

  ZER=(0.0d0,0.0d0)
  A1=0.0d0
  B1=0.0d0
  B2=1.0d0+ZER
  BEALF=SQRT((1.0d0-2.0d0*ANU)/2.0d0/(1.0d0-ANU))
  BA2=BEALF*BEALF
  CB=CRB(BEALF)
  CA=CB*BEALF
  AL=AQB/CA+ZER

  CALL FAFB(CA,CB,FA,FB)

  AK=AL*FB
  AM=AL*FA
  A2=2.0d0*FB/(FB*FB-1.0d0)*B2
  B2=B2/(AL*A2+AK)
  A2=A2*B2

END SUBROUTINE ONDASR

FUNCTION CRB(BEALF)

  implicit none

  include "constants.h"

  double precision U3,BA2,P,Q,FIND,F1,F2,F12,FACT,CRB,BEALF

  U3=1.0d0/3.0d0
  BA2=BEALF*BEALF
  P=8.0d0/3.0d0-16.0d0*BA2
  Q=272.0d0/27.0d0-80.0d0/3.0d0*BA2
  FIND=Q*Q/4.0d0+P*P*P/27.0d0
  IF (FIND>=0.0d0) then
     F1=SQRT(FIND)-Q/2.0d0
     IF (F1>0.0d0) then
        F1=F1**U3
     else
        F1=-(-F1)**U3
     endif
     F2=-SQRT(FIND)-Q/2.0d0
     IF (F2>0.0d0) then
        F2=F2**U3
     else
        F2=-(-F2)**U3
     endif
     FACT=F1+F2+8.0d0/3.0d0
     CRB=SQRT(FACT)
  else
     F1=-27.0d0*Q*Q/(4.0d0*P*P*P)
     F1=SQRT(F1)
     IF (Q<0.0d0) then
        F1=COS((PI-ACOS(F1))/3.0d0)
     else
        F1=COS(ACOS(F1)/3.0d0)
     endif
     F2=-P/3.0d0
     F2=SQRT(F2)
     F12=-2.0d0*F1*F2+8.0d0/3.0d0
     CRB=SQRT(F12)
  endif

END FUNCTION CRB

