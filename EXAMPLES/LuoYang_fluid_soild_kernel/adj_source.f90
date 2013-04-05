implicit none

double precision :: time,pressure
integer,parameter :: NSTEP = 1600
integer :: i

open(1,file="OUTPUT_FILES/S0001.AA.PRE.semp",status="old")

open(10,file="SEM/S0001.AA.BXX.adj",status="unknown")

open(20,file="SEM/S0001.AA.BXY.adj",status="unknown")

open(30,file="SEM/S0001.AA.BXZ.adj",status="unknown")

do i=1,NSTEP

  read(1,*)time,pressure
  if(time>1.274)pressure=0.d0
  write(10,*)time,pressure
  write(20,*)time,pressure
  write(30,*)time,pressure

enddo

end
