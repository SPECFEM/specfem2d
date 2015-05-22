implicit none

double precision :: time,displ_potential
integer,parameter :: NSTEP = 1600
integer :: i

open(1,file="OUTPUT_FILES/AA.S0001.PRE.semp",status="old")

open(10,file="SEM/AA.S0001.POT.adj",status="unknown")

do i=1,NSTEP

  read(1,*)time,displ_potential
  if(time>1.274)displ_potential=0.d0
  write(10,*)time,displ_potential

enddo

close(10)

end
